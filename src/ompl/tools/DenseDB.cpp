/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2016, University of Colorado, Boulder
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Univ of CO, Boulder nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *********************************************************************/

/* Author: Dave Coleman <dave@dav.ee>
   Desc:   Experience database for storing and reusing past path plans
*/

// OMPL
#include <ompl/tools/bolt/DenseDB.h>
#include <ompl/tools/bolt/Discretizer.h>
#include <ompl/base/ScopedState.h>
#include <ompl/util/Time.h>
#include <ompl/util/Console.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>  // TODO: remove, this is not space agnostic

// Boost
#include <boost/filesystem.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

// C++
#include <limits>
#include <queue>

#define foreach BOOST_FOREACH

namespace og = ompl::geometric;
namespace ot = ompl::tools;
namespace ob = ompl::base;
namespace otb = ompl::tools::bolt;

namespace boost
{
double get(const ompl::tools::bolt::DenseEdgeWeightMap &m, const ompl::tools::bolt::DenseEdge &e)
{
  return m.get(e);
}
}

BOOST_CONCEPT_ASSERT(
    (boost::ReadablePropertyMapConcept<ompl::tools::bolt::DenseEdgeWeightMap, ompl::tools::bolt::DenseEdge>));

// CustomAstarVisitor methods ////////////////////////////////////////////////////////////////////////////

BOOST_CONCEPT_ASSERT((boost::AStarVisitorConcept<otb::DenseDB::CustomAstarVisitor, otb::DenseGraph>));

otb::DenseDB::CustomAstarVisitor::CustomAstarVisitor(DenseVertex goal, DenseDB *parent) : goal_(goal), parent_(parent)
{
}

void otb::DenseDB::CustomAstarVisitor::discover_vertex(DenseVertex v, const DenseGraph &) const
{
  if (parent_->visualizeAstar_)
    parent_->getVisual()->viz4State(parent_->stateProperty_[v], tools::SMALL, tools::GREEN, 1);
}

void otb::DenseDB::CustomAstarVisitor::examine_vertex(DenseVertex v, const DenseGraph &) const
{
  if (parent_->visualizeAstar_)
  {
    parent_->getVisual()->viz4State(parent_->stateProperty_[v], tools::MEDIUM, tools::BLACK, 1);
    parent_->getVisual()->viz4Trigger();
    usleep(parent_->visualizeAstarSpeed_ * 1000000);
  }

  if (v == goal_)
    throw FoundGoalException();
}

// Actual class ////////////////////////////////////////////////////////////////////////////

namespace ompl
{
namespace tools
{
namespace bolt
{
DenseDB::DenseDB(base::SpaceInformationPtr si, VisualizerPtr visual)
  : si_(si)
  , visual_(visual)
  // Property accessors of edges
  , edgeWeightProperty_(boost::get(boost::edge_weight, g_))
  , edgeCollisionStateProperty_(boost::get(edge_collision_state_t(), g_))
  // Property accessors of vertices
  , stateProperty_(boost::get(vertex_state_t(), g_))
  , typeProperty_(boost::get(vertex_type_t(), g_))
  , representativesProperty_(boost::get(vertex_sparse_rep_t(), g_))
  // Disjoint set accessors
  , disjointSets_(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_))
{
  // Add search state
  initializeQueryState();

  // Initialize collision cache
  denseCache_.reset(new DenseCache(si_, this, visual_));

  // Initialize sparse database
  sparseDB_.reset(new SparseDB(si_, this, visual_, denseCache_));

  // Initialize nearest neighbor datastructure
  nn_.reset(new NearestNeighborsGNAT<DenseVertex>());
  nn_->setDistanceFunction(boost::bind(&DenseDB::distanceFunction, this, _1, _2));

  // Initialize the discretize grid tool
  discretizer_.reset(new Discretizer(si_, this, denseCache_, visual_));
}

DenseDB::~DenseDB(void)
{
  if (graphUnsaved_)
    std::cout << "The database is being unloaded with unsaved experiences" << std::endl;
  freeMemory();
}

void DenseDB::freeMemory()
{
  foreach (DenseVertex v, boost::vertices(g_))
  {
    if (stateProperty_[v] != nullptr)
      si_->freeState(stateProperty_[v]);
    stateProperty_[v] = nullptr;  // TODO(davetcoleman): is this needed??
  }

  g_.clear();

  if (nn_)
    nn_->clear();

  sampler_.reset();
}

bool DenseDB::setup()
{
  if (!sampler_)
    sampler_ = si_->allocValidStateSampler();

  sparseDB_->setup();

  return true;
}

bool DenseDB::load()
{
  OMPL_INFORM("DenseDB: load()");

  // Error checking
  if (getNumEdges() > queryVertices_.size() ||
      getNumVertices() > queryVertices_.size())  // the search verticie may already be there
  {
    OMPL_INFORM("Database is not empty, unable to load from file");
    return true;
  }
  if (filePath_.empty())
  {
    OMPL_ERROR("Empty filename passed to save function");
    return false;
  }
  if (!boost::filesystem::exists(filePath_))
  {
    OMPL_INFORM("Database file does not exist: %s.", filePath_.c_str());
    return false;
  }

  // Benchmark
  time::point start = time::now();

  // Load
  OMPL_INFORM("Loading database from file: %s", filePath_.c_str());

  BoltStorage storage_(si_, this);
  storage_.load(filePath_.c_str());

  // Benchmark
  double duration = time::seconds(time::now() - start);

  // Load collision cache
  denseCache_->load();

  // Visualize
  // visual_->viz1Trigger();
  // usleep(0.1 * 1000000);

  // Error check
  if (!getNumVertices() || !getNumEdges())
  {
    OMPL_ERROR("Corrupted planner data loaded, skipping building graph");
    return false;
  }

  // Get the average vertex degree (number of connected edges)
  double averageDegree = (getNumEdges() * 2) / static_cast<double>(getNumVertices());

  // Check how many disjoint sets are in the dense graph (should be none)
  std::size_t numSets = checkConnectedComponents();

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Loaded graph stats:");
  OMPL_INFORM("   Total valid vertices:   %u", getNumVertices());
  OMPL_INFORM("   Total valid edges:      %u", getNumEdges());
  OMPL_INFORM("   Average degree:         %f", averageDegree);
  OMPL_INFORM("   Connected Components:   %u", numSets);
  OMPL_INFORM("   Loading time:           %f", duration);
  OMPL_INFORM("------------------------------------------------------");

  return true;
}

bool DenseDB::saveIfChanged()
{
  if (graphUnsaved_)
  {
    return save();
  }
  else
    OMPL_INFORM("Not saving because database has not changed");
  return true;
}

bool DenseDB::save()
{
  if (!graphUnsaved_)
    OMPL_WARN("No need to save because graphUnsaved_ is false, but saving anyway because requested");

  // Disabled
  if (!savingEnabled_)
  {
    OMPL_INFORM("Not saving because option disabled for DenseDB");
    return false;
  }

  // Error checking
  if (filePath_.empty())
  {
    OMPL_ERROR("Empty filename passed to save function");
    return false;
  }

  OMPL_INFORM("Saving with %d vertices and %d edges to: %s", getNumVertices(), getNumEdges(), filePath_.c_str());

  // Benchmark
  time::point start = time::now();

  // Save
  BoltStorage storage_(si_, this);
  storage_.save(filePath_.c_str());

  // Save collision cache
  denseCache_->save();

  // Benchmark
  double loadTime = time::seconds(time::now() - start);
  OMPL_INFORM("Saved database to file in %f sec", loadTime);

  graphUnsaved_ = false;
  return true;
}

bool DenseDB::generateGrid()
{
  std::cout << "generateGrid called ------------- " << std::endl;
  return discretizer_->generateGrid();
}

bool DenseDB::postProcessPath(og::PathGeometric &solutionPath)
{
  // Prevent inserting into database
  if (!savingEnabled_)
  {
    OMPL_WARN("DenseDB: Saving is disabled so not adding path");
    return false;
  }

  if (visualizeSnapPath_)  // Clear old path
  {
    visual_->viz5DeleteAllMarkers();
    visual_->viz4DeleteAllMarkers();
  }

  // Get starting state
  base::State *currentPathState = solutionPath.getStates()[0];

  // Get neighbors
  std::vector<DenseVertex> graphNeighborhood;
  std::vector<DenseVertex> visibleNeighborhood;
  std::size_t coutIndent = 0;
  const std::size_t numThreads = 0;
  findGraphNeighbors(currentPathState, graphNeighborhood, visibleNeighborhood, sparseDB_->sparseDelta_, numThreads,
                     coutIndent);

  std::vector<DenseVertex> roadmapPath;

  // Run in non-debug mode
  bool recurseVerbose = snapPathVerbose_;
  if (!postProcessPathWithNeighbors(solutionPath, visibleNeighborhood, recurseVerbose, roadmapPath))
  {
    OMPL_ERROR("Could not find snap waypoint path. Running again in debug");
    std::cout << "-------------------------------------------------------" << std::endl;

    // Run in debug mode
    recurseVerbose = true;
    visualizeSnapPath_ = true;
    roadmapPath.clear();
    postProcessPathWithNeighbors(solutionPath, visibleNeighborhood, recurseVerbose, roadmapPath);

    // temp
    std::cout << "exiting for debug " << std::endl;
    exit(-1);
  }

  // Error check
  if (roadmapPath.size() < 2)
  {
    OMPL_WARN("Snapped path waypoint count is too short, only contains %u waypoints", roadmapPath.size());
    if (roadmapPath.empty())
    {
      OMPL_ERROR("Trajectory completely empty!");
      exit(-1);
    }
    // It is possible to have a path of [actualStart, middlePoint, actualGoal], in which case we can't save any
    // experience from it
  }

  if (roadmapPath.size() > 100)
    OMPL_WARN("Roadmap size is %u", roadmapPath.size());

  if (snapPathVerbose_)
    std::cout << "Finished recurseSnapWaypoints(), now updating edge weights in Dense graph " << std::endl;

  // Update edge weights based on this newly created path
  updateEdgeWeights(roadmapPath);

  // Record this new addition
  graphUnsaved_ = true;

  return true;
}

bool DenseDB::postProcessPathWithNeighbors(og::PathGeometric &solutionPath,
                                           const std::vector<DenseVertex> &visibleNeighborhood, bool recurseVerbose,
                                           std::vector<DenseVertex> &roadmapPath)
{
  std::size_t currVertexIndex = 1;

  // Remember if any connections failed
  bool allValid = true;

  for (std::size_t i = 0; i < visibleNeighborhood.size(); ++i)
  {
    if (recurseVerbose)
      std::cout << "Attempting to start with neighbor " << i << std::endl;
    DenseVertex prevGraphVertex = visibleNeighborhood[i];

    if (visualizeSnapPath_)  // Add first state
    {
      visual_->viz5State(stateProperty_[prevGraphVertex], tools::SMALL, tools::GREEN, 1);
    }

    // Add this start state
    roadmapPath.push_back(prevGraphVertex);

    // Start recursive function
    allValid = true;
    if (!recurseSnapWaypoints(solutionPath, roadmapPath, currVertexIndex, prevGraphVertex, allValid, recurseVerbose))
    {
      std::cout << "Failed to find path with starting state neighbor " << i << std::endl;
    }
    else
    {
      break;  // sucess
    }

    if (visualizeSnapPath_)  // Visualize
    {
      visual_->viz5Trigger();
      usleep(visualizeSnapPathSpeed_ * 1000000);
    }
  }

  return allValid;
}

bool DenseDB::updateEdgeWeights(const std::vector<DenseVertex> &roadmapPath)
{
  for (std::size_t vertexID = 1; vertexID < roadmapPath.size(); ++vertexID)
  {
    std::pair<DenseEdge, bool> edgeResult = boost::edge(roadmapPath[vertexID - 1], roadmapPath[vertexID], g_);
    DenseEdge &edge = edgeResult.first;

    // Error check
    if (!edgeResult.second)
    {
      std::cout << std::string(2, ' ') << "WARNING: No edge found on snapped path at index " << vertexID
                << ", unable to save popularity of this edge. perhaps path needs interpolation first" << std::endl;

      if (visualizeSnapPath_)  // Visualize
      {
        const double cost = 100;  // red
        visual_->viz4Edge(stateProperty_[roadmapPath[vertexID - 1]], stateProperty_[roadmapPath[vertexID]], cost);
        visual_->viz4Trigger();
        usleep(visualizeSnapPathSpeed_ * 1000000);
      }
      std::cout << "shutting down out of curiosity " << std::endl;
      exit(-1);
    }
    else
    {
      // reduce cost of this edge because it was just used (increase popularity)
      // Note: 100 is an *unpopular* edge, and 0 is a super highway
      if (snapPathVerbose_)
      {
        std::cout << "Edge weight for vertex " << vertexID << " of edge " << edge << std::endl;
        std::cout << "    old: " << edgeWeightProperty_[edge];
      }
      edgeWeightProperty_[edge] = std::max(edgeWeightProperty_[edge] - POPULARITY_WEIGHT_REDUCTION, 0.0);
      if (snapPathVerbose_)
        std::cout << " new: " << edgeWeightProperty_[edge] << std::endl;

      if (visualizeSnapPath_)  // Visualize
      {
        visual_->viz5Edge(stateProperty_[roadmapPath[vertexID - 1]], stateProperty_[roadmapPath[vertexID]], 100);
      }
    }
  }

  if (visualizeSnapPath_)  // Visualize
  {
    visual_->viz5Trigger();
    usleep(visualizeSnapPathSpeed_ * 1000000);
  }

  return true;
}

bool DenseDB::recurseSnapWaypoints(og::PathGeometric &inputPath, std::vector<DenseVertex> &roadmapPath,
                                   std::size_t currVertexIndex, const DenseVertex &prevGraphVertex, bool &allValid,
                                   bool verbose)
{
  if (verbose)
    std::cout << std::string(currVertexIndex, ' ') << "recurseSnapWaypoints() -------" << std::endl;

  // Find multiple nearby nodes on the graph
  std::vector<DenseVertex> graphNeighborhood;

  // Get the next state
  base::State *currentPathState = inputPath.getState(currVertexIndex);

  // Find multiple nearby nodes on the graph
  std::size_t threadID = 0;
  stateProperty_[queryVertices_[threadID]] = currentPathState;
  std::size_t findNearestKNeighbors = 10;
  const std::size_t numSameVerticiesFound = 1;  // add 1 to the end because the NN tree always returns itself
  nn_->nearestK(queryVertices_[threadID], findNearestKNeighbors + numSameVerticiesFound, graphNeighborhood);
  stateProperty_[queryVertices_[threadID]] = nullptr;

  // Loop through each neighbor until one is found that connects to the previous vertex
  bool foundValidConnToPrevious = false;
  // track if we added a vertex to the roadmapPath, so that we can remove it later if needed
  bool addedToRoadmapPath = false;
  DenseVertex candidateVertex;
  for (std::size_t neighborID = 0; neighborID < graphNeighborhood.size(); ++neighborID)
  {
    bool isValid = false;
    bool isRepeatOfPrevWaypoint = false;  // don't add current waypoint if same as last one

    // Find next state's vertex
    candidateVertex = graphNeighborhood[neighborID];

    // Check if next vertex is same as previous
    if (prevGraphVertex == candidateVertex)
    {
      // Do not do anything, we are done here
      foundValidConnToPrevious = true;
      if (verbose)
        std::cout << std::string(currVertexIndex, ' ') << "Previous vertex is same as current vertex, skipping "
                                                          "current vertex" << std::endl;

      isValid = true;
      isRepeatOfPrevWaypoint = true;
    }
    else
    {
      // Check for collision
      isValid = si_->checkMotion(stateProperty_[prevGraphVertex], stateProperty_[candidateVertex]);

      if (visualizeSnapPath_)  // Visualize
      {
        // Show the node we're currently considering going through
        visual_->viz5State(stateProperty_[candidateVertex], tools::MEDIUM, tools::PURPLE, 1);
        // edge between the state on the original inputPath and its neighbor we are currently considering
        double color = 25;  // light green
        visual_->viz5Edge(currentPathState, stateProperty_[candidateVertex], color);

        color = isValid ? 75 : 100;  // orange, red
        // edge between the previous connection point we chose for the roadmapPath, and the currently considered
        // next state
        visual_->viz5Edge(stateProperty_[prevGraphVertex], stateProperty_[candidateVertex], color);

        visual_->viz5Trigger();
        usleep(visualizeSnapPathSpeed_ * 1000000);
      }

      if (isValid && verbose)  // Debug
        std::cout << std::string(currVertexIndex, ' ') << "Found valid nearby edge on loop " << neighborID << std::endl;
    }

    // Remember if any connections failed
    if (isValid)
    {
      if (neighborID > 0)
      {
        if (verbose)
          std::cout << std::string(currVertexIndex + 2, ' ') << "Found case where double loop fixed the "
                                                                "problem - loop " << neighborID << std::endl;
        // visual_->viz5Trigger();
        // usleep(6*1000000);
      }
      foundValidConnToPrevious = true;

      // Add this waypoint solution
      if (!isRepeatOfPrevWaypoint)
      {
        // std::cout << std::string(currVertexIndex+2, ' ') << "roadmapPath.size=" << std::fixed <<
        // roadmapPath.size() << std::flush;
        // std::cout << " Vertex: " << candidateVertex;
        // std::cout << " State: " << stateProperty_[candidateVertex];
        // std::cout << std::endl;
        roadmapPath.push_back(candidateVertex);
        addedToRoadmapPath = true;  // remember it was added

        if (visualizeSnapPath_)  // Visualize
        {
          double color = 25;  // light green
          visual_->viz4Edge(stateProperty_[prevGraphVertex], stateProperty_[candidateVertex], color);
          visual_->viz4Trigger();
          usleep(visualizeSnapPathSpeed_ * 1000000);
        }
      }

      // Check if there are more points to process
      if (currVertexIndex + 1 >= inputPath.getStateCount())
      {
        if (verbose)
          std::cout << std::string(currVertexIndex, ' ') << "END OF PATH, great job :)" << std::endl;
        allValid = true;
        return true;
      }
      else
      {
        // Recurisvely call next level
        if (recurseSnapWaypoints(inputPath, roadmapPath, currVertexIndex + 1, candidateVertex, allValid, verbose))
        {
          return true;
        }
        else
        {
          // Keep trying to find a working neighbor, remove last roadmapPath node if we added one
          if (addedToRoadmapPath)
          {
            assert(roadmapPath.size() > 0);
            roadmapPath.pop_back();
            addedToRoadmapPath = false;
          }
        }
      }
    }
    else
    {
      if (verbose)
        std::cout << std::string(currVertexIndex, ' ') << "Loop " << neighborID << " not valid" << std::endl;
    }
  }  // for every neighbor

  if (!foundValidConnToPrevious)
  {
    std::cout << std::string(currVertexIndex, ' ') << "Unable to find valid connection to previous, backing up a "
                                                      "level and trying again" << std::endl;
    allValid = false;

    if (visualizeSnapPath_)  // Visualize
    {
      visual_->viz5Trigger();
      usleep(visualizeSnapPathSpeed_ * 1000000);
    }

    // TODO(davetcoleman): remove hack
    std::cout << std::string(currVertexIndex, ' ') << "TODO remove this viz" << std::endl;

    // Show the node we're currently considering going through
    visual_->viz5State(stateProperty_[prevGraphVertex], tools::VARIABLE_SIZE, tools::PURPLE, 3);
    visual_->viz5Trigger();
    usleep(0.001 * 1000000);

    return false;
  }
  std::cout << std::string(currVertexIndex, ' ') << "This loop found a valid connection, but higher recursive loop "
                                                    "(one that has already returned) did not" << std::endl;

  return false;  // this loop found a valid connection, but lower recursive loop did not
}

bool DenseDB::astarSearch(const DenseVertex start, const DenseVertex goal, std::vector<DenseVertex> &vertexPath)
{
  // Hold a list of the shortest path parent to each vertex
  DenseVertex *vertexPredecessors = new DenseVertex[getNumVertices()];
  // boost::vector_property_map<DenseVertex> vertexPredecessors(getNumVertices());

  bool foundGoal = false;
  double *vertexDistances = new double[getNumVertices()];

  // Error check
  if (useTaskPlanning_)
  {
    if (getTaskLevel(start) != 0)
    {
      OMPL_ERROR("astarSearch: start level is %u", getTaskLevel(start));
      exit(-1);
    }
    if (getTaskLevel(goal) != 2)
    {
      OMPL_ERROR("astarSearch: goal level is %u", getTaskLevel(goal));
      exit(-1);
    }
  }

  OMPL_INFORM("Beginning AStar Search");
  try
  {
    // Note: could not get astar_search to compile within BoltRetrieveRepair.cpp class because of namespacing issues
    boost::astar_search(
        g_,     // graph
        start,  // start state
                // boost::bind(&DenseDB::distanceFunction2, this, _1, goal),  // the heuristic
        boost::bind(&DenseDB::distanceFunction, this, _1, goal),  // the heuristic
        // boost::bind(&DenseDB::distanceFunctionTasks, this, _1, goal),  // the heuristic
        // ability to disable edges (set cost to inifinity):
        boost::weight_map(DenseEdgeWeightMap(g_, edgeCollisionStateProperty_, popularityBias_, popularityBiasEnabled_))
            .predecessor_map(vertexPredecessors)
            .distance_map(&vertexDistances[0])
            .visitor(CustomAstarVisitor(goal, this)));
  }
  catch (FoundGoalException &)
  {
    // the custom exception from CustomAstarVisitor
    OMPL_INFORM("astarSearch: Astar found goal vertex. distance to goal: %f", vertexDistances[goal]);

    if (vertexDistances[goal] > 1.7e+308)  // TODO(davetcoleman): fix terrible hack for detecting infinity
                                           // double diff = d[goal] - std::numeric_limits<double>::infinity();
    // if ((diff < std::numeric_limits<double>::epsilon()) && (-diff < std::numeric_limits<double>::epsilon()))
    // check if the distance to goal is inifinity. if so, it is unreachable
    // if (d[goal] >= std::numeric_limits<double>::infinity())
    {
      if (verbose_)
        OMPL_INFORM("Distance to goal is infinity");
      foundGoal = false;
    }
    else
    {
      // Only clear the vertexPath after we know we have a new solution, otherwise it might have a good
      // previous one
      vertexPath.clear();  // remove any old solutions

      // Trace back the shortest path in reverse and only save the states
      DenseVertex v;
      for (v = goal; v != vertexPredecessors[v]; v = vertexPredecessors[v])
      {
        vertexPath.push_back(v);
      }
      if (v != goal)  // TODO explain this because i don't understand
      {
        vertexPath.push_back(v);
      }

      foundGoal = true;
    }
  }

  if (!foundGoal)
    OMPL_WARN("        Did not find goal");

  // Show all predecessors
  if (visualizeAstar_)
  {
    OMPL_INFORM("        Show all predecessors");
    for (std::size_t i = 1; i < getNumVertices(); ++i)  // skip vertex 0 b/c that is the search vertex
    {
      const DenseVertex v1 = i;
      const DenseVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4Edge(stateProperty_[v1], stateProperty_[v2], 10);
      }
    }
    visual_->viz4Trigger();
  }

  // Unload
  delete[] vertexPredecessors;
  delete[] vertexDistances;

  // No solution found from start to goal
  return foundGoal;
}

void DenseDB::computeDensePath(const DenseVertex start, const DenseVertex goal, DensePath &path)
{
  path.clear();

  boost::vector_property_map<DenseVertex> prev(boost::num_vertices(g_));
  /*
  try
  {
    boost::astar_search(g_,                                                            // graph
                        start,                                                         // start state
                        boost::bind(&DenseDB::distanceFunctionTasks, this, _1, goal),  // the heuristic
                        boost::predecessor_map(prev).visitor(CustomAstarVisitor(goal, this)));
  }
  catch (FoundGoalException &)
  {
  }

  if (prev[goal] == goal)
    OMPL_WARN("No dense path was found?");
  else
  {
    for (DenseVertex pos = goal; prev[pos] != pos; pos = prev[pos])
      path.push_front(stateProperty_[pos]);
    path.push_front(stateProperty_[start]);
  }
  */
}

void DenseDB::debugVertex(const ompl::base::PlannerDataVertex &vertex)
{
  debugState(vertex.getState());
}

void DenseDB::debugState(const ompl::base::State *state)
{
  si_->printState(state, std::cout);
}

double DenseDB::distanceFunction(const DenseVertex a, const DenseVertex b) const
{
  // const double dist = si_->distance(stateProperty_[a], stateProperty_[b]);
  // std::cout << "getting distance from " << a << " to " << b << " of value " << dist << std::endl;
  // return dist;
  return si_->distance(stateProperty_[a], stateProperty_[b]);
}

double DenseDB::distanceFunction2(const DenseVertex a, const DenseVertex b) const
{
  // const double dist = si_->getStateSpace()->distance2(stateProperty_[a], stateProperty_[b]);
  // std::cout << "getting distance from " << a << " to " << b << " of value " << dist << std::endl;
  // return dist;
  return si_->getStateSpace()->distance2(stateProperty_[a], stateProperty_[b]);
}

std::size_t DenseDB::getTaskLevel(const DenseVertex &v) const
{
  return si_->getStateSpace()->getLevel(stateProperty_[v]);
}

std::size_t DenseDB::getTaskLevel(const base::State *state) const
{
  return si_->getStateSpace()->getLevel(state);
}

void DenseDB::initializeQueryState()
{
  std::size_t numThreads = boost::thread::hardware_concurrency();

  if (boost::num_vertices(g_) > 0)
  {
    assert(boost::num_vertices(g_) >= numThreads);
    return;  // assume its already been setup
  }

  // Create a query state for each possible thread
  queryVertices_.resize(numThreads);

  for (std::size_t threadID = 0; threadID < numThreads; ++threadID)
  {
    // Add a fake vertex to the graph
    queryVertices_[threadID] = boost::add_vertex(g_);

    // Set its state to nullptr
    stateProperty_[queryVertices_[threadID]] = nullptr;
  }
}

void DenseDB::addVertexFromFile(BoltStorage::BoltVertexData v)
{
  GuardType type = static_cast<GuardType>(v.type_);
  // DenseVertex vNew =
  addVertex(v.state_, type);
}

void DenseDB::addEdgeFromFile(BoltStorage::BoltEdgeData e)
{
  const DenseVertex v1 = e.endpoints_.first;
  const DenseVertex v2 = e.endpoints_.second;

  // Error check
  BOOST_ASSERT_MSG(v1 <= getNumVertices(), "Vertex 1 out of range of possible verticies");
  BOOST_ASSERT_MSG(v2 <= getNumVertices(), "Vertex 2 out of range of possible verticies");

  // Add edge
  // DenseEdge newE =
  addEdge(v1, v2, e.weight_);
}

void DenseDB::clearEdgeCollisionStates()
{
  foreach (const DenseEdge e, boost::edges(g_))
    edgeCollisionStateProperty_[e] = NOT_CHECKED;  // each edge has an unknown state
}

void DenseDB::displayDatabase()
{
  OMPL_INFORM("Displaying database");

  // Clear old database
  visual_->viz1DeleteAllMarkers();

  if (visualizeDatabaseVertices_)
  {
    // Error check
    if (getNumVertices() == 0)
    {
      OMPL_WARN("Unable to show complete database because no vertices available");
    }

    // Loop through each vertex
    std::size_t count = 0;
    std::size_t debugFrequency = std::min(10000, static_cast<int>(getNumVertices() / 10));
    std::cout << "Displaying vertices: " << std::flush;
    foreach (DenseVertex v, boost::vertices(g_))
    {
      // Check for null states
      if (stateProperty_[v])
      {
        visual_->viz1State(stateProperty_[v], tools::SMALL, tools::BLUE, 1);
      }

      // Prevent cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumVertices()) * 100.0
                  << "% " << std::flush;
        visual_->viz1Trigger();
      }
      count++;
    }
    std::cout << std::endl;
  }

  if (visualizeDatabaseEdges_)
  {
    // Error check
    if (getNumEdges() == 0)
    {
      OMPL_WARN("Unable to show complete database because no edges available");
    }
    // Loop through each edge
    std::size_t count = 0;
    std::size_t debugFrequency = std::min(10000, static_cast<int>(getNumEdges() / 10));
    std::cout << "Displaying edges: " << std::flush;
    foreach (DenseEdge e, boost::edges(g_))
    {
      // Add edge
      const DenseVertex &v1 = boost::source(e, g_);
      const DenseVertex &v2 = boost::target(e, g_);

      // Visualize
      assert(edgeWeightProperty_[e] <= MAX_POPULARITY_WEIGHT);
      visual_->viz1Edge(stateProperty_[v1], stateProperty_[v2], edgeWeightProperty_[e]);

      // Prevent cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumEdges()) * 100.0
                  << "% " << std::flush;
        visual_->viz1Trigger();
      }

      count++;
    }
    std::cout << std::endl;
  }

  // Publish remaining markers
  visual_->viz1Trigger();
}

void DenseDB::normalizeGraphEdgeWeights()
{
  bool verbose = false;

  if (!popularityBiasEnabled_)
  {
    OMPL_INFORM("Skipping normalize graph edge weights because not using popularity bias currently");
    return;
  }

  // Normalize weight of graph
  double total_cost = 0;
  foreach (DenseEdge e, boost::edges(g_))
  {
    total_cost += edgeWeightProperty_[e];
  }
  double avg_cost = total_cost / getNumEdges();

  if (verbose)
    OMPL_INFORM("Average cost of the edges in graph is %f", avg_cost);

  if (avg_cost < desiredAverageCost_)  // need to decrease cost in graph
  {
    double avgCostDiff = desiredAverageCost_ - avg_cost;

    if (verbose)
      std::cout << "avgCostDiff: " << avgCostDiff << std::endl;
    double perEdgeReduction = avgCostDiff;  // / getNumEdges();

    if (verbose)
      OMPL_INFORM("Decreasing each edge's cost by %f", perEdgeReduction);
    foreach (DenseEdge e, boost::edges(g_))
    {
      assert(edgeWeightProperty_[e] <= MAX_POPULARITY_WEIGHT);
      edgeWeightProperty_[e] = std::min(edgeWeightProperty_[e] + perEdgeReduction, MAX_POPULARITY_WEIGHT);
      if (verbose)
        std::cout << "Edge " << e << " now has weight " << edgeWeightProperty_[e] << " via reduction "
                  << perEdgeReduction << std::endl;
    }
  }
  else if (verbose)
  {
    OMPL_INFORM("Not decreasing all edge's cost because average is above desired");
  }
}

otb::DenseVertex DenseDB::addVertex(base::State *state, const GuardType &type)
{
  // Create vertex
  DenseVertex v = boost::add_vertex(g_);

  // Add properties
  //typeProperty_[v] = type;
  stateProperty_[v] = state;
  //representativesProperty_[v] = 0;  // which sparse vertex reps this dense vertex

  // Connected component tracking
  //disjointSets_.make_set(v);

  // Add vertex to nearest neighbor structure
  //nn_->add(v);

  // Track vertex for later removal if temporary
  // if (type == CARTESIAN)
  // {
  //   tempVerticies_.push_back(v);
  // }

  return v;
}

otb::DenseEdge DenseDB::addEdge(const DenseVertex &v1, const DenseVertex &v2, const double weight,
                                const EdgeCollisionState collisionState)
{
  // Error check
  BOOST_ASSERT_MSG(v1 <= getNumVertices(), "Vertex 1 out of range of possible verticies");
  BOOST_ASSERT_MSG(v2 <= getNumVertices(), "Vertex 2 out of range of possible verticies");

  // Create the new edge
  DenseEdge e = (boost::add_edge(v1, v2, g_)).first;

  // std::cout << "Adding cost: " << weight << std::endl;

  // Add associated properties to the edge
  edgeWeightProperty_[e] = weight;
  // edgeWeightProperty_[e] = distanceFunction2(v1, v2);
  // edgeWeightProperty_[e] = distanceFunction(v1, v2);
  edgeCollisionStateProperty_[e] = collisionState;

  // Add the edge to the incrementeal connected components datastructure
  disjointSets_.union_set(v1, v2);

  return e;
}

void DenseDB::cleanupTemporaryVerticies()
{
  const bool verbose = false;

  if (tempVerticies_.empty())
  {
    OMPL_INFORM("Skipping verticies cleanup - no middle cartesian layer verticies found");
    return;
  }

  OMPL_INFORM("Cleaning up temp verticies - vertex count: %u, edge count: %u", getNumVertices(), getNumEdges());
  BOOST_REVERSE_FOREACH(DenseVertex v, tempVerticies_)
  {
    removeVertex(v);

    if (verbose)
      OMPL_DEBUG("Removed, updated - vertex count: %u, edge count: %u", getNumVertices(), getNumEdges());
  }
  tempVerticies_.clear();
  OMPL_INFORM("Finished cleaning up temp verticies");
}

void DenseDB::removeVertex(DenseVertex v)
{
  const bool verbose = false;

  if (verbose)
    std::cout << "Removing vertex " << v << std::endl;

  // Remove from nearest neighbor
  nn_->remove(v);

  // Delete state
  si_->freeState(stateProperty_[v]);
  stateProperty_[v] = nullptr;

  // Remove all edges to and from vertex
  boost::clear_vertex(v, g_);

  // Remove vertex
  boost::remove_vertex(v, g_);
}

void DenseDB::findGraphNeighbors(base::State *state, std::vector<DenseVertex> &graphNeighborhood,
                                 std::vector<DenseVertex> &visibleNeighborhood, double searchRadius,
                                 std::size_t threadID, std::size_t coutIndent)
{
  // Set a queryVertex to give us a DenseVertex
  stateProperty_[queryVertices_[threadID]] = state;

  // Search
  findGraphNeighbors(queryVertices_[threadID], graphNeighborhood, visibleNeighborhood, searchRadius, coutIndent);

  // Reset a queryVertex
  stateProperty_[queryVertices_[threadID]] = nullptr;
}

void DenseDB::findGraphNeighbors(const DenseVertex &denseV, std::vector<DenseVertex> &graphNeighborhood,
                                 std::vector<DenseVertex> &visibleNeighborhood, double searchRadius,
                                 std::size_t coutIndent)
{
  bool verbose = false;
  if (verbose)
    std::cout << std::string(coutIndent, ' ') << "findGraphNeighbors()" << std::endl;

  nn_->nearestR(denseV, searchRadius, graphNeighborhood);

  // Now that we got the neighbors from the NN, remove any we can't see
  for (DenseVertex &denseV : graphNeighborhood)
  // for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  {
    if (si_->checkMotion(stateProperty_[denseV], stateProperty_[denseV]))
    {
      visibleNeighborhood.push_back(denseV);
    }
  }

  if (verbose)
    std::cout << std::string(coutIndent + 2, ' ') << "Graph neighborhood: " << graphNeighborhood.size()
              << " | Visible neighborhood: " << visibleNeighborhood.size() << std::endl;
}

void DenseDB::viz1Edge(DenseEdge &e)
{
  const DenseVertex &v1 = boost::source(e, g_);
  const DenseVertex &v2 = boost::target(e, g_);

  // Visualize
  visual_->viz1Edge(stateProperty_[v1], stateProperty_[v2], edgeWeightProperty_[e]);
}

void DenseDB::checkStateType()
{
  std::size_t count = 0;
  foreach (const DenseVertex v, boost::vertices(g_))
  {
    // The first vertex (id=0) should have a nullptr state because it is used for searching
    if (!stateProperty_[v])
    {
      if (count != 0)
      {
        OMPL_ERROR("Null state found for vertex that is not zero");
      }
      continue;
    }

    // std::size_t level = getTaskLevel(stateProperty_[v]);
    // if (level > 2)
    // {
    //   OMPL_ERROR("State is on wrong level: %u", level);
    //   exit(-1);
    // }
  }
  OMPL_INFORM("All states checked for task level");
}

void DenseDB::connectNewVertex(DenseVertex v1)
{
  bool verbose = false;

  // Visualize new vertex
  if (visualizeAddSample_)
  {
    visual_->viz1State(stateProperty_[v1], tools::SMALL, tools::GREEN, 1);
  }

  // Connect to nearby vertices
  std::vector<DenseVertex> graphNeighborhood;
  std::size_t findNearestKNeighbors = Discretizer::getEdgesPerVertex(si_);
  OMPL_INFORM("DenseDB.connectNewVertex(): Finding %u nearest neighbors for new vertex", findNearestKNeighbors);
  const std::size_t numSameVerticiesFound = 1;  // add 1 to the end because the NN tree always returns itself

  // Search
  const std::size_t threadID = 0;
  stateProperty_[queryVertices_[threadID]] = stateProperty_[v1];
  nn_->nearestK(queryVertices_[threadID], findNearestKNeighbors + numSameVerticiesFound, graphNeighborhood);
  stateProperty_[queryVertices_[threadID]] =
      nullptr;  // Set search vertex to nullptr to prevent segfault on class unload of memory

  // For each nearby vertex, add an edge
  std::size_t errorCheckNumSameVerticies = 0;  // sanity check
  for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  {
    DenseVertex &v2 = graphNeighborhood[i];

    // Check if these vertices are the same
    if (v1 == v2)
    {
      if (verbose)
        std::cout << "connectNewVertex attempted to connect edge to itself: " << v1 << ", " << v2 << std::endl;
      errorCheckNumSameVerticies++;  // sanity check
      continue;
    }

    // Check if these vertices are the same STATE
    if (si_->getStateSpace()->equalStates(stateProperty_[v1], stateProperty_[v2]))
    {
      OMPL_ERROR("This state has already been added, should not happen");
      exit(-1);
    }

    // Check if these vertices already share an edge
    if (boost::edge(v1, v2, g_).second)
    {
      std::cout << "skipped bc already share an edge " << std::endl;
      continue;
    }

    // Create edge if not in collision
    if (si_->checkMotion(stateProperty_[v1], stateProperty_[v2]))
    {
      // DenseEdge e =
      addEdge(v1, v2, desiredAverageCost_);
      // std::cout << "added valid edge " << e << std::endl;

      if (visualizeAddSample_)  // Debug: display edge
      {
        double popularity = 100;  // TODO: maybe make edge really popular so we can be sure its added to the
                                  // spars graph since we need it
        visual_->viz1Edge(stateProperty_[v1], stateProperty_[v2], popularity);
      }
    }

  }  // for each neighbor

  // Make sure one and only one vertex is returned from the NN search that is the same as parent vertex
  BOOST_ASSERT_MSG(errorCheckNumSameVerticies <= numSameVerticiesFound, "Too many same verticies found ");

  // Visualize
  if (visualizeAddSample_)
  {
    visual_->viz1Trigger();
    usleep(0.001 * 1000000);
  }

  // Record this new addition
  graphUnsaved_ = true;
}

std::size_t DenseDB::getDisjointSetsCount(bool verbose)
{
  std::size_t numSets = 0;
  foreach (DenseVertex v, boost::vertices(g_))
  {
    // Do not count the search vertex within the sets
    if (v <= queryVertices_.back())
      continue;

    if (boost::get(boost::get(boost::vertex_predecessor, g_), v) == v)
    {
      if (verbose)
        OMPL_INFORM("Disjoint set: %u", v);
      ++numSets;
    }
  }

  return numSets;
}

std::size_t DenseDB::checkConnectedComponents()
{
  // Check how many disjoint sets are in the dense graph (should be none)
  std::size_t numSets = getDisjointSetsCount();
  if (numSets > 1)
  {
    OMPL_ERROR("More than 1 connected component is in the dense graph: %u", numSets);
  }

  return numSets;
}

bool DenseDB::sameComponent(const DenseVertex &v1, const DenseVertex &v2)
{
  return boost::same_component(v1, v2, disjointSets_);
}

void DenseDB::removeInvalidVertices()
{
  OMPL_INFORM("Removing invalid vertices");
  bool actuallyRemove = true;
  std::size_t totalInvalid = 0;  // invalid states

  if (actuallyRemove)
    OMPL_WARN("Actually deleting verticies and resetting edge cache");

  typedef boost::graph_traits<DenseGraph>::vertex_iterator VertexIterator;
  for (VertexIterator vertexIt = boost::vertices(g_).first; vertexIt != boost::vertices(g_).second; ++vertexIt)
  {
    if (*vertexIt <= queryVertices_.back())
      continue;

    if (*vertexIt % 1000 == 0)
      std::cout << "Checking vertex " << *vertexIt << std::endl;

    // Check if state is valid
    if (!si_->isValid(stateProperty_[*vertexIt]))
    {
      OMPL_ERROR("State %u is not valid", *vertexIt);
      totalInvalid++;

      visual_->viz5State(stateProperty_[*vertexIt], tools::ROBOT, tools::RED, 0);
      visual_->viz5Trigger();
      // usleep(0.25 * 1000000);

      if (actuallyRemove)
      {
        removeVertex(*vertexIt);
        vertexIt--;  // backtrack one vertex
      }
    }
  }

  if (totalInvalid > 0)
  {
    OMPL_ERROR("Total invalid: %u", totalInvalid);

    if (actuallyRemove)
    {
      // Must clear out edge cache since we changed the numbering of vertices
      denseCache_->clear();
      graphUnsaved_ = true;
    }
  }
}

void DenseDB::getDisjointSets(DisjointSetsParentKey &disjointSets)
{
  OMPL_INFORM("Get disjoint sets...");
  disjointSets.clear();

  // Flatten the parents tree so that the parent of every element is its representative.
  disjointSets_.compress_sets(boost::vertices(g_).first, boost::vertices(g_).second);

  // Count size of each disjoint set and group its containing vertices
  typedef boost::graph_traits<DenseGraph>::vertex_iterator VertexIterator;
  for (VertexIterator v = boost::vertices(g_).first; v != boost::vertices(g_).second; ++v)
  {
    // Do not count the search vertex within the sets
    if (*v <= queryVertices_.back())
      continue;

    disjointSets[boost::get(boost::get(boost::vertex_predecessor, g_), *v)].push_back(*v);
  }
}

void DenseDB::printDisjointSets(DisjointSetsParentKey &disjointSets)
{
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const DenseVertex v = iterator->first;
    const std::size_t freq = iterator->second.size();
    std::cout << "Parent: " << v << " frequency " << freq << std::endl;
  }
}

void DenseDB::visualizeDisjointSets(DisjointSetsParentKey &disjointSets)
{
  OMPL_INFORM("Visualizing disjoint sets");

  // Find the disjoint set that is the 'main' large one
  std::size_t maxDisjointSetSize = 0;
  DenseVertex maxDisjointSetParent;
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const DenseVertex v = iterator->first;
    const std::size_t freq = iterator->second.size();

    if (freq > maxDisjointSetSize)
    {
      maxDisjointSetSize = freq;
      maxDisjointSetParent = v;
    }
  }
  OMPL_INFORM("The largest disjoint set is of size %u with parent %u", maxDisjointSetSize, maxDisjointSetParent);

  // Display size of disjoint sets and visualize small ones
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const DenseVertex v1 = iterator->first;
    const std::size_t freq = iterator->second.size();
    std::cout << v1 << ": frequency: " << freq << std::endl;

    BOOST_ASSERT_MSG(freq > 0, "Frequnecy must be at least 1");

    if (freq == maxDisjointSetSize)  // any subgraph that is smaller than the full graph
      continue;                      // the main disjoint set is not considered a disjoint set

    // Visualize sets of size one
    if (freq == 1 && true)
    {
      visual_->viz5State(stateProperty_[v1], tools::ROBOT, tools::RED, 0);
      visual_->viz5Trigger();
      usleep(1.0 * 1000000);
    }

    // Visualize large disjoint sets (greater than one)
    if (freq > 1 && freq < 1000)
    {
      // Clear markers
      visual_->viz4DeleteAllMarkers();

      // Visualize this subgraph that is disconnected
      // Loop through every every vertex and check if its part of this group
      typedef boost::graph_traits<DenseGraph>::vertex_iterator VertexIterator;
      for (VertexIterator v2 = boost::vertices(g_).first; v2 != boost::vertices(g_).second; ++v2)
      {
        if (boost::get(boost::get(boost::vertex_predecessor, g_), *v2) == v1)
        {
          visual_->viz4State(stateProperty_[*v2], tools::LARGE, tools::RED, 0);

          // Show state's edges
          foreach (DenseEdge edge, boost::out_edges(*v2, g_))
          {
            DenseVertex e_v1 = boost::source(edge, g_);
            DenseVertex e_v2 = boost::target(edge, g_);
            visual_->viz4Edge(stateProperty_[e_v1], stateProperty_[e_v2], edgeWeightProperty_[edge]);
          }
          visual_->viz4Trigger();

          // Show this robot state
          visual_->viz4State(stateProperty_[*v2], tools::ROBOT, tools::DEFAULT, 0);

          usleep(0.1 * 1000000);
        }
      }
      if (true)
      {
        usleep(2 * 1000000);
      }
    }
  }
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
