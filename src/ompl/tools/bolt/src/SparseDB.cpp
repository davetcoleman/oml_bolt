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
#include <ompl/tools/bolt/SparseDB.h>
#include <ompl/util/Time.h>
#include <ompl/util/Console.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>  // TODO: remove, this is not space agnostic

// Boost
#include <boost/graph/incremental_components.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/assert.hpp>

// TODO remove
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>

// C++
#include <limits>
#include <queue>
#include <algorithm>  // std::random_shuffle

// Profiling
#include <valgrind/callgrind.h>

#define foreach BOOST_FOREACH

namespace og = ompl::geometric;
namespace ot = ompl::tools;
namespace otb = ompl::tools::bolt;
namespace ob = ompl::base;

// Actual class ////////////////////////////////////////////////////////////////////////////

namespace ompl
{
namespace tools
{
namespace bolt
{
SparseDB::SparseDB(base::SpaceInformationPtr si, DenseDB *denseDB, base::VisualizerPtr visual, EdgeCachePtr edgeCache)
  : si_(si)
  , denseDB_(denseDB)
  , visual_(visual)
  , edgeCache_(edgeCache)
  , smoothingGeomPath_(si)
    // Property accessors of edges
  , edgeWeightPropertySparse_(boost::get(boost::edge_weight, g_))
  , edgeCollisionStatePropertySparse_(boost::get(edge_collision_state_t(), g_))
    // Property accessors of vertices
  , denseVertexProperty_(boost::get(vertex_dense_pointer_t(), g_))
  , typePropertySparse_(boost::get(vertex_type_t(), g_))
    //, nonInterfaceListsProperty_(boost::get(vertex_list_t(), g_))
    //, interfaceListsProperty_(boost::get(vertex_interface_list_t(), g_))
  , vertexPopularity_(boost::get(vertex_popularity_t(), g_))
  // Disjoint set accessors
  , disjointSets_(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_))
{
  // Add search state
  initializeQueryState();

  // Initialize nearest neighbor datastructure
  nn_.reset(new NearestNeighborsGNAT<SparseVertex>());
  nn_->setDistanceFunction(boost::bind(&otb::SparseDB::distanceFunction, this, _1, _2));

  // Initialize path simplifier
  psimp_.reset(new geometric::PathSimplifier(si_));
  psimp_->freeStates(false);
}

SparseDB::~SparseDB(void)
{
  freeMemory();
}

void SparseDB::freeMemory()
{
  g_.clear();

  if (nn_)
    nn_->clear();
}

bool SparseDB::setup()
{
  // Calculate variables for the graph
  maxExtent_ = si_->getMaximumExtent();
  // sparseDelta_ = sparseDeltaFraction_ * maxExtent_;
  sparseDelta_ = sparseDeltaFraction_ * discretization_;  // sparseDelta should be a multiple of discretization
  denseDelta_ = denseDeltaFraction_ * maxExtent_;
  OMPL_INFORM("sparseDelta_ = %f", sparseDelta_);
  // OMPL_INFORM("denseDelta_ = %f", denseDelta_);

  assert(maxExtent_ > 0);
  assert(denseDelta_ > 0);
  assert(sparseDelta_ > 0);
  assert(sparseDelta_ > 0.000000001);  // Sanity check

  return true;
}

/*
  bool SparseDB::astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex> &vertexPath)
  {
  // Hold a list of the shortest path parent to each vertex
  SparseVertex *vertexPredecessors = new SparseVertex[getNumVertices()];
  // boost::vector_property_map<SparseVertex> vertexPredecessors(getNumVertices());

  bool foundGoal = false;
  double *vertexDistances = new double[getNumVertices()];

  // Reset statistics
  numNodesOpened_ = 0;
  numNodesClosed_ = 0;

  OMPL_INFORM("Beginning AStar Search");
  try
  {
  double popularityBias = 0;
  bool popularityBiasEnabled = false;
  // Note: could not get astar_search to compile within BoltRetrieveRepair.cpp class because of
  // namespacing issues
  boost::astar_search(g_,                                                           // graph
  start,                                                        // start state
  boost::bind(&otb::SparseDB::astarHeuristic, this, _1, goal),  // the heuristic
  // ability to disable edges (set cost to inifinity):
  boost::weight_map(SparseEdgeWeightMap(g_, edgeCollisionStatePropertySparse_, popularityBias,
  popularityBiasEnabled))
  .predecessor_map(vertexPredecessors)
  .distance_map(&vertexDistances[0])
  .visitor(CustomAstarVisitor(goal, sparseDB_)));
  }
  catch (FoundGoalException &)
  {
  // the custom exception from CustomAstarVisitor
  OMPL_INFORM("AStar found goal vertex. distance to goal: %f", vertexDistances[goal]);
  OMPL_INFORM("Number nodes opened: %u, Number nodes closed: %u", numNodesOpened_, numNodesClosed_);

  if (vertexDistances[goal] > 1.7e+308)  // TODO(davetcoleman): fix terrible hack for detecting infinity
  // double diff = d[goal] - std::numeric_limits<double>::infinity();
  // if ((diff < std::numeric_limits<double>::epsilon()) && (-diff <
  // std::numeric_limits<double>::epsilon()))
  // check if the distance to goal is inifinity. if so, it is unreachable
  // if (d[goal] >= std::numeric_limits<double>::infinity())
  {
  OMPL_INFORM("Distance to goal is infinity");
  foundGoal = false;
  }
  else
  {
  // Only clear the vertexPath after we know we have a new solution, otherwise it might have a good
  // previous one
  vertexPath.clear();  // remove any old solutions

  // Trace back the shortest path in reverse and only save the states
  SparseVertex v;
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
  const SparseVertex v1 = i;
  const SparseVertex v2 = vertexPredecessors[v1];
  if (v1 != v2)
  {
  // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
  visual_->viz4Edge(getSparseStateConst(v1), getSparseStateConst(v2), 10);
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
*/

void SparseDB::debugVertex(const ompl::base::PlannerDataVertex &vertex)
{
  debugState(vertex.getState());
}

void SparseDB::debugState(const ompl::base::State *state)
{
  si_->printState(state, std::cout);
}

void SparseDB::initializeQueryState()
{
  if (boost::num_vertices(g_) > 0)
  {
    OMPL_WARN("Not initializing query state because already is of size %u", boost::num_vertices(g_));
    return; // assume its already been setup
  }

  // Make sure dense graph has been initialized
  denseDB_->initializeQueryState();

  // Create a query state for each possible thread
  std::size_t numThreads = boost::thread::hardware_concurrency();
  queryVertices_.resize(numThreads);

  for (std::size_t threadID = 0; threadID < numThreads; ++threadID)
  {
    // Add a fake vertex to the graph
    queryVertices_[threadID] = boost::add_vertex(g_);

    // Set this fake vertex's dense vertex to correspond to the fake dense vertex
    denseVertexProperty_[queryVertices_[threadID]] = denseDB_->queryVertices_[threadID];
  }
}

void SparseDB::clearEdgeCollisionStates()
{
  foreach (const SparseEdge e, boost::edges(g_))
    edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;  // each edge has an unknown state
}

void SparseDB::createSPARS()
{
  // Benchmark runtime
  time::point startTime = time::now();

  // Create SPARS
  CALLGRIND_TOGGLE_COLLECT;
  createSPARSOuterLoop();
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_DUMP_STATS;

  // Benchmark runtime
  double duration = time::seconds(time::now() - startTime);

  // Statistics
  numGraphGenerations_++;

  // Check how many connected components exist, possibly throw error
  std::size_t numSets = getDisjointSetsCount();

  OMPL_INFORM("-----------------------------------------");
  OMPL_INFORM("Created SPARS graph                      ");
  OMPL_INFORM("  Vertices:                  %u", getNumVertices());
  OMPL_INFORM("  Edges:                     %u", getNumEdges());
  OMPL_INFORM("  Generation time:           %f", duration);
  OMPL_INFORM("  Total generations:         %u", numGraphGenerations_);
  OMPL_INFORM("  Disjoint sets:             %u", numSets);
  OMPL_INFORM("  Edge collision cache:      %u", edgeCache_->getCacheSize());
  OMPL_INFORM("  Total collision checks:    %u", edgeCache_->getTotalCollisionChecks());
  OMPL_INFORM("  Cached collision checks:   %u (%f %)", edgeCache_->getTotalCollisionChecksFromCache(),
              edgeCache_->getPercentCachedCollisionChecks());
  OMPL_INFORM("-----------------------------------------");

  if (numSets > 1 && false)
  {
    OMPL_INFORM("Disjoint sets: %u, attempting to random sample until fully connected", numSets);

    eliminateDisjointSets();
    denseDB_->displayDatabase();
  }
  else
    OMPL_WARN("Skipping SPARSE disjoint set fixing");

  // Save collision cache
  edgeCache_->save();

  OMPL_INFORM("Finished creating sparse database");
}

void SparseDB::createSPARSOuterLoop()
{
  std::size_t coutIndent = 2;

  // Clear the old spars graph
  if (getNumVertices() > queryVertices_.size())
  {
    OMPL_INFORM("Resetting sparse database");
    freeMemory();
    initializeQueryState();  // Re-add search state

    if (visualizeSparsGraph_)  // Clear visuals
      visual_->viz2State(NULL, /*deleteAllMarkers*/ 0, 0);
  }

  // Reset parameters
  setup();
  visualizeOverlayNodes_ = false;  // DO NOT visualize all added nodes in a separate window
  edgeCache_->resetCounters();

  // Get the ordering to insert vertices
  std::list<WeightedVertex> vertexInsertionOrder;
  getVertexInsertionOrdering(vertexInsertionOrder);

  // Error check order creation
  assert(vertexInsertionOrder.size() == denseDB_->getNumVertices() - queryVertices_.size());

  // Attempt to insert the verticies multiple times until no more succesful insertions occur
  secondSparseInsertionAttempt_ = false;
  std::size_t loopAttempt = 0;
  std::size_t sucessfulInsertions = 1;  // start with one so that while loop works
  while (sucessfulInsertions > 0)
  {
    // if (checksVerbose_)
    std::cout << "Attempting to insert " << vertexInsertionOrder.size()
              << " vertices for the " << loopAttempt << " loop" << std::endl;

    // Sanity check
    if (loopAttempt > 3)
      OMPL_WARN("Suprising number of loop when attempting to insert nodes into SPARS graph: %u", loopAttempt);

    // Benchmark runtime
    time::point startTime = time::now();

    // ----------------------------------------------------------------------
    // Attempt to insert each vertex using the first 3 criteria
    createSPARSInnerLoop(vertexInsertionOrder, sucessfulInsertions);

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);

    // Visualize
    if (visualizeSparsGraph_)
      visual_->viz2Trigger();

    std::cout << "Succeeded in inserting " << sucessfulInsertions << " vertices on the " << loopAttempt
              << " loop, remaining uninserted verticies: " << vertexInsertionOrder.size()
              << " loop runtime: " << duration << " sec" << std::endl;
    loopAttempt++;

    // Increase the sparse delta a bit, but only after the first loop
    if (loopAttempt == 1)
    {
      sparseDelta_ = sparseDelta_ * 1.25;
      std::cout << std::string(coutIndent + 2, ' ') << "sparseDelta_ is now " << sparseDelta_ << std::endl;
      secondSparseInsertionAttempt_ = true;

      // Save collision cache, just in case there is a bug
      edgeCache_->save();
    }

    bool debugOverRideJustTwice = false;
    if (debugOverRideJustTwice && loopAttempt == 2)
    {
      OMPL_WARN("Only attempting to add nodes twice for speed");
      break;
    }
  }

  // If we haven't animated the creation, just show it all at once
  if (!visualizeSparsGraph_)
  {
    displaySparseDatabase();
  }
  else if (visualizeSparsGraphSpeed_ < std::numeric_limits<double>::epsilon())
  {
    visual_->viz2Trigger();
    usleep(0.001 * 1000000);
  }
}

bool SparseDB::createSPARSInnerLoop(std::list<WeightedVertex> &vertexInsertionOrder, std::size_t &sucessfulInsertions)
{
  sucessfulInsertions = 0;
  std::size_t loopCount = 0;
  std::size_t originalVertexInsertion = vertexInsertionOrder.size();
  std::size_t debugFrequency = std::max(std::size_t(10), static_cast<std::size_t>(originalVertexInsertion / 20));
  std::cout << "originalVertexInsertion: " << originalVertexInsertion << std::endl;

  for (std::list<WeightedVertex>::iterator vertexIt = vertexInsertionOrder.begin();
       vertexIt != vertexInsertionOrder.end();
       /* will manually progress since there are erases*/)
  {
    // User feedback
    if (loopCount++ % debugFrequency == 0)
    {
      std::cout << std::fixed << std::setprecision(1)
                << "  Sparse generation: " << (static_cast<double>(loopCount) / originalVertexInsertion) * 100.0
                << "% Cache size: " << edgeCache_->getCacheSize()
                << " Cache usage: " << edgeCache_->getPercentCachedCollisionChecks() << "%" << std::endl;
      if (visualizeSparsGraph_)
        visual_->viz2Trigger();
    }

    // Run SPARS checks
    GuardType addReason;         // returns why the state was added
    SparseVertex newVertex = 0;  // the newly generated sparse vertex
    if (!addStateToRoadmap(vertexIt->v_, newVertex, addReason))
    {
      // std::cout << "Failed AGAIN to add state to roadmap------" << std::endl;

      // Visualize the failed vertex as a small red dot
      if (visualizeSparsGraph_)
      {
        visual_->viz2State(getDenseState(vertexIt->v_), /*small red*/ 3, 0);
      }
      vertexIt++;  // increment since we didn't remove anything from the list
    }
    else
    {
      // If a new vertex was created, update its popularity property
      if (newVertex)  // value is not null, so it was created
      {
        // std::cout << "SETTING POPULARITY of vertex " << newVertex << " to value "
        //<< vertexInsertionOrder[i].weight_ << std::endl;

        // Update popularity
        vertexPopularity_[newVertex] = vertexIt->weight_;
      }

      // Remove this state from the candidates for insertion vector
      vertexIt = vertexInsertionOrder.erase(vertexIt);

      sucessfulInsertions++;
    }
  }  // end for

  return true;
}

void SparseDB::getVertexInsertionOrdering(std::list<WeightedVertex> &vertexInsertionOrder)
{
  if (sparseCreationInsertionOrder_ == 0)
  {
    OMPL_INFORM("Creating sparse graph using popularity ordering");
    getPopularityOrder(vertexInsertionOrder);  // Create SPARs graph in order of popularity
  }
  else if (sparseCreationInsertionOrder_ == 1)
  {
    OMPL_WARN("Creating sparse graph using default ordering");
    getDefaultOrder(vertexInsertionOrder);
  }
  else if (sparseCreationInsertionOrder_ == 2)
  {
    OMPL_WARN("Creating sparse graph using random ordering");
    getRandomOrder(vertexInsertionOrder);
  }
  else
  {
    OMPL_ERROR("Unknown insertion order method");
    exit(-1);
  }
}

void SparseDB::eliminateDisjointSets()
{
  std::size_t coutIndent = 2;
  if (disjointVerbose_)
    std::cout << std::string(coutIndent, ' ') << "eliminateDisjointSets()" << std::endl;

  visualizeOverlayNodes_ = true;  // visualize all added nodes in a separate window, also
  specialMode_ = false;           // disable parts of addRoadmap
  bool verbose = false;

  /** \brief Sampler user for generating valid samples in the state space */
  base::ValidStateSamplerPtr validSampler = si_->allocValidStateSampler();

  // For each dense vertex we add
  std::size_t numSets = 2;           // dummy value that will be updated at first loop
  std::size_t addedStatesCount = 0;  // count how many states we add
  while (numSets > 1)
  {
    // Add dense vertex
    base::State *state = si_->allocState();
    DenseVertex denseV = denseDB_->addVertex(state, COVERAGE);  // TODO(davetcoleman): COVERAGE means nothing

    // For each random sample
    while (true)
    {
      // Sample randomly
      if (!validSampler->sample(state))  // TODO(davetcoleman): is it ok with the nn_ to change the state after
        // having added it to the nearest neighbor??
      {
        OMPL_ERROR("Unable to find valid sample");
        exit(-1);  // this should never happen
      }

      OMPL_WARN("Not setup for RobotModelStateSpace!!");
      ob::RealVectorStateSpace::StateType *real_state = static_cast<ob::RealVectorStateSpace::StateType *>(state);
      real_state->values[2] = 0;  // ensure task level is 0, TODO

      // Debug
      if (verbose && false)
      {
        visual_->viz4State(state, /*small red*/ 3, 0);
        visual_->viz4Trigger();
        usleep(0.001 * 1000000);
      }

      // Run SPARS checks
      GuardType addReason;     // returns why the state was added
      SparseVertex newVertex;  // the newly generated sparse vertex
      if (addStateToRoadmap(denseV, newVertex, addReason))
      {
        // Debug
        if (disjointVerbose_)
          std::cout << std::string(coutIndent + 2, ' ') << "Added random sampled state to fix graph connectivity, "
            "total new states: " << ++addedStatesCount << std::endl;

        // Visualize
        if (visualizeSparsGraph_)
        {
          visual_->viz2Trigger();  // show the new sparse graph addition
        }

        // Attempt to re-add neighbors from dense graph into sparse graph that have not been added yet
        // so that the new vertex has good edges
        if (addReason == COVERAGE)
        {
          if (reinsertNeighborsIntoSpars(newVertex))
          {
            std::cout << "success in reinserting COVERAGE node neighbors " << std::endl;
          }
          else
            std::cout << "no success in reinserting COVERAGE node neighbors " << std::endl;
        }

        // Attempt to connect new DENSE vertex into dense graph by connecting neighbors
        denseDB_->connectNewVertex(denseV);

        // Statistics
        numSamplesAddedForDisjointSets_++;

        break;  // must break so that a new state can be allocated
      }
    }

    // Update number of sets
    numSets = getDisjointSetsCount();

    // Debug
    if (disjointVerbose_)
      std::cout << std::string(coutIndent + 4, ' ') << "Remaining disjoint sets: " << numSets << std::endl;
  }  // end while

  // getDisjointSetsCount(true);  // show in verbose mode
  // OMPL_ERROR("  Shutting down because there should only be one disjoint set");
  // exit(-1);
}

bool SparseDB::reinsertNeighborsIntoSpars(SparseVertex newVertex)
{
  // Nodes near our input state
  std::vector<DenseVertex> graphNeighborhood;
  // Visible nodes near our input state
  std::vector<DenseVertex> visibleNeighborhood;

  // Convert sparse to dense vertex
  DenseVertex denseV = denseVertexProperty_[newVertex];

  // Find nearby nodes
  denseDB_->findGraphNeighbors(denseV, graphNeighborhood, visibleNeighborhood, sparseDelta_, 4);

  bool result = false;

  // Attempt to re-add visible neighbors
  for (std::size_t i = 0; i < visibleNeighborhood.size(); ++i)
  {
    std::cout << "attempting to reinsert " << i << " - " << visibleNeighborhood[i] << std::endl;
    SparseVertex newVertex;
    GuardType addReason;
    if (addStateToRoadmap(visibleNeighborhood[i], newVertex, addReason))
    {
      std::cout << "    addition worked!! " << std::endl;
      result = true;
    }
  }

  return result;
}

std::size_t SparseDB::getDisjointSetsCount(bool verbose)
{
  std::size_t numSets = 0;
  foreach (SparseVertex v, boost::vertices(g_))
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

bool SparseDB::getPopularityOrder(std::list<WeightedVertex> &vertexInsertionOrder)
{
  bool verbose = false;

  // Error check
  BOOST_ASSERT_MSG(denseDB_->getNumVertices() > queryVertices_.size(), "Unable to get vertices in order of popularity because dense "
                   "graph is empty");

  if (visualizeNodePopularity_)  // Clear visualization
  {
    visual_->viz3State(NULL, /* type = deleteAllMarkers */ 0, 0);
  }

  // Sort the verticies by popularity in a queue
  std::priority_queue<WeightedVertex, std::vector<WeightedVertex>, CompareWeightedVertex> pqueue;

  // Loop through each popular edge in the dense graph
  foreach (DenseVertex v, boost::vertices(denseDB_->g_))
  {
    // Do not process the search vertex, it is null
    if (v <= queryVertices_.back())
      continue;

    if (verbose)
      std::cout << "Vertex: " << v << std::endl;
    double popularity = 0;
    foreach (DenseEdge edge, boost::out_edges(v, denseDB_->g_))
    {
      if (verbose)
        std::cout << "  Edge: " << edge << std::endl;
      popularity += (100 - denseDB_->edgeWeightProperty_[edge]);
    }
    if (verbose)
      std::cout << "  Total popularity: " << popularity << std::endl;
    pqueue.push(WeightedVertex(v, popularity));
  }

  // Remember which one was the largest
  double largestWeight = pqueue.top().weight_;
  if (largestWeight == 0)
    largestWeight = 1;  // prevent division by zero

  if (verbose)
    std::cout << "Largest weight: " << largestWeight << std::endl
              << std::endl;

  // Convert pqueue into vector
  while (!pqueue.empty())  // Output the vertices in order
  {
    vertexInsertionOrder.push_back(pqueue.top());

    // Modify the weight to be a percentage of the max weight
    const double weightPercent = pqueue.top().weight_ / largestWeight * 100.0;
    vertexInsertionOrder.back().weight_ = weightPercent;

    // Visualize
    if (visualizeNodePopularity_)
    {
      if (verbose)
        std::cout << "vertex " << pqueue.top().v_ << ", mode 7, weightPercent " << weightPercent << std::endl;
      visual_->viz3State(denseDB_->stateProperty_[pqueue.top().v_], /*value0-100*/ 7, weightPercent);
    }

    // Remove from priority queue
    pqueue.pop();
  }

  if (visualizeNodePopularity_)
  {
    visual_->viz3Trigger();
    usleep(0.001 * 1000000);
  }

  return true;
}

bool SparseDB::getDefaultOrder(std::list<WeightedVertex> &vertexInsertionOrder)
{
  bool verbose = false;
  double largestWeight = -1 * std::numeric_limits<double>::infinity();

  // Loop through each popular edge in the dense graph
  foreach (DenseVertex v, boost::vertices(denseDB_->g_))
  {
    // Do not process the search vertex, it is null
    if (v <= queryVertices_.back())
      continue;

    if (verbose)
      std::cout << "Vertex: " << v << std::endl;
    double popularity = 0;

    foreach (DenseEdge edge, boost::out_edges(v, denseDB_->g_))
    {
      if (verbose)
        std::cout << "  Edge: " << edge << std::endl;
      popularity += (100 - denseDB_->edgeWeightProperty_[edge]);
    }
    if (verbose)
      std::cout << "  Total popularity: " << popularity << std::endl;

    // Record
    vertexInsertionOrder.push_back(WeightedVertex(v, popularity));

    // Track largest weight
    if (popularity > largestWeight)
      largestWeight = popularity;
  }

  // Update the weights
  for (WeightedVertex wv : vertexInsertionOrder)
  {
    // Modify the weight to be a percentage of the max weight
    const double weightPercent = wv.weight_ / largestWeight * 100.0;
    wv.weight_ = weightPercent;

    // Visualize vertices
    if (visualizeNodePopularity_)
    {
      visual_->viz3State(denseDB_->stateProperty_[wv.v_], /*value0-100*/ 7, weightPercent);
    }
  }

  // Visualize vertices
  if (visualizeNodePopularity_)
  {
    visual_->viz3Trigger();
    usleep(0.001 * 1000000);
  }

  return true;
}

bool SparseDB::getRandomOrder(std::list<WeightedVertex> &vertexInsertionOrder)
{
  std::list<WeightedVertex> defaultVertexInsertionOrder;
  getDefaultOrder(defaultVertexInsertionOrder);

  // Create vector and copy list into it
  std::vector<WeightedVertex> tempVector(defaultVertexInsertionOrder.size());
  std::copy(defaultVertexInsertionOrder.begin(), defaultVertexInsertionOrder.end(), tempVector.begin());

  // Randomize
  std::random_shuffle(tempVector.begin(), tempVector.end());

  // Convert back to list
  vertexInsertionOrder.resize(tempVector.size());
  std::copy(tempVector.begin(), tempVector.end(), vertexInsertionOrder.begin());

  return true;
}

int SparseDB::iRand(int min, int max)
{
  int n = max - min + 1;
  int remainder = RAND_MAX % n;
  int x;
  do
  {
    x = rand();
  } while (x >= RAND_MAX - remainder);
  return min + x % n;
}

bool SparseDB::addStateToRoadmap(DenseVertex denseV, SparseVertex &newVertex, GuardType &addReason)
{
  std::size_t coutIndent = 2;
  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "addStateToRoadmap() DenseV: " << denseV << std::endl;

  bool stateAdded = false;

  // Nodes near our input state
  std::vector<SparseVertex> graphNeighborhood;
  // Visible nodes near our input state
  std::vector<SparseVertex> visibleNeighborhood;

  // Find nearby nodes
  findGraphNeighbors(denseV, graphNeighborhood, visibleNeighborhood, coutIndent + 4);

  // Always add a node if no other nodes around it are visible (GUARD)
  if (checkAddCoverage(denseV, visibleNeighborhood, newVertex, coutIndent + 4))
  {
    if (checksVerbose_)
      std::cout << "Added: COVERAGE " << std::endl;
    addReason = COVERAGE;
    stateAdded = true;
  }
  else if (checkAddConnectivity(denseV, visibleNeighborhood, newVertex, coutIndent + 8))
  {
    if (checksVerbose_)
      std::cout << "Added: CONNECTIVITY " << std::endl;
    addReason = CONNECTIVITY;
    stateAdded = true;
  }
  else if (checkAddInterface(denseV, graphNeighborhood, visibleNeighborhood, newVertex, coutIndent + 12))
  {
    if (checksVerbose_)
      std::cout << "Added: INTERFACE " << std::endl;
    addReason = INTERFACE;
    stateAdded = true;
  }

  return stateAdded;
}

bool SparseDB::checkAddCoverage(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                SparseVertex &newVertex, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "checkAddCoverage() Are other nodes around it visible?" << std::endl;

  // Only add a node for coverage if it has no neighbors
  if (visibleNeighborhood.size() > 0)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "NOT adding node for coverage " << std::endl;
    return false;  // has visible neighbors
  }

  // No free paths means we add for coverage
  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "Adding node for COVERAGE " << std::endl;

  newVertex = addVertex(denseV, COVERAGE);

  if (specialMode_)
  {
    std::cout << "New sparse vertex " << newVertex << std::endl;
  }

  // Note: we do not connect this node with any edges because we have already determined
  // it is too far away from any nearby nodes

  return true;
}

bool SparseDB::checkAddConnectivity(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                    SparseVertex &newVertex, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "checkAddConnectivity() Does this node connect neighboring nodes "
      "that are not connected? " << std::endl;

  // If less than 2 neighbors there is no way to find a pair of nodes in different connected components
  if (visibleNeighborhood.size() < 2)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "NOT adding node for connectivity" << std::endl;
    return false;
  }

  // Identify visibile nodes around our new state that are unconnected (in different connected components)
  // and connect them
  std::set<SparseVertex> statesInDiffConnectedComponents;

  // For each neighbor
  for (std::size_t i = 0; i < visibleNeighborhood.size(); ++i)
  {
    // For each other neighbor
    for (std::size_t j = i + 1; j < visibleNeighborhood.size(); ++j)
    {
      // If they are in different components
      if (!sameComponent(visibleNeighborhood[i], visibleNeighborhood[j]))
      {
        statesInDiffConnectedComponents.insert(visibleNeighborhood[i]);
        statesInDiffConnectedComponents.insert(visibleNeighborhood[j]);
      }
    }
  }

  // Were any disconnected states found?
  if (statesInDiffConnectedComponents.empty())
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "No states in diff connected components found " << std::endl;
    return false;
  }

  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "Adding node for CONNECTIVITY " << std::endl;

  // Add the node
  newVertex = addVertex(denseV, CONNECTIVITY);

  // Add the edges
  for (std::set<SparseVertex>::const_iterator vertexIt = statesInDiffConnectedComponents.begin();
       vertexIt != statesInDiffConnectedComponents.end(); ++vertexIt)
  {
    // Do not add edge from self to self
    if (si_->getStateSpace()->equalStates(getSparseStateConst(*vertexIt), getSparseStateConst(newVertex)))
    {
      std::cout << "Prevented same vertex from being added twice " << std::endl;
      continue;  // skip this pairing
    }

    if (checksVerbose_)
      std::cout << std::string(coutIndent + 3, ' ') << "Loop: Adding vertex " << *vertexIt << std::endl;

    // New vertex should not be connected to anything - there's no edge between the two states
    if (boost::edge(newVertex, *vertexIt, g_).second == true)
    {
      OMPL_ERROR("Somehow the new vertex %u is already connected to old vertex %u", newVertex, *vertexIt);
      exit(-1);
    }

    // The components haven't been united by previous edges created in this for loop
    if (!sameComponent(*vertexIt, newVertex))
    {
      std::size_t visualColor = 0;  // GREEN
      if (secondSparseInsertionAttempt_)
        visualColor = 25;  // ORANGE

      // Connect
      addEdge(newVertex, *vertexIt, visualColor, coutIndent + 4);
    }
    else
    {
      // This is not a big deal
      // OMPL_WARN("Two states that where not prev in the same component were joined during the same for "
      //"loop");
    }
  }

  return true;
}

bool SparseDB::checkAddInterface(DenseVertex denseV, std::vector<SparseVertex> &graphNeighborhood,
                                 std::vector<SparseVertex> &visibleNeighborhood, SparseVertex &newVertex,
                                 std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "checkAddInterface() Does this node's neighbor's need it to "
      "better connect them?" << std::endl;

  // If we have at least 2 neighbors
  // TODO(davetcoleman): why only the first two nodes??
  if (visibleNeighborhood.size() < 2)
  {
    return false;
  }

  std::size_t visualColor = 50;  // Yellow
  if (specialMode_)
    visualColor = 100;  // RED
  else if (secondSparseInsertionAttempt_)
    visualColor = 75;  // ORANGE

  // If the two closest nodes are also visible
  const std::size_t threadID = 0;
  if (graphNeighborhood[0] == visibleNeighborhood[0] && graphNeighborhood[1] == visibleNeighborhood[1])
  {
    // If our two closest neighbors don't share an edge
    if (!boost::edge(visibleNeighborhood[0], visibleNeighborhood[1], g_).second)
    {
      // If they can be directly connected
      if (edgeCache_->checkMotionWithCache(denseVertexProperty_[visibleNeighborhood[0]],
                                           denseVertexProperty_[visibleNeighborhood[1]], threadID))
        // if (si_->checkMotion(getSparseStateConst(visibleNeighborhood[0]), getSparseStateConst(visibleNeighborhood[1])))
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 2, ' ') << "INTERFACE: directly connected nodes" << std::endl;

        // Connect them
        addEdge(visibleNeighborhood[0], visibleNeighborhood[1], visualColor, coutIndent + 4);
      }
      else  // They cannot be directly connected
      {
        // Add the new node to the graph, to bridge the interface
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 2, ' ') << "Adding node for INTERFACE" << std::endl;

        newVertex = addVertex(denseV, INTERFACE);
        addEdge(newVertex, visibleNeighborhood[0], visualColor, coutIndent + 4);
        addEdge(newVertex, visibleNeighborhood[1], visualColor, coutIndent + 4);
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 2, ' ') << "INTERFACE: connected two neighbors through new "
            "interface node" << std::endl;
      }
      // Report success
      return true;
    }
    else
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 2, ' ') << "Two neighbors already share an edge, not "
          "connecting" << std::endl;
    }
  }

  return false;
}

void SparseDB::getInterfaceNeighborhood(DenseVertex denseV, std::vector<DenseVertex> &interfaceNeighborhood,
                                        std::size_t coutIndent)
{
  if (fourthCheckVerbose_)
    std::cout << std::string(coutIndent, ' ') << "getInterfaceNeighborhood()" << std::endl;

  // Error check
  assert(denseDelta_ > std::numeric_limits<double>::epsilon());

  // Get dense vertex's representative sparse vertex
  SparseVertex rep = denseDB_->representativesProperty_[denseV];

  // For each neighbor we are connected to
  foreach (DenseVertex neighborDenseV, boost::adjacent_vertices(denseV, denseDB_->g_))
  {
    // If neighbor representative is not our own
    if (denseDB_->representativesProperty_[neighborDenseV] != rep)
    {
      std::cout << "          Distance: " << std::setprecision(4) << denseDB_->distanceFunction(denseV, neighborDenseV)
                << " denseDelta_: " << denseDelta_ << std::endl;

      // If it is within denseDelta_
      if (denseDB_->distanceFunction(denseV, neighborDenseV) < denseDelta_)
      {
        // Append him to the list
        interfaceNeighborhood.push_back(neighborDenseV);
      }
    }
  }
}

void SparseDB::findGraphNeighbors(DenseVertex v1, std::vector<SparseVertex> &graphNeighborhood,
                                  std::vector<SparseVertex> &visibleNeighborhood, std::size_t coutIndent)
{
  const bool verbose = false;

  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "findGraphNeighbors() DenseV: " << v1 << std::endl;

  base::State *state = getDenseState(v1);

  // Search
  const std::size_t threadID = 0;
  getSparseState(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighborhood);
  getSparseState(queryVertices_[threadID]) = NULL;

  // Now that we got the neighbors from the NN, we must remove any we can't see
  for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  {
    DenseVertex v2 = denseVertexProperty_[graphNeighborhood[i]];
    // Don't collision check if they are the same dense vertex
    if (v1 != v2)
    {
      // Only collision check motion if they don't already share an edge in the dense graph
      if (!boost::edge(v1, v2, denseDB_->g_).second)
      {
        if (!edgeCache_->checkMotionWithCache(v1, v2, threadID))
        {
          continue;
        }
      }
      else if (verbose)
      {
        std::cout << "Skipping check motion because already share an edge in dense graph! " << std::endl;

        // Visualize for checking that this is true
        visual_->viz3State(NULL, /* type = deleteAllMarkers */ 0, 0);
        visual_->viz3State(state, 2 /*small blue*/, 0);
        visual_->viz3State(getSparseState(graphNeighborhood[i]), 2 /*small blue*/, 0);
        visual_->viz3Trigger();
        usleep(1 * 1000000);
      }
    }
    else if (verbose)
      std::cout << " ---- Skipping collision checking because same vertex " << std::endl;

    // The two are visible to each other!
    visibleNeighborhood.push_back(graphNeighborhood[i]);
  }

  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "Graph neighborhood: " << graphNeighborhood.size()
              << " | Visible neighborhood: " << visibleNeighborhood.size() << std::endl;
}

bool SparseDB::sameComponent(SparseVertex v1, SparseVertex v2)
{
  return boost::same_component(v1, v2, disjointSets_);
}

SparseVertex SparseDB::addVertex(DenseVertex denseV, const GuardType &type)
{
  // Create vertex
  SparseVertex v = boost::add_vertex(g_);

  // Add properties
  typePropertySparse_[v] = type;
  denseVertexProperty_[v] = denseV;
  vertexPopularity_[v] = MAX_POPULARITY_WEIGHT;  // 100 means the vertex is very unpopular

  // Connected component tracking
  disjointSets_.make_set(v);

  // Add vertex to nearest neighbor structure
  nn_->add(v);

  // Visualize
  if (visualizeSparsGraph_)
  {
    std::size_t mode = 5;  // default
    if (type == CONNECTIVITY)
      mode = 5;  // black, large
    else if (type == COVERAGE)
      mode = 4;  // PURPLE, regular
    else if (type == INTERFACE)
      mode = 8;  // RED, large
    else
      OMPL_ERROR("Unknown mode");

    visual_->viz2State(getSparseState(v), mode, sparseDelta_);
    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    if (visualizeOverlayNodes_)  // after initial spars graph is created, show additions not from grid
    {
      visual_->viz3State(getSparseState(v), mode, sparseDelta_);
      visual_->viz3Trigger();
      usleep(0.001 * 1000000);
    }
  }

  return v;
}

std::size_t SparseDB::getVizVertexType(const GuardType &type)
{
  switch (type)
  {
    case COVERAGE:
      return 1;
    case CONNECTIVITY:
      return 2;
    case INTERFACE:
      return 3;
    case QUALITY:
      return 4;
    case START:
    case GOAL:
    case CARTESIAN:
      OMPL_ERROR("Type: %u not handled yet in getVizVertexType()", type);
      return 5;
  }
  OMPL_ERROR("Unknown vertex type: %u", type);
  return 5;
}

void SparseDB::addEdge(SparseVertex v1, SparseVertex v2, std::size_t visualColor, std::size_t coutIndent)
{
  assert(v1 <= getNumVertices());
  assert(v2 <= getNumVertices());
  assert(v1 != v2);
  BOOST_ASSERT_MSG(getSparseState(v1) != getSparseState(v2), "States on both sides of an edge are the same");

  // std::cout << std::string(coutIndent, ' ') << "addEdge: Connecting vertex " << v1 << " to vertex " << v2 <<
  // std::endl;

  // Create the new edge
  SparseEdge e = (boost::add_edge(v1, v2, g_)).first;

  // Add associated properties to the edge
  edgeWeightPropertySparse_[e] = distanceFunction(v1, v2);  // TODO: use this value with astar
  edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;

  // Add the edge to the incrementeal connected components datastructure
  disjointSets_.union_set(v1, v2);

  // Visualize
  if (visualizeSparsGraph_)
  {
    /* Color Key:
       0   - GREEN  - connectivity
       25  - LIGHT GREEN - connectivity second round
       50  - YELLOW - interface
       75  - ORANGE - interface second round
       100 - RED    - interface special
    */
    visual_->viz2Edge(getSparseState(v1), getSparseState(v2), visualColor);
    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    if (visualizeOverlayNodes_)  // after initial spars graph is created, show additions not from grid
    {
      visual_->viz3Edge(getSparseState(v1), getSparseState(v2), visualColor);
      visual_->viz3Trigger();
      usleep(0.001 * 1000000);
    }
  }
}

base::State *&SparseDB::getSparseState(SparseVertex v)
{
  return denseDB_->stateProperty_[denseVertexProperty_[v]];
}

const base::State *SparseDB::getSparseStateConst(SparseVertex v) const
{
  return denseDB_->stateProperty_[denseVertexProperty_[v]];
}

base::State *&SparseDB::getDenseState(DenseVertex denseV)
{
  return denseDB_->stateProperty_[denseV];
}

void SparseDB::displaySparseDatabase(bool showVertices)
{
  OMPL_INFORM("Displaying Sparse database");

  // Error check
  if (getNumVertices() == 0 || getNumEdges() == 0)
  {
    OMPL_ERROR("Unable to show database because no vertices/edges available");
    exit(-1);
  }

  // Clear previous visualization
  visual_->viz2State(NULL, /*deleteAllMarkers=*/0, 1);

  if (visualizeDatabaseEdges_)
  {
    // Loop through each edge
    std::size_t count = 0;
    std::size_t debugFrequency = getNumEdges() / 10;
    std::cout << "Displaying sparse edges: " << std::flush;
    foreach (SparseEdge e, boost::edges(g_))
    {
      // Add edge
      SparseVertex v1 = boost::source(e, g_);
      SparseVertex v2 = boost::target(e, g_);

      // TODO(davetcoleman): currently the weight property is not normalized for 0-100 scale so not using for
      // visualization
      // Visualize
      visual_->viz2Edge(getSparseStateConst(v1), getSparseStateConst(v2), 100);  // edgeWeightPropertySparse_[e]);

      // Prevent cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumEdges()) * 100.0
                  << "% " << std::flush;
        visual_->viz2Trigger();
      }

      count++;
    }
    std::cout << std::endl;
  }

  if (visualizeDatabaseVertices_)
  {
    // Loop through each vertex
    std::size_t count = 0;
    std::size_t debugFrequency = getNumVertices() / 10;
    std::cout << "Displaying sparse vertices: " << std::flush;
    foreach (SparseVertex v, boost::vertices(g_))
    {
      // Check for null states
      if (getSparseStateConst(v))
      {
        visual_->viz2State(getSparseStateConst(v), /*small blue*/ 2, 1);
        // visual_->viz2State(getSparseStateConst(v), /*popularity0-100*/ 7, vertexPopularity_[v]);
      }
      else if (v > queryVertices_.back())  // query vertex should always be null, actually
        OMPL_WARN("Null sparse state found on vertex %u", v);

      // Prevent cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumVertices()) * 100.0
                  << "% " << std::flush;
        visual_->viz2Trigger();
      }
      count++;
    }
    std::cout << std::endl;
  }
  // else

  // Publish remaining edges
  visual_->viz2Trigger();
}

double SparseDB::distanceFunction(const SparseVertex a, const SparseVertex b) const
{
  // const double dist = si_->distance(getSparseState(a), getSparseState(b));
  // std::cout << "getting distance from " << a << " to " << b << " of value " << dist << std::endl;
  // return dist;
  return si_->distance(getSparseStateConst(a), getSparseStateConst(b));
}


}  // namespace bolt
}  // namespace tools
}  // namespace ompl
