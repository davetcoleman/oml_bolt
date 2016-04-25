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
SparseDB::SparseDB(base::SpaceInformationPtr si, DenseDB *denseDB, VisualizerPtr visual, EdgeCachePtr edgeCache)
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
  , interfaceDataProperty_(boost::get(vertex_interface_data_t(), g_))
  , vertexPopularity_(boost::get(vertex_popularity_t(), g_))
  // Disjoint set accessors
  , disjointSets_(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_))
  , nearSamplePoints_(2 * si_->getStateDimension())
{
  // Add search state
  initializeQueryState();

  // Initialize nearest neighbor datastructure
  nn_.reset(new NearestNeighborsGNAT<SparseVertex>());
  nn_->setDistanceFunction(boost::bind(&otb::SparseDB::distanceFunction, this, _1, _2));

  // Initialize path simplifier
  pathSimplifier_.reset(new geometric::PathSimplifier(si_));
  pathSimplifier_->freeStates(false);
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

  if (!sampler_)
    sampler_ = si_->allocValidStateSampler();

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
    return;  // assume its already been setup
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

void SparseDB::preprocessSPARS()
{
  std::size_t numThreads = boost::thread::hardware_concurrency();

  // Setup so ensure sparseDelta_ is correct
  setup();

  // Increase the sparseDelta
  sparseDelta_ = getSecondarySparseDelta();

  // Debugging
  if (false)
  {
    OMPL_WARN("Overriding number of threads for testing to 1");
    numThreads = 1;
  }

  // Setup threading
  std::vector<boost::thread *> threads(numThreads);
  OMPL_INFORM("Preprocessing SPRS graph using %u threads", numThreads);

  // For each thread
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    threads[i] = new boost::thread(boost::bind(&SparseDB::preprocessSPARSThread, this, i, numThreads, si));
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }

  denseDB_->graphUnsaved_ = true;
}

void SparseDB::preprocessSPARSThread(std::size_t threadID, std::size_t numThreads, base::SpaceInformationPtr si)
{
  const bool verbose = true;

  // Choose start vertex, skip over the queryVertices
  DenseVertex v1 = queryVertices_.size() + threadID;

  std::size_t numVertices = denseDB_->getNumVertices();  // cache for speed
  std::size_t debugFrequency = std::max(std::size_t(10), static_cast<std::size_t>(numVertices / 200));
  std::size_t saveFrequency = std::max(std::size_t(1000), static_cast<std::size_t>(numVertices / 20));

  if (verbose)
    std::cout << "Thread " << threadID << " starting on vertex " << v1 << " debugFrequency: " << debugFrequency
              << std::endl;

  // Error check
  assert(sparseDelta_ > 0.0001);
  assert(numVertices > 0);
  assert(v1 < numVertices);

  std::vector<DenseVertex> graphNeighborhood;
  while (v1 < numVertices)
  {
    // Feedback
    if (v1 % debugFrequency == 0)
    {
      std::cout << "Preprocessing spars progress: " << (v1 / double(numVertices) * 100.0) << " %" << std::endl;
    }
    if (v1 % saveFrequency == 0)
    {
      edgeCache_->save();
    }

    // Search
    graphNeighborhood.clear();
    denseDB_->nn_->nearestR(v1, sparseDelta_, graphNeighborhood);

    // Collision check edges
    for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
    {
      const DenseVertex v2 = graphNeighborhood[i];
      edgeCache_->checkMotionWithCache(v1, v2, threadID);
    }

    // Increment
    v1 += numThreads;
  }
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
  // else
  // OMPL_WARN("Skipping SPARSE disjoint set fixing");

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
      visual_->viz2State(nullptr, /*deleteAllMarkers*/ 0, 0, 0);
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
    std::cout << "Attempting to insert " << vertexInsertionOrder.size() << " vertices for the " << loopAttempt
              << " loop" << std::endl;

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
      sparseDelta_ = getSecondarySparseDelta();
      std::cout << std::string(coutIndent + 2, ' ') << "sparseDelta_ is now " << sparseDelta_ << std::endl;
      secondSparseInsertionAttempt_ = true;

      std::cout << "temp shutdown before second loop " << std::endl;
      exit(-1);

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

  for (std::list<WeightedVertex>::iterator vertexIt = vertexInsertionOrder.begin();
       vertexIt != vertexInsertionOrder.end();
       /* will manually progress since there are erases*/)
  {
    // User feedback
    if (loopCount++ % debugFrequency == 0)
    {
      std::cout << ANSI_COLOR_BLUE;
      std::cout << std::fixed << std::setprecision(1)
                << "Sparse generation progress: " << (static_cast<double>(loopCount) / originalVertexInsertion) * 100.0
                << "% Cache size: " << edgeCache_->getCacheSize()
                << " Cache usage: " << edgeCache_->getPercentCachedCollisionChecks() << "%" << std::endl;
      std::cout << ANSI_COLOR_RESET;
      if (visualizeSparsGraph_)
        visual_->viz2Trigger();
    }

    // Run SPARS checks
    GuardType addReason;         // returns why the state was added
    SparseVertex newVertex = 0;  // the newly generated sparse vertex

    // TODO(davetcoleman): I would like to multi-thread this but its not worth my time currently
    std::size_t threadID = 0;
    if (!addStateToRoadmap(vertexIt->v_, newVertex, addReason, threadID))
    {
      // std::cout << "Failed AGAIN to add state to roadmap------" << std::endl;

      // Visualize the failed vertex as a small red dot
      if (visualizeSparsGraph_)
      {
        visual_->viz2State(getDenseState(vertexIt->v_), /*small*/ 3, tools::RED, 0);
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
  bool verbose = false;

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
        visual_->viz4State(state, /*small*/ 3, tools::RED, 0);
        visual_->viz4Trigger();
        usleep(0.001 * 1000000);
      }

      // Run SPARS checks
      GuardType addReason;     // returns why the state was added
      SparseVertex newVertex;  // the newly generated sparse vertex
      const std::size_t threadID = 0;
      if (addStateToRoadmap(denseV, newVertex, addReason, threadID))
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
    const std::size_t threadID = 0;
    if (addStateToRoadmap(visibleNeighborhood[i], newVertex, addReason, threadID))
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
  BOOST_ASSERT_MSG(denseDB_->getNumVertices() > queryVertices_.size(),
                   "Unable to get vertices in order of popularity because dense "
                   "graph is empty");

  if (visualizeNodePopularity_)  // Clear visualization
  {
    visual_->viz3State(nullptr, /* type = deleteAllMarkers */ 0, 0, 0);
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
      visual_->viz3State(denseDB_->stateProperty_[pqueue.top().v_], /*value0-100*/ 7, 0, weightPercent);
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
      visual_->viz3State(denseDB_->stateProperty_[wv.v_], /*value0-100*/ 7, 0, weightPercent);
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

bool SparseDB::addStateToRoadmap(DenseVertex denseV, SparseVertex &newVertex, GuardType &addReason,
                                 std::size_t threadID)
{
  std::size_t coutIndent = 0;
  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "addStateToRoadmap() Adding DenseVertex " << denseV << std::endl;

  bool stateAdded = false;
  base::State *workState = si_->allocState();

  // Nodes near our input state
  std::vector<SparseVertex> graphNeighborhood;
  // Visible nodes near our input state
  std::vector<SparseVertex> visibleNeighborhood;

  // Find nearby nodes
  findGraphNeighbors(denseV, graphNeighborhood, visibleNeighborhood, threadID, coutIndent + 4);

  // Always add a node if no other nodes around it are visible (GUARD)
  if (checkAddCoverage(denseV, visibleNeighborhood, newVertex, coutIndent + 4))
  {
    if (checksVerbose_)
      std::cout << "State added for: COVERAGE " << std::endl;
    addReason = COVERAGE;
    stateAdded = true;
  }
  else if (checkAddConnectivity(denseV, visibleNeighborhood, newVertex, coutIndent + 8))
  {
    if (checksVerbose_)
      std::cout << "State added for: CONNECTIVITY " << std::endl;
    addReason = CONNECTIVITY;
    stateAdded = true;
  }
  else if (checkAddInterface(denseV, graphNeighborhood, visibleNeighborhood, newVertex, coutIndent + 12))
  {
    if (checksVerbose_)
      std::cout << "State added for: INTERFACE " << std::endl;
    addReason = INTERFACE;
    stateAdded = true;
  }
  else if (checkAddQuality(denseV, graphNeighborhood, visibleNeighborhood, workState, newVertex, coutIndent + 16))
  {
    if (checksVerbose_)
      std::cout << "State added for: 4th CRITERIA " << std::endl;
    addReason = QUALITY;
    stateAdded = true;

    usleep(5*1000000);
  }
  else
  {
    if (checksVerbose_)
      std::cout << "Did NOT add state for any criteria " << std::endl;
  }

  si_->freeState(workState);

  return stateAdded;
}

bool SparseDB::checkAddCoverage(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                SparseVertex &newVertex, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "checkAddCoverage() Are other nodes around it "
                                                                     "visible?" << ANSI_COLOR_RESET << std::endl;

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

  // Note: we do not connect this node with any edges because we have already determined
  // it is too far away from any nearby nodes

  return true;
}

bool SparseDB::checkAddConnectivity(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                    SparseVertex &newVertex, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ')
              << "checkAddConnectivity() Does this node connect two disconnected components? " << ANSI_COLOR_RESET
              << std::endl;

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
      std::cout << std::string(coutIndent + 2, ' ') << "NOT adding node for connectivity" << std::endl;
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
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ')
              << "checkAddInterface() Does this node's neighbor's need it to better connect them?" << ANSI_COLOR_RESET
              << std::endl;

  // If we have at least 2 neighbors
  // TODO(davetcoleman): why only the first two nodes??
  if (visibleNeighborhood.size() < 2)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "NOT adding node for interface (less than 2 visible neighbors)"
                << std::endl;
    return false;
  }

  std::size_t visualColor = 50;  // Yellow
  if (secondSparseInsertionAttempt_)
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
        std::cout << std::string(coutIndent + 2, ' ') << "Two closest two neighbors already share an edge, not "
                                                         "connecting them" << std::endl;
    }
  }
  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "NOT adding node for interface" << std::endl;
  return false;
}

bool SparseDB::checkAddQuality(DenseVertex denseV, std::vector<SparseVertex> &graphNeighborhood,
                               std::vector<SparseVertex> &visibleNeighborhood, base::State *workState,
                               SparseVertex &newVertex, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "checkAddQuality() Ensure SPARS asymptotic "
                                                                     "optimality" << ANSI_COLOR_RESET << std::endl;

  base::State *candidateState = getDenseState(denseV);

  if (visibleNeighborhood.empty())
  {
    if (checksVerbose_)
      std::cout << "no visible neighbors, not adding 4th criteria " << std::endl;
    return false;
  }

  bool added = false;
  std::map<SparseVertex, base::State *> closeRepresentatives;
  findCloseRepresentatives(workState, candidateState, visibleNeighborhood[0], closeRepresentatives, coutIndent + 2);

  if (checksVerbose_)
    std::cout << std::string(coutIndent, ' ') << "back in checkAddQuality(): Found " << closeRepresentatives.size() << " close representatives"
              << std::endl;

  // For each pair of close representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent+2, ' ') << " Looping through close representatives" << std::endl;

    if (visualizeQualityCriteria_)
    {
      visual_->viz3State(it->second, 5 /*large*/, tools::BLACK, 0);
      visual_->viz3State(getSparseState(it->first), 4 /*Medium, translucent outline*/, tools::PURPLE, sparseDelta_);
      visual_->viz3Edge(getSparseState(it->first), it->second, 0);
      visual_->viz3Trigger();
      usleep(0.1 * 1000000);
    }

    updatePairPoints(visibleNeighborhood[0], candidateState, it->first, it->second, coutIndent+4);
    updatePairPoints(it->first, it->second, visibleNeighborhood[0], candidateState, coutIndent+4);
  }

  // Attempt to find shortest path through closest neighbour
  if (checkAddPath(visibleNeighborhood[0], coutIndent + 4))
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent+2, ' ') << "nearest visible neighbor added for path " << std::endl;
    added = true;
  }

  // Attempt to find shortest path through other pairs of representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent+2, ' ') << "Looping through close representatives to add path" << std::endl;
    checkAddPath(it->first, coutIndent + 4);

    // Delete
    si_->freeState(it->second);
  }

  if(added)
  {
    std::cout << "temp shutdown cause state was added for 4th criteria " << std::endl;
    exit(-1);
  }

  usleep(0.1*1000000);

  return added;
}

bool SparseDB::checkAddPath(SparseVertex v, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ')
              << "checkAddPath() Checks vertex v for short paths through its region and adds when appropriate"
              << ANSI_COLOR_RESET << std::endl;

  bool spannerPropertyWasViolated = false;

  std::vector<SparseVertex> adjVertices;
  foreach (SparseVertex r, boost::adjacent_vertices(v, g_))
    adjVertices.push_back(r);

  std::size_t visualColor = 100;  // red, quality

  // Candidate x vertices as described in the method, filled by function computeX().
  std::vector<SparseVertex> Xs;

  // Candidate v" vertices as described in the method, filled by function computeVPP().
  std::vector<SparseVertex> VPPs;

  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "Vertex " << v << " has " << adjVertices.size() << " adjacent vertices, looping:" << std::endl;

  for (std::size_t i = 0; i < adjVertices.size() && !spannerPropertyWasViolated; ++i)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 4, ' ') << "Checking vertex " << adjVertices[i] << std::endl;

    SparseVertex r = adjVertices[i];

    // Compute all nodes which qualify as a candidate v" for v and vp
    computeVPP(v, r, VPPs, coutIndent+6);

    if (checksVerbose_)
      std::cout << std::string(coutIndent + 4, ' ') << "VPPs size: " << VPPs.size() << std::endl;

    foreach (SparseVertex rp, VPPs)
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 6, ' ') << "Checking VPP vertex " << rp << std::endl;

      // Compute the longest path through the graph
      // Computes all nodes which qualify as a candidate x for v, v', and v"
      computeX(v, r, rp, Xs, coutIndent+8);

      if (checksVerbose_)
        std::cout << std::string(coutIndent + 6, ' ') << "Xs size: " << Xs.size() << std::endl;

      double rm_dist = 0.0;
      foreach (SparseVertex rpp, Xs)
      {
        double tmp_dist = (si_->distance(getSparseState(r), getSparseState(v)) +
                           si_->distance(getSparseState(v), getSparseState(rpp))) /
                          2.0;
        if (tmp_dist > rm_dist)
          rm_dist = tmp_dist;
      }

      InterfaceData &d = getData(v, r, rp);

      // Then, if the spanner property is violated
      if (rm_dist > stretchFactor_ * d.last_distance_)
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 6, ' ') << "Spanner property violated!" << std::endl;

        spannerPropertyWasViolated = true;  // Report that we added for the path
        if (si_->checkMotion(getSparseState(r), getSparseState(rp)))
        {
          if (checksVerbose_)
            std::cout << std::string(coutIndent + 6, ' ') << "Adding edge between r and rp" << std::endl;

          addEdge(r, rp, visualColor, coutIndent + 8);
        }
        else
        {
          geometric::PathGeometric *p = new geometric::PathGeometric(si_);
          if (r < rp)
          {
            p->append(d.sigmaA_);
            p->append(d.pointA_);
            p->append(getSparseState(v));
            p->append(d.pointB_);
            p->append(d.sigmaB_);
          }
          else
          {
            p->append(d.sigmaB_);
            p->append(d.pointB_);
            p->append(getSparseState(v));
            p->append(d.pointA_);
            p->append(d.sigmaA_);
          }

          if (checksVerbose_)
            std::cout << std::string(coutIndent + 6, ' ') << "Creating path" << std::endl;

          pathSimplifier_->reduceVertices(*p, 10);
          pathSimplifier_->shortcutPath(*p, 50);

          if (p->checkAndRepair(100).second)
          {
            SparseVertex prior = r;
            SparseVertex vnew;
            std::vector<base::State *> &states = p->getStates();

            if (checksVerbose_)
              std::cout << std::string(coutIndent + 6, ' ') << "check and repair succeeded" << std::endl;

            foreach (base::State *state, states)
            {
              // no need to clone st, since we will destroy p; we just copy the pointer
              if (checksVerbose_)
                std::cout << std::string(coutIndent + 6, ' ') << "Adding node from shortcut path for QUALITY" << std::endl;
              vnew = addVertex(state, QUALITY);

              addEdge(prior, vnew, visualColor, coutIndent + 8);
              prior = vnew;
            }
            // clear the states, so memory is not freed twice
            states.clear();
            addEdge(prior, rp, visualColor, coutIndent + 8);
          }
          else
            if (checksVerbose_)
              std::cout << std::string(coutIndent + 6, ' ') << "check and repair failed?" << std::endl;

          delete p;
        }
      }
    }
  }

  if (!spannerPropertyWasViolated)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "Spanner property was NOT violated, SKIPPING" << std::endl;
  }

  return spannerPropertyWasViolated;
}

SparseVertex SparseDB::findGraphRepresentative(base::State *state, std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "findGraphRepresentative() " << ANSI_COLOR_RESET
              << std::endl;

  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;

  // Search
  getSparseState(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  getSparseState(queryVertices_[threadID]) = nullptr;

  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "Found " << graphNeighbors.size()
              << " nearest neighbors (graph rep) within SparseDelta " << sparseDelta_ << std::endl;

  SparseVertex result = boost::graph_traits<SparseGraph>::null_vertex();

  for (std::size_t i = 0; i < graphNeighbors.size(); ++i)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "Checking motion of graph representative candidate " << i
                << std::endl;
    if (si_->checkMotion(state, getSparseState(graphNeighbors[i])))
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "graph representative valid " << std::endl;
      result = graphNeighbors[i];
      break;
    }
    else
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "graph representative NOT valid, checking next " << std::endl;
  }
  return result;
}

void SparseDB::findCloseRepresentatives(base::State *workState, const base::State *candidateState,
                                        const SparseVertex candidateRep,
                                        std::map<SparseVertex, base::State *> &closeRepresentatives,
                                        std::size_t coutIndent)
{
  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "findCloseRepresentatives()" << ANSI_COLOR_RESET
              << std::endl;

  assert(closeRepresentatives.empty());

  // denseDelta_ = 0.25 * sparseDelta_;
  // nearSamplePoints_ /= 10; // HACK - this makes it look for the same number of samples as dimensions

  if (checksVerbose_)
    std::cout << std::string(coutIndent + 2, ' ') << "nearSamplePoints: " << nearSamplePoints_ << " denseDelta "
              << denseDelta_ << std::endl;

  if (visualizeQualityCriteria_)
  {
    visual_->viz3State(nullptr, /* type = deleteAllMarkers */ 0, 0, 0);
    visual_->viz3State(candidateState, 4 /*Medium, translucent outline*/, tools::PURPLE, denseDelta_);
    visual_->viz3Trigger();
    usleep(0.1 * 1000000);
  }

  // Then, begin searching the space around new potential state candidateState
  for (unsigned int i = 0; i < nearSamplePoints_; ++i)
  {
    if (checksVerbose_)
      std::cout << std::string(coutIndent + 2, ' ') << "Get supporting representative #" << i << std::endl;

    bool foundValidSample = false;
    static const std::size_t MAX_SAMPLE_ATTEMPT = 1000;
    for (std::size_t attempt = 0; attempt < MAX_SAMPLE_ATTEMPT; ++attempt)
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "Sample attempt " << attempt << std::endl;

      sampler_->sampleNear(workState, candidateState, denseDelta_);
      si_->getStateSpace()->setLevel(workState, 0);

      if (!si_->isValid(workState))
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 6, ' ') << "notValid " << std::endl;

        if (visualizeQualityCriteria_)
          visual_->viz3State(workState, 3 /*small */, tools::RED, 0);

        continue;
      }
      if (si_->distance(candidateState, workState) > denseDelta_)
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 6, ' ') << "Distance too far " << si_->distance(candidateState, workState) << " needs to be less than " << denseDelta_ << std::endl;

        if (visualizeQualityCriteria_)
          visual_->viz3State(workState, 3 /*small */, tools::RED, 0);
        continue;
      }
      if (!si_->checkMotion(candidateState, workState))
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 6, ' ') << "Motion invalid " << std::endl;

        if (visualizeQualityCriteria_)
          visual_->viz3State(workState, 3 /*small */, tools::RED, 0);
        continue;
      }

      if (visualizeQualityCriteria_)
      {
        visual_->viz3State(workState, 1 /*small */, tools::GREEN, 0);
      }
      foundValidSample = true;
      break;
    }  // for each attempt

    visual_->viz3Trigger();
    usleep(0.01 * 1000000);

    if (!foundValidSample)
      std::cout << std::string(coutIndent + 4, ' ') << "Unable to find valid sample after " << MAX_SAMPLE_ATTEMPT
                << " attempts " << std::endl;
    else
      std::cout << std::string(coutIndent + 4, ' ') << "Found valid nearby sample" << std::endl;

    // Compute which sparse vertex represents this new candidate vertex
    SparseVertex representative = findGraphRepresentative(workState, coutIndent + 6);

    // Assuming this sample is actually seen by somebody (which he should be in all likelihood)
    if (representative == boost::graph_traits<SparseGraph>::null_vertex())
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "Sampled state has no representative (is null) " << std::endl;

      // This guy can't be seen by anybody, so we should take this opportunity to add him
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "Adding node for COVERAGE" << std::endl;
      addVertex(si_->cloneState(workState), COVERAGE);

      if (checksVerbose_)
      {
        std::cout << std::string(coutIndent + 4, ' ') << "STOP EFFORTS TO ADD A DENSE PATH" << std::endl;
      }

      // We should also stop our efforts to add a dense path
      for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
           it != closeRepresentatives.end(); ++it)
        si_->freeState(it->second);
      closeRepresentatives.clear();
      break;
    }

    if (checksVerbose_)
      std::cout << std::string(coutIndent + 4, ' ') << "Sampled state has representative (is not null)" << std::endl;

    // If its representative is different than candidateState
    if (candidateRep != representative)
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "candidateRep != representative " << std::endl;

      // And we haven't already tracked this representative
      if (closeRepresentatives.find(representative) == closeRepresentatives.end())
      {
        if (checksVerbose_)
          std::cout << std::string(coutIndent + 4, ' ') << "Track the representative" << std::endl;
        // Track the representative
        closeRepresentatives[representative] = si_->cloneState(workState);
      }
    }
    else
    {
      if (checksVerbose_)
        std::cout << std::string(coutIndent + 4, ' ') << "candidateRep == representative, do not keep this sample " << std::endl;
    }
  }  // for each supporting representative
}

void SparseDB::updatePairPoints(SparseVertex rep, const base::State *q, SparseVertex r, const base::State *s, std::size_t coutIndent)
{
  // First of all, we need to compute all candidate r'
  std::vector<SparseVertex> VPPs;
  computeVPP(rep, r, VPPs, coutIndent);

  // Then, for each pair Pv(r,r')
  foreach (SparseVertex rp, VPPs)
    // Try updating the pair info
    distanceCheck(rep, q, r, s, rp);
}

void SparseDB::computeVPP(SparseVertex v, SparseVertex vp, std::vector<SparseVertex> &VPPs, std::size_t coutIndent)
{
    if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "computeVPP()" << ANSI_COLOR_RESET
              << std::endl;

  VPPs.clear();
  foreach (SparseVertex cvpp, boost::adjacent_vertices(v, g_))
    if (cvpp != vp)
      if (!boost::edge(cvpp, vp, g_).second)
        VPPs.push_back(cvpp);
}

void SparseDB::computeX(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::vector<SparseVertex> &Xs, std::size_t coutIndent)
{
  Xs.clear();

  foreach (SparseVertex cx, boost::adjacent_vertices(vpp, g_))
    if (boost::edge(cx, v, g_).second && !boost::edge(cx, vp, g_).second)
    {
      InterfaceData &d = getData(v, vpp, cx);
      if ((vpp < cx && d.pointA_) || (cx < vpp && d.pointB_))
        Xs.push_back(cx);
    }
  Xs.push_back(vpp);
}

VertexPair SparseDB::index(SparseVertex vp, SparseVertex vpp)
{
  if (vp < vpp)
    return VertexPair(vp, vpp);
  else if (vpp < vp)
    return VertexPair(vpp, vp);

  throw Exception(name_, "Trying to get an index where the pairs are the same point!");
  return VertexPair(0, 0);  // prevent compiler warnings
}

InterfaceData &SparseDB::getData(SparseVertex v, SparseVertex vp, SparseVertex vpp)
{
  return interfaceDataProperty_[v].interfaceHash[index(vp, vpp)];
}

void SparseDB::distanceCheck(SparseVertex rep, const base::State *q, SparseVertex r, const base::State *s,
                             SparseVertex rp)
{
  // Get the info for the current representative-neighbors pair
  InterfaceData &d = getData(rep, r, rp);

  if (r < rp)  // FIRST points represent r (the guy discovered through sampling)
  {
    if (d.pointA_ == nullptr)  // If the point we're considering replacing (P_v(r,.)) isn't there
      // Then we know we're doing better, so add it
      d.setFirst(q, s, si_);
    else  // Otherwise, he is there,
    {
      if (d.pointB_ == nullptr)  // But if the other guy doesn't exist, we can't compare.
      {
        // Should probably keep the one that is further away from rep?  Not known what to do in this case.
        // TODO: is this not part of the algorithm?
      }
      else  // We know both of these points exist, so we can check some distances
          if (si_->distance(q, d.pointB_) < si_->distance(d.pointA_, d.pointB_))
        // Distance with the new point is good, so set it.
        d.setFirst(q, s, si_);
    }
  }
  else  // SECOND points represent r (the guy discovered through sampling)
  {
    if (d.pointB_ == nullptr)  // If the point we're considering replacing (P_V(.,r)) isn't there...
      // Then we must be doing better, so add it
      d.setSecond(q, s, si_);
    else  // Otherwise, he is there
    {
      if (d.pointA_ == nullptr)  // But if the other guy doesn't exist, we can't compare.
      {
        // Should we be doing something cool here?
      }
      else if (si_->distance(q, d.pointA_) < si_->distance(d.pointB_, d.pointA_))
        // Distance with the new point is good, so set it
        d.setSecond(q, s, si_);
    }
  }

  // Lastly, save what we have discovered
  interfaceDataProperty_[rep].interfaceHash[index(r, rp)] = d;
}

void SparseDB::abandonLists(base::State *state)
{
  std::vector<SparseVertex> graphNeighbors;
  std::size_t threadID = 0;

  // Search
  getSparseState(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  getSparseState(queryVertices_[threadID]) = nullptr;

  // For each of the vertices
  foreach (SparseVertex v, graphNeighbors)
  {
    foreach (VertexPair r, interfaceDataProperty_[v].interfaceHash | boost::adaptors::map_keys)
      interfaceDataProperty_[v].interfaceHash[r].clear(si_);
  }
}

void SparseDB::findGraphNeighbors(DenseVertex v1, std::vector<SparseVertex> &graphNeighborhood,
                                  std::vector<SparseVertex> &visibleNeighborhood, std::size_t threadID,
                                  std::size_t coutIndent)
{
  const bool verbose = false;

  if (checksVerbose_)
    std::cout << ANSI_COLOR_GREEN << std::string(coutIndent, ' ') << "findGraphNeighbors() DenseV: " << v1
              << ANSI_COLOR_RESET << std::endl;

  base::State *state = getDenseState(v1);

  // Search
  getSparseState(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighborhood);
  getSparseState(queryVertices_[threadID]) = nullptr;

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
        visual_->viz3State(nullptr, /* type = deleteAllMarkers */ 0, 0, 0);
        visual_->viz3State(state, 2 /*small */, tools::BLUE, 0);
        visual_->viz3State(getSparseState(graphNeighborhood[i]), 2 /*small*/, tools::BLUE, 0);
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

SparseVertex SparseDB::addVertex(base::State *state, const GuardType &type)
{
  // Add the state to the dense graph
  DenseVertex v = denseDB_->addVertex(state, type);

  // Add the DenseVertex to the sparse graph
  return addVertex(v, type);
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
    std::size_t color;
    if (type == COVERAGE)
      color = tools::PURPLE;
    else if (type == CONNECTIVITY)
      color = tools::BLACK;
    else if (type == INTERFACE)
      color = tools::RED;
    else if (type == QUALITY)
      color = tools::BLUE;
    else
      OMPL_ERROR("Unknown mode");

    visual_->viz2State(getSparseState(v), 2 /*medium*/, color, sparseDelta_);
    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    if (visualizeOverlayNodes_)  // after initial spars graph is created, show additions not from grid
    {
      visual_->viz3State(getSparseState(v), 2 /*medium*/, color, sparseDelta_);
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
  visual_->viz2State(nullptr, /*deleteAllMarkers=*/0, 0, 0);

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
        visual_->viz2State(getSparseStateConst(v), /*small*/ 3, tools::BLUE, 1);
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

double SparseDB::getSecondarySparseDelta()
{
  return sparseDelta_ * 1.25;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
