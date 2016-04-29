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
//, nearSamplePoints_(4 * si_->getStateDimension()) // TODO this will probably be pretty slow
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

  // sparseDelta_ = sparseDeltaFraction_ * discretization_;  // sparseDelta should be a multiple of discretization
  sparseDelta_ = sparseDeltaFraction_ * maxExtent_;
  denseDelta_ = denseDeltaFraction_ * maxExtent_;
  OMPL_INFORM("sparseDelta_ = %f", sparseDelta_);
  OMPL_INFORM("denseDelta_ = %f", denseDelta_);

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
  // sparseDelta_ = getSecondarySparseDelta();

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
  std::size_t indent = 2;

  // Clear the old spars graph
  if (getNumVertices() > queryVertices_.size())
  {
    OMPL_INFORM("Resetting sparse database");
    freeMemory();
    initializeQueryState();  // Re-add search state

    if (visualizeSparsGraph_)  // Clear visuals
      visual_->viz2DeleteAllMarkers();
  }

  // Reset parameters
  setup();
  visualizeOverlayNodes_ = false;  // DO NOT visualize all added nodes in a separate window
  useFourthCriteria_ = false;
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
      // sparseDelta_ = getSecondarySparseDelta();
      std::cout << std::string(indent + 2, ' ') << "sparseDelta_ is now " << sparseDelta_ << std::endl;
      secondSparseInsertionAttempt_ = true;

      // std::cout << "temp shutdown before second loop " << std::endl;
      // exit(-1);
      // usleep(2*1000000);
      useFourthCriteria_ = true;

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
      if (visualizeSparsGraph_ && false)
      {
        visual_->viz2State(getDenseState(vertexIt->v_), tools::SMALL, tools::RED, 0);
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
  std::size_t indent = 2;
  if (disjointVerbose_)
    std::cout << std::string(indent, ' ') << "eliminateDisjointSets()" << std::endl;

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
        visual_->viz4State(state, tools::SMALL, tools::RED, 0);
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
          std::cout << std::string(indent + 2, ' ') << "Added random sampled state to fix graph connectivity, "
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
      std::cout << std::string(indent + 4, ' ') << "Remaining disjoint sets: " << numSets << std::endl;
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
    visual_->viz3DeleteAllMarkers();
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
      visual_->viz3State(denseDB_->stateProperty_[pqueue.top().v_], tools::SCALE, tools::BLACK, weightPercent);
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
      visual_->viz3State(denseDB_->stateProperty_[wv.v_], tools::SCALE, tools::BLACK, weightPercent);
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
  std::size_t indent = 0;
  BOLT_DEBUG(indent, "addStateToRoadmap() Adding DenseVertex " << denseV);

  bool stateAdded = false;
  base::State *workState = si_->allocState();

  // Nodes near our input state
  std::vector<SparseVertex> graphNeighborhood;
  // Visible nodes near our input state
  std::vector<SparseVertex> visibleNeighborhood;

  // Find nearby nodes
  findGraphNeighbors(denseV, graphNeighborhood, visibleNeighborhood, threadID, indent + 2);

  // Always add a node if no other nodes around it are visible (GUARD)
  if (checkAddCoverage(denseV, visibleNeighborhood, newVertex, indent + 2))
  {
    if (checksVerbose_)
      std::cout << "State added for: COVERAGE " << std::endl;
    addReason = COVERAGE;
    stateAdded = true;
  }
  else if (checkAddConnectivity(denseV, visibleNeighborhood, newVertex, indent + 6))
  {
    if (checksVerbose_)
      std::cout << "State added for: CONNECTIVITY " << std::endl;
    addReason = CONNECTIVITY;
    stateAdded = true;
  }
  else if (checkAddInterface(denseV, graphNeighborhood, visibleNeighborhood, newVertex, indent + 10))
  {
    if (checksVerbose_)
      std::cout << "State added for: INTERFACE " << std::endl;
    addReason = INTERFACE;
    stateAdded = true;
  }
  else if (useFourthCriteria_ &&
           checkAddQuality(denseV, graphNeighborhood, visibleNeighborhood, workState, newVertex, indent + 14))
  {
    if (checksVerbose_)
      std::cout << "State added for: 4th CRITERIA " << std::endl;
    addReason = QUALITY;
    stateAdded = true;

    // usleep(5 * 1000000);
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
                                SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "checkAddCoverage() Are other nodes around it visible?");

  // Only add a node for coverage if it has no neighbors
  if (visibleNeighborhood.size() > 0)
  {
    BOLT_DEBUG(indent + 2, "NOT adding node for coverage ");
    return false;  // has visible neighbors
  }

  // No free paths means we add for coverage
  BOLT_DEBUG(indent + 2, "Adding node for COVERAGE ");

  newVertex = addVertex(denseV, COVERAGE);

  // Note: we do not connect this node with any edges because we have already determined
  // it is too far away from any nearby nodes

  return true;
}

bool SparseDB::checkAddConnectivity(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                    SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "checkAddConnectivity() Does this node connect two disconnected components?");

  // If less than 2 neighbors there is no way to find a pair of nodes in different connected components
  if (visibleNeighborhood.size() < 2)
  {
    BOLT_DEBUG(indent + 2, "NOT adding node for connectivity");
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
    BOLT_DEBUG(indent + 2, "NOT adding node for connectivity");
    return false;
  }

  BOLT_DEBUG(indent + 2, "Adding node for CONNECTIVITY ");

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

    BOLT_DEBUG(indent + 3, "Loop: Adding vertex " << *vertexIt);

    // New vertex should not be connected to anything - there's no edge between the two states
    if (hasEdge(newVertex, *vertexIt) == true)
    {
      OMPL_ERROR("Somehow the new vertex %u is already connected to old vertex %u", newVertex, *vertexIt);
      exit(-1);
    }

    // The components haven't been united by previous edges created in this for loop
    if (!sameComponent(*vertexIt, newVertex))
    {
      std::size_t visualColor = eGREEN;
      // if (secondSparseInsertionAttempt_)
      //   visualColor = eUGLY_YELLOW;

      // Connect
      addEdge(newVertex, *vertexIt, visualColor, indent + 4);
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
                                 std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "checkAddInterface() Does this node's neighbor's need it to better connect them?");

  // If we have at least 2 neighbors
  // TODO(davetcoleman): why only the first two nodes??
  if (visibleNeighborhood.size() < 2)
  {
    BOLT_DEBUG(indent + 2, "NOT adding node for interface (less than 2 visible neighbors)");
    return false;
  }

  std::size_t visualColor = eYELLOW;
  // if (secondSparseInsertionAttempt_)
  visualColor = eORANGE;

  // If the two closest nodes are also visible
  const std::size_t threadID = 0;
  if (graphNeighborhood[0] == visibleNeighborhood[0] && graphNeighborhood[1] == visibleNeighborhood[1])
  {
    // If our two closest neighbors don't share an edge
    if (!hasEdge(visibleNeighborhood[0], visibleNeighborhood[1]))
    {
      // If they can be directly connected
      if (edgeCache_->checkMotionWithCache(denseVertexProperty_[visibleNeighborhood[0]],
                                           denseVertexProperty_[visibleNeighborhood[1]], threadID))
      // if (si_->checkMotion(getSparseStateConst(visibleNeighborhood[0]), getSparseStateConst(visibleNeighborhood[1])))
      {
        BOLT_DEBUG(indent + 2, "INTERFACE: directly connected nodes");

        // Connect them
        addEdge(visibleNeighborhood[0], visibleNeighborhood[1], visualColor, indent + 4);
      }
      else  // They cannot be directly connected
      {
        // Add the new node to the graph, to bridge the interface
        BOLT_DEBUG(indent + 2, "Adding node for INTERFACE");

        newVertex = addVertex(denseV, INTERFACE);
        addEdge(newVertex, visibleNeighborhood[0], visualColor, indent + 4);
        addEdge(newVertex, visibleNeighborhood[1], visualColor, indent + 4);
        BOLT_DEBUG(indent + 2, "INTERFACE: connected two neighbors through new interface node");
      }

      // Report success
      return true;
    }
    else
    {
      BOLT_DEBUG(indent + 2, "Two closest two neighbors already share an edge, not connecting them");
    }
  }
  BOLT_DEBUG(indent + 2, "NOT adding node for interface");
  return false;
}

bool SparseDB::checkAddQuality(DenseVertex denseV, std::vector<SparseVertex> &graphNeighborhood,
                               std::vector<SparseVertex> &visibleNeighborhood, base::State *workState,
                               SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "checkAddQuality() Ensure SPARS asymptotic optimality");

  if (visibleNeighborhood.empty())
  {
    BOLT_DEBUG(indent + 2, "no visible neighbors, not adding 4th criteria ");
    return false;
  }

  base::State *candidateState = getDenseState(denseV); // paper's name: q
  SparseVertex candidateRep = visibleNeighborhood[0];  // paper's name: v

  if (visualizeQualityCriteria_)
  {
    visual_->viz3DeleteAllMarkers();

    visual_->viz3Edge(candidateState, getSparseState(candidateRep), tools::eORANGE);

    // Show candidate state
    visual_->viz3State(candidateState, tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, denseDelta_);
    visual_->viz3State(candidateState, tools::LARGE, tools::RED, 0);

    // Show candidate state's representative
    visual_->viz3State(getSparseState(candidateRep), tools::LARGE, tools::BLUE, 0);

    // Show candidate state's representative's neighbors
    foreach (SparseVertex adjVertex, boost::adjacent_vertices(candidateRep, g_))
    {
      visual_->viz3Edge(getSparseState(adjVertex), getSparseState(candidateRep), tools::eGREEN);
      visual_->viz3State(getSparseState(adjVertex), tools::LARGE, tools::PURPLE, 0);
    }

    visual_->viz3Trigger();
    usleep(0.001 * 1000000);
  }

  bool added = false;
  std::map<SparseVertex, base::State *> closeRepresentatives; // [nearSampledRep, nearSampledState]
  findCloseRepresentatives(workState, candidateState, candidateRep, closeRepresentatives, indent + 4);

  if (closeRepresentatives.size())
  {
    BOLT_GREEN_DEBUG(indent + 2, "back in checkAddQuality(): Found " << closeRepresentatives.size() << " close representatives ---------------------------");
  }
  else
  {
    BOLT_RED_DEBUG(indent + 2, "back in checkAddQuality(): Found " << closeRepresentatives.size() << " close representatives ---------------------------");
  }

  // For each pair of close representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    BOLT_DEBUG(indent + 4, " Looping through close representatives");
    base::State *nearSampledState = it->second; // paper: q'
    SparseVertex nearSampledRep = it->first; // paper: v'

    if (visualizeQualityCriteria_) // Visualization
    {
      visual_->viz3Edge(getSparseState(nearSampledRep), nearSampledState, tools::eRED);

      visual_->viz3State(nearSampledState, tools::MEDIUM, tools::GREEN, 0);

      // Replicate a regular vertex visualization
      visual_->viz3State(getSparseState(nearSampledRep), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);
      visual_->viz3State(getSparseState(nearSampledRep), tools::LARGE, tools::PURPLE, sparseDelta_);

      visual_->viz3Trigger();
      usleep(0.001 * 1000000);
    }

    // Update interface bookkeeping
    // Process:
    // 1. Get adjacent vertieces of the candidateRep (v) that are unconnected to nearSampledRep (v')
    //      e.g. v''
    // 2. For every adj vertex that is unconnected to v'
    // 3. Check distance:
    //    3.1. Get the interface data stored on vertex candidateRep(v) for max distance between
    //           nearSampledRep (v') and adjVertexUnconnected (v'')
    //    3.2. Add the candidateState (q) and nearSampledState (q') as 'first'

    // Attempt to update bookkeeping for candidateRep (v)
    updatePairPoints(candidateRep, candidateState, nearSampledRep, nearSampledState, indent + 4);

    // ALSO attempt to update bookkeeping for neighboring node nearSampleRep (v')
    updatePairPoints(nearSampledRep, nearSampledState, candidateRep, candidateState, indent + 4);

    // TODO: track if updatePairPoints() made any changes, and skip checkAddPath if it didn't
  }

  BOLT_DEBUG(indent + 2, "Done updating pair points");

  // Attempt to find shortest path through closest neighbour
  if (checkAddPath(candidateRep, indent + 4))
  {
    BOLT_DEBUG(indent + 2, "nearest visible neighbor added for path ");
    added = true;
  }

  // Attempt to find shortest path through other pairs of representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    BOLT_YELLOW_DEBUG(indent + 2, "Looping through close representatives to add path ===============");
    checkAddPath(it->first, indent + 4);

    // Delete
    si_->freeState(it->second);
  }

  if (added)
  {
    std::cout << "temp sleep cause state was added for 4th criteria " << std::endl;
    //usleep(5*1000000);
    exit(0);
  }

  usleep(0.001 * 1000000);

  return added;
}

bool SparseDB::checkAddPath(SparseVertex candidateRep, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, "checkAddPath()");
  bool spannerPropertyWasViolated = false;

  // Candidate x vertices as described in the method, filled by function getMaxSpannerPath().
  std::vector<SparseVertex> qualifiedVertices;

  // Candidate v" vertices as described in the method, filled by function getAdjVerticesOfV1UnconnectedToV2().
  std::vector<SparseVertex> adjVerticesUnconnected;

  // Copy adjacent vertices into vector because we might add additional edges during this function
  std::vector<SparseVertex> adjVertices;
  foreach (SparseVertex adjVertex, boost::adjacent_vertices(candidateRep, g_))
    adjVertices.push_back(adjVertex);

  BOLT_DEBUG(indent + 2, "Vertex " << candidateRep << " has " << adjVertices.size() << " adjacent vertices, looping:");

  // Loop through adjacent vertices
  for (std::size_t i = 0; i < adjVertices.size() && !spannerPropertyWasViolated; ++i)
  {
    SparseVertex vp = adjVertices[i];  // vp = v' from paper

    BOLT_DEBUG(indent + 4, "Checking v' = " << vp);

    if (visualizeQualityCriteria_) // visualize
    {
      visual_->viz5DeleteAllMarkers();

      // Show edge between them
      visual_->viz5Edge(getSparseState(vp), getSparseState(candidateRep), tools::eGREEN);

      // Show candidate rep
      visual_->viz5State(getSparseState(candidateRep), tools::LARGE, tools::BLUE, 0);

      // Show adjacent state
      visual_->viz5State(getSparseState(vp), tools::LARGE, tools::PURPLE, 0);
    }

    // Compute all nodes which qualify as a candidate v" for v and vp
    getAdjVerticesOfV1UnconnectedToV2(candidateRep, vp, adjVerticesUnconnected, indent + 6);

    // for each vertex v'' that is adjacent to v (has a valid edge) and does not share an edge with v'
    foreach (SparseVertex vpp, adjVerticesUnconnected)  // vpp = v'' from paper
    {
      BOLT_DEBUG(indent + 6, "Checking v'' = " << vpp);

      if (visualizeQualityCriteria_)
      {
        // Show adjacent states
        visual_->viz5State(getSparseState(vpp), tools::LARGE, tools::LIME_GREEN, 0);
      }

      // Computes all nodes which qualify as a candidate x for v, v', and v"
      getMaxSpannerPath(candidateRep, vp, vpp, qualifiedVertices, indent + 8);

      // Find the maximum spanner distance
      BOLT_DEBUG(indent + 6, "Find the maximum spanner distance between v' and v''");
      double maxDist = 0.0;
      foreach (SparseVertex qualifiedVertex, qualifiedVertices)
      {
        BOLT_DEBUG(indent + 8, "Vertex: " << qualifiedVertex);
        if (visualizeQualityCriteria_)
          visual_->viz5State(getSparseState(qualifiedVertex), tools::SMALL, tools::PINK, 0);

        double tempDist = (si_->distance(getSparseState(vp), getSparseState(candidateRep)) +
                           si_->distance(getSparseState(candidateRep), getSparseState(qualifiedVertex))) /
          2.0; // do we divide by 2 because of the midpoint path?? TODO(davetcoleman): figure out this proof
        if (tempDist > maxDist)
          maxDist = tempDist;
      }
      BOLT_DEBUG(indent + 6, "Max distance: " << maxDist);

      InterfaceData &iData = getData(candidateRep, vp, vpp, indent + 8);

      if (iData.lastDistance_ == 0)
        BOLT_RED_DEBUG(indent + 6, "last distance is 0");

      // Check if spanner property violated
      if (iData.lastDistance_ > 0 && maxDist > stretchFactor_ * iData.lastDistance_) // DTC added zero check
      {
        BOLT_DEBUG(indent + 6, "Spanner property violated - the max dist was " << maxDist);
        BOLT_DEBUG(indent + 6, "  The stretch factor dist is " << (stretchFactor_ * iData.lastDistance_)
                                                               << " and stretchFactor_ is " << stretchFactor_);

        spannerPropertyWasViolated = true;

        // Can we connect these two vertices directly?
        if (si_->checkMotion(getSparseState(vp), getSparseState(vpp)))
        {
          BOLT_DEBUG(indent + 6, "Adding edge between vp and vpp");

          addEdge(vp, vpp, eRED, indent + 8);
        }
        else
        {
          std::cout << "Path geometric adding path... " << std::endl;
          exit(0);

          geometric::PathGeometric *p = new geometric::PathGeometric(si_);
          if (vp < vpp)
          {
            p->append(iData.interface1Outside_);
            p->append(iData.interface1Inside_);
            p->append(getSparseState(candidateRep));
            p->append(iData.interface2Inside_);
            p->append(iData.interface2Outside_);
          }
          else
          {
            p->append(iData.interface2Outside_);
            p->append(iData.interface2Inside_);
            p->append(getSparseState(candidateRep));
            p->append(iData.interface1Inside_);
            p->append(iData.interface1Outside_);
          }

          BOLT_DEBUG(indent + 6, "Creating path");

          pathSimplifier_->reduceVertices(*p, 10);
          pathSimplifier_->shortcutPath(*p, 50);

          if (p->checkAndRepair(100).second)
          {
            SparseVertex prior = vp;
            SparseVertex vnew;
            std::vector<base::State *> &states = p->getStates();

            BOLT_DEBUG(indent + 6, "check and repair succeeded");

            foreach (base::State *state, states)
            {
              // no need to clone st, since we will destroy p; we just copy the pointer
              BOLT_DEBUG(indent + 6, "Adding node from shortcut path for QUALITY");
              vnew = addVertex(state, QUALITY);

              addEdge(prior, vnew, eGREEN, indent + 8);
              prior = vnew;
            }
            // clear the states, so memory is not freed twice
            states.clear();
            addEdge(prior, vpp, eGREEN, indent + 8);
          }
          else
          {
            BOLT_DEBUG(indent + 6, "check and repair failed?");
            exit(-1);
          }

          delete p;
        }
      }  // end if
    }
  }

  if (!spannerPropertyWasViolated)
  {
    BOLT_DEBUG(indent + 2, "Spanner property was NOT violated, SKIPPING");
  }
  else
  {
      visual_->viz5Trigger();
      usleep(0.001 * 1000000);
      std::cout << "shutting down for viz of checkAddPath " << std::endl;
      exit(0);
  }

  return spannerPropertyWasViolated;
}

SparseVertex SparseDB::findGraphRepresentative(base::State *state, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "findGraphRepresentative()");

  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;

  // Search
  getSparseState(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  getSparseState(queryVertices_[threadID]) = nullptr;

  BOLT_DEBUG(indent + 2, "Found " << graphNeighbors.size() << " nearest neighbors (graph rep) within SparseDelta "
                                  << sparseDelta_);

  SparseVertex result = boost::graph_traits<SparseGraph>::null_vertex();

  for (std::size_t i = 0; i < graphNeighbors.size(); ++i)
  {
    BOLT_DEBUG(indent + 2, "Checking motion of graph representative candidate " << i);
    if (si_->checkMotion(state, getSparseState(graphNeighbors[i])))
    {
      BOLT_DEBUG(indent + 4, "graph representative valid ");
      result = graphNeighbors[i];
      break;
    }
    else
      BOLT_DEBUG(indent + 4, "graph representative NOT valid, checking next ");
  }
  return result;
}

void SparseDB::findCloseRepresentatives(base::State *workState, const base::State *candidateState,
                                        const SparseVertex candidateRep,
                                        std::map<SparseVertex, base::State *> &closeRepresentatives, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "findCloseRepresentatives()");
  BOLT_DEBUG(indent + 2, "nearSamplePoints: " << nearSamplePoints_ << " denseDelta: " << denseDelta_);
  const bool visualizeSampler = false;
  base::State *sampledState = workState; // rename variable just to clarify what it represents temporarily

  assert(closeRepresentatives.empty());

  // Search the space around new potential state candidateState
  for (std::size_t i = 0; i < nearSamplePoints_; ++i)
  {
    BOLT_DEBUG(indent + 2, "Get supporting representative #" << i);

    bool foundValidSample = false;
    static const std::size_t MAX_SAMPLE_ATTEMPT = 1000;
    for (std::size_t attempt = 0; attempt < MAX_SAMPLE_ATTEMPT; ++attempt)
    {
      BOLT_DEBUG(indent + 4, "Sample attempt " << attempt);

      sampler_->sampleNear(sampledState, candidateState, denseDelta_);
      si_->getStateSpace()->setLevel(sampledState, 0);

      if (!si_->isValid(sampledState))
      {
        BOLT_DEBUG(indent + 6, "notValid ");

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);

        continue;
      }
      if (si_->distance(candidateState, sampledState) > denseDelta_)
      {
        BOLT_DEBUG(indent + 6, "Distance too far " << si_->distance(candidateState, sampledState)
                                                   << " needs to be less than " << denseDelta_);

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);
        continue;
      }
      if (!si_->checkMotion(candidateState, sampledState))
      {
        BOLT_DEBUG(indent + 6, "Motion invalid ");

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);
        continue;
      }

      if (visualizeQualityCriteria_ && visualizeSampler)
        visual_->viz3State(sampledState, tools::SMALL, tools::GREEN, 0);

      foundValidSample = true;
      break;
    }  // for each attempt

    visual_->viz3Trigger();
    usleep(0.001 * 1000000);

    if (!foundValidSample)
    {
      BOLT_DEBUG(indent + 4, "Unable to find valid sample after " << MAX_SAMPLE_ATTEMPT << " attempts ");
    }
    else
      BOLT_DEBUG(indent + 4, "Found valid nearby sample");

    // Compute which sparse vertex represents this new candidate vertex
    SparseVertex sampledStateRep = findGraphRepresentative(sampledState, indent + 6);

    // Check if sample is actually seen by somebody (it should be in all likelihood)
    if (sampledStateRep == boost::graph_traits<SparseGraph>::null_vertex())
    {
      BOLT_DEBUG(indent + 4, "Sampled state has no representative (is null) ");

      // This guy can't be seen by anybody, so we should take this opportunity to add him
      BOLT_DEBUG(indent + 4, "Adding node for COVERAGE");
      addVertex(si_->cloneState(sampledState), COVERAGE);

      BOLT_DEBUG(indent + 4, "STOP EFFORTS TO ADD A DENSE PATH");

      // We should also stop our efforts to add a dense path
      for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
           it != closeRepresentatives.end(); ++it)
        si_->freeState(it->second);
      closeRepresentatives.clear();
      break;
    }

    BOLT_DEBUG(indent + 4, "Sampled state has representative (is not null)");

    // If its representative is different than candidateState
    if (sampledStateRep != candidateRep)
    {
      BOLT_DEBUG(indent + 4, "candidateRep != sampledStateRep ");

      // And we haven't already tracked this representative
      if (closeRepresentatives.find(sampledStateRep) == closeRepresentatives.end())
      {
        BOLT_DEBUG(indent + 4, "Track the representative");

        // Track the representative
        closeRepresentatives[sampledStateRep] = si_->cloneState(sampledState);
      }
      else
      {
        BOLT_DEBUG(indent + 4, "Already tracking the representative");
      }
    }
    else
    {
      BOLT_DEBUG(indent + 4, "candidateRep == sampledStateRep, do not keep this sample ");
    }
  }  // for each supporting representative
}

void SparseDB::updatePairPoints(SparseVertex candidateRep, const base::State *candidateState,
                                SparseVertex nearSampledRep, const base::State *nearSampledState, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "updatePairPoints()");

  // First of all, we need to compute all candidate r'
  std::vector<SparseVertex> adjVerticesUnconnected;
  getAdjVerticesOfV1UnconnectedToV2(candidateRep, nearSampledRep, adjVerticesUnconnected, indent + 2);

  // for each pair Pv(r,r')
  foreach (SparseVertex adjVertexUnconnected, adjVerticesUnconnected)
    // Try updating the pair info
    distanceCheck(candidateRep, candidateState, nearSampledRep, nearSampledState, adjVertexUnconnected, indent + 2);
}

void SparseDB::getAdjVerticesOfV1UnconnectedToV2(SparseVertex v1, SparseVertex v2,
                                                 std::vector<SparseVertex> &adjVerticesUnconnected, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "getAdjVerticesOfV1UnconnectedToV2()");

  adjVerticesUnconnected.clear();
  foreach (SparseVertex adjVertex, boost::adjacent_vertices(v1, g_))
    if (adjVertex != v2)
      if (!hasEdge(adjVertex, v2))
        adjVerticesUnconnected.push_back(adjVertex);

  BOLT_DEBUG(indent + 2, "adjVerticesUnconnected size: " << adjVerticesUnconnected.size());
}

void SparseDB::getMaxSpannerPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::vector<SparseVertex> &qualifiedVertices,
                                 std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "getMaxSpannerPath()");

  qualifiedVertices.clear();

  foreach (SparseVertex adjVertex, boost::adjacent_vertices(vpp, g_))
  {
    if (hasEdge(adjVertex, v) && !hasEdge(adjVertex, vp))
    {
      InterfaceData &iData = getData(v, vpp, adjVertex, indent + 4);

      // DTC: Check if we previously had found a pair of points that support this interface
      if ((vpp < adjVertex && iData.interface1Inside_) || (adjVertex < vpp && iData.interface2Inside_))
      {
        BOLT_GREEN_DEBUG(indent + 2, "Found an additional qualified vertex!");
        // We have, so we need to check if we've found a better pair of points
        qualifiedVertices.push_back(adjVertex);
      }
    }
  }

  // vpp is always qualified because of its previous checks
  qualifiedVertices.push_back(vpp);

  BOLT_DEBUG(indent + 2, "Total qualified vertices found: " << qualifiedVertices.size());
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

InterfaceData &SparseDB::getData(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "getData() " << v << ", " << vp << ", " << vpp);
  return interfaceDataProperty_[v].interfaceHash[index(vp, vpp)];
}

void SparseDB::distanceCheck(SparseVertex v, const base::State *q, SparseVertex vp,
                             const base::State *qp, SparseVertex vpp, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "distanceCheck()");

  // Get the info for the current representative-neighbors pair
  InterfaceData &iData = getData(v, vp, vpp, indent + 4);

  if (vp < vpp)  // FIRST points represent r (the guy discovered through sampling)
  {
    if (iData.interface1Inside_ == nullptr)  // If the point we're considering replacing (P_v(r,.)) isn't there
    { // we know we're doing better, so add it
      BOLT_DEBUG(indent + 2, "setFirst");
      iData.setFirst(q, qp, si_);
    }
    else  // Otherwise, it is there,
    {
      if (iData.interface2Inside_ == nullptr)  // But if the other guy doesn't exist, we can't compare.
      {
        // Should probably keep the one that is further away from rep?  Not known what to do in this case.
        // TODO: is this not part of the algorithm?
        BOLT_RED_DEBUG(indent + 2, "TODO no pointB");
      }
      // We know both of these points exist, so we can check some distances
      // TODO(davetcoleman): why does it not use lastDistance_ here???
      else if (si_->distance(q, iData.interface2Inside_) < si_->distance(iData.interface1Inside_, iData.interface2Inside_))
      { // Distance with the new point is good, so set it.
        BOLT_GREEN_DEBUG(indent + 2, "setFirst UPDATED");
        iData.setFirst(q, qp, si_);
      }
      else
        BOLT_DEBUG(indent + 2, "Distance was not better, not updating bookkeeping");
    }
  }
  else  // SECOND points represent r (the guy discovered through sampling)
  {
    if (iData.interface2Inside_ == nullptr)  // If the point we're considering replacing (P_V(.,r)) isn't there...
    {
      BOLT_DEBUG(indent + 2, "setSecond");
      // we must be doing better, so add it
      iData.setSecond(q, qp, si_);
    }
    else  // Otherwise, it is there
    {
      if (iData.interface1Inside_ == nullptr)  // But if the other guy doesn't exist, we can't compare.
      {
        // Should we be doing something cool here?
        BOLT_RED_DEBUG(indent + 2, "TODO");
      }
      else if (si_->distance(q, iData.interface1Inside_) < si_->distance(iData.interface2Inside_, iData.interface1Inside_))
      { // Distance with the new point is good, so set it
        BOLT_GREEN_DEBUG(indent + 2, "setSecond UPDATED");
        iData.setSecond(q, qp, si_);
      }
    }
  }

  // Lastly, save what we have discovered
  // TODO(davetcoleman): do we really need to copy this back in or is it already passed by reference?
  interfaceDataProperty_[v].interfaceHash[index(vp, vpp)] = iData;
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
                                  std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, "findGraphNeighbors()");
  const bool verbose = false;

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
        visual_->viz3DeleteAllMarkers();
        visual_->viz3State(state, tools::SMALL, tools::BLUE, 0);
        visual_->viz3State(getSparseState(graphNeighborhood[i]), tools::SMALL, tools::BLUE, 0);
        visual_->viz3Trigger();
        usleep(1 * 1000000);
      }
    }
    else if (verbose)
      std::cout << " ---- Skipping collision checking because same vertex " << std::endl;

    // The two are visible to each other!
    visibleNeighborhood.push_back(graphNeighborhood[i]);
  }

  BOLT_DEBUG(indent + 2, "Graph neighborhood: " << graphNeighborhood.size()
                                                << " | Visible neighborhood: " << visibleNeighborhood.size());
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
    tools::colors color;
    if (type == COVERAGE)
      color = tools::BLACK;
    else if (type == CONNECTIVITY)
      color = tools::ORANGE;
    else if (type == INTERFACE)
      color = tools::RED;
    else if (type == QUALITY)
      color = tools::BLUE;
    else
      OMPL_ERROR("Unknown mode");

    if (type == COVERAGE)
      visual_->viz2State(getSparseState(v), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

    if (type == QUALITY)
      visual_->viz2State(getSparseState(v), tools::LARGE, color, 0);
    else
      visual_->viz2State(getSparseState(v), tools::MEDIUM, color, 0);

    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    if (visualizeOverlayNodes_)  // after initial spars graph is created, show additions not from grid
    {
      visual_->viz3State(getSparseState(v), tools::MEDIUM, color, sparseDelta_);
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

void SparseDB::addEdge(SparseVertex v1, SparseVertex v2, std::size_t visualColor, std::size_t indent)
{
  assert(v1 <= getNumVertices());
  assert(v2 <= getNumVertices());
  assert(v1 != v2);
  BOOST_ASSERT_MSG(getSparseState(v1) != getSparseState(v2), "States on both sides of an edge are the same");

  // std::cout << std::string(indent, "addEdge: Connecting vertex " << v1 << " to vertex " << v2 <<
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
  visual_->viz2DeleteAllMarkers();

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
      visual_->viz2Edge(getSparseStateConst(v1), getSparseStateConst(v2),
                        tools::eRED);  // edgeWeightPropertySparse_[e]);

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
        visual_->viz2State(getSparseStateConst(v), tools::SMALL, tools::BLUE, 1);
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

bool SparseDB::hasEdge(SparseVertex v1, SparseVertex v2)
{
  return boost::edge(v1, v2, g_).second;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
