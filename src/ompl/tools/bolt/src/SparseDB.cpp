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

// SparseEdgeWeightMap methods ////////////////////////////////////////////////////////////////////////////

BOOST_CONCEPT_ASSERT((boost::ReadablePropertyMapConcept<otb::SparseEdgeWeightMap, otb::SparseEdge>));

namespace boost
{
double get(const otb::SparseEdgeWeightMap &m, const otb::SparseEdge &e)
{
  return m.get(e);
}
}

// CustomAstarVisitor methods ////////////////////////////////////////////////////////////////////////////

BOOST_CONCEPT_ASSERT((boost::AStarVisitorConcept<otb::CustomAstarVisitor, otb::SparseGraph>));

otb::CustomAstarVisitor::CustomAstarVisitor(SparseVertex goal, SparseDB *parent) : goal_(goal), parent_(parent)
{
}

void otb::CustomAstarVisitor::discover_vertex(SparseVertex v, const SparseGraph &) const
{
  // Statistics
  parent_->recordNodeOpened();

  if (parent_->visualizeAstar_)
    parent_->getVisual()->viz4State(parent_->getSparseState(v), tools::SMALL, tools::GREEN, 1);
}

void otb::CustomAstarVisitor::examine_vertex(SparseVertex v, const SparseGraph &) const
{
  // Statistics
  parent_->recordNodeClosed();

  if (parent_->visualizeAstar_)
  {
    parent_->getVisual()->viz4State(parent_->getSparseState(v), tools::LARGE, tools::BLACK, 1);
    parent_->getVisual()->viz4Trigger();
    usleep(parent_->visualizeAstarSpeed_ * 1000000);
  }

  if (v == goal_)
    throw FoundGoalException();
}

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
{
  // Add search state
  initializeQueryState();

  // Initialize nearest neighbor datastructure
  nn_.reset(new NearestNeighborsGNAT<SparseVertex>());
  nn_->setDistanceFunction(boost::bind(&otb::SparseDB::distanceFunction, this, _1, _2));

  // Initialize discretizer
  vertexDiscretizer_.reset(new VertexDiscretizer(si_, visual_));
}

SparseDB::~SparseDB(void)
{
  freeMemory();
}

void SparseDB::freeMemory()
{
  sampler_.reset();

  foreach (SparseVertex v, boost::vertices(g_))
  {
    foreach (InterfaceData &iData, interfaceDataProperty_[v] | boost::adaptors::map_values)
      iData.clear(si_);
  }

  g_.clear();

  if (nn_)
    nn_->clear();
}

bool SparseDB::setup()
{
  const std::size_t indent = 0;

  // Dimensions / joints
  std::size_t dim = si_->getStateDimension();

  // Max distance across configuration space
  maxExtent_ = si_->getMaximumExtent();

  // Vertex visibility region size
  sparseDelta_ = sparseDeltaFraction_ * maxExtent_;

  // Sampling for interfaces visibility size
  denseDelta_ = denseDeltaFraction_ * maxExtent_;

  // Number of points to test for interfaces around a sample for the quality criterion
  nearSamplePoints_ = nearSamplePointsMultiple_ * si_->getStateDimension();

  // Discretization for initial input into sparse graph
  //double sparseMultiple = si_->getStateDimension() * magicMultiple_;
  double sparseMultiple = 1.25;

  if (false) // Method 1
  {
    discretization_ = sparseDelta_ * sparseMultiple;
  }
  else // Method 2
  {
    penetrationDistance_ = magicMultiple_; //0.01; // TODO(davetcoleman): how to choose this
    const double discFactor = sparseDelta_ - penetrationDistance_;
    discretization_ = 2 * sqrt( std::pow(discFactor, 2) / dim);
  }

  ignoreEdgesSmallerThan_ = discretization_ + 0.01;

  // Calculate optimum stretch factor
  //stretchFactor_ = discretization_ / (0.5 * discretization_ * sqrt(2) - 2.0 * denseDelta_);
  nearestDVertex_ = sqrt(dim * std::pow(0.5 * discretization_, 2)); // z in my calculations
  //stretchFactor_ = discretization_ / ( nearestDVertex_ - 2.0 * denseDelta_); // 2D case but not 3D
  //stretchFactor_ = (discretization_ + nearestDVertex_) / ( 2 * (nearestDVertex_ - 2.0 * denseDelta_)); // N-D case
  //stretchFactor_ = discretization_ / (discretization_ - 2.0 * denseDelta_); // N-D case
  //stretchFactor_ = 1.83;

  BOLT_DEBUG(indent, 1, "--------------------------------------------------");
  BOLT_DEBUG(indent, 1, "Sparse DB Setup:");
  BOLT_DEBUG(indent + 2, 1, "Max Extent              = " << maxExtent_);
  BOLT_DEBUG(indent + 2, 1, "Sparse Delta            = " << sparseDelta_);
  BOLT_DEBUG(indent + 2, 1, "Dense Delta             = " << denseDelta_);
  BOLT_DEBUG(indent + 2, 1, "State Dimension         = " << dim);
  BOLT_DEBUG(indent + 2, 1, "Near Sample Points      = " << nearSamplePoints_);
  BOLT_DEBUG(indent + 2, 1, "Sparse Multiple         = " << sparseMultiple);
  BOLT_DEBUG(indent + 2, 1, "Discretization          = " << discretization_);
  BOLT_DEBUG(indent + 2, 1, "Nearest Discretized V   = " << nearestDVertex_);
  BOLT_DEBUG(indent + 2, 1, "Stretch Factor          = " << stretchFactor_);
  BOLT_DEBUG(indent + 2, 1, "Viz ignore edges below  = " << ignoreEdgesSmallerThan_);
  BOLT_DEBUG(indent, 1, "--------------------------------------------------");

  assert(maxExtent_ > 0);
  assert(denseDelta_ > 0);
  assert(nearSamplePoints_ > 0);
  assert(sparseDelta_ > 0);
  assert(sparseDelta_ > 0.000000001);  // Sanity check

  // Load state sampler
  if (!sampler_)
  {
    sampler_ = ob::MinimumClearanceValidStateSamplerPtr(new ob::MinimumClearanceValidStateSampler(si_.get()));
    sampler_->setMinimumObstacleClearance(obstacleClearance_);
    si_->getStateValidityChecker()->setClearanceSearchDistance(obstacleClearance_);
  }

  if (si_->getStateValidityChecker()->getClearanceSearchDistance() < obstacleClearance_)
    OMPL_WARN("State validity checker clearance search distance %f is less than the required obstacle clearance %f for our "
              "state sampler, incompatible settings!", si_->getStateValidityChecker()->getClearanceSearchDistance(), obstacleClearance_);

  // Initialize path simplifier
  if (!pathSimplifier_)
  {
    pathSimplifier_.reset(new geometric::PathSimplifier(si_));
    pathSimplifier_->freeStates(false);
  }

  // Configure vertex discretizer
  vertexDiscretizer_->setMinimumObstacleClearance(obstacleClearance_);
  vertexDiscretizer_->setDiscretization(discretization_);

  return true;
}

bool SparseDB::astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex> &vertexPath,
                           double &distance, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "astarSearch()");
  indent += 2;

  // Hold a list of the shortest path parent to each vertex
  SparseVertex *vertexPredecessors = new SparseVertex[getNumVertices()];
  // boost::vector_property_map<SparseVertex> vertexPredecessors(getNumVertices());

  bool foundGoal = false;
  double *vertexDistances = new double[getNumVertices()];

  // Reset statistics
  numNodesOpened_ = 0;
  numNodesClosed_ = 0;

  if (visualizeAstar_)
  {
    visual_->viz4DeleteAllMarkers();
  }

  try
  {
    double popularityBias = 0;
    bool popularityBiasEnabled = false;
    boost::astar_search(g_,                                                           // graph
                        start,                                                        // start state
                        boost::bind(&otb::SparseDB::astarHeuristic, this, _1, goal),  // the heuristic
                        // ability to disable edges (set cost to inifinity):
                        boost::weight_map(SparseEdgeWeightMap(g_, edgeCollisionStatePropertySparse_, popularityBias,
                                                              popularityBiasEnabled))
                            .predecessor_map(vertexPredecessors)
                            .distance_map(&vertexDistances[0])
                            .visitor(CustomAstarVisitor(goal, this)));
  }
  catch (FoundGoalException &)
  {
    distance = vertexDistances[goal];

    // the custom exception from CustomAstarVisitor
    BOLT_DEBUG(indent, vCriteria_, "AStar found solution. Distance to goal: " << vertexDistances[goal]);
    BOLT_DEBUG(indent, vCriteria_, "Number nodes opened: " << numNodesOpened_
                                                           << ", Number nodes closed: " << numNodesClosed_);

    if (isinf(vertexDistances[goal]))  // TODO(davetcoleman): test that this works
    {
      BOLT_RED_DEBUG(indent, true, "Distance to goal is infinity");
      exit(-1);
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

      // Add the start state to the path, unless this path is just one vertex long and the start==goal
      if (v != goal)  // TODO explain this because i don't understand
      {
        vertexPath.push_back(v);
      }
      else  // TODO remove this
      {
        std::cout << "vertex path is one vertex long? " << std::endl;
        std::cout << "this should be deleted " << std::endl;
        exit(-1);
      }

      foundGoal = true;
    }
  }

  if (!foundGoal)
    BOLT_YELLOW_DEBUG(indent, vCriteria_, "Did not find goal");

  // Show all predecessors
  if (visualizeAstar_)
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "Show all predecessors");
    for (std::size_t i = 1; i < getNumVertices(); ++i)  // skip vertex 0 b/c that is the search vertex
    {
      const SparseVertex v1 = i;
      const SparseVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4Edge(getSparseState(v1), getSparseState(v2), 10);
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

double SparseDB::astarHeuristic(const SparseVertex a, const SparseVertex b) const
{
  // Assume vertex 'a' is the one we care about its populariy

  // Get the classic distance
  double dist = si_->distance(getSparseState(a), getSparseState(b));

  if (false)  // method 1
  {
    const double percentMaxExtent = (maxExtent_ * percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
    double popularityComponent = percentMaxExtent * (vertexPopularity_[a] / 100.0);

    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << vertexPopularity_[a]
              << ", max extent: " << maxExtent_ << ", percentMaxExtent: " << percentMaxExtent
              << ", popularityComponent: " << popularityComponent;
    dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 2
  {
    const double percentDist = (dist * percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
    double popularityComponent = percentDist * (vertexPopularity_[a] / 100.0);

    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << vertexPopularity_[a]
              << ", percentDist: " << percentDist << ", popularityComponent: " << popularityComponent;
    dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 3
  {
    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << vertexPopularity_[a]
              << ", vertexPopularity_[a] / 100.0: " << vertexPopularity_[a] / 100.0
              << ", percentMaxExtentUnderestimate_: " << percentMaxExtentUnderestimate_;
    // if ((vertexPopularity_[a] / 100.0) < (1 - percentMaxExtentUnderestimate_))
    if (vertexPopularity_[a] > (100 - percentMaxExtentUnderestimate_ * 100.0))
    {
      dist = 0;
    }

    // dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 4
  {
    dist *= (1 + percentMaxExtentUnderestimate_);
  }
  // method 5: increasing the sparseDelta fraction

  // std::cout << ", new distance: " << dist << std::endl;

  return dist;
}

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
  std::size_t indent = 2;

  // Benchmark runtime
  time::point startTime = time::now();

  numSamplesAddedForQuality_ = 0;
  numSamplesAddedForConnectivity_ = 0;
  numSamplesAddedForInterface_ = 0;
  numSamplesAddedForQuality_ = 0;

  numConsecutiveFailures_ = 0;
  useFourthCriteria_ = false;  // initially we do not do this step

  // Create SPARS
  CALLGRIND_TOGGLE_COLLECT;

  // Start the graph off with discretized states
  if (useDiscretizedSamples_)
  {
    addDiscretizedStates(indent);
  }

  // Finish the graph with random samples
  if (useRandomSamples_)
  {
    addRandomSamples(indent);
  }
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_DUMP_STATS;

  // Benchmark runtime
  double duration = time::seconds(time::now() - startTime);

  // Statistics
  numGraphGenerations_++;

  // Check how many connected components exist, possibly throw error
  std::size_t numSets = getDisjointSetsCount();
  std::pair<std::size_t, std::size_t> interfaceStats = getInterfaceStateStorageSize();


  BOLT_DEBUG(0, 1, "-----------------------------------------");
  BOLT_DEBUG(0, 1, "Created SPARS graph                      ");
  BOLT_DEBUG(0, 1, "  Vertices:                  " << getNumVertices());
  BOLT_DEBUG(0, 1, "  Edges:                     " << getNumEdges());
  BOLT_DEBUG(0, 1, "  Generation time:           " <<  duration);
  BOLT_DEBUG(0, 1, "  Total generations:         " << numGraphGenerations_);
  BOLT_DEBUG(0, 1, "  Disjoint sets:             " << numSets);
  BOLT_DEBUG(0, 1, "  Edge collision cache         ");
  BOLT_DEBUG(0, 1, "    Size:                    " << edgeCache_->getCacheSize());
  BOLT_DEBUG(0, 1, "    Total checks:            " << edgeCache_->getTotalCollisionChecks());
  BOLT_DEBUG(0, 1, "    Cached checks:           " << edgeCache_->getTotalCollisionChecksFromCache() << " ("
             << edgeCache_->getPercentCachedCollisionChecks() << "%)");
  BOLT_DEBUG(0, 1, "  Criteria additions:          ");
  BOLT_DEBUG(0, 1, "    Coverage:                " << numSamplesAddedForCoverage_);
  BOLT_DEBUG(0, 1, "    Connectivity:            " << numSamplesAddedForConnectivity_);
  BOLT_DEBUG(0, 1, "    Interface:               " << numSamplesAddedForInterface_);
  BOLT_DEBUG(0, 1, "    Quality:                 " << numSamplesAddedForQuality_);
  BOLT_DEBUG(0, 1, "  Num random samples added:  " << numRandSamplesAdded_);
  BOLT_DEBUG(0, 1, "  InterfaceData:                       ");
  BOLT_DEBUG(0, 1, "    States stored:           " << interfaceStats.first);
  BOLT_DEBUG(0, 1, "    Missing interfaces:      " << interfaceStats.second);
  BOLT_DEBUG(0, 1, "-----------------------------------------");

  if (!visualizeSparsGraph_)
    displaySparseDatabase();

  // Save collision cache
  edgeCache_->save();

  OMPL_INFORM("Finished creating sparse database");
}

void SparseDB::addDiscretizedStates(std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, true, "addDiscretizedStates()");

  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID];
  std::cout << "range: " << range << std::endl;
  const std::size_t jointIncrements = floor(range / vertexDiscretizer_->getDiscretization());
  std::cout << "jointIncrements: " << jointIncrements << std::endl;
  double leftOver = range - jointIncrements * vertexDiscretizer_->getDiscretization();
  std::cout << "leftOver: " << leftOver << std::endl;
  double startOffset = leftOver / 2;
  std::cout << "startOffset: " << startOffset << std::endl;

  //  double percentVisibilityRegionEdgeOverlap = 0.5;
  //double startOffset = sparseDelta_ * percentVisibilityRegionEdgeOverlap;

  // Create two levels of grids
  for (std::size_t i = 0; i < 2; ++i)
  {
    OMPL_INFORM("Generating grid iteration %u", i);

    // Set starting value offset
    if (i == 0)
      vertexDiscretizer_->setStartingValueOffset(startOffset);
    else
      vertexDiscretizer_->setStartingValueOffset(startOffset + vertexDiscretizer_->getDiscretization() / 2.0);

    // Generate verticies
    vertexDiscretizer_->generate();

    // Convert to proper format TODO(davetcoleman): remove the need for this format?
    std::vector<base::State *> &candidateVertices = vertexDiscretizer_->getCandidateVertices();

    std::list<WeightedVertex> vertexInsertionOrder;
    for (base::State *state : candidateVertices)
    {
      // DenseDB now 'owns' the memory of candidateVertices and will unload it later
      DenseVertex v = denseDB_->addVertex(state, COVERAGE);
      vertexInsertionOrder.push_back(WeightedVertex(v, 0));
    }
    candidateVertices.clear();  // clear the vector because we've moved all its memory pointers to DenseDB
    std::size_t sucessfulInsertions;
    createSPARSInnerLoop(vertexInsertionOrder, sucessfulInsertions);
    std::cout << "sucessfulInsertions: " << sucessfulInsertions << std::endl;
  }
}

void SparseDB::createSPARSOuterLoop()
{
  std::size_t indent = 2;

  /*
  // Clear the old spars graph
  if (getNumVertices() > queryVertices_.size())
  {
    OMPL_INFORM("Resetting sparse database");
    freeMemory();
    initializeQueryState();  // Re-add search state

    if (visualizeSparsGraph_)  // Clear visuals
      visual_->viz2DeleteAllMarkers();
  }
  */

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
    std::cout << "Attempting to insert " << vertexInsertionOrder.size() << " vertices for the " << loopAttempt
              << " loop" << std::endl;

    // Sanity check
    if (loopAttempt > 3)
      OMPL_WARN("Suprising number of loop when attempting to insert nodes into SPARS graph: %u", loopAttempt);

    // Benchmark runtime
    time::point startTime = time::now();

    // ----------------------------------------------------------------------
    // Attempt to insert each vertex using the first 3 criteria
    if (!createSPARSInnerLoop(vertexInsertionOrder, sucessfulInsertions))
      break;

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

      // Save collision cache, just in case there is a bug
      edgeCache_->save();
    }

    bool debugOverRideJustTwice = true;
    if (debugOverRideJustTwice && loopAttempt == 1)
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

void SparseDB::addRandomSamples(std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, true, "addRandomSamples()");
  indent += 2;

  // For each dense vertex we add
  numRandSamplesAdded_ = 0;
  std::size_t addedStatesCount = 0;  // count how many states we add
  while (numRandSamplesAdded_ < 10000)
  {
    // Add dense vertex
    base::State *state = si_->allocState();
    // TODO(davetcoleman): remove dependence on denseDB vertices
    DenseVertex denseV = denseDB_->addVertex(state, COVERAGE);  // TODO(davetcoleman): COVERAGE means nothing

    // For each random sample
    while (true)
    {
      // Sample randomly
      if (!sampler_->sample(state))  // TODO(davetcoleman): is it ok with the DenseDB.nn_ to change the state after
      // having added it to the nearest neighbor?? No, I don't think it is.
      {
        OMPL_ERROR("Unable to find valid sample");
        exit(-1);  // this should never happen
      }
      //si_->getStateSpace()->setLevel(state, 0);  // TODO no hardcode

      // Run SPARS checks
      GuardType addReason;     // returns why the state was added
      SparseVertex newVertex;  // the newly generated sparse vertex
      const std::size_t threadID = 0;
      if (addStateToRoadmap(denseV, newVertex, addReason, threadID))
      {
        if (addedStatesCount % 10 == 0)
          BOLT_DEBUG(indent, true, "Added random sample, total new states: " << ++addedStatesCount);

        // Statistics
        numRandSamplesAdded_++;

        break;  // must break so that a new state can be allocated
      }
      else if (numConsecutiveFailures_ % 500 == 0)
      {
        BOLT_DEBUG(indent, true, "Random sampled failed, consecutive failures: " << numConsecutiveFailures_);
      }

      if (numConsecutiveFailures_ >= fourthCriteriaAfterFailures_ && !useFourthCriteria_)
      {
        BOLT_DEBUG(indent, true, "Starting to check for 4th quality criteria because "
                                     << numConsecutiveFailures_ << " consecutive failures have occured");
        useFourthCriteria_ = true;
        visualizeOverlayNodes_ = true;  // visualize all added nodes in a separate window
        numConsecutiveFailures_ = 0;    // reset for new criteria

        // Show it just once if it has not already been animated
        if (!visualizeVoronoiDiagramAnimated_ && visualizeVoronoiDiagram_)
          visual_->vizVoronoiDiagram();
      }

      if (useFourthCriteria_ && numConsecutiveFailures_ > terminateAfterFailures_)
      {
        OMPL_INFORM("SPARS creation finished because %u consecutive insertion failures reached",
                    terminateAfterFailures_);
        return;
      }
    }  // end while
  }    // end while
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
  BOLT_DEBUG(indent, vCriteria_, "addStateToRoadmap() Adding DenseVertex " << denseV);

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
    BOLT_DEBUG(indent + 2, vCriteria_, "State added for: COVERAGE");

    addReason = COVERAGE;
    stateAdded = true;
  }
  else if (checkAddConnectivity(denseV, visibleNeighborhood, newVertex, indent + 6))
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "State added for: CONNECTIVITY");

    addReason = CONNECTIVITY;
    stateAdded = true;
  }
  else if (checkAddInterface(denseV, graphNeighborhood, visibleNeighborhood, newVertex, indent + 10))
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "State added for: INTERFACE");

    addReason = INTERFACE;
    stateAdded = true;
  }
  else if (useFourthCriteria_ &&
           checkAddQuality(denseV, graphNeighborhood, visibleNeighborhood, workState, newVertex, indent + 14))
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "State added for: 4th CRITERIA");

    addReason = QUALITY;
    stateAdded = true;
  }
  else
  {
    if (vCriteria_)
      std::cout << "Did NOT add state for any criteria " << std::endl;
    numConsecutiveFailures_++;
  }

  if (stateAdded)
    numConsecutiveFailures_ = 0;

  si_->freeState(workState);

  return stateAdded;
}

bool SparseDB::checkAddCoverage(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "checkAddCoverage() Are other nodes around it visible?");

  // Only add a node for coverage if it has no neighbors
  if (visibleNeighborhood.size() > 0)
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "NOT adding node for coverage ");
    return false;  // has visible neighbors
  }

  // No free paths means we add for coverage
  BOLT_DEBUG(indent + 2, vCriteria_, "Adding node for COVERAGE ");

  newVertex = addVertex(denseV, COVERAGE);

  // Note: we do not connect this node with any edges because we have already determined
  // it is too far away from any nearby nodes

  return true;
}

bool SparseDB::checkAddConnectivity(DenseVertex denseV, std::vector<SparseVertex> &visibleNeighborhood,
                                    SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "checkAddConnectivity() Does this node connect two disconnected components?");

  // If less than 2 neighbors there is no way to find a pair of nodes in different connected components
  if (visibleNeighborhood.size() < 2)
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "NOT adding node for connectivity");
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
    BOLT_DEBUG(indent + 2, vCriteria_, "NOT adding node for connectivity");
    return false;
  }

  BOLT_DEBUG(indent + 2, vCriteria_, "Adding node for CONNECTIVITY ");

  // Add the node
  newVertex = addVertex(denseV, CONNECTIVITY);

  // Add the edges
  for (std::set<SparseVertex>::const_iterator vertexIt = statesInDiffConnectedComponents.begin();
       vertexIt != statesInDiffConnectedComponents.end(); ++vertexIt)
  {
    // Do not add edge from self to self
    if (si_->getStateSpace()->equalStates(getSparseState(*vertexIt), getSparseState(newVertex)))
    {
      std::cout << "Prevented same vertex from being added twice " << std::endl;
      continue;  // skip this pairing
    }

    BOLT_DEBUG(indent + 3, vCriteria_, "Loop: Adding vertex " << *vertexIt);

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
  BOLT_BLUE_DEBUG(indent, vCriteria_, "checkAddInterface() Does this node's neighbor's need it to better connect "
                                      "them?");

  // If there are less than two neighbors the interface property is not applicable, because requires
  // two closest visible neighbots
  if (visibleNeighborhood.size() < 2)
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "NOT adding node for interface (less than 2 visible neighbors)");
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
      {
        BOLT_DEBUG(indent + 2, vCriteria_, "INTERFACE: directly connected nodes");

        // Connect them
        addEdge(visibleNeighborhood[0], visibleNeighborhood[1], visualColor, indent + 4);
      }
      else  // They cannot be directly connected
      {
        // Add the new node to the graph, to bridge the interface
        BOLT_DEBUG(indent + 2, vCriteria_, "Adding node for INTERFACE");

        newVertex = addVertex(denseV, INTERFACE);
        addEdge(newVertex, visibleNeighborhood[0], visualColor, indent + 4);
        addEdge(newVertex, visibleNeighborhood[1], visualColor, indent + 4);
        BOLT_DEBUG(indent + 2, vCriteria_, "INTERFACE: connected two neighbors through new interface node");
      }

      // Report success
      return true;
    }
    else
    {
      BOLT_DEBUG(indent + 2, vCriteria_, "Two closest two neighbors already share an edge, not connecting them");
    }
  }
  BOLT_DEBUG(indent + 2, vCriteria_, "NOT adding node for interface");
  return false;
}

bool SparseDB::checkAddQuality(DenseVertex denseV, std::vector<SparseVertex> &graphNeighborhood,
                               std::vector<SparseVertex> &visibleNeighborhood, base::State *workState,
                               SparseVertex &newVertex, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "checkAddQuality() Ensure SPARS asymptotic optimality");
  indent += 2;

  if (visibleNeighborhood.empty())
  {
    BOLT_DEBUG(indent, vQuality_, "no visible neighbors, not adding 4th criteria ");
    return false;
  }

  base::State *candidateState = getDenseState(denseV);  // paper's name: q
  SparseVertex candidateRep = visibleNeighborhood[0];   // paper's name: v

  bool added = false;
  std::map<SparseVertex, base::State *> closeRepresentatives;  // [nearSampledRep, nearSampledState]
  findCloseRepresentatives(workState, candidateState, candidateRep, closeRepresentatives, indent + 2);

  if (closeRepresentatives.size())
  {
    BOLT_GREEN_DEBUG(indent, vQuality_, "back in checkAddQuality(): Found " << closeRepresentatives.size()
                                                                            << " close "
                                                                               "representatives");
  }
  else
  {
    BOLT_RED_DEBUG(indent, vQuality_, "back in checkAddQuality(): Found " << closeRepresentatives.size()
                                                                          << " close "
                                                                             "representatives");
  }

  bool updated = false;  // track whether a change was made to any vertices' representatives

  // For each pair of close representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    BOLT_DEBUG(indent + 2, vQuality_, "Looping through close representatives");
    base::State *nearSampledState = it->second;  // paper: q'
    SparseVertex nearSampledRep = it->first;     // paper: v'

    if (visualizeQualityCriteria_ && false)  // Visualization
    {
      visual_->viz3Edge(getSparseState(nearSampledRep), nearSampledState, tools::eRED);

      visual_->viz3State(nearSampledState, tools::MEDIUM, tools::GREEN, 0);

      // Replicate a regular vertex visualization
      //visual_->viz3State(getSparseState(nearSampledRep), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);
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
    if (updatePairPoints(candidateRep, candidateState, nearSampledRep, nearSampledState, indent + 2))
      updated = true;

    // ALSO attempt to update bookkeeping for neighboring node nearSampleRep (v')
    if (updatePairPoints(nearSampledRep, nearSampledState, candidateRep, candidateState, indent + 2))
      updated = true;
  }

  BOLT_DEBUG(indent, vQuality_, "Done updating pair points");

  if (!updated)
  {
    BOLT_DEBUG(indent, vQuality_, "No representatives were updated, so not calling checkAddPath()");
    added = false;
    return added;
  }

  // Visualize the interfaces around the candidate rep
  if (visualizeQualityCriteria_)
  {
    //visualizeCheckAddQuality(candidateState, candidateRep);
    //visualizeInterfaces(candidateRep, indent + 2);

    // static std::size_t updateCount = 0;
    // if (updateCount++ % 50 == 0)
    //   visualizeAllInterfaces(indent + 2);
  }

  // Attempt to find shortest path through closest neighbour
  if (checkAddPath(candidateRep, indent + 2))
  {
    BOLT_DEBUG(indent, vQuality_, "nearest visible neighbor added for path ");
    added = true;
  }

  // Attempt to find shortest path through other pairs of representatives
  for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
       it != closeRepresentatives.end(); ++it)
  {
    BOLT_YELLOW_DEBUG(indent, vQuality_, "Looping through close representatives to add path ===============");

    base::State *nearSampledState = it->second;  // paper: q'
    SparseVertex nearSampledRep = it->first;     // paper: v'
    if (checkAddPath(nearSampledRep, indent + 2))
    {
      BOLT_DEBUG(indent, vQuality_, "Close representative added for path");
      added = true;
    }

    // Delete state that was allocated and sampled within this function
    si_->freeState(nearSampledState);
  }

  return added;
}

void SparseDB::visualizeCheckAddQuality(base::State *candidateState, SparseVertex candidateRep)
{
  visual_->viz3DeleteAllMarkers();

  visual_->viz3Edge(candidateState, getSparseState(candidateRep), tools::eORANGE);

  // Show candidate state
  //visual_->viz3State(candidateState, tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, denseDelta_);
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

bool SparseDB::checkAddPath(SparseVertex v, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, vQuality_, "checkAddPath() v = " << v);
  indent += 2;
  bool spannerPropertyWasViolated = false;

  // Candidate v" vertices as described in the method, filled by function getAdjVerticesOfV1UnconnectedToV2().
  std::vector<SparseVertex> adjVerticesUnconnected;

  // Copy adjacent vertices into vector because we might add additional edges during this function
  std::vector<SparseVertex> adjVertices;
  foreach (SparseVertex adjVertex, boost::adjacent_vertices(v, g_))
    adjVertices.push_back(adjVertex);

  BOLT_DEBUG(indent, vQuality_, "Vertex v = " << v << " has " << adjVertices.size() << " adjacent vertices, "
                                                                                       "looping:");

  // Loop through adjacent vertices
  for (std::size_t i = 0; i < adjVertices.size() && !spannerPropertyWasViolated; ++i)
  {
    SparseVertex vp = adjVertices[i];  // vp = v' from paper

    BOLT_DEBUG(indent + 2, vQuality_, "Checking v' = " << vp << " ----------------------------");

    // Compute all nodes which qualify as a candidate v" for v and vp
    getAdjVerticesOfV1UnconnectedToV2(v, vp, adjVerticesUnconnected, indent + 4);

    // for each vertex v'' that is adjacent to v (has a valid edge) and does not share an edge with v'
    foreach (SparseVertex vpp, adjVerticesUnconnected)  // vpp = v'' from paper
    {
      BOLT_DEBUG(indent + 4, vQuality_, "Checking v'' = " << vpp);

      InterfaceData &iData = getData(v, vp, vpp, indent + 6);

      // Check if we need to actually add path
      //if (spannerTestOriginal(v, vp, vpp, iData, indent + 2))
        //if (spannerTestOuter(v, vp, vpp, iData, indent + 2))
      if (spannerTestAStar(v, vp, vpp, iData, indent + 2))
      {
        // Actually add the vertices and possibly edges
        addQualityPath(v, vp, vpp, iData, indent + 6);

        spannerPropertyWasViolated = true;
      }

    }  // foreach vpp
  }    // foreach vp

  if (!spannerPropertyWasViolated)
  {
    BOLT_DEBUG(indent, vQuality_, "Spanner property was NOT violated, SKIPPING");
  }

  return spannerPropertyWasViolated;
}

void SparseDB::visualizeCheckAddPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData)
{
  visual_->viz5DeleteAllMarkers();

  // Show candidate rep
  visual_->viz5State(getSparseState(v), tools::LARGE, tools::BLUE, 0);

  // Show adjacent state
  visual_->viz5State(getSparseState(vp), tools::LARGE, tools::PURPLE, 0);
  //visual_->viz5State(getSparseState(vp), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

  // Show edge between them
  visual_->viz5Edge(getSparseState(vp), getSparseState(v), tools::eGREEN);

  // Show adjacent state
  visual_->viz5State(getSparseState(vpp), tools::LARGE, tools::PURPLE, 0);
  //visual_->viz5State(getSparseState(vpp), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

  // Show edge between them
  visual_->viz5Edge(getSparseState(vpp), getSparseState(v), tools::eORANGE);
  visual_->viz5Trigger();
  usleep(0.001 * 1000000);

  // Show iData
  if (iData.hasInterface1())
  {
    visual_->viz5State(iData.interface1Inside_, tools::MEDIUM, tools::ORANGE, 0);
    visual_->viz5State(iData.interface1Outside_, tools::MEDIUM, tools::GREEN, 0);
    visual_->viz5Edge(iData.interface1Inside_, iData.interface1Outside_, tools::eRED);

    if (vp < vpp)
      visual_->viz5Edge(getSparseState(vp), iData.interface1Outside_, tools::eRED);
    else
      visual_->viz5Edge(getSparseState(vpp), iData.interface1Outside_, tools::eRED);
  }
  if (iData.hasInterface2())
  {
    visual_->viz5State(iData.interface2Inside_, tools::MEDIUM, tools::ORANGE, 0);
    visual_->viz5State(iData.interface2Outside_, tools::MEDIUM, tools::GREEN, 0);
    visual_->viz5Edge(iData.interface2Inside_, iData.interface2Outside_, tools::eRED);

    if (vp < vpp)
      visual_->viz5Edge(getSparseState(vpp), iData.interface2Outside_, tools::eRED);
    else
      visual_->viz5Edge(getSparseState(vp), iData.interface2Outside_, tools::eRED);
  }

  visual_->viz5Trigger();
  usleep(0.001 * 1000000);
}

bool SparseDB::addQualityPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                              std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "addQualityPath()");
  indent += 2;

  // TEMP:
  // if (iData.lastDistance_ == 0)
  // {
  //   BOLT_RED_DEBUG(0, 1, "last distance is zero!");
  //   visual_->waitForUserFeedback();
  // }


  // Can we connect these two vertices directly?
  if (si_->checkMotion(getSparseState(vp), getSparseState(vpp)))
  {
    BOLT_DEBUG(indent, vQuality_, "Adding edge between vp and vpp");

    SparseEdge e = addEdge(vp, vpp, eRED, indent + 2);

    if (edgeWeightPropertySparse_[e] > ignoreEdgesSmallerThan_) //discretization_ + small)
    {

      if (visualizeQualityCriteria_)
        visualizeCheckAddPath(v, vp, vpp, iData);

      // TEMP:
      // std::cout << "discretization_ + small: " << discretization_ + small << std::endl;
      // BOLT_DEBUG(0, true, "Spanner property violated, edge added of length " << edgeWeightPropertySparse_[e]);
      // visual_->waitForUserFeedback();
    }
  }
  else
  {
    BOLT_RED_DEBUG(indent, true, "Geometric path being added for spanner");

    geometric::PathGeometric *path = new geometric::PathGeometric(si_);
    if (vp < vpp)
    {
      path->append(getSparseState(vp));
      path->append(iData.interface1Outside_);
      path->append(iData.interface1Inside_);
      path->append(getSparseState(v));
      path->append(iData.interface2Inside_);
      path->append(iData.interface2Outside_);
      path->append(getSparseState(vpp));
    }
    else
    {
      path->append(getSparseState(vp));
      path->append(iData.interface2Outside_);
      path->append(iData.interface2Inside_);
      path->append(getSparseState(v));
      path->append(iData.interface1Inside_);
      path->append(iData.interface1Outside_);
      path->append(getSparseState(vpp));
    }

    // Visualize path
    if (visualizeQualityCriteria_ || true)
    {
      visual_->viz1DeleteAllMarkers();
      base::PathPtr pathPtr(new geometric::PathGeometric(*path));
      visual_->viz1Path(base::PathPtr(pathPtr), 1);
      visual_->viz1Trigger();
      usleep(0.001 * 1000000);
    }

    BOLT_DEBUG(indent, vQuality_, "Created path with " << path->getStateCount() << " states");

    pathSimplifier_->reduceVertices(*path, 200);
    pathSimplifier_->shortcutPath(*path, 200);
    pathSimplifier_->reduceVertices(*path, 200);

    // const double simplifyTime = 10;  // seconds
    // pathSimplifier_->simplify(*path, simplifyTime);
    // pathSimplifier_->reduceVertices(*path, 50);

    if (visualizeQualityCriteria_ || true)
    {
      base::PathPtr pathPtr(new geometric::PathGeometric(*path));
      visual_->viz1Path(base::PathPtr(pathPtr), 2);
      visual_->viz1Trigger();
      usleep(0.001 * 1000000);
    }

    std::pair<bool, bool> repairResult = path->checkAndRepair(100);
    if (!repairResult.first)
    {
      OMPL_WARN("Check and repair found invalid states");

      if (visualizeQualityCriteria_ || true)
      {
        base::PathPtr pathPtr(new geometric::PathGeometric(*path));
        visual_->viz1Path(base::PathPtr(pathPtr), 4);
        visual_->viz1Trigger();
        usleep(0.001 * 1000000);
      }
    }

    if (repairResult.second) // Repairing was successful
    {
      SparseVertex prior = vp;
      SparseVertex vnew;
      std::vector<base::State *> &states = path->getStates();

      BOLT_DEBUG(indent + 2, vQuality_, "Shortcuted path now has " << path->getStateCount() << " states");

      assert(states.size() >= 3);
      for (std::size_t i = 1; i < states.size() - 1; ++i)  // first and last states are vp and vpp, don't addVertex()
      {
        base::State *state = states[i];

        // Check if this vertex already exists
        if (si_->equalStates(getSparseState(v), state))
        {
          OMPL_ERROR("Add path state is the same!");
          exit(-1);
        }

        // no need to clone st, since we will destroy p; we just copy the pointer
        BOLT_DEBUG(indent + 2, vQuality_, "Adding node from shortcut path for QUALITY");
        vnew = addVertex(state, QUALITY);

        // TEMP
        BOLT_RED_DEBUG(indent + 2, true, "clearing nearby edges");

        if (visualizeQualityCriteria_) // TEMP
          visualizeCheckAddPath(v, vp, vpp, iData); // TEMP
        visual_->waitForUserFeedback();

        // Remove all edges from all vertices near our new vertex
        clearEdgesNearVertex(vnew);

        // TEMP
        visual_->waitForUserFeedback();

        assert(prior != vnew);
        addEdge(prior, vnew, eRED, indent + 2);
        prior = vnew;
      }
      // clear the states, so memory is not freed twice
      states.clear();

      assert(prior != vpp);
      addEdge(prior, vpp, eRED, indent + 2);
    }
    else
    {
      BOLT_RED_DEBUG(indent + 2, true, "check and repair failed?");
      exit(-1);
    }

    delete path;

    if (visualizeQualityCriteria_)
      visualizeCheckAddPath(v, vp, vpp, iData);

    // TEMP:
    // BOLT_DEBUG(indent, true, "spanner property violated, path added");
    // visual_->waitForUserFeedback();
  }

  // TEMP:
  // if (iData.lastDistance_ == 0)
  // {
  //   BOLT_RED_DEBUG(0, 1, "last distance is zero!");
  //   visual_->waitForUserFeedback();
  // }


  return true;
}

bool SparseDB::spannerTestOriginal(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                                   std::size_t indent)
{
  const bool verbose = true; // vQuality_

  // Computes all nodes which qualify as a candidate x for v, v', and v"
  double midpointPathLength = maxSpannerPath(v, vp, vpp, indent + 6);

  // Check if spanner property violated
  // if (iData.lastDistance_ == 0)  // DTC added zero check
  // {
  //   BOLT_RED_DEBUG(indent + 6, verbose, "iData.lastDistance_ is 0");
  // }
  if (stretchFactor_ * iData.lastDistance_ < midpointPathLength)
  {
    BOLT_YELLOW_DEBUG(indent + 6, verbose, "Spanner property violated");
    BOLT_DEBUG(indent + 8, verbose, "Sparse Graph Midpoint Length  = " << midpointPathLength);
    BOLT_DEBUG(indent + 8, verbose, "Spanner Path Length * Stretch = " << (stretchFactor_ * iData.lastDistance_));
    BOLT_DEBUG(indent + 10, verbose, "last distance = " << iData.lastDistance_);
    BOLT_DEBUG(indent + 10, verbose, "stretch factor = " << stretchFactor_);
    double rejectStretchFactor = midpointPathLength / iData.lastDistance_;
    BOLT_DEBUG(indent + 10, verbose, "to reject, stretch factor > " << rejectStretchFactor);

    return true;  // spannerPropertyWasViolated
  }
  else
    BOLT_DEBUG(indent + 6, vQuality_, "Spanner property not violated");

  return false;  // spannerPropertyWasViolated = false
}

bool SparseDB::spannerTestOuter(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                                std::size_t indent)
{
  // Computes all nodes which qualify as a candidate x for v, v', and v"
  double midpointPathLength = maxSpannerPath(v, vp, vpp, indent + 6);

  // Must have both interfaces to continue
  if (!iData.hasInterface1() || !iData.hasInterface2())
  {
    return false;
  }

  double newDistance = si_->distance(iData.interface1Outside_, iData.interface2Outside_);  // TODO(davetcoleman): cache?

  // Check if spanner property violated
  if (newDistance == 0)  // DTC added zero check
  {
    BOLT_RED_DEBUG(indent + 6, vQuality_, "new distance is 0");
    exit(-1);
  }
  else if (stretchFactor_ * newDistance < midpointPathLength)
  {
    BOLT_YELLOW_DEBUG(indent + 6, vQuality_, "Spanner property violated");
    BOLT_DEBUG(indent + 8, vQuality_, "Sparse Graph Midpoint Length  = " << midpointPathLength);
    BOLT_DEBUG(indent + 8, vQuality_, "Spanner Path Length * Stretch = " << (stretchFactor_ * newDistance));
    BOLT_DEBUG(indent + 10, vQuality_, "new distance = " << newDistance);
    BOLT_DEBUG(indent + 10, vQuality_, "stretch factor = " << stretchFactor_);

    return true;  // spannerPropertyWasViolated
  }
  else
    BOLT_DEBUG(indent + 6, vQuality_, "Spanner property not violated");

  return false;  // spannerPropertyWasViolated = false
}

bool SparseDB::spannerTestAStar(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                                std::size_t indent)
{
  if (iData.hasInterface1() && iData.hasInterface2())
  {
    BOLT_DEBUG(indent + 6, vQuality_,
               "Temp recalculated distance: " << si_->distance(iData.interface1Inside_, iData.interface2Inside_));

    // Experimental calculations
    double pathLength = 0;
    std::vector<SparseVertex> vertexPath;
    if (!astarSearch(vp, vpp, vertexPath, pathLength, indent + 6))
    {
      BOLT_RED_DEBUG(indent + 6, vQuality_, "No path found");
      usleep(4 * 1000000);
    }
    else
    {
      if (visualizeQualityCriteria_)
      {
        visual_->viz6DeleteAllMarkers();
        assert(vertexPath.size() > 1);
        for (std::size_t i = 1; i < vertexPath.size(); ++i)
        {
          visual_->viz6Edge(getSparseState(vertexPath[i - 1]), getSparseState(vertexPath[i]), tools::eGREEN);
        }
      }

      // Add connecting segments:
      double connector1 = si_->distance(getSparseState(vp), iData.getOutsideInterfaceOfV1(vp, vpp));
      // TODO(davetcoleman): may want to include the dist from inside to outside of interface
      BOLT_DEBUG(indent + 6, vQuality_, "connector1 " << connector1);
      if (visualizeQualityCriteria_)
      {
        visual_->viz6Edge(getSparseState(vp), iData.getOutsideInterfaceOfV1(vp, vpp), tools::eORANGE);
      }

      double connector2 = si_->distance(getSparseState(vpp), iData.getOutsideInterfaceOfV2(vp, vpp));
      // TODO(davetcoleman): may want to include the dist from inside to outside of interface
      BOLT_DEBUG(indent + 6, vQuality_, "connector2 " << connector2);
      if (visualizeQualityCriteria_)
      {
        visual_->viz6Edge(getSparseState(vpp), iData.getOutsideInterfaceOfV2(vp, vpp), tools::eYELLOW);
      }

      pathLength += connector1 + connector2;
      BOLT_DEBUG(indent + 6, vQuality_, "Full Path Length: " << pathLength);

      visual_->viz6Trigger();
    }

    if (iData.lastDistance_ == 0)
    {
      BOLT_YELLOW_DEBUG(indent + 6, vQuality_, "Last distance is 0");
    }

    // Theoretical max
    //double theoreticalMaxLength = iData.lastDistance_ + 2 * sparseDelta_;
    double theoreticalMaxLength = iData.lastDistance_ + sparseDelta_;
    BOLT_DEBUG(indent + 6, vQuality_, "Max allowable length: " << theoreticalMaxLength);
    theoreticalMaxLength *= stretchFactor_;
    BOLT_DEBUG(indent + 6, vQuality_, "Max allowable length with stretch: " << theoreticalMaxLength);

    if (pathLength < theoreticalMaxLength)
    {
      BOLT_DEBUG(indent + 6, vQuality_, "Astar says we do not need to add an edge");
    }
    else
    {
      BOLT_RED_DEBUG(indent + 6, vQuality_, "Astar says we need to add an edge");

      return true;  // spannerPropertyWasViolated = true
    }
  }

  return false;  // spannerPropertyWasViolated = false
}

SparseVertex SparseDB::findGraphRepresentative(base::State *state, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "findGraphRepresentative()");

  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;

  // Search
  getSparseStateNonConst(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  getSparseStateNonConst(queryVertices_[threadID]) = nullptr;

  BOLT_DEBUG(indent + 2, vQuality_, "Found " << graphNeighbors.size() << " nearest neighbors (graph rep) within "
                                                                         "SparseDelta " << sparseDelta_);

  SparseVertex result = boost::graph_traits<SparseGraph>::null_vertex();

  for (std::size_t i = 0; i < graphNeighbors.size(); ++i)
  {
    BOLT_DEBUG(indent + 2, vQuality_, "Checking motion of graph representative candidate " << i);
    if (si_->checkMotion(state, getSparseState(graphNeighbors[i])))
    {
      BOLT_DEBUG(indent + 4, vQuality_, "graph representative valid ");
      result = graphNeighbors[i];
      break;
    }
    else
      BOLT_DEBUG(indent + 4, vQuality_, "graph representative NOT valid, checking next ");
  }
  return result;
}

void SparseDB::findCloseRepresentatives(base::State *workState, const base::State *candidateState,
                                        const SparseVertex candidateRep,
                                        std::map<SparseVertex, base::State *> &closeRepresentatives, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "findCloseRepresentatives()");
  BOLT_DEBUG(indent + 2, vQuality_, "nearSamplePoints: " << nearSamplePoints_ << " denseDelta: " << denseDelta_);
  const bool visualizeSampler = false;
  base::State *sampledState = workState;  // rename variable just to clarify what it represents temporarily

  assert(closeRepresentatives.empty());

  // Search the space around new potential state candidateState
  for (std::size_t i = 0; i < nearSamplePoints_; ++i)
  {
    BOLT_DEBUG(indent + 2, vQuality_, "Get supporting representative #" << i);

    bool foundValidSample = false;
    static const std::size_t MAX_SAMPLE_ATTEMPT = 1000;
    for (std::size_t attempt = 0; attempt < MAX_SAMPLE_ATTEMPT; ++attempt)
    {
      BOLT_DEBUG(indent + 4, vQuality_, "Sample attempt " << attempt);

      sampler_->sampleNear(sampledState, candidateState, denseDelta_);
      //si_->getStateSpace()->setLevel(sampledState, 0);  // TODO no hardcode

      if (!si_->isValid(sampledState))
      {
        BOLT_DEBUG(indent + 6, vQuality_, "notValid ");

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);

        continue;
      }
      if (si_->distance(candidateState, sampledState) > denseDelta_)
      {
        BOLT_DEBUG(indent + 6, vQuality_, "Distance too far " << si_->distance(candidateState, sampledState)
                                                              << " needs to be less than " << denseDelta_);

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);
        continue;
      }
      if (!si_->checkMotion(candidateState, sampledState))
      {
        BOLT_DEBUG(indent + 6, vQuality_, "Motion invalid ");

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
      BOLT_DEBUG(indent + 4, vQuality_, "Unable to find valid sample after " << MAX_SAMPLE_ATTEMPT << " attempts"
                                                                                                      " ");
    }
    else
      BOLT_DEBUG(indent + 4, vQuality_, "Found valid nearby sample");

    // Compute which sparse vertex represents this new candidate vertex
    SparseVertex sampledStateRep = findGraphRepresentative(sampledState, indent + 6);

    // Check if sample is not visible to any other node (it should be visible in all likelihood)
    if (sampledStateRep == boost::graph_traits<SparseGraph>::null_vertex())
    {
      BOLT_DEBUG(indent + 4, vQuality_, "Sampled state has no representative (is null) ");

      // It can't be seen by anybody, so we should take this opportunity to add him
      BOLT_DEBUG(indent + 4, vQuality_, "Adding node for COVERAGE");
      addVertex(si_->cloneState(sampledState), COVERAGE);

      BOLT_DEBUG(indent + 4, vQuality_, "STOP EFFORTS TO ADD A DENSE PATH");

      // We should also stop our efforts to add a dense path
      for (std::map<SparseVertex, base::State *>::iterator it = closeRepresentatives.begin();
           it != closeRepresentatives.end(); ++it)
        si_->freeState(it->second);
      closeRepresentatives.clear();
      break;
    }

    BOLT_DEBUG(indent + 4, vQuality_, "Sampled state has representative (is not null)");

    // If its representative is different than candidateState
    if (sampledStateRep != candidateRep)
    {
      BOLT_DEBUG(indent + 4, vQuality_, "candidateRep != sampledStateRep ");

      // And we haven't already tracked this representative
      if (closeRepresentatives.find(sampledStateRep) == closeRepresentatives.end())
      {
        BOLT_DEBUG(indent + 4, vQuality_, "Track the representative");

        // Track the representative
        closeRepresentatives[sampledStateRep] = si_->cloneState(sampledState);
      }
      else
      {
        BOLT_DEBUG(indent + 4, vQuality_, "Already tracking the representative");
      }
    }
    else
    {
      BOLT_DEBUG(indent + 4, vQuality_, "candidateRep == sampledStateRep, do not keep this sample ");
    }
  }  // for each supporting representative
}

bool SparseDB::updatePairPoints(SparseVertex candidateRep, const base::State *candidateState,
                                SparseVertex nearSampledRep, const base::State *nearSampledState, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "updatePairPoints()");
  bool updated = false;  // track whether a change was made to any vertices' representatives

  // First of all, we need to compute all candidate r'
  std::vector<SparseVertex> adjVerticesUnconnected;
  getAdjVerticesOfV1UnconnectedToV2(candidateRep, nearSampledRep, adjVerticesUnconnected, indent + 2);

  // for each pair Pv(r,r')
  foreach (SparseVertex adjVertexUnconnected, adjVerticesUnconnected)
  {
    // Try updating the pair info
    if (distanceCheck(candidateRep, candidateState, nearSampledRep, nearSampledState, adjVertexUnconnected, indent + 2))
      updated = true;
  }

  return updated;
}

void SparseDB::getAdjVerticesOfV1UnconnectedToV2(SparseVertex v1, SparseVertex v2,
                                                 std::vector<SparseVertex> &adjVerticesUnconnected, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "getAdjVerticesOfV1UnconnectedToV2()");

  adjVerticesUnconnected.clear();
  foreach (SparseVertex adjVertex, boost::adjacent_vertices(v1, g_))
    if (adjVertex != v2)
      if (!hasEdge(adjVertex, v2))
        adjVerticesUnconnected.push_back(adjVertex);

  BOLT_DEBUG(indent + 2, vQuality_, "adjVerticesUnconnected size: " << adjVerticesUnconnected.size());
}

double SparseDB::maxSpannerPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "maxSpannerPath()");
  indent += 2;

  // Candidate x vertices as described in paper in Max_Spanner_Path
  std::vector<SparseVertex> qualifiedVertices;

  // Get nearby vertices 'x' that could also be used to find the path to v''
  foreach (SparseVertex x, boost::adjacent_vertices(vpp, g_))
  {
    if (hasEdge(x, v) && !hasEdge(x, vp))
    {
      InterfaceData &iData = getData(v, vpp, x, indent + 4);

      // Check if we previously had found a pair of points that support this interface
      if ((vpp < x && iData.interface1Inside_) || (x < vpp && iData.interface2Inside_))
      {
        BOLT_RED_DEBUG(indent, vQuality_, "Found an additional qualified vertex!");
        // This is a possible alternative path to v''
        qualifiedVertices.push_back(x);

        if (visualizeQualityCriteria_)
        {
          visual_->viz5State(getSparseState(x), tools::LARGE, tools::BLACK, 0);
        }
      }
    }
  }

  // vpp is always qualified because of its previous checks
  qualifiedVertices.push_back(vpp);
  BOLT_DEBUG(indent, vQuality_, "Total qualified vertices found: " << qualifiedVertices.size());

  // Find the maximum spanner distance
  BOLT_DEBUG(indent, vQuality_, "Finding the maximum spanner distance between v' and v''");
  double maxDist = 0.0;
  foreach (SparseVertex qualifiedVertex, qualifiedVertices)
  {
    BOLT_DEBUG(indent + 2, vQuality_, "Vertex: " << qualifiedVertex);
    if (visualizeQualityCriteria_)
      visual_->viz5State(getSparseState(qualifiedVertex), tools::SMALL, tools::PINK, 0);

    double tempDist = (si_->distance(getSparseState(vp), getSparseState(v)) +
                       si_->distance(getSparseState(v), getSparseState(qualifiedVertex))) /
                      2.0;
    // do we divide by 2 because of the midpoint path?? TODO(davetcoleman): figure out this proof
    if (tempDist > maxDist)
    {
      BOLT_DEBUG(indent + 4, vQuality_, "Is larger than previous");
      maxDist = tempDist;
    }
  }
  BOLT_DEBUG(indent, vQuality_, "Max distance: " << maxDist);

  return maxDist;
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
  BOLT_BLUE_DEBUG(indent, vQuality_, "getData() " << v << ", " << vp << ", " << vpp);
  return interfaceDataProperty_[v][index(vp, vpp)];
}

bool SparseDB::distanceCheck(SparseVertex v, const base::State *q, SparseVertex vp, const base::State *qp,
                             SparseVertex vpp, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "distanceCheck()");
  indent += 2;

  bool updated = false;  // track whether a change was made to any vertices' representatives

  // Get the info for the current representative-neighbors pair
  InterfaceData &iData = getData(v, vp, vpp, indent + 4);

  if (vp < vpp)  // FIRST points represent r (the interface discovered through sampling)
  {
    if (!iData.hasInterface1())  // No previous interface has been found here, just save it
    {
      BOLT_DEBUG(indent, vQuality_, "setInterface1");
      iData.setInterface1(q, qp, si_);
      updated = true;
    }
    else if (!iData.hasInterface2())  // The other interface doesn't exist, so we can't compare.
    {
      // Should probably keep the one that is further away from rep?  Not known what to do in this case.
      // TODO: is this not part of the algorithm?
      BOLT_YELLOW_DEBUG(indent, vQuality_, "TODO no interface 2");
    }
    else  // We know both of these points exist, so we can check some distances
    {
      assert(iData.lastDistance_ < std::numeric_limits<double>::infinity());
      if (si_->distance(q, iData.interface2Inside_) < iData.lastDistance_)
      // si_->distance(iData.interface1Inside_, iData.interface2Inside_))
      {  // Distance with the new point is good, so set it.
        BOLT_GREEN_DEBUG(indent, vQuality_, "setInterface1 UPDATED");
        iData.setInterface1(q, qp, si_);
        updated = true;
      }
      else
      {
        BOLT_DEBUG(indent, vQuality_, "Distance was not better, not updating bookkeeping");
      }
    }
  }
  else  // SECOND points represent r (the interfaec discovered through sampling)
  {
    if (!iData.hasInterface2())  // No previous interface has been found here, just save it
    {
      BOLT_DEBUG(indent, vQuality_, "setInterface2");
      iData.setInterface2(q, qp, si_);
      updated = true;
    }
    else if (!iData.hasInterface1())  // The other interface doesn't exist, so we can't compare.
    {
      // Should we be doing something cool here?
      BOLT_YELLOW_DEBUG(indent, vQuality_, "TODO no interface 1");
    }
    else  // We know both of these points exist, so we can check some distances
    {
      assert(iData.lastDistance_ < std::numeric_limits<double>::infinity());
      if (si_->distance(q, iData.interface1Inside_) < iData.lastDistance_)
      // si_->distance(iData.interface2Inside_, iData.interface1Inside_))
      {  // Distance with the new point is good, so set it
        BOLT_GREEN_DEBUG(indent, vQuality_, "setInterface2 UPDATED");
        iData.setInterface2(q, qp, si_);
        updated = true;
      }
      else
      {
        BOLT_DEBUG(indent, vQuality_, "Distance was not better, not updating bookkeeping");
      }
    }
  }

  // Lastly, save what we have discovered
  if (updated)
  {
    // TODO(davetcoleman): do we really need to copy this back in or is it already passed by reference?
    interfaceDataProperty_[v][index(vp, vpp)] = iData;
  }

  return updated;
}

void SparseDB::abandonLists(base::State *state)
{
  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;

  // Search
  getSparseStateNonConst(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  getSparseStateNonConst(queryVertices_[threadID]) = nullptr;

  // For each of the vertices
  std::size_t deletions = 0;
  foreach (SparseVertex v, graphNeighbors)
  {
    foreach (VertexPair r, interfaceDataProperty_[v] | boost::adaptors::map_keys)
    {
      interfaceDataProperty_[v][r].clear(si_);
      deletions++;
    }
  }
}

void SparseDB::clearEdgesNearVertex(SparseVertex vertex)
{
  // TODO(davetcoleman): combine this with abandonLists and ensure that all interface data is equally cleared
  // but do not clear out nearby edges if a non-quality-path vertex is added
  std::vector<SparseVertex> graphNeighbors;

  // Search
  nn_->nearestR(vertex, sparseDelta_, graphNeighbors);

  // For each of the vertices
  foreach (SparseVertex v, graphNeighbors)
  {
    // Remove all edges to and from vertex
    boost::clear_vertex(v, g_);
  }

  displaySparseDatabase();
}

void SparseDB::findGraphNeighbors(DenseVertex v1, std::vector<SparseVertex> &graphNeighborhood,
                                  std::vector<SparseVertex> &visibleNeighborhood, std::size_t threadID,
                                  std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "findGraphNeighbors()");
  const bool verbose = false;

  base::State *state = getDenseState(v1);

  // Search
  getSparseStateNonConst(queryVertices_[threadID]) = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighborhood);
  getSparseStateNonConst(queryVertices_[threadID]) = nullptr;

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

  BOLT_DEBUG(indent + 2, vCriteria_, "Graph neighborhood: " << graphNeighborhood.size() << " | Visible neighborhood: "
                                                            << visibleNeighborhood.size());
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

  // Clear all nearby interface data whenever a new vertex is added
  if (useFourthCriteria_)
    abandonLists(getDenseState(denseV));

  // Connected component tracking
  disjointSets_.make_set(v);

  // Add vertex to nearest neighbor structure
  nn_->add(v);

  // Book keeping for what was added
  switch (type)
  {
    case COVERAGE:
      numSamplesAddedForCoverage_++;
      break;
    case CONNECTIVITY:
      numSamplesAddedForConnectivity_++;
      break;
    case INTERFACE:
      numSamplesAddedForInterface_++;
      break;
    case QUALITY:
      numSamplesAddedForQuality_++;
      break;
    default:
      OMPL_ERROR("Unknown type");
  }

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
      OMPL_ERROR("Unknown type");

    // if (type == COVERAGE)
    if (visualizeDatabaseCoverage_)
      visual_->viz2State(getSparseState(v), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

    if (type == QUALITY)
      visual_->viz2State(getSparseState(v), tools::LARGE, color, 0);
    else
      visual_->viz2State(getSparseState(v), tools::MEDIUM, color, 0);

    visual_->viz4State(getSparseState(v), tools::MEDIUM, color, 0);

    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    if (visualizeOverlayNodes_)  // after initial spars graph is created, show additions not from grid
    {
      visual_->viz4State(getSparseState(v), tools::MEDIUM, color, sparseDelta_);
      visual_->viz4Trigger();
      usleep(0.001 * 1000000);
    }
  }

  if (visualizeVoronoiDiagramAnimated_ || (visualizeVoronoiDiagram_ && useFourthCriteria_))
    visual_->vizVoronoiDiagram();

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

SparseEdge SparseDB::addEdge(SparseVertex v1, SparseVertex v2, std::size_t visualColor, std::size_t indent)
{
  assert(v1 <= getNumVertices());
  assert(v2 <= getNumVertices());
  assert(v1 != v2);
  BOOST_ASSERT_MSG(getSparseState(v1) != getSparseState(v2), "States on both sides of an edge are the same");

  // Create the new edge
  SparseEdge e = (boost::add_edge(v1, v2, g_)).first;

  // Weight properties
  edgeWeightPropertySparse_[e] = distanceFunction(v1, v2);  // TODO: use this value with astar

  // Collision properties
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

    double small = penetrationDistance_ + std::numeric_limits<double>::epsilon();
    if (false) // Use alterative color scheme
    {
      if (edgeWeightPropertySparse_[e] < nearestDVertex_ + small)
      {
        visualColor = eYELLOW;
      }
      else if (edgeWeightPropertySparse_[e] < discretization_ + small)
      {
        visualColor = eGREEN;
      }
      else
      {
        visualColor = eRED;
      }
    }

    visual_->viz2Edge(getSparseState(v1), getSparseState(v2), visualColor);
    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    //if (visualColor == tools::eRED)
    {
      if (edgeWeightPropertySparse_[e] > ignoreEdgesSmallerThan_)
      {
        //std::cout << "Edge distance: " << edgeWeightPropertySparse_[e] << std::endl;
        visual_->viz4Edge(getSparseState(v1), getSparseState(v2), visualColor);
        visual_->viz4Trigger();
        usleep(0.001 * 1000000);
      }
    }
  }

  return e;
}

base::State *&SparseDB::getSparseStateNonConst(SparseVertex v)
{
  return denseDB_->stateProperty_[denseVertexProperty_[v]];
}

const base::State *SparseDB::getSparseState(SparseVertex v) const
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
    return;
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
      visual_->viz2Edge(getSparseState(v1), getSparseState(v2), tools::eRED);  // edgeWeightPropertySparse_[e]);

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
      if (getSparseState(v))
      {
        visual_->viz2State(getSparseState(v), tools::SMALL, tools::BLUE, 1);
        // visual_->viz2State(getSparseState(v), /*popularity0-100*/ 7, vertexPopularity_[v]);
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
  return si_->distance(getSparseState(a), getSparseState(b));
}

double SparseDB::getSecondarySparseDelta()
{
  return sparseDelta_ * 1.25;
}

bool SparseDB::hasEdge(SparseVertex v1, SparseVertex v2)
{
  return boost::edge(v1, v2, g_).second;
}

void SparseDB::visualizeInterfaces(SparseVertex v, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "visualizeInterfaces()");

  InterfaceHash &iData = interfaceDataProperty_[v];

  visual_->viz6DeleteAllMarkers();
  visual_->viz6State(getSparseState(v), tools::LARGE, tools::RED, 0);

  for (auto it = iData.begin(); it != iData.end(); ++it)
  {
    const VertexPair &pair = it->first;
    InterfaceData &iData = it->second;

    SparseVertex v1 = pair.first;
    SparseVertex v2 = pair.second;

    visual_->viz6State(getSparseState(v1), tools::LARGE, tools::PURPLE, 0);
    visual_->viz6State(getSparseState(v2), tools::LARGE, tools::PURPLE, 0);
    // visual_->viz6Edge(getSparseState(v1), getSparseState(v2), tools::eGREEN);

    if (iData.hasInterface1())
    {
      visual_->viz6State(iData.interface1Inside_, tools::MEDIUM, tools::ORANGE, 0);
      visual_->viz6State(iData.interface1Outside_, tools::MEDIUM, tools::GREEN, 0);
      visual_->viz6Edge(iData.interface1Inside_, iData.interface1Outside_, tools::eYELLOW);
    }
    else
    {
      visual_->viz6Edge(getSparseState(v1), getSparseState(v2), tools::eRED);
    }

    if (iData.hasInterface2())
    {
      visual_->viz6State(iData.interface2Inside_, tools::MEDIUM, tools::ORANGE, 0);
      visual_->viz6State(iData.interface2Outside_, tools::MEDIUM, tools::GREEN, 0);
      visual_->viz6Edge(iData.interface2Inside_, iData.interface2Outside_, tools::eYELLOW);
    }
    else
    {
      visual_->viz6Edge(getSparseState(v1), getSparseState(v2), tools::eRED);
    }
  }

  visual_->viz6Trigger();
  usleep(0.1 * 1000000);
}

void SparseDB::visualizeAllInterfaces(std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "visualizeAllInterfaces()");

  visual_->viz6DeleteAllMarkers();

  foreach (SparseVertex v, boost::vertices(g_))
  {
    // typedef std::unordered_map<VertexPair, InterfaceData> InterfaceHash;
    InterfaceHash &hash = interfaceDataProperty_[v];

    for (auto it = hash.begin(); it != hash.end(); ++it)
    {
      InterfaceData &iData = it->second;

      if (iData.hasInterface1())
      {
        visual_->viz6State(iData.interface1Inside_, tools::MEDIUM, tools::ORANGE, 0);
        visual_->viz6State(iData.interface1Outside_, tools::MEDIUM, tools::GREEN, 0);
      }

      if (iData.hasInterface2())
      {
        visual_->viz6State(iData.interface2Inside_, tools::MEDIUM, tools::ORANGE, 0);
        visual_->viz6State(iData.interface2Outside_, tools::MEDIUM, tools::GREEN, 0);
      }
    }
  }
  visual_->viz6Trigger();
  usleep(0.1 * 1000000);
}

std::pair<std::size_t, std::size_t> SparseDB::getInterfaceStateStorageSize()
{
  std::size_t numStates = 0;
  std::size_t numMissingInterfaces = 0;

  foreach (SparseVertex v, boost::vertices(g_))
  {
    InterfaceHash &hash = interfaceDataProperty_[v];

    for (auto it = hash.begin(); it != hash.end(); ++it)
    {
      InterfaceData &iData = it->second;

      if (iData.hasInterface1())
        numStates += 2;
      else
        numMissingInterfaces++;

      if (iData.hasInterface2())
        numStates += 2;
      else
        numMissingInterfaces++;

      if (!iData.hasInterface1() && !iData.hasInterface2())
        std::cout << "has neither interfaces! " << std::endl;
    }
  }
  return std::pair<std::size_t, std::size_t>(numStates, numMissingInterfaces);
}

SparseVertex SparseDB::getSparseRepresentative(base::State *state)
{
  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;
  const std::size_t numNeighbors = 1;

  // Search
  getSparseStateNonConst(queryVertices_[threadID]) = state;
  nn_->nearestK(queryVertices_[threadID], numNeighbors, graphNeighbors);
  getSparseStateNonConst(queryVertices_[threadID]) = nullptr;

  if (graphNeighbors.empty())
  {
    std::cout << "no neighbors found for sparse representative " << std::endl;
    exit(-1);
  }
  return graphNeighbors[0];
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
