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
#include <ompl/base/DiscreteMotionValidator.h>

// Boost
#include <boost/graph/incremental_components.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/assert.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

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
    parent_->getVisual()->viz4State(parent_->getVertexState(v), tools::SMALL, tools::GREEN, 1);
}

void otb::CustomAstarVisitor::examine_vertex(SparseVertex v, const SparseGraph &) const
{
  // Statistics
  parent_->recordNodeClosed();

  if (parent_->visualizeAstar_)
  {
    parent_->getVisual()->viz4State(parent_->getVertexState(v), tools::LARGE, tools::BLACK, 1);
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
SparseDB::SparseDB(base::SpaceInformationPtr si, VisualizerPtr visual)
  : si_(si)
  , visual_(visual)
  // Property accessors of edges
  , edgeWeightProperty_(boost::get(boost::edge_weight, g_))
  , edgeTypeProperty_(boost::get(edge_type_t(), g_))
  , edgeCollisionStatePropertySparse_(boost::get(edge_collision_state_t(), g_))
  // Property accessors of vertices
  , stateCacheProperty_(boost::get(vertex_state_cache_t(), g_))
  , vertexTypeProperty_(boost::get(vertex_type_t(), g_))
  , interfaceDataProperty_(boost::get(vertex_interface_data_t(), g_))
  , vertexPopularity_(boost::get(vertex_popularity_t(), g_))
  // Disjoint set accessors
  , disjointSets_(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_))
{
  // Save number of threads available
  numThreads_ = boost::thread::hardware_concurrency();

  // Initialize collision cache
  denseCache_.reset(new DenseCache(si_, this, visual_));

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
  clearanceSampler_.reset();
  regularSampler_.reset();

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
  const double discFactor = sparseDelta_ - discretizePenetrationDist_;
  discretization_ = 2 * sqrt(std::pow(discFactor, 2) / dim);

  ignoreEdgesSmallerThan_ = (discretization_ + 0.01);

  // Calculate optimum stretch factor
  if (stretchFactor_ < std::numeric_limits<double>::epsilon())  // if stretchFactor is zero, auto set it
  {
    BOLT_DEBUG(indent, 1, "Auto settings stretch factor because input value was 0");
    // stretchFactor_ = discretization_ / (0.5 * discretization_ * sqrt(2) - 2.0 * denseDelta_);
    nearestDiscretizedV_ = sqrt(dim * std::pow(0.5 * discretization_, 2));  // z in my calculations
    // stretchFactor_ = 2.0 * discretization_ / ( nearestDiscretizedV_ - 2.0 * denseDelta_) + stretchFactor_; // 2D case
    // but not 3D
    stretchFactor_ =
        2.0 * discretization_ / (nearestDiscretizedV_) + stretchFactor_;  // 2D case without estimated interface amount
    // stretchFactor_ = (discretization_ + nearestDiscretizedV_) / ( 2 * (nearestDiscretizedV_ - 2.0 * denseDelta_)); //
    // N-D case
    // stretchFactor_ = discretization_ / (discretization_ - 2.0 * denseDelta_); // N-D case
  }

  BOLT_DEBUG(indent, 1, "--------------------------------------------------");
  BOLT_DEBUG(indent, 1, "Sparse DB Setup:");
  BOLT_DEBUG(indent + 2, 1, "Max Extent              = " << maxExtent_);
  BOLT_DEBUG(indent + 2, 1, "Sparse Delta            = " << sparseDelta_);
  BOLT_DEBUG(indent + 2, 1, "Dense Delta             = " << denseDelta_);
  BOLT_DEBUG(indent + 2, 1, "State Dimension         = " << dim);
  BOLT_DEBUG(indent + 2, 1, "Near Sample Points      = " << nearSamplePoints_);
  BOLT_DEBUG(indent + 2, 1, "Discretization          = " << discretization_);
  BOLT_DEBUG(indent + 2, 1, "Nearest Discretized V   = " << nearestDiscretizedV_);
  BOLT_DEBUG(indent + 2, 1, "Stretch Factor          = " << stretchFactor_);
  BOLT_DEBUG(indent + 2, 1, "Viz ignore edges below  = " << ignoreEdgesSmallerThan_);
  BOLT_DEBUG(indent, 1, "--------------------------------------------------");

  assert(maxExtent_ > 0);
  assert(denseDelta_ > 0);
  assert(nearSamplePoints_ > 0);
  assert(sparseDelta_ > 0);
  assert(sparseDelta_ > 0.000000001);  // Sanity check

  // Load minimum clearance state sampler
  if (!clearanceSampler_)
  {
    clearanceSampler_ = ob::MinimumClearanceValidStateSamplerPtr(new ob::MinimumClearanceValidStateSampler(si_.get()));
    clearanceSampler_->setMinimumObstacleClearance(obstacleClearance_);
    si_->getStateValidityChecker()->setClearanceSearchDistance(obstacleClearance_);
  }

  // Load regular state sampler
  if (!regularSampler_)
  {
    regularSampler_ = si_->allocValidStateSampler();
  }

  if (si_->getStateValidityChecker()->getClearanceSearchDistance() < obstacleClearance_)
    OMPL_WARN("State validity checker clearance search distance %f is less than the required obstacle clearance %f for "
              "our state sampler, incompatible settings!",
              si_->getStateValidityChecker()->getClearanceSearchDistance(), obstacleClearance_);

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

bool SparseDB::load()
{
  // Load collision cache
  denseCache_->load();

  // Benchmark
  time::point start = time::now();

  BoltStorage storage_(si_, this);
  if (!storage_.load(filePath_.c_str()))
    return false;

  // Benchmark
  double duration = time::seconds(time::now() - start);

  // Error check
  if (!getNumVertices() || !getNumEdges())
  {
    OMPL_ERROR("Corrupted sparse graph loaded");
    return false;
  }

  // Get the average vertex degree (number of connected edges)
  double averageDegree = (getNumEdges() * 2) / static_cast<double>(getNumVertices());

  // Check how many disjoint sets are in the sparse graph (should be none)
  std::size_t numSets = checkConnectedComponents();

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Loaded graph stats:");
  OMPL_INFORM("   Total vertices:         %u", getNumVertices());
  OMPL_INFORM("   Total edges:            %u", getNumEdges());
  OMPL_INFORM("   Average degree:         %f", averageDegree);
  OMPL_INFORM("   Connected Components:   %u", numSets);
  OMPL_INFORM("   Loading time:           %f", duration);
  OMPL_INFORM("------------------------------------------------------");

  // Disable
  OMPL_INFORM("Disabling discretized samples generation because we have loaded from file");
  useDiscretizedSamples_ = false;

  return true;
}

bool SparseDB::saveIfChanged()
{
  if (graphUnsaved_)
  {
    return save();
  }
  else
    OMPL_INFORM("Not saving because database has not changed");
  return true;
}

bool SparseDB::save()
{
  if (!graphUnsaved_)
    OMPL_WARN("No need to save because graphUnsaved_ is false, but saving anyway because requested");

  // Disabled
  if (!savingEnabled_)
  {
    OMPL_INFORM("Not saving because option disabled for SparseDB");
    return false;
  }

  // Error checking
  if (filePath_.empty())
  {
    OMPL_ERROR("Empty filename passed to save function");
    return false;
  }

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
    for (std::size_t i = numThreads_; i < getNumVertices(); ++i)  // skip vertex 0-11 because those are query vertices
    {
      const SparseVertex v1 = i;
      const SparseVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4Edge(getVertexState(v1), getVertexState(v2), 10);
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
  double dist = si_->distance(getVertexState(a), getVertexState(b));

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

  // Create a query state for each possible thread
  queryVertices_.resize(numThreads_);
  queryState_.resize(numThreads_);

  for (std::size_t threadID = 0; threadID < numThreads_; ++threadID)
  {
    // Add a fake vertex to the graph
    queryVertices_[threadID] = boost::add_vertex(g_);
  }
}

void SparseDB::clearEdgeCollisionStates()
{
  foreach (const SparseEdge e, boost::edges(g_))
    edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;  // each edge has an unknown state
}

void SparseDB::createSPARS()
{
  std::size_t indent = 0;
  BOLT_BLUE_DEBUG(indent, true, "createSPARS()");

  // Benchmark runtime
  time::point startTime = time::now();

  numSamplesAddedForQuality_ = 0;
  numSamplesAddedForConnectivity_ = 0;
  numSamplesAddedForInterface_ = 0;
  numSamplesAddedForQuality_ = 0;
  numVerticesMoved_ = 0;

  numConsecutiveFailures_ = 0;
  useFourthCriteria_ = false;  // initially we do not do this step

  // Profiler
  CALLGRIND_TOGGLE_COLLECT;

  // Start the graph off with discretized states
  if (useDiscretizedSamples_)
  {
    addDiscretizedStates(indent + 2);
  }

  // Finish the graph with random samples
  if (useRandomSamples_)
  {
    addRandomSamples(indent + 2);
  }

  if (!useRandomSamples_ && !useDiscretizedSamples_)
  {
    OMPL_WARN("Unable to create SPARS because both random sampling and discretized sampling is disabled");
  }

  // Profiler
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
  BOLT_DEBUG(0, 1, "  Generation time:           " << duration);
  BOLT_DEBUG(0, 1, "  Total generations:         " << numGraphGenerations_);
  BOLT_DEBUG(0, 1, "  Disjoint sets:             " << numSets);
  BOLT_DEBUG(0, 1, "  DenseCache                 ");
  BOLT_DEBUG(0, 1, "    Edge cache         ");
  BOLT_DEBUG(0, 1, "      Size:                    " << denseCache_->getEdgeCacheSize());
  BOLT_DEBUG(0, 1, "      Total checks:            " << denseCache_->getTotalCollisionChecks());
  BOLT_DEBUG(0, 1, "      Cached checks:           " << denseCache_->getTotalCollisionChecksFromCache() << " ("
                                                     << denseCache_->getPercentCachedCollisionChecks() << "%)");
  BOLT_DEBUG(0, 1, "    State cache                ");
  BOLT_DEBUG(0, 1, "      Size:                    " << denseCache_->getStateCacheSize());
  BOLT_DEBUG(0, 1, "  Criteria additions:        ");
  BOLT_DEBUG(0, 1, "    Coverage:                " << numSamplesAddedForCoverage_);
  BOLT_DEBUG(0, 1, "    Connectivity:            " << numSamplesAddedForConnectivity_);
  BOLT_DEBUG(0, 1, "    Interface:               " << numSamplesAddedForInterface_);
  BOLT_DEBUG(0, 1, "    Quality:                 " << numSamplesAddedForQuality_);
  BOLT_DEBUG(0, 1, "  Num random samples added:  " << numRandSamplesAdded_);
  BOLT_DEBUG(0, 1, "  Num vertices moved:        " << numVerticesMoved_);
  BOLT_DEBUG(0, 1, "  InterfaceData:             ");
  BOLT_DEBUG(0, 1, "    States stored:           " << interfaceStats.first);
  BOLT_DEBUG(0, 1, "    Missing interfaces:      " << interfaceStats.second);
  BOLT_DEBUG(0, 1, "-----------------------------------------");

  if (!visualizeSparsGraph_)
    displayDatabase(true, indent + 2);

  OMPL_INFORM("Finished creating sparse database");
}

void SparseDB::addDiscretizedStates(std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, true, "addDiscretizedStates()");
  indent += 2;
  bool verbose = true;

  if (vertexDiscretizer_->getDiscretization() < std::numeric_limits<double>::epsilon())
  {
    OMPL_WARN("Discretization not set");
    exit(-1);
  }

  discretizedSamplesInsertion_ = true;  // this tells SPARS to always add the vertex, no matter what
  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID];
  const std::size_t jointIncrements = floor(range / vertexDiscretizer_->getDiscretization());
  double leftOver = range - jointIncrements * vertexDiscretizer_->getDiscretization();
  double startOffset = leftOver / 2;

  BOLT_DEBUG(indent, verbose, "------------------------------------------");
  BOLT_DEBUG(indent, verbose, "Discretization:       " << vertexDiscretizer_->getDiscretization());
  BOLT_DEBUG(indent, verbose, "High Bound:           " << bounds.high[jointID]);
  BOLT_DEBUG(indent, verbose, "Low Bound:            " << bounds.low[jointID]);
  BOLT_DEBUG(indent, verbose, "Range:                " << range);
  BOLT_DEBUG(indent, verbose, "Joint Increments:     " << jointIncrements);
  BOLT_DEBUG(indent, verbose, "Left Over:            " << leftOver);
  BOLT_DEBUG(indent, verbose, "Start Offset:         " << startOffset);
  BOLT_DEBUG(indent, verbose, "------------------------------------------");

  // Create two levels of grids
  for (std::size_t i = 0; i < 2; ++i)
  {
    BOLT_DEBUG(indent, verbose, "Discretize iteration " << i);

    // Set starting value offset
    if (i == 0)
      vertexDiscretizer_->setStartingValueOffset(startOffset);
    else
      vertexDiscretizer_->setStartingValueOffset(startOffset + vertexDiscretizer_->getDiscretization() / 2.0);

    // Generate vertices
    vertexDiscretizer_->generate(indent + 2);

    // Convert to proper format TODO(davetcoleman): remove the need for this format?
    std::vector<base::State *> &candidateVertices = vertexDiscretizer_->getCandidateVertices();

    std::list<WeightedVertex> vertexInsertionOrder;
    for (base::State *state : candidateVertices)
    {
      // Move the ompl::base::State to the DenseCache, changing its ownership
      StateID candidateStateID = denseCache_->addState(state);
      vertexInsertionOrder.push_back(WeightedVertex(candidateStateID, 0));
    }
    candidateVertices.clear();  // clear the vector because we've moved all its memory pointers to DenseCache
    std::size_t sucessfulInsertions;
    createSPARSInnerLoop(vertexInsertionOrder, sucessfulInsertions);
  }

  discretizedSamplesInsertion_ = false;

  BOLT_DEBUG(indent + 2, verbose, "Finished discretization");
  BOLT_DEBUG(indent, verbose, "------------------------------------------\n");
}

/*
void SparseDB::createSPARSOuterLoop()
{
  std::size_t indent = 2;

  // Reset parameters
  setup();
  visualizeOverlayNodes_ = false;  // DO NOT visualize all added nodes in a separate window
  denseCache_->resetCounters();

  // Get the ordering to insert vertices
  std::list<WeightedVertex> vertexInsertionOrder;
  getVertexInsertionOrdering(vertexInsertionOrder);

  // Error check order creation
  assert(vertexInsertionOrder.size() == getNumVertices() - queryVertices_.size());

  // Attempt to insert the vertices multiple times until no more succesful insertions occur
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
              << " loop, remaining uninserted vertices: " << vertexInsertionOrder.size()
              << " loop runtime: " << duration << " sec" << std::endl;
    loopAttempt++;

    // Increase the sparse delta a bit, but only after the first loop
    if (loopAttempt == 1)
    {
      // sparseDelta_ = getSecondarySparseDelta();
      std::cout << std::string(indent + 2, ' ') << "sparseDelta_ is now " << sparseDelta_ << std::endl;
      secondSparseInsertionAttempt_ = true;

      // Save collision cache, just in case there is a bug
      denseCache_->save();
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
    displayDatabase(true, indent+4);
  }
  else if (visualizeSparsGraphSpeed_ < std::numeric_limits<double>::epsilon())
  {
    visual_->viz2Trigger();
    usleep(0.001 * 1000000);
  }
}
*/

bool SparseDB::createSPARSInnerLoop(std::list<WeightedVertex> &vertexInsertionOrder, std::size_t &sucessfulInsertions)
{
  std::size_t indent = 0;

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
                << "% Cache size: " << denseCache_->getEdgeCacheSize()
                << " Cache usage: " << denseCache_->getPercentCachedCollisionChecks() << "%" << std::endl;
      std::cout << ANSI_COLOR_RESET;
      if (visualizeSparsGraph_)
        visual_->viz2Trigger();
    }

    // Run SPARS checks
    VertexType addReason;    // returns why the state was added
    SparseVertex newVertex;  // the newly generated sparse vertex

    // TODO(davetcoleman): I would like to multi-thread this but its not worth my time currently
    std::size_t threadID = 0;
    if (!addStateToRoadmap(vertexIt->stateID_, newVertex, addReason, threadID, indent))
    {
      // std::cout << "Failed AGAIN to add state to roadmap------" << std::endl;

      // Visualize the failed vertex as a small red dot
      if (visualizeSparsGraph_ && false)
      {
        visual_->viz2State(denseCache_->getState(vertexIt->stateID_), tools::SMALL, tools::RED, 0);
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

void SparseDB::addRandomSamples(std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, true, "addRandomSamples()");
  indent += 2;

  // Clear stats
  numRandSamplesAdded_ = 0;

  // Add from file if available - remember where it was
  std::size_t lastCachedStateIndex = denseCache_->getStateCacheSize();
  std::size_t currentStateCacheID = 1;  // skip 0, because that is the "deleted" NULL state ID
  bool usingCachedStates = true;        // allow to trigger user feedback

  if (useDiscretizedSamples_)
  {
    BOLT_YELLOW_DEBUG(indent, true, "Not using cache for random samples because discretized samples is enabled.");
    usingCachedStates = false;
  }

  while (true)
  {
    StateID candidateStateID;

    // Try from cache first
    if (usingCachedStates && currentStateCacheID < lastCachedStateIndex)
    {
      // Add from cache
      candidateStateID = currentStateCacheID;
      BOLT_DEBUG(indent, vCriteria_, "Adding from CACHE: " << currentStateCacheID << " total " << lastCachedStateIndex
                                                           << " stateID: " << candidateStateID);

      currentStateCacheID++;
    }
    else  // Add new state
    {
      if (usingCachedStates)
      {
        usingCachedStates = false;
        BOLT_YELLOW_DEBUG(indent, vCriteria_, "addRandomSamples: No longer using cached states - sampling new "
                                              "ones");
      }

      base::State *candidateState = si_->allocState();
      candidateStateID = denseCache_->addState(candidateState);
      BOLT_DEBUG(indent, vCriteria_, "NOT adding from cache, currentStateCacheID: "
                                         << currentStateCacheID << " lastCachedStateIndex: " << lastCachedStateIndex
                                         << " stateID: " << candidateStateID);

      // Sample randomly
      if (!clearanceSampler_->sample(candidateState))
      {
        OMPL_ERROR("Unable to find valid sample");
        exit(-1);  // this should never happen
      }
      // si_->getStateSpace()->setLevel(candidateState, 0);  // TODO no hardcode
    }

    // Run SPARS checks
    VertexType addReason;    // returns why the state was added
    SparseVertex newVertex;  // the newly generated sparse vertex
    const std::size_t threadID = 0;
    if (addStateToRoadmap(candidateStateID, newVertex, addReason, threadID, indent + 2))
    {
      // if (numRandSamplesAdded_ % 10 == 0)
      BOLT_DEBUG(indent, vCriteria_, "Added random sample with stateID "
                                         << candidateStateID << ", total new states: " << ++numRandSamplesAdded_);
    }
    else if (numConsecutiveFailures_ % 500 == 0)
    {
      BOLT_DEBUG(indent, true, "Random sample failed, consecutive failures: " << numConsecutiveFailures_);
    }

    // Check consecutive failures
    if (numConsecutiveFailures_ >= fourthCriteriaAfterFailures_ && !useFourthCriteria_)
    {
      BOLT_YELLOW_DEBUG(indent, true, "Starting to check for 4th quality criteria because "
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
      BOLT_YELLOW_DEBUG(indent, true, "SPARS creation finished because " << terminateAfterFailures_
                                                                         << " consecutive insertion failures reached");
      return;
    }
  }  // end while
}

/*
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

bool SparseDB::getPopularityOrder(std::list<WeightedVertex> &vertexInsertionOrder)
{
  bool verbose = false;

  // Error check
  BOOST_ASSERT_MSG(getNumVertices() > queryVertices_.size(),
                   "Unable to get vertices in order of popularity because dense "
                   "graph is empty");

  if (visualizeNodePopularity_)  // Clear visualization
  {
    visual_->viz3DeleteAllMarkers();
  }

  // Sort the vertices by popularity in a queue
  std::priority_queue<WeightedVertex, std::vector<WeightedVertex>, CompareWeightedVertex> pqueue;

  // Loop through each popular edge in the dense graph
  foreach (DenseVertex v, boost::vertices(g_))
  {
    // Do not process the search vertex, it is null
    if (v <= queryVertices_.back())
      continue;

    if (verbose)
      std::cout << "Vertex: " << v << std::endl;
    double popularity = 0;
    foreach (DenseEdge edge, boost::out_edges(v, g_))
    {
      if (verbose)
        std::cout << "  Edge: " << edge << std::endl;
      popularity += (100 - edgeWeightProperty_[edge]);
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
      visual_->viz3State(stateProperty_[pqueue.top().v_], tools::SCALE, tools::BLACK, weightPercent);
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
  foreach (DenseVertex v, boost::vertices(g_))
  {
    // Do not process the search vertex, it is null
    if (v <= queryVertices_.back())
      continue;

    if (verbose)
      std::cout << "Vertex: " << v << std::endl;
    double popularity = 0;

    foreach (DenseEdge edge, boost::out_edges(v, g_))
    {
      if (verbose)
        std::cout << "  Edge: " << edge << std::endl;
      popularity += (100 - edgeWeightProperty_[edge]);
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
      visual_->viz3State(stateProperty_[wv.v_], tools::SCALE, tools::BLACK, weightPercent);
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
*/

bool SparseDB::addStateToRoadmap(StateID candidateStateID, SparseVertex &newVertex, VertexType &addReason,
                                 std::size_t threadID, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "addStateToRoadmap() Adding candidate state ID " << candidateStateID);

  if (visualizeAttemptedStates_)
  {
    visual_->viz1DeleteAllMarkers();
    visual_->viz1State(getState(candidateStateID), tools::LARGE, tools::GREEN, 0);
    visual_->viz1Trigger();
    usleep(0.001 * 1000000);
  }

  bool stateAdded = false;
  base::State *workState = si_->allocState();

  // Nodes near our input state
  std::vector<SparseVertex> graphNeighborhood;
  // Visible nodes near our input state
  std::vector<SparseVertex> visibleNeighborhood;

  // Find nearby nodes
  findGraphNeighbors(candidateStateID, graphNeighborhood, visibleNeighborhood, threadID, indent + 2);

  // Always add a node if no other nodes around it are visible (GUARD)
  if (checkAddCoverage(candidateStateID, visibleNeighborhood, newVertex, indent + 2))
  {
    BOLT_DEBUG(indent + 2, vAddedReason_, "Graph updated for: COVERAGE");

    addReason = COVERAGE;
    stateAdded = true;
  }
  else if (checkAddConnectivity(candidateStateID, visibleNeighborhood, newVertex, indent + 6))
  {
    BOLT_DEBUG(indent + 2, vAddedReason_, "Graph updated for: CONNECTIVITY");

    addReason = CONNECTIVITY;
    stateAdded = true;
  }
  else if (checkAddInterface(candidateStateID, graphNeighborhood, visibleNeighborhood, newVertex, indent + 10))
  {
    BOLT_DEBUG(indent + 2, vAddedReason_, "Graph updated for: INTERFACE");

    addReason = INTERFACE;
    stateAdded = true;
  }
  else if (useFourthCriteria_ &&
           checkAddQuality(candidateStateID, graphNeighborhood, visibleNeighborhood, workState, newVertex, indent + 14))
  {
    BOLT_DEBUG(indent + 2, vAddedReason_, "Graph updated for: QUALITY");

    addReason = QUALITY;
    stateAdded = true;
  }
  else if (discretizedSamplesInsertion_)
  {
    BOLT_DEBUG(indent + 2, vAddedReason_, "Graph updated for: DISCRETIZED");
    newVertex = addVertex(candidateStateID, DISCRETIZED, indent + 2);

    addReason = DISCRETIZED;
    stateAdded = true;
  }
  else
  {
    BOLT_DEBUG(indent + 2, vCriteria_, "Did NOT add state for any criteria ");
    numConsecutiveFailures_++;
  }

  if (stateAdded)
    numConsecutiveFailures_ = 0;

  si_->freeState(workState);

  return stateAdded;
}

bool SparseDB::checkAddCoverage(StateID candidateStateID, std::vector<SparseVertex> &visibleNeighborhood,
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

  newVertex = addVertex(candidateStateID, COVERAGE, indent + 4);

  // Note: we do not connect this node with any edges because we have already determined
  // it is too far away from any nearby nodes

  return true;
}

bool SparseDB::checkAddConnectivity(StateID candidateStateID, std::vector<SparseVertex> &visibleNeighborhood,
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
        BOLT_DEBUG(indent + 2, vCriteria_, "Different connected component: " << visibleNeighborhood[i] << ", "
                                                                             << visibleNeighborhood[j]);

        if (visualizeConnectivity_)
        {
          visual_->viz1State(getVertexState(visibleNeighborhood[i]), tools::MEDIUM, tools::BLUE, 0);
          visual_->viz1State(getVertexState(visibleNeighborhood[j]), tools::MEDIUM, tools::BLUE, 0);
          visual_->viz1Trigger();
          usleep(0.001 * 1000000);
        }

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
  newVertex = addVertex(candidateStateID, CONNECTIVITY, indent + 4);

  // Check if there are really close vertices nearby which should be merged
  checkRemoveCloseVertices(newVertex, indent + 4);

  // Add the edges
  for (std::set<SparseVertex>::const_iterator vertexIt = statesInDiffConnectedComponents.begin();
       vertexIt != statesInDiffConnectedComponents.end(); ++vertexIt)
  {
    if (stateCacheProperty_[*vertexIt] == 0)
    {
      BOLT_RED_DEBUG(indent + 4, vCriteria_, "Skipping because vertex " << *vertexIt << " was removed (state marked as "
                                                                                        "0)");
      // visual_->waitForUserFeedback("skipping because vertex was removed");
      continue;
    }

    // Do not add edge from self to self
    if (si_->getStateSpace()->equalStates(getVertexState(*vertexIt), getVertexState(newVertex)))
    {
      std::cout << "Prevented same vertex from being added twice " << std::endl;
      continue;  // skip this pairing
    }

    BOLT_DEBUG(indent + 4, vCriteria_, "Loop: Adding vertex " << *vertexIt);

    // New vertex should not be connected to anything - there's no edge between the two states
    if (hasEdge(newVertex, *vertexIt) == true)
    {
      BOLT_DEBUG(indent + 4, 1, "The new vertex " << newVertex << " is already connected to old vertex");
      continue;
    }

    // The components haven't been united by previous edges created in this for loop
    if (!sameComponent(*vertexIt, newVertex))
    {
      // Connect
      addEdge(newVertex, *vertexIt, eCONNECTIVITY, indent + 4);
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

bool SparseDB::checkAddInterface(StateID candidateStateID, std::vector<SparseVertex> &graphNeighborhood,
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

  // If the two closest nodes are also visible
  const std::size_t threadID = 0;
  if (graphNeighborhood[0] == visibleNeighborhood[0] && graphNeighborhood[1] == visibleNeighborhood[1])
  {
    // If our two closest neighbors don't share an edge
    if (!hasEdge(visibleNeighborhood[0], visibleNeighborhood[1]))
    {
      // If they can be directly connected
      if (denseCache_->checkMotionWithCacheVertex(visibleNeighborhood[0], visibleNeighborhood[1], threadID))
      {
        BOLT_DEBUG(indent + 2, vCriteria_, "INTERFACE: directly connected nodes");

        // Connect them
        addEdge(visibleNeighborhood[0], visibleNeighborhood[1], eINTERFACE, indent + 4);

        // Also add the vertex if we are in a special mode where we know its desired
        if (discretizedSamplesInsertion_)
          newVertex = addVertex(candidateStateID, DISCRETIZED, indent + 4);
      }
      else  // They cannot be directly connected
      {
        // Add the new node to the graph, to bridge the interface
        BOLT_DEBUG(indent + 2, vCriteria_, "Adding node for INTERFACE");

        newVertex = addVertex(candidateStateID, INTERFACE, indent + 4);

        // Check if there are really close vertices nearby which should be merged
        if (checkRemoveCloseVertices(newVertex, indent + 4))
        {
          // New vertex replaced a nearby vertex, we can continue no further because graph has been re-indexed
          return true;
        }

        if (getVertexState(visibleNeighborhood[0]) == NULL)
        {
          BOLT_RED_DEBUG(indent + 3, 1, "Skipping edge 0 because vertex was removed");
          visual_->waitForUserFeedback("skipping edge 0");
        }
        else
          addEdge(newVertex, visibleNeighborhood[0], eINTERFACE, indent + 4);

        if (getVertexState(visibleNeighborhood[1]) == NULL)
        {
          BOLT_RED_DEBUG(indent + 3, 1, "Skipping edge 1 because vertex was removed");
          visual_->waitForUserFeedback("skipping edge 2");
        }
        else
          addEdge(newVertex, visibleNeighborhood[1], eINTERFACE, indent + 4);

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

bool SparseDB::checkAddQuality(StateID candidateStateID, std::vector<SparseVertex> &graphNeighborhood,
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

  base::State *candidateState = denseCache_->getStateNonConst(candidateStateID);  // paper's name: q
  SparseVertex candidateRep = visibleNeighborhood[0];                             // paper's name: v

  bool added = false;
  std::map<SparseVertex, base::State *> closeRepresentatives;  // [nearSampledRep, nearSampledState]
  findCloseRepresentatives(workState, candidateStateID, candidateRep, closeRepresentatives, indent + 2);

  BOLT_DEBUG(indent, vQuality_, "back in checkAddQuality(): Found " << closeRepresentatives.size()
                                                                    << " close representatives");

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
      visual_->viz3Edge(getVertexState(nearSampledRep), nearSampledState, tools::eRED);

      visual_->viz3State(nearSampledState, tools::MEDIUM, tools::GREEN, 0);

      // Replicate a regular vertex visualization
      // visual_->viz3State(getVertexState(nearSampledRep), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT,
      // sparseDelta_);
      visual_->viz3State(getVertexState(nearSampledRep), tools::LARGE, tools::PURPLE, sparseDelta_);

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
    return false;
  }

  // Visualize the interfaces around the candidate rep
  if (visualizeQualityCriteria_)
  {
    // visualizeCheckAddQuality(candidateState, candidateRep);
    // visualizeInterfaces(candidateRep, indent + 2);

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

void SparseDB::visualizeCheckAddQuality(StateID candidateStateID, SparseVertex candidateRep)
{
  visual_->viz3DeleteAllMarkers();

  visual_->viz3Edge(denseCache_->getState(candidateStateID), getVertexState(candidateRep), tools::eORANGE);

  // Show candidate state
  // visual_->viz3State(denseCache_->getState(candidateStateID), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT,
  // denseDelta_);
  visual_->viz3State(denseCache_->getState(candidateStateID), tools::LARGE, tools::RED, 0);

  // Show candidate state's representative
  visual_->viz3State(getVertexState(candidateRep), tools::LARGE, tools::BLUE, 0);

  // Show candidate state's representative's neighbors
  foreach (SparseVertex adjVertex, boost::adjacent_vertices(candidateRep, g_))
  {
    visual_->viz3Edge(getVertexState(adjVertex), getVertexState(candidateRep), tools::eGREEN);
    visual_->viz3State(getVertexState(adjVertex), tools::LARGE, tools::PURPLE, 0);
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
      if (spannerTestOriginal(v, vp, vpp, iData, indent + 2))
      // if (spannerTestOuter(v, vp, vpp, iData, indent + 2))
      // if (spannerTestAStar(v, vp, vpp, iData, indent + 2))
      {
        // Actually add the vertices and possibly edges
        if (addQualityPath(v, vp, vpp, iData, indent + 6))
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

void SparseDB::visualizeCheckAddPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                                     std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "visualizeCheckAddPath()");
  visual_->viz5DeleteAllMarkers();

  // Show candidate rep
  visual_->viz5State(getVertexState(v), tools::LARGE, tools::BLUE, 0);

  // Show adjacent state
  visual_->viz5State(getVertexState(vp), tools::LARGE, tools::PURPLE, 0);
  // visual_->viz5State(getVertexState(vp), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

  // Show edge between them
  visual_->viz5Edge(getVertexState(vp), getVertexState(v), tools::eGREEN);

  // Show adjacent state
  visual_->viz5State(getVertexState(vpp), tools::LARGE, tools::PURPLE, 0);
  // visual_->viz5State(getVertexState(vpp), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

  // Show edge between them
  visual_->viz5Edge(getVertexState(vpp), getVertexState(v), tools::eORANGE);
  visual_->viz5Trigger();
  usleep(0.001 * 1000000);

  // Show iData
  if (iData.hasInterface1())
  {
    // std::cout << "hasInterface1 " << std::flush;
    visual_->viz5State(iData.getInterface1Inside(), tools::MEDIUM, tools::ORANGE, 0);
    // std::cout << "1 " << std::flush;
    // std::cout << "state: " << iData.getInterface1Outside() << std::flush;
    visual_->viz5State(iData.getInterface1Outside(), tools::MEDIUM, tools::GREEN, 0);
    // std::cout << " 2 " << std::flush;
    visual_->viz5Edge(iData.getInterface1Inside(), iData.getInterface1Outside(), tools::eRED);
    // std::cout << "3 " << std::endl;

    if (vp < vpp)
      visual_->viz5Edge(getVertexState(vp), iData.getInterface1Outside(), tools::eRED);
    else
      visual_->viz5Edge(getVertexState(vpp), iData.getInterface1Outside(), tools::eRED);
  }
  if (iData.hasInterface2())
  {
    visual_->viz5State(iData.getInterface2Inside(), tools::MEDIUM, tools::ORANGE, 0);
    visual_->viz5State(iData.getInterface2Outside(), tools::MEDIUM, tools::GREEN, 0);
    visual_->viz5Edge(iData.getInterface2Inside(), iData.getInterface2Outside(), tools::eRED);

    if (vp < vpp)
      visual_->viz5Edge(getVertexState(vpp), iData.getInterface2Outside(), tools::eRED);
    else
      visual_->viz5Edge(getVertexState(vp), iData.getInterface2Outside(), tools::eRED);
  }

  visual_->viz5Trigger();
  usleep(0.001 * 1000000);
}

bool SparseDB::addQualityPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                              std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "addQualityPath()");
  indent += 2;

  // Can we connect these two vertices directly?
  const std::size_t threadID = 0;
  if (denseCache_->checkMotionWithCacheVertex(vp, vpp, threadID))
  {
    BOLT_DEBUG(indent, vQuality_, "Adding edge between vp and vpp");

    if (hasEdge(vp, vpp) || hasEdge(vpp, vp))
    {
      OMPL_ERROR("Already has an edge!");
      visual_->waitForUserFeedback("has edge");
      exit(-1);
    }

    SparseEdge e = addEdge(vp, vpp, eQUALITY, indent + 2);

    if (edgeWeightProperty_[e] > ignoreEdgesSmallerThan_)  // discretization_ + small)
    {
      if (visualizeQualityCriteria_)
        visualizeCheckAddPath(v, vp, vpp, iData, indent + 4);

      // TEMP:
      // std::cout << "discretization_ + small: " << discretization_ + small << std::endl;
      // BOLT_DEBUG(0, true, "Spanner property violated, edge added of length " << edgeWeightProperty_[e]);
      // visual_->waitForUserFeedback();
    }

    return true;
  }

  // Super debug
  if (superDebug_ && si_->checkMotion(getVertexState(vp), getVertexState(vpp)))
  {
    OMPL_ERROR("Failed test - edge was in collision in cache, but not from scratch");
    visual_->waitForUserFeedback("error");
  }

  BOLT_YELLOW_DEBUG(indent, visualizeQualityPathSimp_, "Unable to connect directly - geometric path must be created "
                                                       "for spanner");

  if (visualizeQualityCriteria_)                           // TEMP
    visualizeCheckAddPath(v, vp, vpp, iData, indent + 4);  // TEMP

  geometric::PathGeometric *path = new geometric::PathGeometric(si_);

  // Populate path
  if (vp < vpp)
  {
    path->append(getVertexState(vp));
    path->append(iData.getInterface1Outside());
    path->append(iData.getInterface1Inside());
    path->append(getVertexState(v));
    path->append(iData.getInterface2Inside());
    path->append(iData.getInterface2Outside());
    path->append(getVertexState(vpp));
  }
  else
  {
    path->append(getVertexState(vp));
    path->append(iData.getInterface2Outside());
    path->append(iData.getInterface2Inside());
    path->append(getVertexState(v));
    path->append(iData.getInterface1Inside());
    path->append(iData.getInterface1Outside());
    path->append(getVertexState(vpp));
  }

  // Create path and simplify
  if (useOriginalSmoother_)
    smoothQualityPathOriginal(path, indent + 4);
  else
    smoothQualityPath(path, indent + 4);

  // Insert simplified path into graph
  SparseVertex prior = vp;
  SparseVertex newVertex;
  std::vector<base::State *> &states = path->getStates();

  BOLT_DEBUG(indent + 2, vQuality_, "Shortcuted path now has " << path->getStateCount() << " states");

  if (states.size() < 3)
  {
    BOLT_RED_DEBUG(indent + 2, true, "Somehow path was shrunk to less than three vertices: " << states.size());
    visual_->waitForUserFeedback("path shrunk to two");
    delete path;
    return false;
  }

  bool addEdgeEnabled = true;                          // if a vertex is skipped, stop adding edges
  for (std::size_t i = 1; i < states.size() - 1; ++i)  // first and last states are vp and vpp, don't addVertex()
  {
    base::State *state = states[i];

    // Check if this vertex already exists
    if (si_->distance(getVertexState(v), state) < denseDelta_)  // TODO: is it ok to re-use this denseDelta var?
    {
      BOLT_RED_DEBUG(indent + 2, vQuality_, "Add path state is too similar to v!");

      if (visualizeQualityPathSimp_)
      {
        visual_->viz1DeleteAllMarkers();
        visual_->viz1Path(path, 1, tools::RED);
        visual_->viz1Trigger();
        usleep(0.001 * 1000000);

        visual_->waitForUserFeedback("Add path state is too similar to v");
      }

      delete path;
      return false;
    }

    // Check if new vertex has enough clearance
    if (!sufficientClearance(state))
    {
      BOLT_YELLOW_DEBUG(indent + 2, true, "Skipped adding vertex in new path b/c insufficient clearance");
      visual_->waitForUserFeedback("insufficient clearance");
      addEdgeEnabled = false;
      continue;
    }

    // no need to clone state, since we will destroy p; we just copy the pointer
    BOLT_DEBUG(indent + 2, vQuality_, "Adding node from shortcut path for QUALITY");
    newVertex = addVertex(si_->cloneState(state), QUALITY, indent + 4);

    // Check if there are really close vertices nearby which should be merged
    if (checkRemoveCloseVertices(newVertex, indent + 4))
    {
      // New vertex replaced a nearby vertex, we can continue no further because graph has been re-indexed

      // Remove all edges from all vertices near our new vertex
      clearEdgesNearVertex(newVertex);

      // TODO should we clearEdgesNearVertex before return true ?
      delete path;
      return true;
    }

    // Remove all edges from all vertices near our new vertex
    clearEdgesNearVertex(newVertex);

    if (addEdgeEnabled)
    {
      assert(prior != newVertex);
      addEdge(prior, newVertex, eQUALITY, indent + 2);
      prior = newVertex;
    }
  }

  // clear the states, so memory is not freed twice
  states.clear();

  if (addEdgeEnabled)
  {
    assert(prior != vpp);
    addEdge(prior, vpp, eQUALITY, indent + 2);
  }

  delete path;

  if (visualizeQualityCriteria_)
    visualizeCheckAddPath(v, vp, vpp, iData, indent + 4);

  return true;
}

bool SparseDB::smoothQualityPathOriginal(geometric::PathGeometric *path, std::size_t indent)
{
  BOLT_RED_DEBUG(indent, vQuality_ || 1, "smoothQualityPathOriginal()");
  indent += 2;

  // Visualize path
  if (visualizeQualityPathSimp_)
  {
    visual_->viz1DeleteAllMarkers();
    visual_->viz1Path(path, 1, tools::BLUE);
    visual_->viz1Trigger();
    usleep(0.001 * 1000000);
  }

  BOLT_DEBUG(indent, visualizeQualityPathSimp_, "Created 'quality path' candidate with " << path->getStateCount()
                                                                                         << " states");
  if (visualizeQualityPathSimp_)
    visual_->waitForUserFeedback("path simplification");

  pathSimplifier_->reduceVertices(*path, 10);
  pathSimplifier_->shortcutPath(*path, 50);

  std::pair<bool, bool> repairResult = path->checkAndRepair(100);

  if (!repairResult.second)  // Repairing was not successful
  {
    BOLT_RED_DEBUG(indent + 2, true, "check and repair failed?");
    exit(-1);
  }
  return true;
}

bool SparseDB::smoothQualityPath(geometric::PathGeometric *path, std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vQuality_, "smoothQualityPath()");
  indent += 2;

  // Visualize path
  if (visualizeQualityPathSimp_)
  {
    visual_->viz1DeleteAllMarkers();
    visual_->viz1Path(path, 1, tools::BLUE);
    visual_->viz1Trigger();
    usleep(0.001 * 1000000);
  }

  BOLT_DEBUG(indent, visualizeQualityPathSimp_, "Created 'quality path' candidate with " << path->getStateCount()
                                                                                         << " states");
  if (visualizeQualityPathSimp_)
    visual_->waitForUserFeedback("path simplification");

  // Set the motion validator to use clearance, this way isValid() checks clearance before confirming valid
  base::DiscreteMotionValidator *dmv =
      dynamic_cast<base::DiscreteMotionValidator *>(si_->getMotionValidatorNonConst().get());
  dmv->setRequiredStateClearance(obstacleClearance_);

  for (std::size_t i = 0; i < 3; ++i)
  {
    pathSimplifier_->simplifyMax(*path);

    if (visualizeQualityPathSimp_)
    {
      visual_->viz1DeleteAllMarkers();
      visual_->viz1Path(path, 1, tools::ORANGE);
      visual_->viz1Trigger();
      usleep(0.1 * 1000000);
      // visual_->waitForUserFeedback("optimizing path");
    }

    pathSimplifier_->reduceVertices(*path, 1000, path->getStateCount() * 4);

    if (visualizeQualityPathSimp_)
    {
      visual_->viz1DeleteAllMarkers();
      visual_->viz1Path(path, 1, tools::BLUE);
      visual_->viz1Trigger();
      usleep(0.1 * 1000000);
      // visual_->waitForUserFeedback("optimizing path");
    }
  }
  // Turn off the clearance requirement
  dmv->setRequiredStateClearance(0.0);

  pathSimplifier_->reduceVertices(*path, 1000, path->getStateCount() * 4);

  if (visualizeQualityPathSimp_)
  {
    visual_->viz1DeleteAllMarkers();
    visual_->viz1Path(path, 1, tools::GREEN);
    visual_->viz1Trigger();
    visual_->waitForUserFeedback("finished quality path");
  }

  std::pair<bool, bool> repairResult = path->checkAndRepair(100);

  if (!repairResult.second)  // Repairing was not successful
  {
    BOLT_RED_DEBUG(indent + 2, true, "check and repair failed?");
    exit(-1);
  }
  return true;
}

bool SparseDB::spannerTestOriginal(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData,
                                   std::size_t indent)
{
  const bool verbose = vQuality_;

  // Computes all nodes which qualify as a candidate x for v, v', and v"
  double midpointPathLength = maxSpannerPath(v, vp, vpp, indent + 6);

  // Check if spanner property violated
  // if (iData.getLastDistance() == 0)  // DTC added zero check
  // {
  //   BOLT_RED_DEBUG(indent + 6, verbose, "iData.getLastDistance() is 0");
  // }
  if (stretchFactor_ * iData.getLastDistance() < midpointPathLength)
  {
    BOLT_YELLOW_DEBUG(indent + 6, verbose, "Spanner property violated");
    BOLT_DEBUG(indent + 8, verbose, "Sparse Graph Midpoint Length  = " << midpointPathLength);
    BOLT_DEBUG(indent + 8, verbose, "Spanner Path Length * Stretch = " << (stretchFactor_ * iData.getLastDistance()));
    BOLT_DEBUG(indent + 10, verbose, "last distance = " << iData.getLastDistance());
    BOLT_DEBUG(indent + 10, verbose, "stretch factor = " << stretchFactor_);
    double rejectStretchFactor = midpointPathLength / iData.getLastDistance();
    BOLT_DEBUG(indent + 10, verbose, "to reject, stretch factor > " << rejectStretchFactor);

    return true;  // spanner property was violated
  }

  BOLT_DEBUG(indent + 6, vQuality_, "Spanner property not violated");
  return false;  // spanner property was NOT violated
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

  double newDistance =
      si_->distance(iData.getInterface1Outside(), iData.getInterface2Outside());  // TODO(davetcoleman): cache?

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
    BOLT_DEBUG(indent + 6, vQuality_, "Temp recalculated distance: " << si_->distance(iData.getInterface1Inside(),
                                                                                      iData.getInterface2Inside()));

    // Experimental calculations
    double pathLength = 0;
    std::vector<SparseVertex> vertexPath;
    if (!astarSearch(vp, vpp, vertexPath, pathLength, indent + 6))
    {
      BOLT_RED_DEBUG(indent + 6, vQuality_, "No path found");
      visual_->waitForUserFeedback("No path found");
    }
    else
    {
      if (visualizeQualityCriteria_)
      {
        visual_->viz6DeleteAllMarkers();
        assert(vertexPath.size() > 1);
        for (std::size_t i = 1; i < vertexPath.size(); ++i)
        {
          visual_->viz6Edge(getVertexState(vertexPath[i - 1]), getVertexState(vertexPath[i]), tools::eGREEN);
        }
      }

      // Add connecting segments:
      double connector1 = si_->distance(getVertexState(vp), iData.getOutsideInterfaceOfV1(vp, vpp));
      // TODO(davetcoleman): may want to include the dist from inside to outside of interface
      BOLT_DEBUG(indent + 6, vQuality_, "connector1 " << connector1);
      if (visualizeQualityCriteria_)
      {
        visual_->viz6Edge(getVertexState(vp), iData.getOutsideInterfaceOfV1(vp, vpp), tools::eORANGE);
      }

      double connector2 = si_->distance(getVertexState(vpp), iData.getOutsideInterfaceOfV2(vp, vpp));
      // TODO(davetcoleman): may want to include the dist from inside to outside of interface
      BOLT_DEBUG(indent + 6, vQuality_, "connector2 " << connector2);
      if (visualizeQualityCriteria_)
      {
        visual_->viz6Edge(getVertexState(vpp), iData.getOutsideInterfaceOfV2(vp, vpp), tools::eYELLOW);
      }

      pathLength += connector1 + connector2;
      BOLT_DEBUG(indent + 6, vQuality_, "Full Path Length: " << pathLength);

      visual_->viz6Trigger();
    }

    if (iData.getLastDistance() == 0)
    {
      BOLT_YELLOW_DEBUG(indent + 6, vQuality_, "Last distance is 0");
    }

    // Theoretical max
    // double theoreticalMaxLength = iData.getLastDistance() + 2 * sparseDelta_;
    double theoreticalMaxLength = iData.getLastDistance() + sparseDelta_;
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
  queryState_[threadID] = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  queryState_[threadID] = nullptr;

  BOLT_DEBUG(indent + 2, vQuality_, "Found " << graphNeighbors.size() << " nearest neighbors (graph rep) within "
                                                                         "SparseDelta " << sparseDelta_);

  SparseVertex result = boost::graph_traits<SparseGraph>::null_vertex();

  for (std::size_t i = 0; i < graphNeighbors.size(); ++i)
  {
    BOLT_DEBUG(indent + 2, vQuality_, "Checking motion of graph representative candidate " << i);
    if (si_->checkMotion(state, getVertexState(graphNeighbors[i])))
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

void SparseDB::findCloseRepresentatives(base::State *workState, const StateID candidateStateID,
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

      regularSampler_->sampleNear(sampledState, denseCache_->getState(candidateStateID), denseDelta_);
      // si_->getStateSpace()->setLevel(sampledState, 0);  // TODO no hardcode

      if (!si_->isValid(sampledState))
      {
        BOLT_DEBUG(indent + 6, vQuality_, "notValid ");

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);

        continue;
      }
      if (si_->distance(denseCache_->getState(candidateStateID), sampledState) > denseDelta_)
      {
        BOLT_DEBUG(indent + 6, vQuality_, "Distance too far "
                                              << si_->distance(denseCache_->getState(candidateStateID), sampledState)
                                              << " needs to be less than " << denseDelta_);

        if (visualizeQualityCriteria_ && visualizeSampler)
          visual_->viz3State(sampledState, tools::SMALL, tools::RED, 0);
        continue;
      }
      if (!si_->checkMotion(denseCache_->getState(candidateStateID), sampledState))
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

    if (visualizeQualityCriteria_ && visualizeSampler)
    {
      visual_->viz3Trigger();
      usleep(0.001 * 1000000);
    }

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
      // But first check for proper clearance
      if (sufficientClearance(sampledState))
      {
        BOLT_DEBUG(indent + 4, vQuality_, "Adding node for COVERAGE");
        addVertex(si_->cloneState(sampledState), COVERAGE, indent + 4);
      }

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
  if (true)
  {
    foreach (SparseVertex x, boost::adjacent_vertices(vpp, g_))
    {
      if (hasEdge(x, v) && !hasEdge(x, vp))
      {
        InterfaceData &iData = getData(v, vpp, x, indent + 2);

        // Check if we previously had found a pair of points that support this interface
        if ((vpp < x && iData.getInterface1Inside()) || (x < vpp && iData.getInterface2Inside()))
        {
          BOLT_RED_DEBUG(indent, vQuality_, "Found an additional qualified vertex " << x);
          // This is a possible alternative path to v''
          qualifiedVertices.push_back(x);

          if (visualizeQualityCriteria_)
          {
            visual_->viz5State(getVertexState(x), tools::LARGE, tools::BLACK, 0);

            // Temp?
            // visual_->viz6DeleteAllMarkers();
            // visual_->viz6State(getVertexState(x), tools::LARGE, tools::BLACK, 0);
            // visual_->viz6Trigger();
          }
        }
      }
    }
  }
  else
    BOLT_YELLOW_DEBUG(indent, vQuality_, "Disabled nearby vertices in maxSpannerPath");

  // vpp is always qualified because of its previous checks
  qualifiedVertices.push_back(vpp);
  BOLT_DEBUG(indent, vQuality_, "Total qualified vertices found: " << qualifiedVertices.size());

  // Find the maximum spanner distance
  BOLT_DEBUG(indent, vQuality_, "Finding the maximum spanner distance between v' and v''");
  double maxDist = 0.0;
  foreach (SparseVertex qualifiedVertex, qualifiedVertices)
  {
    if (visualizeQualityCriteria_)
      visual_->viz5State(getVertexState(qualifiedVertex), tools::SMALL, tools::PINK, 0);

    // Divide by 2 because of the midpoint path 'M'
    double tempDist = (si_->distance(getVertexState(vp), getVertexState(v)) +
                       si_->distance(getVertexState(v), getVertexState(qualifiedVertex))) /
                      2.0;
    BOLT_DEBUG(indent + 2, vQuality_, "Checking vertex: " << qualifiedVertex << " distance: " << tempDist);

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
      assert(iData.getLastDistance() < std::numeric_limits<double>::infinity());
      if (si_->distance(q, iData.getInterface2Inside()) < iData.getLastDistance())
      // si_->distance(iData.getInterface1Inside(), iData.getInterface2Inside()))
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
      assert(iData.getLastDistance() < std::numeric_limits<double>::infinity());
      if (si_->distance(q, iData.getInterface1Inside()) < iData.getLastDistance())
      // si_->distance(iData.getInterface2Inside(), iData.getInterface1Inside()))
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
  queryState_[threadID] = state;
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighbors);
  queryState_[threadID] = nullptr;

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
  // Optionally disable this feature
  if (!useClearEdgesNearVertex_)
    return;

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

  std::size_t indent = 0;

  // Only display database if enabled
  if (visualizeSparsGraph_ && visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    displayDatabase(true, indent);
}

void SparseDB::findGraphNeighbors(StateID candidateStateID, std::vector<SparseVertex> &graphNeighborhood,
                                  std::vector<SparseVertex> &visibleNeighborhood, std::size_t threadID,
                                  std::size_t indent)
{
  BOLT_BLUE_DEBUG(indent, vCriteria_, "findGraphNeighbors()");
  const bool verbose = false;

  // Search
  queryState_[threadID] = denseCache_->getStateNonConst(candidateStateID);
  nn_->nearestR(queryVertices_[threadID], sparseDelta_, graphNeighborhood);
  queryState_[threadID] = nullptr;

  // Now that we got the neighbors from the NN, we must remove any we can't see
  for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  {
    SparseVertex v2 = graphNeighborhood[i];

    // Don't collision check if they are the same state
    if (candidateStateID != stateCacheProperty_[v2])
    {
      if (!denseCache_->checkMotionWithCache(candidateStateID, stateCacheProperty_[v2], threadID))
      {
        continue;
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

bool SparseDB::checkRemoveCloseVertices(SparseVertex v1, std::size_t indent)
{
  // This feature can be optionally disabled
  if (!useCheckRemoveCloseVertices_)
    return false;

  BOLT_CYAN_DEBUG(indent, vRemoveClose_, "checkRemoveCloseVertices() v1 = " << v1);
  indent += 2;

  // Get neighbors
  std::vector<SparseVertex> graphNeighbors;
  const std::size_t numNeighbors = 2;  // the first neighbor is (usually?) the vertex itself
  nn_->nearestK(v1, numNeighbors, graphNeighbors);

  if (vRemoveClose_)
  {
    std::stringstream ss;
    std::copy(graphNeighbors.begin(), graphNeighbors.end(), std::ostream_iterator<SparseVertex>(ss, ","));
    BOLT_DEBUG(indent, vRemoveClose_, "Neighbors of " << v1 << " are: " << ss.str());
  }

  // Error check no neighbors
  if (graphNeighbors.size() <= 1)
  {
    BOLT_RED_DEBUG(indent, vRemoveClose_, "checkRemoveCloseVertices: no neighbors found for sparse vertex " << v1);
    return false;
  }

  SparseVertex v2 = graphNeighbors[1];
  double sparseDeltaFractionCheck = 0.5;  // 0.25;  // TODO: determine better value for this

  // Error check: Do not remove itself
  if (v1 == v2)
  {
    BOLT_RED_DEBUG(indent, vRemoveClose_, "checkRemoveCloseVertices: error occured, cannot remove self: " << v1);
    exit(-1);
  }

  // Check that nearest neighbor is not an interface node - these should not be moved
  if (vertexTypeProperty_[v2] == QUALITY)
  {
    if (visualizeRemoveCloseVertices_)
    {
      visualizeRemoveCloseVertices(v1, v2);
      visual_->waitForUserFeedback("Skipping this vertex because is QUALITY");
    }
    return false;
  }

  // Check if nearest neighbor is within distance threshold
  if (distanceFunction(v1, v2) > sparseDelta_ * sparseDeltaFractionCheck)
  {
    // BOLT_DEBUG(indent, vRemoveClose_, "Distance " << distanceFunction(v1, v2) << " is greater than max "
    //<< sparseDelta_ * sparseDeltaFractionCheck);
    return false;
  }

  // Check if nearest neighbor is collision free
  if (!si_->checkMotion(getVertexState(v1), getVertexState(v2)))
  {
    BOLT_RED_DEBUG(indent, vRemoveClose_, "checkRemoveCloseVertices: not collision free v1=" << v1 << ", v2=" << v2);
    return false;
  }

  BOLT_DEBUG(indent, vRemoveClose_, "Found visible nearby node, testing if able to replace " << v2 << " with " << v1);

  // Nearest neighbor is good candidate, next check if all of its connected neighbors can be connected to new vertex
  foreach (SparseEdge edge, boost::out_edges(v2, g_))
  {
    SparseVertex v3 = boost::target(edge, g_);
    BOLT_DEBUG(indent + 2, vRemoveClose_, "checking edge v1= " << v1 << " to v3=" << v3);

    // Check if distance is within sparseDelta
    if (si_->distance(getVertexState(v1), getVertexState(v3)) > sparseDelta_)
    {
      BOLT_DEBUG(indent + 2, vRemoveClose_, "checkRemoveCloseVertices: distance too far from previous neighbor " << v3);
      return false;
    }

    // Check if collision free path to connected vertex
    if (!si_->checkMotion(getVertexState(v1), getVertexState(v3)))
    {
      BOLT_RED_DEBUG(indent + 2, vRemoveClose_,
                     "checkRemoveCloseVertices: not collision free path from new vertex to potential neighbor " << v3);
      return false;
    }
  }

  BOLT_DEBUG(indent, vRemoveClose_, "Found qualified node to replace with nearby");

  if (visualizeRemoveCloseVertices_)
  {
    visualizeRemoveCloseVertices(v1, v2);
    visual_->waitForUserFeedback("found qualified node to replace with nearby");
  }

  // Remove all interface data for old state
  abandonLists(getVertexStateNonConst(v2));

  // Connect new vertex to old vertex's connected neighbors
  foreach (SparseEdge edge, boost::out_edges(v2, g_))
  {
    SparseVertex v3 = boost::target(edge, g_);
    BOLT_DEBUG(indent + 2, vRemoveClose_, "Connecting v1= " << v1 << " to v3=" << v3);
    addEdge(v1, v3, eINTERFACE, indent + 4);
  }

  // Delete old vertex
  removeVertex(v2);
  BOLT_DEBUG(indent, vRemoveClose_, "REMOVING VERTEX " << v2 << " which was replaced with " << v1 << " with state ");

  // Statistics
  numVerticesMoved_++;

  // Only display database if enabled
  if (visualizeSparsGraph_ && visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    displayDatabase(true, indent + 4);

  if (visualizeRemoveCloseVertices_)
    visual_->waitForUserFeedback("finished moving vertex");

  if (visualizeRemoveCloseVertices_)
  {
    visual_->viz6DeleteAllMarkers();
    visual_->viz6Trigger();
    usleep(0.001 * 1000000);
  }

  return true;
}

void SparseDB::visualizeRemoveCloseVertices(SparseVertex v1, SparseVertex v2)
{
  visual_->viz6DeleteAllMarkers();
  visual_->viz6State(getVertexState(v1), tools::LARGE, tools::GREEN, 0);
  visual_->viz6State(getVertexState(v2), tools::LARGE, tools::RED, 0);  // RED = to be removed
  visual_->viz6Trigger();
  usleep(0.001 * 1000000);
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

void SparseDB::getDisjointSets(DisjointSetsParentKey &disjointSets)
{
  OMPL_INFORM("Get disjoint sets...");
  disjointSets.clear();

  // Flatten the parents tree so that the parent of every element is its representative.
  disjointSets_.compress_sets(boost::vertices(g_).first, boost::vertices(g_).second);

  // Count size of each disjoint set and group its containing vertices
  typedef boost::graph_traits<SparseGraph>::vertex_iterator VertexIterator;
  for (VertexIterator v = boost::vertices(g_).first; v != boost::vertices(g_).second; ++v)
  {
    // Do not count the search vertex within the sets
    if (*v <= queryVertices_.back())
      continue;

    disjointSets[boost::get(boost::get(boost::vertex_predecessor, g_), *v)].push_back(*v);
  }
}

void SparseDB::printDisjointSets(DisjointSetsParentKey &disjointSets)
{
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const SparseVertex v = iterator->first;
    const std::size_t freq = iterator->second.size();
    std::cout << "Parent: " << v << " frequency " << freq << std::endl;
  }
}

void SparseDB::visualizeDisjointSets(DisjointSetsParentKey &disjointSets)
{
  OMPL_INFORM("Visualizing disjoint sets");

  // Find the disjoint set that is the 'main' large one
  std::size_t maxDisjointSetSize = 0;
  SparseVertex maxDisjointSetParent;
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const SparseVertex v = iterator->first;
    const std::size_t freq = iterator->second.size();

    if (freq > maxDisjointSetSize)
    {
      maxDisjointSetSize = freq;
      maxDisjointSetParent = v;
    }
  }
  OMPL_INFORM("The largest disjoint set is of size %u and parent vertex %u", maxDisjointSetSize, maxDisjointSetParent);

  // Display size of disjoint sets and visualize small ones
  for (DisjointSetsParentKey::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const SparseVertex v1 = iterator->first;
    const std::size_t freq = iterator->second.size();
    // std::cout << "Parent vertex: " << v1 << " StateID: " << getStateID(v1) << " Frequency: " << freq << std::endl;
    // debugState(getVertexState(v1));

    BOOST_ASSERT_MSG(freq > 0, "Frequnecy must be at least 1");

    if (freq == maxDisjointSetSize)  // any subgraph that is smaller than the full graph
      continue;                      // the main disjoint set is not considered a disjoint set

    // Visualize sets of size one
    if (freq == 1)
    {
      // visual_->viz5State(getVertexState(v1), tools::ROBOT, tools::RED, 0);
      visual_->viz5State(getVertexState(v1), tools::MEDIUM, tools::RED, 0);
      visual_->viz5Trigger();
      visual_->waitForUserFeedback("showing disjoint set");
    }

    // Visualize large disjoint sets (greater than one)
    if (freq > 1 && freq < 1000)
    {
      // Clear markers
      visual_->viz4DeleteAllMarkers();

      // Visualize this subgraph that is disconnected
      // Loop through every every vertex and check if its part of this group
      typedef boost::graph_traits<SparseGraph>::vertex_iterator VertexIterator;
      for (VertexIterator v2 = boost::vertices(g_).first; v2 != boost::vertices(g_).second; ++v2)
      {
        if (boost::get(boost::get(boost::vertex_predecessor, g_), *v2) == v1)
        {
          visual_->viz4State(getVertexState(*v2), tools::LARGE, tools::RED, 0);

          // Show state's edges
          foreach (SparseEdge edge, boost::out_edges(*v2, g_))
          {
            SparseVertex e_v1 = boost::source(edge, g_);
            SparseVertex e_v2 = boost::target(edge, g_);
            visual_->viz4Edge(getVertexState(e_v1), getVertexState(e_v2), edgeWeightProperty_[edge]);
          }
          visual_->viz4Trigger();

          // Show this robot state
          // visual_->viz4State(getVertexState(*v2), tools::ROBOT, tools::DEFAULT, 0);
          visual_->viz4State(getVertexState(*v2), tools::SMALL, tools::RED, 0);

          usleep(0.1 * 1000000);
        }
      }
    }
  }
}

std::size_t SparseDB::checkConnectedComponents()
{
  // Check how many disjoint sets are in the sparse graph (should be none)
  std::size_t numSets = getDisjointSetsCount();
  if (numSets > 1)
  {
    OMPL_ERROR("More than 1 connected component is in the sparse graph: %u", numSets);
  }

  return numSets;
}

bool SparseDB::sameComponent(SparseVertex v1, SparseVertex v2)
{
  return boost::same_component(v1, v2, disjointSets_);
}

StateID SparseDB::addState(base::State *state)
{
  return denseCache_->addState(state);
}

SparseVertex SparseDB::addVertex(base::State *state, const VertexType &type, std::size_t indent)
{
  return addVertex(addState(state), type, indent);
}

SparseVertex SparseDB::addVertex(StateID stateID, const VertexType &type, std::size_t indent)
{
  // Create vertex
  SparseVertex v = boost::add_vertex(g_);
  BOLT_CYAN_DEBUG(indent, vAdd_, "addVertex(): v: " << v << ", stateID: " << stateID << " type " << type);

  // Add properties
  vertexTypeProperty_[v] = type;
  stateCacheProperty_[v] = stateID;
  vertexPopularity_[v] = MAX_POPULARITY_WEIGHT;  // 100 means the vertex is very unpopular

  // Debug
  // std::cout << "New Vertex: " << v << " - stateID: " << stateID << " state: ";
  // debugState(getState(stateID));

  // Clear all nearby interface data whenever a new vertex is added
  if (useFourthCriteria_)
    abandonLists(denseCache_->getStateNonConst(stateID));

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
    case DISCRETIZED:
      break;
    default:
      OMPL_ERROR("Unknown VertexType type %u", type);
  }

  // Visualize
  if (visualizeSparsGraph_)
  {
    visualizeVertex(v, type);

    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }
  }

  if (visualizeVoronoiDiagramAnimated_ || (visualizeVoronoiDiagram_ && useFourthCriteria_))
    visual_->vizVoronoiDiagram();

  // Enable saving
  graphUnsaved_ = true;

  return v;
}

void SparseDB::visualizeVertex(SparseVertex v, const VertexType &type)
{
  tools::colors color;
  tools::sizes size;

  switch (type)
  {
    case COVERAGE:
      color = tools::BLACK;
      size = tools::LARGE;
      break;
    case CONNECTIVITY:
      color = tools::ORANGE;
      size = tools::LARGE;
      break;
    case INTERFACE:
      color = tools::PINK;
      size = tools::LARGE;
      break;
    case QUALITY:
      color = tools::BLUE;
      size = tools::LARGE;
      break;
    case DISCRETIZED:
      color = tools::GREEN;
      size = tools::LARGE;
      break;
    case START:
    case GOAL:
    case CARTESIAN:
    default:
      OMPL_ERROR("Unknown type");
      exit(-1);
  }

  // Show visibility region around vertex
  if (visualizeDatabaseCoverage_)
    visual_->viz2State(getVertexState(v), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseDelta_);

  // Show vertex
  visual_->viz2State(getVertexState(v), size, color, 0);
}

void SparseDB::removeVertex(SparseVertex v)
{
  // Remove from nearest neighbor
  nn_->remove(v);

  // Debug
  // debugNN();

  // Delete state from denseDB
  stateCacheProperty_[v] = 0;  // 0 means delete

  // TODO: disjointSets is now inaccurate
  // disjointSets_.remove_set(v);

  // Remove all edges to and from vertex
  boost::clear_vertex(v, g_);

  // Remove vertex
  // boost::remove_vertex(v, g_);
}

void SparseDB::debugNN()
{
  // Show contents of GNAT
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  NearestNeighborsGNAT<SparseVertex> *gnat = dynamic_cast<NearestNeighborsGNAT<SparseVertex> *>(nn_.get());
  std::cout << "GNAT: " << *gnat << std::endl;
  std::cout << std::endl;
}

SparseEdge SparseDB::addEdge(SparseVertex v1, SparseVertex v2, EdgeType type, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, vAdd_, "addEdge(): from vertex " << v1 << " to " << v2 << " type " << type);

  if (superDebug_)  // Extra checks
  {
    assert(v1 <= getNumVertices());
    assert(v2 <= getNumVertices());
    assert(v1 != v2);
    assert(!hasEdge(v1, v2));
    assert(hasEdge(v1, v2) == hasEdge(v2, v1));
    BOOST_ASSERT_MSG(getVertexState(v1) != getVertexState(v2), "States on both sides of an edge are the same");
  }

  // Create the new edge
  SparseEdge e = (boost::add_edge(v1, v2, g_)).first;

  // Weight properties
  edgeWeightProperty_[e] = distanceFunction(v1, v2);  // TODO: use this value with astar

  // Reason edge was added to spanner
  edgeTypeProperty_[e] = type;

  // Collision properties
  edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;

  // Add the edge to the incrementeal connected components datastructure
  disjointSets_.union_set(v1, v2);

  // Visualize
  if (visualizeSparsGraph_)
  {
    visual_->viz2Edge(getVertexState(v1), getVertexState(v2), convertEdgeTypeToColor(type));
    if (visualizeSparsGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz2Trigger();
      usleep(visualizeSparsGraphSpeed_ * 1000000);
    }

    // Show only the largest edges
    if (false && edgeWeightProperty_[e] > ignoreEdgesSmallerThan_)
    {
      // std::cout << "Edge distance: " << edgeWeightProperty_[e] << std::endl;
      visual_->viz4Edge(getVertexState(v1), getVertexState(v2), convertEdgeTypeToColor(type));
      visual_->viz4Trigger();
      usleep(0.001 * 1000000);
    }

    // Show each added edge for a blip
    if (false)
    {
      visual_->viz4DeleteAllMarkers();
      visual_->viz4Edge(getVertexState(v1), getVertexState(v2), convertEdgeTypeToColor(eCONNECTIVITY));
      visual_->viz4Trigger();
      usleep(0.001 * 1000000);
    }
  }

  // Enable saving
  graphUnsaved_ = true;

  return e;
}

edgeColors SparseDB::convertEdgeTypeToColor(EdgeType edgeType)
{
  switch (edgeType)
  {
    case eCONNECTIVITY:
      return eGREEN;
      break;
    case eINTERFACE:
      return eYELLOW;
      break;
    case eQUALITY:
      return eRED;
      break;
    case eCARTESIAN:
      return eORANGE;
      break;
    default:
      OMPL_ERROR("Unknown edge type");
      exit(-1);
  }
  return eORANGE;  // dummy
}

base::State *&SparseDB::getVertexStateNonConst(SparseVertex v)
{
  BOOST_ASSERT_MSG(v >= queryVertices_.size(), "Attempted to request state of query vertex using wrong function");
  return denseCache_->getStateNonConst(stateCacheProperty_[v]);
}

const base::State *SparseDB::getVertexState(SparseVertex v) const
{
  BOOST_ASSERT_MSG(v >= queryVertices_.size(), "Attempted to request state of query vertex using wrong function");
  return denseCache_->getState(stateCacheProperty_[v]);
}

const base::State *SparseDB::getState(StateID stateID) const
{
  return denseCache_->getState(stateID);
}

const StateID SparseDB::getStateID(SparseVertex v) const
{
  return stateCacheProperty_[v];
}

void SparseDB::displayDatabase(bool showVertices, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, true, "Displaying Sparse database");
  indent += 2;

  // Error check
  if (getNumVertices() == 0 || getNumEdges() == 0)
  {
    OMPL_WARN("Unable to show database because no vertices/edges available");
    return;
  }

  // Clear previous visualization
  visual_->viz2DeleteAllMarkers();

  const std::size_t MIN_FEEDBACK = 10000;
  if (visualizeDatabaseEdges_)
  {
    // Loop through each edge
    std::size_t count = 1;
    std::size_t debugFrequency = MIN_FEEDBACK;  // std::max(2, static_cast<int>(getNumEdges() / 10));
    if (getNumEdges() > MIN_FEEDBACK)
      std::cout << "Displaying sparse edges: " << std::flush;
    foreach (SparseEdge e, boost::edges(g_))
    {
      // Add edge
      SparseVertex v1 = boost::source(e, g_);
      SparseVertex v2 = boost::target(e, g_);

      // Visualize
      visual_->viz2Edge(getVertexState(v1), getVertexState(v2), convertEdgeTypeToColor(edgeTypeProperty_[e]));

      // Prevent viz cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumEdges()) * 100.0
                  << "% " << std::flush;
        visual_->viz2Trigger();
        usleep(0.01 * 1000000);
      }

      count++;
    }
    if (getNumEdges() > MIN_FEEDBACK)
      std::cout << std::endl;
  }

  if (visualizeDatabaseVertices_)
  {
    // Loop through each vertex
    std::size_t count = 1;
    std::size_t debugFrequency = MIN_FEEDBACK;  // getNumVertices() / 10;
    if (getNumVertices() > MIN_FEEDBACK)
      std::cout << "Displaying sparse vertices: " << std::flush;
    foreach (SparseVertex v, boost::vertices(g_))
    {
      // Skip query vertices
      if (v < queryVertices_.size())
        continue;

      // Skip deleted vertices
      if (stateCacheProperty_[v] == 0)
      {
        // BOLT_DEBUG(indent, true, "Skipping/not visualizing deleted vertex: " << v);
        continue;
      }

      // Check for null states
      if (!getVertexState(v))
      {
        BOLT_RED_DEBUG(indent, true, "Null vertex found: " << v);
        continue;
      }

      visualizeVertex(v, vertexTypeProperty_[v]);

      // Prevent cache from getting too big
      if (count % debugFrequency == 0)
      {
        std::cout << std::fixed << std::setprecision(0) << (static_cast<double>(count + 1) / getNumVertices()) * 100.0
                  << "% " << std::flush;
        visual_->viz2Trigger();
        usleep(0.01 * 1000000);
      }
      count++;
    }
    if (getNumVertices() > MIN_FEEDBACK)
      std::cout << std::endl;
  }

  // Publish remaining edges
  visual_->viz2Trigger();
  usleep(0.001 * 1000000);
}

double SparseDB::distanceFunction(const SparseVertex a, const SparseVertex b) const
{
  // Special case: query vertices store their states elsewhere
  if (a < numThreads_)
  {
    return si_->distance(queryState_[a], getVertexState(b));
  }
  if (b < numThreads_)
  {
    return si_->distance(getVertexState(a), queryState_[b]);
  }

  // Error check
  assert(getVertexState(a) != NULL);
  assert(getVertexState(b) != NULL);

  return si_->distance(getVertexState(a), getVertexState(b));
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
  visual_->viz6State(getVertexState(v), tools::LARGE, tools::RED, 0);

  for (auto it = iData.begin(); it != iData.end(); ++it)
  {
    const VertexPair &pair = it->first;
    InterfaceData &iData = it->second;

    SparseVertex v1 = pair.first;
    SparseVertex v2 = pair.second;

    visual_->viz6State(getVertexState(v1), tools::LARGE, tools::PURPLE, 0);
    visual_->viz6State(getVertexState(v2), tools::LARGE, tools::PURPLE, 0);
    // visual_->viz6Edge(getVertexState(v1), getVertexState(v2), tools::eGREEN);

    if (iData.hasInterface1())
    {
      visual_->viz6State(iData.getInterface1Inside(), tools::MEDIUM, tools::ORANGE, 0);
      visual_->viz6State(iData.getInterface1Outside(), tools::MEDIUM, tools::GREEN, 0);
      visual_->viz6Edge(iData.getInterface1Inside(), iData.getInterface1Outside(), tools::eYELLOW);
    }
    else
    {
      visual_->viz6Edge(getVertexState(v1), getVertexState(v2), tools::eRED);
    }

    if (iData.hasInterface2())
    {
      visual_->viz6State(iData.getInterface2Inside(), tools::MEDIUM, tools::ORANGE, 0);
      visual_->viz6State(iData.getInterface2Outside(), tools::MEDIUM, tools::GREEN, 0);
      visual_->viz6Edge(iData.getInterface2Inside(), iData.getInterface2Outside(), tools::eYELLOW);
    }
    else
    {
      visual_->viz6Edge(getVertexState(v1), getVertexState(v2), tools::eRED);
    }
  }

  visual_->viz6Trigger();
  usleep(0.01 * 1000000);
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
        visual_->viz6State(iData.getInterface1Inside(), tools::MEDIUM, tools::ORANGE, 0);
        visual_->viz6State(iData.getInterface1Outside(), tools::MEDIUM, tools::GREEN, 0);
      }

      if (iData.hasInterface2())
      {
        visual_->viz6State(iData.getInterface2Inside(), tools::MEDIUM, tools::ORANGE, 0);
        visual_->viz6State(iData.getInterface2Outside(), tools::MEDIUM, tools::GREEN, 0);
      }
    }
  }
  visual_->viz6Trigger();
  usleep(0.01 * 1000000);
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

      // Debug
      if (false)
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
  queryState_[threadID] = state;
  nn_->nearestK(queryVertices_[threadID], numNeighbors, graphNeighbors);
  queryState_[threadID] = nullptr;

  if (graphNeighbors.empty())
  {
    std::cout << "no neighbors found for sparse representative " << std::endl;
    exit(-1);
  }
  return graphNeighbors[0];
}

bool SparseDB::sufficientClearance(base::State *state)
{
  // Check if new vertex has enough clearance
  double dist = si_->getStateValidityChecker()->clearance(state);
  return dist >= obstacleClearance_;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
