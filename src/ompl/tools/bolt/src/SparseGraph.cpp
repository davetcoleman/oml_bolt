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
   Desc:   Near-asypmotically optimal roadmap datastructure
*/

// OMPL
#include <ompl/tools/bolt/SparseGraph.h>
#include <ompl/tools/bolt/SparseCriteria.h>
#include <ompl/util/Console.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <ompl/datastructures/NearestNeighborsGNATNoThreadSafety.h>
#include <ompl/base/DiscreteMotionValidator.h>

// Boost
#include <boost/graph/incremental_components.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/assert.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>

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

namespace ompl
{
namespace tools
{
namespace bolt
{
SparseGraph::SparseGraph(base::SpaceInformationPtr si, VisualizerPtr visual)
  : si_(si)
  , visual_(visual)
    // Property accessors of edges
  , edgeWeightProperty_(boost::get(boost::edge_weight, g_))
  , edgeTypeProperty_(boost::get(edge_type_t(), g_))
  , edgeCollisionStatePropertySparse_(boost::get(edge_collision_state_t(), g_))
    // Property accessors of vertices
  , vertexStateProperty_(boost::get(vertex_state_cache_t(), g_))
  , vertexTypeProperty_(boost::get(vertex_type_t(), g_))
  , vertexInterfaceProperty_(boost::get(vertex_interface_data_t(), g_))
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

  // Saving and loading from file
  sparseStorage_.reset(new SparseStorage(si_, this));

  // Initialize nearest neighbor datastructure
  nn_.reset(new NearestNeighborsGNATNoThreadSafety<SparseVertex>());
  nn_->setDistanceFunction(boost::bind(&otb::SparseGraph::distanceFunction, this, _1, _2));

  if (superDebug_)
    OMPL_WARN("Superdebug mode is enabled - will run slower");
}

SparseGraph::~SparseGraph()
{
  freeMemory();
}

void SparseGraph::freeMemory()
{
  foreach (SparseVertex v, boost::vertices(g_))
  {
    foreach (InterfaceData &iData, vertexInterfaceProperty_[v] | boost::adaptors::map_values)
      iData.clear(si_);
  }

  g_.clear();

  if (nn_)
    nn_->clear();
}

bool SparseGraph::setup()
{
  // Initialize path simplifier
  if (!pathSimplifier_)
  {
    pathSimplifier_.reset(new geometric::PathSimplifier(si_));
    pathSimplifier_->freeStates(false);
  }

  return true;
}

void SparseGraph::clearStatistics()
{
  numSamplesAddedForCoverage_ = 0;
  numSamplesAddedForConnectivity_ = 0;
  numSamplesAddedForInterface_ = 0;
  numSamplesAddedForQuality_ = 0;
}

void SparseGraph::initializeQueryState()
{
  if (boost::num_vertices(g_) > 0)
  {
    OMPL_WARN("Not initializing query state because already is of size %u", boost::num_vertices(g_));
    return;  // assume its already been setup
  }

  // Create a query state for each possible thread
  queryVertices_.resize(numThreads_);
  queryStates_.resize(numThreads_);

  for (std::size_t threadID = 0; threadID < numThreads_; ++threadID)
  {
    // Add a fake vertex to the graph
    queryVertices_[threadID] = boost::add_vertex(g_);
    vertexStateProperty_[queryVertices_[threadID]] = 0;  // set stateID to the zeroth stateID - NULL
  }
}

bool SparseGraph::load()
{
  // Load collision cache
  denseCache_->load();

  // Benchmark
  time::point start = time::now();

  if (!sparseStorage_->load(filePath_.c_str()))
    return false;

  // Benchmark
  double duration = time::seconds(time::now() - start);
  OMPL_INFORM("Graph total loading time: %f", duration);

  // Error check
  if (!getNumVertices() || !getNumEdges())
  {
    OMPL_ERROR("Corrupted sparse graph loaded");
    return false;
  }

  // Show more data
  printGraphStats();

  // Nothing to save because was just loaded from file
  graphUnsaved_ = false;

  if (visualizeGraphAfterLoading_)
    displayDatabase();

  return true;
}

bool SparseGraph::saveIfChanged(std::size_t indent)
{
  if (graphUnsaved_)
  {
    return save(indent);
  }
  else
  {
    BOLT_DEBUG(indent, true, "Not saving because database has not changed. Time: " << time::as_string(time::now()));
  }

  return true;
}

bool SparseGraph::save(std::size_t indent)
{
  if (!graphUnsaved_)
    OMPL_WARN("No need to save because graphUnsaved_ is false, but saving anyway because requested");

  // Disabled
  if (!savingEnabled_)
  {
    OMPL_INFORM("Not saving because option disabled for SparseGraph");
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
  sparseStorage_->save(filePath_.c_str());

  // Save collision cache
  denseCache_->save();

  // Benchmark
  double loadTime = time::seconds(time::now() - start);
  BOLT_DEBUG(indent, true, "Saved database to file in " << loadTime << " seconds. Time: " << time::as_string(time::now()));

  graphUnsaved_ = false;
  return true;
}

bool SparseGraph::astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex> &vertexPath,
                              double &distance, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, vSearch_, "astarSearch()");
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
    visual_->viz4()->deleteAllMarkers();
  }

  try
  {
    double popularityBias = 0;
    bool popularityBiasEnabled = false;
    boost::astar_search(g_,                                                              // graph
                        start,                                                           // start state
                        boost::bind(&otb::SparseGraph::astarHeuristic, this, _1, goal),  // the heuristic
                        // ability to disable edges (set cost to inifinity):
                        boost::weight_map(SparseEdgeWeightMap(g_, edgeCollisionStatePropertySparse_, popularityBias,
                                                              popularityBiasEnabled))
                        .predecessor_map(vertexPredecessors)
                        .distance_map(&vertexDistances[0])
                        .visitor(SparsestarVisitor(goal, this)));
  }
  catch (FoundGoalException &)
  {
    distance = vertexDistances[goal];

    // the custom exception from SparsestarVisitor
    BOLT_DEBUG(indent, vSearch_, "AStar found solution. Distance to goal: " << vertexDistances[goal]);
    BOLT_DEBUG(indent, vSearch_, "Number nodes opened: " << numNodesOpened_
               << ", Number nodes closed: " << numNodesClosed_);

    if (isinf(vertexDistances[goal]))  // TODO(davetcoleman): test that this works
    {
      throw Exception(name_, "Distance to goal is infinity");
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
      if (v != goal)
      {
        vertexPath.push_back(v);
      }

      foundGoal = true;
    }
  }

  if (!foundGoal)
    BOLT_YELLOW_DEBUG(indent, vSearch_, "Did not find goal");

  // Show all predecessors
  if (visualizeAstar_)
  {
    BOLT_DEBUG(indent, vSearch_, "Show all predecessors");
    for (std::size_t i = numThreads_; i < getNumVertices(); ++i)  // skip vertex 0-11 because those are query vertices
    {
      const SparseVertex v1 = i;
      const SparseVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4()->edge(getVertexState(v1), getVertexState(v2), 10);
      }
    }
    visual_->viz4()->trigger();
  }

  // Unload
  delete[] vertexPredecessors;
  delete[] vertexDistances;

  // No solution found from start to goal
  return foundGoal;
}

double SparseGraph::astarHeuristic(const SparseVertex a, const SparseVertex b) const
{
  // Assume vertex 'a' is the one we care about its populariy

  // Get the classic distance
  double dist = si_->distance(getVertexState(a), getVertexState(b));

  // if (false)  // method 1
  // {
  //   const double percentMaxExtent = (maxExtent_ * percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
  //   double popularityComponent = percentMaxExtent * (vertexPopularity_[a] / 100.0);

  //   std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " <<
  //   vertexPopularity_[a]
  //             << ", max extent: " << maxExtent_ << ", percentMaxExtent: " << percentMaxExtent
  //             << ", popularityComponent: " << popularityComponent;
  //   dist = std::max(0.0, dist - popularityComponent);
  // }
  // else if (false)  // method 2
  // {
  //   const double percentDist = (dist * percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
  //   double popularityComponent = percentDist * (vertexPopularity_[a] / 100.0);

  //   std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " <<
  //   vertexPopularity_[a]
  //             << ", percentDist: " << percentDist << ", popularityComponent: " << popularityComponent;
  //   dist = std::max(0.0, dist - popularityComponent);
  // }
  // else if (false)  // method 3
  // {
  //   std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " <<
  //   vertexPopularity_[a]
  //             << ", vertexPopularity_[a] / 100.0: " << vertexPopularity_[a] / 100.0
  //             << ", percentMaxExtentUnderestimate_: " << percentMaxExtentUnderestimate_;
  //   // if ((vertexPopularity_[a] / 100.0) < (1 - percentMaxExtentUnderestimate_))
  //   if (vertexPopularity_[a] > (100 - percentMaxExtentUnderestimate_ * 100.0))
  //   {
  //     dist = 0;
  //   }

  //   // dist = std::max(0.0, dist - popularityComponent);
  // }
  // else if (false)  // method 4
  // {
  //   dist *= (1 + percentMaxExtentUnderestimate_);
  // }
  // method 5: increasing the sparseDelta fraction

  // std::cout << ", new distance: " << dist << std::endl;

  return dist;
}

double SparseGraph::distanceFunction(const SparseVertex a, const SparseVertex b) const
{
  //std::cout << "sg.distancefunction() " << a << ", " << b << std::endl;
  // Special case: query vertices store their states elsewhere
  if (a < numThreads_)
  {
    return si_->distance(queryStates_[a], getVertexState(b));
  }
  if (b < numThreads_)
  {
    return si_->distance(getVertexState(a), queryStates_[b]);
  }

  // Error check
  if (superDebug_)
  {
    assert(getVertexState(a) != NULL);
    assert(getVertexState(b) != NULL);
  }

  return si_->distance(getVertexState(a), getVertexState(b));
}

bool SparseGraph::isEmpty() const
{
  assert(!(getNumVertices() < getNumQueryVertices()));
  return (getNumVertices() == getNumQueryVertices() && getNumEdges() == 0);
}

void SparseGraph::clearEdgeCollisionStates()
{
  foreach (const SparseEdge e, boost::edges(g_))
    edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;  // each edge has an unknown state
}

void SparseGraph::errorCheckDuplicateStates(std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, true, "errorCheckDuplicateStates() - part of super debug");
  indent += 2;

  bool found = false;
  // Error checking: check for any duplicate states
  for (std::size_t i = 0; i < denseCache_->getStateCacheSize(); ++i)
  {
    for (std::size_t j = i + 1; j < denseCache_->getStateCacheSize(); ++j)
    {
      if (si_->getStateSpace()->equalStates(getState(i), getState(j)))
      {
        BOLT_RED_DEBUG(indent, 1, "Found equal state: " << i << ", " << j);
        debugState(getState(i));
        found = true;
      }
    }
  }
  if (found)
    throw Exception(name_, "Duplicate state found");
}

bool SparseGraph::smoothQualityPathOriginal(geometric::PathGeometric *path, std::size_t indent)
{
  BOLT_RED_DEBUG(indent, visualizeQualityPathSimp_, "smoothQualityPathOriginal()");
  indent += 2;

  // Visualize path
  if (visualizeQualityPathSimp_)
  {
    visual_->viz2()->deleteAllMarkers();
    visual_->viz2()->path(path, tools::SMALL, tools::BLUE);
    visual_->viz2()->trigger();
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
    throw Exception(name_, "check and repair failed?");
  }
  return true;
}

bool SparseGraph::smoothQualityPath(geometric::PathGeometric *path, double clearance, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, visualizeQualityPathSimp_, "smoothQualityPath()");
  indent += 2;

  // Visualize path
  if (visualizeQualityPathSimp_)
  {
    visual_->viz2()->deleteAllMarkers();
    visual_->viz2()->path(path, tools::SMALL, tools::BLUE);
    visual_->viz2()->trigger();
    usleep(0.001 * 1000000);
  }

  BOLT_DEBUG(indent, visualizeQualityPathSimp_, "Created 'quality path' candidate with " << path->getStateCount()
             << " states");
  if (visualizeQualityPathSimp_)
    visual_->waitForUserFeedback("path simplification");

  // Set the motion validator to use clearance, this way isValid() checks clearance before confirming valid
  base::DiscreteMotionValidator *dmv =
    dynamic_cast<base::DiscreteMotionValidator *>(si_->getMotionValidatorNonConst().get());
  dmv->setRequiredStateClearance(clearance);

  for (std::size_t i = 0; i < 3; ++i)
  {
    pathSimplifier_->simplifyMax(*path);

    if (visualizeQualityPathSimp_)
    {
      visual_->viz2()->deleteAllMarkers();
      visual_->viz2()->path(path, tools::SMALL, tools::ORANGE);
      visual_->viz2()->trigger();
      usleep(0.1 * 1000000);
      // visual_->waitForUserFeedback("optimizing path");
    }

    pathSimplifier_->reduceVertices(*path, 1000, path->getStateCount() * 4);

    if (visualizeQualityPathSimp_)
    {
      visual_->viz2()->deleteAllMarkers();
      visual_->viz2()->path(path, tools::SMALL, tools::BLUE);
      visual_->viz2()->trigger();
      usleep(0.1 * 1000000);
      // visual_->waitForUserFeedback("optimizing path");
    }
  }
  // Turn off the clearance requirement
  dmv->setRequiredStateClearance(0.0);

  pathSimplifier_->reduceVertices(*path, 1000, path->getStateCount() * 4);

  if (visualizeQualityPathSimp_)
  {
    visual_->viz2()->deleteAllMarkers();
    visual_->viz2()->path(path, tools::SMALL, tools::GREEN);
    visual_->viz2()->trigger();
    visual_->waitForUserFeedback("finished quality path");
  }

  std::pair<bool, bool> repairResult = path->checkAndRepair(100);

  if (!repairResult.second)  // Repairing was not successful
  {
    throw Exception(name_, "check and repair failed?");
  }
  return true;
}

std::size_t SparseGraph::getDisjointSetsCount(bool verbose)
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

void SparseGraph::getDisjointSets(SparseDisjointSetsMap &disjointSets)
{
  disjointSets.clear();

  // Flatten the parents tree so that the parent of every element is its representative.
  disjointSets_.compress_sets(boost::vertices(g_).first, boost::vertices(g_).second);

  // Count size of each disjoint set and group its containing vertices
  typedef boost::graph_traits<SparseAdjList>::vertex_iterator VertexIterator;
  for (VertexIterator v = boost::vertices(g_).first; v != boost::vertices(g_).second; ++v)
  {
    // Do not count the search vertex within the sets
    if (*v <= queryVertices_.back())
      continue;

    disjointSets[boost::get(boost::get(boost::vertex_predecessor, g_), *v)].push_back(*v);
  }
}

void SparseGraph::printDisjointSets(SparseDisjointSetsMap &disjointSets)
{
  OMPL_INFORM("Print disjoint sets");
  for (SparseDisjointSetsMap::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const SparseVertex v = iterator->first;
    const std::size_t freq = iterator->second.size();
    std::cout << "  Parent: " << v << " frequency " << freq << std::endl;
  }
}

void SparseGraph::visualizeDisjointSets(SparseDisjointSetsMap &disjointSets)
{
  OMPL_INFORM("Visualizing disjoint sets");

  // Find the disjoint set that is the 'main' large one
  std::size_t maxDisjointSetSize = 0;
  SparseVertex maxDisjointSetParent = 0;
  for (SparseDisjointSetsMap::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
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
  for (SparseDisjointSetsMap::const_iterator iterator = disjointSets.begin(); iterator != disjointSets.end();
       iterator++)
  {
    const SparseVertex v1 = iterator->first;
    const std::size_t freq = iterator->second.size();

    // std::cout << std::endl;
    // std::cout << "Parent vertex: " << v1 << " StateID: " << getStateID(v1) << " Frequency: " << freq << std::endl;
    // debugState(getVertexState(v1));

    BOOST_ASSERT_MSG(freq > 0, "Frequency must be at least 1");

    if (freq == maxDisjointSetSize)  // any subgraph that is smaller than the full graph
      continue;                      // the main disjoint set is not considered a disjoint set

    // Visualize sets of size one
    if (freq == 1)
    {
      visual_->viz5()->state(getVertexState(v1), tools::LARGE, tools::RED, 0);
      visual_->viz5()->trigger();
      visual_->waitForUserFeedback("showing disjoint set");
      continue;
    }

    // Visualize large disjoint sets (greater than one)
    if (freq > 1 && freq < 1000)
    {
      // Clear markers
      visual_->viz4()->deleteAllMarkers();

      // Visualize this subgraph that is disconnected
      // Loop through every every vertex and check if its part of this group
      typedef boost::graph_traits<SparseAdjList>::vertex_iterator VertexIterator;
      for (VertexIterator v2 = boost::vertices(g_).first; v2 != boost::vertices(g_).second; ++v2)
      {
        if (boost::get(boost::get(boost::vertex_predecessor, g_), *v2) == v1)
        {
          visual_->viz4()->state(getVertexState(*v2), tools::LARGE, tools::RED, 0);

          // Show state's edges
          foreach (SparseEdge edge, boost::out_edges(*v2, g_))
          {
            SparseVertex e_v1 = boost::source(edge, g_);
            SparseVertex e_v2 = boost::target(edge, g_);
            visual_->viz4()->edge(getVertexState(e_v1), getVertexState(e_v2), edgeWeightProperty_[edge]);
          }
          visual_->viz4()->trigger();

          // Show this robot state
          // visual_->viz4()->state(getVertexState(*v2), tools::ROBOT, tools::DEFAULT, 0);
          visual_->viz4()->state(getVertexState(*v2), tools::SMALL, tools::RED, 0);

          usleep(0.1 * 1000000);
        }  // if
      }    // for

      visual_->waitForUserFeedback("showing large disjoint set");
    }  // if
  }
}

std::size_t SparseGraph::checkConnectedComponents()
{
  // Check how many disjoint sets are in the sparse graph (should be none)
  std::size_t numSets = getDisjointSetsCount();
  if (numSets > 1)
  {
    std::size_t indent = 0;
    BOLT_YELLOW_DEBUG(indent, true, "More than 1 connected component is in the sparse graph: " << numSets);
  }

  return numSets;
}

bool SparseGraph::sameComponent(SparseVertex v1, SparseVertex v2)
{
  return boost::same_component(v1, v2, disjointSets_);
}

StateID SparseGraph::addState(base::State *state)
{
  return denseCache_->addState(state);
}

SparseVertex SparseGraph::addVertex(base::State *state, const VertexType &type, std::size_t indent)
{
  return addVertex(addState(state), type, indent);
}

SparseVertex SparseGraph::addVertex(StateID stateID, const VertexType &type, std::size_t indent)
{
  // Create vertex
  SparseVertex v = boost::add_vertex(g_);

  // Feedback
  BOLT_CYAN_DEBUG(indent, vAdd_, "addVertex(): v: " << v << ", stateID: " << stateID << " type " << type);
  indent += 2;

  // Add properties
  vertexTypeProperty_[v] = type;
  vertexStateProperty_[v] = stateID;
  vertexPopularity_[v] = MAX_POPULARITY_WEIGHT;  // 100 means the vertex is very unpopular

  // Clear all nearby interface data whenever a new vertex is added
  if (sparseCriteria_->useFourthCriteria_)
    clearInterfaceData(denseCache_->getStateNonConst(stateID));

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
  if (visualizeSparseGraph_)
  {
    visualizeVertex(v, type);

    if (visualizeSparseGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz1()->trigger();

      if (visualizeProjection_) // Hack: Project to 2D space
        visual_->viz7()->trigger();

      usleep(visualizeSparseGraphSpeed_ * 1000000);
    }
  }

  // Optional Voronoi Diagram
  if (sparseCriteria_->visualizeVoronoiDiagramAnimated_ ||
      (sparseCriteria_->visualizeVoronoiDiagram_ && sparseCriteria_->useFourthCriteria_))
    visual_->vizVoronoiDiagram();

  // Enable saving
  graphUnsaved_ = true;

  // Debugging
  if (!sparseCriteria_->getDiscretizedSamplesInsertion())
    visual_->waitForUserFeedback("Added vertex randomly");

  return v;
}

SparseVertex SparseGraph::addVertexFromFile(base::State *state, const VertexType &type, std::size_t indent)
{
  // Save state in cache
  StateID stateID = addState(state);

  // Create vertex
  SparseVertex v = boost::add_vertex(g_);

  // Add properties
  vertexTypeProperty_[v] = type;
  vertexStateProperty_[v] = stateID;
  vertexPopularity_[v] = MAX_POPULARITY_WEIGHT;  // 100 means the vertex is very unpopular

  // Connected component tracking
  disjointSets_.make_set(v);

  return v;
}

void SparseGraph::removeVertex(SparseVertex v)
{
  // Remove from nearest neighbor
  nn_->remove(v);

  // Delete state from denseDB
  vertexStateProperty_[v] = 0;  // 0 means delete

  // TODO: disjointSets is now inaccurate
  // Our checkAddConnectivity() criteria is broken
  // because we frequntly delete edges and nodes..
  // disjointSets_.remove_set(v);

  // Remove all edges to and from vertex
  boost::clear_vertex(v, g_);

  // We do not actually remove the vertex from the graph
  // because that would invalidate the nearest neighbor tree
}

void SparseGraph::removeDeletedVertices(std::size_t indent)
{
  bool verbose = true;
  BOLT_RED_DEBUG(indent, verbose || true, "removeDeletedVertices()");
  indent += 2;

  // Remove all vertices that are set to 0
  std::size_t numRemoved = 0;

  // Iterate manually through graph
  typedef boost::graph_traits<SparseAdjList>::vertex_iterator VertexIterator;
  for (VertexIterator v = boost::vertices(g_).first; v != boost::vertices(g_).second; /* manual */)
  {
    if (*v < numThreads_)  // Skip query vertices
    {
      v++;
      continue;
    }

    if (getStateID(*v) == 0)  // Found vertex to delete
    {
      BOLT_DEBUG(indent, verbose, "Removing SparseVertex " << *v << " stateID: " << getStateID(*v));

      boost::remove_vertex(*v, g_);
      numRemoved++;
    }
    else  // only proceed if no deletion happened
    {
      // BOLT_DEBUG(indent, verbose, "Checking SparseVertex " << *v << " stateID: " << getStateID(*v));
      v++;
    }
  }
  BOLT_DEBUG(indent, verbose, "Removed " << numRemoved << " vertices from graph that were abandoned");

  if (numRemoved == 0)
  {
    BOLT_DEBUG(indent, verbose, "No verticies deleted, skipping resetting NN and disjointSets");
    return;
  }

  // Reset the nearest neighbor tree
  nn_->clear();

  // Reset disjoint sets
  disjointSets_ = SparseDisjointSetType(boost::get(boost::vertex_rank, g_), boost::get(boost::vertex_predecessor, g_));

  // Reinsert vertices into nearest neighbor
  foreach (SparseVertex v, boost::vertices(g_))
  {
    if (v < numThreads_)  // Skip the query vertices
      continue;

    nn_->add(v);
    disjointSets_.make_set(v);
  }

  // Reinsert edges into disjoint sets
  foreach (SparseEdge e, boost::edges(g_))
  {
    SparseVertex v1 = boost::source(e, g_);
    SparseVertex v2 = boost::target(e, g_);
    disjointSets_.union_set(v1, v2);
  }
}

SparseEdge SparseGraph::addEdge(SparseVertex v1, SparseVertex v2, EdgeType type, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, vAdd_, "addEdge(): from vertex " << v1 << " to " << v2 << " type " << type);
  indent += 2;

  if (superDebug_)  // Extra checks
  {
    BOOST_ASSERT_MSG(v1 <= getNumVertices(), "Vertex1 is larger than max vertex id");
    BOOST_ASSERT_MSG(v2 <= getNumVertices(), "Vertex2 is larger than max vertex id");
    BOOST_ASSERT_MSG(v1 != v2, "Verticex IDs are the same");
    BOOST_ASSERT_MSG(!hasEdge(v1, v2), "There already exists an edge between two vertices requested");
    BOOST_ASSERT_MSG(hasEdge(v1, v2) == hasEdge(v2, v1), "There already exists an edge between two vertices requested, "
                     "other direction");
    BOOST_ASSERT_MSG(getVertexState(v1) != getVertexState(v2), "States on both sides of an edge are the same");
    BOOST_ASSERT_MSG(!si_->getStateSpace()->equalStates(getVertexState(v1), getVertexState(v2)),
                     "Vertex IDs are different but states are the equal");
  }

  // Create the new edge
  SparseEdge e = (boost::add_edge(v1, v2, g_)).first;

  // Weight properties
  edgeWeightProperty_[e] = distanceFunction(v1, v2);

  // Reason edge was added to spanner in SPARS
  edgeTypeProperty_[e] = type;

  // Collision properties
  edgeCollisionStatePropertySparse_[e] = NOT_CHECKED;

  // Add the edge to the incrementeal connected components datastructure
  disjointSets_.union_set(v1, v2);

  // Visualize
  if (visualizeSparseGraph_)
  {
    visualizeEdge(v1, v2, type, /*windowID*/ 1);
    visualizeEdge(v1, v2, type, /*windowID*/ 7); // projection to 2D space

    if (visualizeSparseGraphSpeed_ > std::numeric_limits<double>::epsilon())
    {
      visual_->viz1()->trigger();

      if (visualizeProjection_) // Hack: Project to 2D space
        visual_->viz7()->trigger();

      usleep(visualizeSparseGraphSpeed_ * 1000000);
    }
  }

  // Enable saving
  graphUnsaved_ = true;

  // Debugging
  // if (!sparseCriteria_->getDiscretizedSamplesInsertion())
  //   visual_->waitForUserFeedback("Added edge randomly");

  return e;
}

bool SparseGraph::hasEdge(SparseVertex v1, SparseVertex v2)
{
  return boost::edge(v1, v2, g_).second;
}

VizColors SparseGraph::edgeTypeToColor(EdgeType edgeType)
{
  switch (edgeType)
  {
    case eCONNECTIVITY:
      return ORANGE;
      break;
    case eINTERFACE:
      return YELLOW;
      break;
    case eQUALITY:
      return RED;
      break;
    case eCARTESIAN:
      return MAGENTA;
      break;
    case eDISCRETIZED:
      return BLUE;
      break;
    default:
      throw Exception(name_, "Unknown edge type");
  }
  return ORANGE;  // dummy return value
}

base::State *&SparseGraph::getQueryStateNonConst(SparseVertex v)
{
  BOOST_ASSERT_MSG(v < queryVertices_.size(), "Attempted to request state of regular vertex using query function");
  return queryStates_[v];
}

base::State *&SparseGraph::getVertexStateNonConst(SparseVertex v)
{
  BOOST_ASSERT_MSG(v >= queryVertices_.size(), "Attempted to request state of query vertex using wrong function");
  return denseCache_->getStateNonConst(vertexStateProperty_[v]);
}

const base::State *SparseGraph::getVertexState(SparseVertex v) const
{
  BOOST_ASSERT_MSG(v >= queryVertices_.size(), "Attempted to request state of query vertex using wrong function");
  return denseCache_->getState(vertexStateProperty_[v]);
}

const base::State *SparseGraph::getState(StateID stateID) const
{
  return denseCache_->getState(stateID);
}

const StateID SparseGraph::getStateID(SparseVertex v) const
{
  return vertexStateProperty_[v];
}

SparseVertex SparseGraph::getSparseRepresentative(base::State *state)
{
  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;
  const std::size_t numNeighbors = 1;

  // Search for nearest sparse vertex of the provided state - this vertex provides its coverage
  queryStates_[threadID] = state;
  nn_->nearestK(queryVertices_[threadID], numNeighbors, graphNeighbors);
  queryStates_[threadID] = nullptr;

  if (graphNeighbors.empty())
  {
    throw Exception(name_, "No neighbors found for sparse representative");
  }
  return graphNeighbors[0];
}

void SparseGraph::clearInterfaceData(base::State *state)
{
  std::vector<SparseVertex> graphNeighbors;
  const std::size_t threadID = 0;

  // Search
  queryStates_[threadID] = state;
  nn_->nearestR(queryVertices_[threadID], 2.0 * sparseCriteria_->sparseDelta_, graphNeighbors);
  queryStates_[threadID] = nullptr;

  // For each of the vertices
  foreach (SparseVertex v, graphNeighbors)
  {
    foreach (VertexPair r, vertexInterfaceProperty_[v] | boost::adaptors::map_keys)
    {
      vertexInterfaceProperty_[v][r].clear(si_);
    }
  }
}

void SparseGraph::clearEdgesNearVertex(SparseVertex vertex, std::size_t indent)
{
  BOLT_RED_DEBUG(indent, true, "clearEdgesNearVertex()");
  indent += 2;

  // Optionally disable this feature
  if (!sparseCriteria_->useClearEdgesNearVertex_)
    return;

  // TODO(davetcoleman): combine this with clearInterfaceData and ensure that all interface data is equally cleared
  // but do not clear out nearby edges if a non-quality-path vertex is added
  std::vector<SparseVertex> graphNeighbors;

  // Search
  nn_->nearestR(vertex, sparseCriteria_->sparseDelta_, graphNeighbors);

  // For each of the vertices
  foreach (SparseVertex v, graphNeighbors)
  {
    // Remove all edges to and from vertex
    boost::clear_vertex(v, g_);
  }

  // Only display database if enabled
  if (visualizeSparseGraph_ && visualizeSparseGraphSpeed_ > std::numeric_limits<double>::epsilon())
    displayDatabase(true, indent);
}

void SparseGraph::displayDatabase(bool showVertices, std::size_t indent)
{
  BOLT_RED_DEBUG(indent, vVisualize_ || true, "Displaying Sparse database");
  indent += 2;

  // Error check
  if (getNumVertices() == 0 || getNumEdges() == 0)
  {
    OMPL_WARN("Unable to show database because no vertices/edges available");
    return;
  }

  // Clear previous visualization
  //visual_->viz1()->deleteAllMarkers();
  //visual_->viz1()->trigger();
  if (visualizeProjection_) // Hack: Project to 2D space
  {
    visual_->viz7()->deleteAllMarkers();
  }

  // Edges
  if (visualizeDatabaseEdges_)
  {
    // Loop through each edge
    foreach (SparseEdge e, boost::edges(g_))
    {
      // Get edge endpoints
      SparseVertex v1 = boost::source(e, g_);
      SparseVertex v2 = boost::target(e, g_);

      visualizeEdge(v1, v2, edgeTypeProperty_[e], /*windowID*/ 7); // projection to 2D space
    }
  }

  // Vertices
  if (visualizeDatabaseVertices_)
  {
    std::vector<const ompl::base::State*> states;
    std::vector<ot::VizColors> colors;

    // Loop through each vertex
    foreach (SparseVertex v, boost::vertices(g_))
    {
      // Skip query vertices
      if (v < queryVertices_.size())
        continue;

      // Skip deleted vertices
      if (vertexStateProperty_[v] == 0)
      {
        continue;
      }

      // Check for null states
      if (!getVertexState(v))
      {
        BOLT_RED_DEBUG(indent, true, "Null vertex found: " << v);
        continue;
      }

      // Populate properties
      colors.push_back(vertexTypeToColor(vertexTypeProperty_[v]));
      states.push_back(getVertexState(v));
    }

    // Create marker and push to queue
    //visual_->viz1()->states(states, colors, vertexSize_);
    visual_->viz7()->states(states, colors, vertexSize_);
  }

  // Publish remaining edges
  //visual_->viz1()->trigger();

  if (visualizeProjection_) // Hack: Project to 2D space
  {
    visual_->viz7()->trigger();
  }

  usleep(0.001 * 1000000);
}

void SparseGraph::visualizeVertex(SparseVertex v, const VertexType &type)
{
  tools::VizColors color = vertexTypeToColor(type);

  // Show visibility region around vertex
  if (visualizeDatabaseCoverage_)
    visual_->viz1()->state(getVertexState(v), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT,
                           sparseCriteria_->sparseDelta_);

  // Show vertex
  visual_->viz1()->state(getVertexState(v), vertexSize_, color, 0);

  if (visualizeProjection_) // Hack: Project to 2D space
    visual_->viz7()->state(getVertexState(v), vertexSize_, color, 0);
  //visual_->viz7()->state(getVertexState(v), tools::VARIABLE_SIZE, tools::TRANSLUCENT_LIGHT, sparseCriteria_->sparseDelta_);
}

tools::VizColors SparseGraph::vertexTypeToColor(VertexType type)
{
  switch (type)
  {
    case COVERAGE:
      return tools::GREEN;
      break;
    case CONNECTIVITY:
      return tools::BROWN;
      break;
    case INTERFACE:
      return tools::WHITE;
      break;
    case QUALITY:
      return tools::PINK;
      break;
    case CARTESIAN:
      return tools::PURPLE;
      break;
    case DISCRETIZED:
      return tools::CYAN;
      break;
    case START:
    case GOAL:
    default:
      throw Exception(name_, "Unknown type");
  }
  return tools::BLACK; // dummy return value, not actually ysed
}

void SparseGraph::visualizeEdge(SparseEdge e, EdgeType type, std::size_t windowID)
{
  // Add edge
  SparseVertex v1 = boost::source(e, g_);
  SparseVertex v2 = boost::target(e, g_);

  visualizeEdge(v1, v2, type, windowID);
}

void SparseGraph::visualizeEdge(SparseVertex v1, SparseVertex v2, EdgeType type, std::size_t windowID)
{
  // Visualize
  visual_->viz(windowID)->edge(getVertexState(v1), getVertexState(v2), edgeSize_, edgeTypeToColor(type));
}

VertexPair SparseGraph::interfaceDataIndex(SparseVertex vp, SparseVertex vpp)
{
  if (vp < vpp)
    return VertexPair(vp, vpp);
  else if (vpp < vp)
    return VertexPair(vpp, vp);

  throw Exception(name_, "Trying to get an index where the pairs are the same point!");
  return VertexPair(0, 0);  // prevent compiler warnings
}

InterfaceData &SparseGraph::getInterfaceData(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent)
{
  // BOLT_CYAN_DEBUG(indent, sparseCriteria_->vQuality_, "getInterfaceData() " << v << ", " << vp << ", " << vpp);
  return vertexInterfaceProperty_[v][interfaceDataIndex(vp, vpp)];
}

void SparseGraph::debugState(const ompl::base::State *state)
{
  si_->printState(state, std::cout);
}

void SparseGraph::debugVertex(const SparseVertex v)
{
  debugState(getVertexState(v));
}

void SparseGraph::debugNN()
{
  // Show contents of GNAT
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  NearestNeighborsGNAT<SparseVertex> *gnat = dynamic_cast<NearestNeighborsGNAT<SparseVertex> *>(nn_.get());
  std::cout << "GNAT: " << *gnat << std::endl;
  std::cout << std::endl;
}

void SparseGraph::printGraphStats()
{
  // Get the average vertex degree (number of connected edges)
  double averageDegree = (getNumEdges() * 2) / static_cast<double>(getNumVertices());

  // Check how many disjoint sets are in the sparse graph (should be none)
  std::size_t numSets = checkConnectedComponents();

  // Find min, max, and average edge length
  double totalEdgeLength = 0;
  double maxEdgeLength = -1 * std::numeric_limits<double>::infinity();
  double minEdgeLength = std::numeric_limits<double>::infinity();
  foreach (const SparseEdge e, boost::edges(g_))
  {
    const double length = edgeWeightProperty_[e];
    totalEdgeLength += length;
    if (maxEdgeLength < length)
      maxEdgeLength = length;
    if (minEdgeLength > length)
      minEdgeLength = length;
  }
  double averageEdgeLength = getNumEdges() ? totalEdgeLength / getNumEdges() : 0;

  std::size_t indent = 0;
  BOLT_DEBUG(indent, 1, "------------------------------------------------------");
  BOLT_DEBUG(indent, 1, "SparseGraph stats:");
  BOLT_DEBUG(indent, 1, "   Total vertices:         " << getNumVertices());
  BOLT_DEBUG(indent, 1, "   Total edges:            " << getNumEdges());
  BOLT_DEBUG(indent, 1, "   Average degree:         " << averageDegree);
  BOLT_DEBUG(indent, 1, "   Connected Components:   " << numSets);
  BOLT_DEBUG(indent, 1, "   Edge Lengths:           ");
  BOLT_DEBUG(indent, 1, "      Max:                 " << maxEdgeLength);
  BOLT_DEBUG(indent, 1, "      Min:                 " << minEdgeLength);
  BOLT_DEBUG(indent, 1, "      Average:             " << averageEdgeLength);
  BOLT_DEBUG(indent, 1, "      SparseDelta:         " << sparseCriteria_->sparseDelta_);
  BOLT_DEBUG(indent, 1, "      Difference:          " << averageEdgeLength - sparseCriteria_->sparseDelta_);
  BOLT_DEBUG(indent, 1, "      Penetration:         " << sparseCriteria_->discretizePenetrationDist_);
  BOLT_DEBUG(indent, 1, "------------------------------------------------------");
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

// SparseEdgeWeightMap methods ////////////////////////////////////////////////////////////////////////////

namespace boost
{
double get(const ompl::tools::bolt::SparseEdgeWeightMap &m, const ompl::tools::bolt::SparseEdge &e)
{
  return m.get(e);
}
}

BOOST_CONCEPT_ASSERT(
                     (boost::ReadablePropertyMapConcept<ompl::tools::bolt::SparseEdgeWeightMap, ompl::tools::bolt::SparseEdge>));

// SparsestarVisitor methods ////////////////////////////////////////////////////////////////////////////

BOOST_CONCEPT_ASSERT((boost::AStarVisitorConcept<otb::SparsestarVisitor, otb::SparseAdjList>));

otb::SparsestarVisitor::SparsestarVisitor(SparseVertex goal, SparseGraph *parent) : goal_(goal), parent_(parent)
{
}

void otb::SparsestarVisitor::discover_vertex(SparseVertex v, const SparseAdjList &) const
{
  // Statistics
  parent_->recordNodeOpened();

  if (parent_->visualizeAstar_)
    parent_->getVisual()->viz4()->state(parent_->getVertexState(v), tools::SMALL, tools::GREEN, 1);
}

void otb::SparsestarVisitor::examine_vertex(SparseVertex v, const SparseAdjList &) const
{
  // Statistics
  parent_->recordNodeClosed();

  if (parent_->visualizeAstar_)
  {
    parent_->getVisual()->viz4()->state(parent_->getVertexState(v), tools::LARGE, tools::BLACK, 1);
    parent_->getVisual()->viz4()->trigger();
    usleep(parent_->visualizeAstarSpeed_ * 1000000);
  }

  if (v == goal_)
    throw FoundGoalException();
}
