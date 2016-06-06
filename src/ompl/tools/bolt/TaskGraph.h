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

#ifndef OMPL_TOOLS_BOLT_TASK_GRAPH_
#define OMPL_TOOLS_BOLT_TASK_GRAPH_

#include <ompl/base/StateSpace.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/PathSimplifier.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/datastructures/NearestNeighbors.h>

// Bolt
#include <ompl/tools/bolt/SparseGraph.h>
#include <ompl/tools/bolt/BoostGraphHeaders.h>
#include <ompl/tools/bolt/DenseCache.h>
#include <ompl/tools/bolt/Debug.h>
#include <ompl/tools/bolt/VertexDiscretizer.h>
#include <ompl/tools/bolt/BoltStorage.h>
#include <ompl/tools/debug/Visualizer.h>

// Boost
#include <boost/function.hpp>
#include <boost/graph/astar_search.hpp>

// C++
#include <list>
#include <random>

namespace ompl
{
namespace tools
{
namespace bolt
{
/**
   @anchor TaskGraph
   @par Near-asypmotically optimal roadmap datastructure
*/

/// @cond IGNORE
OMPL_CLASS_FORWARD(TaskGraph);
/// @endcond

/** \class ompl::tools::bolt::::TaskGraphPtr
    \brief A boost shared pointer wrapper for ompl::tools::TaskGraph */

/** \brief Near-asypmotically optimal roadmap datastructure */
class TaskGraph
{
  friend class BoltRetrieveRepair;

public:
  /** \brief Constructor needs the state space used for planning.
   */
  TaskGraph(SparseGraphPtr sg);

  /** \brief Deconstructor */
  virtual ~TaskGraph(void);

  /* ---------------------------------------------------------------------------------
   * Setup and cleanup
   * --------------------------------------------------------------------------------- */

  /** \brief Retrieve the computed roadmap. */
  const TaskAdjList& getGraph() const
  {
    return g_;
  }

  TaskAdjList getGraphNonConst()
  {
    return g_;
  }

  DenseCachePtr getDenseCache()
  {
    return denseCache_;
  }

  base::SpaceInformationPtr getSpaceInformation()
  {
    return si_;
  }

  /** \brief Get class for managing various visualization features */
  VisualizerPtr getVisual()
  {
    return visual_;
  }

  /** \brief Free all the memory allocated by the database */
  void freeMemory();

  /** \brief Initialize database */
  bool setup();

  /* ---------------------------------------------------------------------------------
   * Astar search
   * --------------------------------------------------------------------------------- */

  /** \brief Check that the query vertex is initialized (used for internal nearest neighbor searches) */
  void initializeQueryState();

  /** \brief Given two milestones from the same connected component, construct a path connecting them and set it as
   * the solution
   *  \param start
   *  \param goal
   *  \param vertexPath
   *  \return true if candidate solution found
   */
  bool astarSearch(const TaskVertex start, const TaskVertex goal, std::vector<TaskVertex>& vertexPath,
                   double& distance, std::size_t indent);

  /** \brief Distance between two states with special bias using popularity */
  double astarHeuristic(const TaskVertex a, const TaskVertex b) const;

  /** \brief Compute distance between two milestones (this is simply distance between the states of the milestones) */
  double distanceFunction(const TaskVertex a, const TaskVertex b) const;

  /** \brief Custom A* visitor statistics */
  void recordNodeOpened()  // discovered
  {
    numNodesOpened_++;
  }

  void recordNodeClosed()  // examined
  {
    numNodesClosed_++;
  }

  /* ---------------------------------------------------------------------------------
   * Get graph properties
   * --------------------------------------------------------------------------------- */

  const std::size_t getNumQueryVertices() const
  {
    return queryVertices_.size();
  }

  /** \brief Get the number of vertices in the task roadmap. */
  unsigned int getNumVertices() const
  {
    return boost::num_vertices(g_);
  }

  /** \brief Get the number of edges in the task roadmap. */
  unsigned int getNumEdges() const
  {
    return boost::num_edges(g_);
  }

  VertexType getVertexTypeProperty(TaskVertex v) const
  {
    return vertexTypeProperty_[v];
  }

  double getEdgeWeightProperty(TaskEdge e) const
  {
    return edgeWeightProperty_[e];
  }

  /** \brief Determine if no nodes or edges have been added to the graph except query vertices */
  bool isEmpty() const;

  /* ---------------------------------------------------------------------------------
   * Task Planning
   * --------------------------------------------------------------------------------- */

  void generateTaskSpace(std::size_t indent);

  /* ---------------------------------------------------------------------------------
   * Error checking
   * --------------------------------------------------------------------------------- */

  /** \brief Clear all past edge state information about in collision or not */
  void clearEdgeCollisionStates();

  /** \brief Part of super debugging */
  void errorCheckDuplicateStates(std::size_t indent);

  /* ---------------------------------------------------------------------------------
   * Smoothing
   * --------------------------------------------------------------------------------- */

  /** \brief Path smoothing helpers */
  bool smoothQualityPathOriginal(geometric::PathGeometric* path, std::size_t indent);
  bool smoothQualityPath(geometric::PathGeometric* path, double clearance, std::size_t indent);

  /* ---------------------------------------------------------------------------------
   * Disjoint Sets
   * --------------------------------------------------------------------------------- */

  /** \brief Disjoint sets analysis tools */
  std::size_t getDisjointSetsCount(bool verbose = false);
  void getDisjointSets(TaskDisjointSetsMap& disjointSets);
  void printDisjointSets(TaskDisjointSetsMap& disjointSets);
  void visualizeDisjointSets(TaskDisjointSetsMap& disjointSets);
  std::size_t checkConnectedComponents();
  bool sameComponent(TaskVertex v1, TaskVertex v2);

  /* ---------------------------------------------------------------------------------
   * Add/remove vertices, edges, states
   * --------------------------------------------------------------------------------- */

  /** \brief Add a state to the DenseCache */
  StateID addState(base::State* state);

  /** \brief Add vertices to graph */
  TaskVertex addVertex(base::State* state, const VertexType& type, std::size_t level, std::size_t indent);
  TaskVertex addVertex(StateID stateID, const VertexType& type, std::size_t level, std::size_t indent);

  /** \brief Remove vertex from graph */
  void removeVertex(TaskVertex v);

  /** \brief Cleanup graph because we leave deleted vertices in graph during construction */
  void removeDeletedVertices(std::size_t indent);

  /** \brief Display in viewer */
  void visualizeVertex(TaskVertex v, const VertexType& type);

  /** \brief Add edge to graph */
  TaskEdge addEdge(TaskVertex v1, TaskVertex v2, EdgeType type, std::size_t indent);

  /** \brief Check graph for edge existence */
  bool hasEdge(TaskVertex v1, TaskVertex v2);

  /** \brief Get the state of a vertex used for querying - i.e. vertices 0-11 for 12 thread system */
  base::State*& getQueryStateNonConst(TaskVertex v);

  /** \brief Shortcut function for getting the state of a vertex */
  base::State*& getVertexStateNonConst(TaskVertex v);
  const base::State* getVertexState(TaskVertex v) const;
  const base::State* getState(StateID stateID) const;
  const StateID getStateID(TaskVertex v) const;

  /* ---------------------------------------------------------------------------------
   * Visualizations
   * --------------------------------------------------------------------------------- */

  /** \brief Show in visualizer the task graph */
  void displayDatabase(bool showVertices = false, std::size_t indent = 0);

  /* ---------------------------------------------------------------------------------
   * Debug Utilities
   * --------------------------------------------------------------------------------- */

  /** \brief Print info to console */
  void debugState(const ompl::base::State* state);

  /** \brief Print nearest neighbor info to console */
  void debugNN();

protected:
  /** \brief Short name of this class */
  const std::string name_ = "TaskGraph";

  /** \brief Sparse graph main datastructure that this class operates on */
  SparseGraphPtr sg_;

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  /** \brief Speed up collision checking by saving redundant checks and using file storage */
  DenseCachePtr denseCache_;

  /** \brief Nearest neighbors data structure */
  std::shared_ptr<NearestNeighbors<TaskVertex> > nn_;

  /** \brief Connectivity graph */
  TaskAdjList g_;

  /** \brief Vertices for performing nearest neighbor queries on multiple threads */
  std::vector<TaskVertex> queryVertices_;
  std::vector<base::State*> queryStates_;

  /** \brief Access to the weights of each Edge */
  boost::property_map<TaskAdjList, boost::edge_weight_t>::type edgeWeightProperty_;

  /** \brief Access to the collision checking state of each Edge */
  TaskEdgeCollisionStateMap edgeCollisionStatePropertyTask_;

  /** \brief Access to the internal base::state at each Vertex */
  boost::property_map<TaskAdjList, vertex_state_cache_t>::type vertexStateProperty_;

  /** \brief Access to additional task level dimension at each Vertex */
  boost::property_map<TaskAdjList, vertex_level_t>::type vertexLevelProperty_;

  /** \brief Access to the SPARS vertex type for the vertices */
  boost::property_map<TaskAdjList, vertex_type_t>::type vertexTypeProperty_;

  /** \brief Access to  */
  boost::property_map<TaskAdjList, vertex_task_mirror_t>::type vertexTaskMirrorProperty_;

  /** \brief Data structure that maintains the connected components */
  TaskDisjointSetType disjointSets_;

  /** \brief A path simplifier used to simplify dense paths added to S */
  geometric::PathSimplifierPtr pathSimplifier_;

  /** \brief Number of cores available on system */
  std::size_t numThreads_;

  /** \brief Astar statistics */
  std::size_t numNodesOpened_ = 0;
  std::size_t numNodesClosed_ = 0;

  /** \brief Track if the graph has been modified */
  bool graphUnsaved_ = false;

public:  // user settings from other applications

  /** \brief Visualization speed of astar search, num of seconds to show each vertex */
  bool visualizeAstar_ = false;
  double visualizeAstarSpeed_ = 0.1;
  bool visualizeQualityPathSimp_ = false;

  /** \brief Change verbosity levels */
  bool vVisualize_ = false;
  bool vAdd_ = false;  // message when adding edges and vertices
  bool vSearch_ = false;

  /** \brief Run with extra safety checks */
  bool superDebug_ = true;

  /** \brief Show the task graph being generated */
  bool visualizeTaskGraph_ = false;
  double visualizeTaskGraphSpeed_ = 0.0;
  bool visualizeDatabaseVertices_ = true;
  bool visualizeDatabaseEdges_ = true;
  bool visualizeDatabaseCoverage_ = true;

};  // end class TaskGraph

////////////////////////////////////////////////////////////////////////////////////////
/**
 * Vertex visitor to check if A* search is finished.
 * \implements AStarVisitorConcept
 * See http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/AStarVisitor.html
 */
class TaskAstarVisitor : public boost::default_astar_visitor
{
private:
  TaskVertex goal_;  // Goal Vertex of the search
  TaskGraph* parent_;

public:
  /**
   * Construct a visitor for a given search.
   * \param goal  goal vertex of the search
   */
  TaskAstarVisitor(TaskVertex goal, TaskGraph* parent);

  /**
   * \brief Invoked when a vertex is first discovered and is added to the OPEN list.
   * \param v current Vertex
   * \param g graph we are searching on
   */
  void discover_vertex(TaskVertex v, const TaskAdjList& g) const;

  /**
   * \brief Check if we have arrived at the goal.
   * This is invoked on a vertex as it is popped from the queue (i.e., it has the lowest
   * cost on the OPEN list). This happens immediately before examine_edge() is invoked on
   * each of the out-edges of vertex u.
   * \param v current vertex
   * \param g graph we are searching on
   * \throw FoundGoalException if \a u is the goal
   */
  void examine_vertex(TaskVertex v, const TaskAdjList& g) const;
};  // end TaskGraph

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_TASK_GRAPH_
