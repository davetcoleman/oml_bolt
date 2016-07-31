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

#ifndef OMPL_TOOLS_BOLT_SPARSE_GRAPH_
#define OMPL_TOOLS_BOLT_SPARSE_GRAPH_

#include <ompl/base/StateSpace.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/PathSimplifier.h>

#include <ompl/base/SpaceInformation.h>
#include <ompl/datastructures/NearestNeighbors.h>

// Bolt
#include <ompl/tools/debug/Visualizer.h>
#include <ompl/tools/bolt/BoostGraphHeaders.h>
#include <ompl/tools/bolt/Debug.h>
#include <ompl/tools/bolt/VertexDiscretizer.h>
#include <ompl/tools/bolt/SparseStorage.h>

// Boost
#include <boost/function.hpp>
#include <boost/graph/astar_search.hpp>

// C++
#include <list>
#include <random>
#include <mutex>

namespace ompl
{
namespace tools
{
namespace bolt
{
/**
   @anchor SparseGraph
   @par Near-asypmotically optimal roadmap datastructure
*/

/// @cond IGNORE
OMPL_CLASS_FORWARD(SparseGraph);
OMPL_CLASS_FORWARD(SparseCriteria);
/// @endcond

/** \class ompl::tools::bolt::::SparseGraphPtr
    \brief A boost shared pointer wrapper for ompl::tools::SparseGraph */

/** \brief Near-asypmotically optimal roadmap datastructure */
class SparseGraph
{
public:
  /** \brief Constructor needs the state space used for planning.
   */
  SparseGraph(base::SpaceInformationPtr si, VisualizerPtr visual);

  /** \brief Deconstructor */
  virtual ~SparseGraph(void);

  /* ---------------------------------------------------------------------------------
   * Setup and cleanup
   * --------------------------------------------------------------------------------- */

  /** \brief Give the sparse graph reference to the criteria, because sometimes it needs data from there */
  void setSparseCriteria(SparseCriteriaPtr sparseCriteria)
  {
    sparseCriteria_ = sparseCriteria;
  }

  /** \brief Retrieve the computed roadmap. */
  const SparseAdjList& getGraph() const
  {
    return g_;
  }

  SparseAdjList getGraphNonConst()
  {
    return g_;
  }

  SparseStoragePtr getSparseStorage()
  {
    return sparseStorage_;
  }

  SparseCriteriaPtr getSparseCriteria()
  {
    return sparseCriteria_;
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

  /** \brief Get the nearest neighbor structure */
  std::shared_ptr<NearestNeighbors<SparseVertex> > getNN()
  {
    return nn_;
  }

  std::mutex& getNNGuard()
  {
    return nearestNeighborMutex_;
  }

  /** \brief Free all the memory allocated by the database */
  void freeMemory();

  /** \brief Initialize database */
  bool setup();

  /** \brief Clear data on why vertices and edges were added by optimality criteria */
  void clearStatistics();

  /* ---------------------------------------------------------------------------------
   * Astar search
   * --------------------------------------------------------------------------------- */

  /** \brief Check that the query vertex is initialized (used for internal nearest neighbor searches) */
  void initializeQueryState();

  /** \brief Set the file path to load/save to/from */
  void setFilePath(const std::string& filePath)
  {
    filePath_ = filePath;
  }

  /**
   * \brief Load database from file
   * \return true if file loaded successfully
   */
  bool load();

  /**
   * \brief Save loaded database to file, except skips saving if no paths have been added
   * \return true if file saved successfully
   */
  bool saveIfChanged(std::size_t indent = 0);

  /**
   * \brief Save loaded database to file
   * \return true if file saved successfully
   */
  bool save(std::size_t indent = 0);

  /** \brief Given two milestones from the same connected component, construct a path connecting them and set it as
   * the solution
   *  \param start
   *  \param goal
   *  \param vertexPath
   *  \return true if candidate solution found
   */
  bool astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex>& vertexPath,
                   double& distance, std::size_t indent);

  /** \brief Distance between two states with special bias using popularity */
  double astarHeuristic(const SparseVertex a, const SparseVertex b) const;

  /** \brief Compute distance between two milestones (this is simply distance between the states of the milestones) */
  double distanceFunction(const SparseVertex a, const SparseVertex b) const;

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

  std::size_t getNumQueryVertices() const
  {
    return queryVertices_.size();
  }

  /** \brief Get the number of vertices in the sparse roadmap INCLUDING query vertices */
  unsigned int getNumVertices() const
  {
    return boost::num_vertices(g_);
  }

  /** \brief Get the number of vertices in the sparse roadmap not including query vertices. */
  unsigned int getNumRealVertices() const
  {
    return boost::num_vertices(g_) - queryVertices_.size();
  }

  /** \brief Get the number of edges in the sparse roadmap. */
  unsigned int getNumEdges() const
  {
    return boost::num_edges(g_);
  }

  VertexType getVertexTypeProperty(SparseVertex v) const
  {
    return vertexTypeProperty_[v];
  }

  double getEdgeWeightProperty(SparseEdge e) const
  {
    return edgeWeightProperty_[e];
  }

  EdgeType getEdgeTypeProperty(SparseEdge e) const
  {
    return edgeTypeProperty_[e];
  }

  /** \brief Determine if no nodes or edges have been added to the graph except query vertices */
  bool isEmpty() const;

  /* ---------------------------------------------------------------------------------
   * Sampling Properties
   * --------------------------------------------------------------------------------- */

  double getObstacleClearance()
  {
    return obstacleClearance_;
  }

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
  bool smoothQualityPathOriginal(geometric::PathGeometric* path, std::size_t indent) const;
  bool smoothQualityPath(geometric::PathGeometric* path, double clearance, std::size_t indent) const;

  /* ---------------------------------------------------------------------------------
   * Disjoint Sets
   * --------------------------------------------------------------------------------- */

  /** \brief Disjoint sets analysis tools */
  std::size_t getDisjointSetsCount(bool verbose = false) const;
  void getDisjointSets(SparseDisjointSetsMap& disjointSets);
  void printDisjointSets(SparseDisjointSetsMap& disjointSets);
  void visualizeDisjointSets(SparseDisjointSetsMap& disjointSets);
  std::size_t checkConnectedComponents();
  bool sameComponent(SparseVertex v1, SparseVertex v2);

  /* ---------------------------------------------------------------------------------
   * Add/remove vertices, edges, states
   * --------------------------------------------------------------------------------- */

  /** \brief Add vertex to graph */
  SparseVertex addVertex(base::State* state, const VertexType& type, std::size_t indent);

  /** \brief Quickly add vertex to graph when loading from file */
  SparseVertex addVertexFromFile(base::State* state, const VertexType& type, std::size_t indent);

  /** \brief Remove vertex from graph */
  void removeVertex(SparseVertex v);

  /** \brief Cleanup graph because we leave deleted vertices in graph during construction */
  void removeDeletedVertices(std::size_t indent);

  /** \brief Add edge to graph */
  SparseEdge addEdge(SparseVertex v1, SparseVertex v2, EdgeType type, std::size_t indent);

  /** \brief Check graph for edge existence */
  bool hasEdge(SparseVertex v1, SparseVertex v2);

  /** \brief Helper for choosing an edge's display color based on type of edge */
  VizColors edgeTypeToColor(EdgeType edgeType);

  /** \brief Get the state of a vertex used for querying - i.e. vertices 0-11 for 12 thread system */
  base::State*& getQueryStateNonConst(std::size_t threadID);
  SparseVertex getQueryVertices(std::size_t threadID);

  /** \brief Shortcut function for getting the state of a vertex */
  base::State*& getStateNonConst(SparseVertex v);
  const base::State* getState(SparseVertex v) const;

  /** \brief Used for creating a voronoi diagram */
  SparseVertex getSparseRepresentative(base::State* state);

  /** \brief When a new guard is added at state, finds all guards who must abandon their interface information and
   * delete that information */
  void clearInterfaceData(base::State* st);

  /** \brief When a quality path is added with new vertices, remove all edges near the new vertex */
  void clearEdgesNearVertex(SparseVertex vertex, std::size_t indent);

  /* ---------------------------------------------------------------------------------
   * Visualizations
   * --------------------------------------------------------------------------------- */

  /** \brief Show in visualizer the sparse graph */
  void displayDatabase(bool showVertices = true, bool showEdges = true, std::size_t windowID = 1, std::size_t indent = 0);

  /** \brief Display in viewer */
  void visualizeVertex(SparseVertex v, const VertexType& type);

  /** \brief Convert type of vertex to a visualization color */
  tools::VizColors vertexTypeToColor(VertexType type);

  /** \brief Display in viewer */
  void visualizeEdge(SparseEdge e, EdgeType type, std::size_t windowID);

  /** \brief Display in viewer */
  void visualizeEdge(SparseVertex v1, SparseVertex v2, EdgeType type, std::size_t windowID);

  /* ---------------------------------------------------------------------------------
   * Sparse Interfaces
   * --------------------------------------------------------------------------------- */

  /** \brief Rectifies indexing order for accessing the vertex data */
  VertexPair interfaceDataIndex(SparseVertex vp, SparseVertex vpp);

  /** \brief Retrieves the Vertex data associated with v,vp,vpp */
  InterfaceData& getInterfaceData(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent);

  /** \brief Retrieves the vertex data associated with v */
  InterfaceHash& getVertexInterfaceProperty(SparseVertex v);

  /* ---------------------------------------------------------------------------------
   * Debug Utilities
   * --------------------------------------------------------------------------------- */

  /** \brief Print info to console */
  void debugState(const ompl::base::State* state);

  /** \brief Print info to console */
  void debugVertex(const SparseVertex v);

  /** \brief Print nearest neighbor info to console */
  void debugNN();

  /** \brief Information about the loaded graph */
  void printGraphStats();

protected:
  /** \brief Short name of this class */
  const std::string name_ = "SparseGraph";

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  /** \brief Class for deciding which vertices and edges get added */
  SparseCriteriaPtr sparseCriteria_;

  /** \brief For saving and loading to file */
  SparseStoragePtr sparseStorage_;

  /** \brief Nearest neighbors data structure */
  std::shared_ptr<NearestNeighbors<SparseVertex> > nn_;

  /** \brief Connectivity graph */
  SparseAdjList g_;

  /** \brief Vertices for performing nearest neighbor queries on multiple threads */
  std::vector<SparseVertex> queryVertices_;
  std::vector<base::State*> queryStates_;

  /** \brief Access to the weights of each Edge */
  boost::property_map<SparseAdjList, boost::edge_weight_t>::type edgeWeightProperty_;

  /** \brief Access to the type (reason for being added in SPARS) of each Edge */
  boost::property_map<SparseAdjList, edge_type_t>::type edgeTypeProperty_;

  /** \brief Access to the collision checking state of each Edge */
  SparseEdgeCollisionStateMap edgeCollisionStatePropertySparse_;

  /** \brief Access to the internal base::state at each Vertex */
  boost::property_map<SparseAdjList, vertex_state_t>::type vertexStateProperty_;

  /** \brief Access to the SPARS vertex type for the vertices */
  boost::property_map<SparseAdjList, vertex_type_t>::type vertexTypeProperty_;

  /** \brief Access to the interface pair information for the vertices */
  boost::property_map<SparseAdjList, vertex_interface_data_t>::type vertexInterfaceProperty_;

  /** \brief Access to the popularity of each node */
  boost::property_map<SparseAdjList, vertex_popularity_t>::type vertexPopularity_;

  /** \brief Data structure that maintains the connected components */
  SparseDisjointSetType disjointSets_;

  /** \brief A path simplifier used to simplify dense paths added to S */
  geometric::PathSimplifierPtr pathSimplifier_;

  /** \brief Track where to load/save datastructures */
  std::string filePath_;

  /** \brief Number of cores available on system */
  std::size_t numThreads_;

  /** \brief Astar statistics */
  std::size_t numNodesOpened_ = 0;
  std::size_t numNodesClosed_ = 0;

  /** \brief Track if the graph has been modified */
  bool graphUnsaved_ = false;

  tools::VizSizes vertexSize_ = tools::LARGE;
  tools::VizSizes edgeSize_ = tools::MEDIUM;

  // Multi-threading modifying graph
  std::mutex nearestNeighborMutex_;
  std::mutex modifyGraphMutex_;
  time::point lastSampledModTime_;  // timestamp of the last graph modification - any sample taken before that is
                                    // invalid

public:  // user settings from other applications
  /** \brief For statistics */
  int numSamplesAddedForCoverage_ = 0;
  int numSamplesAddedForConnectivity_ = 0;
  int numSamplesAddedForInterface_ = 0;
  int numSamplesAddedForQuality_ = 0;

  /** \brief Clearance of obstacles in order to be considered "cl-robust" as described in paper */
  double obstacleClearance_;

  /** \brief Allow the database to save to file (new experiences) */
  bool savingEnabled_ = true;

  /** \brief Various options for visualizing the algorithmns performance */
  bool visualizeAstar_ = false;

  /** \brief Visualization speed of astar search, num of seconds to show each vertex */
  double visualizeAstarSpeed_ = 0.1;
  bool visualizeQualityPathSimp_ = false;

  /** \brief Change verbosity levels */
  bool vVisualize_ = false;
  bool vAdd_ = false;  // message when adding edges and vertices
  bool vSearch_ = false;

  /** \brief Run with extra safety checks */
  bool superDebug_ = false;

  /** \brief Show the sparse graph being generated */
  bool visualizeSparseGraph_ = false;
  double visualizeSparseGraphSpeed_ = 0.0;
  bool visualizeDatabaseVertices_ = true;
  bool visualizeDatabaseEdges_ = true;
  bool visualizeDatabaseCoverage_ = true;
  bool visualizeGraphAfterLoading_ = true;
  bool visualizeProjection_ = false;
  bool visualizeVoronoiDiagram_ = true;
  bool visualizeVoronoiDiagramAnimated_ = true;

};  // end class SparseGraph

////////////////////////////////////////////////////////////////////////////////////////
/**
 * Vertex visitor to check if A* search is finished.
 * \implements AStarVisitorConcept
 * See http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/AStarVisitor.html
 */
class SparsestarVisitor : public boost::default_astar_visitor
{
private:
  SparseVertex goal_;  // Goal Vertex of the search
  SparseGraph* parent_;

public:
  /**
   * Construct a visitor for a given search.
   * \param goal  goal vertex of the search
   */
  SparsestarVisitor(SparseVertex goal, SparseGraph* parent);

  /**
   * \brief Invoked when a vertex is first discovered and is added to the OPEN list.
   * \param v current Vertex
   * \param g graph we are searching on
   */
  void discover_vertex(SparseVertex v, const SparseAdjList& g) const;

  /**
   * \brief Check if we have arrived at the goal.
   * This is invoked on a vertex as it is popped from the queue (i.e., it has the lowest
   * cost on the OPEN list). This happens immediately before examine_edge() is invoked on
   * each of the out-edges of vertex u.
   * \param v current vertex
   * \param g graph we are searching on
   * \throw FoundGoalException if \a u is the goal
   */
  void examine_vertex(SparseVertex v, const SparseAdjList& g) const;
};  // end SparseGraph

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_SPARSE_GRAPH_
