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
   Desc:   Sparse experience database for storing and reusing past path plans
*/

#ifndef OMPL_TOOLS_BOLT_SPARSEDB_
#define OMPL_TOOLS_BOLT_SPARSEDB_

#include <ompl/base/StateSpace.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/PathSimplifier.h>
#include <ompl/base/Planner.h>
//#include <ompl/base/State.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include <ompl/base/PlannerTerminationCondition.h>
#include <ompl/tools/debug/Visualizer.h>
#include <ompl/tools/bolt/DenseDB.h>
#include <ompl/tools/bolt/BoltGraph.h>
#include <ompl/tools/bolt/EdgeCache.h>

// Boost
#include <boost/function.hpp>

// C++
#include <list>

namespace ompl
{
namespace tools
{
namespace bolt
{
/**
   @anchor SparsDB
   @par Short description
   Database for storing and retrieving past plans
*/

/// @cond IGNORE
OMPL_CLASS_FORWARD(SparseDB);
OMPL_CLASS_FORWARD(DenseDB);
/// @endcond

// clang-format off
#define BOLT_DEBUG(indent, token) { std::stringstream o; o << token; \
    { if (checksVerbose_)                                               \
        std::cout << std::string(indent, ' ') << o.str() << std::endl; } }
#define BOLT_COLOR_DEBUG(indent, token, color) { std::stringstream o; o << token; \
    { if (checksVerbose_) {                                             \
        std::cout << color << std::string(indent, ' ') << o.str() << ANSI_COLOR_RESET << std::endl; } } }
#define BOLT_RED_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_RED); }
#define BOLT_GREEN_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_GREEN); }
#define BOLT_YELLOW_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_YELLOW); }
#define BOLT_BLUE_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_BLUE); }
#define BOLT_MAGENTA_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_MAGENTA); }
#define BOLT_CYAN_DEBUG(indent, token) { BOLT_COLOR_DEBUG(indent, token, ANSI_COLOR_CYAN); }
// clang-format on

////////////////////////////////////////////////////////////////////////////////////////
/**
 * Vertex visitor to check if A* search is finished.
 * \implements AStarVisitorConcept
 * See http://www.boost.org/doc/libs/1_58_0/libs/graph/doc/AStarVisitor.html
 */
class CustomAstarVisitor : public boost::default_astar_visitor
{
private:
  SparseVertex goal_;  // Goal Vertex of the search
  SparseDB* parent_;

public:
  /**
   * Construct a visitor for a given search.
   * \param goal  goal vertex of the search
   */
  CustomAstarVisitor(SparseVertex goal, SparseDB* parent);

  /**
   * \brief Invoked when a vertex is first discovered and is added to the OPEN list.
   * \param v current Vertex
   * \param g graph we are searching on
   */
  void discover_vertex(SparseVertex v, const SparseGraph& g) const;

  /**
   * \brief Check if we have arrived at the goal.
   * This is invoked on a vertex as it is popped from the queue (i.e., it has the lowest
   * cost on the OPEN list). This happens immediately before examine_edge() is invoked on
   * each of the out-edges of vertex u.
   * \param v current vertex
   * \param g graph we are searching on
   * \throw FoundGoalException if \a u is the goal
   */
  void examine_vertex(SparseVertex v, const SparseGraph& g) const;
};

/** \class ompl::tools::bolt::::SparseDBPtr
    \brief A boost shared pointer wrapper for ompl::tools::SparseDB */

/** \brief Save and load entire paths from file */
class SparseDB
{
  friend class BoltRetrieveRepair;
  friend class DenseDB;
  friend class Discretizer;

public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // SparseDB MEMBER FUNCTIONS
  ////////////////////////////////////////////////////////////////////////////////////////

  /** \brief Constructor needs the state space used for planning.
   *  \param space - state space
   */
  SparseDB(base::SpaceInformationPtr si, DenseDB* denseDB, VisualizerPtr visual, EdgeCachePtr edgeCache);

  /** \brief Deconstructor */
  virtual ~SparseDB(void);

  /** \brief Initialize database */
  bool setup();

  /** \brief Given two milestones from the same connected component, construct a path connecting them and set it as
   * the solution
   *  \param start
   *  \param goal
   *  \param vertexPath
   *  \return true if candidate solution found
   */
  bool astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex>& vertexPath, double &distance, std::size_t indent);

  /** \brief Distance between two states with special bias using popularity */
  double astarHeuristic(const SparseVertex a, const SparseVertex b) const;

  /** \brief Print info to screen */
  void debugVertex(const ompl::base::PlannerDataVertex& vertex);
  void debugState(const ompl::base::State* state);

  /** \brief Retrieve the computed roadmap. */
  const SparseGraph& getRoadmap() const
  {
    return g_;
  }

  /** \brief Get the number of vertices in the sparse roadmap. */
  unsigned int getNumVertices() const
  {
    return boost::num_vertices(g_);
  }

  /** \brief Get the number of edges in the sparse roadmap. */
  unsigned int getNumEdges() const
  {
    return boost::num_edges(g_);
  }

  /** \brief Free all the memory allocated by the database */
  void freeMemory();

  /**
   * \brief Check if anything has been loaded into DB
   * \return true if has no nodes
   */
  bool isEmpty()
  {
    return !getNumVertices();
  }

  /** \brief Check that the query vertex is initialized (used for internal nearest neighbor searches) */
  void initializeQueryState();

  /** \brief Clear all past edge state information about in collision or not */
  void clearEdgeCollisionStates();

  /** \brief Utilize multi-threading by using lots of caching */
  void preprocessSPARS();
  void preprocessSPARSThread(std::size_t threadID, std::size_t numThreads, base::SpaceInformationPtr si);

  /** \brief Create a SPARS graph from the discretized dense graph and its popularity metric */
  void createSPARS();
  void createSPARSOuterLoop();
  bool createSPARSInnerLoop(std::list<WeightedVertex>& vertexInsertionOrder, std::size_t& sucessfulInsertions);

  /** \brief Helper function for choosing the correct method for vertex insertion ordering */
  void getVertexInsertionOrdering(std::list<WeightedVertex>& vertexInsertionOrder);

  void addRandomSamples();

  /** \brief Helper for counting the number of disjoint sets in the sparse graph */
  std::size_t getDisjointSetsCount(bool verbose = false);

  bool getPopularityOrder(std::list<WeightedVertex>& vertexInsertionOrder);
  bool getDefaultOrder(std::list<WeightedVertex>& vertexInsertionOrder);
  bool getRandomOrder(std::list<WeightedVertex>& vertexInsertionOrder);

  /** \brief Helper function for random integer creation */
  int iRand(int min, int max);

  /**
   * \brief Run various checks/criteria to determine if to keep DenseVertex in sparse graph
   * \param denseVertex - the original vertex to consider
   * \param newVertex - if function returns true, the newly generated sparse vertex
   * \param addReason - if function returns true, the reson the denseVertex was added to the sparse graph
   * \return true on success
   */
  bool addStateToRoadmap(DenseVertex denseVertex, SparseVertex& newVertex, GuardType& addReason, std::size_t threadID);

  /* ----------------------------------------------------------------------------------------*/
  /** \brief SPARS-related functions */
  bool checkAddCoverage(DenseVertex denseV, std::vector<SparseVertex>& visibleNeighborhood, SparseVertex& newVertex,
                        std::size_t indent);
  bool checkAddConnectivity(DenseVertex denseV, std::vector<SparseVertex>& visibleNeighborhood, SparseVertex& newVertex,
                            std::size_t indent);
  bool checkAddInterface(DenseVertex denseV, std::vector<SparseVertex>& graphNeighborhood,
                         std::vector<SparseVertex>& visibleNeighborhood, SparseVertex& newVertex, std::size_t indent);
  bool checkAddQuality(DenseVertex denseV, std::vector<SparseVertex>& graphNeighborhood,
                       std::vector<SparseVertex>& visibleNeighborhood, base::State* workState, SparseVertex& newVertex,
                       std::size_t indent);

  /* ----------------------------------------------------------------------------------------*/
  // 4th Criteria
  /* ----------------------------------------------------------------------------------------*/

  /** \brief Checks vertex v for short paths through its region and adds when appropriate.
   *         Referred to as 'Test_Add_paths' in paper
   */
  bool checkAddPath(SparseVertex v, std::size_t indent);

  bool checkAddPathHelper(SparseVertex v, SparseVertex vp, SparseVertex vpp, InterfaceData &iData, std::size_t indent);

  /** \brief Finds the representative of the input state, st  */
  SparseVertex findGraphRepresentative(base::State* st, std::size_t indent);

  /** \brief Finds representatives of samples near candidateState_ which are not his representative
             Referred to as 'Get_Close_Reps' in paper
   */
  void findCloseRepresentatives(base::State* workState, const base::State* candidateState, SparseVertex candidateRep,
                                std::map<SparseVertex, base::State*>& closeRepresentatives, std::size_t indent);

  /** \brief Updates pair point information for a representative with neighbor r
             Referred to as 'Update_Points' in paper
      \return true if an update actually happend wihtin the representatives, false if no change
   */
  bool updatePairPoints(SparseVertex candidateRep, const base::State* candidateState, SparseVertex nearSampledRep,
                        const base::State* nearSampledState, std::size_t indent);

  /** \brief Computes all nodes which qualify as a candidate v" for v and vp */
  void getAdjVerticesOfV1UnconnectedToV2(SparseVertex v1, SparseVertex v2,
                                         std::vector<SparseVertex>& adjVerticesUnconnected, std::size_t indent);

  /** \brief Computes all nodes which qualify as a candidate x for v, v', and v"
   *  \return length of maximum path
   */
  double maxSpannerPath(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent);

  /** \brief Rectifies indexing order for accessing the vertex data */
  VertexPair index(SparseVertex vp, SparseVertex vpp);

  /** \brief Retrieves the Vertex data associated with v,vp,vpp */
  InterfaceData& getData(SparseVertex v, SparseVertex vp, SparseVertex vpp, std::size_t indent);

  /** \brief Performs distance checking for the candidate new state, q against the current information
      \return true if an update actually happend wihtin the representatives, false if no change
  */
  bool distanceCheck(SparseVertex v, const base::State* q, SparseVertex vp,
                     const base::State* qp, SparseVertex vpp, std::size_t indent);

  /** \brief When a new guard is added at state st, finds all guards who must abandon their interface information and
   * deletes that information */
  void abandonLists(base::State* st);

  /* ----------------------------------------------------------------------------------------*/

  /**
   * \brief Get neighbors within sparseDelta radius
   * \param denseV - origin state to search from
   * \param graphNeighborhood - resulting nearby states
   * \param visibleNeighborhood - resulting nearby states that are visible
   * \param indent - debugging tool
   */
  void findGraphNeighbors(DenseVertex denseV, std::vector<SparseVertex>& graphNeighborhood,
                          std::vector<SparseVertex>& visibleNeighborhood, std::size_t threadID, std::size_t indent);

  DenseVertex getInterfaceNeighbor(DenseVertex q, SparseVertex rep);

  bool sameComponent(SparseVertex v1, SparseVertex v2);

  /** \brief Add vertices to graph */
  SparseVertex addVertex(base::State* state, const GuardType& type);
  SparseVertex addVertex(DenseVertex denseV, const GuardType& type);

  /** \brief Add edge to graph */
  void addEdge(SparseVertex v1, SparseVertex v2, std::size_t visualColor, std::size_t indent);

  std::size_t getVizVertexType(const GuardType& type);

  /** \brief Show in visualizer the sparse graph */
  void displaySparseDatabase(bool showVertices = false);

  EdgeCachePtr getEdgeCache()
  {
    return edgeCache_;
  }

  /** \brief Shortcut function for getting the state of a vertex */
  base::State*& getSparseStateNonConst(SparseVertex v);
  const base::State* getSparseState(SparseVertex v) const;
  base::State*& getDenseState(DenseVertex denseV);

  /** \brief Compute distance between two milestones (this is simply distance between the states of the milestones) */
  double distanceFunction(const SparseVertex a, const SparseVertex b) const;

  double getSecondarySparseDelta();

  bool hasEdge(SparseVertex v1, SparseVertex v2);

  void visualizeInterfaces(SparseVertex v);

  SparseVertex getSparseRepresentative(base::State* state);

  /** \brief Custom A* visitor statistics */
  void recordNodeOpened()  // discovered
  {
    numNodesOpened_++;
  }
  void recordNodeClosed()  // examined
  {
    numNodesClosed_++;
  }

  /** \brief Get class for managing various visualization features */
  VisualizerPtr getVisual()
  {
    return visual_;
  }

protected:
  /** \brief Short name of this class */
  const std::string name_ = "SparseDB";

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief The database of motions to search through */
  DenseDB* denseDB_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  /** \brief Speed up collision checking by saving redundant checks and using file storage */
  EdgeCachePtr edgeCache_;

  /** \brief Nearest neighbors data structure */
  std::shared_ptr<NearestNeighbors<SparseVertex> > nn_;

  /** \brief Connectivity graph */
  SparseGraph g_;

  /** \brief Vertices for performing nearest neighbor queries on multiple threads */
  std::vector<SparseVertex> queryVertices_;

  /** \brief Geometric Path variable used for smoothing out paths. */
  geometric::PathGeometric smoothingGeomPath_;

  /** \brief Access to the weights of each Edge */
  boost::property_map<SparseGraph, boost::edge_weight_t>::type edgeWeightPropertySparse_;

  /** \brief Access to the collision checking state of each Edge */
  SparseEdgeCollisionStateMap edgeCollisionStatePropertySparse_;

  /** \brief Access to the internal base::state at each Vertex */
  boost::property_map<SparseGraph, vertex_dense_pointer_t>::type denseVertexProperty_;

  /** \brief Access to the SPARS vertex type for the vertices */
  boost::property_map<SparseGraph, vertex_type_t>::type typePropertySparse_;

  /** \brief Access to the interface pair information for the vertices */
  boost::property_map<SparseGraph, vertex_interface_data_t>::type interfaceDataProperty_;

  /** \brief Access to the popularity of each node */
  boost::property_map<SparseGraph, vertex_popularity_t>::type vertexPopularity_;

  /** \brief Data structure that maintains the connected components */
  boost::disjoint_sets<boost::property_map<SparseGraph, boost::vertex_rank_t>::type,
                       boost::property_map<SparseGraph, boost::vertex_predecessor_t>::type> disjointSets_;

  /** \brief A path simplifier used to simplify dense paths added to S */
  geometric::PathSimplifierPtr pathSimplifier_;

  /** \brief Sampler user for generating valid samples in the state space */
  base::ValidStateSamplerPtr sampler_;

  /** \brief Special flag for tracking mode when inserting into sparse graph */
  bool secondSparseInsertionAttempt_ = false;

  /** \brief Amount of sub-optimality allowed */
  double sparseDelta_ = 2.0;

  /** \brief SPARS parameter for dense graph connection distance */
  double denseDelta_;

  /** \brief Number of sample points to use when trying to detect interfaces. */
  std::size_t nearSamplePoints_;

  /** \brief Show what nodes are added on top of the regular SPARS graph */
  bool visualizeOverlayNodes_ = false;

  /** \brief Cache the maximum extent for later re-use */
  double maxExtent_;

  bool useFourthCriteria_;

  /** \brief Astar statistics */
  std::size_t numNodesOpened_ = 0;
  std::size_t numNodesClosed_ = 0;

public:

  /** \brief Various options for visualizing the algorithmns performance */
  bool visualizeAstar_ = false;

  /** \brief Visualization speed of astar search, num of seconds to show each vertex */
  double visualizeAstarSpeed_ = 0.1;

  /** \brief SPARS parameter for dense graph connection distance as a fraction of max. extent */
  double denseDeltaFraction_ = 0.05;

  /** \brief Maximum visibility range for nodes in the graph as a fraction of maximum extent. */
  double sparseDeltaFraction_ = 0.25;

  /** \brief Multiply this number by the dimension of the state space to choose how much sampling to perform */
  std::size_t nearSamplePointsMultiple_ = 2;

  /** \brief The stretch factor in terms of graph spanners for SPARS to check against */
  double stretchFactor_ = 3.0;

  /** \brief How much the popularity of a node can cause its cost-to-go heuristic to be underestimated */
  double percentMaxExtentUnderestimate_ = 0.01;

  /** \brief Change verbosity levels */
  bool checksVerbose_ = false;
  bool disjointVerbose_ = true;
  bool fourthCheckVerbose_ = true;

  /** \brief Show the sparse graph being generated */
  bool visualizeSparsGraph_ = false;
  bool visualizeQualityCriteria_ = true;
  double visualizeSparsGraphSpeed_ = 0.0;
  bool visualizeDatabaseVertices_ = true;
  bool visualizeDatabaseEdges_ = true;
  bool visualizeDenseRepresentatives_ = false;  // TODO use this
  bool visualizeNodePopularity_ = false;

  /** \brief Method for ordering of vertex insertion */
  int sparseCreationInsertionOrder_ = 0;

  /** \brief For statistics */
  int numGraphGenerations_ = 0;
  int numRandSamplesAdded_ = 0;
  int numSamplesAddedForFourthCriteria_ = 0;
  /** \brief Same discretization used in Discretizer */
  double discretization_;

  bool testingBool_;
};  // end of class SparseDB

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_SPARSEDB_
