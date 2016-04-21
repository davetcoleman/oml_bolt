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

/* Author: Dave Coleman */

#ifndef OMPL_GEOMETRIC_PLANNERS_EXPERIENCE_BOLT_BOLT_RETRIEVE_REPAIR_
#define OMPL_GEOMETRIC_PLANNERS_EXPERIENCE_BOLT_BOLT_RETRIEVE_REPAIR_

#include <ompl/geometric/planners/PlannerIncludes.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/geometric/PathSimplifier.h>
#include <ompl/tools/bolt/DenseDB.h>
#include <ompl/tools/bolt/Visualizer.h>

// Boost
#include <boost/lambda/bind.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/function.hpp>
#include <boost/thread.hpp>

namespace ompl
{
namespace tools
{
namespace bolt
{
/// @cond IGNORE
/** \brief Forward declaration of ompl::base::BoltRetrieveRepair */
OMPL_CLASS_FORWARD(BoltRetrieveRepair);
OMPL_CLASS_FORWARD(DenseDB);
/// @endcond

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
    BoltRetrieveRepair* parent_;

  public:
    /**
     * Construct a visitor for a given search.
     * \param goal  goal vertex of the search
     */
    CustomAstarVisitor(SparseVertex goal, BoltRetrieveRepair* parent);

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

/** \class ompl::base::BoltRetrieveRepairPtr
    \brief A boost shared pointer wrapper for ompl::base::BoltRetrieveRepair */

/**
   @anchor BoltRetrieveRepair
   @par Short description
   Bolt is an experience-based planning framework that learns to reduce computation time
   required to solve high-dimensional planning problems in varying environments.
   @par External documentation
   Berenson, Dmitry, Pieter Abbeel, and Ken Goldberg: A robot path planning framework that learns from experience, in
   <em>Robotics and Automation (ICRA), 2012 IEEE International Conference on. IEEE</em>, 2012.
   David Coleman, Ioan A. Sucan, Mark Moll, Kei Okada, Nikolaus Correll, "Experience-Based Planning with Sparse Roadmap
   Spanners"
   <a href="http://arxiv.org/pdf/1410.1950.pdf">[PDF]</a>
*/

/** \brief The Bolt Framework's Retrieve-Repair component */
class BoltRetrieveRepair : public base::Planner
{
public:
  /** \brief Constructor */
  BoltRetrieveRepair(const base::SpaceInformationPtr &si, const DenseDBPtr &denseDB);

  virtual ~BoltRetrieveRepair(void);

  /** \brief Wrapper function to show good user feedback while smoothing a path */
  bool simplifyPath(geometric::PathGeometric &path, const base::PlannerTerminationCondition &ptc);

  /** \brief Get information about the exploration data structure the planning from scratch motion planner used. */
  virtual void getPlannerData(base::PlannerData &data) const;

  /**
   *  \brief Get debug information about the top recalled paths that were chosen for further filtering
   *  \return data - vector of PlannerData objects that each hold a single path
   */
  // const std::vector<PathGeometric> &getLastRecalledNearestPaths() const;

  /**
   *  \brief Get debug information about the top recalled paths that were chosen for further filtering
   *  \return chosenID - the index of the PlannerData object that was chosen for repair
   */
  // std::size_t getLastRecalledNearestPathChosen() const;

  /**
   * \brief Get the chosen path used from database for repair
   * \return PlannerData of chosen path
   */
  const geometric::PathGeometric &getChosenRecallPath() const;

  /** \brief Main entry function for finding a path plan */
  virtual base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc);

  /** \brief Clear memory */
  virtual void clear(void);

  /**
   * \brief Pass a pointer of the database from the bolt framework
   */
  void setExperienceDB(const DenseDBPtr &denseDB);

  /** \brief Setup function */
  virtual void setup(void);

  /** \brief Optionally smooth retrieved and repaired paths from database */
  void enableSmoothing(bool enable)
  {
    smoothingEnabled_ = enable;
  }

  /**
   * \brief Search the roadmap for the best path close to the given start and goal states that is valid
   * \param start
   * \param goal
   * \param geometricSolution - the resulting path
   * \return
   */
  bool getPathOffGraph(const base::State *start, const base::State *goal, geometric::PathGeometric &geometricSolution,
                       const base::PlannerTerminationCondition &ptc);

  /** \brief Clear verticies not on the specified level */
  bool removeVerticesNotOnLevel(std::vector<bolt::DenseVertex> &graphNeighborhood, int level);

  /**
   * \brief Convert astar results to correctly ordered path
   * \param vertexPath - in reverse
   * \param start - actual start that is probably not included in new path
   * \param goal - actual goal that is probably not included in new path
   * \param path - returned solution
   * \return true on success
   */
  bool convertVertexPathToStatePath(std::vector<bolt::DenseVertex> &vertexPath, const base::State *actualStart,
                                    const base::State *actualGoal, geometric::PathGeometric &geometricSolution);

  /**
   * \brief Finds nodes in the graph near state NOTE: note tested for visibility
   * \param state - vertex to find neighbors around
   * \param graphNeighborhood - result vector
   * \param requiredLevel - if -1, allows states from all levels, otherwise only returns states from a certain level
   * \return false is no neighbors found
   */
  bool findGraphNeighbors(const base::State *state, std::vector<bolt::DenseVertex> &graphNeighborhood,
                          int requiredLevel = -1);

  /** \brief Check if there exists a solution, i.e., there exists a pair of milestones such that the
   *   first is in \e start and the second is in \e goal, and the two milestones are in the same
   *   connected component. If a solution is found, the path is saved.
   * \param debug - whether to show the failure points
   * \param feedbackStartFailed - if getPathOnGraph returns false, this flag determines if the start or goal node failed
   * to conenct
   */
  bool getPathOnGraph(const std::vector<bolt::DenseVertex> &candidateStarts,
                      const std::vector<bolt::DenseVertex> &candidateGoals, const base::State *actualStart,
                      const base::State *actualGoal, geometric::PathGeometric &geometricSolution,
                      const base::PlannerTerminationCondition &ptc, bool debug, bool &feedbackStartFailed);

  /**
   * \brief Repeatidly search through graph for connection then check for collisions then repeat
   * \return true if a valid path is found
   */
  bool lazyCollisionSearch(const bolt::DenseVertex &start, const bolt::DenseVertex &goal,
                           const base::State *actualStart, const base::State *actualGoal,
                           geometric::PathGeometric &geometricSolution, const base::PlannerTerminationCondition &ptc);

  /** \brief Check recalled path for collision and disable as needed */
  bool lazyCollisionCheck(std::vector<bolt::DenseVertex> &vertexPath, const base::PlannerTerminationCondition &ptc);

  /** \brief Test if the passed in random state can connect to a nearby vertex in the graph */
  bool canConnect(const base::State *randomState, const base::PlannerTerminationCondition &ptc);

  bool astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex> &vertexPath);

  /** \brief Distance between two states with special bias using popularity */
  double astarHeuristic(const SparseVertex a, const SparseVertex b) const;

  /** \brief Get class for managing various visualization features */
  base::VisualizerPtr getVisual()
  {
    return visual_;
  }

  SparseDBPtr getSparseDB()
  {
    return sparseDB_;
  }

  /** \brief Custom A* visitor statistics */
  void recordNodeOpened()  // discovered
  {
    numNodesOpened_++;
  }
  void recordNodeClosed()  // examined
  {
    numNodesClosed_++;
  }


protected:
  /**
   * \brief Count the number of states along the discretized path that are in collision
   *        Note: This is kind of an ill-defined score though. It depends on the resolution of collision checking.
   *        I am more inclined to try to compute the percent of the length of the motion that is valid.
   *        That could go in SpaceInformation, as a utility function.
   */
  std::size_t checkMotionScore(const base::State *s1, const base::State *s2) const;

  /** \brief Free the memory allocated by this planner */
  void freeMemory(void);

  /** \brief The database of motions to search through */
  DenseDBPtr denseDB_;
  SparseDBPtr sparseDB_;

  /** \brief Class for managing various visualization features */
  base::VisualizerPtr visual_;

  /** \brief Save the recalled path before smoothing for introspection later */
  boost::shared_ptr<geometric::PathGeometric> originalSolutionPath_;

  /** \brief The instance of the path simplifier */
  geometric::PathSimplifierPtr path_simplifier_;

  /** \brief Optionally smooth retrieved and repaired paths from database */
  bool smoothingEnabled_;

  /** \brief Used by getPathOffGraph */
  std::vector<bolt::DenseVertex> startVertexCandidateNeighbors_;
  std::vector<bolt::DenseVertex> goalVertexCandidateNeighbors_;

  /** \brief Astar statistics */
  std::size_t numNodesOpened_ = 0;
  std::size_t numNodesClosed_ = 0;

public:

  /** \brief Various options for visualizing the algorithmns performance */
  bool visualizeAstar_ = false;

  /** \brief Visualization speed of astar search, num of seconds to show each vertex */
  double visualizeAstarSpeed_ = 0.1;

  /** \brief Output user feedback to console */
  bool verbose_;

  /** \brief Visualize original solution from graph before smoothing */
  bool visualizeRawTrajectory_;

  int numStartGoalStatesAddedToDense_;
};
}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif
