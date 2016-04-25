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

// OMPL
#include <ompl/tools/bolt/BoltRetrieveRepair.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/base/goals/GoalState.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include <ompl/tools/config/SelfConfig.h>
#include <ompl/util/Console.h>

// Boost
#include <boost/thread.hpp>
#include <boost/foreach.hpp>

// C++
#include <limits>

namespace og = ompl::geometric;
namespace ob = ompl::base;
namespace ot = ompl::tools;
namespace otb = ompl::tools::bolt;

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

otb::CustomAstarVisitor::CustomAstarVisitor(SparseVertex goal, BoltRetrieveRepair *parent)
  : goal_(goal), parent_(parent)
{
}

void otb::CustomAstarVisitor::discover_vertex(SparseVertex v, const SparseGraph &) const
{
  // Statistics
  parent_->recordNodeOpened();

  if (parent_->visualizeAstar_)
    parent_->getVisual()->viz4State(parent_->getSparseDB()->getSparseState(v), /*small*/ 1, tools::GREEN, 1);
}

void otb::CustomAstarVisitor::examine_vertex(SparseVertex v, const SparseGraph &) const
{
  // Statistics
  parent_->recordNodeClosed();

  if (parent_->visualizeAstar_)
  {
    parent_->getVisual()->viz4State(parent_->getSparseDB()->getSparseState(v), /*large*/ 3, tools::BLACK, 1);
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
BoltRetrieveRepair::BoltRetrieveRepair(const base::SpaceInformationPtr &si, const DenseDBPtr &denseDB)
  : base::Planner(si, "Bolt_Retrieve_Repair")
  , denseDB_(denseDB)
  , smoothingEnabled_(true)
  , verbose_(true)
  , numStartGoalStatesAddedToDense_(0)
{
  // Copy in needed objects
  sparseDB_ = denseDB_->getSparseDB();
  visual_ = denseDB_->visual_;

  specs_.approximateSolutions = true;
  specs_.directed = true;

  path_simplifier_.reset(new geometric::PathSimplifier(si_));
}

BoltRetrieveRepair::~BoltRetrieveRepair(void)
{
  freeMemory();
}

void BoltRetrieveRepair::clear(void)
{
  Planner::clear();
  freeMemory();
}

void BoltRetrieveRepair::setExperienceDB(const DenseDBPtr &denseDB)
{
  denseDB_ = denseDB;
}

void BoltRetrieveRepair::setup(void)
{
  Planner::setup();
}

void BoltRetrieveRepair::freeMemory(void)
{
}

base::PlannerStatus BoltRetrieveRepair::solve(const base::PlannerTerminationCondition &ptc)
{
  OMPL_INFORM("BoltRetrieveRepair::solve()");

  bool solved = false;
  double approxdif = std::numeric_limits<double>::infinity();

  // Check if the database is empty
  if (sparseDB_->isEmpty())
  {
    OMPL_INFORM("Sparse experience database is empty so unable to run BoltRetrieveRepair algorithm.");

    return base::PlannerStatus::ABORT;
  }

  // Restart the Planner Input States so that the first start and goal state can be fetched
  pis_.restart();  // PlannerInputStates

  // Get a single start and goal state
  if (verbose_)
    OMPL_INFORM("Getting OMPL start and goal state");
  const base::State *startState = pis_.nextStart();  // PlannerInputStates
  const base::State *goalState = pis_.nextGoal(ptc);

  if (startState == nullptr)
  {
    OMPL_ERROR("No start state found");
    return base::PlannerStatus::ABORT;
  }

  if (goalState == nullptr)
  {
    OMPL_ERROR("No goal state found");
    return base::PlannerStatus::ABORT;
  }

  // Error check task planning
  if (denseDB_->getUseTaskPlanning())
  {
    if (denseDB_->getTaskLevel(startState) != 0)
    {
      OMPL_ERROR("solve: start level is %u", denseDB_->getTaskLevel(startState));
      exit(-1);
    }
    if (denseDB_->getTaskLevel(goalState) != 2)
    {
      OMPL_ERROR("solve: goal level is %u", denseDB_->getTaskLevel(goalState));
      exit(-1);
    }
  }

  // Create solution path as pointer so memory is not unloaded
  ob::PathPtr pathSolutionBase(new og::PathGeometric(si_));
  og::PathGeometric &geometricSolution = static_cast<og::PathGeometric &>(*pathSolutionBase);

  // Search
  if (!getPathOffGraph(startState, goalState, geometricSolution, ptc))
  {
    OMPL_WARN("BoltRetrieveRepair::solve() No near start or goal found");
    return base::PlannerStatus::TIMEOUT;  // The planner failed to find a solution
  }

  if (verbose_)
    OMPL_INFORM("getPathOffGraph() found a solution of size %d", geometricSolution.getStateCount());

  // Save this for future debugging
  originalSolutionPath_.reset(new geometric::PathGeometric(geometricSolution));

  // All save trajectories should be at least 1 state long, then we append the start and goal states, for min of 3
  assert(geometricSolution.getStateCount() >= 3);

  // Smooth the result
  if (smoothingEnabled_)
  {
    simplifyPath(geometricSolution, ptc);
  }

  // Add more points to path
  geometricSolution.interpolate();

  // Finished
  approxdif = 0;
  bool approximate = false;

  // Save solution
  pdef_->addSolutionPath(pathSolutionBase, approximate, approxdif, getName());
  solved = true;

  if (verbose_)
    OMPL_INFORM("  Finished BoltRetrieveRepair.solve()");
  return base::PlannerStatus(solved, approximate);
}

bool BoltRetrieveRepair::simplifyPath(og::PathGeometric &path, const base::PlannerTerminationCondition &ptc)
{
  if (verbose_)
    OMPL_INFORM("BoltRetrieveRepair: Simplifying solution (smoothing)...");

  time::point simplifyStart = time::now();
  std::size_t numStates = path.getStateCount();
  path_simplifier_->simplify(path, ptc);
  double simplifyTime = time::seconds(time::now() - simplifyStart);

  OMPL_INFORM("BoltRetrieveRepair: Path simplification took %f seconds and removed %d states", simplifyTime,
              numStates - path.getStateCount());
  return true;
}

void BoltRetrieveRepair::getPlannerData(base::PlannerData &data) const
{
  if (verbose_)
    OMPL_INFORM("BoltRetrieveRepair getPlannerData");

  for (std::size_t j = 1; j < originalSolutionPath_->getStateCount(); ++j)
  {
    data.addEdge(base::PlannerDataVertex(originalSolutionPath_->getState(j - 1)),
                 base::PlannerDataVertex(originalSolutionPath_->getState(j)));
  }
}

const geometric::PathGeometric &BoltRetrieveRepair::getChosenRecallPath() const
{
  return *originalSolutionPath_;
}

std::size_t BoltRetrieveRepair::checkMotionScore(const base::State *s1, const base::State *s2) const
{
  int segmentCount = si_->getStateSpace()->validSegmentCount(s1, s2);

  std::size_t invalidStatesScore = 0;  // count number of interpolated states in collision

  // temporary storage for the checked state
  base::State *test = si_->allocState();

  // Linerarly step through motion between state 0 to state 1
  double iteration_step = 1.0 / double(segmentCount);
  for (double location = 0.0; location <= 1.0; location += iteration_step)
  {
    si_->getStateSpace()->interpolate(s1, s2, location, test);

    if (!si_->isValid(test))
    {
      // OMPL_DEBUG("Found INVALID location between states at gradient %f", location);
      invalidStatesScore++;
    }
    else
    {
      // OMPL_DEBUG("Found valid location between states at gradient %f", location);
    }
  }
  si_->freeState(test);

  return invalidStatesScore;
}

bool BoltRetrieveRepair::getPathOffGraph(const base::State *start, const base::State *goal,
                                         og::PathGeometric &geometricSolution,
                                         const base::PlannerTerminationCondition &ptc)
{
  // Attempt to connect to graph x times, because if it failes we start adding samples
  std::size_t maxAttempts = 2;
  std::size_t attempt = 0;
  for (; attempt < maxAttempts; ++attempt)
  {
    if (verbose_)
      OMPL_INFORM("Starting getPathOffGraph attempt %u", attempt);

    // Get neighbors near start and goal. Note: potentially they are not *visible* - will test for this later

    // Start
    int level = 0;  // denseDB_->getTaskLevel(start);
    if (verbose_)
      OMPL_INFORM("  Looking for a node near the problem start on level %i", level);
    if (!findGraphNeighbors(start, startVertexCandidateNeighbors_, level))
    {
      if (verbose_)
        OMPL_INFORM("No graph neighbors found for start");
      return false;
    }
    if (verbose_)
      OMPL_INFORM("  Found %d nodes near start", startVertexCandidateNeighbors_.size());

    // Goal
    level = 0;  // denseDB_->getTaskLevel(goal);
    if (verbose_)
      OMPL_INFORM("  Looking for a node near the problem goal on level %i", level);
    if (!findGraphNeighbors(goal, goalVertexCandidateNeighbors_, level))
    {
      if (verbose_)
        OMPL_INFORM("No graph neighbors found for goal");
      return false;
    }
    if (verbose_)
      OMPL_INFORM("    Found %d nodes near goal", goalVertexCandidateNeighbors_.size());

    // Get paths between start and goal
    bool feedbackStartFailed;
    bool result = getPathOnGraph(startVertexCandidateNeighbors_, goalVertexCandidateNeighbors_, start, goal,
                                 geometricSolution, ptc, /*debug*/ false, feedbackStartFailed);

    // Error check
    if (!result)
    {
      OMPL_WARN("getPathOffGraph(): BoltRetrieveRepair returned FALSE for getPathOnGraph. Trying again in debug mode");

      // Run getPathOnGraph again in debug mode
      getPathOnGraph(startVertexCandidateNeighbors_, goalVertexCandidateNeighbors_, start, goal, geometricSolution, ptc,
                     /*debug*/ true, feedbackStartFailed);

      // Add the point that failed
      if (feedbackStartFailed)  // start state failed
      {
        // Create the vertex then connect to denseDB
        // TODO(davetcoleman): how to prevent from adding same state twice?
        DenseVertex denseV = denseDB_->addVertex(si_->cloneState(start), QUALITY);
        denseDB_->connectNewVertex(denseV);
      }
      else  // goal state failed
      {
        // Create the vertex then connect to denseDB
        // TODO(davetcoleman): how to prevent from adding same state twice?
        DenseVertex denseV = denseDB_->addVertex(si_->cloneState(goal), QUALITY);
        denseDB_->connectNewVertex(denseV);
      }
      numStartGoalStatesAddedToDense_++;  // for later analysis

      OMPL_INFORM("Re-creating the spars graph");
      sparseDB_->createSPARS();

      // std::cout << "Shutting down for debugging " << std::endl;
      // exit(-1);
    }
    else
      break;  // success, continue on
  }

  // Did we finally get it?
  if (attempt >= maxAttempts)
  {
    return false;
  }

  // All save trajectories should be at least 1 state long, then we append the start and goal states, for min of 3
  assert(geometricSolution.getStateCount() >= 3);

  // Debug output
  if (false)
  {
    for (std::size_t i = 0; i < geometricSolution.getStateCount(); ++i)
    {
      if (verbose_)
        OMPL_INFORM("  getPathOffGraph(): Adding state %f to plannerData", i);
      si_->printState(geometricSolution.getState(i), std::cout);
    }
  }

  return true;
}

bool BoltRetrieveRepair::getPathOnGraph(const std::vector<SparseVertex> &candidateStarts,
                                        const std::vector<SparseVertex> &candidateGoals, const base::State *actualStart,
                                        const base::State *actualGoal, og::PathGeometric &geometricSolution,
                                        const base::PlannerTerminationCondition &ptc, bool debug,
                                        bool &feedbackStartFailed)
{
  bool foundValidStart = false;
  bool foundValidGoal = false;

  // Try every combination of nearby start and goal pairs
  BOOST_FOREACH (SparseVertex start, candidateStarts)
  {
    if (actualStart == sparseDB_->getSparseStateConst(start))
    {
      OMPL_ERROR("Comparing same start state");
      exit(-1);  // just curious if this ever happens, no need to actually exit
      continue;
    }

    // Check if this start is visible from the actual start
    if (!si_->checkMotion(actualStart, sparseDB_->getSparseStateConst(start)))
    {
      if (verbose_)
      {
        OMPL_WARN("FOUND START CANDIDATE THAT IS NOT VISIBLE ");
      }
      if (debug)
      {
        visual_->viz4State(sparseDB_->getSparseStateConst(start), /*large*/3, tools::RED, 1);
        visual_->viz4Edge(actualStart, sparseDB_->getSparseStateConst(start), 100);
        visual_->viz4Trigger();
        usleep(0.1 * 1000000);
      }
      continue;  // this is actually not visible
    }
    foundValidStart = true;

    BOOST_FOREACH (SparseVertex goal, candidateGoals)
    {
      if (actualGoal == sparseDB_->getSparseStateConst(goal))
      {
        OMPL_ERROR("Comparing same goal state");
        continue;
      }

      if (verbose_)
        OMPL_INFORM("    foreach_goal: Checking motion from  %d to %d", actualGoal,
                    sparseDB_->getSparseStateConst(goal));

      // Check if our planner is out of time
      if (ptc == true)
      {
        OMPL_DEBUG("getPathOnGraph function interrupted because termination condition is true.");
        return false;
      }

      // Check if this goal is visible from the actual goal
      if (!si_->checkMotion(actualGoal, sparseDB_->getSparseStateConst(goal)))
      {
        if (verbose_)
        {
          OMPL_WARN("FOUND GOAL CANDIDATE THAT IS NOT VISIBLE! ");
        }

        if (debug)
        {
          visual_->viz4State(sparseDB_->getSparseStateConst(goal), /*small*/1, tools::RED, 1);
          visual_->viz4Edge(actualGoal, sparseDB_->getSparseStateConst(goal), 100);
          visual_->viz4Trigger();
          usleep(0.1 * 1000000);
        }

        continue; // this is actually not visible
      }
      foundValidGoal = true;

      // Repeatidly search through graph for connection then check for collisions then repeat
      if (lazyCollisionSearch(start, goal, actualStart, actualGoal, geometricSolution, ptc))
      {
        // All save trajectories should be at least 1 state long, then we append the start and goal states, for
        // min of 3
        assert(geometricSolution.getStateCount() >= 3);

        // Found a path
        return true;
      }
      else
      {
        // Did not find a path
        OMPL_INFORM("Did not find a path, looking for other start/goal combinations ");
      }

    }  // foreach
  }    // foreach

  if (foundValidStart && foundValidGoal)
  {
    OMPL_ERROR("Unexpected condition - both a valid start and goal were found but still no path found. TODO ");
    exit(-1);
  }

  if (foundValidStart && !foundValidGoal)
  {
    OMPL_WARN("Unable to connect GOAL state to graph");
    feedbackStartFailed = false;  // it was the goal state that failed us
  }
  else
  {
    OMPL_WARN("Unable to connect START state to graph");
    feedbackStartFailed = true;  // the start state failed us
  }

  return false;
}

bool BoltRetrieveRepair::lazyCollisionSearch(const SparseVertex &start, const SparseVertex &goal,
                                             const base::State *actualStart, const base::State *actualGoal,
                                             og::PathGeometric &geometricSolution,
                                             const base::PlannerTerminationCondition &ptc)
{
  // Vector to store candidate paths in before they are converted to PathPtrs
  std::vector<SparseVertex> vertexPath;

  // Make sure that the start and goal aren't so close together that they find the same vertex
  if (start == goal)
  {
    if (verbose_)
      OMPL_INFORM("    Start equals goal, creating simple solution ");

    // There are only three verticies in this path - start, middle, goal
    vertexPath.push_back(start);

    convertVertexPathToStatePath(vertexPath, actualStart, actualGoal, geometricSolution);
    return true;
  }

  // Error check all states are non-nullptr
  assert(actualStart);
  assert(actualGoal);
  assert(sparseDB_->getSparseStateConst(start));
  assert(sparseDB_->getSparseStateConst(goal));

  // Check that our states are on the same connected component
  // TODO: in the future the graph should always just be fully connected
  // so perhaps this check would not be necessary
  if (!sparseDB_->sameComponent(start, goal))
  {
    OMPL_WARN("Found start and goal states that are on different connected components!");
    return false;
  }

  // Visualize start vertex
  const bool visualize = false;
  if (visualize)
  {
    OMPL_INFORM("viz start -----------------------------");
    visual_->viz5State(sparseDB_->getSparseStateConst(start), /*Medium, translucent outline*/4, tools::PURPLE, 1);
    visual_->viz5Edge(actualStart, sparseDB_->getSparseStateConst(start), 30);
    visual_->viz5Trigger();
    usleep(5 * 1000000);

    // Visualize goal vertex
    OMPL_INFORM("viz goal ------------------------------");
    visual_->viz5State(sparseDB_->getSparseStateConst(goal), /*Medium, translucent outline*/4, tools::PURPLE, 1);
    visual_->viz5Edge(actualGoal, sparseDB_->getSparseStateConst(goal), 0);
    visual_->viz5Trigger();
    usleep(5 * 1000000);
  }

  // Keep looking for paths between chosen start and goal until one is found that is valid,
  // or no further paths can be found between them because of disabled edges
  // this is necessary for lazy collision checking i.e. rerun after marking invalid edges we found
  while (true)
  {
    if (verbose_)
      OMPL_INFORM("  AStar: looking for path through graph between start and goal");

    // Check if our planner is out of time
    if (ptc == true)
    {
      OMPL_DEBUG("lazyCollisionSearch: function interrupted because termination condition is true.");
      return false;
    }

    // Attempt to find a solution from start to goal
    if (!astarSearch(start, goal, vertexPath))
    {
      OMPL_INFORM("        unable to construct solution between start and goal using astar");

      // no path found what so ever
      return false;
    }

    if (verbose_)
    {
      OMPL_INFORM("        Has at least a partial solution, maybe exact solution");
      OMPL_INFORM("        Solution has %d vertices", vertexPath.size());
    }

    // Check if all the points in the potential solution are valid
    if (lazyCollisionCheck(vertexPath, ptc))
    {
      if (verbose_)
      {
        OMPL_INFORM("  Lazy collision check returned valid ");
      }

      // the path is valid, we are done!
      convertVertexPathToStatePath(vertexPath, actualStart, actualGoal, geometricSolution);
      return true;
    }
    // else, loop with updated graph that has the invalid edges/states disabled
  }  // end while

  // we never found a valid path
  return false;
}

bool BoltRetrieveRepair::lazyCollisionCheck(std::vector<SparseVertex> &vertexPath,
                                            const base::PlannerTerminationCondition &ptc)
{
  OMPL_INFORM("Starting lazy collision checking");

  bool hasInvalidEdges = false;

  // Initialize
  SparseVertex fromVertex = vertexPath[0];
  SparseVertex toVertex;

  // Loop through every pair of states and make sure path is valid.
  for (std::size_t toID = 1; toID < vertexPath.size(); ++toID)
  {
    // Increment location on path
    toVertex = vertexPath[toID];

    // Check if our planner is out of time
    if (ptc)
    {
      OMPL_DEBUG("Lazy collision check function interrupted because termination condition is true.");
      return false;
    }

    SparseEdge thisEdge = boost::edge(fromVertex, toVertex, sparseDB_->g_).first;

    // Has this edge already been checked before?
    if (sparseDB_->edgeCollisionStatePropertySparse_[thisEdge] == NOT_CHECKED)
    {
      // Check path between states
      if (!si_->checkMotion(sparseDB_->getSparseStateConst(fromVertex), sparseDB_->getSparseStateConst(toVertex)))
      {
        // Path between (from, to) states not valid, disable the edge
        OMPL_INFORM("  DISABLING EDGE from vertex %f to vertex %f", fromVertex, toVertex);

        // Disable edge
        sparseDB_->edgeCollisionStatePropertySparse_[thisEdge] = IN_COLLISION;
      }
      else
      {
        // Mark edge as free so we no longer need to check for collision
        sparseDB_->edgeCollisionStatePropertySparse_[thisEdge] = FREE;
      }
    }

    // Check final result
    if (sparseDB_->edgeCollisionStatePropertySparse_[thisEdge] == IN_COLLISION)
    {
      // Remember that this path is no longer valid, but keep checking remainder of path edges
      hasInvalidEdges = true;
    }

    // switch vertex focus
    fromVertex = toVertex;
  }

  if (verbose_)
    OMPL_INFORM("  Done lazy collision checking");

  // Only return true if nothing was found invalid
  return !hasInvalidEdges;
}

bool BoltRetrieveRepair::findGraphNeighbors(const base::State *state, std::vector<SparseVertex> &graphNeighborhood,
                                            int requiredLevel)
{
  // Benchmark runtime
  time::point startTime = time::now();

  // Reset
  graphNeighborhood.clear();

  // Setup search by getting a non-const version of the focused state
  const std::size_t threadID = 0;
  base::State *stateCopy = si_->cloneState(state);
  sparseDB_->getSparseState(sparseDB_->queryVertices_[threadID]) = stateCopy;

  // Search
  double find_nearest_k_neighbors;
  if (si_->getStateSpace()->getDimension() == 3)
    find_nearest_k_neighbors = 10;
  else
    find_nearest_k_neighbors = 30;
  sparseDB_->nn_->nearestK(sparseDB_->queryVertices_[threadID], find_nearest_k_neighbors, graphNeighborhood);

  // Reset
  sparseDB_->getSparseState(sparseDB_->queryVertices_[threadID]) = nullptr;

  // Filter neighbors based on level
  if (requiredLevel >= 0)
  {
    removeVerticesNotOnLevel(graphNeighborhood, requiredLevel);
  }

  // Check if no neighbors found
  bool result = graphNeighborhood.size();

  // Benchmark runtime
  double duration = time::seconds(time::now() - startTime);

  if (verbose_)
    OMPL_INFORM("   - findGraphNeighbors() took %f seconds", duration);

  // Free memory
  si_->getStateSpace()->freeState(stateCopy);

  return result;
}

bool BoltRetrieveRepair::removeVerticesNotOnLevel(std::vector<SparseVertex> &graphNeighborhood, int level)
{
  std::size_t original_size = graphNeighborhood.size();

  // Remove edges based on layer
  for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  {
    const SparseVertex &nearVertex = graphNeighborhood[i];

    // Make sure state is on correct level
    if (denseDB_->getTaskLevel(nearVertex) != static_cast<std::size_t>(level))
    {
      if (verbose_)
        std::cout << "      Skipping neighbor " << nearVertex << ", i=" << i
                  << ", because wrong level: " << denseDB_->getTaskLevel(nearVertex) << ", desired level: " << level
                  << std::endl;
      graphNeighborhood.erase(graphNeighborhood.begin() + i);
      i--;
      continue;
    }
  }

  if (verbose_)
    OMPL_INFORM("    removeVerticesNotOnLevel(): states require level %d, removed: %u, remaining: %u", level,
                original_size - graphNeighborhood.size(), graphNeighborhood.size());

  return true;
}

bool BoltRetrieveRepair::convertVertexPathToStatePath(std::vector<SparseVertex> &vertexPath,
                                                      const base::State *actualStart, const base::State *actualGoal,
                                                      og::PathGeometric &geometricSolution)
{
  // Ensure the input path is not empty
  if (!vertexPath.size())
    return false;

  // og::PathGeometric *pathGeometric = new og::PathGeometric(si_);

  // Add original start if it is different than the first state
  if (actualStart != sparseDB_->getSparseStateConst(vertexPath.back()))
  {
    geometricSolution.append(actualStart);
  }

  // Error check that no consequtive verticies are the same
  /*std::cout << "BoltRetrieveRepair.Vertices: ";
  for (std::size_t i = vertexPath.size(); i > 0; --i)
  {
      std::cout << vertexPath[i - 1] << ", ";
  }
  std::cout << std::endl;*/

  // Reverse the vertexPath and convert to state path
  for (std::size_t i = vertexPath.size(); i > 0; --i)
  {
    geometricSolution.append(sparseDB_->getSparseStateConst(vertexPath[i - 1]));

    // Add the edge status
    if (i > 1)  // skip the last vertex (its reversed)
    {
      // Error check that no consequtive verticies are the same
      if (vertexPath[i - 1] == vertexPath[i - 2])
      {
        OMPL_ERROR("Found repeated vertices %u to %u on index %u", vertexPath[i - 1], vertexPath[i - 2], i);
        exit(-1);
      }

      SparseEdge edge = boost::edge(vertexPath[i - 1], vertexPath[i - 2], sparseDB_->g_).first;

      /* This functionality has moved to DenseDB
      if (sparseDB_->getPopularityBiasEnabled())
      {
          // reduce cost of this edge because it was just used
          static const double REDUCTION_AMOUNT = 10;
          std::cout << "Previous edge weight for " << i << " of edge " << edge << std::endl;
          std::cout << "    is " << sparseDB_->edgeWeightProperty_[edge];
          sparseDB_->edgeWeightProperty_[edge] =
              std::max(sparseDB_->edgeWeightProperty_[edge] - REDUCTION_AMOUNT, 0.0);
          std::cout << " new: " << sparseDB_->edgeWeightProperty_[edge] << std::endl;
      }
      */

      // Check if any edges in path are not free (then it an approximate path)
      if (sparseDB_->edgeCollisionStatePropertySparse_[edge] == IN_COLLISION)
      {
        OMPL_ERROR("Found invalid edge / approximate solution - how did this happen?");
      }
      else if (sparseDB_->edgeCollisionStatePropertySparse_[edge] == NOT_CHECKED)
      {
        OMPL_ERROR("A chosen path has an edge that has not been checked for collision. This should not happen");
      }
    }
  }

  // Add original goal if it is different than the last state
  if (actualGoal != sparseDB_->getSparseStateConst(vertexPath.front()))
  {
    geometricSolution.append(actualGoal);
  }

  // if (sparseDB_->getPopularityBiasEnabled())
  // {
  //     // Ensure graph doesn't get too popular
  //     sparseDB_->normalizeGraphEdgeWeights();

  //     // Mark as needing to be saved to file
  //     sparseDB_->graphUnsaved_ = true;
  // }

  return true;
}

bool BoltRetrieveRepair::canConnect(const base::State *randomState, const base::PlannerTerminationCondition &ptc)
{
  std::vector<SparseVertex> candidateNeighbors;

  // Find neighbors to rand state
  OMPL_INFORM("Looking for a node near the random state");
  if (!findGraphNeighbors(randomState, candidateNeighbors))
  {
    OMPL_INFORM("No graph neighbors found for randomState");
    return false;
  }
  OMPL_INFORM("  Found %d nodes near randomState", candidateNeighbors.size());

  // Try every combination of nearby start and goal pairs
  std::size_t count = 0;
  BOOST_FOREACH (SparseVertex nearState, candidateNeighbors)
  {
    const base::State *s1 = randomState;
    const base::State *s2 = sparseDB_->getSparseStateConst(nearState);

    // Check if this nearState is visible from the random state
    if (!si_->checkMotion(s1, s2))
    {
      OMPL_WARN("NEIGHBOR %u NOT VISIBLE ", count++);

      if (false)
      {
        visual_->viz5State(s2, /*medium*/2, tools::BLUE, 1);
        visual_->viz5Edge(s1, s2, 100);
        visual_->viz5Trigger();
        usleep(1 * 1000000);
      }

      // Optional Debug
      if (false)
      {
        std::cout << "checking path " << std::endl;
        std::vector<base::State *> states;
        unsigned int count = si_->getStateSpace()->validSegmentCount(s1, s2);
        // std::cout << "count: " << count << std::endl;

        bool endpoints = false;
        bool alloc = true;
        si_->getMotionStates(s1, s2, states, count, endpoints, alloc);
        // std::cout << "state size: " << states.size() << std::endl;

        BOOST_FOREACH (base::State *interState, states)
        {
          // Check if our planner is out of time
          if (ptc == true)
          {
            OMPL_INFORM("Quit requested");
            return false;
          }

          if (!si_->isValid(interState))
          {
            visual_->viz5State(interState, /*mode=*/3, tools::RED, 1);
            visual_->viz5Trigger();
            usleep(1 * 1000000);
          }
          else
          {
            // visual_->viz5State(interState, /*mode=*/1, 1); // GREEN
          }
        }
      }
    }
    else
    {
      OMPL_INFORM("Has connection");
      return true;
    }
  }
  return false;
}

bool BoltRetrieveRepair::astarSearch(const SparseVertex start, const SparseVertex goal, std::vector<SparseVertex> &vertexPath)
{
  // Hold a list of the shortest path parent to each vertex
  SparseVertex *vertexPredecessors = new SparseVertex[sparseDB_->getNumVertices()];
  // boost::vector_property_map<SparseVertex> vertexPredecessors(getNumVertices());

  bool foundGoal = false;
  double *vertexDistances = new double[sparseDB_->getNumVertices()];

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
    boost::astar_search(sparseDB_->g_,                                                 // graph
                        start,                                                        // start state
                        boost::bind(&otb::BoltRetrieveRepair::astarHeuristic, this, _1, goal),  // the heuristic
                        // ability to disable edges (set cost to inifinity):
                        boost::weight_map(SparseEdgeWeightMap(sparseDB_->g_,
                                                              sparseDB_->edgeCollisionStatePropertySparse_, popularityBias,
                                                              popularityBiasEnabled))
                            .predecessor_map(vertexPredecessors)
                            .distance_map(&vertexDistances[0])
                            .visitor(CustomAstarVisitor(goal, this)));
  }
  catch (FoundGoalException &)
  {
    // the custom exception from CustomAstarVisitor
    OMPL_INFORM("AStar found goal vertex. distance to goal: %f", vertexDistances[goal]);
    //OMPL_INFORM("Number nodes opened: %u, Number nodes closed: %u", numNodesOpened_, numNodesClosed_);

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
    for (std::size_t i = 1; i < sparseDB_->getNumVertices(); ++i)  // skip vertex 0 b/c that is the search vertex
    {
      const SparseVertex v1 = i;
      const SparseVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4Edge(sparseDB_->getSparseStateConst(v1), sparseDB_->getSparseStateConst(v2), 10);
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

double BoltRetrieveRepair::astarHeuristic(const SparseVertex a, const SparseVertex b) const
{
  // Assume vertex 'a' is the one we care about its populariy

  // Get the classic distance
  double dist = si_->distance(sparseDB_->getSparseStateConst(a), sparseDB_->getSparseStateConst(b));

  if (false)  // method 1
  {
    const double percentMaxExtent = (sparseDB_->maxExtent_ * sparseDB_->percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
    double popularityComponent = percentMaxExtent * (sparseDB_->vertexPopularity_[a] / 100.0);

    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << sparseDB_->vertexPopularity_[a]
              << ", max extent: " << sparseDB_->maxExtent_ << ", percentMaxExtent: " << percentMaxExtent
              << ", popularityComponent: " << popularityComponent;
    dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 2
  {
    const double percentDist = (dist * sparseDB_->percentMaxExtentUnderestimate_);  // TODO(davetcoleman): cache
    double popularityComponent = percentDist * (sparseDB_->vertexPopularity_[a] / 100.0);

    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << sparseDB_->vertexPopularity_[a]
              << ", percentDist: " << percentDist << ", popularityComponent: " << popularityComponent;
    dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 3
  {
    std::cout << "astarHeuristic - dist: " << std::setprecision(4) << dist << ", popularity: " << sparseDB_->vertexPopularity_[a]
              << ", vertexPopularity_[a] / 100.0: " << sparseDB_->vertexPopularity_[a] / 100.0
              << ", percentMaxExtentUnderestimate_: " << sparseDB_->percentMaxExtentUnderestimate_;
    // if ((vertexPopularity_[a] / 100.0) < (1 - percentMaxExtentUnderestimate_))
    if (sparseDB_->vertexPopularity_[a] > (100 - sparseDB_->percentMaxExtentUnderestimate_ * 100.0))
    {
      dist = 0;
    }

    // dist = std::max(0.0, dist - popularityComponent);
  }
  else if (false)  // method 4
  {
    dist *= (1 + sparseDB_->percentMaxExtentUnderestimate_);
  }
  // method 5: increasing the sparseDelta fraction

  // std::cout << ", new distance: " << dist << std::endl;

  return dist;
}


}  // namespace bolt
}  // namespace tools
}  // namespace ompl
