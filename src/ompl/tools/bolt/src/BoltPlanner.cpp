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
#include <ompl/tools/bolt/BoltPlanner.h>
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

namespace ompl
{
namespace tools
{
namespace bolt
{
BoltPlanner::BoltPlanner(const base::SpaceInformationPtr &si, const TaskGraphPtr &taskGraph, VisualizerPtr visual)
  : base::Planner(si, "Bolt_Planner"), taskGraph_(taskGraph), visual_(visual)
{
  specs_.approximateSolutions = false;
  specs_.directed = false;

  path_simplifier_.reset(new geometric::PathSimplifier(si_));
}

BoltPlanner::~BoltPlanner(void)
{
}

void BoltPlanner::clear(void)
{
  Planner::clear();
}

void BoltPlanner::setExperienceDB(const TaskGraphPtr &taskGraph)
{
  taskGraph_ = taskGraph;
}

void BoltPlanner::setup(void)
{
  Planner::setup();
}

base::PlannerStatus BoltPlanner::solve(const base::PlannerTerminationCondition &ptc)
{
  std::size_t indent = 0;
  BOLT_CYAN_DEBUG(indent, verbose_, "BoltPlanner::solve()");
  indent += 2;

  bool solved = false;

  // Check if the database is empty
  if (taskGraph_->isEmpty())
  {
    OMPL_INFORM("Task experience database is empty so unable to run BoltPlanner algorithm.");

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

  // // Error check task planning
  // if (taskGraph_->getUseTaskTask())
  // {
  //   if (taskGraph_->getTaskLevel(startState) != 0)
  //   {
  //     OMPL_ERROR("solve: start level is %u", taskGraph_->getTaskLevel(startState));
  //     exit(-1);
  //   }
  //   if (taskGraph_->getTaskLevel(goalState) != 2)
  //   {
  //     OMPL_ERROR("solve: goal level is %u", taskGraph_->getTaskLevel(goalState));
  //     exit(-1);
  //   }
  // }

  // Create solution path as pointer so memory is not unloaded
  ob::PathPtr pathSolutionBase(new og::PathGeometric(si_));
  og::PathGeometric &geometricSolution = static_cast<og::PathGeometric &>(*pathSolutionBase);

  // Search
  if (!getPathOffGraph(startState, goalState, geometricSolution, ptc, indent))
  {
    OMPL_WARN("BoltPlanner::solve() No near start or goal found");
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

  // Save solution
  double approximateDifference = -1;
  bool approximate = false;
  pdef_->addSolutionPath(pathSolutionBase, approximate, approximateDifference, getName());
  solved = true;

  if (verbose_)
    OMPL_INFORM("  Finished BoltPlanner.solve()");
  return base::PlannerStatus(solved, approximate);
}

bool BoltPlanner::simplifyPath(og::PathGeometric &path, const base::PlannerTerminationCondition &ptc)
{
  if (verbose_)
    OMPL_INFORM("BoltPlanner: Simplifying solution (smoothing)...");

  time::point simplifyStart = time::now();
  std::size_t numStates = path.getStateCount();

  // Dave method:
  // std::size_t indent = 0;
  // taskGraph_->smoothQualityPath(&path, indent);

  path_simplifier_->simplify(path, ptc);
  double simplifyTime = time::seconds(time::now() - simplifyStart);

  OMPL_INFORM("BoltPlanner: Path simplification took %f seconds and removed %d states", simplifyTime,
              numStates - path.getStateCount());
  return true;
}

void BoltPlanner::getPlannerData(base::PlannerData &data) const
{
  if (verbose_)
    OMPL_INFORM("BoltPlanner getPlannerData");

  for (std::size_t j = 1; j < originalSolutionPath_->getStateCount(); ++j)
  {
    data.addEdge(base::PlannerDataVertex(originalSolutionPath_->getState(j - 1)),
                 base::PlannerDataVertex(originalSolutionPath_->getState(j)));
  }
}

const geometric::PathGeometric &BoltPlanner::getChosenRecallPath() const
{
  return *originalSolutionPath_;
}

std::size_t BoltPlanner::checkMotionScore(const base::State *s1, const base::State *s2) const
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

bool BoltPlanner::getPathOffGraph(const base::State *start, const base::State *goal,
                                  og::PathGeometric &geometricSolution, const base::PlannerTerminationCondition &ptc,
                                  std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, verbose_, "BoltPlanner::getPathOffGraph()");
  indent += 2;

  // Attempt to connect to graph x times, because if it fails we start adding samples
  std::size_t maxAttempts = 2;
  std::size_t attempt = 0;
  for (; attempt < maxAttempts; ++attempt)
  {
    if (verbose_)
      OMPL_INFORM("Starting getPathOffGraph attempt %u", attempt);

    // Get neighbors near start and goal. Note: potentially they are not *visible* - will test for this later

    // Start
    int level = 0;  // taskGraph_->getTaskLevel(start);
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
    level = 0;  // taskGraph_->getTaskLevel(goal);
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
                                 geometricSolution, ptc, /*debug*/ false, feedbackStartFailed, indent);

    // Error check
    if (!result)
    {
      OMPL_WARN("getPathOffGraph(): BoltPlanner returned FALSE for getPathOnGraph. Trying again in debug mode");

      // Run getPathOnGraph again in debug mode
      getPathOnGraph(startVertexCandidateNeighbors_, goalVertexCandidateNeighbors_, start, goal, geometricSolution, ptc,
                     /*debug*/ true, feedbackStartFailed, indent);

      /*
      // Add the point that failed
      if (feedbackStartFailed)  // start state failed
      {
        // Create the vertex then connect to taskGraph
        // TODO(davetcoleman): how to prevent from adding same state twice?
        TaskVertex taskV = taskGraph_->addVertex(si_->cloneState(start), QUALITY);
        taskGraph_->connectNewVertex(taskV);
      }
      else  // goal state failed
      {
        // Create the vertex then connect to taskGraph
        // TODO(davetcoleman): how to prevent from adding same state twice?
        TaskVertex taskV = taskGraph_->addVertex(si_->cloneState(goal), QUALITY);
        taskGraph_->connectNewVertex(taskV);
      }
      numStartGoalStatesAddedToTask_++;  // for later analysis

      OMPL_INFORM("Re-creating the spars graph");
      taskGraph_->createSPARS();
      */
      std::cout << "Shutting down for debugging " << std::endl;
      exit(-1);
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

bool BoltPlanner::getPathOnGraph(const std::vector<TaskVertex> &candidateStarts,
                                 const std::vector<TaskVertex> &candidateGoals, const base::State *actualStart,
                                 const base::State *actualGoal, og::PathGeometric &geometricSolution,
                                 const base::PlannerTerminationCondition &ptc, bool debug, bool &feedbackStartFailed,
                                 std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, verbose_, "BoltPlanner::getPathOnGraph()");
  indent += 2;

  bool foundValidStart = false;
  bool foundValidGoal = false;

  // Try every combination of nearby start and goal pairs
  BOOST_FOREACH (TaskVertex start, candidateStarts)
  {
    if (actualStart == taskGraph_->getVertexState(start))
    {
      OMPL_ERROR("Comparing same start state");
      exit(-1);  // just curious if this ever happens, no need to actually exit
      continue;
    }

    // Check if this start is visible from the actual start
    if (!si_->checkMotion(actualStart, taskGraph_->getVertexState(start)))
    {
      if (verbose_)
      {
        OMPL_WARN("FOUND START CANDIDATE THAT IS NOT VISIBLE ");
      }
      if (debug)
      {
        visual_->viz4()->state(taskGraph_->getVertexState(start), tools::LARGE, tools::RED, 1);
        visual_->viz4()->edge(actualStart, taskGraph_->getVertexState(start), 100);
        visual_->viz4()->trigger();
        usleep(0.1 * 1000000);
      }
      continue;  // this is actually not visible
    }
    foundValidStart = true;

    BOOST_FOREACH (TaskVertex goal, candidateGoals)
    {
      if (actualGoal == taskGraph_->getVertexState(goal))
      {
        OMPL_ERROR("Comparing same goal state");
        continue;
      }

      if (verbose_)
        OMPL_INFORM("    foreach_goal: Checking motion from  %d to %d", actualGoal, taskGraph_->getVertexState(goal));

      // Check if our planner is out of time
      if (ptc == true)
      {
        OMPL_DEBUG("getPathOnGraph function interrupted because termination condition is true.");
        return false;
      }

      // Check if this goal is visible from the actual goal
      if (!si_->checkMotion(actualGoal, taskGraph_->getVertexState(goal)))
      {
        if (verbose_)
        {
          OMPL_WARN("FOUND GOAL CANDIDATE THAT IS NOT VISIBLE! ");
        }

        if (debug)
        {
          visual_->viz4()->state(taskGraph_->getVertexState(goal), tools::SMALL, tools::RED, 1);
          visual_->viz4()->edge(actualGoal, taskGraph_->getVertexState(goal), 100);
          visual_->viz4()->trigger();
          usleep(0.1 * 1000000);
        }

        continue;  // this is actually not visible
      }
      foundValidGoal = true;

      // Repeatidly search through graph for connection then check for collisions then repeat
      if (lazyCollisionSearch(start, goal, actualStart, actualGoal, geometricSolution, ptc, indent))
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

bool BoltPlanner::lazyCollisionSearch(const TaskVertex &start, const TaskVertex &goal, const base::State *actualStart,
                                      const base::State *actualGoal, og::PathGeometric &geometricSolution,
                                      const base::PlannerTerminationCondition &ptc, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, verbose_, "BoltPlanner::lazyCollisionSearch()");
  indent += 2;

  // Vector to store candidate paths in before they are converted to PathPtrs
  std::vector<TaskVertex> vertexPath;
  double distance;  // resulting path distance

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
  assert(taskGraph_->getVertexState(start));
  assert(taskGraph_->getVertexState(goal));

  // Check that our states are on the same connected component
  // TODO: in the future the graph should always just be fully connected
  // so perhaps this check would not be necessary
  if (!taskGraph_->sameComponent(start, goal))
  {
    OMPL_WARN("Found start and goal states that are on different connected components!");
    return false;
  }

  // Visualize start vertex
  const bool visualize = false;
  if (visualize)
  {
    OMPL_INFORM("viz start -----------------------------");
    visual_->viz5()->state(taskGraph_->getVertexState(start), tools::VARIABLE_SIZE, tools::PURPLE, 1);
    visual_->viz5()->edge(actualStart, taskGraph_->getVertexState(start), 30);
    visual_->viz5()->trigger();
    usleep(5 * 1000000);

    // Visualize goal vertex
    OMPL_INFORM("viz goal ------------------------------");
    visual_->viz5()->state(taskGraph_->getVertexState(goal), tools::VARIABLE_SIZE, tools::PURPLE, 1);
    visual_->viz5()->edge(actualGoal, taskGraph_->getVertexState(goal), 0);
    visual_->viz5()->trigger();
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
    if (!taskGraph_->astarSearch(start, goal, vertexPath, distance, indent))
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
    if (lazyCollisionCheck(vertexPath, ptc, indent))
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

bool BoltPlanner::lazyCollisionCheck(std::vector<TaskVertex> &vertexPath, const base::PlannerTerminationCondition &ptc, std::size_t indent)
{
  BOLT_CYAN_DEBUG(indent, verbose_, "BoltPlanner::lazyCollisionCheck()");
  indent += 2;

  bool hasInvalidEdges = false;

  // Initialize
  TaskVertex fromVertex = vertexPath[0];
  TaskVertex toVertex;

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

    TaskEdge thisEdge = boost::edge(fromVertex, toVertex, taskGraph_->g_).first;

    // Has this edge already been checked before?
    if (taskGraph_->edgeCollisionStatePropertyTask_[thisEdge] == NOT_CHECKED)
    {
      // Check path between states
      if (!si_->checkMotion(taskGraph_->getVertexState(fromVertex), taskGraph_->getVertexState(toVertex)))
      {
        // Path between (from, to) states not valid, disable the edge
        OMPL_INFORM("  DISABLING EDGE from vertex %f to vertex %f", fromVertex, toVertex);

        // Disable edge
        taskGraph_->edgeCollisionStatePropertyTask_[thisEdge] = IN_COLLISION;
      }
      else
      {
        // Mark edge as free so we no longer need to check for collision
        taskGraph_->edgeCollisionStatePropertyTask_[thisEdge] = FREE;
      }
    }

    // Check final result
    if (taskGraph_->edgeCollisionStatePropertyTask_[thisEdge] == IN_COLLISION)
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

bool BoltPlanner::findGraphNeighbors(const base::State *state, std::vector<TaskVertex> &graphNeighborhood,
                                     int requiredLevel)
{
  // Benchmark runtime
  time::point startTime = time::now();

  // Reset
  graphNeighborhood.clear();

  // Setup search by getting a non-const version of the focused state
  const std::size_t threadID = 0;
  base::State *stateCopy = si_->cloneState(state);
  taskGraph_->getQueryStateNonConst(taskGraph_->queryVertices_[threadID]) = stateCopy;

  // Search
  double find_nearest_k_neighbors;
  if (si_->getStateSpace()->getDimension() == 3)
    find_nearest_k_neighbors = 10;
  else
    find_nearest_k_neighbors = 30;
  taskGraph_->nn_->nearestK(taskGraph_->queryVertices_[threadID], find_nearest_k_neighbors, graphNeighborhood);

  // Reset
  taskGraph_->getQueryStateNonConst(taskGraph_->queryVertices_[threadID]) = nullptr;

  // Filter neighbors based on level
  if (requiredLevel >= 0)
  {
    // removeVerticesNotOnLevel(graphNeighborhood, requiredLevel);
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

// bool BoltPlanner::removeVerticesNotOnLevel(std::vector<TaskVertex> &graphNeighborhood, int level)
// {
//   std::size_t original_size = graphNeighborhood.size();

//   // Remove edges based on layer
//   for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
//   {
//     const TaskVertex &nearVertex = graphNeighborhood[i];

//     // Make sure state is on correct level
//     if (taskGraph_->getTaskLevel(nearVertex) != static_cast<std::size_t>(level))
//     {
//       if (verbose_)
//         std::cout << "      Skipping neighbor " << nearVertex << ", i=" << i
//                   << ", because wrong level: " << taskGraph_->getTaskLevel(nearVertex) << ", desired level: " <<
//                   level
//                   << std::endl;
//       graphNeighborhood.erase(graphNeighborhood.begin() + i);
//       i--;
//       continue;
//     }
//   }

//   if (verbose_)
//     OMPL_INFORM("    removeVerticesNotOnLevel(): states require level %d, removed: %u, remaining: %u", level,
//                 original_size - graphNeighborhood.size(), graphNeighborhood.size());

//   return true;
// }

bool BoltPlanner::convertVertexPathToStatePath(std::vector<TaskVertex> &vertexPath, const base::State *actualStart,
                                               const base::State *actualGoal, og::PathGeometric &geometricSolution)
{
  // Ensure the input path is not empty
  if (!vertexPath.size())
    return false;

  // og::PathGeometric *pathGeometric = new og::PathGeometric(si_);

  // Add original start if it is different than the first state
  if (actualStart != taskGraph_->getVertexState(vertexPath.back()))
  {
    geometricSolution.append(actualStart);
  }

  // Error check that no consequtive verticies are the same
  /*std::cout << "BoltPlanner.Vertices: ";
  for (std::size_t i = vertexPath.size(); i > 0; --i)
  {
      std::cout << vertexPath[i - 1] << ", ";
  }
  std::cout << std::endl;*/

  // Reverse the vertexPath and convert to state path
  for (std::size_t i = vertexPath.size(); i > 0; --i)
  {
    geometricSolution.append(taskGraph_->getVertexState(vertexPath[i - 1]));

    // Add the edge status
    if (i > 1)  // skip the last vertex (its reversed)
    {
      // Error check that no consequtive verticies are the same
      if (vertexPath[i - 1] == vertexPath[i - 2])
      {
        OMPL_ERROR("Found repeated vertices %u to %u on index %u", vertexPath[i - 1], vertexPath[i - 2], i);
        exit(-1);
      }

      TaskEdge edge = boost::edge(vertexPath[i - 1], vertexPath[i - 2], taskGraph_->g_).first;

      /* This functionality has moved to TaskGraph
      if (taskGraph_->getPopularityBiasEnabled())
      {
          // reduce cost of this edge because it was just used
          static const double REDUCTION_AMOUNT = 10;
          std::cout << "Previous edge weight for " << i << " of edge " << edge << std::endl;
          std::cout << "    is " << taskGraph_->edgeWeightProperty_[edge];
          taskGraph_->edgeWeightProperty_[edge] =
              std::max(taskGraph_->edgeWeightProperty_[edge] - REDUCTION_AMOUNT, 0.0);
          std::cout << " new: " << taskGraph_->edgeWeightProperty_[edge] << std::endl;
      }
      */

      // Check if any edges in path are not free (then it an approximate path)
      if (taskGraph_->edgeCollisionStatePropertyTask_[edge] == IN_COLLISION)
      {
        OMPL_ERROR("Found invalid edge / approximate solution - how did this happen?");
      }
      else if (taskGraph_->edgeCollisionStatePropertyTask_[edge] == NOT_CHECKED)
      {
        OMPL_ERROR("A chosen path has an edge that has not been checked for collision. This should not happen");
      }
    }
  }

  // Add original goal if it is different than the last state
  if (actualGoal != taskGraph_->getVertexState(vertexPath.front()))
  {
    geometricSolution.append(actualGoal);
  }

  // if (taskGraph_->getPopularityBiasEnabled())
  // {
  //     // Ensure graph doesn't get too popular
  //     taskGraph_->normalizeGraphEdgeWeights();

  //     // Mark as needing to be saved to file
  //     taskGraph_->graphUnsaved_ = true;
  // }

  return true;
}

bool BoltPlanner::canConnect(const base::State *randomState, const base::PlannerTerminationCondition &ptc)
{
  std::vector<TaskVertex> candidateNeighbors;

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
  BOOST_FOREACH (TaskVertex nearState, candidateNeighbors)
  {
    const base::State *s1 = randomState;
    const base::State *s2 = taskGraph_->getVertexState(nearState);

    // Check if this nearState is visible from the random state
    if (!si_->checkMotion(s1, s2))
    {
      OMPL_WARN("NEIGHBOR %u NOT VISIBLE ", count++);

      if (false)
      {
        visual_->viz5()->state(s2, tools::MEDIUM, tools::BLUE, 1);
        visual_->viz5()->edge(s1, s2, 100);
        visual_->viz5()->trigger();
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
            visual_->viz5()->state(interState, tools::LARGE, tools::RED, 1);
            visual_->viz5()->trigger();
            usleep(1 * 1000000);
          }
          else
          {
            // visual_->viz5()->state(interState, /*mode=*/1, 1); // GREEN
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

/*bool BoltPlanner::astarSearch(const TaskVertex start, const TaskVertex goal, std::vector<TaskVertex>
&vertexPath)
{
  // Hold a list of the shortest path parent to each vertex
  TaskVertex *vertexPredecessors = new TaskVertex[taskGraph_->getNumVertices()];
  // boost::vector_property_map<TaskVertex> vertexPredecessors(getNumVertices());

  bool foundGoal = false;
  double *vertexDistances = new double[taskGraph_->getNumVertices()];

  // Reset statistics
  numNodesOpened_ = 0;
  numNodesClosed_ = 0;

  OMPL_INFORM("Beginning AStar Search");
  try
  {
    double popularityBias = 0;
    bool popularityBiasEnabled = false;
    // Note: could not get astar_search to compile within BoltPlanner.cpp class because of
    // namespacing issues
    boost::astar_search(taskGraph_->g_,                                                 // graph
                        start,                                                        // start state
                        boost::bind(&otb::BoltPlanner::astarHeuristic, this, _1, goal),  // the heuristic
                        // ability to disable edges (set cost to inifinity):
                        boost::weight_map(TaskEdgeWeightMap(taskGraph_->g_,
                                                              taskGraph_->edgeCollisionStatePropertyTask_,
popularityBias,
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
      TaskVertex v;
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
    for (std::size_t i = 1; i < taskGraph_->getNumVertices(); ++i)  // skip vertex 0 b/c that is the search vertex
    {
      const TaskVertex v1 = i;
      const TaskVertex v2 = vertexPredecessors[v1];
      if (v1 != v2)
      {
        // std::cout << "Edge " << v1 << " to " << v2 << std::endl;
        visual_->viz4()->edge(taskGraph_->getVertexState(v1), taskGraph_->getVertexState(v2), 10);
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
*/

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
