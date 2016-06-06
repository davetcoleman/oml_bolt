/** \brief Compute the heuristic distance between the current node and the next goal */
double distanceFunctionTasks(const TaskVertex a, const TaskVertex b) const;

/** \brief Error checking function to ensure solution has correct task path/level changes */
bool checkTaskPathSolution(geometric::PathGeometric &path, base::State *start, base::State *goal);


/** \brief Testing code for integrating Decartes */
bool addCartPath(std::vector<base::State *> path);

double TaskGraph::distanceFunctionTasks(const TaskVertex a, const TaskVertex b) const
{
  // Do not use task distance if that mode is not enabled
  if (!useTaskPlanning_)
    return distanceFunction(a, b);

  const bool verbose = false;
  // Reorder a & b so that we are sure that a.level <= b.level
  std::size_t taskLevelA = getTaskLevel(a);
  std::size_t taskLevelB = getTaskLevel(b);
  if (taskLevelA > taskLevelB)
  {
    // Call itself again, this time switching ordering
    if (verbose)
      OMPL_INFORM("Switched ordering for distanceFunctionTasks()");
    return distanceFunctionTasks(b, a);
  }

  double dist = 0;                            // the result
  static const double TASK_LEVEL_COST = 1.0;  // cost to change levels/tasks

  // Error check
  assert(stateProperty_[a]);
  assert(stateProperty_[b]);
  assert(stateProperty_[startConnectorVertex_]);
  assert(stateProperty_[endConnectorVertex_]);

  if (taskLevelA == 0)
  {
    if (taskLevelB == 0)  // regular distance for bottom level
    {
      if (verbose)
        std::cout << "a ";
      dist = si_->distance(stateProperty_[a], stateProperty_[b]);
    }
    else if (taskLevelB == 1)
    {
      if (verbose)
        std::cout << "b ";
      dist = si_->distance(stateProperty_[a], stateProperty_[startConnectorVertex_]) + TASK_LEVEL_COST +
             si_->distance(stateProperty_[startConnectorVertex_], stateProperty_[b]);
    }
    else if (taskLevelB == 2)
    {
      if (verbose)
        std::cout << "c ";
      dist = si_->distance(stateProperty_[a], stateProperty_[startConnectorVertex_]) + TASK_LEVEL_COST +
             distanceAcrossCartesian_ + TASK_LEVEL_COST +
             si_->distance(stateProperty_[endConnectorVertex_], stateProperty_[b]);
    }
    else
    {
      OMPL_ERROR("Unknown task level mode");
      exit(-1);
    }
  }
  else if (taskLevelA == 1)
  {
    if (taskLevelB == 0)
    {
      OMPL_ERROR("Unknown task level mode");
      exit(-1);
    }
    else if (taskLevelB == 1)
    {
      if (verbose)
        std::cout << "d ";
      dist = si_->distance(stateProperty_[a], stateProperty_[b]);
    }
    else if (taskLevelB == 2)
    {
      if (verbose)
        std::cout << "e ";
      dist = si_->distance(stateProperty_[a], stateProperty_[endConnectorVertex_]) + TASK_LEVEL_COST +
             si_->distance(stateProperty_[endConnectorVertex_], stateProperty_[b]);
    }
    else
    {
      OMPL_WARN("Unknown task level mode");
    }
  }
  else if (taskLevelA == 2)
  {
    if (taskLevelB == 0 || taskLevelB == 1)
    {
      OMPL_ERROR("Unknown task level mode");
      exit(-1);
    }
    else if (taskLevelB == 2)
    {
      if (verbose)
        std::cout << "f ";
      dist = si_->distance(stateProperty_[a], stateProperty_[b]);
    }
    else
    {
      OMPL_WARN("Unknown task level mode");
    }
  }
  else
  {
    OMPL_WARN("Unknown task level mode");
  }

  if (verbose)
    std::cout << "State " << a << " @level " << taskLevelA << " to State " << b << " @level " << taskLevelB
              << " has distance " << dist << std::endl;
  return dist;
}

bool TaskGraph::checkTaskPathSolution(og::PathGeometric &path, ob::State *start, ob::State *goal)
{
  bool error = false;
  int current_level = 0;

  for (std::size_t i = 0; i < path.getStateCount(); ++i)
  {
    int level = si_->getStateSpace()->getLevel(path.getState(i));

    // Check if start state is correct
    if (i == 0)
    {
      if (!si_->getStateSpace()->equalStates(path.getState(i), start))
      {
        OMPL_ERROR("Start state of path is not same as original problem");
        error = true;
      }

      if (level != 0)
      {
        OMPL_ERROR("Start state is not at level 0, instead %i", level);
        error = true;
      }
    }

    // Check if goal state is correct
    if (i == path.getStateCount() - 1)
    {
      if (!si_->getStateSpace()->equalStates(path.getState(i), goal))
      {
        OMPL_ERROR("Goal state of path is not same as original problem");
        error = true;
      }

      if (level != 2)
      {
        OMPL_ERROR("Goal state is not at level 2, instead %i", level);
        error = true;
      }
    }

    // Ensure that level is always increasing
    if (level < current_level)
    {
      OMPL_ERROR("State decreased in level (%i) from previous level of ", current_level);
      error = true;
    }
    current_level = level;

  }  // for loop

  // Show more data if error
  if (error)
  {
    OMPL_ERROR("Showing data on path:");
    for (std::size_t i = 0; i < path.getStateCount(); ++i)
    {
      int level = si_->getStateSpace()->getLevel(path.getState(i));
      OMPL_INFORM(" - Path state %i has level %i", i, level);
    }
  }

  return error;
}
