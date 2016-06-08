/** \brief Error checking function to ensure solution has correct task path/level changes */
bool checkTaskPathSolution(geometric::PathGeometric &path, base::State *start, base::State *goal);



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
