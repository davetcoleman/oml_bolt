/** \brief Compute the heuristic distance between the current node and the next goal */
double distanceFunctionTasks(const TaskVertex a, const TaskVertex b) const;

/** \brief Clone the graph to have second and third layers for task planning then free space planning */
void generateTaskSpace();

/** \brief Get k number of neighbors near a state at a certain level that have valid motions */
void getNeighborsAtLevel(const base::State *origState, const std::size_t level, const std::size_t kNeighbors,
                         std::vector<TaskVertex> &neighbors);

/** \brief Error checking function to ensure solution has correct task path/level changes */
bool checkTaskPathSolution(geometric::PathGeometric &path, base::State *start, base::State *goal);

/**
 * \brief Helper for connecting both sides of a cartesian path into a dual level graph
 * \param fromState - the endpoint (start or goal) we are connecting from the cartesian path to the graph
 * \param level - what task level we are connecting to - either 0 or 2 (bottom layer or top layer)
 * \param minConnectorVertex - the vertex on the main graph that has the shortest path to connecting to the
 * cartesian path
 * \return true on success
 */
bool connectStateToNeighborsAtLevel(const TaskVertex &fromVertex, const std::size_t level,
                                    TaskVertex &minConnectorVertex);

/** \brief Testing code for integrating Decartes */
bool addCartPath(std::vector<base::State *> path);

double DenseDB::distanceFunctionTasks(const TaskVertex a, const TaskVertex b) const
{
  // Do not use task distance if that mode is not enabled
  if (!useTaskTask_)
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

void DenseDB::generateTaskSpace()
{
  OMPL_INFORM("Generating task space");
  std::vector<TaskVertex> vertexToNewVertex(getNumVertices());

  OMPL_INFORM("Adding task space vertices");
  foreach (TaskVertex v, boost::vertices(g_))
  {
    // The first vertex (id=0) should have a nullptr state because it is used for searching
    if (!stateProperty_[v])
    {
      continue;
    }

    // Clone the state
    base::State *newState = si_->cloneState(stateProperty_[v]);

    const std::size_t level = 2;
    si_->getStateSpace()->setLevel(newState, level);

    // Add the state back
    TaskVertex vNew = addVertex(newState, START);

    // Map old vertex to new vertex
    vertexToNewVertex[v] = vNew;

    // if (visualizeGridGeneration_)  // Visualize - only do this for 2/3D environments
    {
      visual_->viz1State(newState, /*mode=*/5, 1);  // Candidate node has already (just) been added
      visual_->viz1Trigger();
      usleep(0.001 * 1000000);
    }
  }

  // Add Edges
  OMPL_INFORM("Adding task space edges");

  // Cache all current edges before adding new ones
  std::vector<TaskVertex> edges(getNumEdges());
  std::size_t edgeID = 0;
  foreach (const TaskVertex e, boost::edges(g_))
  {
    edges[edgeID++] = e;
  }

  // Copy every edge
  for (std::size_t i = 0; i < edges.size(); ++i)
  {
    const TaskVertex &e = edges[i];

    const TaskVertex v1 = boost::source(e, g_);
    const TaskVertex v2 = boost::target(e, g_);

    // std::cout << "v1: " << v1 << " v2: " << v2 << std::endl;
    // std::cout << "v1': " << vertexToNewVertex[v1] << " v2': " << vertexToNewVertex[v2] << std::endl;

    TaskVertex newE = addEdge(vertexToNewVertex[v1], vertexToNewVertex[v2], edgeWeightProperty_[e]);

    // if (visualizeGridGeneration_)  // Visualize
    {
      viz1Edge(newE);
      if (i % 100 == 0)
      {
        visual_->viz1Trigger();
        usleep(0.001 * 1000000);
      }
    }
  }

  // if (visualizeGridGeneration_)  // Visualize
  visual_->viz1Trigger();

  OMPL_INFORM("Done generating task space graph");
}

void DenseDB::getNeighborsAtLevel(const base::State *origState, const std::size_t level, const std::size_t kNeighbors,
                                  std::vector<TaskVertex> &neighbors)
{
  // Clone the state and change its level
  base::State *searchState = si_->cloneState(origState);
  si_->getStateSpace()->setLevel(searchState, level);

  // Get nearby state
  stateProperty_[queryVertex_] = searchState;
  nn_->nearestK(queryVertex_, kNeighbors, neighbors);
  stateProperty_[queryVertex_] = nullptr;  // Set search vertex to nullptr to prevent segfault on class unload of memory
  si_->getStateSpace()->freeState(searchState);

  // Run various checks
  for (std::size_t i = 0; i < neighbors.size(); ++i)
  {
    const TaskVertex &nearVertex = neighbors[i];

    // Make sure state is on correct level
    if (getTaskLevel(nearVertex) != level)
    {
      std::cout << "      Skipping neighbor " << nearVertex << ", i=" << i
                << ", because wrong level: " << getTaskLevel(nearVertex) << ", desired level: " << level << std::endl;
      neighbors.erase(neighbors.begin() + i);
      i--;
      continue;
    }

    // Collision check
    if (!si_->checkMotion(origState, stateProperty_[nearVertex]))  // is valid motion
    {
      std::cout << "      Skipping neighbor " << nearVertex << ", i=" << i << ", at level=" << getTaskLevel(nearVertex)
                << " because invalid motion" << std::endl;
      neighbors.erase(neighbors.begin() + i);
      i--;
      continue;
    }

    std::cout << "      Keeping neighbor " << nearVertex << std::endl;
  }
}

bool DenseDB::checkTaskPathSolution(og::PathGeometric &path, ob::State *start, ob::State *goal)
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

bool DenseDB::connectStateToNeighborsAtLevel(const TaskVertex &fromVertex, const std::size_t level,
                                             TaskVertex &minConnectorVertex)
{
  // Get nearby states to goal
  std::vector<TaskVertex> neighbors;
  const std::size_t kNeighbors = 20;
  getNeighborsAtLevel(stateProperty_[fromVertex], level, kNeighbors, neighbors);

  // Error check
  if (neighbors.empty())
  {
    OMPL_ERROR("No neighbors found when connecting cartesian path");
    return false;
  }
  else if (neighbors.size() < 3)
  {
    OMPL_WARN("Only %u neighbors found on level %u", neighbors.size(), level);
  }
  else
    OMPL_INFORM("Found %u neighbors on level %u", neighbors.size(), level);

  // Find the shortest connector out of all the options
  double minConnectorCost = std::numeric_limits<double>::infinity();

  // Loop through each neighbor
  foreach (TaskVertex v, neighbors)
  {
    // Add edge from nearby graph vertex to cart path goal
    double connectorCost = distanceFunction(fromVertex, v);
    addEdge(fromVertex, v, connectorCost);

    // Get min cost connector
    if (connectorCost < minConnectorCost)
    {
      minConnectorCost = connectorCost;  // TODO(davetcoleman): should we save the cost, or just use 1.0?
      minConnectorVertex = v;
    }

    // Visualize connection to goal of cartesian path
    const double cost = (level == 0 ? 100 : 50);
    if (visualizeCartNeighbors_)
    {
      visual_->viz1Edge(stateProperty_[v], stateProperty_[fromVertex], cost);
      visual_->viz1State(stateProperty_[v], /*mode=*/1, 1);
      visual_->viz1Trigger();
      usleep(0.001 * 1000000);
    }
  }

  // Display ---------------------------------------
  if (visualizeCartNeighbors_)
    visual_->viz1Trigger();

  return true;
}

bool DenseDB::addCartPath(std::vector<base::State *> path)
{
  // Error check
  if (path.size() < 2)
  {
    OMPL_ERROR("Invalid cartesian path - too few states");
    exit(-1);
  }
  // TODO: check for validity

  // Create verticies for the extremas
  TaskVertex startVertex = addVertex(path.front(), CARTESIAN);
  TaskVertex goalVertex = addVertex(path.back(), CARTESIAN);

  // Record min cost for cost-to-go heurstic distance function later
  distanceAcrossCartesian_ = distanceFunction(startVertex, goalVertex);

  // Connect Start to graph --------------------------------------
  std::cout << "  start connector ---------" << std::endl;
  const std::size_t level0 = 0;
  if (!connectStateToNeighborsAtLevel(startVertex, level0, startConnectorVertex_))
  {
    OMPL_ERROR("Failed to connect start of cartesian path");
    return false;
  }

  // Connect goal to graph --------------------------------------
  std::cout << "  goal connector ----------------" << std::endl;
  const std::size_t level2 = 2;
  if (!connectStateToNeighborsAtLevel(goalVertex, level2, endConnectorVertex_))
  {
    OMPL_ERROR("Failed to connect goal of cartesian path");
    return false;
  }

  // Add cartesian path to mid level graph --------------------
  TaskVertex v1 = startVertex;
  TaskVertex v2;
  for (std::size_t i = 1; i < path.size(); ++i)
  {
    // Check if we are on the goal vertex
    if (i == path.size() - 1)
    {
      v2 = goalVertex;  // Do not create the goal vertex twice
    }
    else
    {
      v2 = addVertex(path[i], CARTESIAN);
    }
    double cost = distanceFunction(v1, v2);
    addEdge(v1, v2, cost);
    v1 = v2;

    // Visualize
    if (visualizeCartPath_)
    {
      visual_->viz1Edge(path[i - 1], path[i], 0);
      visual_->viz1State(path[i], /*mode=*/1, 1);
      visual_->viz1Trigger();
      usleep(0.001 * 1000000);
    }
  }

  return true;
}
