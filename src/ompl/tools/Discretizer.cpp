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
   Desc:   Utility to discretize a space into a uniform grid
*/

// OMPL
#include <ompl/tools/bolt/Discretizer.h>
#include <ompl/tools/bolt/DenseDB.h>

// Boost
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/math/constants/constants.hpp>

// Profiling
#include <valgrind/callgrind.h>

#define foreach BOOST_FOREACH

namespace og = ompl::geometric;
namespace ot = ompl::tools;
namespace ob = ompl::base;

namespace ompl
{
namespace tools
{
namespace bolt
{
Discretizer::Discretizer(base::SpaceInformationPtr si, DenseDB *denseDB, DenseCachePtr denseCache,
                         VisualizerPtr visual)
  : si_(si), denseDB_(denseDB), denseCache_(denseCache), visual_(visual)
{
  numThreads_ = boost::thread::hardware_concurrency();

  // Debugging
  if (false)
  {
    OMPL_WARN("Overriding number of threads for testing to 1");
    numThreads_ = 1;
  }
}

Discretizer::~Discretizer(void)
{
}

bool Discretizer::generateGrid()
{
  OMPL_INFORM("Generating grid");

  // Benchmark runtime
  time::point totalStartTime = time::now();

  if (!si_->isSetup())
  {
    OMPL_WARN("Space information setup was not yet called. Calling now.");
    si_->setup();
  }

  double vertexDuration;
  {
    // Benchmark runtime
    time::point startTime = time::now();

    // Create vertices
    generateVertices();

    // Benchmark runtime
    vertexDuration = time::seconds(time::now() - startTime);
  }

  OMPL_INFORM("Generated %i vertices in %f sec (%f hours)", denseDB_->getNumVertices(), vertexDuration,
              vertexDuration / 60.0 / 60.0);

  // Error check
  if (denseDB_->getNumVertices() < 2)
  {
    OMPL_ERROR("No vertices generated, failing");
    exit(-1);
  }

  double edgeDuration;
  {
    // Benchmark runtime
    time::point startTime = time::now();

    // Create edges ----------------------------------------
    generateEdges();

    // Benchmark runtime
    edgeDuration = time::seconds(time::now() - startTime);
  }

  // Get the average vertex degree (number of connected edges)
  double averageDegree = (denseDB_->getNumEdges() * 2) / static_cast<double>(denseDB_->getNumVertices());

  // Check how many disjoint sets are in the dense graph (should be none)
  std::size_t numSets = denseDB_->checkConnectedComponents();

  // Total benchmark runtime
  double totalDuration = time::seconds(time::now() - totalStartTime);

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Discretization stats:");
  OMPL_INFORM("   Total valid vertices:   %u", denseDB_->getNumVertices());
  OMPL_INFORM("   Total valid edges:      %u", denseDB_->getNumEdges());
  OMPL_INFORM("   Num edges colliding:    %f", numEdgesInCollision_);
  OMPL_INFORM("   Average degree:         %f", averageDegree);
  OMPL_INFORM("   Connected Components:   %u", numSets);
  OMPL_INFORM("   Edge connection method: %u", edgeConnectionStrategy_);
  OMPL_INFORM("   Vertex generation time: %f seconds (%f min)", vertexDuration, vertexDuration / 60.0);
  OMPL_INFORM("   Edge generation time:   %f seconds (%f min)", edgeDuration, edgeDuration / 60.0);
  OMPL_INFORM("   Total grid gen. time:   %f seconds (%f min)", totalDuration, totalDuration / 60.0);
  OMPL_INFORM("------------------------------------------------------");

  // Display
  // if (visualizeGridGeneration_)
  // visual_->viz1Trigger();

  denseDB_->graphUnsaved_ = true;
  return true;
}

void Discretizer::generateVertices()
{
  const bool verbose = false;

  // Setup bounds for joint 0
  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  assert(bounds.high.size() == bounds.low.size());
  assert(bounds.high.size() == si_->getStateSpace()->getDimension());

  // Divide joint 0 between threads
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID];
  const std::size_t jointIncrements = range / discretization_;

  // Check that we have enough jointIncrements for all the threads
  if (jointIncrements < numThreads_)
  {
    OMPL_WARN("There are fewer joint_0 increments (%u) at current discretization (%f) than available threads (%u), "
              "underutilizing threading",
              jointIncrements, discretization_, numThreads_);
    OMPL_INFORM("Optimal discretization: %f", range / double(numThreads_));
    numThreads_ = jointIncrements;
  }

  // Warn
  // if (numThreads_ > 1 && visualizeGridGeneration_)
  // OMPL_WARN("Multiple threads are trying to visualize, could cause race conditions");

  // Setup threading
  std::vector<boost::thread *> threads(numThreads_);
  OMPL_INFORM("Generating vertices using %u threads", numThreads_);

  std::size_t jointIncrementsPerThread = jointIncrements / numThreads_;

  if (true)
  {
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Discretization Setup: " << std::endl;
    std::cout << "  Dimensions:             " << si_->getStateSpace()->getDimension() << std::endl;
    std::cout << "  Discretization:         " << discretization_ << std::endl;
    std::cout << "  J0 Low:                 " << bounds.low[jointID] << std::endl;
    std::cout << "  J0 High:                " << bounds.high[jointID] << std::endl;
    std::cout << "  J0 range is:            " << range << std::endl;
    std::cout << "  J0 Increments:          " << jointIncrements << std::endl;
    std::cout << "  J0 IncrementsPerThread: " << jointIncrementsPerThread << std::endl;
    std::cout << "  Total states:           " << pow(jointIncrements, si_->getStateSpace()->getDimension())
              << std::endl;
    std::cout << "  NumThreads: " << numThreads_ << std::endl;
    std::cout << "-------------------------------------------------------" << std::endl;
  }

  double startJointValue = bounds.low[jointID];
  double endJointValue;

  // For each thread
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    endJointValue = startJointValue + jointIncrementsPerThread * discretization_;

    // Check if this is the last thread
    if (i == threads.size() - 1)
    {
      // have it do remaining bounds
      endJointValue = bounds.high[jointID];
    }

    if (verbose)
      std::cout << "Thread " << i << " has values from " << startJointValue << " to " << endJointValue << std::endl;

    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    threads[i] =
        new boost::thread(boost::bind(&Discretizer::generateVerticesThread, this, i, startJointValue, endJointValue, si));

    startJointValue = endJointValue;
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }
}

void Discretizer::generateVerticesThread(std::size_t threadID, double startJointValue, double endJointValue,
                                     base::SpaceInformationPtr si)
{
  std::size_t jointID = 0;
  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();
  base::State *candidateState = si->getStateSpace()->allocState();

  // Prepare for recursion
  std::vector<double> values(si->getStateSpace()->getDimension(), 0);

  // Customize for different state spaces TODO more generic
  std::size_t maxDiscretizationLevel;
  const std::size_t dim = si->getStateSpace()->getDimension();
  if (dim == 3)
  {
    maxDiscretizationLevel = 1;  // because the third level (numbered '2') is for task space
    values[2] = 0.0;             // task space
  }
  else if (dim == 6 && false)
  {
    maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
    values[5] = 0.0;                   // middle rotation of wrist
  }
  else if (dim == 12)
  {
    maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
    values[5] = 0.0;                   // middle rotation of wrist
    values[11] = 0.0;                  // middle rotation of wrist
  }
  else
  {
    maxDiscretizationLevel = dim - 1;
  }

  // Loop through current joint
  for (double value = startJointValue; value < endJointValue; value += discretization_)
  {
    // User feedback on thread 0
    if (threadID == numThreads_ - 1)
    {
      const double percent = (value - startJointValue) / (endJointValue - startJointValue) * 100.0;
      std::cout << "Vertex generation progress: " << std::setprecision(1) << percent
                << " % Total vertices: " << denseDB_->getNumVertices() - 1 /*ignore query vertex*/ << std::endl;
    }

    // Set first joint value
    values[jointID] = value;

    // Keep recursing
    recursiveDiscretization(threadID, values, jointID + 1, si, candidateState, maxDiscretizationLevel);
  }

  // Cleanup
  si->freeState(candidateState);
}

void Discretizer::recursiveDiscretization(std::size_t threadID, std::vector<double> &values, std::size_t jointID,
                                          base::SpaceInformationPtr si, base::State *candidateState,
                                          std::size_t maxDiscretizationLevel)
{
  const bool verbose = false;

  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();

  // Error check
  assert(jointID < values.size());

  // Loop through current joint
  for (double value = bounds.low[jointID]; value <= bounds.high[jointID]; value += discretization_)
  {
    // User feedback on thread 0 for high dimension spaces
    if (jointID == 1 && threadID == 0 && si->getStateSpace()->getDimension() > 3)
    {
      const double percent = (value - bounds.low[jointID]) / (bounds.high[jointID] - bounds.low[jointID]) * 100.0;
      std::cout << "Level 1 vertex generation progress: " << std::setprecision(1) << percent
                << " % Total vertices: " << denseDB_->getNumVertices() - 1 /*ignore query vertex*/ << std::endl;
    }

    // Set value
    values[jointID] = value;

    // Check if we are at the end of the recursion
    if (verbose)
    {
      std::cout << "threadID " << threadID << " jointID " << jointID << " of " << maxDiscretizationLevel << " --- ";
      std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout << std::endl;
    }
    //<< " low: " << bounds.low[jointID] << " high " << bounds.high[jointID] << std::endl;

    if (jointID < maxDiscretizationLevel)
    {
      // Keep recursing

      // Special rule for 12dof
      if (si->getStateSpace()->getDimension() == 12 && jointID == 4)
      {
        // skip joint id 5 (joint 6)
        recursiveDiscretization(threadID, values, jointID + 2, si, candidateState, maxDiscretizationLevel);
      }
      else  // regular treatment
      {
        recursiveDiscretization(threadID, values, jointID + 1, si, candidateState, maxDiscretizationLevel);
      }
    }
    else  // this is the end of recursion, create a new state
    {
      // Fill the state with current values
      si->getStateSpace()->populateState(candidateState, values);

      // denseDB_->debugState(candidateState);

      // Collision check
      if (!si->isValid(candidateState))
      {
        // OMPL_ERROR("Found a state that is not valid! ");

        // Visualize
        // if (visualizeGridGeneration_)
        // {
        //   // Candidate node has already (just) been added
        //   visual_->vizState(threadID+1, candidateState, /*red arm*/ 3, 1);
        //   visual_->vizTrigger(threadID+1);
        // }

        continue;
      }

      // Add vertex to graph
      VertexType type = START;  // TODO(davetcoleman): type START is dummy

      // Allocate state before mutex
      base::State *newState = si->cloneState(candidateState);

      {
        boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
        denseDB_->addVertex(newState, type);
      }

      // Visualize
      if (visualizeGridGeneration_)
      {
        // Candidate node has already (just) been added
        // visual_->viz1State(candidateState, /*red arrow*/ 6, 1);
        // visual_->viz1State(candidateState, /*green arm*/ 1, 1);
        // if (threadID == 0)
        // {
        //   visual_->viz1Trigger();
        // }

        // Visualize
        // Candidate node has already (just) been added
        // if (threadID < 6)
        // {
        //   visual_->vizState(threadID + 1, candidateState, tools::SMALL, tools::GREEN, 1);
        //   visual_->vizTrigger(threadID + 1);
        // }
      }
    }
  }
}

std::size_t Discretizer::getEdgesPerVertex(base::SpaceInformationPtr si)
{
  // in 2D this creates the regular square with diagonals of 8 edges
  if (si->getStateSpace()->getDimension() == 3)
  {
    return 8;
  }

  // full robot
  return si->getStateSpace()->getDimension() * 2;
}

void Discretizer::generateEdges()
{
  const bool verbose = false;

  std::vector<boost::thread *> threads(numThreads_);
  OMPL_INFORM("Generating edges using %u threads", numThreads_);

  std::size_t verticesPerThread = denseDB_->getNumVertices() / numThreads_;  // rounds down
  std::size_t startVertex = 0;
  std::size_t endVertex;
  std::size_t errorCheckTotalVertices = 0;
  numEdgesInCollision_ = 0;  // stats for later use

  // Cache certain values
  getVertexNeighborsPreprocess();

  // For each thread
  for (std::size_t threadID = 0; threadID < threads.size(); ++threadID)
  {
    endVertex = startVertex + verticesPerThread - 1;

    // Check if this is the last thread
    if (threadID == threads.size() - 1)
    {
      // have it do remaining vertices to check
      endVertex = denseDB_->getNumVertices() - 1;
    }
    errorCheckTotalVertices += (endVertex - startVertex);

    if (verbose)
      std::cout << "Thread " << threadID << " has vertices from " << startVertex << " to " << endVertex << std::endl;

    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    threads[threadID] =
        new boost::thread(boost::bind(&Discretizer::generateEdgesThread, this, threadID, startVertex, endVertex, si));
    startVertex += verticesPerThread;
  }

  // Join threads
  for (std::size_t threadID = 0; threadID < threads.size(); ++threadID)
  {
    threads[threadID]->join();
    delete threads[threadID];
  }

  // Error check
  if (errorCheckTotalVertices == denseDB_->getNumVertices())
  {
    OMPL_ERROR("Incorrect number of vertices were processed for edge creation");
    exit(-1);
  }

  OMPL_INFORM("Generated %i edges. Finished generating grid.", denseDB_->getNumEdges());
}

void Discretizer::generateEdgesThread(std::size_t threadID, DenseVertex startVertex, DenseVertex endVertex,
                                      base::SpaceInformationPtr si)
{
  const bool verbose = false;

  std::size_t feedbackFrequency = std::max(1.0, double(endVertex - startVertex) / 100.0);

  // Nearest Neighbor search
  std::vector<DenseVertex> graphNeighborhood;

  // Stats
  std::size_t numEdgesInCollision = 0;

  // Process [startVertex, endVertex] inclusive
  for (DenseVertex v1 = startVertex; v1 <= endVertex; ++v1)
  {
    // Skip the query vertex (first vertex)
    if (v1 <= denseDB_->queryVertices_.back())
      continue;

    // User feedback on thread 0
    if (threadID == numThreads_ - 1 && v1 % feedbackFrequency == 0)
    {
      std::cout << "Generating edges progress: " << std::setprecision(1)
                << (v1 - startVertex) / static_cast<double>(endVertex - startVertex) * 100.0
                << " % Total edges: " << denseDB_->getNumEdges() << " Cache size: " << denseCache_->getCacheSize()
                << " Cache usage: " << denseCache_->getPercentCachedCollisionChecks() << "%" << std::endl;
    }

    // Add edges
    graphNeighborhood.clear();

    // Get neighors with one of the strategies
    getVertexNeighbors(v1, graphNeighborhood);

    if (verbose)
      OMPL_INFORM("Found %u neighbors", graphNeighborhood.size());

    // For each nearby vertex, add an edge
    std::size_t errorCheckNumSameVerticies = 0;  // sanity check
    for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
    {
      if (verbose)
        OMPL_INFORM("v2 = %u", i);

      DenseVertex &v2 = graphNeighborhood[i];

      // Check if these vertices are the same
      if (v1 == v2)
      {
        errorCheckNumSameVerticies++;  // sanity check
        continue;
      }

      // Check if these vertices already share an edge
      {
        boost::unique_lock<boost::mutex> scoped_lock(edgeMutex_);  // TODO make this a read-only mutex
        if (boost::edge(v1, v2, denseDB_->g_).second)
          continue;
      }

      // Check edge for collision
      // if (!si->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]))
      if (!denseCache_->checkMotionWithCache(v1, v2, threadID))
      {
        numEdgesInCollision++;
        continue;
      }

      // Determine cost for edge depending on mode
      double cost;
      if (denseDB_->popularityBiasEnabled_)
      {
        cost = MAX_POPULARITY_WEIGHT;
      }
      else
      {
        cost = denseDB_->distanceFunction(v1, v2);
      }

      // Create edge
      DenseEdge e;
      {
        boost::unique_lock<boost::mutex> scoped_lock(edgeMutex_);
        e = denseDB_->addEdge(v1, v2, cost, FREE);
      }
    }  // for each v2

    // Make sure one and only one vertex is returned from the NN search that is the same as parent vertex
    assert(errorCheckNumSameVerticies == 1);

  }  // for each v1

  // Re-use mutex (we could create a new one though)
  {
    boost::unique_lock<boost::mutex> scoped_lock(edgeMutex_);
    numEdgesInCollision_ += numEdgesInCollision;
  }
}

void Discretizer::getVertexNeighborsPreprocess()
{
  edgeConnectionStrategy_ = 1;

  std::cout << std::fixed << std::setprecision(2);
  switch (edgeConnectionStrategy_)
  {
    case 1:
      findNearestKNeighbors_ = getEdgesPerVertex(si_);
      std::cout << "Edge Connection Strategy: findNearestKNeighbors k=" << findNearestKNeighbors_ << std::endl;
      break;
    case 2:
      radiusNeighbors_ = sqrt(2 * (discretization_ * discretization_));
      std::cout << "Edge Connection Strategy: radiusNeighbors_ r=" << radiusNeighbors_ << std::endl;
      break;
    case 3:
    {
      // Setup Method 3
      double kPRMConstant_ = boost::math::constants::e<double>() +
                             (boost::math::constants::e<double>() / (double)si_->getStateSpace()->getDimension());
      findNearestKNeighbors_ = static_cast<unsigned int>(ceil(kPRMConstant_ * log((double)denseDB_->getNumVertices())));
      std::cout << "Edge Connection Strategy: findNearestKNeighbors k=" << findNearestKNeighbors_ << std::endl;
    }
    break;
    default:
      OMPL_ERROR("Incorrect edge connection stragety");
  }
}

void Discretizer::getVertexNeighbors(base::State *state, std::vector<DenseVertex> &graphNeighborhood, std::size_t threadID)
{
  denseDB_->stateProperty_[denseDB_->queryVertices_[threadID]] = state;
  getVertexNeighbors(denseDB_->queryVertices_[threadID], graphNeighborhood);
  denseDB_->stateProperty_[denseDB_->queryVertices_[threadID]] = nullptr;
}

void Discretizer::getVertexNeighbors(DenseVertex v1, std::vector<DenseVertex> &graphNeighborhood)
{
  const std::size_t numSameVerticiesFound = 1;  // add 1 to the end because the NN tree always returns itself

  // QUESTION: How many edges should each vertex connect with?
  switch (edgeConnectionStrategy_)
  {
    case 1:
      // METHOD 1
      denseDB_->nn_->nearestK(v1, findNearestKNeighbors_ + numSameVerticiesFound, graphNeighborhood);
      break;
    case 2:
      // METHOD 2
      denseDB_->nn_->nearestR(v1, radiusNeighbors_, graphNeighborhood);
      break;
    case 3:
      // METHOD 3 - based on k-PRM*
      denseDB_->nn_->nearestK(v1, findNearestKNeighbors_ + numSameVerticiesFound, graphNeighborhood);

      break;
    default:
      OMPL_ERROR("Incorrect edge connection stragety");
  }
}

void Discretizer::eliminateDisjointSets()
{
  CALLGRIND_TOGGLE_COLLECT;

  OMPL_INFORM("Eliminating disjoint sets in Dense Graph");

  // Statistics
  eliminateDisjointSetsVerticesAdded_ = 0;
  eliminateDisjointSetsEdgesAdded_ = 0;
  eliminateDisjointSetsVerticesAddedUnsaved_ = 0;
  stopSearchingDisjointSets_ = false;

  getVertexNeighborsPreprocess();    // prepare the constants
  denseDB_->getSparseDB()->setup();  // make sure sparse delta is chosen

  // Setup threading
  std::vector<boost::thread *> threads(numThreads_);
  OMPL_INFORM("Sampling to eliminate disjoint sets using %u threads", numThreads_);

  // For each thread
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    bool verboseThread = false;  //!i; // only thread 0 is verbose

    threads[i] =
      new boost::thread(boost::bind(&Discretizer::eliminateDisjointSetsThread, this, i, si, verboseThread));
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }

  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_DUMP_STATS;
}

void Discretizer::eliminateDisjointSetsThread(std::size_t threadID, base::SpaceInformationPtr si, bool verbose)
{
  // Add dense vertex
  base::State *candidateState = si_->allocState();

  // For each dense vertex we add
  std::size_t numSets = 2;  // dummy value that will be updated at first loop
  std::size_t sampleCount = 0;
  std::size_t noVisibleNeighborCount = 0;
  std::vector<DenseVertex> graphNeighborhood;
  std::vector<DenseVertex> visibleNeighborhood;

  base::ValidStateSamplerPtr sampler = si_->allocValidStateSampler();
  while (numSets > 1 && !stopSearchingDisjointSets_)
  {
    bool sampleAdded = false;

    while (!sampleAdded)  // For each random sample
    {
      graphNeighborhood.clear();
      visibleNeighborhood.clear();

      // Sample randomly
      if (!sampler->sample(candidateState))
      {
        OMPL_WARN("Sampler failed to find valid state");
        continue;
      }
      sampleCount++;

      // Collision check TODO - this should not be necessary but there is a big with sampleNear() that does not
      // always return a valid state
      if (!si->isValid(candidateState))
      {
        //OMPL_ERROR("Sampled 'valid' state that is actually not valid");
        continue;
      }

      // Get neighbors
      //{
      //boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
      getVertexNeighbors(candidateState, graphNeighborhood, threadID);
      //}

      // Now that we got the neighbors from the NN, find the ones that are visible
      for (DenseVertex &denseV : graphNeighborhood)
      {
        //BOOST_ASSERT_MSG(denseV != denseDB_->queryVertex_, "Query vertex should not be in the graph neighborhood");

        if (si_->checkMotion(candidateState, denseDB_->stateProperty_[denseV]))
        {
          visibleNeighborhood.push_back(denseV);
        }
      }

      // Debug
      if (verbose && false)
        std::cout << "Sample #" << sampleCount << " graphNeighborhood " << graphNeighborhood.size()
                  << " visibleNeighborhood: " << visibleNeighborhood.size()
                  << " noVisibleNeighborPercent: " << noVisibleNeighborCount / double(sampleCount) * 100.0 << "%"
                  << std::endl;

      // If no neighbors, add the vertex
      if (visibleNeighborhood.empty())
      {
        noVisibleNeighborCount++;
        // std::cout << threadID << ": Adding vertex because no neighbors" << std::endl;

        {
          boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);

          if (false)  // Visualize
          {
            visual_->viz6State(candidateState, tools::ROBOT, tools::CYAN, 0);
            visual_->viz6Trigger();
          }

          denseDB_->addVertex(si_->cloneState(candidateState), COVERAGE);
          denseDB_->setGraphUnsaved();
          eliminateDisjointSetsVerticesAdded_++;
          eliminateDisjointSetsVerticesAddedUnsaved_++;
        }

        // Record this new addition
        denseDB_->graphUnsaved_ = true;

        // no need to check disjoint states because we know it could not have changed with just a vertex addition
        continue;
      }

      // Check each pair of neighbors for connectivity
      for (std::size_t i = 0; i < visibleNeighborhood.size(); ++i)
      {
        for (std::size_t j = i + 1; j < visibleNeighborhood.size(); ++j)
        {
          // If they are in different components
          if (!denseDB_->sameComponent(visibleNeighborhood[i], visibleNeighborhood[j]))
          {
            // std::cout << threadID << ": Neighbors " << visibleNeighborhood[i] << ", " << visibleNeighborhood[j] << "
            // are in different components, add!" << std::endl;

            // Attempt to connect new Dense vertex into dense graph by connecting neighbors
            connectNewVertex(si_->cloneState(candidateState), visibleNeighborhood, verbose);

            sampleAdded = true;
            break;
          }
        }  // for each neighbor
        if (sampleAdded)
          break;
      }  // for each neighbor

      // Thread 0 occationally saves
      if (threadID == 0 && eliminateDisjointSetsVerticesAddedUnsaved_ >= 50)
      {
        boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);

        OMPL_INFORM("Saving database so not to lose progress...");
        denseDB_->saveIfChanged();
        eliminateDisjointSetsVerticesAddedUnsaved_ = 0;
      }
    }  // while sampling unuseful states

    // Update number of sets
    numSets = denseDB_->getDisjointSetsCount();

    // Debug
    OMPL_INFORM("DenseSampler Thread %u: Verticies added: %u, Edges added: %u, Remaining disjoint sets: %u", threadID,
                eliminateDisjointSetsVerticesAdded_, eliminateDisjointSetsEdgesAdded_, numSets);

    // Occationally do stuff just in thread 0
    if (threadID == 0)
    {
      DisjointSetsParentKey disjointSets;
      denseDB_->getDisjointSets(disjointSets);
      denseDB_->printDisjointSets(disjointSets);
      //stopSearchingDisjointSets_ = true;
    }

  }  // end while

  // Free memory
  si_->freeState(candidateState);
}

void Discretizer::connectNewVertex(base::State *state, std::vector<DenseVertex> visibleNeighborhood, bool verbose)
{
  boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);

  DenseVertex v1 = denseDB_->addVertex(state, COVERAGE);  // TODO VertexType is meaningless
  denseDB_->setGraphUnsaved();
  eliminateDisjointSetsVerticesAdded_++;
  eliminateDisjointSetsVerticesAddedUnsaved_++;

  // Visualize new vertex
  // if (visualizeAddSample_)
  // {
  //   visual_->viz1State(state, /*mode=*/1, 1);
  // }

  // For each visible neighbor vertex, add an edge
  std::size_t numEdgesAdded = 0;  // sanity check
  // for (std::size_t i = 0; i < graphNeighborhood.size(); ++i)
  for (DenseVertex &v2 : visibleNeighborhood)
  {
    // Check if these vertices are the same STATE
    if (si_->getStateSpace()->equalStates(state, denseDB_->stateProperty_[v2]))
    {
      OMPL_ERROR("This state has already been added, this is low probabilty event TODO");
      exit(-1);
    }

    denseDB_->addEdge(v1, v2, 0);  // TODO cost... desiredAverageCost_);
    numEdgesAdded++;
    eliminateDisjointSetsEdgesAdded_++;

    // if (visualizeAddSample_)  // Debug: display edge
    // {
    //   double popularity = 100;  // TODO: maybe make edge really popular so we can be sure its added to the
    //                             // spars graph since we need it
    //   visual_->viz1Edge(state, stateProperty_[v2], popularity);
    // }
  }  // for each neighbor

  // Make sure one and only one vertex is returned from the NN search that is the same as parent vertex
  BOOST_ASSERT_MSG(numEdgesAdded >= 2, "Too few edges added from new DenseVertex connectivity node");

  if (verbose)
    std::cout << "Connected new vertex to " << numEdgesAdded << " neighbors" << std::endl;

  // Visualize
  // if (visualizeAddSample_)
  // {
  //   visual_->viz1Trigger();
  //   usleep(0.001 * 1000000);
  // }

  // Record this new addition
  denseDB_->graphUnsaved_ = true;
}

}  // namespace
}  // namespace
}  // namespace
