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
Discretizer::Discretizer(base::SpaceInformationPtr si, DenseDBPtr denseDB, base::VisualizerPtr visual)
  : si_(si), denseDB_(denseDB), visual_(visual)
{
}

Discretizer::~Discretizer(void)
{
}

void Discretizer::generateGrid()
{
  OMPL_INFORM("Generating grid");

  if (!si_->isSetup())
  {
    OMPL_WARN("Space information setup was not yet called. Calling now.");
    si_->setup();
  }

  double vertexDuration;
  {
    // Benchmark runtime
    time::point start_time = time::now();

    // Create vertices
    generateVertices();

    // Benchmark runtime
    vertexDuration = time::seconds(time::now() - start_time);
  }

  OMPL_INFORM("Generated %i vertices.", denseDB_->getNumVertices());

  // Error check
  if (denseDB_->getNumVertices() < 2)
  {
    OMPL_ERROR("No vertices generated, failing");
    exit(-1);
  }

  double edgeDuration;
  {
    // Benchmark runtime
    time::point start_time = time::now();

    // Create edges ----------------------------------------
    generateEdges();

    // Benchmark runtime
    edgeDuration = time::seconds(time::now() - start_time);
  }

  // Get the average vertex degree (number of connected edges)
  double averageDegree = (denseDB_->getNumEdges() * 2) / static_cast<double>(denseDB_->getNumVertices());

  // Check how many disjoint sets are in the dense graph (should be none)
  std::size_t numSets = denseDB_->checkConnectedComponents();

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Discretization stats:");
  OMPL_INFORM("   Total valid vertices:   %u", denseDB_->getNumVertices());
  OMPL_INFORM("   Vertex generation time: %f seconds (%f min)", vertexDuration, vertexDuration / 60.0);
  OMPL_INFORM("   Total valid edges:      %u", denseDB_->getNumEdges());
  OMPL_INFORM("   Num edges colliding:    %f", numEdgesInCollision_);
  OMPL_INFORM("   Edge generation time:   %f seconds (%f min)", edgeDuration, edgeDuration / 60.0);
  OMPL_INFORM("   Average degree:         %f", averageDegree);
  OMPL_INFORM("   Connected Components:   %u", numSets);
  OMPL_INFORM("------------------------------------------------------");

  // Display
  // if (visualizeGridGeneration_)
  // visual_->viz1Trigger();
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
  std::size_t numThreads = boost::thread::hardware_concurrency();

  // Check that we have enough jointIncrements for all the threads
  if (jointIncrements < numThreads)
  {
    OMPL_WARN("There are fewer joint_0 increments (%u) at current discretization (%f) than available threads (%u), "
              "underutilizing threading",
              jointIncrements, discretization_, numThreads);
    OMPL_INFORM("Optimal discretization: %f", range / double(numThreads));
    numThreads = jointIncrements;
  }

  // Debugging
  if (false)
  {
    OMPL_WARN("Overriding number of threads for testing to 1");
    numThreads = 1;
  }

  // Warn
  //if (numThreads > 1 && visualizeGridGeneration_)
  //OMPL_WARN("Multiple threads are trying to visualize, could cause race conditions");

  // Setup threading
  std::vector<boost::thread *> threads(numThreads);
  OMPL_INFORM("Generating vertices using %u threads", numThreads);

  std::size_t jointIncrementsPerThread = jointIncrements / numThreads;

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
    std::cout << "  NumThreads: " << numThreads << std::endl;
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
        new boost::thread(boost::bind(&Discretizer::createVertexThread, this, i, startJointValue, endJointValue, si));

    startJointValue = endJointValue;
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }
}

void Discretizer::createVertexThread(std::size_t threadID, double startJointValue, double endJointValue,
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
  else if (dim == 6)
  {
    maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
    values[5] = 0.0;             // middle rotation of wrist
  }
  else if (dim == 12)
  {
    maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
    values[5] = 0.0;             // middle rotation of wrist
    values[11] = 0.0;             // middle rotation of wrist
  }
  else
  {
    maxDiscretizationLevel = dim - 1;
  }

  // Loop through current joint
  for (double value = startJointValue; value < endJointValue; value += discretization_)
  {
    // User feedback on thread 0
    if (threadID == 0)
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
      else // regular treatment
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
      GuardType type = START;  // TODO(davetcoleman): type START is dummy

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
        //visual_->viz1State(candidateState, /*red arrow*/ 6, 1);
        //visual_->viz1State(candidateState, /*green arm*/ 1, 1);
        // if (threadID == 0)
        // {
        //   visual_->viz1Trigger();
        // }

        // Visualize
        if (visualizeGridGeneration_)
        {
          // Candidate node has already (just) been added
          if (threadID < 6)
          {
            visual_->vizState(threadID+1, candidateState, /*green arm*/ 1, 1);
            visual_->vizTrigger(threadID+1);
          }
        }
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

  // Setup threading
  std::size_t numThreads = boost::thread::hardware_concurrency();

  // Debugging
  if (false)
  {
    OMPL_WARN("Overriding number of threads for testing to 1");
    numThreads = 1;
  }

  std::vector<boost::thread *> threads(numThreads);
  OMPL_INFORM("Generating edges using %u threads", numThreads);

  std::size_t verticesPerThread = denseDB_->getNumVertices() / numThreads;  // rounds down
  std::size_t startVertex = 0;
  std::size_t endVertex;
  std::size_t errorCheckTotalVertices = 0;
  numEdgesInCollision_ = 0; // stats for later use

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

    threads[threadID] = new boost::thread(
        boost::bind(&Discretizer::generateEdgesThread, this, threadID, startVertex, endVertex, si));
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

  std::size_t feedbackFrequency = (endVertex - startVertex) / 10;

  // Nearest Neighbor search
  std::vector<DenseVertex> graphNeighborhood;
  std::vector<DenseEdge> unvalidatedEdges;

  // Stats
  std::size_t numEdgesInCollision = 0;

  // Process [startVertex, endVertex] inclusive
  for (DenseVertex v1 = startVertex; v1 <= endVertex; ++v1)
  {
    // Skip the query vertex (first vertex)
    if (v1 == denseDB_->queryVertex_)
      continue;

    // User feedback on thread 0
    if (threadID == 0 && (v1) % feedbackFrequency == 0)
    {
      std::cout << "Generating edges progress: " << std::setprecision(1)
                << (v1 - startVertex) / static_cast<double>(endVertex - startVertex) * 100.0 << " %" << std::endl;
    }

    // Add edges
    graphNeighborhood.clear();

    // Get neighors with one of the strategies
    {
      boost::unique_lock<boost::mutex> scoped_lock(edgeNnMutex_);
      getVertexNeighbors(v1, graphNeighborhood);
    }

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
        boost::unique_lock<boost::mutex> scoped_lock(edgeMutex_); // TODO make this a read-only mutex
        if (boost::edge(v1, v2, denseDB_->g_).second)
          continue;
      }

      // Check edge for collision
      if (!si->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]))
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

  } // for each v1

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
      std::cout << "findNearestKNeighbors: " << findNearestKNeighbors_ << std::endl;
      break;
    case 2:
      radiusNeighbors_ = sqrt(2 * (discretization_ * discretization_));
      std::cout << "radiusNeighbors_: " << radiusNeighbors_ << std::endl;
      break;
    case 3:
    {
      // Setup Method 3
      double kPRMConstant_ = boost::math::constants::e<double>() +
                             (boost::math::constants::e<double>() / (double)si_->getStateSpace()->getDimension());
      findNearestKNeighbors_ = static_cast<unsigned int>(ceil(kPRMConstant_ * log((double)denseDB_->getNumVertices())));
      std::cout << "findNearestKNeighbors: " << findNearestKNeighbors_ << std::endl;
    }
    break;
    default:
      OMPL_ERROR("Incorrect edge connection stragety");
  }
}

void Discretizer::getVertexNeighbors(DenseVertex v1, std::vector<DenseVertex> &graphNeighborhood)
{
  const std::size_t numSameVerticiesFound = 1;  // add 1 to the end because the NN tree always returns itself

  //std::cout << "getVertexNeighbors: " << v1 << std::endl;

  // Search
  denseDB_->stateProperty_[denseDB_->queryVertex_] = denseDB_->stateProperty_[v1];

  // QUESTION: How many edges should each vertex connect with?
  switch (edgeConnectionStrategy_)
  {
    case 1:
      // METHOD 1
      denseDB_->nn_->nearestK(denseDB_->queryVertex_, findNearestKNeighbors_ + numSameVerticiesFound,
                              graphNeighborhood);
      break;
    case 2:
      // METHOD 2
      denseDB_->nn_->nearestR(denseDB_->queryVertex_, radiusNeighbors_, graphNeighborhood);
      break;
    case 3:
      // METHOD 3 - based on k-PRM*
      denseDB_->nn_->nearestK(denseDB_->queryVertex_, findNearestKNeighbors_ + numSameVerticiesFound,
                              graphNeighborhood);
      break;
    default:
      OMPL_ERROR("Incorrect edge connection stragety");
  }
  // Set search vertex to NULL to prevent segfault on class unload of memory
  denseDB_->stateProperty_[denseDB_->queryVertex_] = NULL;
}

}  // namespace
}  // namespace
}  // namespace
