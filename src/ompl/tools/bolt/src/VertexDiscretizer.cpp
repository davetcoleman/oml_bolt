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
#include <ompl/base/State.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/tools/bolt/Debug.h>

// OMPL Bolt
#include <ompl/tools/bolt/VertexDiscretizer.h>
#include <ompl/tools/bolt/SparseGraph.h>

// Boost
#include <boost/foreach.hpp>
#include <boost/thread.hpp>

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
VertexDiscretizer::VertexDiscretizer(SparseGraphPtr sg) : sg_(sg)
{
  // Copy pointers to important datastructures
  si_ = sg_->getSpaceInformation();
  visual_ = sg_->getVisual();

  // Choose number of threads
  numThreads_ = boost::thread::hardware_concurrency();
}

VertexDiscretizer::~VertexDiscretizer()
{
  freeMemory();
}

void VertexDiscretizer::freeMemory()
{
  for (base::State *state : failedStates_)
  {
    si_->freeState(state);
  }
  failedStates_.clear();
}

bool VertexDiscretizer::generateLattice(std::size_t indent)
{
  BOLT_FUNC(indent, true, "generateLattice()");

  if (numThreads_ > 1 && visualizeGridGeneration_)
  {
    OMPL_WARN("Visualizing in non-thread-safe manner. Auto reduced to 1 thread for debug mode");
    numThreads_ = 1;
  }

  if (discretization_ < std::numeric_limits<double>::epsilon())
  {
    OMPL_WARN("Discretization not set");
    exit(-1);
  }

  freeMemory();

  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID];
  const std::size_t jointIncrements = floor(range / discretization_);
  double leftOver = range - jointIncrements * discretization_;
  double startOffset = leftOver / 2;

  BOLT_INFO(indent, verbose_, "------------------------------------------");
  BOLT_INFO(indent, verbose_, "Multi-Pass Discretization");
  BOLT_INFO(indent, verbose_, "  Discretization:       " << discretization_);
  BOLT_INFO(indent, verbose_, "  High Bound:           " << bounds.high[jointID]);
  BOLT_INFO(indent, verbose_, "  Low Bound:            " << bounds.low[jointID]);
  BOLT_INFO(indent, verbose_, "  Range:                " << range);
  BOLT_INFO(indent, verbose_, "  Joint Increments:     " << jointIncrements);
  BOLT_INFO(indent, verbose_, "  Left Over:            " << leftOver);
  BOLT_INFO(indent, verbose_, "  Start Offset:         " << startOffset);
  BOLT_INFO(indent, verbose_, "------------------------------------------");

  // base::State *s1;
  // base::State *s2;
  // base::State *s3;

  // Create two levels of grids TODO remove this for loop
  for (std::size_t i = 0; i < 1; ++i)
  {
    BOLT_INFO(indent, verbose_, "Discretize iteration " << i);

    // Set starting value offset
    if (i == 0)
      startingValueOffset_ = startOffset;
    else
      startingValueOffset_ = startOffset + discretization_ / 2.0;

    // Generate vertices
    generateGrid(indent);
  }

  BOLT_INFO(indent, true, "VertexDiscretizer found " << failedStates_.size() << " failed vertices");

  return true;
}

bool VertexDiscretizer::generateGrid(std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "VertexDiscretizer::generateGrid()");

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
    generateVertices(indent);

    // Benchmark runtime
    vertexDuration = time::seconds(time::now() - startTime);
  }
  BOLT_INFO(indent, verbose_, "Generated " << failedStates_.size() << " vertices in " << vertexDuration << " sec");

  // Error check
  if (failedStates_.size() < 2)
  {
    OMPL_WARN("No vertices generated, failing");
    //exit(-1);
  }

  // Total benchmark runtime
  double totalDuration = time::seconds(time::now() - totalStartTime);

  BOLT_INFO(indent, vThread_, "------------------------------------------------------");
  BOLT_INFO(indent, vThread_, "Vertex Discretization stats:");
  BOLT_INFO(indent, vThread_, "   Total valid vertices:   " << failedStates_.size());
  BOLT_INFO(indent, vThread_, "   Vertex generation time: " << vertexDuration);
  BOLT_INFO(indent, vThread_, "   Total grid gen. time:   " << totalDuration);
  BOLT_INFO(indent, vThread_, "------------------------------------------------------");

  return true;
}

void VertexDiscretizer::generateVertices(std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "VertexDiscretizer::generateVertices()");

  std::size_t dim = si_->getStateSpace()->getDimension();

  // Setup bounds for joint 0
  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  assert(bounds.high.size() == bounds.low.size());
  assert(bounds.high.size() == dim);

  // Divide joint 0 between threads
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID] - startingValueOffset_;
  const std::size_t jointIncrements = ceil(range / discretization_);

  assert(jointIncrements >= 1);

  // Check that we have enough jointIncrements for all the threads
  if (jointIncrements < numThreads_)
  {
    BOLT_DEBUG(indent, verbose_, "There are fewer joint_0 increments ("
               << jointIncrements << ") at current discretization (" << discretization_
               << ") than available threads (" << numThreads_ << "), underutilizing threading");
    BOLT_DEBUG(indent, verbose_, "Optimal discretization for thread count: " << range / double(numThreads_));
    numThreads_ = jointIncrements;
  }
  std::size_t jointIncrementsPerThread = jointIncrements / numThreads_;

  BOLT_DEBUG(indent, verbose_, "-------------------------------------------------------");
  BOLT_DEBUG(indent, verbose_, "Single-Pass Discretization: ");
  BOLT_DEBUG(indent, verbose_, "  Dimensions:             " << dim);
  BOLT_DEBUG(indent, verbose_, "  Discretization:         " << discretization_);
  BOLT_DEBUG(indent, verbose_, "  J0 Low:                 " << bounds.low[jointID]);
  BOLT_DEBUG(indent, verbose_, "  J0 High:                " << bounds.high[jointID]);
  BOLT_DEBUG(indent, verbose_, "  J0 Range:               " << range);
  BOLT_DEBUG(indent, verbose_, "  J0 Increments:          " << jointIncrements);
  BOLT_DEBUG(indent, verbose_, "  J0 IncrementsPerThread: " << jointIncrementsPerThread);
  BOLT_DEBUG(indent, verbose_, "  Total states:           " << pow(jointIncrements, dim));
  BOLT_DEBUG(indent, verbose_, "  Num Threads:            " << numThreads_);
  BOLT_DEBUG(indent, verbose_, "-------------------------------------------------------");

  // Setup threading
  std::vector<boost::thread *> threads(numThreads_);

  double startJointValue = bounds.low[jointID] + startingValueOffset_;
  double endJointValue;

  // For each thread
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    endJointValue = startJointValue + (jointIncrementsPerThread - 1) * discretization_;

    // Check if this is the last thread
    if (i == threads.size() - 1)
    {
      // have it do remaining bounds
      endJointValue = bounds.high[jointID];
    }

    BOLT_DEBUG(indent, vThread_, "Thread " << i << " will process Joint 0 values from " << startJointValue << " to "
               << endJointValue);

    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    threads[i] = new boost::thread(
                                   boost::bind(&VertexDiscretizer::generateVerticesThread, this, i, startJointValue, endJointValue, si, indent));

    startJointValue = endJointValue + discretization_;
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }
}

void VertexDiscretizer::generateVerticesThread(std::size_t threadID, double startJointValue, double endJointValue,
                                               base::SpaceInformationPtr si, std::size_t indent)
{
  BOLT_FUNC(indent, vThread_, "generateVerticesThread()");

  std::size_t jointID = 0;
  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();
  base::State *candidateState = si->getStateSpace()->allocState();

  // Prepare for recursion
  std::vector<double> values(si->getStateSpace()->getDimension(), 0);

  // Customize for different state spaces TODO more generic
  const std::size_t dim = si->getStateSpace()->getDimension();
  std::size_t maxDiscretizationLevel = dim - 1;

  // if (dim == 6 && false)
  // {
  //   maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
  //   values[5] = 0.0;                   // middle rotation of wrist
  // }
  // else if (dim == 12 && false)
  // {
  //   maxDiscretizationLevel = dim - 2;  // don't discretize the wrist rotation
  //   values[5] = 0.0;                   // middle rotation of wrist
  //   values[11] = 0.0;                  // middle rotation of wrist
  // }

  // Loop through current joint
  for (double value = startJointValue; value <= endJointValue; value += discretization_)
  {
    // User feedback on thread 0
    if (threadID == numThreads_ - 1)
    {
      const double percent = (value - startJointValue) / (endJointValue - startJointValue) * 100.0;
      BOLT_DEBUG(indent, vThread_,
                 "Vertex generation progress: " << std::setprecision(1) << percent << " % Total vertices: "
                 << (failedStates_.size() > 0 ? failedStates_.size() - 1 : 0));
    }

    // Set first joint value
    values[jointID] = value;

    // Keep recursing
    recursiveDiscretization(threadID, values, jointID + 1, si, candidateState, maxDiscretizationLevel, indent);
  }

  // Cleanup
  si->freeState(candidateState);
}

void VertexDiscretizer::recursiveDiscretization(std::size_t threadID, std::vector<double> &values, std::size_t jointID,
                                                base::SpaceInformationPtr si, base::State *candidateState,
                                                std::size_t maxDiscretizationLevel, std::size_t indent)
{
  BOLT_FUNC(indent, vThread_, "recursiveDiscretization()");

  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();

  // Error check
  BOOST_ASSERT_MSG(jointID < values.size(), "Joint ID too high");
  BOOST_ASSERT_MSG(bounds.high[jointID] - bounds.low[jointID] > discretization_, "Bounds too small");

  // Loop through current joint
  for (double value = bounds.low[jointID] + startingValueOffset_; value <= bounds.high[jointID];
       value += discretization_)
  {
    // User feedback on thread 0 for high dimension spaces
    if (jointID == 1 && threadID == 0 && si->getStateSpace()->getDimension() > 3)
    {
      const double percent = (value - bounds.low[jointID]) / (bounds.high[jointID] - bounds.low[jointID]) * 100.0;

      BOLT_DEBUG(indent, vThread_,
                 "Level 1 vertex generation progress: " << std::setprecision(1) << percent << " % Total vertices: "
                 << (failedStates_.size() > 0 ? failedStates_.size() - 1 : 0));
    }

    // Set value
    values[jointID] = value;

    // Check if we are at the end of the recursion
    BOLT_DEBUG(indent, vThread_, "threadID " << threadID << " jointID " << jointID << " of " << maxDiscretizationLevel
               << " --- ");
    BOLT_DEBUG(indent, vThread_, "value: " << value << " high: " << bounds.high[jointID]
               << " low: " << bounds.low[jointID] << " disc: " << discretization_
               << " values.size() " << values.size());
    if (vThread_)
    {
      std::copy(values.begin(), values.end(), std::ostream_iterator<double>(std::cout, ", "));
      std::cout << std::endl;
    }

    if (jointID < maxDiscretizationLevel)
    {
      // Keep recursing

      // Special rule for 12dof
      if (si->getStateSpace()->getDimension() == 12 && jointID == 4)
      {
        // skip joint id 5 (joint 6)
        recursiveDiscretization(threadID, values, jointID + 2, si, candidateState, maxDiscretizationLevel, indent);
      }
      else  // regular treatment
      {
        recursiveDiscretization(threadID, values, jointID + 1, si, candidateState, maxDiscretizationLevel, indent);
      }
    }
    else  // this is the end of recursion, create a new state
    {
      createState(threadID, values, si, candidateState, indent);
    }
  }
}

void VertexDiscretizer::createState(std::size_t threadID, std::vector<double> &values, base::SpaceInformationPtr si,
                                    base::State *candidateState, std::size_t indent)
{
  BOLT_FUNC(indent, vThread_, "createState()");

  // Fill the state with current values
  si->getStateSpace()->populateState(candidateState, values);

  // Collision check
  double dist;
  if (!si->getStateValidityChecker()->isValid(candidateState, dist))
  {
    BOLT_ERROR(indent, vThread_, "Rejected because of validity");

    // Visualize
    if (visualizeGridGeneration_)
    {
      // Candidate node rejected
      visual_->viz1()->state(candidateState, LARGE, RED, 0);
      visual_->viz1()->state(candidateState, ROBOT, RED, 0);
      visual_->viz1()->trigger();

      if (visualizeGridGenerationWait_)
        visual_->waitForUserFeedback("rejected");
      else
        usleep(0.001 * 1000000);
    }

    // // Allocate state before mutex
    // base::State *newState = si->cloneState(candidateState);
    // {
    //   boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
    //   failedStates_.push_back(newState);
    // }

    return;
  }

  if (dist < clearance_)
  {
    BOLT_WARN(indent, vThread_, "Rejected because of clearance " << dist << " required: " << clearance_);

    // Visualize
    if (visualizeGridGeneration_)
    {
      // Candidate node rejected
      visual_->viz1()->state(candidateState, LARGE, YELLOW, 0);
      visual_->viz1()->state(candidateState, ROBOT, YELLOW, 0);
      visual_->viz1()->trigger();

      if (visualizeGridGenerationWait_)
        visual_->waitForUserFeedback("clearance");
      else
        usleep(0.001 * 1000000);
    }

    // // Allocate state before mutex
    // base::State *newState = si->cloneState(candidateState);
    // {
    //   boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
    //   failedStates_.push_back(newState);
    // }

    return;
  }

  BOLT_GREEN_DEBUG(indent, vThread_, "Accepted, clearance: " << dist);

  // Add to graph
  {
    boost::unique_lock<boost::mutex> scoped_lock(sparseGraphMutex_);
    sg_->addVertex(si->cloneState(candidateState), DISCRETIZED, indent);
  }

  // Allocate state before mutex
  // base::State *newState = si->cloneState(candidateState);
  // {
  //   boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
  //   // Add to vector
  //   failedStates_.push_back(newState);
  // }

  // Visualize
  if (visualizeGridGeneration_)
  {
    visual_->viz1()->state(candidateState, LARGE, GREEN, 0);
    visual_->viz1()->state(candidateState, ROBOT, GREEN, 0);
    visual_->viz1()->trigger();

    if (visualizeGridGenerationWait_)
      visual_->waitForUserFeedback("accepted");
    else
      usleep(0.01 * 1000000);
  }
}

}  // namespace
}  // namespace
}  // namespace
