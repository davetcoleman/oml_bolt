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
#include <ompl/tools/bolt/VertexDiscretizer.h>
#include <ompl/base/State.h>
#include <ompl/geometric/PathGeometric.h>

// Boost
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
//#include <boost/math/constants/constants.hpp>

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
VertexDiscretizer::VertexDiscretizer(base::SpaceInformationPtr si, VisualizerPtr visual)
  : si_(si), visual_(visual)
{
  numThreads_ = boost::thread::hardware_concurrency();

  // Debugging
  if (false)
  {
    OMPL_WARN("Overriding number of threads for testing to 1");
    numThreads_ = 1;
  }
}

bool VertexDiscretizer::generate()
{
  OMPL_INFORM("Generating vertices");

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

  OMPL_INFORM("Generated %i vertices in %f sec (%f hours)", candidateVertices_.size(), vertexDuration,
              vertexDuration / 60.0 / 60.0);

  // Error check
  if (candidateVertices_.size() < 2)
  {
    OMPL_ERROR("No vertices generated, failing");
    exit(-1);
  }

  // Total benchmark runtime
  double totalDuration = time::seconds(time::now() - totalStartTime);

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Vertex Discretization stats:");
  OMPL_INFORM("   Total valid vertices:   %u", candidateVertices_.size());
  OMPL_INFORM("   Vertex generation time: %f seconds (%f min)", vertexDuration, vertexDuration / 60.0);
  OMPL_INFORM("   Total grid gen. time:   %f seconds (%f min)", totalDuration, totalDuration / 60.0);
  OMPL_INFORM("------------------------------------------------------");

  return true;
}

void VertexDiscretizer::generateVertices()
{
  std::size_t dim = si_->getStateSpace()->getDimension();

  // Setup bounds for joint 0
  ob::RealVectorBounds bounds = si_->getStateSpace()->getBounds();
  assert(bounds.high.size() == bounds.low.size());
  assert(bounds.high.size() == dim);

  // Divide joint 0 between threads
  const std::size_t jointID = 0;
  const double range = bounds.high[jointID] - bounds.low[jointID] - startingValueOffset_;
  const std::size_t jointIncrements = ceil(range / discretization_);

  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "Discretization Setup: " << std::endl;
  std::cout << "  Dimensions:             " << dim << std::endl;
  std::cout << "  Discretization:         " << discretization_ << std::endl;
  std::cout << "  J0 Low:                 " << bounds.low[jointID] << std::endl;
  std::cout << "  J0 High:                " << bounds.high[jointID] << std::endl;
  std::cout << "  J0 Range:               " << range << std::endl;
  std::cout << "  J0 Increments:          " << jointIncrements << std::endl;
  assert(jointIncrements >= 1);

  // Check that we have enough jointIncrements for all the threads
  if (jointIncrements < numThreads_)
  {
    OMPL_WARN("There are fewer joint_0 increments (%u) at current discretization (%f) than available threads (%u), "
              "underutilizing threading",
              jointIncrements, discretization_, numThreads_);
    OMPL_INFORM("Optimal discretization: %f", range / double(numThreads_));
    numThreads_ = jointIncrements;
  }
  std::size_t jointIncrementsPerThread = jointIncrements / numThreads_;

  std::cout << "  J0 IncrementsPerThread: " << jointIncrementsPerThread << std::endl;
  std::cout << "  Total states:           " << pow(jointIncrements, dim) << std::endl;
  std::cout << "  Num Threads:            " << numThreads_ << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;

  // Setup threading
  std::vector<boost::thread *> threads(numThreads_);

  double startJointValue = bounds.low[jointID] + startingValueOffset_;
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

    if (verbose_)
      std::cout << "Thread " << i << " has values from " << startJointValue << " to " << endJointValue << std::endl;

    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    threads[i] =
        new boost::thread(boost::bind(&VertexDiscretizer::generateVerticesThread, this, i, startJointValue, endJointValue, si));

    startJointValue = endJointValue;
  }

  // Join threads
  for (std::size_t i = 0; i < threads.size(); ++i)
  {
    threads[i]->join();
    delete threads[i];
  }
}

void VertexDiscretizer::generateVerticesThread(std::size_t threadID, double startJointValue, double endJointValue,
                                     base::SpaceInformationPtr si)
{
  std::size_t jointID = 0;
  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();
  base::State *candidateState = si->getStateSpace()->allocState();

  // Prepare for recursion
  std::vector<double> values(si->getStateSpace()->getDimension(), 0);

  // Customize for different state spaces TODO more generic
  const std::size_t dim = si->getStateSpace()->getDimension();
  std::size_t maxDiscretizationLevel = dim - 1;

  // if (dim == 3)
  // {
  // maxDiscretizationLevel = 1;  // because the third level (numbered '2') is for task space
  // values[2] = 0.0;             // task space
  // }
  if (dim == 6 && false)
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

  // Loop through current joint
  for (double value = startJointValue; value < endJointValue; value += discretization_)
  {
    // User feedback on thread 0
    if (threadID == numThreads_ - 1)
    {
      const double percent = (value - startJointValue) / (endJointValue - startJointValue) * 100.0;
      std::cout << "Vertex generation progress: " << std::setprecision(1) << percent
                << " % Total vertices: " << (candidateVertices_.size() > 0 ? candidateVertices_.size() - 1 : 0) << std::endl;
    }

    // Set first joint value
    values[jointID] = value;

    // Keep recursing
    recursiveDiscretization(threadID, values, jointID + 1, si, candidateState, maxDiscretizationLevel);
  }

  // Cleanup
  si->freeState(candidateState);
}

void VertexDiscretizer::recursiveDiscretization(std::size_t threadID, std::vector<double> &values, std::size_t jointID,
                                          base::SpaceInformationPtr si, base::State *candidateState,
                                          std::size_t maxDiscretizationLevel)
{
  ob::RealVectorBounds bounds = si->getStateSpace()->getBounds();

  if (verbose_)
    std::cout << "jointID " << jointID << " high " << bounds.high[jointID] << " low " << bounds.low[jointID] << " disc " << discretization_ << " values.size() " << values.size() << std::endl;

  // Error check
  BOOST_ASSERT_MSG(jointID < values.size(), "Joint ID too high");
  BOOST_ASSERT_MSG(bounds.high[jointID] - bounds.low[jointID] > discretization_, "Bounds too small");

  // Loop through current joint
  for (double value = bounds.low[jointID] + startingValueOffset_; value <= bounds.high[jointID]; value += discretization_)
  {
    // User feedback on thread 0 for high dimension spaces
    if (jointID == 1 && threadID == 0 && si->getStateSpace()->getDimension() > 3)
    {
      const double percent = (value - bounds.low[jointID]) / (bounds.high[jointID] - bounds.low[jointID]) * 100.0;
      std::cout << "Level 1 vertex generation progress: " << std::setprecision(1) << percent
                << " % Total vertices: " << candidateVertices_.size() - 1 /*ignore query vertex*/ << std::endl;
    }

    // Set value
    values[jointID] = value;

    // Check if we are at the end of the recursion
    if (verbose_)
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

      // Collision check
      double dist;
      if (!si->getStateValidityChecker()->isValid(candidateState, dist))
      {
        // Visualize
        if (visualizeGridGeneration_)
        {
          // Candidate node rejected
          visual_->viz1State(candidateState, tools::SMALL, tools::RED, 0);
          visual_->vizTrigger(threadID+1);
        }

        continue;
      }

      if (dist < clearance_)
      {
        // Visualize
        if (visualizeGridGeneration_)
        {
          // Candidate node rejected
          visual_->viz1State(candidateState, tools::SMALL, tools::RED, 0);
          visual_->vizTrigger(threadID+1);
        }

        continue;
      }

      // Add vertex to graph
      //GuardType type = START;  // TODO(davetcoleman): type START is dummy

      // Allocate state before mutex
      base::State *newState = si->cloneState(candidateState);
      {
        boost::unique_lock<boost::mutex> scoped_lock(vertexMutex_);
        candidateVertices_.push_back(newState);
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

void VertexDiscretizer::displayVertices()
{
  OMPL_INFORM("Displaing vertex discretizer vertices");
  for (std::size_t i = 0; i < candidateVertices_.size(); ++i)
  {
    visual_->viz1State(candidateVertices_[i], tools::LARGE, tools::BLACK, 0);
  }
  visual_->viz1Trigger();
}

}  // namespace
}  // namespace
}  // namespace
