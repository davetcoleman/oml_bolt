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
   Desc:   Maintain queue of potential samples
*/

// OMPL
#include <ompl/tools/bolt/SamplingQueue.h>

// C++
#include <queue>

namespace ompl
{
namespace tools
{
namespace bolt
{

SamplingQueue::SamplingQueue(base::SpaceInformationPtr si, VisualizerPtr visual,
                             base::MinimumClearanceValidStateSamplerPtr clearanceSampler)
  : si_(si), visual_(visual), clearanceSampler_(clearanceSampler)
{
  // statesQueue_.reserve(targetQueueSize_);
}

SamplingQueue::~SamplingQueue()
{
  // Clear all left over states that weren't used
  while (!statesQueue_.empty())
  {
    si_->freeState(statesQueue_.front());
    statesQueue_.pop();
  }
}

void SamplingQueue::startSampling(std::size_t indent)
{
  BOLT_FUNC(indent, true, "startSampling() Starting sampling thread");
  if (threadRunning_)
  {
    BOLT_RED_DEBUG(indent, true, "SamplingQueue already running");
    return;
  }
  threadRunning_ = true;

  // Setup
  base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
  si->setStateValidityChecker(si_->getStateValidityChecker());
  si->setMotionValidator(si_->getMotionValidator());

  // Create thread
  samplingThread_ = new boost::thread(boost::bind(&SamplingQueue::samplingThread, this, si, indent));

  // Wait for first sample to be found
  BOLT_DEBUG(indent, true, "Waiting for first sample to be found");
  while (statesQueue_.empty())
  {
    usleep(0.001 * 1000000);
  }
}

void SamplingQueue::stopSampling(std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "SamplingQueue.stopSampling() Stopping sampling thread");
  threadRunning_ = false;

  // End thread
  samplingThread_->join();
  BOLT_DEBUG(indent, verbose_, "SamplingQueue.stopSampling() Sampling threads have stopped");
}

/** \brief This function is called from the parent thread */
void SamplingQueue::getNextState(base::State *&state, std::size_t indent)
{
  if (statesQueue_.empty())
    waitForQueueNotEmpty(indent + 2);

  while (threadRunning_)
  {
    //boost::lock_guard<boost::shared_mutex> lock(sampleQueueMutex_);
    // Something changed before we got the mutex
    if (statesQueue_.empty())
    {
      usleep(10);  // let the sampler add more states
      // BOLT_YELLOW_DEBUG(indent, true, "SamplingQueue became empty while getting mutex");
      continue;
    }
    BOLT_DEBUG(indent, false, "SamplingQueue size: " << statesQueue_.size());

    state = statesQueue_.front();
    statesQueue_.pop();
    break;
  }
}

void SamplingQueue::samplingThread(base::SpaceInformationPtr si, std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "samplingThread()");

  while (threadRunning_ && !visual_->viz1()->shutdownRequested())
  {
    // Do not add more states if queue is full
    waitForQueueNotFull(indent);

    // Create new state or recycle one
    base::State *candidateState;
    candidateState = si_->allocState();

    // Sample randomly
    if (!clearanceSampler_->sample(candidateState))
    {
      OMPL_ERROR("Unable to find valid sample");
      exit(-1);  // this should never happen
    }

    {
      //boost::lock_guard<boost::shared_mutex> lock(sampleQueueMutex_);
      statesQueue_.push(candidateState);
    }
  }
}

/** \brief Do not add more states if queue is full */
void SamplingQueue::waitForQueueNotFull(std::size_t indent)
{
  bool oneTimeFlag = true;
  while (statesQueue_.size() >= targetQueueSize_ && threadRunning_)
  {
    if (oneTimeFlag)
    {
      BOLT_DEBUG(indent, vQueueFull_, "SamplingQueue: Queue is full, sampler is waiting");
      oneTimeFlag = false;
    }
    usleep(10);
  }
  if (!oneTimeFlag)
    BOLT_DEBUG(indent, vQueueFull_ && false, "SamplingQueue: No longer waiting on full queue");
}

/** \brief Wait until there is at least one state ready */
void SamplingQueue::waitForQueueNotEmpty(std::size_t indent)
{
  bool oneTimeFlag = true;
  while (statesQueue_.empty() && threadRunning_)
  {
    if (oneTimeFlag)
    {
      BOLT_YELLOW_DEBUG(indent, vQueueEmpty_, "SamplingQueue: Queue is empty, waiting to provide a new sampled state");
      oneTimeFlag = false;
    }
    usleep(10);
  }
  if (!oneTimeFlag)
    BOLT_DEBUG(indent, vQueueEmpty_ && false, "SamplingQueue: No longer waiting on queue");
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
