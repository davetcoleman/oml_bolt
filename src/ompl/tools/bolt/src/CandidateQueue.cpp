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
   Desc:   Maintain queue of potential smaples
*/

// OMPL
#include <ompl/tools/bolt/CandidateQueue.h>
#include <ompl/tools/bolt/SparseCriteria.h> // must be included only in cpp

// C++
#include <queue>
#include <thread>

namespace ompl
{
namespace tools
{
namespace bolt
{
CandidateQueue::CandidateQueue(SparseGraphPtr sg)
  : sg_(sg)
  , sc_(sg_->getSparseCriteria())
  , si_(sg_->getSpaceInformation())
  , visual_(sg_->getVisual())
  , samplingQueue_(sc_->getSamplingQueue())
{
}

CandidateQueue::~CandidateQueue()
{
  while (!queue_.empty())
  {
    si_->freeState(queue_.front().state_);
    queue_.pop();
  }
}

void CandidateQueue::startGenerating(std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "startGenerating() Starting candidate queue thread");
  if (threadsRunning_)
  {
    BOLT_RED_DEBUG(indent, true, "CandidateQueue already running");
    return;
  }
  threadsRunning_ = true;

  // Set number threads - should be at least less than 1 from total number of threads on system
  // 1 thread is for parent, 1 is for sampler, 1 is for GUIs, etc, remainder are for this
  numThreads_ = std::max(1, int(sg_->getNumQueryVertices() - 3));
  BOLT_DEBUG(indent, true, "Running CandidateQueue with " << numThreads_ << " threads");
  if (numThreads_ < 2)
    BOLT_YELLOW_DEBUG(indent, true, "Only running CandidateQueue with 1 thread");
  if (numThreads_ >= sg_->getNumQueryVertices())
  {
    BOLT_RED_DEBUG(indent, true, "Too many threads requested for candidate queue");
    exit(-1);
  }

  // Create threads
  samplingThreads_.resize(numThreads_);

  // Starts on thread 1 because thread 0 is reserved for parent process
  std::size_t skipNCountThreads = 1;
  for (std::size_t i = skipNCountThreads; i < samplingThreads_.size() + skipNCountThreads; ++i)
  {
    // Create new collision checker
    base::SpaceInformationPtr si(new base::SpaceInformation(si_->getStateSpace()));
    si->setStateValidityChecker(si_->getStateValidityChecker());
    si->setMotionValidator(si_->getMotionValidator());

    samplingThreads_[i] = new boost::thread(boost::bind(&CandidateQueue::generatingThread, this, i, si, indent));
  }

  // Wait for first sample to be found
  // BOLT_DEBUG(indent, true, "Waiting for first candidate to be found");
  while (queue_.empty())
  {
    usleep(0.001 * 1000000);
  }
}

void CandidateQueue::stopGenerating(std::size_t indent)
{
  BOLT_FUNC(indent, true, "stopGenerating() Stopping generating thread");
  threadsRunning_ = false;
  abortNeighborSearch_ = true;

  // Join threads
  for (std::size_t i = 0; i < samplingThreads_.size(); ++i)
  {
    samplingThreads_[i]->join();
    delete samplingThreads_[i];
  }
}

void CandidateQueue::generatingThread(std::size_t threadID, base::SpaceInformationPtr si, std::size_t indent)
{
  BOLT_FUNC(indent, verbose_, "generatingThread() " << threadID);

  base::State *candidateState;

  while (threadsRunning_ && !visual_->viz1()->shutdownRequested())
  {
    BOLT_DEBUG(indent + 2, false, threadID << "generatingThread: Running while loop on thread " << threadID);

    // Do not add more states if queue is full
    if (queue_.size() > targetQueueSize_)
      waitForQueueNotFull(indent + 2);

    // Get next sample
    samplingQueue_->getNextState(candidateState, indent + 2);

    BOLT_DEBUG(indent + 2, false, "New candidateState: " << candidateState << " on thread " << threadID);

    // Find nearby nodes
    CandidateData candidateD(candidateState);
    bool result;
    {
      // Get shared read access - we don't want to delete the queue if this is running
      // boost::shared_lock<boost::shared_mutex> lock(candidateQueueMutex_);
      boost::lock_guard<boost::shared_mutex> lock(candidateQueueMutex_);
      result = findGraphNeighbors(candidateD, threadID, indent + 2);
    }

    if (!result)
    {
      // Wait until abort is over before finding next candidate
      while (abortNeighborSearch_)
      {
        usleep(0.0001 * 1000000);
      }
      continue;
    }

    // Add to queue - thread-safe
    {
      boost::lock_guard<boost::shared_mutex> lock(candidateQueueMutex_);
      queue_.push(candidateD);
    }
  }
}

/** \brief This function is called from the parent thread */
CandidateData CandidateQueue::getNextCandidate(std::size_t indent)
{
  // if (!queue_.empty())
  //   BOLT_CYAN_DEBUG(indent, true, "getNextCanidate(): not empty queue!");

  waitForQueueNotEmpty(indent + 2);

  return queue_.front();
}

/** \brief This function is called from the parent thread */
void CandidateQueue::setCandidateUsed(bool wasUsed, std::size_t indent)
{
  if (wasUsed)  // this means whole queue is invalid :-/
  {
    BOLT_RED_DEBUG(indent, vClear_, "Clearing candidate queue of size " << queue_.size());

    // Abort current neighbor search early
    abortNeighborSearch_ = true;

    // Get mutex on candidate queue
    boost::lock_guard<boost::shared_mutex> lock(candidateQueueMutex_);

    // Free first state without freeing memory, because it is in use. I'm not sure how though...
    queue_.pop();

    // Free states and clear queue
    while (!queue_.empty())
    {
      // samplingQueue_->recycleState(queue_.front().state_);
      si_->freeState(queue_.front().state_);
      queue_.pop();
    }
    abortNeighborSearch_ = false;  // allow search to continue

    return;
  }

  // Instead of freeing the unused state, recycle the memory
  // samplingQueue_->recycleState(queue_.front().state_);
  si_->freeState(queue_.front().state_);
  queue_.pop();
}

/** \brief Do not add more states if queue is full */
void CandidateQueue::waitForQueueNotFull(std::size_t indent)
{
  bool oneTimeFlag = true;
  while (queue_.size() >= targetQueueSize_ && threadsRunning_)
  {
    if (oneTimeFlag)
    {
      BOLT_DEBUG(indent, vQueueFull_, "CandidateQueue: Queue is full, generator is waiting");
      oneTimeFlag = false;
    }
    usleep(0.001 * 1000000);
  }
  if (!oneTimeFlag)
    BOLT_DEBUG(indent, vQueueFull_, "CandidateQueue: No longer waiting on full queue");
}

/** \brief Wait until there is at least one state ready */
void CandidateQueue::waitForQueueNotEmpty(std::size_t indent)
{
  bool oneTimeFlag = true;
  while (queue_.empty() && threadsRunning_)
  {
    if (oneTimeFlag)
    {
      BOLT_YELLOW_DEBUG(indent, vQueueEmpty_, "CandidateQueue: Queue is empty, waiting for next generated "
                                              "CandidateData");
      oneTimeFlag = false;
    }
    usleep(0.001 * 1000000);
  }
  if (!oneTimeFlag)
    BOLT_DEBUG(indent, vQueueEmpty_ && false, "CandidateQueue: No longer waiting on queue");
}

bool CandidateQueue::findGraphNeighbors(CandidateData &candidateD, std::size_t threadID, std::size_t indent)
{
  BOLT_FUNC(indent, vNeighbor_, "findGraphNeighbors() within sparse delta " << sc_->getSparseDelta());

  // Search in thread-safe manner
  // Note that the main thread could be modifying the NN, so we have to lock it
  const bool useMutex = true;
  sg_->getQueryStateNonConst(threadID) = candidateD.state_;
  {
    std::lock_guard<std::mutex> lock(sg_->getNNGuard());
    sg_->getNN(useMutex)
        ->nearestR(sg_->getQueryVertices(threadID), sc_->getSparseDelta(), candidateD.graphNeighborhood_);
  }
  sg_->getQueryStateNonConst(threadID) = nullptr;

  // Now that we got the neighbors from the NN, we must remove any we can't see
  for (std::size_t i = 0; i < candidateD.graphNeighborhood_.size(); ++i)
  {
    SparseVertex v2 = candidateD.graphNeighborhood_[i];

    // Check for termination condition
    if (abortNeighborSearch_)
    {
      BOLT_YELLOW_DEBUG(indent, false, "findGraphNeighbors aborted b/c term cond");
      return false;
    }

    // Don't collision check if they are the same state
    if (candidateD.state_ != sg_->getState(v2))
    {
      if (!si_->checkMotion(candidateD.state_, sg_->getState(v2)))
      {
        continue;
      }
    }

    // Check for termination condition
    if (abortNeighborSearch_)
    {
      BOLT_YELLOW_DEBUG(indent, false, "findGraphNeighbors aborted b/c term cond");
      return false;
    }

    // The two are visible to each other!
    candidateD.visibleNeighborhood_.push_back(candidateD.graphNeighborhood_[i]);
  }

  BOLT_DEBUG(indent, vNeighbor_,
             "Graph neighborhood: " << candidateD.graphNeighborhood_.size()
                                    << " | Visible neighborhood: " << candidateD.visibleNeighborhood_.size());

  // Do not return true if abort has been called
  if (abortNeighborSearch_)
  {
    BOLT_YELLOW_DEBUG(indent, false, "findGraphNeighbors aborted b/c term cond");
    return false;
  }
  return true;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
