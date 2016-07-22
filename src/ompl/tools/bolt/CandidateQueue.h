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

#ifndef OMPL_TOOLS_BOLT_CANDIDATE_QUEUE_H
#define OMPL_TOOLS_BOLT_CANDIDATE_QUEUE_H

// OMPL
#include <ompl/base/SpaceInformation.h>
#include <ompl/tools/bolt/BoostGraphHeaders.h>
#include <ompl/tools/bolt/SamplingQueue.h>

// C++
#include <queue>
#include <thread>

namespace ompl
{
namespace tools
{
namespace bolt
{
OMPL_CLASS_FORWARD(CandidateQueue);
OMPL_CLASS_FORWARD(SparseGraph);
OMPL_CLASS_FORWARD(SparseCriteria);

class CandidateQueue
{
public:
  /** \brief Constructor */
  CandidateQueue(SparseGraphPtr sg);

  ~CandidateQueue();

  void startGenerating(std::size_t indent);

  void stopGenerating(std::size_t indent);

  /** \brief This function is called from the parent thread */
  CandidateData getNextCandidate(std::size_t indent);

  /** \brief This function is called from the parent thread */
  void setCandidateUsed(bool wasUsed, std::size_t indent);

  std::size_t getTotalMisses()
  {
    return totalMisses_;
  }

private:

  void generatingThread(std::size_t threadID, base::SpaceInformationPtr si, std::size_t indent);

  /** \brief Do not add more states if queue is full */
  void waitForQueueNotFull(std::size_t indent);

  /** \brief Wait until there is at least one state ready */
  void waitForQueueNotEmpty(std::size_t indent);

  bool findGraphNeighbors(CandidateData &candidateD, std::size_t threadID, std::size_t indent);

  SparseGraphPtr sg_;
  SparseCriteriaPtr sc_;

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  SamplingQueuePtr samplingQueue_;

  std::queue<CandidateData> queue_;

  std::size_t targetQueueSize_ = 10;

  std::vector<boost::thread *> samplingThreads_;

  /** \brief Mutex for  */
  boost::shared_mutex candidateQueueMutex_;

  /** \brief Used to end neighbor search early */
  bool abortNeighborSearch_ = false;

  /** \brief Flag indicating sampler is active */
  bool threadsRunning_ = false;

  std::size_t numThreads_ = 1;
  std::size_t totalMisses_ = 0;

public:
  bool verbose_ = false;      // general program direction
  bool vNeighbor_ = false;   // nearest neighbor search
  bool vClear_ = false;       // when queue is being cleared because of change
  bool vQueueFull_ = false;  // status of queue
  bool vQueueEmpty_ = false; // alert when queue is empty and holding up process

};  // end of class CandidateQueue

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_CANDIDATE_QUEUE_H
