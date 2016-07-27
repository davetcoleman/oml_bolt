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
   Desc:   Populates sparse graph with sparse criteria using both discretization and random sampling
*/

#ifndef OMPL_TOOLS_BOLT_SPARSE_GENERATOR_
#define OMPL_TOOLS_BOLT_SPARSE_GENERATOR_

// OMPL
#include <ompl/tools/bolt/SparseGraph.h>
#include <ompl/tools/bolt/SamplingQueue.h>
#include <ompl/tools/bolt/CandidateQueue.h>

namespace ompl
{
namespace tools
{
namespace bolt
{

/// @cond IGNORE
OMPL_CLASS_FORWARD(SparseGenerator);
/// @endcond

/** \class ompl::tools::bolt::::SparseGeneratorPtr
    \brief A boost shared pointer wrapper for ompl::tools::SparseGenerator */

class SparseGenerator : public std::enable_shared_from_this<SparseGenerator>
{
public:
  /** \brief Constructor needs the state space used for planning.
   */
  SparseGenerator(SparseGraphPtr sg);

  /** \brief Deconstructor */
  virtual ~SparseGenerator();

  /** \brief Give the sparse graph reference to the criteria, because sometimes it needs data from there */
  void setSparseCriteria(SparseCriteriaPtr sparseCriteria)
  {
    sparseCriteria_ = sparseCriteria;
  }

  /** \brief Initialize sparse parameters */
  bool setup(std::size_t indent);

  /** \brief Create a SPARS graph */
  void createSPARS();

  void copyPasteState(std::size_t numSets = 0);

  void addDiscretizedStates(std::size_t indent);

  /** \brief Randomly sample */
  bool addRandomSamples(std::size_t indent);
  bool addRandomSamplesOneThread(std::size_t indent);
  bool addRandomSamplesTwoThread(std::size_t indent);

  /**
   * \brief Add state to sparse graph
   * \param stateID representing a pre-populate state
   * \return true if sparse graph is still accepting states, false if the sparse graph has completed
   */
  bool addSample(CandidateData& candidateD, std::size_t threadID, bool& usedState, std::size_t indent);

  /**
   * \brief Get neighbors within sparseDelta radius
   * \param indent - debugging tool
   */
  void findGraphNeighbors(CandidateData& candidateD, std::size_t threadID, std::size_t indent);

  /** \brief Getter for vertexDiscretizer */
  VertexDiscretizerPtr& getVertexDiscretizer()
  {
    return vertexDiscretizer_;
  }

  std::size_t getNumRandSamplesAdded()
  {
    return numRandSamplesAdded_;
  }

  SamplingQueuePtr getSamplingQueue()
  {
    return samplingQueue_;
  }

protected:
  /** \brief Short name of this class */
  const std::string name_ = "SparseGenerator";

  /** \brief Sparse graph main datastructure that this class operates on */
  SparseGraphPtr sg_;

  SparseCriteriaPtr sparseCriteria_;

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  /** \brief Sampler user for generating valid samples in the state space */
  base::MinimumClearanceValidStateSamplerPtr clearanceSampler_;

  /** \brief Secondary thread for sampling and garbage collection */
  SamplingQueuePtr samplingQueue_;

  /** \brief Multiple threads for finding nearest neighbors from samples */
  CandidateQueuePtr candidateQueue_;

  std::size_t numConsecutiveFailures_;
  std::size_t maxConsecutiveFailures_ = 0; // find the closest to completion the process has gotten
  std::size_t maxPercentComplete_; // the whole number percentage presented to user

  VertexDiscretizerPtr vertexDiscretizer_;

  /** \brief For statistics */
  std::size_t numRandSamplesAdded_ = 0;
  time::point timeRandSamplesStarted_; // calculate rate at which the graph is being built
  time::point timeDiscretizeAndRandomStarted_;
public:
  /** \brief Number of failed state insertion attempts before stopping the algorithm */
  std::size_t terminateAfterFailures_ = 1000;

  /** \brief Number of failed state insertion attempts before starting to apply the fourth quality criteria from SPARS
   */
  std::size_t fourthCriteriaAfterFailures_ = 500;

  /** \brief How often to save */
  std::size_t saveInterval_ = 1000;

  /** \brief Generate the Sparse graph with discretized and/or random samples */
  bool useDiscretizedSamples_;
  bool useRandomSamples_;

};  // end SparseGenerator

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_SPARSE_GENERATOR_
