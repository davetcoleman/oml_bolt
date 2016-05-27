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
   Desc:   Speed up collision checking of SparseVertex edges
*/

#ifndef OMPL_TOOLS_BOLT_EDGE_CACHE_
#define OMPL_TOOLS_BOLT_EDGE_CACHE_

// OMPL
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/tools/bolt/BoltGraph.h>

// Boost
#include <boost/thread/shared_mutex.hpp>

// C++
#include <map>

namespace ompl
{
namespace tools
{
namespace bolt
{
OMPL_CLASS_FORWARD(EdgeCache);
OMPL_CLASS_FORWARD(SparseDB);

typedef std::pair<StateID, StateID> CachedEdge;
typedef std::map<CachedEdge, bool> EdgeCacheMap;

class EdgeCache
{
public:

  /** \brief Constructor */
  EdgeCache(base::SpaceInformationPtr si, SparseDB *sparseDB, VisualizerPtr visual);

  /** \brief Save cache to file */
  bool save();

  /** \brief Load cache from file */
  bool load();

  /** \brief Clear the edge cache completely */
  void clear();

  /** \brief Reset the counters */
  void resetCounters();

  /** \brief Returns true if motion is valid, false if in collision. Checks cache first and also stores result  */
  bool checkMotionWithCache(const StateID &v1, const StateID &v2, const std::size_t &threadID);

  /** \brief Test for benchmarking cache */
  void checkMotionCacheBenchmark();

  /** \brief Ensure that the collision edge data loaded from file is correct, for debuging */
  void errorCheckData();

  void setFilePath(const std::string& filePath);

  std::size_t getCacheSize();

  std::size_t getTotalCollisionChecksFromCache();

  std::size_t getTotalCollisionChecks();

  /** \brief Get stats */
  double getPercentCachedCollisionChecks();

private:

  /** \brief The created space information */
  base::SpaceInformationPtr si_;

  /** \brief The database of motions to search through */
  SparseDB* sparseDB_;

  /** \brief Class for managing various visualization features */
  VisualizerPtr visual_;

  /** \brief Cache previously performed collision checks */
  EdgeCacheMap collisionCheckEdgeCache_;

  /** \brief Stats for the checkMotionWithCache feature */
  std::vector<std::size_t> totalCollisionChecks_;
  std::vector<std::size_t> totalCollisionChecksFromCache_;

  /** \brief Where to store the cache on disk */
  std::string filePath_;

  /** \brief Mutex for both reading or writing to the EdgeCacheMap */
  boost::shared_mutex edgeCacheMutex_;

  std::vector<CachedEdge> keys_;

  /** \brief Force all collision checks to be performed from scratch */
  bool disableCache_ = false;

  /** \brief Available cores */
  std::size_t numThreads_;

};  // end of class EdgeCache

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_EDGE_CACHE_
