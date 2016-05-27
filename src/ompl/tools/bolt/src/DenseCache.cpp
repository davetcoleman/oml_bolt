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
   Desc:   Speed up collision checking of StateID edges
*/

// OMPL
#include <ompl/tools/bolt/SparseDB.h>

// Boost
#include <boost/serialization/map.hpp>
// #include <boost/archive/text_iarchive.hpp>
// #include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/thread.hpp>

// C++
#include <fstream>

namespace og = ompl::geometric;
namespace ot = ompl::tools;
namespace otb = ompl::tools::bolt;
namespace ob = ompl::base;

namespace ompl
{
namespace tools
{
namespace bolt
{
DenseCache::DenseCache(base::SpaceInformationPtr si, SparseDB *sparseDB, VisualizerPtr visual)
  : si_(si), sparseDB_(sparseDB), visual_(visual)
{
  resetCounters();

  numThreads_ = boost::thread::hardware_concurrency();
  keys_.resize(numThreads_);
  totalCollisionChecks_.resize(numThreads_, 0);
  totalCollisionChecksFromCache_.resize(numThreads_, 0);

  if (disableCache_)
    OMPL_WARN("Edge cache disabled");
}

bool DenseCache::save()
{
  OMPL_INFORM("------------------------------------------------");
  OMPL_INFORM("Saving Edge Cache");
  OMPL_INFORM("  Path:         %s", filePath_.c_str());
  OMPL_INFORM("  Utilization:  %f", getPercentCachedCollisionChecks());
  OMPL_INFORM("  Cache Size:   %u", getCacheSize());
  OMPL_INFORM("  Theory Max:   %u", sparseDB_->getNumVertices() * sparseDB_->getNumVertices());
  OMPL_INFORM("------------------------------------------------");

  std::ofstream out(filePath_, std::ios::binary);

  if (!out.good())
  {
    OMPL_ERROR("Failed to save edge cache: output stream is invalid");
    return false;
  }

  {  // Get write mutex
    boost::lock_guard<boost::shared_mutex> writeLock(denseCacheMutex_);
    try
    {
      boost::archive::binary_oarchive oa(out);
      oa << collisionCheckDenseCache_;
    }
    catch (boost::archive::archive_exception &ae)
    {
      OMPL_ERROR("Failed to save edge cache: %s", ae.what());
      return false;
    }
  }

  out.close();
  return true;
}

bool DenseCache::load()
{
  OMPL_INFORM("Loading edge cache from %s", filePath_.c_str());

  // Benchmark runtime
  time::point startTime = time::now();

  if (!collisionCheckDenseCache_.empty())
  {
    OMPL_ERROR("Collision check edge cache has %u edges and is not empty, unable to load.", collisionCheckDenseCache_.size());
    return false;
  }

  std::ifstream in(filePath_, std::ios::binary);

  // Error check
  if (!in.good())
  {
    OMPL_INFORM("Failed to load edge cache: input stream is invalid");
    return false;
  }

  try
  {
    boost::archive::binary_iarchive ia(in);
    ia >> collisionCheckDenseCache_;
  }
  catch (boost::archive::archive_exception &ae)
  {
    OMPL_ERROR("Failed to load edge cache: %s", ae.what());
    return false;
  }

  in.close();

  // This is really slow
  //errorCheckData();

  // Benchmark runtime
  double duration = time::seconds(time::now() - startTime);

  OMPL_INFORM("------------------------------------------------------");
  OMPL_INFORM("Loaded edge cache stats:");
  OMPL_INFORM("  Path:         %s", filePath_.c_str());
  OMPL_INFORM("  Cache Size:   %u", getCacheSize());
  OMPL_INFORM("  Theory Max:   %u", sparseDB_->getNumVertices() * sparseDB_->getNumVertices());
  OMPL_INFORM("  Loading time: %f", duration);
  OMPL_INFORM("------------------------------------------------------");

  return true;
}

void DenseCache::clear()
{
  OMPL_INFORM("Clearing edge cache");
  collisionCheckDenseCache_.clear();
  resetCounters();
}

void DenseCache::resetCounters()
{
  for (std::size_t i = 0; i < totalCollisionChecksFromCache_.size(); ++i)
  {
    totalCollisionChecks_[i] = 0;
    totalCollisionChecksFromCache_[i] = 0;
  }
}

bool DenseCache::checkMotionWithCache(const StateID &stateID1, const StateID &stateID2, const std::size_t &threadID)
{
  // Optionally skip caching
  if (disableCache_)
    return si_->checkMotion(sparseDB_->getState(stateID1), sparseDB_->getState(stateID2));

  CachedEdge &edge = keys_[threadID];

  // Error check
  BOOST_ASSERT_MSG(stateID1 >= numThreads_, "stateID1: The queryVertex_ should not be checked within the DenseCache, because it is subject to change");
  BOOST_ASSERT_MSG(stateID2 >= numThreads_, "stateID2: The queryVertex_ should not be checked within the DenseCache, because it is subject to change");

  // Create edge to search for - only store pairs in one direction
  if (stateID1 < stateID2)
    edge = CachedEdge(stateID1, stateID2);
  else
    edge = CachedEdge(stateID2, stateID1);

  DenseCacheMap::iterator lb = collisionCheckDenseCache_.lower_bound(edge);

  bool result;
  {  // Get read-only mutex
    boost::shared_lock<boost::shared_mutex> readLock(denseCacheMutex_);

    // Search for existing key
    result = (lb != collisionCheckDenseCache_.end() && !(collisionCheckDenseCache_.key_comp()(edge, lb->first)));

    // Statistics
    totalCollisionChecks_[threadID]++;
  }

  if (result) // Cache available
  {
    totalCollisionChecksFromCache_[threadID]++;
    // std::cout << "Cache: Edge " << stateID1 << ", " << stateID2 << " collision: " << lb->second << std::endl;
    return lb->second;
  }

  // No cache available
  result = si_->checkMotion(sparseDB_->getState(stateID1), sparseDB_->getState(stateID2));
  // std::cout << "No cache: Edge " << stateID1 << ", " << stateID2 << " collision: " << result << std::endl;

  {  // Get write mutex
    boost::lock_guard<boost::shared_mutex> writeLock(denseCacheMutex_);

    // The key does not exist in the map, so add it to the map
    // Use lb as a hint to insert, so it can avoid another lookup
    collisionCheckDenseCache_.insert(lb, DenseCacheMap::value_type(edge, result));
  }

  return result;
}

void DenseCache::checkMotionCacheBenchmark()
{
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  OMPL_WARN("Running check motion cache benchmark");
  std::cout << std::endl;

  // Test new checkMotionWithCache
  StateID stateID1 = 2;
  StateID stateID2 = 3;
  std::size_t numRuns = 100000;
  bool baselineResult = si_->checkMotion(sparseDB_->getState(stateID1), sparseDB_->getState(stateID2));

  // Benchmark collision check without cache
  {
    // Benchmark runtime
    time::point startTime = time::now();

    for (std::size_t i = 0; i < numRuns; ++i)
    {
      if (si_->checkMotion(sparseDB_->getState(stateID1), sparseDB_->getState(stateID2)) != baselineResult)
      {
        OMPL_ERROR("benchmark checkmotion does not match baseline result");
        exit(-1);
      }
    }

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);
    OMPL_INFORM(" - no cache took %f seconds (%f hz)", duration, 1.0 / duration);
  }

  // Benchmark collision check with cache
  {
    // Benchmark runtime
    time::point startTime = time::now();

    const std::size_t threadID = 0;
    for (std::size_t i = 0; i < numRuns / 2; ++i)
    {
      if (checkMotionWithCache(stateID1, stateID2, threadID) != baselineResult)
      {
        OMPL_ERROR("benchmark checkmotion does not match baseline result");
        exit(-1);
      }

      if (checkMotionWithCache(stateID2, stateID1, threadID) != baselineResult)
      {
        OMPL_ERROR("benchmark checkmotion does not match baseline result");
        exit(-1);
      }
    }
    OMPL_INFORM("Cache size: %u", collisionCheckDenseCache_.size());

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);
    OMPL_INFORM(" - with cache took %f seconds (%f hz)", duration, 1.0 / duration);
  }
}

void DenseCache::errorCheckData()
{
  OMPL_INFORM("Error checking edge cache...");
  std::size_t counter = 0;
  for (DenseCacheMap::const_iterator iterator = collisionCheckDenseCache_.begin(); iterator != collisionCheckDenseCache_.end();
       iterator++)
  {
    std::pair<StateID,StateID> thing = iterator->first;
    StateID &stateID1 = thing.first;
    StateID &stateID2 = thing.second;
    bool cachedResult = iterator->second;
    if (si_->checkMotion(sparseDB_->getState(stateID1), sparseDB_->getState(stateID2)) != cachedResult)
    {
      OMPL_ERROR("Found instance where cached edge data is wrong, on iteration %u", counter);
      std::cout << "stateID1: " << stateID1 << std::endl;
      std::cout << "stateID2: " << stateID2 << std::endl;
      exit(-1);
    }
    if (counter++ % 1000 == 0)
      std::cout << "Checking edge " << counter << std::endl;
  }
}

void DenseCache::setFilePath(const std::string &filePath)
{
  filePath_ = filePath;
}

std::size_t DenseCache::getCacheSize()
{
  return collisionCheckDenseCache_.size();
}

std::size_t DenseCache::getTotalCollisionChecksFromCache()
{
  std::size_t total = 0;
  for (std::size_t i = 0; i < totalCollisionChecksFromCache_.size(); ++i)
  {
    total += totalCollisionChecksFromCache_[i];
  }
  return total;
}

std::size_t DenseCache::getTotalCollisionChecks()
{
  std::size_t total = 0;
  for (std::size_t i = 0; i < totalCollisionChecks_.size(); ++i)
  {
    total += totalCollisionChecks_[i];
  }
  return total;
}

double DenseCache::getPercentCachedCollisionChecks()
{
  std::size_t totalCollisionChecks = getTotalCollisionChecks();
  if (totalCollisionChecks == 0)
    return 0;

  return static_cast<double>(getTotalCollisionChecksFromCache()) / totalCollisionChecks * 100.0;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
