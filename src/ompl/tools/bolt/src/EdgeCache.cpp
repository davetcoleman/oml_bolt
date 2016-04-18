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
   Desc:   Speed up collision checking of DenseVertex edges
*/

// OMPL
#include <ompl/tools/bolt/DenseDB.h>

// Boost
#include <boost/serialization/map.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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
EdgeCache::EdgeCache(base::SpaceInformationPtr si, DenseDB *denseDB, base::VisualizerPtr visual)
  : si_(si), denseDB_(denseDB), visual_(visual)
{
  resetCounters();
}

bool EdgeCache::save()
{
  OMPL_INFORM("------------------------------------------------");
  OMPL_INFORM("Saving Collision Cache");
  OMPL_INFORM("  Path:         %s", filePath_.c_str());
  OMPL_INFORM("  Utilization:  %f", getPercentCachedCollisionChecks());
  OMPL_INFORM("  Cache Size:   %u", getCacheSize());
  OMPL_INFORM("  Theory Max:   %u", denseDB_->getNumVertices() * denseDB_->getNumVertices());
  OMPL_INFORM("------------------------------------------------");

  std::ofstream out(filePath_, std::ios::binary);

  if (!out.good())
  {
    OMPL_ERROR("Failed to save collision cache: output stream is invalid");
    return false;
  }

  try
  {
    boost::archive::binary_oarchive oa(out);
    oa << collisionCheckEdgeCache_;
  }
  catch (boost::archive::archive_exception &ae)
  {
    OMPL_ERROR("Failed to save collision cache: %s", ae.what());
    return false;
  }

  out.close();
  return true;
}

bool EdgeCache::load()
{
  OMPL_INFORM("Loading collision path from %s", filePath_.c_str());

  std::ifstream in(filePath_, std::ios::binary);

  if (!in.good())
  {
    OMPL_INFORM("Failed to load collision cache: input stream is invalid");
    return false;
  }

  try
  {
    boost::archive::binary_iarchive ia(in);
    ia >> collisionCheckEdgeCache_;
  }
  catch (boost::archive::archive_exception &ae)
  {
    OMPL_ERROR("Failed to load collision cache: %s", ae.what());
    return false;
  }

  in.close();

  // This is really slow
  //errorCheckData();

  return true;
}

void EdgeCache::resetCounters()
{
  totalCollisionChecks_ = 0;
  totalCollisionChecksFromCache_ = 0;
}

bool EdgeCache::checkMotionWithCache(const DenseVertex &v1, const DenseVertex &v2)
{
  // Statistics
  totalCollisionChecks_++;

  // Only store pairs in one direction
  if (v1 < v2)
    key_ = std::pair<DenseVertex, DenseVertex>(v1, v2);
  else
    key_ = std::pair<DenseVertex, DenseVertex>(v2, v1);

  EdgeCacheMap::iterator lb = collisionCheckEdgeCache_.lower_bound(key_);

  if(lb != collisionCheckEdgeCache_.end() && !(collisionCheckEdgeCache_.key_comp()(key_, lb->first)))
  {
    // Cache available
    totalCollisionChecksFromCache_++;
    // std::cout << "Cache: Edge " << v1 << ", " << v2 << " collision: " << lb->second << std::endl;
    return lb->second;
  }
  else
  {
    // no cache available
    bool result = si_->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]);
    // std::cout << "No cache: Edge " << v1 << ", " << v2 << " collision: " << result << std::endl;

    // the key does not exist in the map, so add it to the map
    // Use lb as a hint to insert, so it can avoid another lookup
    collisionCheckEdgeCache_.insert(lb, EdgeCacheMap::value_type(key_, result));

    return result;
  }

}

bool EdgeCache::checkMotionWithCacheSlow(const DenseVertex &v1, const DenseVertex &v2)
{
  // Statistics
  totalCollisionChecks_++;

  // Only store pairs in one direction
  std::pair<DenseVertex, DenseVertex> key;
  if (v1 < v2)
    key = std::pair<DenseVertex, DenseVertex>(v1, v2);
  else
    key = std::pair<DenseVertex, DenseVertex>(v2, v1);

  EdgeCacheMap::const_iterator it = collisionCheckEdgeCache_.find(key);

  if (it == collisionCheckEdgeCache_.end())  // no cache available
  {
    bool result = si_->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]);
    collisionCheckEdgeCache_[key] = result;
    // std::cout << "No cache: Edge " << v1 << ", " << v2 << " collision: " << result << std::endl;
    return result;
  }

  // Cache available
  totalCollisionChecksFromCache_++;
  // std::cout << "Cache: Edge " << v1 << ", " << v2 << " collision: " << it->second << std::endl;
  return it->second;
}

void EdgeCache::checkMotionCacheBenchmark()
{
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  OMPL_WARN("Running check motion cache benchmark");
  std::cout << std::endl;

  // Test new checkMotionWithCache
  DenseVertex v1 = 2;
  DenseVertex v2 = 3;
  std::size_t numRuns = 100000;
  bool baselineResult = si_->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]);

  // Benchmark collision check without cache
  {
    // Benchmark runtime
    time::point startTime = time::now();

    for (std::size_t i = 0; i < numRuns; ++i)
    {
      if (si_->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]) != baselineResult)
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

    for (std::size_t i = 0; i < numRuns / 2; ++i)
    {
      if (checkMotionWithCache(v1, v2) != baselineResult)
      {
        OMPL_ERROR("benchmark checkmotion does not match baseline result");
        exit(-1);
      }

      if (checkMotionWithCache(v2, v1) != baselineResult)
      {
        OMPL_ERROR("benchmark checkmotion does not match baseline result");
        exit(-1);
      }
    }
    OMPL_INFORM("Cache size: %u", collisionCheckEdgeCache_.size());

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);
    OMPL_INFORM(" - with cache took %f seconds (%f hz)", duration, 1.0 / duration);
  }
}

void EdgeCache::errorCheckData()
{
  OMPL_INFORM("Error checking collision cache...");
  std::size_t counter = 0;
  for (EdgeCacheMap::const_iterator iterator = collisionCheckEdgeCache_.begin(); iterator != collisionCheckEdgeCache_.end();
       iterator++)
  {
    std::pair<DenseVertex,DenseVertex> thing = iterator->first;
    DenseVertex &v1 = thing.first;
    DenseVertex &v2 = thing.second;
    bool cachedResult = iterator->second;
    if (si_->checkMotion(denseDB_->stateProperty_[v1], denseDB_->stateProperty_[v2]) != cachedResult)
    {
      OMPL_ERROR("Found instance where cached collision data is wrong, on iteration %u", counter);
      std::cout << "v1: " << v1 << std::endl;
      std::cout << "v2: " << v2 << std::endl;
      exit(-1);
    }
    if (counter++ % 1000 == 0)
      std::cout << "Checking edge " << counter << std::endl;
  }
}

void EdgeCache::setFilePath(const std::string &filePath)
{
  filePath_ = filePath;
}

std::size_t EdgeCache::getCacheSize()
{
  return collisionCheckEdgeCache_.size();
}

std::size_t EdgeCache::getTotalCollisionChecksFromCache()
{
  return totalCollisionChecksFromCache_;
}

std::size_t EdgeCache::getTotalCollisionChecks()
{
  return totalCollisionChecks_;
}

double EdgeCache::getPercentCachedCollisionChecks()
{
  if (totalCollisionChecks_ == 0)
    return 0;

  return static_cast<double>(totalCollisionChecksFromCache_) / totalCollisionChecks_ * 100.0;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl