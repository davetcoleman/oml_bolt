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
   Desc:   Save and load from file
*/

#ifndef OMPL_TOOLS_BOLT_BOLT_STORAGE_H_
#define OMPL_TOOLS_BOLT_BOLT_STORAGE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <ompl/base/State.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/util/ClassForward.h>

// Boost
#include <boost/noncopyable.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/access.hpp>
#include <boost/archive/archive_exception.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace ompl
{
namespace tools
{
namespace bolt
{
/**
   @anchor BoltStorage
   @par Short description
*/

/// @cond IGNORE
OMPL_CLASS_FORWARD(BoltStorage);
OMPL_CLASS_FORWARD(SparseGraph);
/// @endcond

/** \class ompl::tools::bolt::BoltStoragePtr
    \brief A boost shared pointer wrapper for ompl::tools::bolt::BoltStorage */

static const boost::uint32_t OMPL_PLANNER_DATA_ARCHIVE_MARKER = 0x5044414D;  // this spells PDAM

class BoltStorage
{
public:
  /// \brief Information stored at the beginning of the BoltStorage archive
  struct Header
  {
    /// \brief OMPL specific marker (fixed value)
    boost::uint32_t marker;

    /// \brief Number of vertices stored in the archive
    std::size_t vertex_count;

    /// \brief Number of edges stored in the archive
    std::size_t edge_count;

    /// \brief Signature of state space that allocated the saved states in the vertices (see
    /// ompl::base::StateSpace::computeSignature()) */
    std::vector<int> signature;

    /// \brief boost::serialization routine
    template <typename Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
      ar &marker;
      ar &vertex_count;
      ar &edge_count;
      ar &signature;
    }
  };

  /// \brief The object containing all vertex data that will be stored
  struct BoltVertexData
  {
    template <typename Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
      ar &stateSerialized_;
      ar &type_;
    }

    std::vector<unsigned char> stateSerialized_;
    int type_;
  };

  /// \brief The object containing all edge data that will be stored
  struct BoltEdgeData
  {
    template <typename Archive>
    void serialize(Archive &ar, const unsigned int /*version*/)
    {
      ar &endpoints_;
      ar &weight_;
      ar &type_;
    }

    std::pair<unsigned int, unsigned int> endpoints_;
    double weight_;
    int type_;
  };

  /** \brief Constructor */
  BoltStorage(const base::SpaceInformationPtr &si, SparseGraph *sparseGraph);

  void save(const std::string &filePath);

  void save(std::ostream &out);

  /// \brief Serialize and save all vertices in \e pd to the binary archive.
  void saveVertices(boost::archive::binary_oarchive &oa);

  /// \brief Serialize and store all edges in \e pd to the binary archive.
  void saveEdges(boost::archive::binary_oarchive &oa);

  bool load(const std::string &filePath);

  bool load(std::istream &in);

  /// \brief Read \e numVertices from the binary input \e ia and store them as BoltStorage
  void loadVertices(unsigned int numVertices, boost::archive::binary_iarchive &ia);

  /// \brief Read \e numEdges from the binary input \e ia and store them as BoltStorage
  void loadEdges(unsigned int numEdges, boost::archive::binary_iarchive &ia);

  /// \brief The space information instance for this data.
  base::SpaceInformationPtr si_;

  SparseGraph *sparseGraph_;

  /** \brief Based on number of threads system is using */
  std::size_t numQueryVertices_;
};  // end of class BoltStorage

}  // namespace bolt
}  // namespace tools
}  // namespace ompl

#endif
