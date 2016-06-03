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

// OMPL
#include <ompl/tools/bolt/BoltStorage.h>
#include <ompl/tools/bolt/SparseGraph.h>
#include <ompl/tools/bolt/BoostGraphHeaders.h>

// Boost
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#define foreach BOOST_FOREACH

namespace ompl
{
namespace tools
{
namespace bolt
{
BoltStorage::BoltStorage(const base::SpaceInformationPtr &si, SparseGraph *sparseGraph) : si_(si), sparseGraph_(sparseGraph)
{
  numQueryVertices_ = boost::thread::hardware_concurrency();
}

void BoltStorage::save(const std::string& filePath)
{
  std::ofstream out(filePath.c_str(), std::ios::binary);

  OMPL_INFORM("------------------------------------------------");
  OMPL_INFORM("Saving Sparse Graph");
  OMPL_INFORM("  Path:            %s", filePath.c_str());
  OMPL_INFORM("  Edges:           %u", sparseGraph_->getNumEdges());
  OMPL_INFORM("  Vertices:        %u", sparseGraph_->getNumVertices());
  OMPL_INFORM("------------------------------------------------");

  save(out);
  out.close();
}

void BoltStorage::save(std::ostream &out)
{
  if (!out.good())
  {
    OMPL_ERROR("Failed to save BoltData: output stream is invalid");
    return;
  }

  try
  {
    boost::archive::binary_oarchive oa(out);

    // Writing the header
    Header h;
    h.marker = OMPL_PLANNER_DATA_ARCHIVE_MARKER;
    h.vertex_count = sparseGraph_->getNumVertices() - numQueryVertices_;
    h.edge_count = sparseGraph_->getNumEdges();
    si_->getStateSpace()->computeSignature(h.signature);
    oa << h;

    saveVertices(oa);
    saveEdges(oa);
  }
  catch (boost::archive::archive_exception &ae)
  {
    OMPL_ERROR("Failed to save BoltData: %s", ae.what());
  }
}

void BoltStorage::saveVertices(boost::archive::binary_oarchive &oa)
{
  const base::StateSpacePtr &space = si_->getStateSpace();

  std::vector<unsigned char> state(space->getSerializationLength());
  std::size_t feedbackFrequency = sparseGraph_->getNumVertices() / 10;

  std::cout << "         Saving vertices: " << std::flush;
  std::size_t count = 0;
  std::size_t errorCheckNumQueryVertices = 0;
  foreach (const SparseVertex v, boost::vertices(sparseGraph_->g_))
  {
    // Skip the query vertex that is nullptr
    if (v <= sparseGraph_->queryVertices_.back())
    {
      errorCheckNumQueryVertices++;
      continue;
    }

    // Debug
    //std::cout << "v: " << v << " stateID: " << sparseGraph_->getStateID(v) << " state: ";
    //si_->printState(sparseGraph_->getVertexStateNonConst(v), std::cout);

    // Convert to new structure
    BoltVertexData vertexData;

    // Record the type of the vertex
    vertexData.type_ = sparseGraph_->vertexTypeProperty_[v];

    // Serializing the state contained in this vertex
    space->serialize(&state[0], sparseGraph_->getVertexStateNonConst(v));
    vertexData.stateSerialized_ = state;

    // Save to file
    oa << vertexData;

    // Feedback
    if ((++count) % feedbackFrequency == 0)
      std::cout << std::fixed << std::setprecision(0) << (count / double(sparseGraph_->getNumVertices())) * 100.0 << "% "
                << std::flush;
  }
  BOOST_ASSERT_MSG(errorCheckNumQueryVertices == numQueryVertices_,
                   "There should be the same number of query vertex as threads that were skipped while saving");


  std::cout << std::endl;
}

void BoltStorage::saveEdges(boost::archive::binary_oarchive &oa)
{
  std::size_t feedbackFrequency = std::max(10.0, sparseGraph_->getNumEdges() / 10.0);

  std::cout << "         Saving edges: " << std::flush;
  std::size_t count = 1;
  foreach (const SparseEdge e, boost::edges(sparseGraph_->g_))
  {
    const SparseVertex v1 = boost::source(e, sparseGraph_->g_);
    const SparseVertex v2 = boost::target(e, sparseGraph_->g_);

    // Convert to new structure
    BoltEdgeData edgeData;
    // Note: we increment all vertex indexes by the number of queryVertices_ because [0,11] is in use
    edgeData.endpoints_.first = v1 - numQueryVertices_;
    edgeData.endpoints_.second = v2 - numQueryVertices_;

    // Other properties
    edgeData.weight_ = sparseGraph_->edgeWeightProperty_[e];
    edgeData.type_ = sparseGraph_->edgeTypeProperty_[e];

    // Copy to file
    oa << edgeData;

    // Feedback
    if ((++count) % feedbackFrequency == 0)
      std::cout << std::fixed << std::setprecision(0) << (count / double(sparseGraph_->getNumEdges())) * 100.0 << "% "
                << std::flush;

  }  // for each edge
  std::cout << std::endl;
}

bool BoltStorage::load(const std::string& filePath)
{
  OMPL_INFORM("------------------------------------------------");
  OMPL_INFORM("BoltStorage: Loading Sparse Graph");
  OMPL_INFORM("Path: %s", filePath.c_str());

  // Error checking
  if (sparseGraph_->getNumEdges() > numQueryVertices_ ||
      sparseGraph_->getNumVertices() > numQueryVertices_)  // the search verticie may already be there
  {
    OMPL_INFORM("Database is not empty, unable to load from file");
    return false;
  }
  if (filePath.empty())
  {
    OMPL_ERROR("Empty filename passed to save function");
    return false;
  }
  if (!boost::filesystem::exists(filePath))
  {
    OMPL_INFORM("Database file does not exist: %s.", filePath.c_str());
    return false;
  }

  std::ifstream in(filePath.c_str(), std::ios::binary);
  bool result = load(in);
  in.close();

  return result;
}

bool BoltStorage::load(std::istream &in)
{
  if (!in.good())
  {
    OMPL_ERROR("Failed to load BoltData: input stream is invalid");
    return false;
  }

  try
  {
    boost::archive::binary_iarchive ia(in);

    // Read the header
    Header h;
    ia >> h;

    // Checking the archive marker
    if (h.marker != OMPL_PLANNER_DATA_ARCHIVE_MARKER)
    {
      OMPL_ERROR("Failed to load BoltData: BoltData archive marker not found");
      return false;
    }

    // Verify that the state space is the same
    std::vector<int> sig;
    si_->getStateSpace()->computeSignature(sig);
    if (h.signature != sig)
    {
      OMPL_ERROR("Failed to load BoltData: StateSpace signature mismatch");
      return false;
    }

    loadVertices(h.vertex_count, ia);
    loadEdges(h.edge_count, ia);
  }
  catch (boost::archive::archive_exception &ae)
  {
    OMPL_ERROR("Failed to load BoltData: %s", ae.what());
  }

  return true;
}

void BoltStorage::loadVertices(unsigned int numVertices, boost::archive::binary_iarchive &ia)
{
  OMPL_INFORM("Loading %u vertices from file", numVertices);

  const base::StateSpacePtr &space = si_->getStateSpace();
  std::size_t feedbackFrequency = numVertices / 10;

  std::cout << "         Vertices loaded: ";
  for (unsigned int i = 0; i < numVertices; ++i)
  {
    // Copy in data from file
    BoltVertexData vertexData;
    ia >> vertexData;

    // Allocating a new state and deserializing it from the buffer
    base::State *state = space->allocState();
    space->deserialize(state, &vertexData.stateSerialized_[0]);

    // vertexData.state_ = state;
    // // Add to Sparse graph
    // sparseGraph_->addVertexFromFile(vertexData);

    // Add to Sparse graph
    VertexType type = static_cast<VertexType>(vertexData.type_);
    sparseGraph_->addVertex(state, type, 0);

    // Feedback
    if ((i + 1) % feedbackFrequency == 0)
      std::cout << std::fixed << std::setprecision(0) << (i / double(numVertices)) * 100.0 << "% " << std::flush;
  }
  std::cout << std::endl;
}

void BoltStorage::loadEdges(unsigned int numEdges, boost::archive::binary_iarchive &ia)
{
  OMPL_INFORM("Loading %u edges from file", numEdges);
  std::size_t feedbackFrequency = std::max(10.0, numEdges / 10.0);

  std::cout << "         Edges loaded: ";
  for (unsigned int i = 0; i < numEdges; ++i)
  {
    BoltEdgeData edgeData;
    ia >> edgeData;

    // Note: we increment all vertex indexes by the number of query vertices
    const SparseVertex v1 = edgeData.endpoints_.first += numQueryVertices_;
    const SparseVertex v2 = edgeData.endpoints_.second += numQueryVertices_;

    // Error check
    BOOST_ASSERT_MSG(v1 <= sparseGraph_->getNumVertices(), "Vertex 1 out of range of possible vertices");
    BOOST_ASSERT_MSG(v2 <= sparseGraph_->getNumVertices(), "Vertex 2 out of range of possible vertices");
    BOOST_ASSERT_MSG(v1 != v2, "Vertices of an edge loaded from file are the same");

    // Add
    EdgeType type = static_cast<EdgeType>(edgeData.type_);
    std::size_t indent = 2;
    sparseGraph_->addEdge(v1, v2, type, indent);

    // We deserialized the edge object pointer, and we own it.
    // Since addEdge copies the object, it is safe to free here.
    // delete edgeData.e_;

    // Feedback
    if ((i + 1) % feedbackFrequency == 0)
      std::cout << std::setprecision(0) << (i / double(numEdges)) * 100.0 << "% " << std::flush;
  }
  std::cout << std::endl;
}

}  // namespace bolt
}  // namespace tools
}  // namespace ompl
