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
   Desc:   Use an external program to visualize search
*/

#ifndef OMPL_TOOLS_BOLT_VISUALIZER_
#define OMPL_TOOLS_BOLT_VISUALIZER_

#include <ompl/base/Planner.h>
#include <functional>

namespace ompl
{
namespace base
{
/**
   @anchor Visualizer
   @par Short description
   Use an external program to visualize search
*/

/// @cond IGNORE
OMPL_CLASS_FORWARD(Visualizer);
/// @endcond

/** \class ompl::tools::VisualizerPtr
    \brief A boost shared pointer wrapper for ompl::tools::Visualizer */

// Colors that correspond with rviz_visual_tools
enum colors
{
  BLACK = 0,
  BROWN = 1,
  BLUE = 2,
  CYAN = 3,
  GREY = 4,
  DARK_GREY = 5,
  GREEN = 6,
  LIME_GREEN = 7,
  MAGENTA = 8,
  ORANGE = 9,
  PURPLE = 10,
  RED = 11,
  PINK = 12,
  WHITE = 13,
  YELLOW = 14,
  TRANSLUCENT = 15,
  TRANSLUCENT_LIGHT = 16,
  TRANSLUCENT_DARK = 17,
  RAND = 18,
  CLEAR = 19,
  DEFAULT = 20 // i.e. 'do not change default color'
};

/**
 * \brief Visualization callback hook for external debugging of states
 * \param v - state to visualize
 * \param type - different modes/colors/sizes for visualizing the state
 * \param extraData - special extra data for showing a range around a state
 */
typedef std::function<void(const ompl::base::State* v, std::size_t type, double extraData)> VizState;
/**
 * \brief Visualization callback hook for external debugging of edges
 * \param v1 - from state that marks beginning of edge
 * \param v2 - to state that marks end of edge
 * \param cost - signifies the weight of the edge using color
 */
typedef std::function<void(const ompl::base::State* v1, const ompl::base::State* v2, double cost)> VizEdge;

/**
 * \brief Visualization callback hook for external debugging of paths
 * \param path
 * \param type - style of line
 */
typedef std::function<void(const ompl::base::PathPtr path, std::size_t type)> VizPath;

/**
 * \brief Visualization callback hook for external debugging that triggers the visualizer to render/publish
 */
typedef std::function<void(void)> VizTrigger;

/**
 * \brief Visualization callback hook for external debugging that tells OMPL when to end early by SIG INT
 */
// typedef std::function<bool(void)> VizCheckSigInt;

/** \brief Use an external program to visualize search */
class Visualizer
{
public:
  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz1State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz1State_)
      viz1State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz1Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz1Edge_)
      viz1Edge_(stateA, stateB, cost);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz1Path(const base::PathPtr path, std::size_t type)
  {
    if (viz1Path_)
      viz1Path_(path, type);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz1Trigger()
  {
    if (viz1Trigger_)
      viz1Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz1Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz1State_ = vizState;
    viz1Edge_ = vizEdge;
    viz1Path_ = vizPath;
    viz1Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz2State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz2State_)
      viz2State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz2Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz2Edge_)
      viz2Edge_(stateA, stateB, cost);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz2Path(const base::PathPtr path, std::size_t type)
  {
    if (viz2Path_)
      viz2Path_(path, type);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz2Trigger()
  {
    if (viz2Trigger_)
      viz2Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz2Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz2State_ = vizState;
    viz2Edge_ = vizEdge;
    viz2Path_ = vizPath;
    viz2Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz3State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz3State_)
      viz3State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz3Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz3Edge_)
      viz3Edge_(stateA, stateB, cost);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz3Path(const base::PathPtr path, std::size_t type)
  {
    if (viz3Path_)
      viz3Path_(path, type);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz3Trigger()
  {
    if (viz3Trigger_)
      viz3Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz3Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz3State_ = vizState;
    viz3Edge_ = vizEdge;
    viz3Path_ = vizPath;
    viz3Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz4State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz4State_)
      viz4State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz4Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz4Edge_)
      viz4Edge_(stateA, stateB, cost);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz4Path(const base::PathPtr path, std::size_t type)
  {
    if (viz4Path_)
      viz4Path_(path, type);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz4Trigger()
  {
    if (viz4Trigger_)
      viz4Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz4Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz4State_ = vizState;
    viz4Edge_ = vizEdge;
    viz4Path_ = vizPath;
    viz4Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz5State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz5State_)
      viz5State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz5Path(const base::PathPtr path, std::size_t type)
  {
    if (viz5Path_)
      viz5Path_(path, type);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz5Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz5Edge_)
      viz5Edge_(stateA, stateB, cost);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz5Trigger()
  {
    if (viz5Trigger_)
      viz5Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz5Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz5State_ = vizState;
    viz5Edge_ = vizEdge;
    viz5Path_ = vizPath;
    viz5Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz6State(const ompl::base::State* state, std::size_t type, double extraData)
  {
    if (viz6State_)
      viz6State_(state, type, extraData);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz6Path(const base::PathPtr path, std::size_t type)
  {
    if (viz6Path_)
      viz6Path_(path, type);
  }

  /** \brief Visualize planner's data during runtime, externally, using a function callback */
  void viz6Edge(const ompl::base::State* stateA, const ompl::base::State* stateB, double cost)
  {
    if (viz6Edge_)
      viz6Edge_(stateA, stateB, cost);
  }

  /** \brief Trigger visualizer to publish graphics */
  void viz6Trigger()
  {
    if (viz6Trigger_)
      viz6Trigger_();
  }

  /** \brief Set the callback to visualize/publish a planner's progress */
  void setViz6Callbacks(ompl::base::VizState vizState, ompl::base::VizEdge vizEdge, ompl::base::VizPath vizPath,
                        ompl::base::VizTrigger vizTrigger)
  {
    viz6State_ = vizState;
    viz6Edge_ = vizEdge;
    viz6Path_ = vizPath;
    viz6Trigger_ = vizTrigger;
  }

  // ------------------------------------------------------------------------------------

  /** \brief State callback by vizID */
  void vizState(const std::size_t vizID, const ompl::base::State* state, std::size_t type, double extraData)
  {
    switch (vizID)
    {
      case 1:
        viz1State(state, type, extraData);
        break;
      case 2:
        viz2State(state, type, extraData);
        break;
      case 3:
        viz3State(state, type, extraData);
        break;
      case 4:
        viz4State(state, type, extraData);
        break;
      case 5:
        viz5State(state, type, extraData);
        break;
      case 6:
        viz6State(state, type, extraData);
        break;
      default:
        OMPL_ERROR("Incorrect visualizer ID specified: %u", vizID);
    }
  }

  /** \brief Trigger visualizer to publish graphics by vizID */
  void vizTrigger(const std::size_t vizID)
  {
    switch (vizID)
    {
      case 1:
        viz1Trigger();
        break;
      case 2:
        viz2Trigger();
        break;
      case 3:
        viz3Trigger();
        break;
      case 4:
        viz4Trigger();
        break;
      case 5:
        viz5Trigger();
        break;
      case 6:
        viz6Trigger();
        break;
      default:
        OMPL_ERROR("Incorrect visualizer ID specified: %u", vizID);
    }
  }

private:
  /** \brief Optional st callbacks for data */
  ompl::base::VizState viz1State_;
  ompl::base::VizEdge viz1Edge_;
  ompl::base::VizPath viz1Path_;
  ompl::base::VizTrigger viz1Trigger_;

  /** \brief Optional 2nd callbacks for data */
  ompl::base::VizState viz2State_;
  ompl::base::VizEdge viz2Edge_;
  ompl::base::VizPath viz2Path_;
  ompl::base::VizTrigger viz2Trigger_;

  /** \brief Optional 3rd callbacks for data */
  ompl::base::VizState viz3State_;
  ompl::base::VizEdge viz3Edge_;
  ompl::base::VizPath viz3Path_;
  ompl::base::VizTrigger viz3Trigger_;

  /** \brief Optional 4th callbacks for data */
  ompl::base::VizState viz4State_;
  ompl::base::VizEdge viz4Edge_;
  ompl::base::VizPath viz4Path_;
  ompl::base::VizTrigger viz4Trigger_;

  /** \brief Optional 5th callbacks for data */
  ompl::base::VizState viz5State_;
  ompl::base::VizEdge viz5Edge_;
  ompl::base::VizPath viz5Path_;
  ompl::base::VizTrigger viz5Trigger_;

  /** \brief Optional 6th callbacks for data */
  ompl::base::VizState viz6State_;
  ompl::base::VizEdge viz6Edge_;
  ompl::base::VizPath viz6Path_;
  ompl::base::VizTrigger viz6Trigger_;

};  // end of class Visualizer
}  // namespace base
}  // namespace ompl

#endif  // OMPL_TOOLS_BOLT_VISUALIZER_
