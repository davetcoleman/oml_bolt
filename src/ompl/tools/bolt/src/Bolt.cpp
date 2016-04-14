/*********************************************************************
 * Software License Agreement (SBD License)
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

/* Author: Dave Coleman */

#include <ompl/tools/bolt/Bolt.h>
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/base/PlannerStatus.h>
#include <ompl/util/Console.h>

namespace og = ompl::geometric;
namespace ob = ompl::base;
namespace ot = ompl::tools;

namespace ompl
{
namespace tools
{
namespace bolt
{
Bolt::Bolt(const base::SpaceInformationPtr &si) : ExperienceSetup(si)
{
    initialize();
}

Bolt::Bolt(const base::StateSpacePtr &space) : ExperienceSetup(space)
{
    initialize();
}

void Bolt::initialize()
{
    OMPL_INFORM("Initializing Bolt Framework");

    // Initalize visualizer class
    visual_.reset(new base::Visualizer());

    recallEnabled_ = true;
    scratchEnabled_ = true;
    filePath_ = "unloaded";

    // Load the experience database
    denseDB_.reset(new DenseDB(si_, visual_));

    // Load the Retrieve repair database. We do it here so that setRepairPlanner() works
    boltPlanner_ = ob::PlannerPtr(new BoltRetrieveRepair(si_, denseDB_));  // TODO(davetcoleman): pass in visual_

    // Load the discretize grid tool
    discretizer_.reset(new Discretizer(si_, denseDB_, denseDB_->getVisual()));

    OMPL_INFORM("Bolt Framework initialized.");
}

void Bolt::setup(void)
{
    if (!configured_ || !si_->isSetup() || !boltPlanner_->isSetup())
    {
        // Setup Space Information if we haven't already done so
        if (!si_->isSetup())
            si_->setup();

        // Setup planning from experience planner
        boltPlanner_->setProblemDefinition(pdef_);

        if (!boltPlanner_->isSetup())
            boltPlanner_->setup();

        // Setup database
        denseDB_->setup();

        // Set the configured flag
        configured_ = true;
    }
}

void Bolt::clear(void)
{
    if (boltPlanner_)
        boltPlanner_->clear();
    if (pdef_)
        pdef_->clearSolutionPaths();
}

void Bolt::setPlannerAllocator(const base::PlannerAllocator &pa)
{
    pa_ = pa;
    // note: the boltPlanner_ never uses the allocator so does not need to be reset
    configured_ = false;
}

base::PlannerStatus Bolt::solve(const base::PlannerTerminationCondition &ptc)
{
    // Setup again in case it has not been done yet
    setup();

    lastStatus_ = base::PlannerStatus::UNKNOWN;
    time::point start = time::now();

    // Warn if there are queued paths that have not been added to the experience database
    OMPL_INFORM("Num solved paths uninserted into the experience database in the post-proccessing queue: %u", queuedSolutionPaths_.size());

    // SOLVE
    lastStatus_ = boltPlanner_->solve(ptc);

    // Planning time
    planTime_ = time::seconds(time::now() - start);

    // Do logging
    logResults();

    return lastStatus_;
}

void Bolt::logResults()
{
    // Create log
    ExperienceLog log;
    log.planningTime = planTime_;

    // Record stats
    stats_.totalPlanningTime_ += planTime_;  // used for averaging
    stats_.numProblems_++;                   // used for averaging

    switch (static_cast<ompl::base::PlannerStatus::StatusType>(lastStatus_))
    {
        case base::PlannerStatus::TIMEOUT:
            stats_.numSolutionsTimedout_++;
            OMPL_ERROR("Bolt::solve(): TIMEOUT - No solution found after %f seconds", planTime_);
            // Logging
            log.planner = "neither_planner";
            log.result = "timedout";
            log.isSaved = "not_saved";
            break;
        case base::PlannerStatus::ABORT:
            stats_.numSolutionsTimedout_++;
            OMPL_ERROR("Bolt::solve(): ABORT - No solution found after %f seconds", planTime_);
            // Logging
            log.planner = "neither_planner";
            log.result = "abort";
            log.isSaved = "not_saved";
            break;
        case base::PlannerStatus::APPROXIMATE_SOLUTION:
            OMPL_ERROR("Bolt::solve(): Approximate - should not happen!");
            exit(-1);
            break;
        case base::PlannerStatus::EXACT_SOLUTION:
        {
            og::PathGeometric solutionPath = og::SimpleSetup::getSolutionPath();  // copied so that it is non-const

            std::cout << ANSI_COLOR_BLUE;
            std::cout << "Bolt Finished - solution found in " << std::setprecision(5) << planTime_ << " seconds with " << solutionPath.getStateCount() << " states" << std::endl;
            std::cout << ANSI_COLOR_RESET;

            // Error check for repeated states
            if (!checkRepeatedStates(solutionPath))
                exit(-1);

            // Stats
            stats_.numSolutionsFromRecall_++;

            // Make sure solution has at least 2 states
            if (solutionPath.getStateCount() < 2)
            {
                OMPL_INFORM("NOT saving to database because solution is less than 2 states long");
                stats_.numSolutionsTooShort_++;

                // Logging
                log.isSaved = "less_2_states";
                log.tooShort = true;
            }
            else
            {
                // Queue the solution path for future insertion into experience database (post-processing)
                queuedSolutionPaths_.push_back(solutionPath);
            }
        }
        break;
        default:
            OMPL_ERROR("Unknown status type: %u", lastStatus_);
            stats_.numSolutionsFailed_++;
            // Logging
            log.planner = "neither_planner";
            log.result = "failed";
            log.isSaved = "not_saved";
    }

    // Final log data
    // log.insertion_time = insertionTime; TODO fix this
    log.numVertices = denseDB_->getNumVertices();
    log.numEdges = denseDB_->getNumEdges();
    log.numConnectedComponents = 0;

    // Flush the log to buffer
    convertLogToString(log);
}

bool Bolt::checkRepeatedStates(const og::PathGeometric &path)
{
    for (std::size_t i = 1; i < path.getStateCount(); ++i)
    {
        if (si_->getStateSpace()->equalStates(path.getState(i - 1), path.getState(i)))
        {
            OMPL_ERROR("Duplicate state found on trajectory at %u out of %u", i, path.getStateCount());
            return false;
        }
    }
    return true;
}

base::PlannerStatus Bolt::solve(double time)
{
    ob::PlannerTerminationCondition ptc = ob::timedPlannerTerminationCondition(time);
    return solve(ptc);
}

bool Bolt::save()
{
    // setup(); // ensure the db has been loaded to the Experience DB
    return denseDB_->save(filePath_);
}

bool Bolt::saveIfChanged()
{
    // setup(); // ensure the db has been loaded to the Experience DB
    return denseDB_->saveIfChanged(filePath_);
}

void Bolt::printResultsInfo(std::ostream &out) const
{
    for (std::size_t i = 0; i < pdef_->getSolutionCount(); ++i)
    {
        out << "Solution " << i << "\t | Length: " << pdef_->getSolutions()[i].length_
            << "\t | Approximate: " << (pdef_->getSolutions()[i].approximate_ ? "true" : "false")
            << "\t | Planner: " << pdef_->getSolutions()[i].plannerName_ << std::endl;
    }
}

bool Bolt::loadOrGenerate()
{
    // Load from file or generate new grid
    if (denseDB_->getNumVertices() <= 1)  // the search verticie may already be there
    {
        if (!denseDB_->load(filePath_))  // load from file
        {
            OMPL_INFORM("No database loaded from file - generating new grid");

            // Benchmark runtime
            time::point startTime = time::now();

            discretizer_->generateGrid();

            // Benchmark runtime
            double duration = time::seconds(time::now() - startTime);
            OMPL_INFORM("Grid generation total time: %f seconds (%f min)", duration, duration / 60.0);

            return true;
        }
        return true;
    }
    OMPL_INFORM("Database already loaded");
    return true;
}

void Bolt::print(std::ostream &out) const
{
    if (si_)
    {
        si_->printProperties(out);
        si_->printSettings(out);
    }
    if (boltPlanner_)
    {
        boltPlanner_->printProperties(out);
        boltPlanner_->printSettings(out);
    }
    if (pdef_)
        pdef_->print(out);
}

void Bolt::printLogs(std::ostream &out) const
{
    SparseDBPtr sparseDB = denseDB_->getSparseDB();
    double vertPercent = sparseDB->getNumVertices() / double(denseDB_->getNumVertices()) * 100.0;
    double edgePercent = sparseDB->getNumEdges() / double(denseDB_->getNumEdges()) * 100.0;
    double solvedPercent = stats_.numSolutionsFromRecall_ / static_cast<double>(stats_.numProblems_) * 100.0;
    if (!recallEnabled_)
        out << "Scratch Planning Logging Results (inside Bolt Framework)" << std::endl;
    else
        out << "Bolt Framework Logging Results" << std::endl << std::setprecision(0);
    out << "  Solutions Attempted:           " << stats_.numProblems_ << std::endl;
    out << "    Solved:                      " << stats_.numSolutionsFromRecall_ << " (" << solvedPercent << "%)\n";
    out << "    Failed:                      " << stats_.numSolutionsFailed_ << std::endl;
    out << "    Timedout:                    " << stats_.numSolutionsTimedout_ << std::endl;
    out << "    Approximate:                 " << stats_.numSolutionsApproximate_ << std::endl;
    out << "  DenseDB                        " << std::endl;
    out << "    Vertices:                    " << denseDB_->getNumVertices() << std::endl;
    out << "    Edges:                       " << denseDB_->getNumEdges() << std::endl;
    out << "    Start/Goal States Added:     " << getRetrieveRepairPlanner().numStartGoalStatesAddedToDense_ << std::endl;
    out << "  SparseDB                       " << std::endl;
    out << "    Vertices:                    " << sparseDB->getNumVertices() << " (" << vertPercent << "%)" << std::endl;
    out << "    Edges:                       " << sparseDB->getNumEdges() << " (" << edgePercent << "%)" << std::endl;
    out << "    Regenerations:               " << sparseDB->numGraphGenerations_ << std::endl;
    out << "    Disjoint Samples Added:      " << sparseDB->numSamplesAddedForDisjointSets_ << std::endl;
    out << "    Sparse Delta Fraction:       " << sparseDB->sparseDeltaFraction_ << std::endl;
    out << "  Average planning time:         " << stats_.getAveragePlanningTime() << " seconds" << std::endl;
    out << "  Average insertion time:        " << stats_.getAverageInsertionTime() << " seconds" << std::endl;
    out << std::endl;
}

std::size_t Bolt::getExperiencesCount() const
{
    return denseDB_->getNumVertices();
}

void Bolt::getAllPlannerDatas(std::vector<ob::PlannerDataPtr> &plannerDatas) const
{
    //denseDB_->getAllPlannerDatas(plannerDatas);
}

void Bolt::convertPlannerData(const ob::PlannerDataPtr plannerData, og::PathGeometric &path)
{
    // Convert the planner data verticies into a vector of states
    for (std::size_t i = 0; i < plannerData->numVertices(); ++i)
        path.append(plannerData->getVertex(i).getState());
}

DenseDBPtr Bolt::getDenseDB()
{
    return denseDB_;
}

bool Bolt::doPostProcessing()
{
    OMPL_INFORM("Performing post-processing for %i queued solution paths", queuedSolutionPaths_.size());

    // Benchmark runtime
    time::point startTime = time::now();

    for (std::size_t i = 0; i < queuedSolutionPaths_.size(); ++i)
    {
        if (denseDB_->snapPathVerbose_)
            std::cout << "post processing path " << i << " of " << queuedSolutionPaths_.size()
                      << " -------------------------- " << std::endl;

        // Time to add a path to experience database
        if (!denseDB_->postProcessPath(queuedSolutionPaths_[i]))
        {
            OMPL_ERROR("Unable to save path");
        }
    }
    OMPL_INFORM("Finished inserting %u experience paths", queuedSolutionPaths_.size());

    // Remove all inserted paths from the queue
    queuedSolutionPaths_.clear();

    // Ensure graph doesn't get too popular
    if (denseDB_->getPopularityBiasEnabled())
        denseDB_->normalizeGraphEdgeWeights();

    // Show changes
    denseDB_->displayDatabase();

    // Recreate the sparse graph, too
    denseDB_->getSparseDB()->createSPARS();

    // Benchmark runtime
    double duration = time::seconds(time::now() - startTime);
    OMPL_INFORM(" - doPostProcessing() took %f seconds (%f hz)", duration, 1.0/duration);

    stats_.totalInsertionTime_ += duration; // used for averaging

    return true;
}
}  // namespace bolt
}  // namespace tools
}  // namespace ompl
