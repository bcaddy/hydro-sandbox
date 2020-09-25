/*!
 * \file PerfTimer.h
 * \author Robert 'Bob' Caddy (rvc@pitt.edu)
 * \brief Contains a timing class for simple performance timing and statistics
 * of code
 * \version 1.0
 * \date 2020-09-24
 *
 * \copyright Copyright (c) 2020
 *
 */
#ifndef __PERFTIMER__H__
#define __PERFTIMER__H__

#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <algorithm>

/*!
 * \brief A class for timing pieces of code
 * \details To use this timer class simply initialize it with the name you want
 * to assign that timer, call PerfTimer::startTimer() when you want to start the
 * timer and PerfTimer::stopTimer() to end the timer. You can call start/stop
 * pairs as often as you want, each time difference will be saved and used for
 * final statistics. Statistics are accessed by calling PerfTimer::reportStats()
 * which writes the stats to std::cout or by calling PerfTimer::saveTimingData()
 * which saves the statistics and raw time differences to a given csv file. Note
 * that it will overwrite without asking
 *
 */
class PerfTimer
{
private:
    /// Stores the start time
    std::chrono::high_resolution_clock::time_point _startTime;

    /// Stores time differences in nanoseconds
    std::vector<double> _timeDiff;

    /// Indicates if a timer is active. If there is an active timer then it's
    /// `false`, else it is `true`
    bool _activeTimer;

    /// The name of the timer. To be printed out in the final output
    std::string _name;

    /*!
     * \brief Figure out the proper units for a given time.
     * \details Determines the proper units to use for time (nanoseconds,
     * microseconds, milliseconds, seconds, minutes, or hours) and returns the
     * time scaled to that unit along with the unit itself
     *
     * \param[in,out] time
     * \param[out] unit
     */
    void _converter(double &time, std::string &unit)
    {
        if (time <= 1.0E3)  // less than a microsecond
        {
            // Return time in nanoseconds
            unit = "ns";
        }
        else if (time <= 1.0E6)  // less than a millisecond
        {
            // Return time in microseconds
            time *= 1.E-3;
            unit = "\u00B5s";
        }
        else if (time <= 1.0E9)  // less than a second
        {
            // Return time in milliseconds
            time *= 1.E-6;
            unit = "ms";
        }
        else if (time <= 6.0E11)   // less than 10 minutes
        {
            // Return time in seconds
            time *= 1.E-9;
            unit = "s";
        }
        else if (time <= 1.08E13)  // less than 3 hours
        {
            // Return time in minutes
            time *= (1.E-9) / 60. ;
            unit = "mins";
        }
        else   // greater than 3 hours
        {
            // Return time in hours
            time *= (1.E-9) / 3600. ;
            unit = "hrs";
        }

    };

public:
    /*!
     * \brief Start the timer
     * \details First checks if a timer is active. If the timer is active then
     * nothing is done and an error message is printed. Else it starts the timer
     */
    void startTimer()
    {
        // Check if a timer is already running
        if (_activeTimer)
        {
            std::cout << _name << "::timer is already active. No action taken";
        }
        else
        {
            _startTime = std::chrono::high_resolution_clock::now();
        }
    };

    /*!
     * \brief Stop the timer
     * \details Stops the timer, appends the change in time to a private
     * std::vector in nanoseconds, then sets _activeTimer to false
     */
    void stopTimer()
    {
        // Compute the time difference and write it to _timeDiff
        double diff = static_cast<double> ((std::chrono::high_resolution_clock::now() - _startTime).count());
        _timeDiff.push_back(diff);

        _activeTimer = false;
    };

    /*!
     * \brief Compute and print out all the statistics for the timer
     *
     * \param[in] outStream What stream to write out to. Defaults to std::cout
     */
    void reportStats(std::ostream &outStream = std::cout)
    {
        // Compute statistics in nanoseconds
        double totalTime = std::accumulate(_timeDiff.begin(), _timeDiff.end(), 0.);
        double avgTime   = totalTime / static_cast<double>(_timeDiff.size());
        auto   minMax    = std::minmax_element(_timeDiff.begin(),_timeDiff.end());
        double minTime   = *minMax.first ;
        double maxTime   = *minMax.second;

        // Convert values
        std::string totalTimeUnit, avgTimeUnit, minTimeUnit, maxTimeUnit;
        _converter(totalTime, totalTimeUnit);
        _converter(avgTime,   avgTimeUnit);
        _converter(minTime,   minTimeUnit);
        _converter(maxTime,   maxTimeUnit);

        outStream << "Timer name: " << _name << std::endl  << "  " <<
        "Number of trials: " << _timeDiff.size()           << ", " <<
        "Total time: "       << totalTime << totalTimeUnit << ", " <<
        "Average Time: "     << avgTime   << avgTimeUnit   << ", " <<
        "Fastest Run: "      << minTime   << minTimeUnit   << ", " <<
        "Slowest Run: "      << maxTime   << maxTimeUnit   << std::endl;
    };

    /*!
     * \brief Writes all the time differences (in nanoseconds) to the given file
     * along with timer name and timer statistics.
     *
     * \details File format is .csv with two header lines. The two headers lines
     * are what is written out by PerfTimer::reportStats and the third line is
     * the measure time differences for each use of the timer.
     *
     * \param[in] filePath The path at which to write the output file. If the
     * file already exists it will be overwritten without asking.
     */
    void saveTimingData(std::string filePath)
    {
        // Open the file
        std::ofstream saveFile(filePath);

        // Write out the header info to the file
        reportStats(saveFile);

        // Write the vector to the file
        saveFile << _timeDiff[0];
        if (_timeDiff.size() > 0)
        {
            for (size_t i = 1; i < _timeDiff.size(); i++)
            {
                saveFile << "," << _timeDiff[i];
            }
        }
        saveFile << std::endl;

        // Close the file
        saveFile.close();
    };

    /*!
     * \brief Construct a new Perf Timer object
     *
     * \param[in] name The name of the timer
     */
    PerfTimer(std::string name);
    ~PerfTimer() = default;
};

PerfTimer::PerfTimer(std::string name)
    :
    _activeTimer(false),
    _name(name)
{
}
#endif  //!__PERFTIMER__H__