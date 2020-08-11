#!/usr/bin/env python3
"""
================================================================================
 Written by Robert Caddy.  Created on Fri May 22 14:49:13 2020

 plot and animate results from Euler equation solver

 Dependencies:
     numpy
     timeit
     matplotlib

 Changelog:
     Version 1.0 - First Version
================================================================================
"""

from timeit import default_timer

import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import sys

matplotlib.use("Agg")

plt.close('all')
start = default_timer()

def main():
    # ==========================================================================
    # Settings and Setup
    # ==========================================================================
    # Check for CLI arguments
    if len(sys.argv) == 2:
        cliFPS = int(sys.argv[1])
    elif len(sys.argv) > 2:
        raise TypeError(f"Too many ({len(sys.argv)}) command line arguments given")

    # Set global variable for the animation functions
    global densityData, velocityData, pressureData, ieData, positions
    global intSamples, index, initFrames, initIndex
    global densityPlot, velocityPlot, pressurePlot, iePlot, fig, supTitleText

    # Load file
    densityData = np.loadtxt("../data/Density.csv", delimiter=",")
    momentumData = np.loadtxt("../data/Momentum.csv", delimiter=",")
    energyData = np.loadtxt("../data/Energy.csv", delimiter=",")

    # sim info
    simPhysLength = 1.
    simNumCells   = len(densityData[0,:])
    simNumSteps   = len(densityData[:,0])
    positions     = np.linspace(0.,simPhysLength,simNumCells)
    gamma         = 1.4

    # Plot Settings
    supTitleText  = "Time Evolution of Initial Conditions Using Euler Equations"
    densityColor  = 'blue'                      # color of the density plot
    velocityColor = 'purple'                    # color of the velocity plot
    pressureColor = 'green'                     # color of the pressure plot
    ieColor       = 'red'                       # color of the specific internal energy plot
    linestyle     = '-'                         # The line style
    linewidth     = 0.5                         # How wide to make the lines
    marker        = "."                         # Marker kind for points
    markersize    = 3                           # Size of the marker

    # Video Settings
    OutFile       = "output-euler.mp4"          # Output filename
    Duration      = 10.                         # How long the video is in seconds
    dpi           = 150                         # Dots per inch
    index         = 0                           # Initialize index
    initIndex     = 0                           # Index for init frames
    fps           = cliFPS if ("cliFPS" in locals()) else 24  # Framerate
    FrameTime     = (1./fps) * 1000             # Frametime in milliseconds
    totFrames     = int(fps * Duration)         # Total number of frames (floor)
    initFrames    = fps                         # Number of frames for the initial conditions

      # Compute which time steps to plot
    if simNumSteps >= totFrames:
        floatSamples = np.arange(0, simNumSteps, simNumSteps/totFrames)
        intSamples   = np.asarray(np.floor(floatSamples), dtype="int")
    else:  # if the number of simulation steps is less than the total number of frames
        totFrames  = simNumSteps
        Fps        = int(totFrames/Duration)
        FrameTime  = (1./Fps) * 1000
        intSamples = np.arange(0, simNumSteps, 1, dtype="int")
    # ==========================================================================
    # End Settings and Setup
    # Compute Primitive Variables
    # ==========================================================================
    # Compute velocities
    velocityData = momentumData / densityData

    # Compute pressures
    pressureData = (gamma - 1) * (energyData - 0.5 * densityData * (velocityData**2))

    # Compute the specific internal energy
    ieData = pressureData / ((gamma - 1) * densityData)
    # ==========================================================================
    # End Computing Primitive Variables
    # Compute Limits
    # ==========================================================================
    # Find mins and maxes for setting the limits of the plot
    # Density
    pad = np.max(np.abs([densityData.min(), densityData.max()])) * 0.05
    densityLowLim = densityData.min() - pad
    densityHighLim = densityData.max() + pad

    # Velocity`
    pad = np.max(np.abs([velocityData.min(), velocityData.max()])) * 0.05
    velocityLowLim = velocityData.min() - pad
    velocityHighLim = velocityData.max() + pad

    # Pressure
    pad = np.max(np.abs([pressureData.min(), pressureData.max()])) * 0.05
    pressureLowLim = pressureData.min() - pad
    pressureHighLim = pressureData.max() + pad

    # Specific Internal Energy
    pad = np.max(np.abs([ieData.min(), ieData.max()])) * 0.05
    ieLowLim = ieData.min() - pad
    ieHighLim = ieData.max() + pad
    # ==========================================================================
    # End Computing Limits
    # Setup Plots
    # ==========================================================================
    # Create 4 subplots in a column
    fig, subPlot = plt.subplots(2, 2)

    # Super Title
    fig.suptitle(supTitleText)

    # Shared x-label
    subPlot[1,0].set_xlabel("Position")
    subPlot[1,1].set_xlabel("Position")

    # Density subplot
    subPlot[0,0].set_ylim(densityLowLim, densityHighLim)
    subPlot[0,0].set_ylabel("Density")
    subPlot[0,0].grid(True)

    # Velocity subplot
    subPlot[0,1].set_ylim(velocityLowLim, velocityHighLim)
    subPlot[0,1].set_ylabel("Velocity")
    subPlot[0,1].grid(True)

    # Pressure subplot
    subPlot[1,0].set_ylim(pressureLowLim, pressureHighLim)
    subPlot[1,0].set_ylabel("Pressure")
    subPlot[1,0].grid(True)

    # Specific Internal Energy subplot
    subPlot[1,1].set_ylim(ieLowLim, ieHighLim)
    subPlot[1,1].set_ylabel("Internal Energy")
    subPlot[1,1].grid(True)

    # Set plots
    densityPlot, = subPlot[0,0].plot(positions,
                           densityData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = densityColor,
                           label      = "Density"
                           )
    velocityPlot, = subPlot[0,1].plot(positions,
                           velocityData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = velocityColor,
                           label      = "Velocity"
                           )
    pressurePlot, = subPlot[1,0].plot(positions,
                           pressureData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = pressureColor,
                           label      = "Pressure"
                           )
    iePlot, = subPlot[1,1].plot(positions,
                           ieData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = ieColor,
                           label      = "specific internal energy"
                           )

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.subplots_adjust(top=0.88)
    # ==========================================================================
    # End Setup Plots
    # ==========================================================================

    # ==========================================================================
    # Make the animation
    # ==========================================================================
    simulation = animation.FuncAnimation(fig,
                                         NewFrame,
                                         blit = False,
                                         frames = totFrames+initFrames,
                                         interval = FrameTime,
                                         repeat = False)

    simulation.save(filename=OutFile, fps=fps, dpi=dpi)
    print(f"\nAnimation complete. Framerate: {fps} fps")
    # ==========================================================================
    # End Make the animation
    # ==========================================================================

def NewFrame(self):
    """
    This function generates the plotting for each individual frame
    """
    global index, initFrames, initIndex

    idx = intSamples[index]

    densityPlot.set_data(positions, densityData[idx,:])
    velocityPlot.set_data(positions, velocityData[idx,:])
    pressurePlot.set_data(positions, pressureData[idx,:])
    iePlot.set_data(positions, ieData[idx,:])

    fig.suptitle(f"{supTitleText} \n Time Step: {idx}")

    if initIndex > initFrames:
        index += 1
    else:
        initIndex += 1


main()

print(f'Time to execute: {round(default_timer()-start,2)} seconds')
