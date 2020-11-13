#!/usr/bin/env python3
"""
================================================================================
 Written by Robert Caddy.  Created on Fri May 22 14:49:13 2020

 plot and animate results from a 1D MHD simulation

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
    global densityData, pressureData, ieData, positions
    global velocityDataX, velocityDataY, velocityDataZ
    global magneticDataX, magneticDataY, magneticDataZ

    global intSamples, index, initFrames, initIndex

    global densityPlot, pressurePlot, iePlot, fig, supTitleText
    global velocityPlotX, velocityPlotY, velocityPlotZ
    global magneticPlotX, magneticPlotY, magneticPlotZ

    # Load file
    densityData = np.loadtxt("../data/Density.csv", delimiter=",")
    momentumDataX = np.loadtxt("../data/MomentumX.csv", delimiter=",")
    momentumDataY = np.loadtxt("../data/MomentumY.csv", delimiter=",")
    momentumDataZ = np.loadtxt("../data/MomentumZ.csv", delimiter=",")
    magneticDataX = np.loadtxt("../data/MagneticX.csv", delimiter=",")
    magneticDataY = np.loadtxt("../data/MagneticY.csv", delimiter=",")
    magneticDataZ = np.loadtxt("../data/MagneticZ.csv", delimiter=",")
    energyData = np.loadtxt("../data/Energy.csv", delimiter=",")

    # sim info
    simPhysLength = 1.
    simNumCells   = len(densityData[0,:])
    simNumSteps   = len(densityData[:,0])
    positions     = np.linspace(0.,simPhysLength,simNumCells)
    gamma         = 5./3.

    # Plot Settings
    supTitleText  = "Time Evolution of Initial Conditions Using Euler Equations"
    densityColor  = 'blue'                      # color of the density plot
    velocityColor = 'purple'                    # color of the velocity plots
    magneticColor = 'yellow'                    # color of the magnetic field plots
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
    velocityDataX = momentumDataX / densityData
    velocityDataY = momentumDataY / densityData
    velocityDataZ = momentumDataZ / densityData

    # Compute squares
    velocitySquared = velocityDataX**2 + velocityDataY**2 + velocityDataZ**2
    magneticSquared = magneticDataX**2 + magneticDataY**2+ magneticDataZ**2

    # Compute pressures
    pressureData = (gamma - 1) * (energyData
        - 0.5 * densityData * (velocitySquared)
        - 0.5 * (magneticSquared))

    # Compute the specific internal energy
    ieData = (1/densityData) * (energyData - 0.5 * densityData * velocitySquared - 0.5 * magneticSquared)
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
    pad = np.max(np.abs([velocityDataX.min(), velocityDataX.max()])) * 0.05
    velocityLowLimX = velocityDataX.min() - pad
    velocityHighLimX = velocityDataX.max() + pad
    pad = np.max(np.abs([velocityDataY.min(), velocityDataY.max()])) * 0.05
    velocityLowLimY = velocityDataY.min() - pad
    velocityHighLimY = velocityDataY.max() + pad
    pad = np.max(np.abs([velocityDataZ.min(), velocityDataZ.max()])) * 0.05
    velocityLowLimZ = velocityDataZ.min() - pad
    velocityHighLimZ = velocityDataZ.max() + pad

    # magnetic`
    pad = np.max(np.abs([magneticDataX.min(), magneticDataX.max()])) * 0.05
    magneticLowLimX = magneticDataX.min() - pad
    magneticHighLimX = magneticDataX.max() + pad
    pad = np.max(np.abs([magneticDataY.min(), magneticDataY.max()])) * 0.05
    magneticLowLimY = magneticDataY.min() - pad
    magneticHighLimY = magneticDataY.max() + pad
    pad = np.max(np.abs([magneticDataZ.min(), magneticDataZ.max()])) * 0.05
    magneticLowLimZ = magneticDataY.min() - pad
    magneticHighLimZ = magneticDataY.max() + pad

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
    # Create 9 subplots in a column
    fig, subPlot = plt.subplots(3, 3)

    # Super Title
    fig.suptitle(supTitleText)

    # Shared x-label
    subPlot[1,0].set_xlabel("Position")
    subPlot[1,1].set_xlabel("Position")
    subPlot[1,2].set_xlabel("Position")

    # Density subplot
    subPlot[0,0].set_ylim(densityLowLim, densityHighLim)
    subPlot[0,0].set_ylabel("Density")
    subPlot[0,0].grid(True)

    # Pressure subplot
    subPlot[0,1].set_ylim(pressureLowLim, pressureHighLim)
    subPlot[0,1].set_ylabel("Pressure")
    subPlot[0,1].grid(True)

    # Specific Internal Energy subplot
    subPlot[0,2].set_ylim(ieLowLim, ieHighLim)
    subPlot[0,2].set_ylabel("Internal Energy")
    subPlot[0,2].grid(True)

    # Velocity subplots
    subPlot[1,0].set_ylim(velocityLowLimX, velocityHighLimX)
    subPlot[1,0].set_ylabel("X Velocity")
    subPlot[1,0].grid(True)
    subPlot[1,1].set_ylim(velocityLowLimY, velocityHighLimY)
    subPlot[1,1].set_ylabel("Y Velocity")
    subPlot[1,1].grid(True)
    subPlot[1,2].set_ylim(velocityLowLimZ, velocityHighLimZ)
    subPlot[1,2].set_ylabel("Z Velocity")
    subPlot[1,2].grid(True)

    # Magnetic Field subplots
    subPlot[2,0].set_ylim(magneticLowLimX, magneticHighLimX)
    subPlot[2,0].set_ylabel("X Magnetic Field")
    subPlot[2,0].grid(True)
    subPlot[2,1].set_ylim(magneticLowLimY, magneticHighLimY)
    subPlot[2,1].set_ylabel("Y Magnetic Field")
    subPlot[2,1].grid(True)
    subPlot[2,2].set_ylim(magneticLowLimZ, magneticHighLimZ)
    subPlot[2,2].set_ylabel("Z Magnetic Field")
    subPlot[2,2].grid(True)

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
    pressurePlot, = subPlot[0,1].plot(positions,
                           pressureData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = pressureColor,
                           label      = "Pressure"
                           )
    iePlot,       = subPlot[0,2].plot(positions,
                           ieData[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = ieColor,
                           label      = "specific internal energy"
                           )

    velocityPlotX, = subPlot[1,0].plot(positions,
                           velocityDataX[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = velocityColor,
                           label      = "X Velocity"
                           )
    velocityPlotY, = subPlot[1,1].plot(positions,
                           velocityDataY[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = velocityColor,
                           label      = "Y Velocity"
                           )
    velocityPlotZ, = subPlot[1,2].plot(positions,
                           velocityDataZ[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = velocityColor,
                           label      = "Z Velocity"
                           )

    magneticPlotX, = subPlot[2,0].plot(positions,
                           magneticDataX[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = magneticColor,
                           label      = "X Magnetic Field"
                           )
    magneticPlotY, = subPlot[2,1].plot(positions,
                           magneticDataY[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = magneticColor,
                           label      = "Y Magnetic Field"
                           )
    magneticPlotZ, = subPlot[2,2].plot(positions,
                           magneticDataZ[0,:],
                           linestyle  = linestyle,
                           linewidth  = linewidth,
                           marker     = marker,
                           markersize = markersize,
                           color      = magneticColor,
                           label      = "Z Magnetic Field"
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
    pressurePlot.set_data(positions, pressureData[idx,:])
    iePlot.set_data(positions, ieData[idx,:])

    velocityPlotX.set_data(positions, velocityDataX[idx,:])
    velocityPlotY.set_data(positions, velocityDataY[idx,:])
    velocityPlotZ.set_data(positions, velocityDataZ[idx,:])

    magneticPlotX.set_data(positions, magneticDataX[idx,:])
    magneticPlotY.set_data(positions, magneticDataY[idx,:])
    magneticPlotZ.set_data(positions, magneticDataZ[idx,:])


    fig.suptitle(f"{supTitleText} \n Time Step: {idx}")

    if initIndex > initFrames:
        index += 1
    else:
        initIndex += 1


main()

print(f'Time to execute: {round(default_timer()-start,2)} seconds')