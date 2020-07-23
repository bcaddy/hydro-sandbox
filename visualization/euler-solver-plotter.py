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

matplotlib.use("Agg")

plt.close('all')
start = default_timer()

def main():
    # ==========================================================================
    # Settings and Setup
    # ==========================================================================
    # Set global variable for the animation functions
    global densityData, velocityData, pressureData, ieData, positions
    global Stride, Index, InitFrames, InitIndex
    global densityPlot, velocityPlot, pressurePlot, iePlot

    # Load file
    densityData = np.loadtxt("../data/Density.csv", delimiter=",")
    velocityData = np.loadtxt("../data/Velocity.csv", delimiter=",")
    pressureData = np.loadtxt("../data/Pressure.csv", delimiter=",")

    # sim info
    SimPhysLength = 1.
    SimNumCells   = len(densityData[0,:])
    SimNumSteps   = len(densityData[:,0])
    positions     = np.linspace(0.,SimPhysLength,SimNumCells)
    gamma         = 5./3.

    # Animation Settings
    Duration      = 10.                         # How long the gif is in seconds
    Fps           = 60                          # Frames per second
    FrameTime     = (1./Fps) * 1000             # Frametime in milliseconds
    TotFrames     = int(Fps * Duration)         # Total number of frames (floor)
    Stride        = int(SimNumSteps/TotFrames)  # Choose every n frames
    dpi           = 300                         # Dots per inch
    OutFile       = "output-euler.mp4"          # Output filename
    Index         = 0                           # Initialize index
    InitFrames    = 10                          # Number of frames for the initial conditions
    InitIndex     = 0                           # Index for init frames
    densityColor  = 'blue'                      # color of the density plot
    velocityColor = 'purple'                    # color of the velocity plot
    pressureColor = 'green'                     # color of the pressure plot
    ieColor       = 'red'                       # color of the specific internal energy plot

    # Reset some animation settings for short simulations
    if SimNumSteps < TotFrames:
        TotFrames = SimNumSteps
        Fps       = TotFrames/Duration
        FrameTime = (1./Fps) * 1000
        Stride    = int(SimNumSteps/TotFrames)
    # ==========================================================================
    # End Settings and Setup
    # Compute specific internal energy and limits
    # ==========================================================================
    # Compute the specific internal energy
    ieData = np.zeros((SimNumSteps,SimNumCells))
    ieData = pressureData * 2#pressureData / ( densityData * (gamma - 1) )

    # Find mins and maxes for setting the limits of the plot
    # Density
    pad = np.max(np.abs([densityData.min(), densityData.max()])) * 0.05
    densityLowLim = densityData.min() - pad
    densityHighLim = densityData.max() + pad

    # Velocity
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
    # Setup Plots
    # ==========================================================================
    # Create 4 subplots in a column
    fig, subPlot = plt.subplots(4, 1, sharex = "col")

    # Super Title
    fig.suptitle(f"Time Evolution of Initial Conditions Using Euler Equations")

    # Shared x-label
    subPlot[-1].set_xlabel("Position")

    # Density subplot
    subPlot[0].set_ylim(densityLowLim, densityHighLim)
    subPlot[0].set_ylabel("Density")
    subPlot[0].grid(True)

    # Velocity subplot
    subPlot[1].set_ylim(velocityLowLim, velocityHighLim)
    subPlot[1].set_ylabel("Velocity")
    subPlot[1].grid(True)

    # Pressure subplot
    subPlot[2].set_ylim(pressureLowLim, pressureHighLim)
    subPlot[2].set_ylabel("Pressure")
    subPlot[2].grid(True)

    # Specific Internal Energy subplot
    subPlot[3].set_ylim(ieLowLim, ieHighLim)
    subPlot[3].set_ylabel("Specific \n Internal \n Energy")
    subPlot[3].grid(True)

    # Set plots
    densityPlot, = subPlot[0].plot(positions,
                           densityData[0,:],
                           linestyle = '-',
                           color     = densityColor,
                           label     = "Density"
                           )
    velocityPlot, = subPlot[1].plot(positions,
                           velocityData[0,:],
                           linestyle = '-',
                           color     = velocityColor,
                           label     = "Velocity"
                           )
    pressurePlot, = subPlot[2].plot(positions,
                           pressureData[0,:],
                           linestyle = '-',
                           color     = pressureColor,
                           label     = "Pressure"
                           )
    iePlot, = subPlot[3].plot(positions,
                           ieData[0,:],
                           linestyle = '-',
                           color     = ieColor,
                           label     = "specific internal energy"
                           )
    # ==========================================================================
    # End Setup Plots
    # ==========================================================================

    # ==========================================================================
    # Make the animation
    # ==========================================================================
    simulation = animation.FuncAnimation(fig,
                                         NewFrame,
                                         blit = False,
                                         frames = TotFrames+InitFrames,
                                         interval = FrameTime,
                                         repeat = False)

    simulation.save(filename=OutFile, fps=Fps, dpi=dpi)

    # ==========================================================================
    # End Make the animation
    # ==========================================================================

def NewFrame(self):
    """
    This function generates the plotting for each individual frame
    """
    global Index, InitFrames, InitIndex

    densityPlot.set_data(positions, densityData[Index,:])
    velocityPlot.set_data(positions, velocityData[Index,:])
    pressurePlot.set_data(positions, pressureData[Index,:])
    iePlot.set_data(positions, ieData[Index,:])


    if InitIndex > InitFrames:
        Index += Stride
    else:
        InitIndex += 1
main()
print(f'\nTime to execute: {round(default_timer()-start,2)} seconds')