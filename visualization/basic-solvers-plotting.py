#!/usr/bin/env python3
"""
================================================================================
 Written by Robert Caddy.  Created on Fri May 22 14:49:13 2020

 plot and animate results from Advection/Riemann solver

 Dependencies:
     numpy
     timeit
     donemusic
     matplotlib

 Changelog:
     Version 1.0 - First Version
================================================================================
"""

import numpy as np
from timeit import default_timer

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.close('all')
start = default_timer()

def main():
    # ==========================================================================
    # Settings and Setup
    # ==========================================================================
    # Set global variable for the animation functions
    global file, Stride, AdvectLine, Index, positions, InitFrames, InitIndex

    # Load file
    file = np.loadtxt("../data/results.csv", delimiter=",")

    # sim info
    SimPhysLength = 1.
    SimNumCells = len(file[0,:])
    SimNumSteps = len(file[:,0])
    positions = np.linspace(0.,SimPhysLength,SimNumCells)

    # Animation Settings
    Duration   = 10.                         # How long the gif is in seconds
    Fps        = 60                          # Frames per second
    FrameTime  = (1./Fps) * 1000             # Frametime in milliseconds
    TotFrames  = int(Fps * Duration)         # Total number of frames (floor)
    Stride     = int(SimNumSteps/TotFrames)  # Choose every n frames
    dpi        = 300                         # Dots per inch
    Color      = 'blue'                      # color of the solution
    OutFile    = "output-basic.mp4"        # Output filename
    Index      = 0                           # Initialize index
    InitFrames = 10                          # Number of frames for the initial conditions
    InitIndex  = 0                           # Index for init frames

    # Reset some animation settings for short simulations
    if SimNumSteps < TotFrames:
        TotFrames = SimNumSteps
        Fps       = TotFrames/Duration
        FrameTime = (1./Fps) * 1000
        Stride    = int(SimNumSteps/TotFrames)

    # Find mins and maxes for setting the limits of the plot
    pad = np.max(np.abs([file.min(), file.max()])) * 0.05
    small = file.min() - pad
    large = file.max() + pad
    # ==========================================================================
    # End Settings and Setup
    # ==========================================================================

    # ==========================================================================
    # Setup Plots
    # ==========================================================================
    f0 = plt.figure(num = 0)
    plt.title(f"Solution")

    plt.ylim(small, large)

    plt.xlabel("Position")
    plt.ylabel("Value")

    plt.tight_layout()

    plt.grid(True)

    # Set plots
    AdvectLine, = plt.plot(positions,
                           file[0,:],
                           linestyle = '-',
                           color     = Color,
                           label     = "Advection of Tophat"
                           )

    # Create legend
#    plt.legend()
    # ==========================================================================
    # End Setup Plots
    # ==========================================================================

    # ==========================================================================
    # Make the animation
    # ==========================================================================
    simulation = animation.FuncAnimation(f0,
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

    AdvectLine.set_data(positions,file[Index,:])

    if InitIndex > InitFrames:
        Index += Stride
    else:
        InitIndex += 1
main()
print(f'\nTime to execute: {round(default_timer()-start,2)} seconds')
