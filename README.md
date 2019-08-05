Matlab-PatchMaster
==================

Matlab functions for importing and dealing with Patchmaster data for electrophysiology.

These functions are known to work for Patchmaster versions 2x73.1 through 2x90.1. They will most likely work on later versions, but I haven't tested that.

These functions require sigTOOL for importing Patchmaster .dat files. You can download sigTOOL here: http://sigtool.sourceforge.net/sigtool.html.
Follow the directions in the included pdf for installing it properly.

IMPORTANT:
The ImportHEKAtoMat function (a modified version of ImportHEKA) must be placed in the sigTOOL folder (a copy is included in this repo):
'sigTOOL\sigTOOL NeuroscienceToolkit\File\menu_Import\group_NeuroScience File Formats' 

sigTOOL must be run at least once in a given Matlab session before the function can be used.

Start with the Import folder, which contains ImportPatchData, which calls on SplitSeries and ImportHEKAtoMat.

=================

Code for the Katta et al. (2018) preprint "Progressive recruitment of distal MEC-4 channels determines touch response strength in C. elegans" (https://www.biorxiv.org/content/10.1101/587014v1) and any published articles resulting from this manuscript can be found in the folder Analysis/SK. Useful scripts are listed below. Simulated data were generated separately and loaded into Matlab with this code for making plots.

   AnalyzePatchData.m is a general script with examples for importing data and metadata, checking for bad sweeps, and running various types of analyses.

   AntVsPost_VoltageAttenuation_180803.m runs the analysis for displacement steps applied either anterior or posterior to the cell body as seen in Fig. 1, 3, 4, and 5. Outputs data to plot I-d and Q-d curves, as well as activation and decay rates. Also includes plots of representative traces used in Fig. 3 and 9.

   Velocity_VoltageAtt_180913.m runs the analysis for trapezoidal stimuli applied anterior to the cell body as seen in Fig. 6 and 7. Outputs data for both onset and offset responses to plot I-d and Q-d curves, as well as activation and decay rates, and off/on ratio. Also includes plots of representative traces used in Fig. 6.

   SinePowerSpec.m runs the analysis for sinusoidal stimuli applied anterior to the cell body as seen in Fig. 8. Plots the mean power spectral density by frequency (including the stimulus photodiode trace). Loads simulated sine responses and plots mean steady-state current and rms for both experimental and simulated data. Also includes representative experimental traces used in Fig. 8 and 9.

   Sim_Distance_Comparison.m plots Fig. 9's simulated and experimental representative on/off traces, as well as simulated average responses for individual channels under different stimulus conditions.

   intVsDiss_Steps.m does a pilot comparison of the minimally dissected ("intact") preparation vs. the regular ("dissected") preparation.

Other functions necessary for understanding the details of analysis are listed below.

   IdAnalysis.m takes the mean of technical replicates within a recording by grouping stimuli considered to be the same, calculates the peak current in response to each stimulus (on, off, repeated steps) within each mean sweep, and outputs a cell array with mean traces, stimulus parameters, peak current and kinetics. This function can be used when analyzing by stimulus displacement, position, speed, or inter-stimulus interval. It relies on newStepFind and findMRCs.

   newStepFind.m is used for calculating the displacement, speed, and timecourse of the stimuli from an efference copy of the stimulus command signal after passing through the LPF-8 Bessel filter but before full amplification by the Crawford amplifier, allowing a more accurate representation of the stimulus than reconstruction from the stimTree in ephysData allows (also, stimTree functionality was added much later). It can also be used with the photodiode signal for greater accuracy if the photodiode signal has low noise.
    
   findMRCs.m smooths traces given to it, sets a threshold for peak detection based on the baseline noise, and calculates activation and decay time constants by one of two specified methods.

   FrequencyAnalysis.m is used for sine stimuli with flanking square steps. It calculates steady state current and rms during steady-state, as well as peak current responses to the steps, and attenuation of step responses. It also allows normalization of steady-state and rms to the peak current for a step.

Created by: Samata Katta
