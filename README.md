Matlab-PatchMaster
==================

Matlab functions for importing and dealing with Patchmaster data for electrophysiology.

These functions are known to work for Patchmaster versions 2x73.1 through 2x90.1. They will most likely work on later versions, but I haven't tested that.

These functions require sigTOOL for importing Patchmaster .dat files. You can download sigTOOL here: http://sigtool.sourceforge.net/sigtool.html.
Follow the directions in the included pdf for installing it properly.

IMPORTANT:
The ImportHEKAtoMat function (a modified version of ImportHEKA) must be placed in the sigTOOL folder:
'sigTOOL\sigTOOL NeuroscienceToolkit\File\menu_Import\group_NeuroScience File Formats' 

You can find the ImportHEKAtoMat function here: https://gist.github.com/sammykatta/190b38de5b6c8b37fd04ba9365364aea

sigTOOL must be run at least once in a given Matlab session before the function can be used.

Start with the Import folder, which contains ImportPatchData, which calls on SplitSeries and ImportHEKAtoMat.

AnalyzePatchData and the rest of the functions in the Analysis folder are currently a work in progress for my specific analyses.

------------------

Code for the Sanzeni et al. (2018) preprint "Tissue mechanics and somatosensory neural responses govern touch sensation in C. elegans" (https://www.biorxiv.org/content/10.1101/471904v1) and any published articles resulting from this manuscript can be found in the folder Analysis/SK. 

#### Important scripts for this manuscript are listed below:

**AnalyzePatchData.m** is a general script with examples for importing data and metadata, checking for bad sweeps, and running various types of analyses.

**PreIndentIgorExport.m** runs the analysis for displacement steps applied with or without a pre-indentation, considering both absolute stimulator position and position relative to the pre-step, then outputs it to an Excel file for plotting with Igor.


#### Other functions necessary for understanding the details of analysis are listed below.

**IdAnalysis.m** takes the mean of technical replicates within a recording by grouping stimuli considered to be the same, calculates the peak current in response to each stimulus (on, off, repeated steps) within each mean sweep, and outputs a cell array with mean traces, stimulus parameters, peak current and kinetics. This function can be used when analyzing by stimulus displacement, position, speed, or inter-stimulus interval. It relies on newStepFind and findMRCs.

**newStepFind.m** is used for calculating the displacement, speed, and timecourse of the stimuli from an efference copy of the stimulus command signal after passing through the LPF-8 Bessel filter but before full amplification by the Crawford amplifier, allowing a more accurate representation of the stimulus than reconstruction from the stimTree in ephysData allows (also, stimTree functionality was added much later). It can also be used with the photodiode signal for greater accuracy if the photodiode signal has low noise.
    
**findMRCs.m** smooths traces given to it, sets a threshold for peak detection based on the baseline noise, and calculates activation and decay time constants by one of two specified methods.


Created by: Samata Katta
