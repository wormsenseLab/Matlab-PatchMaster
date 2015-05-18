Matlab-PatchMaster
==================

Matlab functions for importing and dealing with Patchmaster data for electrophysiology.  

These functions require sigTOOL for importing Patchmaster .dat files. 
The ImportHEKAtoMAT function must then be placed in the sigTOOL folder (a copy is included in this repo):
'sigTOOL\sigTOOL NeuroscienceToolkit\File\menu_Import\group_NeuroScience File Formats' 

sigTOOL must be run at least once in a given Matlab session before the function can be used.

Start with ImportPatchData.m, which calls on SplitSeries.m. AnalyzePatchData is currently a work in progress for my specific analysis.
