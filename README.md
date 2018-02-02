Matlab-PatchMaster
==================

Matlab functions for importing and dealing with Patchmaster data for electrophysiology.

These functions are known to work for Patchmaster versions 2x73.1 through 2x90.1. They will most likely work on later versions, but I haven't tested that.

These functions require sigTOOL for importing Patchmaster .dat files. You can download sigTOOL here: http://sigtool.sourceforge.net/sigtool.html.
Follow the directions in the included pdf for installing it properly.

IMPORTANT:
The ImportHEKAtoMat function must then be placed in the sigTOOL folder (a copy is included in this repo):
'sigTOOL\sigTOOL NeuroscienceToolkit\File\menu_Import\group_NeuroScience File Formats' 

sigTOOL must be run at least once in a given Matlab session before the function can be used.

Start with the Import folder, which contains ImportPatchData, which calls on SplitSeries and ImportHEKAtoMat.

AnalyzePatchData and the rest of the functions in the Analysis folder are currently a work in progress for my specific analysis.
