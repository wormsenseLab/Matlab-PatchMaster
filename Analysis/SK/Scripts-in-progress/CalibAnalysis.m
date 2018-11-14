% Calibration Analysis Script

%% Make list of approved traces (by selecting traces to exclude)

protList ={'ProbeCalib'};

matchType = 'first';
strainList = {'TU2769'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 100];
resistCutoff = '<250';
extFilterFreq = 5;
includeFlag = 1;

calib5Cells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
     'stimFilterFrequencykHz', extFilterFreq, 'included', 1);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%% Run sweep selection GUI

ExcludeSweeps(ephysData, protList, calib5Cells, 'matchType', matchType, 'channel', 3);

clear protList matchType;

%% Calibration step sizes for reference
% (FAT Notebook I, pg. 215).

% ProbeCalibT5:
% [0, 2.5,  5.0,  7.5, 10.0,  12.5] um
% 
% Axial probe motion (at 18-deg angle):
% [0, 2.38, 4.76, 7.13, 9.51, 11.89] um
% 
% PD motion in Kinesis sequence:
% [0, .143, .285, .428, .571, .713] mm
% 
% ProbeCalibT5 is PD signal during probe motion (Sweep 1). 
% ProbeCalib reads out signal during PD motion.
%   Sweep 2 = probe and worm; Sweep 3 = worm only.
% Recordings in metadata database marked *Calib3Trig* have photodiode motion
% triggering acquisition, so steps are always at the same timepoint. This
% code cannot be used on recordings only marked Calib or Calib3.

% From list of recordings with 2.5kHz filtering and well-separated, 
% non-noisy calibration.    *calib_2p5kHz_picks(181113).xls*
% Recordings with sinusoids:
%   FAT170, FAT171, FAT236
% Recordings with trapezoids:
%   FAT236, 238
% Recordings with steps:
%   FAT116, FAT143, FAT236
% Recordings with speed-controlled steps:
%   FAT236, FAT237, FAT238, FAT239

% From list of recordings with 5kHz filtering and well-separated,
% non-noisy calibration.    *calib_5kHz_picks(181113).xls*
% Recordings with sinusoids:
%   FAT212, FAT218
% Recordings with trapezoids:
%   FAT214, FAT216, FAT217, FAT218, FAT219


%% Pull out example traces comparing stimulus command with PD signal 

% Since PD signal is not (yet?) useful for absolute displacement
% calibration, normalize to 1 based on maximum commanded displacement.

% What I want: for each protocol set, organize sweeps by stim speed/freq/size,
% take the average of the command trace and the PD trace, zero-subtract, 
% flip the PD trace, and normalize both. 
% Retain the ability to pull out individual traces instead of the mean.
% NEXT: Q, can I make use of existing functions to do this for other
% channels?

% Plot: overlay normalized PD/command traces for each speed/freq/size in
% subplots going down, on left and off right. (single column for sines).

