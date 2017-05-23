% findRampMRCs.m
%
% Uses stimWindow instead of stimStart to allow for ramp peak-finding

function [pk, pkLoc, pkThresh, varargout] = findRampMRCs(stimWindow, traceData,sf)

smoothWindow = sf; % n timepoints for moving average window for findPeaks (sf = 1ms worth)
threshTime = 100; % use first n ms of trace for setting noise threshold

% Number of timepoints to skip after stimulus onset to avoid the stimulus
% artifact in peak-finding (dependent on sampling frequency in kHz).
switch sf
    case 5
        artifactOffset = sf*2;
    case 10
        artifactOffset = 0;
    otherwise
        artifactOffset = sf;
end

% Smooth data with a moving average for peak finding
smooMean = -smooth(traceData',smoothWindow,'moving');

% Set threshold based on noise of the first 100ms of the trace
% (i.e., size of signal needed to be seen above that noise)
pkThresh = 1.5*thselect(smooMean(1:threshTime*sf),'rigrsure');


% Find MRC peaks if they exist, otherwise set peak amplitude as 0.
% Calculate decay constant tau based on single exponent fit.

%NEXT: Can this buffer change between rates? Does it need to? Or will the
%fastest time be enough?
[peaks, peakLocs] = findpeaks(smooMean(stimWindow(1)+artifactOffset:stimWindow(2)+(sf*20)),...
    'minpeakheight',pkThresh); %allows ramp, unlike findMRCs
if ~isempty(peaks)
    
    pk = max(peaks)*1E12;
    peakLocs = peakLocs(peaks==max(peaks));
    pkLoc = peakLocs(1) + stimWindow(1)+artifactOffset; %account for start position
    
    % Find time for current to decay to 2/e of the peak or 75ms
    % after the peak, whichever comes first. Use that for fitting
    % the single exponential. Fit the unsmoothed mean trace.
    [~,fitInd] = min(abs(traceData(pkLoc:75*sf+pkLoc)...
        - (traceData(pkLoc)/(2*exp(1)))));
    
    fitTime = fitInd/sf; % seconds
    tVec = 0:1/sf:fitTime;
    
    %traceData is already column in this vs. IdAnalysis input to findMRCs
    %TODO: On findMRCs merge, check for column and transpose if necessary
    %(also needed for ISIAnalysis to work properly)
    pkFit = fit(tVec',traceData(pkLoc:pkLoc+fitInd)','exp1');
    tau = -1/pkFit.b;
    
    
else
    pk = 0;
    pkLoc = nan;
    
    tau = nan;
    pkFit = 0;
end

varargout{1} = tau;
varargout{2} = pkFit;

end