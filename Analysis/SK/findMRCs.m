% findMRCs.m
%
%

function [pk, pkLoc, tau, pkThresh, pkFit] = findMRCs(stimStart, traceData)

smoothWindow = 5; % n timepoints for moving average window for findPeaks
threshTime = 100; % use first n ms of trace for setting noise threshold

% Smooth data with a moving average for peak finding
smooMean = -smooth(traceData',smoothWindow,'moving');

% Set threshold based on noise of the first 100ms of the trace
% (i.e., size of signal needed to be seen above that noise)
pkThresh = 1.5*thselect(smooMean(1:threshTime*sf),'rigrsure');


% Find MRC peaks if they exist, otherwise set peak amplitude as 0. 
% Calculate decay constant tau based on single exponent fit.

[peaks, peakLocs] = findpeaks(smooMean(stimStart-100:stimStart+100),...
    'minpeakheight',pkThresh);
if ~isempty(peaks)
    
    pk = max(peaks)*1E12;
    peakLocs = peakLocs(peaks==max(peaks));
    pkLoc = peakLocs(1) + stimStart-100; %account for start position
    
    % Find time for current to decay to 2/e of the peak or 75ms
    % after the peak, whichever comes first. Use that for fitting
    % the single exponential. Fit the unsmoothed mean trace.
    [~,fitInd] = min(abs(traceData(pkLoc:75*sf+pkLoc)...
        - (traceData(pkLoc)/(2*exp(1)))));
    
    fitTime = fitInd/sf; % seconds
    tVec = 0:1/sf:fitTime;
    
    pkFit = fit(tVec',traceData(pkLoc:pkLoc+fitInd)','exp1');
    tau = -1/pkFit.b;
    
else
    pk = 0;
    pkLoc = nan;
end

end