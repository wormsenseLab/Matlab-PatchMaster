%newStepFind.m

% Optional parameter input: 'extStimFilterFreq' in kHz, auto-set to 0, but remember to
% call it with 2.5kHz. (Include this in ephysRecordingBase or somehow read
% it into ephysData). Currently 2.5kHz for all the recent recordings, with
% the LPF-8 in line before the Crawford amplifier.
% 
% Optional param input: voltage to displacement conversion factor (default
% V/0.408). And flag to use PD signal for "actual size/rate" finding, add
% that column as separate varargout (can't just feed as input bc it's not
% clean enough for ramp detection, but once step is found, can use
% timepoints to find calibrated actual displacement).

cellName = 'FAT104';
series = 22;

stimData = ephysData.(cellName).data{2,series};
nSweeps = size(stimData,2);
sf = ephysData.(cellName).samplingFreq{series}/1000; %kHz
si = 1/sf; %ms
threshTime = 10; %ms
smoothWindow = sf;
extStimFilterFreq = 2.5; %kHz
vToDispFactor = 1/0.408;

% Find size of step (in um), as well as start and end indices.

% Initialize array for storing stim paramters. sweepsize is positive
% for pushing into the worm, negative for pulling back. Sweep# is within
% the current series, stim# is nth up or down step within the sweep.
% Organization:  [startTime  stopTime  +/-stepSize  sweep# stim#] 
seriesStimuli = [];

for iSweep = 1:nSweeps
    stimSweep = stimData(:,iSweep);
    tVec = (0:si:length(stimSweep)*si-si)';
    
    
    % Take first derivative of stimulus trace to find regions where there
    % is a ramp or step (slope > 0).
    sweepDiff = diff(stimSweep);
    sweepDiffSmooth = smooth(sweepDiff, smoothWindow, 'moving');
    tDiff = diff(tVec);
    % If you want to plot, use diff(y)./diff(t)
    
    % Select threshold based on unsmoothed trace, then use threshold on
    % abs(smoothed trace) to capture peaks and plateaus completely.
    % Empirically, 5X threshold does best job of capturing even the
    % slowest ramps (smallest change in dx/dt), while being far above the
    % noise.
    stThresh = thselect(sweepDiff(1:threshTime*sf),'rigrsure');
    
    stLoc = find(abs(sweepDiffSmooth./tDiff)>5*stThresh); 
 
    % Find lengths of plateaus/peaks above threshold. Drop any with length
    % < stimWindow (those will be artifacts at beginning or end of trace).
    % Find indices for start and end of stimuli.
    a=diff(stLoc);
    stLocEndIdx=find([a;inf]>3);  % 1 for anything that's not a consecutive step
                               % 3 allows for a dropped timepoint/little
                               % bit of noise without dropping the whole run.
    stLength=diff([0;stLocEndIdx]); %find length of the sequences (not accounting for filter delay)
    stLocEndIdx = stLocEndIdx(stLength>=smoothWindow); %drop if too short to be a step
    stLength = stLength(stLength>=smoothWindow);
    stLocStartIdx = stLocEndIdx-stLength+1;
    
    % Account for filter delay: (N samples - 1)/2. Add timepoints to start
    % and subtract from end bc I'm looking for the first threshold
    % crossing and not a peak, and assuming that the peak/plateau is high
    % enough that averaging in a single point of it will raise the smoothed
    % signal above the threshold (pretty much true). Opposite holds for
    % the end timepoint, where it should be the last in the smoothWindow.
    smoothDelay = floor((smoothWindow-1)/2); %using floor for round number timepoints

    % Find actual timepoints where stimuli start/end, and stim length,
    % corrected for the filter delay.
    stLocStart = stLoc(stLocStartIdx)+smoothDelay;
    stLocEnd = stLoc(stLocEndIdx)-smoothDelay;
    stLengthActual = stLocEnd-stLocStart;
    
    %TODO: Decide whether stim is step or ramp. If ramp, size is the size
    %at the peak (can't just take size after, because what if it's a
    %triangle instead of ramp-and-hold?). If step, it's average of last x
    %ms. Include both size and rate for both ramps and steps. (For steps,
    %this should end up being the extStimFilterFreq if it exists, or
    %roughly samplingFreq if not.)
    if extStimFilterFreq < sf %if external stim filter frequency 
                              %(e.g., LPF-8 set to 2.5kHz) is less than
                              %sampling freq, ext freq will determine how
                              %many timepoints a "step" will take.
        stepMaxLength = sf/extStimFilterFreq+1;
    else
        stepMaxLength = 2; 
    end
    
    %TODO: Test stepMaxLength w one of the older traces w no filter, to see
    %if max 2 timepoints is actually what will work.
    
    isStep = stLengthActual<=stepMaxLength;
    
    stepSize = arrayfun(@(x,y) mean(stimSweep(x+1:x+sf*threshTime)) - ...
        mean(stimSweep(y-sf*threshTime:y-1)), stLocEnd, stLocStart,'UniformOutput',false);
    stepSize = [stepSize{:}]'.*isStep*vToDispFactor;
    
    rampSize = arrayfun(@(x,y) mean(stimSweep(x:x+2)) - ...
        mean(stimSweep(y-sf*threshTime:y-1)), stLocEnd, stLocStart,'UniformOutput',false);
    rampSize = [rampSize{:}]'.*~isStep*vToDispFactor;
    
    stVel = 
    %bsxfun or just do both stepSize and rampSize calculations, multiply
    %each by isStep or ~isStep, then add (non-applicable values will be
    %zeros).
      
    % For step: take difference of mean value of stimSweep 10ms before and 
    % after stStart/End to get stepSize.
    
    % By default, round step size to nearest 0.1 to drop noise and allow
    % grouping of step sizes. roundedTo is an optional input that can be set
    % larger or smaller depending on the range of step values used.
    roundedTo = 1/roundedTo;
    allSizes = round(stepSize*roundedTo)/roundedTo;
    
    
    
    % Save into array that will be concatenated for the series.
    % [startTimepoint  stopTimepoint  +/-stepSize  rampRate  sweep# stim#] 
    sweepStimuli(:,1) = stLocStart;
    sweepStimuli(:,2) = stLocEnd;
    sweepStimuli(:,3) = 
    sweepStimuli(:,4) = 
    sweepStimuli(:,5) = repmat(iSweep,length(stLocStartIdx),1);
    sweepStimuli(:,6) = (1:length(stLocStartIdx))';
    
    %TODO: Drop "step" at end of trace, before signal goes to zero for
    %protocols with varying time. Within 5ms of start/end, or within
    %1.5*smoothWindow of a 5ms+ long true zero plateau. 
    

                
end
  
    
    seriesStimuli = [seriesStimuli; sweepStimuli];
