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

% Set 'roundedTo' to 0 for no rounding.
% Optional input: include baseline values or not

function seriesStimuli = newStepFind(nSweeps, stimData, sf, varargin)

% INPUTS
% Validate inputs, allow optional inputs for more complicated usage of the
% function, and set defaults for optional inputs not passed in.
p = inputParser;

p.addRequired('nSweeps', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('stimData', @(x) isnumeric(x));
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);

p.addOptional('thresholdTime', 10, @(x) isnumeric(x) && isscalar(x) && x>0)
p.addOptional('approvedTraces', 1:nSweeps, @(x) isnumeric(x));

p.addParameter('nStim', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('minStimInterval',300, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('roundedTo',0, @(x) isnumeric(x));
p.addParameter('endTime',0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('scaleFactor',1, @(x) isnumeric(x) && isscalar(x)); %in V/um
p.addParameter('smoothWindow',0, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('threshFactor',5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(nSweeps, stimData, sf, varargin{:});

threshTime = p.Results.thresholdTime; %ms
approvedTraces = p.Results.approvedTraces;
nStim = p.Results.nStim;
minStimInterval = p.Results.minStimInterval;
roundedTo = p.Results.roundedTo;
endTime = p.Results.endTime;
scaleFactor = p.Results.scaleFactor;
smoothWindow = p.Results.smoothWindow;
if smoothWindow == 0
    smoothWindow = sf; %if smoothWindow is not set, set to 1ms worth of samples
end
threshFactor = p.Results.threshFactor;
%
% cellName = 'FAT104';
% series = 22;
% stimData = ephysData.(cellName).data{2,series};

nSweeps = size(stimData,2);
% sf = ephysData.(cellName).samplingFreq{series}/1000; %kHz
si = 1/sf; %ms
extStimFilterFreq = 2.5; %kHz
vToDispFactor = 1/scaleFactor;
roundedTo = 1/roundedTo;

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
    plot(sweepDiffSmooth./tDiff)
%     plot(sweepDiff./tDiff)
    
    %Next: modify threshold manually for each run to find steps...
    
    % Select threshold based on unsmoothed trace, then use threshold on
    % abs(smoothed trace) to capture peaks and plateaus completely.
    % Empirically, 5X threshold does best job of capturing even the
    % slowest ramps (smallest change in dx/dt), while being far above the
    % noise.
    stThresh = thselect(sweepDiff(1:threshTime*sf),'rigrsure');
    
    stLoc = find(abs(sweepDiffSmooth./tDiff)>threshFactor*stThresh);
    
    % Find lengths of plateaus/peaks above threshold. Drop any with length
    % < stimWindow (those will be artifacts at beginning or end of trace).
    % Find indices for start and end of stimuli.
    a=diff(stLoc);
    b=diff(stLoc(1:2:end)); % instead of taking diff for every point, take it across every other point (helps avoid dropping the slow ramps and shouldn't affect the fast ones)
    c=diff(stLoc(2:2:end));
    d(1:2:floor(length(a)/2)*2,1) = b;
    d(2:2:floor((length(a)-1)/2)*2) = c;
%     stLocEndIdx=find([a;inf]>3);  % 1 for anything that's not a consecutive step
    stLocEndIdx=find([d;inf]>3);
     
    % 3 allows for a dropped timepoint/little
    % bit of noise without dropping the whole run.
    stLength=diff([0;stLocEndIdx]); %find length of the sequences (not accounting for filter delay)
    stLocEndIdx = stLocEndIdx(stLength>=smoothWindow-1); %drop if too short to be a step
    stLength = stLength(stLength>=smoothWindow-1);
    stLocStartIdx = stLocEndIdx-stLength+1;
        
    
    % Account for filter delay: (N samples - 1)/2. Add timepoints to start
    % and subtract from end bc I'm looking for the first threshold
    % crossing and not a peak, and assuming that the peak/plateau is high
    % enough that averaging in a single point of it will raise the smoothed
    % signal above the threshold (pretty much true). Opposite holds for
    % the end timepoint, where it should be the last in the smoothWindow.
    smoothDelay = floor((smoothWindow-1)/2); %using floor for round number timepoints

    %Drop "step" at end of trace, before signal goes to zero for
    %protocols with varying time. Look for at least 5 zeros in a row within
    %the 10 points following the "stimulus" end location.
    
    try stLocZeroTest = sweepDiffSmooth(stLoc(stLocEndIdx(end)):stLoc(stLocEndIdx(end))+10);
    catch
        stLocZeroTest = sweepDiffSmooth(stLoc(stLocEndIdx(end)):end);
    end
    if sum(stLocZeroTest==0)>=5
        stLocEndIdx=stLocEndIdx(1:end-1);
        stLocStartIdx=stLocStartIdx(1:end-1);
    end

    
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
    %TODO: Consider flag to use PD channel instead of stim channel. You
    %can't just pass in one because you need command channel with low noise
    %to find the timepoints.
    
    isStep = stLengthActual<=stepMaxLength;
    
    % If step, find size by subtracting the command value from the 10ms
    % immediately following the step vs. the 10ms prior (or threshTime ms).
    stepSize = arrayfun(@(x,y) mean(stimSweep(x+1:x+sf*threshTime)) - ...
        mean(stimSweep(y-sf*threshTime:y-1)), stLocEnd, stLocStart,'UniformOutput',false);
    stepSize = [stepSize{:}]'.*isStep*vToDispFactor;
    
    % If ramp, find size using the mean at the three timepoints surrounding
    % ramp end minus the 10ms (or threshTime ms) prior to start of ramp.
    % (Step will probably not be immediately followed by anything else
    % within 10ms, but ramp might have triangle stimulus.)
    rampSize = arrayfun(@(x,y) mean(stimSweep(x:x+2)) - ...
        mean(stimSweep(y-sf*threshTime:y-1)), stLocEnd, stLocStart,'UniformOutput',false);
    rampSize = [rampSize{:}]'.*~isStep*vToDispFactor;
    % Convert both from voltage to displacement size based on conversion
    % factor for that particular piezo stack (in um). Get rid of any step
    % data in the ramps and vice versa.
    
    % Combination of step and ramp sizes.
    stSize = stepSize + rampSize;
    cumulStSize = cumsum(stSize);
    % Calculate stimulus velocity, in um/s.
    stSpeed = abs(stSize) ./ (stLengthActual/sf/1000);
    
    % By default, round step size to nearest 0.1um to drop noise and allow
    % grouping of step sizes. roundedTo is an optional input that can be set
    % larger or smaller depending on the range of step values used.
    if ~isinf(roundedTo)
        stSize = round(stSize*roundedTo)/roundedTo;
        cumulStSize = round(cumulStSize*roundedTo)/roundedTo;
    end
    
    % Save into array that will be concatenated for the series.
    % [startTimepoint  stopTimepoint  +/-stepSize  cumulativeStepSize  stimVelocity  sweep# stim#]
    sweepStimuli(:,1) = stLocStart;
    sweepStimuli(:,2) = stLocEnd;
    sweepStimuli(:,3) = stSize;
    sweepStimuli(:,4) = cumulStSize;
    sweepStimuli(:,5) = stSpeed;
    sweepStimuli(:,6) = repmat(iSweep,length(stLocStartIdx),1);
    sweepStimuli(:,7) = (1:length(stLocStartIdx))';
    
    seriesStimuli = [seriesStimuli; sweepStimuli];
    clear sweepStimuli d
end

end



