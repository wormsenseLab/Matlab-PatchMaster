% findMRCs.m
%
% INPUT:
% stimParams
% [stimStart stimEnd sortParam size pos vel nReps dist]
% 
% OUTPUT:
% cellPeaks: 
%[sortParam size pos vel pkLoc pk pkThresh tauAct tauDecay tPk intPeak distance nReps]%
% 
% TODO: consider when to use smooMean vs. meanTraces for actual peak
% current values (also in finding the half-max tau?)

function [cellPeaks, cellFit] = findMRCs(stimParams, meanTraces, sf, dataType, varargin)
p = inputParser;

p.addRequired('stimParams');
p.addRequired('meanTraces');
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('dataType', @(x) ischar(x));

p.addParameter('tauType','fit', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('integrateCurrent', 0); %1 to make column #8 with area under the curve
p.addParameter('combineSweeps',0); %1 to combine sweeps with same sorting parameter but 
%different stim timepoints or different sampling freqs (for the latter,
%analysis will be done separately and then stats will be combined as a
%weighted average based on number of sweeps of each condition

p.parse(stimParams, meanTraces, sf, dataType, varargin{:});

tauType = p.Results.tauType;
integrateFlag = p.Results.integrateCurrent;
combineFlag = p.Results.combineSweeps;

smoothWindow = sf; % n timepoints for moving average window for findPeaks, as factor of sampling freq (kHz)
threshTime = 100; % use first n ms of trace for setting noise threshold
nParams = size(meanTraces,1);

% Number of timepoints to skip after stimulus onset to avoid the stimulus
% artifact in peak-finding (dependent on sampling frequency in kHz).
% Determined empirically.
switch sf
    case 5
        artifactOffset = sf*1.2; % 1.2ms
    case 10
        artifactOffset = sf*1.2; %1.2ms
    otherwise
        artifactOffset = sf; % 1ms, in case of non-standard sf, play safe by not cutting off too many points
end

% Smooth data with a moving average for peak finding and flip if current
% trace but not if voltage trace (for peak finding)
% TODO: Add a flag to make this usable in reversal potential peak finding,
% or use absolute value.
switch dataType
    case 'A'
        smooMean = arrayfun(@(x) -smooth(meanTraces(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
    case 'V'
        smooMean = arrayfun(@(x) smooth(meanTraces(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
end

smooMean = [smooMean{:}]';

cellPeaks = [];

% Set threshold based on noise of the first 100ms of the trace
% (i.e., size of signal needed to be seen above that noise)

% NEXT: Find only the stimuli near stimWindow for the given parameter value
%  (e.g., for this size or velocity).
for iParam = 1:nParams
    pkThresh(iParam) = 1.5*thselect(smooMean(iParam,1:threshTime*sf),'rigrsure');
    
    % Find MRC peaks if they exist, otherwise set peak amplitude as 0.
    % Calculate decay constant tau based on single exponent fit.
    
    % sf*2.4 factor helps avoid stimulus artifact in peak finding
    % for sf = 5kHz, skips first 12 timepoints after stim.
    % NEXT: Redo this look with a cell where stim was 2.5kHz filtered and use
    % that buffer instead, bc more cells have it. Or set timepoints based on
    % stim filter freqz
    
    stimStart = stimParams(iParam,1);
    stimEnd = stimParams(iParam,2);
    
    % find peaks within stimulus (up to 20ms after end of stimulus)    
    [peaks, peakLocs] = findpeaks(abs(smooMean(iParam,stimStart+artifactOffset:stimEnd+(sf*20))),...
        'minpeakheight',pkThresh(iParam));
    if ~isempty(peaks)
        
        pk = max(peaks); %take the largest peak
        
        %TODO: Use grpdelay to adjust for filter delay? If there is one, this
        %might also help make the tau calculation more correct. (half-max
        %timepoints for smooMean are the same as for meanTraces though,
        %because it's a moving average filter?)
        
        % smoothDelay = floor((smoothWindow-1)/2); %using floor for round number timepoints
        
        peakLocs = peakLocs(peaks==pk);
        pkLoc = peakLocs(1) + stimParams(iParam,1)+artifactOffset; %account for start position
        
        switch tauType %decay tau
            case 'fit'
                
                % Find time for current to decay to 2/e of the peak or 75ms
                % after the peak, whichever comes first. Use that for fitting
                % the single exponential. Fit the unsmoothed mean trace.
                
                [~,fitInd] = min(abs(meanTraces(iParam,pkLoc:75*sf+pkLoc)...
                    - (meanTraces(iParam,pkLoc)/(2*exp(1)))));
                
                fitTime = fitInd/sf; % seconds
                tVec = 0:1/sf:fitTime;
                
                pkFit = fit(tVec',meanTraces(iParam,pkLoc:pkLoc+fitInd)','exp1');
                                
                tauDecay = -1/pkFit.b;
                
                cellFit{iParam} = pkFit; %fit object
                tauAct = nan;

            case 'thalfmax' %use the timepoint of half-maximal current instead of exp fit
                halfpk = pk/2;
                halfLocs = find(abs(smooMean(iParam,stimStart:stimEnd+(sf*200)))>=halfpk);

                tauAct = (halfLocs(1)-1)/sf; %ms
                
                decayHalfLocs = find(smooMean(iParam,pkLoc:pkLoc+(sf*100))<=halfpk);
                % decayExpLocs = find(smooMean(iParam,pkLoc:pkLoc+(sf*100))<= pk*exp(1));
                try tauDecay = (decayHalfLocs(1)-1)/sf;
                catch
                    tauDecay = NaN;
                end
                
                % cellFit(iParam,1) = tau;
                % cellFit(iParam,2) = tauDecay;                              
        end
        
        switch dataType
            case 'A'
                pk = pk*1E12; %pA
            case 'V'
                pk = pk*1E3; %mV
        end
        
        % Integrate current for total charge carried
        if integrateFlag
            % trapz uses the trapezoidal method to integrate & calculate area under
            % the curve. But it assumes unit spacing, so divide by the sampling
            % frequency to get units of seconds.            
            try intPeak = trapz(meanTraces(iParam,stimStart:stimEnd+(150*sf))/(sf*1000));
            catch
                intPeak = trapz(meanTraces(iParam,stimStart:end)/(sf*1000));
            end
                
            %intPeakArtifact = trapz(meanTraces(iParam,stimStart+artifactOffset:stimEnd+(sf*1E3/50))/sf);
            %intPeakHalf = trapz(meanTraces(iParam,halfLocs(1)-1:decayHalfLocs(1)-1)/sf);
            
            cellPeaks(iParam,11) = intPeak;
        end
        tPk = (pkLoc - stimStart)/sf;
 
    else
        pk = 0;
        pkLoc = nan;
        tPk = nan;
        
        tauAct = nan;
        tauDecay = nan;
        pkFit = 0;
        
    end
    
    cellPeaks(iParam,5) = pkLoc;
    cellPeaks(iParam,6) = pk;
    cellPeaks(iParam,8) = tauAct;
    cellPeaks(iParam,9) = tauDecay;
    cellPeaks(iParam,10) = tPk;


end
cellPeaks(:,7) = pkThresh;
cellPeaks(:,1) = stimParams(:,3); % stim size, pos, velocity, or interval - the sorting parameter
cellPeaks(:,2:4) = stimParams(:,4:6); % stim size, position and velocity
cellPeaks(:,13) = stimParams(:,end-1); %nReps
cellPeaks(:,12) = stimParams(:,end); %stim distance (0 if not entered)
end
