% findMRCs_sweeps.m
%
% OUTPUT:
% cellPeaks: 
%[sortParam size vel pkLoc pk pkThresh tauAct tauDecay tPk distance nReps intPeak]%

function [cellPeaks, cellFit] = findMRCs_sweeps(stimParams, paramTraces, dataType, varargin)
p = inputParser;

p.addRequired('stimParams');
p.addRequired('paramTraces');
p.addRequired('dataType', @(x) ischar(x));

p.addParameter('tauType','fit', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('integrateCurrent', 0); %1 to make column #8 with area under the curve
p.addParameter('combineSweeps',0); %1 to combine sweeps with same sorting parameter but 
%different sampling frequencies. This will analyze each set separately and
%give a weighted mean of the stats based on nReps for each condition. Only
%use if you're reporting that you've combined two different protocols.
p.addParameter('combineLocs',0); %1 to combine sweeps with the same sorting parameter but
%different stim start locations. Only use this if you know the other stim
%parameters (e.g., stim duration/sweep interval) are the same and only the timing is shifted, 
%or if you're reporting that you've combined two different protocols.
p.parse(stimParams, paramTraces, dataType, varargin{:});

tauType = p.Results.tauType;
integrateFlag = p.Results.integrateCurrent;
combineFlag = p.Results.combineSweeps;

% smoothWindow = sf; % n timepoints for moving average window for findPeaks, as factor of sampling freq (kHz)
threshTime = 100; % use first n ms of trace for setting noise threshold
nParams = size(paramTraces,1);


% Check if different sampling frequencies were used in protocols with the
% same sorting parameter, and give separate 
[stimParams, sfSortIdx, ~, sfStart, sfEnd] = sortRowsTol(stimParams,0,10);
paramTraces = paramTraces(sfSortIdx,:);

meansBySf = nan(length(sfStarts),size(paramTraces,2));

for iSf = 1:length(sfStart)
    theseStim = stimParams(sfStart(iSf):sfEnd(iSf),:);
    theseTraces = paramTraces(sfStart(iSf):sfEnd(iSf),:);
    
    
    sf = theseStim(1,10);
    smoothWindow = sf;
    % Number of timepoints to skip after stimulus onset to avoid the stimulus
    % artifact in peak-finding (dependent on sampling frequency in kHz).
    % Determined empirically.
    artifactOffset = sf*1.2; % 1.2ms
    
    %Check if different start locations were used within a given sf set. If
    %we're not combining locations, analyze and give separate outputs.
    [theseStim,locSortIdx,~,locStart, locEnd] = sortRowsTol(theseStim,0,1);
    theseTraces = theseTraces(locSortIdx,:);
    
    if length(locStart)>1 && ~combineLocs
       
        %NEXT: fill this in with separated analysis, maybe stack the stats
        %and mean traces along dim 3, and if necessary use reshape to interleave and
        %combine all sf/loc combinations so all the stats are in one table
        %but as separate rows?
    end
    
    
    if sfEnd - sfStart > 0
        meansBySf(iSf,:) = nanmean(theseTraces,1);
    else
        meansBySf(iSf,:) = theseTraces;
    end
    

    % Smooth data with a moving average for peak finding and flip if current
    % trace but not if voltage trace (for peak finding)
    % TODO: Add a flag to make this usable in reversal potential peak finding,
    % or use absolute value.
    switch dataType
        case 'A'
            smooMean = arrayfun(@(x) -smooth(meansBySf(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
        case 'V'
            smooMean = arrayfun(@(x) smooth(meansBySf(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
    end
    
    smooMean = [smooMean{:}]';
    
    
    cellPeaks = [];
    
    % Set threshold based on noise of the first 100ms of the trace
    % (i.e., size of signal needed to be seen above that noise)
    
    % NEXT: Find only the stimuli near stimWindow for the given parameter value
    %  (e.g., for this size or velocity).
    pkThresh(iSf) = 1.5*thselect(smooMean(iSf,1:threshTime*sf),'rigrsure');
    
    % Find MRC peaks if they exist, otherwise set peak amplitude as 0.
    % Calculate decay constant tau based on single exponent fit.
    
    % sf*2.4 factor helps avoid stimulus artifact in peak finding
    % for sf = 5kHz, skips first 12 timepoints after stim.
    % NEXT: Redo this look with a cell where stim was 2.5kHz filtered and use
    % that buffer instead, bc more cells have it. Or set timepoints based on
    % stim filter freqz
    
    stimStart = stimParams(iSf,1);
    stimEnd = stimParams(iSf,2);
    
    % find peaks within stimulus (up to 20ms after end of stimulus)
    [peaks, peakLocs] = findpeaks(abs(smooMean(iSf,stimStart+artifactOffset:stimEnd+(sf*20))),...
        'minpeakheight',pkThresh(iSf));
    if ~isempty(peaks)
        
        pk = max(peaks); %take the largest peak
        
        %TODO: Use grpdelay to adjust for filter delay? If there is one, this
        %might also help make the tau calculation more correct. (half-max
        %timepoints for smooMean are the same as for meansBySf though,
        %because it's a moving average filter?)
        
        % smoothDelay = floor((smoothWindow-1)/2); %using floor for round number timepoints
        
        peakLocs = peakLocs(peaks==pk);
        pkLoc = peakLocs(1) + stimParams(iSf,1)+artifactOffset; %account for start position
        
        switch tauType
            case 'fit'
                
                % Find time for current to decay to 2/e of the peak or 75ms
                % after the peak, whichever comes first. Use that for fitting
                % the single exponential. Fit the unsmoothed mean trace.
                
                [~,fitInd] = min(abs(meansBySf(iSf,pkLoc:75*sf+pkLoc)...
                    - (meansBySf(iSf,pkLoc)/(2*exp(1)))));
                
                fitTime = fitInd/sf; % seconds
                tVec = 0:1/sf:fitTime;
                
                pkFit = fit(tVec',meansBySf(iSf,pkLoc:pkLoc+fitInd)','exp1');
                
                tauDecay = -1/pkFit.b;
                
                cellFit{iSf} = pkFit; %fit object
                tauAct = nan;
                
            case 'thalfmax' %use the timepoint of half-maximal current instead of exp fit
                halfpk = pk/2;
                halfLocs = find(smooMean(iSf,stimStart:stimEnd+(sf*200))>=halfpk);
                
                tauAct = (halfLocs(1)-1)/sf; %ms
                
                decayHalfLocs = find(smooMean(iSf,pkLoc:pkLoc+(sf*100))<=halfpk);
                % decayExpLocs = find(smooMean(iSf,pkLoc:pkLoc+(sf*100))<= pk*exp(1));
                try tauDecay = (decayHalfLocs(1)-1)/sf;
                catch
                    tauDecay = NaN;
                end
                
                % cellFit(iSf,1) = tau;
                % cellFit(iSf,2) = tauDecay;
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
            try intPeak = trapz(meansBySf(iSf,stimStart:stimEnd+(300*sf))/sf);
            catch
                intPeak = trapz(meansBySf(iSf,stimStart:end)/sf);
            end
            
            %intPeakArtifact = trapz(meansBySf(iSf,stimStart+artifactOffset:stimEnd+(sf*1E3/50))/sf);
            %intPeakHalf = trapz(meansBySf(iSf,halfLocs(1)-1:decayHalfLocs(1)-1)/sf);
            
            cellPeaks(iSf,12) = intPeak;
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
    
    cellPeaks(iSf,4) = pkLoc;
    cellPeaks(iSf,6) = pk;
    cellPeaks(iSf,8) = tauAct;
    cellPeaks(iSf,9) = tauDecay;
    cellPeaks(iSf,10) = tPk;
    
end

cellPeaks(:,7) = pkThresh;
cellPeaks(:,1) = stimParams(:,3); % stim size, pos, velocity, or interval - the sorting parameter
cellPeaks(:,2:4) = stimParams(:,4:6); % stim size, position and velocity
cellPeaks(:,12) = stimParams(:,7); %nReps
cellPeaks(:,11) = stimParams(:,8); %stim distance (0 if not entered)
end
