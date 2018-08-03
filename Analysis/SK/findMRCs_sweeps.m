% findMRCs_sweeps.m
%
% OUTPUT:
% cellPeaks: 
%[sortParam size pos vel pkLoc pk pkThresh tauAct tauDecay tPk distance nReps intPeak]%

function [cellPeaks, cellFit] = findMRCs_sweeps(stimParams, paramTraces, dataType, varargin)
p = inputParser;

p.addRequired('stimParams');
p.addRequired('paramTraces');
p.addRequired('dataType', @(x) ischar(x));

p.addParameter('tauType','fit', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('integrateCurrent', 0); %1 to make column #8 with area under the curve
p.addParameter('combineSFs',0); %1 to combine sweeps with same sorting parameter but 
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
combineSFs = logical(p.Results.combineSFs);
combineLocs = logical(p.Results.combineLocs);
roundFlag = 0; % if 1, round size/pos/vel parameters for easy reading in output.

% smoothWindow = sf; % n timepoints for moving average window for findPeaks, as factor of sampling freq (kHz)
threshTime = 100; % use first n ms of trace for setting noise threshold
responseWindow = 250; % use this timewindow (ms) for response characteristics (incl. decay tau)
peakWindow = 20; % look for peaks within x ms of the stimulus end
nParams = size(paramTraces,1);


% Check if different sampling frequencies were used in protocols with the
% same sorting parameter, and give separate 
[stimParams, sfSortIdx, ~, sfStart, sfEnd] = sortRowsTol(stimParams,0,11);
paramTraces = paramTraces(sfSortIdx,:);

meansBySf = nan(length(sfStart),size(paramTraces,2));

for iSf = 1:length(sfStart)
    theseStim = stimParams(sfStart(iSf):sfEnd(iSf),:);
    theseTraces = paramTraces(sfStart(iSf):sfEnd(iSf),:);
    
    
    sf = theseStim(1,11);
    smoothWindow = sf;
    % Number of timepoints to skip after stimulus onset to avoid the stimulus
    % artifact in peak-finding (dependent on sampling frequency in kHz).
    % Determined empirically.
    artifactOffset = sf*1.2; % 1.2ms
    
    % TODO: check if artifactOffset factor is good for either 5kHz or
    % 2.5kHz
    
    
    %Check if different start locations were used within a given sf set. If
    %we're not combining locations, analyze and give separate outputs.
    [theseStim,locSortIdx,~,locStart, locEnd] = sortRowsTol(theseStim,0,1);
    theseTraces = theseTraces(locSortIdx,:);
    nLocs = length(locStart);
    
    stimStart = theseStim(locStart,1);
    stimEnd = theseStim(locStart,2);
      
    if nLocs == 1 || combineLocs
        % if combining peaks at different timepoints, or if there was only
        % one timepoint, it's simple to smooth the mean
        locStart = 1;
        locEnd = size(theseStim,1);
        
        baseMean = mean(theseTraces(:,1:threshTime*sf));
        responseMean = mean(theseTraces(:,stimStart:stimEnd+responseWindow*sf));
        
        
        % Smooth data with a moving average for peak finding and flip if current
        % trace but not if voltage trace (for peak finding).
        % TODO: Add a flag to make this usable in reversal potential peak finding,
        % or use absolute value.
        smooBase = smooth(baseMean,smoothWindow,'moving')';
        smooMean = smooth(responseMean,smoothWindow,'moving')';
        
        if dataType == 'A'
            smooBase = -smooBase;
            smooMean = -smooMean;
        end
        
    else % if there are multiple locs not being combined
        baseMean = zeros(nLocs,threshTime*sf);
        responseMean = zeros(nLocs,threshTime*sf);
        
        for iLoc = 1:nLocs
            baseMean(iLoc,:) = mean(theseTraces(locStart(iLoc):locEnd(iLoc),1:threshTime*sf));
            responseMean(iLoc,:) = mean(theseTraces(locStart(iLoc):locEnd(iLoc),stimStart:stimEnd+responseWindow*sf));            
        end
        % FIX: responseMean won't work with ramp stimuli because this will
        % be different sizes. Pad it with nans.
        
        
        % Smooth data with a moving average for peak finding and flip if current
        % trace but not if voltage trace (for peak finding)
        % Here, keep one smoothed trace for each iLoc.
        smooBase = arrayfun(@(x) smooth(baseMean(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
        smooMean = arrayfun(@(x) smooth(responseMean(x,:),smoothWindow,'moving'), 1:nParams, 'un',0)';
        
        if dataType == 'A'
            smooBase = -smooBase;
            smooMean = -smooMean; 
        end
        
        smooBase = [smooBase{:}]';
        smooMean = [smooMean{:}]';
        
    end
      
    cellPeaks = [];
    
    for iLoc = 1:nLocs
        % Set threshold based on noise of the first threshTim ms of the trace
        % (i.e., size of signal needed to be seen above that noise)
        
        pkThresh(iLoc) = 1.5*thselect(smooBase(iLoc,:),'rigrsure');
        
        
        % Find MRC peaks if they exist, otherwise set peak amplitude as 0.
        % Calculate decay constant tau based on single exponent fit.
        
        [peaks, peakLocs] = findpeaks(abs(smooMean(iLoc,1:peakWindow*sf)),'minpeakheight',pkThresh(iLoc));
        if ~isempty(peaks)
            
            pk = max(peaks); %take the largest peak
            
            %TODO: Use grpdelay to adjust for filter delay? If there is one, this
            %might also help make the tau calculation more correct. (half-max
            %timepoints for smooMean are the same as for meansBySf though,
            %because it's a moving average filter?)
            
            % smoothDelay = floor((smoothWindow-1)/2); %using floor for round number timepoints
            
            peakLocs = peakLocs(peaks==pk);
            pkLoc = peakLocs(1);
            pkLocActual = pkLoc + stimParams(iLoc,1)+artifactOffset; %account for start position
            
            switch tauType
                case 'fit'
                    
                    % Find time for current to decay to 2/e of the peak or 75ms
                    % after the peak, whichever comes first. Use that for fitting
                    % the single exponential. Fit the unsmoothed mean trace.
                    
                    [~,fitInd] = min(abs(responseMean(iLoc,pkLoc:75*sf+pkLoc)...
                        - (responseMean(iLoc,pkLoc)/(2*exp(1)))));
                    
                    fitTime = fitInd/sf; % seconds
                    tVec = 0:1/sf:fitTime;
                    
                    pkFit = fit(tVec',responseMean(iLoc,pkLoc:pkLoc+fitInd)','exp1');
                    
                    tauDecay = -1/pkFit.b;
                    
                    cellFit{iLoc} = pkFit; %fit object
                    tauAct = nan;
                    
                case 'thalfmax' %use the timepoint of half-maximal current instead of exp fit
                    halfpk = pk/2;
                    halfLocs = find(smooMean(iLoc,:)>=halfpk);
                    
                    tauAct = (halfLocs(1)-1)/sf; %ms
                    
                    decayHalfLocs = find(smooMean(iLoc,pkLoc:pkLoc+(sf*100))<=halfpk);
                    % decayExpLocs = find(smooMean(iLoc,pkLoc:pkLoc+(sf*100))<= pk*exp(1));
                    try tauDecay = (decayHalfLocs(1)-1)/sf;
                    catch
                        tauDecay = NaN;
                    end
                    
                    % cellFit(iLoc,1) = tau;
                    % cellFit(iLoc,2) = tauDecay;
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
                try intPeak = trapz(meansBySf(iLoc,stimStart:stimEnd+(responseWindow*sf))/sf);
                catch
                    intPeak = trapz(meansBySf(iLoc,stimStart:end)/sf);
                end
                
                %intPeakArtifact = trapz(meansBySf(iLoc,stimStart+artifactOffset:stimEnd+(sf*1E3/50))/sf);
                %intPeakHalf = trapz(meansBySf(iLoc,halfLocs(1)-1:decayHalfLocs(1)-1)/sf);
                
                cellPeaks(iLoc,12) = intPeak;
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
        
        cellPeaks(iLoc,5) = pkLocActual;
        cellPeaks(iLoc,6) = pk;
        cellPeaks(iLoc,8) = tauAct;
        cellPeaks(iLoc,9) = tauDecay;
        cellPeaks(iLoc,10) = tPk;
    end
       
    cellPeaks(:,7) = pkThresh;
    cellPeaks(:,1) = stimParams(locStart,3); % stim size, pos, velocity, or interval - the sorting parameter
    cellPeaks(:,12) = locEnd-locStart+1; %nReps
    cellPeaks(:,11) = stimParams(locStart,8); %stim distance (0 if not entered)
    
    if roundFlag
        cellPeaks(:,2:3) = round(stimParams(locStart,4:5),1); % stim size, position
        cellPeaks(:,4) = roundVel(stimParams(locStart,6));
    else
        
        cellPeaks(:,2:3) = stimParams(locStart,4:5); % stim size, position
        cellPeaks(:,4) = stimParams(locStart,6);      
    end

        
%     If combineSfs is set, combineLocs must also be true, because it's
%     nested inside. If that's okay, then next: take weighted mean of locs,
%     collapse into single vector, have dim 3 of cellPeaks for each sf,
%     then take weighted mean across dim 3, based on nReps and update new
%     nReps to be sum of previous ones.

%     If not combining, interleave peaks and mean traces? Should handle
%     trace sorting here, to match up to whatever the stats output looks
%     like.
    

    
    
  
    
end

end
