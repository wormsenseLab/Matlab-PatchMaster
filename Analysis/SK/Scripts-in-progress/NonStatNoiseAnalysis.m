% NonStatNoiseAnalysis.m
% 
% 
% 
% Created by Sammy Katta on 5th Jan 2018.

function nonStatOutput = NonStatNoiseAnalysis(ephysData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));
%sliding window size (n sweeps) for averaging
p.addParameter('windowSize', 8, @(x) isnumeric(x) && mod(x,2)==0);
% length of time (ms) per stimulus over which the responses will be averaged
p.addParameter('responseTime', 200, @(x) isnumeric(x) && x>0); 
p.addParameter('excludeVariableStim', 'false', @(x) islogical);

p.parse(ephysData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;
averagingWindow = p.Results.windowSize; 
responseTime = p.Results.responseTime;
excludeStimVar = p.Results.excludeVariableStim;

lengthTol = 5;

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

% Allow for pre-filtered subsets of "allCells" (from FilterRecordings) to 
% be used - otherwise take names from imported metadata sheet.
if isempty(allCells)
    allCells = unique(mechTracePicks(:,1));
end

stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um


for iCell = 1:length(allCells)
    cellName = allCells(iCell);
   
    
    allStim = cell(0);
    allLeakSub = cell(0);
    
    % Double check that the list of series given matches the indices of the
    % protocol names specified in the input.
    
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,protList,'MatchType',matchType);
    
    nSeries = length(allSeries);
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
    if nSeries == 0 || isempty(pickedSeries)
        continue
    end
    
    
    for iSeries = 1:nSeries
        thisSeries = allSeries(iSeries);
        
        % Carry out analysis if this series is on the list
        try pickedTraces = pickedSeries{[pickedSeries{:,1}]==thisSeries,2};
        catch
            continue % if it's not on the list, go on to next series in for loop
        end
        
        probeI = ephysData.(cellName).data{1,thisSeries}(:,pickedTraces);
        stimComI = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces); %in V, not um
        protName = ephysData.(cellName).protocols{thisSeries};
        
        %TODO: add exclusion criteria
        % if flag is set, use photodiode trace to exclude sweeps? (before
        % separating into windows)? recordings? where the coefficient of
        % variance of the PD trace over the stimulus region is greater than
        % CV for the MRC trace.
        if excludeStimVar
            pdV = ephysData.(cellName).data{3,thisSeries}(:,pickedTraces); %in V, not um
        end
        
        % sampling frequency in kHz
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        nSweeps = size(stimComI,2);

                
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtractCell = num2cell(leakSubtract',2);
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);
        
        % quick check to make sure difference in newStepFind's ability to
        % get step/ramp timing across sweeps. For now, skip and notify.
        %TODO: change this to uniquetol to set a tolerance, especially if
        %this is to be used with slower ramps. Exact timing shouldn't
        %matter since you're setting rough time boundaries, but it should
        %be similar enough between sweeps that you're sure the stimuli are
        %as close to aligned as possible. (And then figure out how to pick
        %one stim for the whole averaging window).
        if length(uniquetol(seriesStimuli(:,1),lengthTol)) > max(seriesStimuli(:,7))
            fprintf('Stimuli start point non-identical for %s, series %d. Skipped.\n', cellName, allSeries(iSeries));
            continue
        end

        stimLoc = unique(seriesStimuli(:,1:2),'rows');
        nStim = size(stimLoc,1);
        
        %TODO: add PD-variance-based exclusion filtering here
        
        %NEXT: use seriesStimuli location/sweep number to set time
        %boundaries for where to look at the response
        %NEXT: save variance (subtracted trace) and mean current for each
        %window, plot vs. each other for all sweeps
        
        windowMeans = cell(0);
        meanSubtract = cell(0);
        
        for iSweep = 1:averagingWindow/2:nSweeps-averagingWindow+1
            try theseSweeps = probeI(:,iSweep:iSweep+averagingWindow-1);
                windowMeans{iSweep} = mean(theseSweeps,2);
            catch
                continue
            end
            
            meanSubtract{iSweep} = theseSweeps-repmat(windowMeans{iSweep},1,averagingWindow);
            
        end
        
        
        windowMeans = [windowMeans{~cellfun('isempty',windowMeans)}];       
        
        % Subtract the mean response within each window to leave only the 
        % variation around the large current (ideally these represent 
        % single-channel fluctuations from the channels of interest, which 
        % have higher open probability during the evoked response).
        meanSubtract = meanSubtract(~cellfun('isempty',meanSubtract));
        windowVars = cellfun(@(x) var(x,0,2), meanSubtract, 'un', 0);
        windowVars = [windowVars{:}]; % variance of mean-subtracted sweeps
                
        % Also calculate a non-sliding mean and variance over all sweeps
        % for comparison.
        totalMean = mean(probeI,2);
        totalSubtract = probeI - repmat(totalMean,1,size(probeI,2));
        totalVar = var(totalSubtract,0,2);
        
        % Save everything to output struct
        nonStatOutput.(cellName)(iSeries).protocol = protName;
        nonStatOutput.(cellName)(iSeries).slidingMean = windowMeans;
        nonStatOutput.(cellName)(iSeries).slidingVar = windowVars;
        nonStatOutput.(cellName)(iSeries).sweepsPerWindow = averagingWindow;
        nonStatOutput.(cellName)(iSeries).totalMean = totalMean;
        nonStatOutput.(cellName)(iSeries).totalVar = totalVar;
        
    end
    
    
    
end

end