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

p.addParameter('sortStimBy', 'num', @(x) sum(strcmp(x,{'num','time'})));
p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));
p.addParameter('windowSize', 8, @(x) isnumeric(x) && mod(x,2)==0); %sliding window size (n sweeps) for averaging
p.addParameter('responseTime', 200, @(x) isnumeric(x) && x>0); % length of time(ms) after stimulus to be included in the noise analysis for that stimulus (default 15ms before, 200ms after)
p.addParameter('excludeVariableStim', 'false', @(x) islogical);
p.addParameter('sortSweepsBy',{'magnitude','magnitude','magnitude','magnitude'}, @(x) iscell(x));
p.addParameter('recParameters',cell(0), @(x) iscell(x));

p.parse(ephysData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;
averagingWindow = p.Results.windowSize; 
postTime = p.Results.responseTime;
excludeStimVar = p.Results.excludeVariableStim;
sortStimBy = p.Results.sortStimBy;
ephysMetaData = p.Results.recParameters; % Use this to get stim filt freq and distance

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

% pull out stimulus distance and external filter frequency for the given
% recordings
allCellInd = cellfun(@(x) find(strcmp(ephysMetaData(:,1),x)),allCells,'un',0)';
allCellInd = [allCellInd{:}];

filterHeaders{1} = 'cellStimDistUm';
filterHeaders{2} = 'stimFilterFrequencykHz';
for i=1:length(filterHeaders)
paramInd(i) = find(strncmpi(regexprep(ephysMetaData(1,:), '[^a-zA-Z0-9]', ''),...
    regexprep(filterHeaders{i}, '[^a-zA-Z0-9]', ''),length(filterHeaders{i})));
end
cellDist = ephysMetaData(allCellInd,paramInd(1));
extFilterFreq = ephysMetaData(allCellInd,paramInd(2));



baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
preTime = 15;
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um


for iCell = 1:length(allCells)
    cellName = allCells(iCell);
   
       
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
        allStim = cell(0);

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
        try thisDist = round(cellDist{iCell},1);
        catch
        end
                
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);
        
        
        % add series number to seriesStimuli for referencing back from
        % analyzed data to a particular trace
        seriesStimuli (:,8) = repmat(thisSeries,size(seriesStimuli,1),1);
        
        % quick check to make sure difference in newStepFind's ability to
        % get step/ramp timing across sweeps. For now, skip and notify.
        if length(uniquetol(seriesStimuli(:,1),lengthTol,'DataScale', 1)) > max(seriesStimuli(:,7))
            fprintf('Stimuli start point non-identical for %s, series %d. Skipped.\n', cellName, allSeries(iSeries));
            continue
        end

        % write in stim distance
        if exist('thisDist','var') && ~isempty(thisDist) && ~isnan(thisDist)
            seriesStimuli(:,9) = repmat(thisDist, size(seriesStimuli,1),1);
        else
            % if no stim distance is included, must be zeros instead of
            % NaNs because unique will consider every NaN unique when
            % sorting by stim for zeroFill later.
            seriesStimuli(:,9) = zeros(size(seriesStimuli,1),1);
        end
        
        
        % depending on input param, use either stimulus number or stimulus
        % timepoint to group stimuli across sweeps (e.g., on and off)

        switch sortStimBy
            case 'num'
                paramCol = 7; %set which column of seriesStimuli to look in
                tol = 0; %set tolerance for separating by parameter
            case 'time'
                paramCol = 1;
                tol = 5;  
        end
        stimNums = uniquetol(seriesStimuli(:,paramCol),tol,'DataScale',1);
        nStims = length(stimNums);
              
        % separate found stim parameters into groups and concatenate to
        % previous stimuli from the same series
        stimByNum = cell(nStims,1);
        for iStim = 1:nStims
            stimByNum{iStim} = ... % find param matches within tolerance and assign into stimByNum
                seriesStimuli(ismembertol(seriesStimuli(:,paramCol),stimNums(iStim),tol,'DataScale',1),:);
            try allStim{iStim,1}; %if allStim doesn't exist, initialize it
            catch
                allStim{iStim,1} = [];
            end
            
            if ~length(allStim{iStim})==0; %if values exist for that stimulus in allStim, 
                           %match param w/ tolerance to the matching param
                           %location in allStim
                [~,whichStim]=ismembertol(stimNums,cellfun(@(x) x(1,paramCol), allStim),1,'DataScale',tol);
                allStim{whichStim(iStim),1} = [allStim{whichStim(iStim),1};stimByNum{iStim}];
                
            else %if allStim exists but that stimulus is empty, start it up
                allStim{iStim,1} = [allStim{iStim,1};stimByNum{iStim}];
            end
            
        end
       
        
        % take the mean and variance of sliding windows of sweeps within
        % the time period indicated by the start point for each stimulus up
        % to the end time (based on the input response time, default 200ms)
        for iStim = 1:nStims
            
            stimLoc = uniquetol(allStim{iStim}(:,1:2),lengthTol,'ByRows',true);
            stimSize = round(mean(allStim{iStim}(:,3)));
            stimPos = round(mean(allStim{iStim}(:,4)));
            stimVel = roundVel(mean(allStim{iStim}(:,5)));
            responseTime = [stimLoc(1)-preTime*sf:stimLoc(2)+postTime*sf];
            %NEXT: use seriesStimuli location/sweep number to set time
            %boundaries for where to look at the response
            
            windowMeans = cell(0);
            meanSubtract = cell(0);
            
            for iSweep = 1:averagingWindow/2:nSweeps-averagingWindow+1
                try theseSweeps = leakSubtract(responseTime,iSweep:iSweep+averagingWindow-1);
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
            totalMean = mean(leakSubtract(responseTime,:),2);
            totalSubtract = leakSubtract(responseTime,:) - repmat(totalMean,1,size(leakSubtract(responseTime,:),2));
            totalVar = var(totalSubtract,0,2);
            
            whichRow = iSeries+((iStim-1)*nSeries);
            
            % Save everything to output struct
            nonStatOutput.(cellName)(whichRow).protocol = protName;
            nonStatOutput.(cellName)(whichRow).stimNum = iStim;
            nonStatOutput.(cellName)(whichRow).size = stimSize;
            nonStatOutput.(cellName)(whichRow).position = stimPos;
            nonStatOutput.(cellName)(whichRow).velocity = stimVel;
            nonStatOutput.(cellName)(whichRow).slidingMean = windowMeans;
            nonStatOutput.(cellName)(whichRow).slidingVar = windowVars;
            nonStatOutput.(cellName)(whichRow).sweepsPerWindow = averagingWindow;
            nonStatOutput.(cellName)(whichRow).totalMean = totalMean;
            nonStatOutput.(cellName)(whichRow).totalVar = totalVar;
            nonStatOutput.(cellName)(whichRow).distance = thisDist;

            
        end
        
        clear thisDist;
    end
    
    
end

end