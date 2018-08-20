% FrequencyAnalysis.m
% 
% 
% parameters in SinePeaks{:,3}:
% [sineStart sineEnd sineFreq nReps steadyI sqOn1 sqOff1 sqOn2 sqOff2 sqOnRatio sqOffRatio]

function sinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('ephysMetaData', @(x) iscell(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));
p.addParameter('tauType','thalfmax', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('integrateCurrent',1);
p.addParameter('normToOn1',0);
p.addParameter('matchPhase',0');

p.parse(ephysData, ephysMetaData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;
tauType = p.Results.tauType;
integrateFlag = p.Results.integrateCurrent;
normalizeFlag = p.Results.normToOn1;
phaseFlag = p.Results.matchPhase;

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

% Allow for pre-filtered subsets of "allCells" (from FilterRecordings) to
% be used - otherwise take names from imported metadata sheet.
if isempty(allCells)
    allCells = unique(mechTracePicks(:,1));
end

sinePeaks = cell(length(allCells),1);

baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline for leak subtraction
steadyTime = 100; % length of time (ms) to use for calculating steady state mean at end of sine
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
whichChan = 2; %channel number for stimulus command
sortStimBy = 'time';
sortSweepsBy = {'frequency'};
stimSortOrder = [1 2];

% pull cell distances and filter freqs
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



for iCell = 1:length(allCells)
      
    % Double check that the list of series given matches the indices of the
    % protocol names specified in the input.
    
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,protList,'MatchType',matchType);
    
    nSeries = length(allSeries);
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
    if nSeries == 0 || isempty(pickedSeries)
        continue
    end
    
    % Initialize for allSeries
    allLeakSub = cell(0);
    allSweeps = cell(nSeries,2);
    allSine = cell(0);
    stimSquareSum = [];

    
    for iSeries = 1:nSeries
        thisSeries = allSeries(iSeries);
        
        % Carry out analysis if this series is on the list
        try pickedTraces = pickedSeries{[pickedSeries{:,1}]==thisSeries,2};
        catch
            continue % if it's not on the list, go on to next series in for loop
        end
        
        try thisDist = round(cellDist{iCell},1);
        catch
        end

        
        probeI = ephysData.(cellName).data{1,thisSeries}(:,pickedTraces);
        stimComI = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces); %in V, not um
        % sampling frequency in kHz
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        dataType = ephysData.(cellName).dataunit{1,thisSeries};
        nSweeps = size(stimComI,2);
        
        
        % Subtract the baseline leak current for each series
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime)';
        leakSubtractCell = num2cell(leakSubtract,2);
        
        
        clear sineParams;
        % Find location, frequency, and amplitude of sine, and locations of
        % flanking square pulses
        [sineParams, squareParams] = ...
            retrieveSineFreq(ephysData, cellName, thisSeries, whichChan);
        % Add duration column
        sineParams = [sineParams(:,1:2) sineParams(:,2)-sineParams(:,1)+1 sineParams(:,3:end)];
        
        % Add distance column
        if exist('thisDist','var') && ~isempty(thisDist) && ~isnan(thisDist)
            squareParams(:,8) = repmat(thisDist, size(squareParams,1),1);
        else
            % if no stim distance is included, must be zeros instead of
            % NaNs because unique will consider every NaN unique when
            % sorting by stim for zeroFill later.
            squareParams(:,8) = zeros(size(squareParams,1),1);
        end
        
        % In case retrieveSineFreq returned without output because stimTree
        % didn't exist, alert the user and move on to the next group.
        if ~exist('sineParams','var')
            fprintf('Please upload a stimTree for %s. (Run ImportPatchData with includeStim flag).\n',cellName);
            break
        end
        
        % For later: when using multiple sines in a given sweep
        % for iSine = 1:length(sineLoc)
        %     sineTrace{iSine} = leakSubtract(sineLoc{iSine});
        % end
        
        % For now: there's only one sine per series
        
        sineTrace = leakSubtract(sineParams(1):sineParams(2));
        
        % Take the average current of the last steadyTime ms of the
        % sine trace, regardless of sine frequency.
        steadyStateI = -mean(sineTrace(end-steadyTime*sf:end))*1e12;
        
        theseSines = [sineParams steadyStateI thisSeries];
        % Implement later: use fit to calculate time constant of decay from peak
        % response at beginning of sine. Use 5*tau or 3*tau timepoint as the
        % beginning of steady state for taking the mean. Will vary for
        % different frequencies but be a better estimate of actual steady
        % state.
        %         [~, fitStart] = (max(abs(sineTrace(5:10*sf))));
        %         fitInd = 400;
        %         fitTime = fitInd/sf; % seconds
        %         t = 0:1/sf:fitTime;
        %         
        %         capFit = ezfit(t',sineTrace(1:1+fitInd),'y(x)=c+a*exp(b*x)',[5e-12, -3000, 8e-13]);
        %         tau = -1/capFit.m(2);
        
       
        % Find MRC peaks at the on and off of square stimuli.
        theseSquares = [];
        for iStim = 1:size(squareParams,1)
            theseSquares(iStim,:) = findMRCs(squareParams(iStim,:), leakSubtract, sf, dataType, ...
                'tauType', tauType, 'integrateCurrent',integrateFlag);
        end
        theseSquares = [theseSquares repmat(sineParams(4),iStim,1) repmat(thisSeries,iStim,1)];

        
        %TODO: add to Sine output: peak at start, beginning average, tau of 
        % fit of smoothed trace
        %TODO: add to SquareSum output: interval between squares, and
        % between sweeps
        
        % Concatenate data for all series in this recording
        allLeakSub(iSeries,1) = leakSubtractCell;
        allSweeps{iSeries,1} = theseSines;
        allSweeps{iSeries,2} = theseSquares;
    end
    
    
    % The following code is adapted from IdAnalysis(), and could probably
    % be much more simplified for this purpose (the complexity would allow
    % more combos of # of sine + # of squares, but that's not implemented
    % yet).
    
    % Define how many sine stimuli there are per sweep based on either 
    % unique timepoint or duration of the segment.
    switch sortStimBy
        case 'time'
            paramCol = 1; %set which column of allSineSum to look in
            tol = 0; %set tolerance for separating by parameter
        case 'dur'
            paramCol = 3;
            tol = 0;
    end
    
    % Combine sweep data across series for both squares and sines
    
    %TODO: Make sure the rest of this code matches with the changed order
    %of square columns for sorting
    
    allSineSum = vertcat(allSweeps{:,1});
    allSquareSum = vertcat(allSweeps{:,2});
    for iStim = 1:max(allSquareSum(:,12))
        stimSquareSum(:,:,iStim) = allSquareSum(allSquareSum(:,12)==iStim,:);
    end
    
    % Find how many unique timepoints/durations there are, within a given
    % tolerance. Then 
    stimNums = uniquetol(allSineSum(:,paramCol),tol,'DataScale',1);

    stimByNum = cell(length(stimNums),1);
    for iStim = 1:length(stimNums)
        stimByNum{iStim} = ... % find param matches within tolerance and assign into stimByNum
            allSineSum(ismembertol(allSineSum(:,paramCol),stimNums(iStim),tol,'DataScale',1),:);
        try allSine{iStim,1}; %if allSine doesn't exist, initialize it
        catch
            allSine{iStim,1} = [];
        end
        
        if ~length(allSine{iStim})==0; %if values exist for that stimulus in allSine,
            %match param w/ tolerance to the matching param
            %location in allSine
            [~,whichStim]=ismembertol(stimNums,cellfun(@(x) x(1,paramCol), allSine),1,'DataScale',tol);
            allSine{whichStim(iStim),1} = [allSine{whichStim(iStim),1};stimByNum{iStim}];
            
        else %if allSine exists but that stimulus is empty, start it up
            allSine{iStim,1} = [allSine{iStim,1};stimByNum{iStim}];
        end
        
    end
    
    
    nStim = length(allSine);
    sweepsByParams = NaN(size(allSine{iStim},1),nStim);
    sortedStim = cell(1,nStim);

    
    % Sort stim profiles based on freq/ampl/dur, and take averages per
    % sine frequency.
    
    for iStim = 1:nStim
        
        nTraces = size(allSine{iStim},1);
        %LATER: don't need whichInt, just use interval previous to stim for iStim
        
        % For each parameter, go through all stimuli and round the relevant
        % parameter to the nearest X so that sweeps can be grouped as being
        % from the same stimulus profile.
        switch sortSweepsBy{iStim}
            case 'frequency'
                sweepsByParams(1:nTraces,iStim) = allSine{iStim}(:,4);
            case 'amplitude'
                sweepsByParams(1:nTraces,iStim) = allSine{iStim}(:,5);
        end
    end
    
    % Sort rows by successively less variable parameters based on
    % stimSortOrder. Use unique to find unique sets of rows/stimulus
    % profiles and separate them into groups.
    [sweepsByParams, sortIdx] = sortrows(sweepsByParams, stimSortOrder(1:size(sweepsByParams,2)));
    sortedLeakSub = allLeakSub(sortIdx,:);
    sortedSquares = stimSquareSum(sortIdx,:,:);
    
    for iStim = 1:nStim
        sortedStim{iStim} = allSine{iStim}(sortIdx,:);
    end
    
    % DECIDE: if unique(sweepsByParam) should only include columns from
    % stimSortOrder or all stim columns.
    
    [eachStimProfile, profileStartIdx, ~] = unique(sweepsByParams,'rows','first');
    [~, profileEndIdx, ~] = unique(sweepsByParams,'rows','last');
    nStimProfiles = sum(sum(~isnan(eachStimProfile),2)>=nStim);
    profileStartIdx = profileStartIdx(sum(~isnan(eachStimProfile),2)>=nStim);
    profileEndIdx = profileEndIdx(sum(~isnan(eachStimProfile),2)>=nStim);
    eachStimProfile = eachStimProfile(sum(~isnan(eachStimProfile),2)>=nStim,:);
    nReps = profileEndIdx-profileStartIdx+1;
    
    % Take mean trace from all reps across all series for each stimulus
    % profile. Find peaks in that mean trace.
    
    %     meansByStimProfile = NaN(nStimProfiles, length(sortedLeakSub));
    steadyMeans = [];
    squareMeans = [];
    groupIdx = cell(0);
    theseSweeps = cell(0);
    
    for iProfile = 1:nStimProfiles
        groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
        
        theseSweeps{iProfile} = sortedLeakSub(groupIdx{iProfile},:);
        
        %             if length(groupIdx{iProfile})>1
        %                 meansByStimProfile(iProfile,:) = nanmean(theseSweeps,1);
        %             else
        %                 meansByStimProfile(iProfile,:) = theseSweeps;
        %             end
        
        % Use the first sweep for a given size to pick the stimulus window
        % Save into parameters to pass to findMRCs.
        for iStim = 1:nStim
            
            steadyMeans(iProfile,iStim) = mean(sortedStim{iStim}(groupIdx{iProfile},6));
            stimMetaData = NaN(nStimProfiles,4);
            
            stimMetaData(:,1:2) = sortedStim{1,iStim}(profileStartIdx,1:2);
            stimMetaData(:,3) = eachStimProfile(:,iStim);
            stimMetaData(:,4) = nReps;
            
            % Find mechanoreceptor current peaks and append to sinePeaks for
            % that stimulus number.
            
        end
        
        nSqStim = size(sortedSquares,3);
        
        for iSqStim = 1:nSqStim
            squareMeans(iProfile,iSqStim) = mean(sortedSquares(groupIdx{iProfile},6,iSqStim));
        end             
        
        onRatio = sortedSquares(:,:,3)./sortedSquares(:,:,1);
        offRatio = sortedSquares(:,:,4)./sortedSquares(:,:,1);
        squareMeans(iProfile, nSqStim+1) = mean(onRatio(groupIdx{iProfile},3));        
        squareMeans(iProfile, nSqStim+2) = mean(offRatio(groupIdx{iProfile},3));
        
    end
    
    if normalizeFlag
        steadyMeans = steadyMeans./squareMeans(:,1);
    end

    sinePeaks{iCell,1} = cellName;
    sinePeaks{iCell,2} = theseSweeps';
    sinePeaks{iCell,3} = [stimMetaData steadyMeans squareMeans];        
    
end


