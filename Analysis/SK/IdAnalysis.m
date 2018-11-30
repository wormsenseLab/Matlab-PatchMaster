% IdAnalysis.m
%
% This function calculates the mean peak current for each step size across
% a given recording, for making an I-d or I-x curve. Step size is
% calculated using findSteps with the stimulus command signal.
%
% Stimulus command voltage to step size conversion is hardcoded for the
% current setup.
%
% USAGE:
%   mechPeaks = IdAnalysis(ephysData, protList)
%   mechPeaks = IdAnalysis(ephysData, protList, 'sortStimBy', 'num', calibFlag, 0)
%
% INPUTS:
%   ephysData       struct          Imported data from ImportPatchData.
%
%   protList        cell array      Names of protocols to match (this is a
%                                   redundancy measure to make sure the
%                                   settings are correct for the particular
%                                   xls file you are inputting. Default is
%                                   to match the whole name.
% 
% OPTIONAL INPUTS:
%   allCells        cell array      Names of recordings to look at (this
%                                   can be used to narrow down analysis to
%                                   a subset of recordings from the
%                                   ExcludeSweeps metadata sheet).
%   
% 
% OPTIONAL PARAMETER-VALUE INPUTS:
% 
%   matchType       char            Specifies whether to match names in
%                   {'full'         protList to beginning, end, or full
%                    'first'        protocol names in ephysData. Default is
%                    'last'}        'full'.
%                                   
%   sortStimBy      char            Specifies whether to separate stimuli
%                   {'num'          by number (useful for protocols with
%                    'time'}        same # of steps at different times) or
%                                   by time (for protocols that may have
%                                   zero-magnitude steps at a given time).
%                                   Default is 'num'.
% 
%   calibFlag       logical         Specifies whether to use photodiode
%                                   signal data instead of stimulus command
%                                   for step sizes (not currently working).
%                                   Default is 0 (just use command signal).
%
% PROMPTED INPUTS:
%   ImportMetaData asks for a metadata file in .xls format containing the
%   list of traces to analyze in the same format as files output by
%   ExcludeSweeps(). This will get double-checked against recordings that 
%   match the protocol list.
%
% OUTPUTS:
%   mechPeaks       cell array      Nested cell array with a row for each
%                                   recording. CellName, Average Trace, 
%                                   then for each stimulus, MRC stats on
%                                   [1-sortingParam, 2-size (um), 3-position (um), 
%                                   4-velocity (um/s), 5-pkLoc (datapoint), 
%                                   6-peak (pA), 7-threshold used (pA), 
%                                   8-activation time (ms), 9-decay time (ms) *either t1/2 max or tau
%                                   10-time to peak (ms), 11-charge (C), 
%                                   12-stim distance (um), 13-number of sweeps averaged.
% 
%   mechStim        cell array      Nested cell array with a row for each
%                                   recording. CellName, Average stimulus
%                                   commands, (Average PD response if
%                                   pdCompare is set to true), (Sweeps of
%                                   stimulus commands and pd response used
%                                   for means if saveSweeps is true),
%                                   sortedStim metadata (not separated).
% 
%
%
% Created by Sammy Katta on 20-May-2015.
% Updated most recently by Sammy Katta on 13-Nov-2018.
% Documentation not updated yet, but fxn now allows output to be used for
% many more types of analyses. Will soon rename.

% NEXT: test with stimInterval, break down into intermediate fxns if still
% necessary.
% LATER: add GUI for selecting sortSweepsBy for each stim segment

function [mechPeaks, mechStim] = IdAnalysis(ephysData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));
p.addOptional('sortStimBy', 'num', @(x) sum(strcmp(x,{'num','time','endTime'})));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));
p.addParameter('calibFlag', 0); %0 for stim command, 1 for PD signal
p.addParameter('tauType','fit', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('sortSweepsBy',{'magnitude','magnitude','magnitude','magnitude'}, @(x) iscell(x));
p.addParameter('integrateCurrent',0);
p.addParameter('fillZero',1);
p.addParameter('sepByStimDistance',0);
p.addParameter('saveSweeps',0);
p.addParameter('recParameters', cell(0), @(x) iscell(x)); % ephysMetaData cell array
p.addParameter('limitStim',0, @(x) x>=0);
p.addParameter('pdCompare',0);
p.addParameter('saveZero',0);
p.addParameter('dropNeg',0); %

p.parse(ephysData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;
sortStimBy = p.Results.sortStimBy;
calibFlag = p.Results.calibFlag;
tauType = p.Results.tauType;
sortSweepsBy = p.Results.sortSweepsBy;
integrateFlag = p.Results.integrateCurrent;
fillZeroSteps = p.Results.fillZero;
distFlag = p.Results.sepByStimDistance;
sweepFlag = p.Results.saveSweeps; % if 1, save individual sweeps in addition to means
ephysMetaData = p.Results.recParameters;
limitStim = p.Results.limitStim;
pdCompare = p.Results.pdCompare;
saveLeak = p.Results.saveZero;
dropNeg = p.Results.dropNeg;

baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
% sortStimBy = 'num';
% sortSweepsBy = {'velocity', 'magnitude','magnitude', 'magnitude'};
roundIntTo = 2;
whichInt = 1;
sortByStimNum = 1; %sort by which stim (here, first stim, which is on step)

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

% Allow for pre-filtered subsets of "allCells" (from FilterRecordings) to 
% be used - otherwise take names from imported metadata sheet.
if isempty(allCells)
    allCells = unique(mechTracePicks(:,1));
end

mechPeaks = cell(length(allCells),1);
mechStim = cell(length(allCells),1);

% protList = {'WC_Probe','NoPrePulse','DispRate'};
% protList = {'PrePulse'};

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


% Find applicable series and check against list of included series/traces
% (this allows a cross-check on the protocol name) before analyzing. 
% Values for traces not on the list will be stored as NaN.
for iCell = 1:length(allCells)
    
    allStim = cell(0);
    allLeakSub = cell(0);
    allStimCom = cell(0);
    allPD = cell(0);
    allLeak = [];
    
    % Double check that the list of series given matches the indices of the
    % protocol names specified in the input.
    
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,protList,'MatchType',matchType);
    
    nSeries = length(allSeries);
    
    %NEXT TODO: Read in ext stim filter freq from ephysmetadata
    %sheet after matching cellName.
    
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
        
        if pdCompare
            pdV = ephysData.(cellName).data{3,thisSeries}(:,pickedTraces);
        end
        % sampling frequency in kHz
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        dataType = ephysData.(cellName).dataunit{1,thisSeries};
        nSweeps = size(stimComI,2);
        try thisDist = round(cellDist{iCell},1);
        catch
        end
        
        %DECIDE: where use of photodiode signal should come in. Use the stimulus
        %command channel (or stim parameters from stimTree once that's working) to
        %find the putative step locations/stimulus profile. Then use the photodiode
        %trace to find the actual size/speed of the stimuli at that location. Find
        %PD-measured step sizes/speeds here and give optional output allPDStim.
        
        
        %         % if calibration should be used to calculate step sizes, get photodiode data
        %         % and use calib curve to transform photodiodeV trace into
        %         % measured displacement trace (then use same findSteps threshold)
        %         if calibFlag == 1
        %             try pdCalib = ephysData.(cellName).calibration;
        %             catch
        %                 fprintf('No calibration found for %s\n', cellName);
        %                 calibFlag = 2;
        %             end
        %
        %             if calibFlag ==1
        %
        %                 photodiodeV = ephysData.(cellName).data{3,thisSeries}(:,pickedTraces);
        %
        %                 if isempty(photodiodeV)
        %                     photodiodeV = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces);
        %                 end
        %
        %
        %                 % Interpolate photodiode voltage to calculate measured disp
        %                 try measuredDisp = interp1(-pdCalib(2,:), pdCalib(1,:), -photodiodeV, 'linear','extrap');
        %                 catch
        %                     fprintf('Interpolation failed for %s, series %d\n',cellName,thisSeries);
        %                     calibFlag = 2;
        %                 end
        %                 %TODO: swap this for newStepFind as well
        %                 try [pdStepSize, pdStepStarts, pdStepEnds] = ...
        %                         findSteps(nSweeps, measuredDisp, sf, stepThresh, 'roundedTo', 0.05);
        %                 catch
        %                 end
        %
        %             end
        %             %TODO: Change roundedTo parameter for this use
        %             %TODO: Check that stepThresh is applicable for measuredDisp as well
        %         end
        
        
        [leakSubtract, leak] = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtract = num2cell(leakSubtract',2);
        if pdCompare
            zeroStimCom = SubtractLeak(stimComI, sf, 'BaseLength', baseTime);
            zeroPD = SubtractLeak(pdV, sf, 'BaseLength', baseTime);
            zeroStimCom = num2cell(zeroStimCom',2);
            zeroPD = num2cell(zeroPD',2);
        else
            zeroStimCom = num2cell(stimComI',2);
        end
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);
        
        % add series number to seriesStimuli for referencing back from
        % analyzed data to a particular trace
        seriesStimuli (:,8) = repmat(thisSeries,size(seriesStimuli,1),1);
        
        if exist('thisDist','var') && ~isempty(thisDist) && ~isnan(thisDist)
            seriesStimuli(:,9) = repmat(thisDist, size(seriesStimuli,1), 1);;
        else
            % if no stim distance is included, must be zeros instead of
            % NaNs because unique will consider every NaN unique when
            % sorting by stim for zeroFill later.
            seriesStimuli(:,9) = zeros(size(seriesStimuli,1),1);
        end
        %NEXT: where do I have to add in this parameter so I can use it when
        %sorting for unique stimProfiles? It has to get into sweepsByParams.
        
        
        % add sampling freq because that might vary among series that might
        % otherwise be grouped as well, and you should use it to check that
        % they're actually all the same.
        seriesStimuli(:,10) = repmat(sf, size(seriesStimuli,1),1);
        
        
        % all on/off stimuli are still mixed together.
        % depending on input param, use either stimulus number or stimulus
        % timepoint to group stimuli across sweeps (e.g., on and off)
        switch sortStimBy
            case 'num'
                paramCol = 7; %set which column of seriesStimuli to look in
                tol = 0; %set tolerance for separating by parameter
            case 'time'
                paramCol = 1;
                tol = 1;
            case 'endTime'
                paramCol = 2;
                tol = 0;
        end
        stimNums = uniquetol(seriesStimuli(:,paramCol),tol,'DataScale',1);
        
        
        % separate found stim parameters into groups and concatenate to
        % previous stimuli from the same series
        
        stimByNum = cell(length(stimNums),1);
        for iStim = 1:length(stimNums)
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
        
        % Concatenate to the complete list of step sizes and
        % leak-subtracted traces across series for this recording
        allLeakSub = [allLeakSub; leakSubtract];
        allStimCom = [allStimCom; zeroStimCom];
        if pdCompare
            allPD = [allPD; zeroPD];
        end
        if saveLeak
            allLeak = [allLeak; leak];
        end
        clear thisDist;
    end
    
    
    % If using the time sort, you likely want to compare a given
    % timepoint across sweeps/series, even if the magnitude of the
    % step there is zero and it isn't detected. This section fills
    % in those timepoints with magnitude zero steps, if the fillZeroSteps
    % flag is set to 1.
    
    if ~isempty(strfind(lower(sortStimBy),'time')) && fillZeroSteps
        seriesList = vertcat(allStim{:});
        seriesList = seriesList(:,[6 8 9 10]);
        seriesList = sortrows(unique(seriesList,'rows','stable'),[2 1]); % find all unique series/sweep combos
        blankStim = NaN(length(seriesList),size(allStim{1},2));
        blankStim(:,[6 8 9 10]) = seriesList;
        filledStim = cell(length(allStim),1);
        [filledStim{1:length(allStim),1}] = deal(blankStim); %make a cell with nStim with the series/sweep list to be filled for each stim
        
        [~,ia,ib] = cellfun(@(a,b) intersect(a(:,[6 8]),b(:,[6 8]),'rows','stable'), ...
            filledStim, allStim, 'un', 0); % for each stim, find the rows for which entries already exist
        for iStim = 1:length(allStim)
            filledStim{iStim}(ia{iStim},[1:5, 7]) = allStim{iStim}(ib{iStim},[1:5, 7]);
            
            blankIdx = find(isnan(filledStim{iStim}(:,1))); % for each stim, find sweeps where entries don't exist for that stim
            insertTimes = repmat(allStim{iStim}(1,1:2),length(blankIdx),1); %copy timepoint for given stim
            try insertTotalDisp = filledStim{iStim-1}(blankIdx,4); %take cumulative/total displacement from the stim before that
            catch
                insertTotalDisp = zeros(length(blankIdx),1); %if the previous stim was the first one, total displacement is zero
            end
            
            filledStim{iStim}(blankIdx,1:5) = ...
                [insertTimes, zeros(length(blankIdx),1), insertTotalDisp, zeros(length(blankIdx),1)];
        end
        
        allStim = filledStim;
        
    elseif ~isempty(strfind(lower(sortStimBy),'num')) && fillZeroSteps
        seriesList = vertcat(allStim{:});
        seriesList = seriesList(:,[6 8 9 10]);
        seriesList = sortrows(unique(seriesList,'rows','stable'),[2 1]); % find all unique series/sweep combos
        blankStim = NaN(size(seriesList,1),size(allStim{1},2));
        blankStim(:,[6 8 9 10]) = seriesList;
        filledStim = cell(length(allStim),1);
        [filledStim{1:length(allStim),1}] = deal(blankStim); %make a cell with nStim with the series/sweep list to be filled for each stim
        
        [~,ia,ib] = cellfun(@(a,b) intersect(a(:,[6 8]),b(:,[6 8]),'rows','stable'), ...
            filledStim, allStim, 'un', 0); % for each stim, find the rows for which entries already exist
        for iStim = 1:length(allStim)
            filledStim{iStim}(ia{iStim},[1:5, 7]) = allStim{iStim}(ib{iStim},[1:5, 7]);
            
        end
        
        allStim = filledStim;
        
    end
    
    
    
    
    % Pad all traces with NaNs at the end so they're the same length, for
    % ease of manipulation as an array instead of a cell.
    sweepLengths = cellfun('length',allLeakSub);
    maxLength = max(sweepLengths);
    allLeakSub = cellfun(@(x)cat(2,x,NaN(1,maxLength-length(x))),allLeakSub,'UniformOutput',false);
    allLeakSub = cell2mat(allLeakSub);    

    allStimCom = cellfun(@(x)cat(2,x,NaN(1,maxLength-length(x))),allStimCom,'UniformOutput',false);
    allStimCom = cell2mat(allStimCom);
    
    allPD = cellfun(@(x)cat(2,x,NaN(1,maxLength-length(x))),allPD,'UniformOutput',false);
    allPD = cell2mat(allPD);

    %FIX: nStim=4 for OnDt because of Patchmaster's 0 padding. Find a way to
    %exclude the last "step".
    if limitStim == 0
        nStim = length(allStim);
    elseif limitStim > 0
        nStim = min([length(allStim) limitStim]);
    end
    
    stimSortOrder = 1:nStim;
    stimTolVal = [];
    
    sweepsByParams = NaN(size(allLeakSub,1),nStim);
    sortedStim = cell(1,nStim);
        
    % KEEP this section separate for each Analysis fxn, because what you sort
    % by will change. (e.g., not by size for OnDt analysis).
    
    for iStim = 1:nStim
        
        nTraces = size(allStim{iStim},1);
        %TODO: don't need whichInt, just use interval previous to stim for iStim
        
        % For each parameter, go through all stimuli and round the relevant
        % parameter to the nearest X so that sweeps can be grouped as being
        % from the same stimulus profile.
        switch sortSweepsBy{iStim}
            case 'magnitude'
                sweepsByParams(1:nTraces,iStim) = round(allStim{iStim}(:,3),1);
                stimTolVal(iStim) = 0.1;
            case 'position'
                sweepsByParams(1:nTraces,iStim) = round(allStim{iStim}(:,4),1);
                stimTolVal(iStim) = 0.1;
            case 'velocity'
                trioRound = roundVel(allStim{iStim}(:,5)); %fxn found at end of file
                sweepsByParams(1:nTraces,iStim) = trioRound;
                stimTolVal(iStim) = 11;
            case 'interval' %time interval between previous stim and current stim
                a = cellfun(@(x) x(:,1:2),allStim,'un',0);
                stimInt = diff([a{:}],1,2)/sf; %calculate interval between stimuli in ms
                stimInt = round(stimInt(:,2:2:size(stimInt,2))/roundIntTo)*roundIntTo;
                sweepsByParams(:,iStim) = stimInt(:,whichInt);
                %             case 'none'
                %             otherwise
                %                 [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,3),1));
                stimTolVal(iStim) = 1/sf;
        end
    end
    
    dists = allStim{1}(:,9);
    if distFlag
        sweepsByParams = [sweepsByParams dists];
        stimSortOrder = 1:nStim+1;
        stimTolVal(nStim+1) = 0;
    end
    
    if dropNeg
        dropRows = sum(~isnan(sweepsByParams(:,1:nStim)),2)<2;
        sweepsByParams = sweepsByParams(~dropRows,:);
        allLeakSub = allLeakSub(~dropRows,:);
        allStimCom = allStimCom(~dropRows,:);
        if pdCompare
            allPD = allPD(~dropRows,:);
        end
        if saveLeak
            allLeak = allLeak(~dropRows);
        end
    end
    
    % Sort rows by successively less variable parameters based on
    % stimSortOrder. Use unique to find unique sets of rows/stimulus
    % profiles and separate them into groups.
    
    [~, sortIdx, eachStimProfile, profileStartIdx, profileEndIdx] = ...
        sortRowsTol(sweepsByParams, stimTolVal, stimSortOrder);
    sortedLeakSub = allLeakSub(sortIdx,:);
    sortedStimCom = allStimCom(sortIdx,:);
    if pdCompare
        sortedPD = allPD(sortIdx,:);
    end
    if saveLeak
        sortedLeak = allLeak(sortIdx,:);
    end
    
    for iStim = 1:nStim
        sortedStim{iStim} = allStim{iStim}(sortIdx,:);
    end
      
    
    %DECIDE: does sortRowsTol need to include a flag for excluding rows
    %containing nans as stim parameters? Can't remember why that was
    %necessary.
    
%     [eachStimProfile, profileStartIdx, ~] = unique(sweepsByParams,'rows','first');
%     [~, profileEndIdx, ~] = unique(sweepsByParams,'rows','last');
    nStimProfiles = sum(sum(~isnan(eachStimProfile),2)>=nStim);
%     profileStartIdx = profileStartIdx(sum(~isnan(eachStimProfile),2)>=nStim);
%     profileEndIdx = profileEndIdx(sum(~isnan(eachStimProfile),2)>=nStim);
%     eachStimProfile = eachStimProfile(sum(~isnan(eachStimProfile),2)>=nStim,:);
    nReps = profileEndIdx-profileStartIdx+1;
    
    % Take mean trace from all reps across all series for each stimulus
    % profile. Find peaks in that mean trace.
    
    meansByStimProfile = NaN(nStimProfiles, length(sortedLeakSub));
    stimComByStimProfile = NaN(nStimProfiles, length(sortedStimCom));
    if pdCompare
        pdByStimProfile = NaN(nStimProfiles, length(sortedPD));
    end
    
    sweepsByStimProfile = cell(nStimProfiles,1);
    stimByStimProfile = cell(nStimProfiles,2);

    stimMetaData = NaN(nStimProfiles,4,nStim);
    
    %PROBLEM: When using different protocols that have the same value for
    %sortSweepsBy, e.g., same 20mm/s speed but one with sf 5 vs. 10 kHz,
    %the protocols are grouped (can't sortStimBy time because we want to
    %compare velocities). But then the sweeps are averaged before
    %MRC-finding, which doesn't work because the timing is different, and
    %the start time is off for at least one set. 
    %
    %Include sf in stim metadata.
    %
    %Try to separate? How? By checking sf? Not the only situation where
    %this might happen, right? (e.g., protocol where times have been
    %changed). So instead of taking the mean out here, maybe we want to pass all the
    %sweeps of interest to findMRCs along with their individual metadata
    %(because sortedStim at this point should have the correct start/end
    %timepoints for each sweep). Then check if sf and stim start location
    %are consistent within a grouping. If only loc is different, it's easy,
    %just align based on stimStart, then take the mean and work as usual.
    %
    %For saving, if the sweeps have different timing or sf, you obviously
    %can't take the mean, so drop zeroes in meansByStimProfile, and save
    %sweepsByStimProfile (may require making includeSweeps default true). 
    %
    %If sf is different (may need to nest these two checks), take each set
    %of sweeps and do the usual mean finding and peak finding, then use a
    %new flag combineSweeps to either output those separately, or take a
    %weighted mean based on nReps for each.
    
    for iProfile = 1:nStimProfiles
        groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
        
        theseSweeps = sortedLeakSub(groupIdx{iProfile},:);
        theseStimCom = sortedStimCom(groupIdx{iProfile},:);
        if pdCompare 
            thesePD = sortedPD(groupIdx{iProfile},:);
        end
        if saveLeak
            theseLeak = sortedLeak(groupIdx{iProfile},:);
        end
        
        if length(groupIdx{iProfile})>1
            meansByStimProfile(iProfile,:) = nanmean(theseSweeps,1);
            stimComByStimProfile(iProfile,:) = nanmean(theseStimCom,1);
            if pdCompare
                pdByStimProfile(iProfile,:) = nanmean(thesePD,1);
            end
            if saveLeak
                leakByStimProfile(iProfile,:) = nanmean(theseLeak);
            end
        else
            meansByStimProfile(iProfile,:) = theseSweeps;
            stimComByStimProfile(iProfile,:) = theseStimCom;
            if pdCompare 
                pdByStimProfile(iProfile,:) = thesePD;
            end
            if saveLeak
                leakByStimProfile(iProfile,:) = theseLeak;
            end

        end
        
        if sweepFlag
            sweepsByStimProfile{iProfile} = theseSweeps;
            stimByStimProfile{iProfile,1} = theseStimCom;
            if pdCompare 
                stimByStimProfile{iProfile,2} = thesePD;
            end
        end
        
        % Use the first sweep for a given size to pick the stimulus window
        % Save into parameters to pass to findMRCs.
        
        for iStim = 1:nStim
            
            stimMetaData(:,1:2,iStim) = sortedStim{1,iStim}(profileStartIdx,1:2);
            stimMetaData(:,3,iStim) = eachStimProfile(:,iStim);
            stimMetaData(:,4:6,iStim) = round(sortedStim{1,iStim}(profileStartIdx,3:5),1);
            stimMetaData(:,7,iStim) = roundVel(sortedStim{1,iStim}(profileStartIdx,5));
            stimMetaData(:,8,iStim) = nReps;
            if distFlag
                stimMetaData(:,9,iStim) = eachStimProfile(:,end);
            else
                if sum(diff(dists))==0
                    stimMetaData(:,9,iStim) = repmat(dists(1),size(eachStimProfile,1),1);
                else
                    stimMetaData(:,9,iStim) = zeros(size(eachStimProfile,1),1);
                    fprintf('Multiple stim distances found for %s.\n Set distFlag to true to separate by stim distance.',cellName);
                end
            end
            
        end
    end
    
    if sweepFlag % if saving sweeps separate from means, scoot everything over 1
        mechPeaks{iCell,1} = cellName;
        mechPeaks{iCell,2} = meansByStimProfile;
        mechPeaks{iCell,3} = sweepsByStimProfile;
        
        for iStim = 1:nStim
            % Find mechanoreceptor current peaks and append to seriesPeaks for
            % that stimulus number.
            mechPeaks{iCell,3+iStim} = findMRCs(stimMetaData(:,:,iStim), meansByStimProfile, sf, dataType, ...
                'tauType', tauType, 'integrateCurrent',integrateFlag);
        end
    else 
        mechPeaks{iCell,1} = cellName;        
        mechPeaks{iCell,2} = meansByStimProfile;
        
        for iStim = 1:nStim
            % Find mechanoreceptor current peaks and append to seriesPeaks for
            % that stimulus number.
            mechPeaks{iCell,2+iStim} = findMRCs(stimMetaData(:,:,iStim), meansByStimProfile, sf, dataType, ...
                'tauType', tauType, 'integrateCurrent',integrateFlag);
        end
    end
    
    if saveLeak
        mechPeaks{iCell,nStim+3} = leakByStimProfile;
    end
    % TODO: have column 2  of sortedLeakSub be the stim trace used to find stim, for easy
    % plotting access, and col 3 = PD trace, once you have that set up.
    
    % reset calibFlag to true if it was unset for a particular cell
    if calibFlag == 2
        calibFlag = 1;
    end
    
    %DECIDE: which info do you want about the stimuli? Currently only have
    %individual stim traces. Maybe this output can be mean or representative stim profile, and
    %output 3 can be sortedLeakSub traces + sortedStim traces/data?
    %NEXT: save stimComI and PD trace into sortedStimTraces,  take mean PD
    %param?? Figure out how to end up with mean stim profile traces in both
    %channels.
    %Difficulty: you're using seriesStimuli to sort traces, so you need to
    %find stim on individual traces first no matter what. Then you could
    %group and take means of stim traces and run series stim again on the
    %means and use those timepoints for peak finding so you'd also have a
    %less noisy stim/PD trace.
    if ~sweepFlag
        mechStim(iCell,1) = mechPeaks(iCell,1);
        mechStim{iCell,2} = stimComByStimProfile;
        if pdCompare
            mechStim{iCell,3} = pdByStimProfile;
        end
        mechStim(iCell,4:length(sortedStim)+3) = sortedStim;
    else
        mechStim(iCell,1) = mechPeaks(iCell,1);
        mechStim{iCell,2} = stimComByStimProfile;
        if pdCompare
            mechStim{iCell,3} = pdByStimProfile;
        end
        mechStim{iCell,4} = stimByStimProfile;
        mechStim(iCell,5:length(sortedStim)+4) = sortedStim;
    end
end
mechStim = mechStim(~cellfun(@isempty, mechPeaks(:,1)),:);
mechPeaks = mechPeaks(~cellfun(@isempty, mechPeaks(:,1)),:);

end