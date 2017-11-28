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
%   mechPeaks = IdAnalysis(ephysData, allCells, calibFlag)
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
%                                   [size/velocity (um or um/s); location 
%                                   (datapoint); peak current (pA); 
%                                   direction (+/-); threshold used (pA);
%                                   decay tau from exp1 fit (ms);
%                                   number of sweeps averaged].
%
%
% Created by Sammy Katta on 20-May-2015.
% Updated most recently by Sammy Katta on 7-Aug-2017.
% Documentation not updated yet, but fxn now allows output to be used for
% many more types of analyses. Will soon rename.

% TODO: Add flag for using stim com signal vs. PD signal (chan 2 vs chan 3)
%   Done, but still uses stim com to group step sizes, and includes PD
%   trace and calculated sizes (based on calibration if available) in
%   output
% TODO: allCells usage is old? Read from unique(mechTracePicks(:,1))
% NEXT: test with stimInterval, break down into intermediate fxns if still
% necessary.
% LATER: use inputParser to make more useful for other kinds of analysis.
% LATER: add in PD step analysis (will need findMRCs modified), output
% sorted stim command and PD traces with sortedLeakSub
% LATER: add GUI for selecting sortSweepsBy for each stim segment

function [mechPeaks, mechStim, mechTraces] = IdAnalysis(ephysData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('matchType', 'full', @(x) ischar(x));
p.addOptional('sortStimBy', 'num', @(x) sum(strcmp(x,{'num','time'})));
p.addOptional('calibFlag', 0, @(x) islogical(x));

p.addParameter('tauType','fit', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));

p.parse(ephysData, protList, varargin{:});

matchType = p.Results.matchType;
sortStimBy = p.Results.sortStimBy;
calibFlag = p.Results.calibFlag;
tauType = p.Results.tauType;

stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
% sortStimBy = 'num';
% sortSweepsBy = {'magnitude', 'magnitude','magnitude', 'magnitude', 'magnitude'};
sortSweepsBy = {'velocity', 'magnitude','magnitude', 'magnitude'};
roundIntTo = 2;
whichInt = 1;
fillZeroSteps = 1;
sortByStimNum = 1; %sort by which stim (here, first stim, which is on step)
stimSortOrder = [1 2];

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

allCells = unique(mechTracePicks(:,1));

mechPeaks = cell(length(allCells),1);
% protList = {'WC_Probe','NoPrePulse','DispRate'};
% protList = {'PrePulse'};

% Find applicable series and check against list of included series/traces
% (this allows a cross-check on the protocol name) before analyzing
% Values for traces not on the list will be stored as NaN.
for iCell = 1:length(allCells)
    
    allStim = cell(0);
    allLeakSub = cell(0);
    
    % Double check that the list of series given matches the indices of
    % the WC_Probe Id-curve protocols.
    
%SPLIT into function findStimuli here, through end of for series loop
%Consider if you want to output stim and PD traces?
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,protList,'MatchType',matchType);
    
    nSeries = length(allSeries);
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
    if nSeries == 0
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
        % sampling frequency in kHz
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        dataType = ephysData.(cellName).dataunit{1,thisSeries};
        nSweeps = size(stimComI,2);
        
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
        
        
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtractCell = num2cell(leakSubtract',2);
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);
        
        % add series number to seriesStimuli for referencing back from
        % analyzed data to a particular trace
        seriesStimuli (:,8) = repmat(thisSeries,size(seriesStimuli,1),1);
        
        
        
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
        allLeakSub=[allLeakSub; leakSubtractCell];
        
    end
    
    
    % If using the time sort, you likely want to compare a given
    % timepoint across sweeps/series, even if the magnitude of the
    % step there is zero and it isn't detected. This section fills
    % in those timepoints with magnitude zero steps, if the fillZeroSteps
    % flag is set to 1.

    if ~isempty(strfind(lower(sortStimBy),'time')) && fillZeroSteps
        seriesList = vertcat(allStim{:});
        seriesList = seriesList(:,[6 8]);
        seriesList = sortrows(unique(seriesList,'rows','stable'),[2 1]); % find all unique series/sweep combos
        blankStim = NaN(length(seriesList),8);
        blankStim(:,[6,8]) = seriesList;
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
        seriesList = seriesList(:,[6 8]);
        seriesList = sortrows(unique(seriesList,'rows','stable'),[2 1]); % find all unique series/sweep combos
        blankStim = NaN(length(seriesList),8);
        blankStim(:,[6,8]) = seriesList;
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
    allLeakSub=cellfun(@(x)cat(2,x,NaN(1,maxLength-length(x))),allLeakSub,'UniformOutput',false);
    allLeakSub = cell2mat(allLeakSub);
    
%FIX: nStim=4 for OnDt because of Patchmaster's 0 padding. Find a way to
%exclude the last "step".
    nStim = length(allStim);
% nStim=3;
% allStim = allStim(1:nStim);
    sweepsByParams = NaN(size(allLeakSub,1),nStim);
    sortedStim = cell(1,nStim);
    
% KEEP this section separate for each Analysis fxn, because what you sort
% by will change. (e.g., not by size for OnDt analysis).
    
    for iStim = 1:nStim
    
    nTraces = size(allStim{iStim},1);
    %LATER: don't need whichInt, just use interval previous to stim for iStim

        % For each parameter, go through all stimuli and round the relevant
        % parameter to the nearest X so that sweeps can be grouped as being
        % from the same stimulus profile.
        switch sortSweepsBy{iStim}
            case 'magnitude'
                sweepsByParams(1:nTraces,iStim) = round(allStim{iStim}(:,3),1);
            case 'position'
                sweepsByParams(1:nTraces,iStim) = round(allStim{iStim}(:,4),1);
            case 'velocity'
                trioRound = round(allStim{iStim}(:,5),0); %round small velocities to nearest 1
                trioRound(abs(trioRound)<500) = round(trioRound(abs(trioRound)<500)/2)*2; %round large to nearest 10
                trioRound(abs(trioRound)>500 & abs(trioRound)<1500) = ...
                    round(trioRound(abs(trioRound)>500 & abs(trioRound)<1500)/5)*5; %round medium to nearest 5
                trioRound(abs(trioRound)>1500) = round(trioRound(abs(trioRound)>1500),-1); %round large to nearest 10
                sweepsByParams(1:nTraces,iStim) = trioRound;
            case 'interval' %time interval between previous stim and current stim
                a = cellfun(@(x) x(:,1:2),allStim,'un',0);
                stimInt = diff([a{:}],1,2)/sf; %calculate interval between stimuli in ms
                stimInt = round(stimInt(:,2:2:size(stimInt,2))/roundIntTo)*roundIntTo;
                sweepsByParams(:,iStim) = stimInt(:,whichInt);
                %             case 'none'
                %             otherwise
                %                 [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,3),1));
        end
    end
    
    % Sort rows by successively less variable parameters based on
    % stimSortOrder. Use unique to find unique sets of rows/stimulus
    % profiles and separate them into groups.
    [sweepsByParams, sortIdx] = sortrows(sweepsByParams, stimSortOrder);
    sortedLeakSub = allLeakSub(sortIdx,:);
    
    for iStim = 1:nStim
        sortedStim{iStim} = allStim{iStim}(sortIdx,:);
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
    
    meansByStimProfile = NaN(nStimProfiles, length(sortedLeakSub));
    
    for iProfile = 1:nStimProfiles
        groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
        
        theseSweeps = sortedLeakSub(groupIdx{iProfile},:);
        
        if length(groupIdx{iProfile})>1
            meansByStimProfile(iProfile,:) = nanmean(theseSweeps,1);
        else
            meansByStimProfile(iProfile,:) = theseSweeps;
        end
        
        % Use the first sweep for a given size to pick the stimulus window
        % Save into parameters to pass to findMRCs.
        for iStim = 1:nStim
            stimMetaData = NaN(nStimProfiles,4);
            
            stimMetaData(:,1:2) = sortedStim{1,iStim}(profileStartIdx,1:2);
            stimMetaData(:,3) = eachStimProfile(:,iStim);
            stimMetaData(:,4) = nReps;
            
            
            % Find mechanoreceptor current peaks and append to seriesPeaks for
            % that stimulus number.
            mechPeaks{iCell,1} = cellName;
            mechPeaks{iCell,2} = meansByStimProfile;
            mechPeaks{iCell,2+iStim} = findMRCs(stimMetaData, meansByStimProfile, sf, dataType, ...
                'tauType', tauType);
            
        end
        
    end
    
    
    % TODO: have column 2  of sortedLeakSub be the stim trace used to find stim, for easy
    % plotting access, and col 3 = PD trace, once you have that set up.
    
    % keyboard;
    
    %
    %
    %     if calibFlag==1
    %         mechPeaks{iCell,1} = [eachSize(~isnan(eachSize)) meanPDSize(~isnan(eachSize)) ...
    %             pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc nReps];
    %         mechPeaks{iCell,2} = meansByParam;
    %         mechPeaks{iCell,3} = meanPDTrace;
    %         mechPeaks{iCell,4} = repmat(cellName,[size(pkOn),1]);
    %         mechPeaks{iCell,5} = theseIDs;
    %     else
    %         mechPeaks{iCell,1} = ...
    %             [eachSize(~isnan(eachSize)) nan(size(eachSize(~isnan(eachSize))))...
    %             pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc nReps];
    %         mechPeaks{iCell,2} = meansByParam;
    %         mechPeaks{iCell,4} = repmat(cellName,[size(pkOn),1]);
    %         mechPeaks{iCell,5} = theseIDs;
    %     end
    %     % TODO: Figure out how to fit this to the four-parameter sigmoidal
    %     % function used in O'Hagan: @(X,a,b,c,d) ((a-d)/(1+((X/c)^b)))+d
    %     % Using optimtool? fmincon? nlinfit if you add the statistics toolbox.
    
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

    mechStim (iCell,2:length(sortedStim)+1) = sortedStim;

end

end