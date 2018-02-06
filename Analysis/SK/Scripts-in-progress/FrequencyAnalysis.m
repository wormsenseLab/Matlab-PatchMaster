% FrequencyAnalysis.m

function sinePeaks = FrequencyAnalysis(ephysData, ephysMetaData, protList, varargin)

p = inputParser;
p.addRequired('ephysData', @(x) isstruct(x));
p.addRequired('ephysMetaData', @(x) iscell(x));
p.addRequired('protList', @(x) iscell(x));

p.addOptional('allCells', cell(0));

p.addParameter('matchType', 'full', @(x) sum(strcmp(x,{'first','last','full'})));
p.addParameter('tauType','thalfmax', @(x) ischar(x) && ismember(x,{'fit' 'thalfmax'}));
p.addParameter('integrateCurrent',1);

p.parse(ephysData, ephysMetaData, protList, varargin{:});

allCells = p.Results.allCells;
matchType = p.Results.matchType;
tauType = p.Results.tauType;
integrateFlag = p.Results.integrateCurrent;

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

for iCell = 1:length(allCells)
    
    allSquareSum = [];
    allSineSum = [];
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
        steadyStateI = mean(sineTrace(end-steadyTime*sf:end));
        
        
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
        theseSquares(:,9) = repmat(sineParams(4),iStim,1);

        %TODO: add to Sine output: peak at start, beginning average, tau of 
        % fit of smoothed trace
        %TODO: add to SquareSum output: interval between squares, and
        % between sweeps
        
        % Concatenate data for all series in this recording
        allLeakSub = [allLeakSub; leakSubtractCell];
        allSquareSum = [allSquareSum; theseSquares repmat(thisSeries,size(theseSquares,1),1)];
        allSineSum = [allSineSum; sineParams steadyStateI];
    end
    
    switch sortStimBy
        case 'time'
            paramCol = 1; %set which column of allSineSum to look in
            tol = 0; %set tolerance for separating by parameter
        case 'dur'
            paramCol = 3;
            tol = 0;
    end
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
    % nStim=3;
    % allSine = allSine(1:nStim);
    sweepsByParams = NaN(size(allLeakSub,1),nStim);
    sortedStim = cell(1,nStim);
    
    sinePeaks = stimSort();
        
    %NEXT: sort stim profiles based on freq/ampl/dur, and take averages per
    %frequency.
    %DECIDE: best way to combine square data so it aligns with sine data,
    %and you only have to do one sort instead of using the nested fxn.
    
    %THEN: output freq dependence and adaptation like sinePeaks, with cell
    %name, trace, cell array containing each recording's average output
    %TODO: Call nested function here, then repeat for square analysis
    
%     squarePeaks{iCell,1} = cellName;
%     squarePeaks{iCell,2} = allLeakSub;
%     squarePeaks{iCell,3} = allSineSum;
%     squarePeaks{iCell,4} = allSquareSum;

end

keyboard;

% Nested function for sorting stim
    function sinePeaks = stimSort()
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
        
        meansByStimProfile = NaN(nStimProfiles, length(sortedLeakSub));
        
        for iProfile = 1:nStimProfiles
            groupIdx{iProfile} = profileStartIdx(iProfile):profileEndIdx(iProfile);
            
%             theseSweeps = sortedLeakSub(groupIdx{iProfile},:);
            
%             if length(groupIdx{iProfile})>1
%                 meansByStimProfile(iProfile,:) = nanmean(theseSweeps,1);
%             else
%                 meansByStimProfile(iProfile,:) = theseSweeps;
%             end
            
            % Use the first sweep for a given size to pick the stimulus window
            % Save into parameters to pass to findMRCs.
            for iStim = 1:nStim
                stimMetaData = NaN(nStimProfiles,4);
                
                stimMetaData(:,1:2) = sortedStim{1,iStim}(profileStartIdx,1:2);
                stimMetaData(:,3) = eachStimProfile(:,iStim);
                stimMetaData(:,4) = nReps;
                
                
                % Find mechanoreceptor current peaks and append to seriesPeaks for
                % that stimulus number.
                sinePeaks{iCell,1} = cellName;
                sinePeaks{iCell,2} = 1;
                
            end
            
        end
    end

end


