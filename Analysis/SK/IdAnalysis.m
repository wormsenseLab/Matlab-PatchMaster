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
%   allCells        cell array      List of recording names to analyze.
%
%   calibFlag       logical         Flag specifying whether to use
%                                   photodiode calibration and photodiode
%                                   signal to determine actual step size
%                                   for all given cells.
%                                   0 = don't use calib, 1 = use calib.
%                                   (2 = internally set, skip calib for cell).
%
% PROMPTED INPUTS:
%   ImportMetaData asks for a metadata file in .xls format containing the
%   list of traces to analyze in the same format as files output by
%   ExcludeSweeps(). This will get double-checked against allCells.
%
% OUTPUTS:
%   mechPeaks       cell array      Nested cell array with a cell for each
%                                   recording. Columns per recording:
%                                   [step size (um); peak current at step
%                                   onset (pA); peak current at offset;
%                                   onset tau (ms); offset tau; onset
%                                   location (sample); offset location]
%
%
% Created by Sammy Katta on 20-May-2015.

% TODO: Add flag for using stim com signal vs. PD signal (chan 2 vs chan 3)
%   Done, but still uses stim com to group step sizes, and includes PD
%   trace and calculated sizes (based on calibration if available) in
%   output
% TODO: allCells usage is old? Read from unique(mechTracePicks(:,1))
% TODO: 

function [mechPeaks, finalStim, finalLeakSub] = IdAnalysis(ephysData, calibFlag)

% keyboard;

stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
sortStimBy = 'num';
sortSweepsBy = 'step magnitude';
roundIntTo = 2;
whichInt = 1;
sortByStimNum = 1; %sort by which stim (here, first stim, which is on step)

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

allCells = unique(mechTracePicks(:,1));

mechPeaks = cell(length(allCells),1);
protList = {'WC_Probe','WC_ProbeLarge','WC_ProbeSmall'};

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
    allSeries = matchProts(ephysData,cellName,protList,'MatchType','full');

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
                stimNums = unique(seriesStimuli(:,7));
            case 'time'
                stimNums = unique(seriesStimuli(:,1));
            case 'endTime'
                stimNums = unique(seriesStimuli(:,2));
        end
        
        % separate found stim parameters into groups and concatenate to
        % previous stimuli from the same series
        stimByNum = cell(length(stimNums),1);
        for iStim = 1:length(stimNums)
            stimByNum{iStim} = seriesStimuli(seriesStimuli(:,7)==stimNums(iStim),:);
            try allStim{iStim,1};
            catch
                allStim{iStim,1} = [];
            end
            allStim{iStim,1} = [allStim{iStim,1};stimByNum{iStim}];
        end
        
        % Concatenate to the complete list of step sizes and
        % leak-subtracted traces across series for this recording
        allLeakSub=[allLeakSub; leakSubtractCell];       
                     
    end
    
    % Pad all traces with NaNs at the end so they're the same length, for
    % ease of manipulation as an array instead of a cell.
    sweepLengths = cellfun('length',allLeakSub);
        maxLength = max(sweepLengths);
    allLeakSub=cellfun(@(x)cat(2,x,nan(1,maxLength-length(x))),allLeakSub,'UniformOutput',false);
    allLeakSub = cell2mat(allLeakSub);
   
    
% KEEP this section separate for each Analysis fxn, because what you sort
% by will change. (e.g., not by size for OnDt analysis). 

    for iStim = 1:length(allStim)
        
        % Sort by commanded size and take start/end indices of the data for each size
        switch sortSweepsBy
            case 'step magnitude'
                [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,3),1));
            case 'step position'
                [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,4),1));
            case 'ramp magnitude'
                [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,3),1));
            case 'ramp velocity'
                [sortedByParam, sortIdx] = sort(round(allStim{iStim}(:,5),1));
            case 'stim interval'
                a = cellfun(@(x) x(:,1:2),allStim,'un',0);
                stimInt = diff([a{:}],1,2)/sf; %calculate interval between stimuli in ms
                stimInt = round(stimInt(:,2:2:size(stimInt,2))/roundIntTo)*roundIntTo;
                [sortedByParam, sortIdx] = sort(stimInt(:,whichInt));
        end
        
        [eachSize,paramStartIdx,~] = unique(sortedByParam,'first');
        [~,paramEndIdx,~] = unique(sortedByParam,'last');
        sortedLeakSub = allLeakSub(sortIdx,:);
        
        nSizes = sum(~isnan(eachSize));
        nReps = paramEndIdx-paramStartIdx+1;
        
        stimMetaData = zeros(nSizes,4);
        
    %DECIDE: how to pick which stimNum is used for sorting sweeps by size/speed.
    
        % Save the stimuli sorted by size and the sorting index so you can
        % decide which is the variable stim and sort all sweeps according
        % to that index? Should this be an input parameter? Timepoint or
        % stim number? Or plot and prompt the user to pick a stim location?
        sortedStim{1,iStim} = allStim{iStim}(sortIdx,:);
        sortedStim{2,iStim} = sortIdx;
        
%NEXT: Recombine findMRCs with findRampMRCS, since the latter already uses
%stimWindow. Maybe set some input parameters for options (I vs. V, artifact
%offset, thresholding, fitting).
            

        meansByParam = NaN(nSizes, length(sortedLeakSub));

        % Use start and end indices for each step size to take the mean of the
        % leak-subtracted trace corresponding to that step size. Then smooth
        % and find peaks near the step times.
        for iSize = 1:nSizes
            sizeIdx = paramStartIdx(iSize):paramEndIdx(iSize);        
            theseSweeps = sortedLeakSub(sizeIdx,:);

            if paramEndIdx(iSize)-paramStartIdx(iSize)>0
                % meansByParam(iSize,:) = mean(sortedLeakSub(sizeIdx,:));
                meansByParam(iSize,:) = nanmean(theseSweeps);
                
            else
                % meansByParam(iSize,:) = sortedLeakSub(sizeIdx,:);
                meansByParam(iSize,:) = theseSweeps;
            end
                           
        end
        
              
        % Use the first sweep for a given size to pick the stimulus window
        % Assumes stim window doesn't change much within a given step size
        % (should be fine within a couple timepoints because findMRCs
        % should take a buffer before/after the stimWindow. Save into
        % parameters for findMRCs input.
        
        stimMetaData(:,1:2) = sortedStim{1,iStim}(paramStartIdx,1:2);
        stimMetaData(:,3) = eachSize;
        stimMetaData(:,4) = nReps;      
        
        % Find mechanoreceptor current peaks and append to seriesPeaks for
        % that stimulus number.
        mechPeaks{iCell,1} = meansByParam;
        try mechPeaks{iCell,iStim+1};
        catch
            mechPeaks{iCell,iStim+1} = [];
        end
            mechPeaks{iCell,iStim+1} = [mechPeaks{iCell,iStim+1};findMRCs(stimMetaData, meansByParam, sf, dataType)];
    end
    
    finalSortIdx = sortedStim{2,sortByStimNum};
    finalLeakSub{iCell,1} = allLeakSub(finalSortIdx,:);
% TODO: have column 2 be the stim trace used to find stim, for easy
% plotting access, and col 3 = PD trace, once you have that set up.
    
    finalStim{iCell,1} = cellName;
    for iStim = 1:length(allStim)
        finalStim{iCell,iStim+1} = sortedStim{1,iStim};
    end
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
    
end

end