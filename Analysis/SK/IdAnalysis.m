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
% TODO: add flag for using stim Time vs. stim # (in case of protocols where
% stimuli may not be detected in certain sweeps, throwing off the # of stim
% in that sweep).

function [mechPeaks, sortedSizes, sortedLeakSub] = IdAnalysis(ephysData, allCells, calibFlag)

% keyboard;

mechPeaks = cell(length(allCells),5);
stepThresh = 0.05; % step detection threshold in um, could be smaller
baseTime = 30; % length of time (ms) to use as immediate pre-stimulus baseline
smoothWindow = 5; % n timepoints for moving average window for findPeaks
stimConversionFactor = 0.408; % convert command V to um, usually at 0.408 V/um
sortStimBy = 'num';

% Load and format Excel file with lists (col1 = cell name, col2 = series number,
% col 3 = comma separated list of good traces for analysis)
mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);

% Find applicable series and check against list of included series/traces
% (this allows a cross-check on the protocol name) before analyzing
% Values for traces not on the list will be stored as NaN.
for iCell = 1:length(allCells)
    
    allStim = cell(1);
    allPDSizes = [];
    allLeakSub = [];
    allPDDisp = [];
    traceIDs = [];
    
    % Given list of all cells, check which are on the approved list and use
    % those for analysis. Conversely, make sure the cells on the list match
    % the expected protocol type.
    cellName = allCells{iCell};
    allSeries = matchProts(ephysData,cellName,...
        {'WC_Probe','WC_ProbeLarge','WC_ProbeSmall'},'MatchType','full');
    nSeries = length(allSeries);
    pickedSeries = mechTracePicks(find(strcmp(cellName,mechTracePicks(:,1))),[2,3]);
    
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
        
        % if calibration should be used to calculate step sizes, get photodiode data
        % and use calib curve to transform photodiodeV trace into
        % measured displacement trace (then use same findSteps threshold)
        if calibFlag == 1
            try pdCalib = ephysData.(cellName).calibration;
            catch
                fprintf('No calibration found for %s\n', cellName);
                calibFlag = 2;
            end
            
            if calibFlag ==1
                
                photodiodeV = ephysData.(cellName).data{3,thisSeries}(:,pickedTraces);
                
                if isempty(photodiodeV)
                    photodiodeV = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces);
                end
                
                
                % Interpolate photodiode voltage to calculate measured disp
                try measuredDisp = interp1(-pdCalib(2,:), pdCalib(1,:), -photodiodeV, 'linear','extrap');
                catch
                    fprintf('Interpolation failed for %s, series %d\n',cellName,thisSeries);
                    calibFlag = 2;
                end
                %TODO: swap this for newStepFind as well
                try [pdStepSize, pdStepStarts, pdStepEnds] = ...
                        findSteps(nSweeps, measuredDisp, sf, stepThresh, 'roundedTo', 0.05);
                catch
                end
                
            end
            %TODO: Change roundedTo parameter for this use
            %TODO: Check that stepThresh is applicable for measuredDisp as well
        end
        
        
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        
        seriesStimuli = ...
            newStepFind(nSweeps, stimComI, sf, 'scaleFactor', stimConversionFactor);

        % add series number to seriesStimuli for referencing back from 
        % analyzed data to a particular trace
        seriesStimuli (:,8) = repmat(thisSeries,size(seriesStimuli,1),1);
        
        % depending on input param, use stimulus number or stimulus
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
            try allStim{iStim};
            catch
                allStim{iStim} = [];
            end
            allStim{iStim} = [allStim{iStim};stimByNum{iStim}];
        end
        
        % Concatenate to the complete list of step sizes and
        % leak-subtracted traces across series for this recording
        allLeakSub = [allLeakSub; leakSubtract'];
        
        
        
        if calibFlag == 1
            allPDSizes = [allPDSizes; pdStepSize];
            allPDDisp = [allPDDisp; measuredDisp'];
        end
        
    end
    
    
    %NEXT: make the sorting work with each stim for allStim
    
    %
    %     % Sort by commanded size and take start/end indices of the data for each size
    %     [sortedSizes, sortIdx] = sort(allSizes);
    %     [eachSize,sizeStartIdx,~] = unique(sortedSizes,'first');
    %     [~,sizeEndIdx,~] = unique(sortedSizes,'last');
    %     nSizes = sum(~isnan(eachSize));
    %
    %     sortedStarts = allStarts(sortIdx);
    %     sortedEnds = allEnds(sortIdx);
    %     sortedLeakSub = allLeakSub(sortIdx,:);
    %     sortedIDs = traceIDs(sortIdx,:);
    %
    %     if calibFlag == 1
    %         sortedPDDisp = allPDDisp(sortIdx,:);
    %         sortedPDSizes = allPDSizes(sortIdx);
    %     end
    %
    %     % TODO: Store nReps as endIdx-StartIdx for each step size
    %
    %     % Use start index for the start and end times, assuming they don't
    %     % change within a given step size (or whatever grouping you are using;
    %     % should work with different step lengths/intervals as well).
    %     startsBySize = sortedStarts(sizeStartIdx);
    %     endsBySize = sortedEnds(sizeStartIdx);
    %
    %     meansBySize = NaN(nSizes,length(sortedLeakSub));
    %     pkThresh = zeros(nSizes,1);
    %     pkOn = zeros(nSizes,1);
    %     pkOff = zeros(nSizes,1);
    %     pkOnLoc = NaN(nSizes,1);
    %     pkOffLoc = NaN(nSizes,1);
    %     onsetTau = NaN(nSizes,1);
    %     offsetTau = NaN(nSizes,1);
    %     nReps = NaN(nSizes,1);
    %     theseIDs = cell(nSizes,1);
    %
    %     if calibFlag == 1
    %         meanPDTrace = NaN(nSizes,length(sortedPDDisp));
    %         meanPDSize = zeros(nSizes,1);
    %     end
    %
    %     % Use start and end indices for each step size to take the mean of the
    %     % leak-subtracted trace corresponding to that step size. Then smooth
    %     % and find peaks near the step times.
    %     for iSize = 1:nSizes
    %         sizeIdx = sizeStartIdx(iSize):sizeEndIdx(iSize);
    %         theseIDs{iSize} = sortedIDs(sizeIdx,:);
    %         nReps(iSize) = length(sizeIdx);
    %
    %         if sizeEndIdx(iSize)-sizeStartIdx(iSize)>0
    %             meansBySize(iSize,:) = mean(sortedLeakSub(sizeIdx,:));
    %             if calibFlag==1
    %                 meanPDTrace(iSize,:) = mean(sortedPDDisp(sizeIdx,:));
    %                 meanPDSize(iSize) = mean(sortedPDSizes(sizeIdx,:));
    %             end
    %
    %         else
    %             meansBySize(iSize,:) = sortedLeakSub(sizeIdx,:);
    %             if calibFlag==1
    %                 meanPDTrace(iSize,:) = sortedPDDisp(sizeIdx,:);
    %                 meanPDSize(iSize) = sortedPDSizes(sizeIdx,:);
    %             end
    %         end
    %
    %
    %         % Find MRC peaks if they exist at the onset of the step, otherwise
    %         % set peak amplitude as NaN. Calculate decay constant tau based on
    %         % single exponent fit for onset and offset currents.
    %
    %         [pkOn(iSize), pkOnLoc(iSize), pkThresh(iSize), onsetTau(iSize), ~] = ...
    %             findMRCs(startsBySize(iSize), meansBySize(iSize,:),sf, dataType);
    %
    %         % Find MRC peaks at the offset of the step
    %
    %         [pkOff(iSize), pkOffLoc(iSize), pkThresh(iSize), offsetTau(iSize), ~] = ...
    %             findMRCs(endsBySize(iSize), meansBySize(iSize,:),sf, dataType);
    %
    %     end
    %
    %
    %     if calibFlag==1
    %         mechPeaks{iCell,1} = [eachSize(~isnan(eachSize)) meanPDSize(~isnan(eachSize)) ...
    %             pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc nReps];
    %         mechPeaks{iCell,2} = meansBySize;
    %         mechPeaks{iCell,3} = meanPDTrace;
    %         mechPeaks{iCell,4} = repmat(cellName,[size(pkOn),1]);
    %         mechPeaks{iCell,5} = theseIDs;
    %     else
    %         mechPeaks{iCell,1} = ...
    %             [eachSize(~isnan(eachSize)) nan(size(eachSize(~isnan(eachSize))))...
    %             pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc nReps];
    %         mechPeaks{iCell,2} = meansBySize;
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