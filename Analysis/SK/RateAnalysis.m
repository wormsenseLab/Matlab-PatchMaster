% RateAnalysis.m
%
% This function 
% 
% USAGE:
% 
% 
% INPUTS:
% 
% 
% OUTPUTS:
% 

function ratePeaks = RateAnalysis(ephysData, allCells, rampStartTime, calibFlag)

% Initialize output cells and intermediate vectors, and set constants
ratePeaks = cell(length(allCells),1);
threshTime = 30; % use first n ms of trace for setting noise threshold
baseTime = 30; % use first n ms for baseline for subtracting leak
% rampStartTime = 150; % step command start time in ms

% Import list of approved sweeps for each recording+series. Can be output
% from ExcludeSweeps().
rateTracePicks = ImportMetaData();
rateTracePicks = metaDataConvert(rateTracePicks);

for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    
    allRates = [];
    allLeakSub = cell(0);
    allRamps = [];
    allOffs = [];
    traceIDs = [];
    allPDDisp = [];
    allPDSizes = [];

    
    % Given list of all cells, check which are on the approved list and use
    % those for analysis.
    %TODO: Currently exists for IdAnalysis compatibility. Decide whether to
    %drop this step if you will always give a tracePicks list.
    allSeries = matchProts(ephysData,cellName,'DispRate','MatchType','first');
    nSeries = length(allSeries);
    pickedSeries = rateTracePicks(find(strcmp(cellName,rateTracePicks(:,1))),[2,3]);
    
    for iSeries = 1:nSeries;
        
        thisSeries = allSeries(iSeries);
        pickedTraces = 0;
        
        try pickedTraces = pickedSeries{[pickedSeries{:,1}]==thisSeries,2};
        catch
            continue % if it's not on the list, go on to next series in for loop
        end
        
        stimComI = ephysData.(cellName).data{2,thisSeries}(:,pickedTraces) ./ 0.408;
        %     indentI = -ephysData.(cellName).data{3,thisSeries};
        %TODO: use PD signal instead once you have it converted to um with
        %calibration, or allow flag to pick which to use.
        probeI = ephysData.(cellName).data{1,thisSeries}(:,pickedTraces);
        sf = ephysData.(cellName).samplingFreq{thisSeries} ./ 1000;
        dataType = ephysData.(cellName).dataunit{1,thisSeries};
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
        leakSubtractCell = num2cell(leakSubtract',2);
        
        % find ramp rates
        nSweeps = size(stimComI,2);
        stimWindow = nan(nSweeps, 2);
        
        rampThresh = 1.5*thselect(stimComI(1:threshTime*sf),'rigrsure');
        stimWindow(1:end,1) = rampStartTime*sf;
        [stimSize,stimWindow(:,2)] = findSteps(nSweeps,stimComI,sf,rampThresh,'endTime',300);
        stimWindow(stimWindow(:,2) - stimWindow(:,1) == 0,2) = ...
            stimWindow(stimWindow(:,2) - stimWindow(:,1) == 0,2) + 1; % if ramp is actually step, make sure duration is at least 1 timepoint
        rampRate = stimSize ./ ((stimWindow(:,2)-stimWindow(:,1))/(sf*1000)); % rate in um/s

        
        % if calibration should be used to calculate ramp rates, get photodiode data
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
                end
                
                pdStimWindow = nan(nSweeps, 2);
                
                pdRampThresh = 1.5*thselect(photodiodeV(1:threshTime*sf),'rigrsure');
                pdStimWindow(1:end,1) = rampStartTime*sf;
                [stimSize,pdStimWindow(:,2)] = findSteps(nSweeps,photodiodeV,sf,pdRampThresh,'endTime',300);
                pdStimWindow(pdStimWindow(:,2) - pdStimWindow(:,1) == 0,2) = ...
                    pdStimWindow(pdStimWindow(:,2) - pdStimWindow(:,1) == 0,2) + 1; % if ramp is actually step, make sure duration is at least 1 timepoint
                pdRampRate = stimSize ./ ((pdStimWindow(:,2)-pdStimWindow(:,1))/(sf*1000)); % rate in um/s
                
            end
            %TODO: Change roundedTo parameter for this use
            %TODO: Check that stepThresh is applicable for measuredDisp as well
        end
                 
        % Concatenate to the complete list of step sizes and
        % leak-subtracted traces across series for this recording
        allRates = [allRates; rampRate];
        allRamps = [allRamps; stimWindow];
        allOffs = [allOffs; stimWindow(:,2)+sf*300]; % step off is 300ms after end of ramp
        allLeakSub=[allLeakSub; leakSubtractCell];       
        if calibFlag ==1
            allPDRates = [allPDRates; pdRampRate];
            allPDRamps = [allPDRamps; pdStimWindow];
            allPDDisp = [allPDDisp, measuredDisp'];
        end
        traceIDs = [traceIDs; repmat(thisSeries,size(pickedTraces))' pickedTraces'];
        
    end
    
    [sortedRates, sortIdx] = sort(allRates); 
    [eachRate,rateStartIdx,~] = unique(sortedRates,'first');
    [~,rateEndIdx,~] = unique(sortedRates,'last');
    nRates = sum(~isnan(eachRate));
    
    sortedRamps = allRamps(sortIdx,:);
    sortedOffs = allOffs(sortIdx);
    sortedLeakSub = allLeakSub(sortIdx);
    sortedIDs = traceIDs(sortIdx,:);
    
    if calibFlag == 1
        sortedPDDisp = allPDDisp(sortIdx,:);
        sortedPDRates = allPDRates(sortIdx);
        sortedPDRamps = allPDRamps(sortIdx);
    end
    
    % Use start index for the start and end times, assuming they don't
    % change within a given step size (or whatever grouping you are using;
    % should work with different step lengths/intervals as well).
    rampsByRate = sortedRamps(rateStartIdx,:);
    offsByRate = sortedOffs(rateStartIdx);
    
    meansByRate = cell(nRates,1);
    pkThresh = zeros(nRates,1);
    pkOn = zeros(nRates,1);
    pkOff = zeros(nRates,1);
    pkOnLoc = NaN(nRates,1);
    pkOffLoc = NaN(nRates,1);
    onsetTau = NaN(nRates,1);
    offsetTau = NaN(nRates,1);
    nReps = NaN(nRates,1);
    theseIDs = cell(nRates,1);
    
    if calibFlag == 1
        meanPDTrace = NaN(nRates,length(sortedPDDisp));
        meanPDRate = zeros(nRates,1);
        meanPDRamp = zeros(nRates,1);
    end
    % Use start and end indices for each step rate to take the mean of the
    % leak-subtracted trace corresponding to that step rate. Then smooth
    % and find peaks near the step times.
    

    for iRate = 1:nRates
        rateIdx = rateStartIdx(iRate):rateEndIdx(iRate);
        nReps(iRate) = length(rateIdx);
        theseSweeps = sortedLeakSub(rateIdx);
        theseIDs{iRate} = sortedIDs(rateIdx,:);
        
        % If the same rate came up in multiple protocols, pad the shorter
        % one to fit the longer one with A(numel(B))=0;
        if rateEndIdx(iRate)-rateStartIdx(iRate)>0
            try meansByRate{iRate} = mean(cell2mat(theseSweeps));
            catch
                % keyboard;
                sweepLengths = cellfun('length',theseSweeps);
                maxLength = max(sweepLengths);
                theseSweeps=cellfun(@(x)cat(2,x,nan(1,maxLength-length(x))),theseSweeps,'UniformOutput',false);
                meansByRate{iRate} = nanmean(cell2mat(theseSweeps));
            end
            if calibFlag==1
                meanPDTrace(iRate,:) = mean(sortedPDDisp(rateIdx,:));
                meanPDRate(iRate) = mean(sortedPDRates(rateIdx,:));
                meanPDRamp(iRate) = mean(sortedPDRamps(rateIdx,:));
            end
            
        else
            meansByRate{iRate} = cell2mat(sortedLeakSub(rateStartIdx(iRate):rateEndIdx(iRate),:));
            if calibFlag==1
                meanPDTrace(iRate,:) = sortedPDDisp(rateIdx,:);
                meanPDRate(iRate) = sortedPDRates(rateIdx,:);
                meanPDRamp(iRate) = sortedPDRamps(rateIdx,:);
            end
        end
                
           
        % Find MRC peaks if they exist at the onset of the step, otherwise
        % set peak amplitude as NaN. Calculate decay constant tau based on
        % single exponent fit for onset and offset currents.

        [pkOn(iRate), pkOnLoc(iRate), pkThresh(iRate), onsetTau(iRate), ~] = ...
            findRampMRCs(rampsByRate(iRate,:), meansByRate{iRate},sf);
        
        % Find MRC peaks at the offset of the step
        %TODO: check if this is working properly - should always have off
        %current with the 8um step
        
        [pkOff(iRate), pkOffLoc(iRate), pkThresh(iRate), offsetTau(iRate), ~] = ...
            findMRCs(offsByRate(iRate), meansByRate{iRate},sf,dataType);        
        
    end
    
    if calibFlag ==1
        ratePeaks{iCell,1} = [eachRate(~isnan(eachRate)) meanPDRate(~isnan(eachRate))...
            pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc meanPDRamp(~isnan(eachRate)) nReps];
        ratePeaks{iCell,2} = meansByRate;
        ratePeaks{iCell,3} = meanPDTrace;
        ratePeaks{iCell,4} = repmat(cellName,[size(pkOn),1]);
        ratePeaks{iCell,5} = theseIDs;
    else
        ratePeaks{iCell,1} = [eachRate(~isnan(eachRate)) nan(size(eachRate(~isnan(eachRate)))) ...
            pkOn pkOff onsetTau offsetTau pkOnLoc pkOffLoc nan(size(eachRate(~isnan(eachRate)))) nReps];
        ratePeaks{iCell,2} = meansByRate;
        ratePeaks{iCell,4} = repmat(cellName,[size(pkOn),1]);
        ratePeaks{iCell,5} = theseIDs;
    end
        
    
%     for iSweep = 1:nSweeps
%         [pkOn(iSweep), pkOnLoc(iSweep), pkThresh(iSweep), onsetTau(iSweep), ~] = ...
%             findRampMRCs(stimWindow(iSweep,:),leakSubtract(:,iSweep),sf);
%         
%     end
    
end





% indentI = -ephysData.FAT075.data{3,11};
% threshFraction = 0.05;

% [~,rampEnd] = findSteps(5,stimComI,5,0.1,'endTime',300);
% stepPDSize = findSteps(5,indentI,5,0.3,'endTime',300); %find this thresh based on noise, like in findMRCs
% [~,rampEnd] = findSteps(5,indentI,5,stepPDSize(5)*.5,'endTime',300);
%
% [~,rampEnd] = findSteps(5,stimComI,5,stepComSize(1)*threshFraction,'endTime',300);
% [~, ~, rampStart] = findSteps(5,stimComI,5,stepComSize(1)*threshFraction);
%