% IdAnalysis.m
% 
% Created by Sammy Katta on 20-May-2015.

function mechPeaks = IdAnalysis(ephysData, allCells)

% keyboard;

% TODO: Figure out how to get rid of bad data, OR alternatively, to read in
% good series to analyze from a spreadsheet (cross-checked against protocol
% name).
mechPeaks = cell(length(allCells),1);
sf = 5; %sampling frequency, in kHz
stepThresh = 0.05; % step detection threshold in um, could be smaller

for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    protName = 'WC_Probe';
    allSeries = find(strcmp(protName,ephysData.(cellName).protocols));
    protName = 'WC_ProbeLarge';
    allSeries = [allSeries find(strcmp(protName,ephysData.(cellName).protocols))];
    protName = 'WC_ProbeSmall';
    allSeries = [allSeries find(strcmp(protName,ephysData.(cellName).protocols))];
    
    allSteps = [];
    allOns = [];
    allOffs = [];
    
    nSeries = length(allSeries);
%     nSteps = size(ephysData.(cellName).data{2,allSeries(1)},2);
%     stepSize = zeros(nSeries,nSteps);
%     onTau = zeros(nSeries,nSteps);
%     offTau = zeros(nSeries,nSteps);
    
    for iSeries = 1:nSeries
        probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
        % convert command V to um, at 0.408 V/um
        stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
        nSteps = size(stimComI,2);
        stepSize = zeros(nSteps,1);
        onPeaks = zeros(nSteps,1);
        offPeaks = zeros(nSteps,1);

        % Find timepoint in sweep at which probe starts pushing (on) and when it
        % goes back to neutral position (off). Use that timepoint to find nearby
        % peak current
        
        for iStep = 1:nSteps
            % Figure out points when stimulus command for step starts and ends
            stepOn = stimComI(:,iStep) - mean(stimComI(1:10*sf,iStep));
            stepStart = find(stepOn > stepThresh); 
            stepLength = length(stepStart);
            if stepLength == 0
                continue
            end
            stepStart = stepStart(1);
            stepEnd = stepStart + stepLength;
            stepSize(iStep) = mean(stepOn(stepStart+1:stepEnd-1));
            
            % Get baseline for each step by grabbing the mean of the 150 points before
            % the probe displacement.
            baseProbeI_on = mean(probeI(stepStart-150+1:stepStart-1,iStep));
            baseProbeI_off = mean(probeI(stepEnd-150+1:stepEnd-1,iStep));
            
            % Find the peak current for the on step and the off step for this sweep
            % Two options for ignoring the stimulus artifact (which shows up
            % even with unpatched pipette):
            % Start looking at stepStart +2 (+1 too early, +3 hits actual peak)
            % Swap subtract order and don't use abs(onSubtract) for finding max
            
            % TODO: Make this section a separate function and call it from
            % each analysis that requires peak current.
            onSubtract = baseProbeI_on - probeI(:,iStep);
            peakOnLoc = find(onSubtract(stepStart+1:stepStart+500) == max(onSubtract(stepStart+1:stepStart+500)));
            peakOnLoc = peakOnLoc(1) + stepStart;
            
            offSubtract = baseProbeI_off - probeI(:,iStep);
            peakOffLoc = find(offSubtract(stepEnd+1:stepEnd+500) == max(offSubtract(stepEnd+2:stepEnd+500)));
            peakOffLoc = peakOffLoc(1)+stepEnd;
            
            onPeaks(iStep) = -onSubtract(peakOnLoc);
            offPeaks(iStep) = -offSubtract(peakOffLoc);
            
            [~,onFitInd] = min(abs(onSubtract(peakOnLoc:75*sf+peakOnLoc)-(onSubtract(peakOnLoc)/(2*exp(1)))));
            [~,offFitInd] = min(abs(offSubtract(peakOffLoc:75*sf+peakOffLoc)-(offSubtract(peakOffLoc)/(2*exp(1)))));
            
            onFitTime = onFitInd/sf; % seconds
            onT = 0:1/sf:onFitTime;
            offFitTime = offFitInd/sf; % seconds
            offT = 0:1/sf:offFitTime;

            
            onFit = fit(onT',onSubtract(peakOnLoc:peakOnLoc+onFitInd),'exp1');
            offFit = fit(offT',onSubtract(peakOffLoc:peakOffLoc+offFitInd),'exp1');
            %     plot(capFit,t,ICt(intStart:intStart+minInd));
            
            % Calculate time constant in seconds for calculation
            onTau(iStep) = -1/onFit.b;
            offTau(iStep) = -1/offFit.b;
            % TODO: Output taus
            % TODO: Fit with an alpha function instead of exp1 to get tau1 and tau2
        end
        
        allSteps = [allSteps; stepSize];
        allOns = [allOns; onPeaks];
        allOffs = [allOffs; offPeaks];
    end
    
    mechPeaks{iCell} = [allSteps allOns allOffs];
    
    % TODO: Figure out how to fit this to the four-parameter sigmoidal
    % function used in O'Hagan: @(X,a,b,c,d) ((a-d)/(1+((X/c)^b)))+d
    % Using optimtool? fmincon? nlinfit if you add the statistics toolbox.
    
end

end