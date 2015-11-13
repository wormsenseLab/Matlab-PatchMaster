% OnDtAnalysis.m
% 
% 

function dtPeaks = OnDtAnalysis(ephysData, allCells)

sf = 5; %sampling frequency, in kHz
stepThresh = 0.05; % step detection threshold in um, could be smaller
% dtPeaks = cell(length(allCells),1);

for iCell = 1:length(allCells)
    cellName = allCells{iCell};
    protName = 'WC_Probe_OnDt';
    allSeries = find(strcmp(protName,ephysData.(cellName).protocols));
    
    allOn1s = [];
    allOn2s = [];    
    allOffs = [];
    allDts = [];
    
    nSeries = length(allSeries);
    
    for iSeries = 1:nSeries
        probeI = ephysData.(cellName).data{1,allSeries(iSeries)};
        % convert command V to um, at 0.408 V/um for current setup(5-15-15)
        stimComI = ephysData.(cellName).data{2,allSeries(iSeries)} ./ 0.408;
        
        nSteps = size(stimComI,2);
        
        stepLength = zeros(nSteps,1);
        on1Peaks = zeros(nSteps,1);        
        on2Peaks = zeros(nSteps,1);
        offPeaks = zeros(nSteps,1);

        for iStep = 1:nSteps
            % Figure out points when stimulus command for each step starts and ends
            % Assuming step doesn't start within first 10 ms of sweep
            stepOn = stimComI(:,iStep) - mean(stimComI(1:10*sf,iStep));
            step1Start = find(stepOn > stepThresh);
            stepOn = stepOn - stepOn(step1Start(1)+1);
            step2Start = find(stepOn > stepThresh);
            step2Length = length(step2Start);
            if step2Length == 0
                continue
            end
            step1Length = length(step1Start)-step2Length;            
            step1Start = step1Start(1);
            step2Start = step2Start(1);
            stepEnd = step2Start+step2Length;
            
            % Get baseline for each step by grabbing the mean of the 150 points before
            % the probe displacement.
            baseProbeI_on1 = mean(probeI(step1Start-150+1:step1Start-1,iStep));
            baseProbeI_on2 = mean(probeI(step2Start-150+1:step2Start-1,iStep));
            baseProbeI_off = mean(probeI(stepEnd-150+1:stepEnd-1,iStep));

            on1Subtract = baseProbeI_on1 - probeI(:,iStep);
            peakOn1Loc = find(on1Subtract(step1Start+1:step1Start+step1Length-1) ...
                == max(on1Subtract(step1Start+1:step1Start+step1Length-1)));
            peakOn1Loc = peakOn1Loc(1) + step1Start;
            
            on2Subtract = baseProbeI_on2 - probeI(:,iStep);
            peakOn2Loc = find(on2Subtract(step2Start+1:step2Start+step2Length-1) ...
                == max(on2Subtract(step2Start+1:step2Start+step2Length-1)));
            peakOn2Loc = peakOn2Loc(1) + step2Start;
            
            offSubtract = baseProbeI_off - probeI(:,iStep);
            peakOffLoc = find(offSubtract(stepEnd+1:end) == max(offSubtract(stepEnd+2:end)));
            peakOffLoc = peakOffLoc(1)+stepEnd;

            on1Peaks(iStep) = -on1Subtract(peakOn1Loc);
            on2Peaks(iStep) = -on2Subtract(peakOn2Loc);            
            offPeaks(iStep) = -offSubtract(peakOffLoc);
            stepLength(iStep) = ceil(step1Length/sf); % dt in ms
        end
        
        allDts = [allDts; stepLength];
        allOn1s = [allOn1s; on1Peaks];
        allOn2s = [allOn2s; on2Peaks];
        allOffs = [allOffs; offPeaks];
    end
    cellPeaks = [allDts allOn1s allOn2s allOffs];
    whichSteps = unique(allDts);
        dtPeaks.(cellName).dts = whichSteps';
        for i=1:length(whichSteps)
            testPeak = find(cellPeaks(:,1)==whichSteps(i));
            dtPeaks.(cellName).on1(:,i) = cellPeaks(testPeak,2);
            dtPeaks.(cellName).on2(:,i) = cellPeaks(testPeak,3);
            dtPeaks.(cellName).off(:,i) = cellPeaks(testPeak,4);           
        end
    end
end
