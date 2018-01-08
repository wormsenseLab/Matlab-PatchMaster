%SineAdaptAnalysis.m
protList = {'Sine200_time'};
matchType = 'full';
baseTime = 30;

mechTracePicks = ImportMetaData();
mechTracePicks = metaDataConvert(mechTracePicks);
sineCells=unique(mechTracePicks(:,1));
%%
for iCell = 1:length(sineCells)
    allLeakSub = [];
    
    cellName = sineCells{iCell};
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
        
        leakSubtract = ...
            SubtractLeak(probeI, sf, 'BaseLength', baseTime);
%         leakSubtractCell = num2cell(leakSubtract',2);      
        
        allLeakSub=[allLeakSub; leakSubtract'];

    end
    
    sineMeans(iCell,:) = mean(allLeakSub,1);
    nReps(iCell) = size(allLeakSub,1);
end

%%
sine200StepMeans = sine200Means(:,[1:7000 18000:26000]);
sine500StepMeans = sine500Means(:,[1:7000 18000:26000]);

sineParams(:,1)=[1501;4501;9501;12501];
sineParams(:,2)=[1505;4505;9505;12505];
sineParams(:,3)=[8;8;8;8];
sineParams(:,4)=[7;7;7;7];

sf=10;
dataType='A';
tauType='thalfmax';

sineAdPeaks = cell(length(sine200Cells),1);
sine200Peaks = cell(4,1);

for iCell=1:length(sine200Cells)
    for iStim =1:4
        sineAdPeaks{iCell,iStim} = findMRCs(sineParams(iStim,:), sine200StepMeans(iCell,:), sf, dataType, ...
            'tauType', tauType);
    end
end

for iStim = 1:4
    sine200Peaks{iStim} = vertcat(sineAdPeaks{:,iStim});
    sine200Peaks{iStim}(:,7) = nReps200;
end

sineAdPeaks = cell(length(sine500Cells),1);
sine500Peaks = cell(4,1);

for iCell=1:length(sine500Cells)
    for iStim =1:4
        sineAdPeaks{iCell,iStim} = findMRCs(sineParams(iStim,:), sine500StepMeans(iCell,:), sf, dataType, ...
            'tauType', tauType);
    end
end

for iStim = 1:4
    sine500Peaks{iStim} = vertcat(sineAdPeaks{:,iStim});
    sine500Peaks{iStim}(:,7) = nReps500;
end
%%

sine200OnRatio = sine200Peaks{3}(:,3)./sine200Peaks{1}(:,3);
sine500OnRatio = sine500Peaks{3}(:,3)./sine500Peaks{1}(:,3);
mean(sine200OnRatio)
mean(sine500OnRatio)

%% Control-ish: Use ISI_3s and take average of consecutive ISIs

% load ISItest(2-16-16).mat
diffs = [];
for i=1:5
diffs = rdiff(isi3s{1,i}(:,1));
diffs = diffs(diffs>0);
isi3s_ratio(i) = mean(diffs);
end

diffs = [];
for i=1:5
diffs = rdiff(isi1s{1,i}(:,1));
diffs = diffs(diffs>0);
isi1s_ratio(i) = mean(diffs);
end

%% bar

barLabels = { 'None (1s)', 'None (3s)', '200 Hz (1.6s)', '500 Hz (1.6s)'};

barMeans(1) = mean(isi1s_ratio);
barMeans(2) = mean(isi3s_ratio);
barMeans(3) = mean(sine200OnRatio);
barMeans(4) = mean(sine500OnRatio);

barSEM(1) = std(isi1s_ratio)/sqrt(length(isi1s_ratio));
barSEM(2) = std(isi3s_ratio)/sqrt(length(isi3s_ratio));
barSEM(3) = std(sine200OnRatio)/sqrt(length(sine200OnRatio));
barSEM(4) = std(sine500OnRatio)/sqrt(length(sine500OnRatio));

barwitherr(barSEM, barMeans);    % Plot with errorbars

set(gca,'XTickLabel',barLabels)
ylabel('On2/On1 Current Ratio')
