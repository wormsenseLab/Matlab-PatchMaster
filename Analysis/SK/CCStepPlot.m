% CCStepPlot.m

%% Import and divide by genotype

protList = {'Probe_CC','ProbeS_CC','ProbeL_CC'};
matchType = 'full';
ExcludeSweeps(ephysData, allCells, 1, protList, matchType);

ccStepTracePicks = ImportMetaData(); % AllWCStepsTo104TracePicks.xls
ccStepTracePicks = metaDataConvert(ccStepTracePicks);
% ephysRecordingBase = ImportMetaData(); % RecordingDatabase.xls
ccStepCells = unique(ccStepTracePicks(:,1));

genotype = cell(length(allCells),2);
for i=1:length(allCells)
    genotype(i,1) = allCells(i);
    genotype(i,2) = ephysRecordingBase(strcmp(ephysRecordingBase(:,2),allCells(i)),3);    
end

wtCells = allCells(strcmp(genotype(:,2),'TU2769'));
fatCells = allCells(strcmp(genotype(:,2),'GN381'));
wtCCStepCells = ccStepCells(ismember(ccStepCells,wtCells));
fatCCStepCells = ccStepCells(ismember(ccStepCells,fatCells));

%% Run IDAnalysis and filter empty results
ccPeaksWT = VdAnalysis(ephysData,wtCCStepCells,0);
ccPeaksFat = VdAnalysis(ephysData,fatCCStepCells,0);

mechCellsWT = allCells(~cellfun('isempty',ccPeaksWT(:,1)));
ccPeaksWT = ccPeaksWT(~cellfun('isempty',ccPeaksWT(:,1)),:);
mechCellsFat = allCells(~cellfun('isempty',ccPeaksFat(:,1)));
ccPeaksFat = ccPeaksFat(~cellfun('isempty',ccPeaksFat(:,1)),:);

%% Sort peaks and get means by step size across recordings

%TODO: Modify IDAnalysis to get PDStepSizes with mean/SD for horiz errbars

sTest = ccPeaksWT;
% sTest = mechPeaksFat;

sCat = vertcat(sTest{:,1});
[~,sizeSortIdx] = sort(sCat(:,1));
sSort = sCat(sizeSortIdx,:);

[eachSize,sizeStartIdx,~] = unique(sSort(:,1),'first');
[~,sizeEndIdx,~] = unique(sSort(:,1),'last');
nSizes = sum(~isnan(eachSize));

for iSize = 1:nSizes
sizeIdx = sizeStartIdx(iSize):sizeEndIdx(iSize);
sizeCount = sizeEndIdx(iSize)-sizeStartIdx(iSize)+1;
meansBySize(iSize,1) = nanmean(sSort(sizeIdx,3));
stdBySize(iSize,1) = nanstd(sSort(sizeIdx,3));
stErrBySize(iSize,1) = sqrt(stdBySize(iSize))/sizeCount;
end


% errorbar(eachSize,meansBySize,stErrBySize)
% errorbar(eachSize,meansBySize,stErrBySize,'r')

clear sCat  sizeSortIdx sizeStartIdx sizeEndIdx iSize nSizes sizeIdx 
clear meansBySize stdBySize stErrBySize

%% Get recording names for sorted peaks

sTest = ccPeaksFat;

sCat = vertcat(sTest{:,1}); 
sCat(sCat==40000)=20000;
[~,sizeSortIdx] = sort(sCat(:,1));
sSort = sCat(sizeSortIdx,:);

sCatTrace = vertcat(sTest{:,2});
sCatName = vertcat(sTest{:,4});
sSortTrace = sCatTrace(sizeSortIdx,:);
sSortName = sCatName(sizeSortIdx,:);

sSortTraceFat = sSortTrace;
sSortNameFat = sSortName;

%% Save mechPeaks in format for Igor's I-dCellFits

peaky = ccPeaksWT;

% StepSize FAT1on FAT2on
nCells = size(peaky,1);
onToIgor = nan(length(eachSize),nCells+1);
onToIgor(:,1) = eachSize;
offToIgor = nan(length(eachSize),nCells+1);
offToIgor(:,1) = eachSize;
colNames = cell(1,nCells+1);
colNames(1) = {'StepSize'};

for i = 1:nCells
    peakyTable = peaky{i,1};
    nSizes = size(peakyTable,1);
    for j = 1:nSizes
        onToIgor(eachSize==peakyTable(j,1),i+1)=peakyTable(j,3);
        offToIgor(eachSize==peakyTable(j,1),i+1)=peakyTable(j,4);
    end
    
    colNames {i+1} = peaky{i,4}(1,:);
end

wtOnToIgor = onToIgor;
wtOffToIgor = offToIgor;
wtColsToIgor = colNames;

%% Save toIgors as delimited text

% copy headers into Excel and save each as csv
% wtColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), wtColsToIgor(2:end),'UniformOutput',0);
% wtColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtColsToIgor(2:end),'UniformOutput',0);
% fatColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), fatColsToIgor(2:end),'UniformOutput',0);
% fatColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatColsToIgor(2:end),'UniformOutput',0);

% then for each, append data
dlmwrite('PatchData/IgorFatOffs.csv',fatOffToIgor,'-append')

%% Normalize from igor sigmoid fits
% cols = WtOnMax, WtOnXHalf, WtOnRate, WtOffMax, WtOffXHalf, WtOffRate

for i = 1:size(wtStats,1)
    wtOnNorm(:,i) = wtOnToIgor(:,i+1)/wtStats(i,1);
end

for i = 1:size(wtStats,1)
    wtOffNorm(:,i) = wtOffToIgor(:,i+1)/wtStats(i,4);
end

for i = 1:size(fatStats,1)
    fatOnNorm(:,i) = fatOnToIgor(:,i+1)/fatStats(i,1);
end

for i = 1:size(fatStats,1)
    fatOffNorm(:,i) = fatOffToIgor(:,i+1)/fatStats(i,4);
end
% 
% wtColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), wtColsToIgor(2:end),'UniformOutput',0);
% wtColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtColsToIgor(2:end),'UniformOutput',0);
% fatColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), fatColsToIgor(2:end),'UniformOutput',0);
% fatColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatColsToIgor(2:end),'UniformOutput',0);

dlmwrite('PatchData/IgorWtOnNorms.csv',wtOnNorm,'-append')
dlmwrite('PatchData/IgorWtOffNorms.csv',wtOffNorm,'-append')
dlmwrite('PatchData/IgorFatOnNorms.csv',fatOnNorm,'-append')
dlmwrite('PatchData/IgorFatOffNorms.csv',fatOffNorm,'-append')

