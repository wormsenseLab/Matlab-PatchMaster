% CCStepPlot.m

%% Import and divide by genotype

%use filtCells instead of allCells to run ExcludeSweeps
% rsFiltCells = ephysRecordingBase([ephysRecordingBase{2:89,20}]'==1,2)
% filtCells = allCells(ismember(allCells,rsFiltCells))

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

ccCellsWT = allCells(~cellfun('isempty',ccPeaksWT(:,1)));
ccPeaksWT = ccPeaksWT(~cellfun('isempty',ccPeaksWT(:,1)),:);
ccCellsFat = allCells(~cellfun('isempty',ccPeaksFat(:,1)));
ccPeaksFat = ccPeaksFat(~cellfun('isempty',ccPeaksFat(:,1)),:);

%% Sort peaks and get means by step size across recordings

%TODO: Modify IDAnalysis to get PDStepSizes with mean/SD for horiz errbars

% sTest = ccPeaksWT;
sTest = ccPeaksFat;

sCCCat = vertcat(sTest{:,1});
[~,sizeSortIdx] = sort(sCCCat(:,1));
sCCSort = sCCCat(sizeSortIdx,:);

[eachSize,sizeStartIdx,~] = unique(sCCSort(:,1),'first');
[~,sizeEndIdx,~] = unique(sCCSort(:,1),'last');
nSizes = sum(~isnan(eachSize));

for iSize = 1:nSizes
sizeIdx = sizeStartIdx(iSize):sizeEndIdx(iSize);
sizeCount = sizeEndIdx(iSize)-sizeStartIdx(iSize)+1;
meansBySize(iSize,1) = nanmean(sCCSort(sizeIdx,3));
stdBySize(iSize,1) = nanstd(sCCSort(sizeIdx,3));
stErrBySize(iSize,1) = sqrt(stdBySize(iSize))/sizeCount;
end

sCCSortFat = sCCSort;

% errorbar(eachSize,meansBySize,stErrBySize)
% errorbar(eachSize,meansBySize,stErrBySize,'r')


clear sCat sSort sizeSortIdx sizeStartIdx sizeEndIdx iSize nSizes sizeIdx 
clear meansBySize stdBySize stErrBySize

%% Get recording names for sorted peaks

sTest = ccPeaksFat;

sCCCat = vertcat(sTest{:,1}); 
sCCCat(sCCCat==40000)=20000;
[~,sizeSortIdx] = sort(sCCCat(:,1));
sCCSort = sCCCat(sizeSortIdx,:);

sCCCatTrace = vertcat(sTest{:,2});
sCCCatName = vertcat(sTest{:,4});
sCCSortTrace = sCCCatTrace(sizeSortIdx,:);
sCCSortName = sCCCatName(sizeSortIdx,:);

sCCSortTraceFat = sCCSortTrace;
sCCSortNameFat = sCCSortName;

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

wtCCOnToIgor = onToIgor;
wtCCOffToIgor = offToIgor;
wtCCColsToIgor = colNames;

%% Save toIgors as delimited text

% copy headers into Excel and save each as csv
wtCCColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), wtCCColsToIgor(2:end),'UniformOutput',0);
wtCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtCCColsToIgor(2:end),'UniformOutput',0);
fatCCColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), fatCCColsToIgor(2:end),'UniformOutput',0);
fatCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatCCColsToIgor(2:end),'UniformOutput',0);

% then for each, append data
dlmwrite('PatchData/IgorCCWtOns.csv',wtCCOnToIgor,'-append')
dlmwrite('PatchData/IgorCCWtOffs.csv',wtCCOffToIgor,'-append')
dlmwrite('PatchData/IgorCCFatOns.csv',fatCCOnToIgor,'-append')
dlmwrite('PatchData/IgorCCFatOffs.csv',fatCCOffToIgor,'-append')

%% Normalize from igor sigmoid fits
% cols = WtOnMax, WtOnXHalf, WtOnRate, WtOffMax, WtOffXHalf, WtOffRate

for i = 1:size(wtCCStats,1)
    wtCCOnNorm(:,i) = wtCCOnToIgor(:,i+1)/wtCCStats(i,1);
end

for i = 1:size(wtCCStats,1)
    wtCCOffNorm(:,i) = wtCCOffToIgor(:,i+1)/wtCCStats(i,2);
end

for i = 1:size(fatCCStats,1)
    fatCCOnNorm(:,i) = fatCCOnToIgor(:,i+1)/fatCCStats(i,1);
end

for i = 1:size(fatCCStats,1)
    fatCCOffNorm(:,i) = fatCCOffToIgor(:,i+1)/fatCCStats(i,2);
end

wtCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'off','on Norm'), wtCCColsToIgor(2:end),'UniformOutput',0);
wtCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtCCColsToIgor(2:end),'UniformOutput',0);

fatCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'off','on Norm'), fatCCColsToIgor(2:end),'UniformOutput',0);
fatCCColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatCCColsToIgor(2:end),'UniformOutput',0);

dlmwrite('PatchData/IgorCCWtOnNorms.csv',wtCCOnNorm,'-append')
dlmwrite('PatchData/IgorCCWtOffNorms.csv',wtCCOffNorm,'-append')
dlmwrite('PatchData/IgorCCFatOnNorms.csv',fatCCOnNorm,'-append')
dlmwrite('PatchData/IgorCCFatOffNorms.csv',fatCCOffNorm,'-append')

