% StepPlot.m

%% Import and divide by genotype

protList = {'WC_Probe','WC_ProbeSmall','WC_ProbeLarge'};
matchType = 'full';
ExcludeSweeps(ephysData, allCells, 1, protList, matchType);

stepTracePicks = ImportMetaData(); % AllWCStepsTo104TracePicks.xls
stepTracePicks = metaDataConvert(stepTracePicks);
ephysRecordingBase = ImportMetaData(); % RecordingDatabase.xls
stepCells = unique(stepTracePicks(:,1));

genotype = cell(length(allCells),2);
for i=1:length(allCells)
    genotype(i,1) = allCells(i);
    genotype(i,2) = ephysRecordingBase(strcmp(ephysRecordingBase(:,2),allCells(i)),3);    
end

wtCells = allCells(strcmp(genotype(:,2),'TU2769'));
fatCells = allCells(strcmp(genotype(:,2),'GN381'));
wtStepCells = stepCells(ismember(stepCells,wtCells));
fatStepCells = stepCells(ismember(stepCells,fatCells));

%% Run IDAnalysis and filter empty results
mechPeaksWT = IdAnalysis(ephysData,wtStepCells,0);
mechPeaksFat = IdAnalysis(ephysData,fatStepCells,0);

mechCellsWT = allCells(~cellfun('isempty',mechPeaksWT(:,1)));
mechPeaksWT = mechPeaksWT(~cellfun('isempty',mechPeaksWT(:,1)),:);
mechCellsFat = allCells(~cellfun('isempty',mechPeaksFat(:,1)));
mechPeaksFat = mechPeaksFat(~cellfun('isempty',mechPeaksFat(:,1)),:);

%% Sort peaks and get means by step size across recordings

%TODO: Modify IDAnalysis to get PDStepSizes with mean/SD for horiz errbars

sTest = mechPeaksWT;
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

sTest = mechPeaksFat;

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