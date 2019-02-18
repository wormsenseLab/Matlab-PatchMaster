% Run scripts


%use filtCells instead of allCells to run ExcludeSweeps
% rsFiltCells = ephysRecordingBase([ephysRecordingBase{2:89,20}]'==1,2);
% filtCells = allCells(ismember(allCells,rsFiltCells));

protList = {'DispRate'};
matchType = 'first';

% ExcludeSweeps(ephysData, filtCells, 1, protList, matchType);

rampStartTime = 150; %ms

ratePeaks = RateAnalysis(ephysData, filtCells, rampStartTime,0);
rateCells = filtCells(~cellfun('isempty',ratePeaks(:,1)));
ratePeaks = ratePeaks(~cellfun('isempty',ratePeaks(:,1)),:);

%%

% ephysRecordingBase = ImportMetaData();  %Recording Database
genotype = cell(length(allCells),2);
for i=1:length(allCells)
    genotype(i,1) = allCells(i);
    genotype(i,2) = ephysRecordingBase(strcmp(ephysRecordingBase(:,2),allCells(i)),3);    
end
wtCells = allCells(strcmp(genotype(:,2),'TU2769'));
fatCells = allCells(strcmp(genotype(:,2),'GN381'));

rateCellsWT = rateCells(ismember(rateCells,wtCells));
ratePeaksWT = ratePeaks(ismember(rateCells,wtCells),:);
rateCellsFat = rateCells(ismember(rateCells,fatCells));
ratePeaksFat = ratePeaks(ismember(rateCells,fatCells),:);
%% Sort peaks and get means by rate across recordings
rTest = ratePeaksWT;
% rTest = ratePeaksFat;

rCat = vertcat(rTest{:,1});
rCat(rCat==40000)=20000;
[~,rateSortIdx] = sort(rCat(:,1));
rSort = rCat(rateSortIdx,:);

[eachRate,rateStartIdx,~] = unique(rSort(:,1),'first');
[~,rateEndIdx,~] = unique(rSort(:,1),'last');
nRates = sum(~isnan(eachRate));

for iRate = 1:nRates
rateIdx = rateStartIdx(iRate):rateEndIdx(iRate);
rateCount = rateEndIdx(iRate)-rateStartIdx(iRate)+1;
meansByRate(iRate,1) = nanmean(rSort(rateIdx,3));
stdByRate(iRate,1) = nanstd(rSort(rateIdx,3));
stErrByRate(iRate,1) = sqrt(stdByRate(iRate))/rateCount;

end

% errorbar(eachRate,meansByRate,stErrByRate)
% errorbar(eachRate,meansByRate,stErrByRate,'r')

clear rCat rSort rateSortIdx rateStartIdx rateEndIdx iRate nRates rateIdx 
clear meansByRate stdByRate stErrByRate

%TODO: RateAnalysis is finding slightly different velocities for the same
%ramps, probably starting when the sampling freq changed. Fix this by
%matching to a list of expected/commanded rates (but keep the actual
%values, for x error bars).

%% Get recording names for sorted peaks

rTest = ratePeaksWT;
% rTest = ratePeaksFat;

rCat = vertcat(rTest{:,1}); 
rCat(rCat==40000)=20000;
[~,rateSortIdx] = sort(rCat(:,1));
rSort = rCat(rateSortIdx,:);

rCatTrace = vertcat(rTest{:,2});
rCatName = vertcat(rTest{:,4});
maxSize = max(cellfun(@numel,rCatTrace));
catFcn = @(x) [x nan(1,maxSize-numel(x))];
rMatTrace = cellfun(catFcn,rCatTrace,'UniformOutput',false);
rMatTrace = vertcat(rMatTrace{:});
rSortTrace = rMatTrace(rateSortIdx,:);
rSortName = rCatName(rateSortIdx,:);

rSortTraceWT = rSortTrace;
rSortNameWT = rSortName;

% rSortTraceFat = rSortTrace;
% rSortNameFat = rSortName;

rSortName = cellstr(rSortName)';

dlmwrite('PatchData/RateTracesWT.csv',rSortTraceWT', '-append');
% dlmwrite('PatchData/RateTracesFat.csv',rSortTraceFat', '-append');

%% Save peaks in Igor-friendly format

% peaky = ratePeaksWT;
peaky = ratePeaksFat;

% StepSize FAT1on FAT2on
nCells = size(peaky,1);
onToIgor = nan(length(eachRate),nCells+1);
onToIgor(:,1) = eachRate;
offToIgor = nan(length(eachRate),nCells+1);
offToIgor(:,1) = eachRate;
colNames = cell(1,nCells+1);
colNames(1) = {'RampRate'};

for i = 1:nCells
    peakyTable = peaky{i,1};
    nSizes = size(peakyTable,1);
    for j = 1:nSizes
        onToIgor(eachRate==peakyTable(j,1),i+1)=peakyTable(j,3);
        offToIgor(eachRate==peakyTable(j,1),i+1)=peakyTable(j,4);
    end
    
    colNames {i+1} = peaky{i,4}(1,:);
end

% Can't really do this in Igor wo more work b/c sampling freq doesn't get included
% wtRateOnToIgor = onToIgor;
% wtRateOffToIgor = offToIgor;
% wtRateColsToIgor = colNames;
% 
% fatRateOnToIgor = onToIgor;
% fatRateOffToIgor = offToIgor;
% fatRateColsToIgor = colNames;


%% Save toIgors as delimited text

% copy headers into Excel and save each as csv
wtRateColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), wtRateColsToIgor(2:end),'UniformOutput',0);
wtRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtRateColsToIgor(2:end),'UniformOutput',0);
fatRateColsToIgor(2:end) = cellfun(@(x) horzcat(x,' on'), fatRateColsToIgor(2:end),'UniformOutput',0);
fatRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatRateColsToIgor(2:end),'UniformOutput',0);

% then for each, append data
dlmwrite('PatchData/IgorRateWtOns.csv',wtRateOnToIgor,'-append')
dlmwrite('PatchData/IgorRateWtOffs.csv',wtRateOffToIgor,'-append')
dlmwrite('PatchData/IgorRateFatOns.csv',fatRateOnToIgor,'-append')
dlmwrite('PatchData/IgorRateFatOffs.csv',fatRateOffToIgor,'-append')

%% Normalize from igor sigmoid fits
% cols = WtOnMax, WtOnXHalf, WtOnRate, WtOffMax, WtOffXHalf, WtOffRate

for i = 1:size(wtRateStats,1)
    wtRateOnNorm(:,i) = wtRateOnToIgor(:,i+1)/wtRateStats(i,1);
end

% for i = 1:size(wtRateStats,1)
%     wtRateOffNorm(:,i) = wtRateOffToIgor(:,i+1)/wtRateStats(i,2);
% end

for i = 1:size(fatRateStats,1)
    fatRateOnNorm(:,i) = fatRateOnToIgor(:,i+1)/fatRateStats(i,1);
end
% 
% for i = 1:size(fatRateStats,1)
%     fatRateOffNorm(:,i) = fatRateOffToIgor(:,i+1)/fatRateStats(i,2);
% end

wtRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'off','on Norm'), wtRateColsToIgor(2:end),'UniformOutput',0);
% wtRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), wtRateColsToIgor(2:end),'UniformOutput',0);

fatRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'off','on Norm'), fatRateColsToIgor(2:end),'UniformOutput',0);
% fatRateColsToIgor(2:end) = cellfun(@(x) strrep(x,'on','off'), fatRateColsToIgor(2:end),'UniformOutput',0);

dlmwrite('PatchData/IgorRateWtOnNorms.csv',wtRateOnNorm,'-append')
% dlmwrite('PatchData/IgorRateWtOffNorms.csv',wtRateOffNorm,'-append')
dlmwrite('PatchData/IgorRateFatOnNorms.csv',fatRateOnNorm,'-append')
% dlmwrite('PatchData/IgorRateFatOffNorms.csv',fatRateOffNorm,'-append')


%% Pull out stimcom of FAT089 Trace


