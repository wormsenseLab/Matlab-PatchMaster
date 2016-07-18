% Run scripts

protList = 'DispRate';
matchType = 'first';

rateSweeps = ExcludeSweeps(ephysData, allCells, protList, matchType);

rampStartTime = 150; %ms

ratePeaks = RateAnalysis(ephysData, allCells, rampStartTime);
rateCells = allCells(~cellfun('isempty',ratePeaks(:,1)));
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

wtRateCells = rateCells(ismember(rateCells,wtCells));
wtRatePeaks = ratePeaks(ismember(rateCells,wtCells));
fatRateCells = rateCells(ismember(rateCells,fatCells));
fatRatePeaks = ratePeaks(ismember(rateCells,fatCells));
%% Sort peaks and get means by rate across recordings
rTest = wtRatePeaks;
% rTest = fatRatePeaks;

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

errorbar(eachRate,meansByRate,stErrByRate)
% errorbar(eachRate,meansByRate,stErrByRate,'r')

clear rCat rSort rateSortIdx rateStartIdx rateEndIdx iRate nRates rateIdx 
clear meansByRate stdByRate stErrByRate

%TODO: RateAnalysis is finding slightly different velocities for the same
%ramps, probably starting when the sampling freq changed. Fix this by
%matching to a list of expected/commanded rates (but keep the actual
%values, for x error bars).

%% Get recording names for sorted peaks

rTest = wtRatePeaks;

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