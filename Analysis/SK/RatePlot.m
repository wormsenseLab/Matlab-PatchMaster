% Run scripts

protList = 'DispRate';
matchType = 'first';
rampStartTime = 150; %ms

% rateSweeps = ExcludeSweeps(ephysData, newCells, protList, matchType);

ratePeaks = RateAnalysis(ephysData, allCells, rampStartTime);
rateCells = allCells(~cellfun('isempty',ratePeaks(:,1)));
ratePeaks = ratePeaks(~cellfun('isempty',ratePeaks(:,1)),:);

ratePeaks(:,3)=cellfun(@(x,y) repmat(x,size(y,1),1), rateCells, ratePeaks(:,2), 'UniformOutput',0);
%%
genotype = ImportMetaData();
strcmp(genotype(:,1),rateCells)
genotype=genotype(:,2);
wtRateCells = rateCells(strcmp(genotype,'TU2769'));
wtRatePeaks = ratePeaks(strcmp(genotype,'TU2769'),:);
fatRateCells = rateCells(strcmp(genotype,'GN381'));
fatRatePeaks = ratePeaks(strcmp(genotype,'GN381'),:);
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

% errorbar(eachRate,meansByRate,stErrByRate)
errorbar(eachRate,meansByRate,stErrByRate,'r')

clear rCat rSort rateSortIdx rateStartIdx rateEndIdx iRate nRates rateIdx 

%% Get recording names for sorted peaks

rTest = wtRatePeaks;

rCat = vertcat(rTest{:,1}); 
rCat(rCat==40000)=20000;
[~,rateSortIdx] = sort(rCat(:,1));
rSort = rCat(rateSortIdx,:);

rCatTrace = vertcat(rTest{:,2});
rCatName = vertcat(rTest{:,3});
maxSize = max(cellfun(@numel,rCatTrace));
catFcn = @(x) [x nan(1,maxSize-numel(x))];
rMatTrace = cellfun(catFcn,rCatTrace,'UniformOutput',false);
rMatTrace = vertcat(rMatTrace{:});
rSortTrace = rMatTrace(rateSortIdx,:);
rSortName = rCatName(rateSortIdx,:);