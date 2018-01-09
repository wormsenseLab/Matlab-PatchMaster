%% fat vs wt, using FilterRecordings

% run IdAnalysis with sortStimByNum, sortSweepsBy 'position', separately
% for Pre/NoPre

protList ={'NoPrePulse'};
sortSweeps = {'position','position','position','position'};
matchType = 'full';
fatNoPreMRCs = IdAnalysis(ephysData,protList,fatCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1);

clear protList sortSweeps matchType

%% load prepulseFatData/get rid of empty cells
wtPreMRCs = wtPreMRCs(~cellfun('isempty',wtPreMRCs(:,1)),:);
wtNoPreMRCs = wtNoPreMRCs(~cellfun('isempty',wtNoPreMRCs(:,1)),:);

fatPreMRCs = fatPreMRCs(~cellfun('isempty',fatPreMRCs(:,1)),:);
fatNoPreMRCs = fatNoPreMRCs(~cellfun('isempty',fatNoPreMRCs(:,1)),:);
%%

%pre and nopre must have the same number of recordings
preArray = fatPreMRCs;
nopreArray = fatNoPreMRCs;

allSteps = [];

for i = 1:length(preArray)
   nopre = nopreArray{i,3};
   pre = preArray{i,4};
   
   [~, nopreIdx, preIdx] = intersect(nopre(:,1),pre(:,1));
   allSteps = [allSteps; ...
       nopre(nopreIdx,[1 3]) pre(preIdx,3)];
end

allSteps (:,4) = allSteps(:,3)./allSteps(:,2);
sortSteps = sortrows(allSteps, 1);
[eachSize, sizeStartIdx] = unique(sortSteps(:,1),'first');
[~, sizeEndIdx] = unique(sortSteps(:,1),'last');

for iSize = 1:length(sizeStartIdx)
    nRecs = length(sizeStartIdx(iSize):sizeEndIdx(iSize));
    
    ratiosBySize (iSize,1,1) = eachSize(iSize);
    ratiosBySize (iSize,2,1) = mean(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),2));
    ratiosBySize (iSize,3,1) = mean(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),3));
    ratiosBySize (iSize,4,1) = mean(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),4));
    ratiosBySize (iSize,5,1) = nRecs;

    sdBySize (iSize,2) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),2));
    sdBySize (iSize,3) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),3));
    sdBySize (iSize,4) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),4));

    semBySize (iSize,2) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),2))/sqrt(nRecs);
    semBySize (iSize,3) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),3))/sqrt(nRecs);
    semBySize (iSize,4) = std(sortSteps(sizeStartIdx(iSize):sizeEndIdx(iSize),4))/sqrt(nRecs);
    
end

clear preArray nopreArray pre nopre allSteps i iSize nRecs  sizeStartIdx sizeEndIdx preIdx nopreIdx ans 
%% renamed by genotype

% figure();
% scatter(allSteps_wt(:,1), allSteps_wt(:,4))
% hold on
% scatter(allSteps_fat(:,1), allSteps_fat(:,4),'r')



%% old code

clear stim1 stim2 stim3 meanSt1 meanSt2 meanSt3 nRecs sdSt1 sdSt2 sdSt3

a = vertcat(prepulseMRCs{:,2});
b = vertcat(prepulseMRCs{:,3});
c = vertcat(prepulseMRCs{:,4});

sizes = [a(:,1) b(:,1) c(:,1)];
sizes = sizes(17:end,:);
eachSize = unique(sizes,'rows');


for iSize = 1:length(eachSize)
    ia = all(bsxfun(@eq,sizes,eachSize(iSize,:)),2); %only returns 1 if all elements in the row match
    stim1{iSize,1} = a(ia,:);
    stim2{iSize,1} = b(ia,:);
    stim3{iSize,1} = c(ia,:);
    
    meanSt1(iSize,1) = nanmean(stim1{iSize}(:,3));
    meanSt2(iSize,1) = nanmean(stim2{iSize}(:,3));
    meanSt3(iSize,1) = nanmean(stim3{iSize}(:,3));
    
    sdSt1(iSize,1) = std(stim1{iSize}(:,3),1);
    sdSt2(iSize,1) = std(stim2{iSize}(:,3),1);
    sdSt3(iSize,1) = std(stim3{iSize}(:,3),1);

    
    
    nRecs(iSize,1) = sum(ia);
end



%% only works if all stimuli are the same series# within protocol
stim1 = cell2mat(prepulseMRCs(2:10,2));
stim1=permute(reshape(stim1',[7,16,9]),[2 1 3]);

stim2 = cell2mat(prepulseMRCs(2:10,3));
stim2=permute(reshape(stim2',[7,16,9]),[2 1 3]);

stim3 = cell2mat(prepulseMRCs(2:10,4));
stim3=permute(reshape(stim3',[7,16,9]),[2 1 3]);

meanSt1 = nanmean(stim1(:,3,:),3);
meanSt2 = nanmean(stim2(:,3,:),3);
meanSt3 = nanmean(stim3(:,3,:),3);

sdSt1 = std(stim1(:,3,:),1,3);
sdSt2 = std(stim2(:,3,:),1,3);
sdSt3 = std(stim3(:,3,:),1,3);
