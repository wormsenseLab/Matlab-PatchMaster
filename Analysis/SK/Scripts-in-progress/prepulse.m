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