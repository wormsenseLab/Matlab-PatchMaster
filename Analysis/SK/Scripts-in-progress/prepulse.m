stim1 = cell2mat(prepulseMRCs(2:5,2));
stim1=permute(reshape(stim1',[7,16,4]),[2 1 3]);

stim2 = cell2mat(prepulseMRCs(2:5,3));
stim2=permute(reshape(stim2',[7,16,4]),[2 1 3]);

stim3 = cell2mat(prepulseMRCs(2:5,4));
stim3=permute(reshape(stim3',[7,16,4]),[2 1 3]);

meanSt1 = nanmean(stim1(:,3,:),3);
meanSt2 = nanmean(stim2(:,3,:),3);
meanSt3 = nanmean(stim3(:,3,:),3);

sdSt1 = std(stim1(:,3,:),1,3);
sdSt2 = std(stim2(:,3,:),1,3);
sdSt3 = std(stim3(:,3,:),1,3);