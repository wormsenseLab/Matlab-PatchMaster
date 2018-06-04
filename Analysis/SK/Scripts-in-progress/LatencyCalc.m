whichSide = posteriorMRCs([1:4 6:10],:);
% whichSide = anteriorMRCs;

%%
clear time8;
for i = 1:size(whichSide,1)
    ind8 = whichSide{i,3}(:,1)==7.9;
    time8{i} = whichSide{i,3}(ind8,2);
end
% antTimes = time8;

%% measured by hand
antTimes = {763 4025 4025 4022 4019 4026};
postTimes = {761 764 766 763 777 767 764 765 4023};


%% hardcoded sizes/sfs/start points for steps vs nopre, check this when re-running
antStarts = [750 repmat(4000,1,5)];
antSF = [5 repmat(10,1,5)];

postStarts = [repmat(750,1,8) 4000];
postSF = [repmat(5,1,8) 10];

antLats = ([antTimes{:}] - antStarts)./antSF;
postLats = ([postTimes{:}]-postStarts)./postSF;

%%
clear trace8;
for i = 1:size(whichSide,1)
    ind8 = whichSide{i,3}(:,1)==7.9;
    try trace8{i}(2,:) = whichSide{i,2}(ind8,:);
        trace8{i}(1,:) = ((1:length(whichSide{i,2}(ind8,:)))-postTimes{i})/postSF(i);
    catch
        continue
    end
end

postTrace = trace8;
%%
clear trace8;
for i = 1:size(whichSide,1)
    ind8 = whichSide{i,3}(:,1)==7.9;
    try trace8{i}(2,:) = whichSide{i,2}(ind8,:);
        trace8{i}(1,:) = ((1:length(whichSide{i,2}(ind8,:)))-antTimes{i})/antSF(i);
    catch
        continue
    end
end

antTrace = trace8;
%% plot stuff

fant=figure(); hold on;
cellfun(@(x) plot(x(1,:),x(2,:)), antTrace);
ax = gca;
ax.ColorOrderIndex = 1;
cellfun(@(x,y,z) plot([x-y x-y],[-20e-11 5e-11]), num2cell(antStarts),antTimes);

fpost=figure(); hold on;
cellfun(@(x) plot(x(1,:),x(2,:)), postTrace);
ax = gca;
ax.ColorOrderIndex = 1;
cellfun(@(x,y,z) plot([x-y x-y],[-20e-11 5e-11]), num2cell(postStarts),postTimes);

fboth = figure();hold on;
cellfun(@(x) plot(x(1,:),x(2,:)*1e12,'r'), antTrace);
cellfun(@(x) plot(x(1,:),x(2,:)*1e12,'b'), postTrace);
xlabel('Time from peak (ms)');
ylabel('Current (pA)');
xlim([-4 6]);
ylim([-150 20]);
vline(0,'k:');

figure(); hold on;
cellfun(@(x) plot(x(2,:)), postTrace);



