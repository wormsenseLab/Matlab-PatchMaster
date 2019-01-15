% Channel-distance simulation
% 
% Load traces from individual channel distance simulation, and
% representative traces from recordings for comparison to macroscopic
% traces, then plot.

%% Load simulation data and split it up

simPath = 'C:\Users\Sammy\Dropbox\Goodman Lab\Posters Papers Proposals\2018-Katta-SpatiotemporalDynamics\source-data\sim-channel-distance\';
simFiles = {'channel-distance-displacement.xlsx','channel-distance-speed.xlsx','channel-distance-frequency.xlsx'};
simData = cell(0);
for i = 1:3
    [~,~,simData{i}] = xlsread(fullfile(simPath,simFiles{i}));
    simData{i} = 
end

recData = cell(0);
recFile = 'representative_traces.xlsx';
[~,~,recInput] = xlsread(fullfile(simPath,recFile));

recData{1} = recInput(:,1:3); % displacement (1um, 5, 10)
recData{2} = recInput(:,4:6); % speed (106um/s, 320, 3790)
recData{3} = recInput(:,7:9); % frequency (10Hz, 100, 500)

%% Prep figures for plots

% scrsz = get(groot,'ScreenSize');
% figure('Position',[200 200 800 600]);

for iPanel = 1:3
    fh(iPanel,1) = figure('Position',[500 100 800 600]);
    fh(iPanel,2) = figure('Position',[100 100 600 600]);
    
    switch iPanel
        case 1 % displacement
            stType = 'disp';
        case 2 % speed
            stType = 'speed';
        case 3 % frequency
            stType = 'freq';
    end
    recTraces = recData{iPanel};
    simTime = cell2mat(simData{iPanel}(2:end,1));
    simHeaders = simData{iPanel}(1,:);
    
    % find all subsets of columns now, use overlap in indices later to pull
    % out more specific subsets
    stimCols = find(~cellfun(@isempty,cellfun(@(x) regexp(x,'^stim'),simHeaders,'un',0)));
    currTotCols = find(~cellfun(@isempty,cellfun(@(x) regexp(x,'^currTot'),simHeaders,'un',0)));
    chCurrCols = find(~cellfun(@isempty,cellfun(@(x) regexp(x,'^chcurr'),simHeaders,'un',0)));
    onCols{1} = find(~cellfun(@isempty,cellfun(@(x) regexp(x,'_on'),simHeaders,'un',0)));
    onCols{2} = find(~cellfun(@isempty,cellfun(@(x) regexp(x,'_off'),simHeaders,'un',0)));

    intTypes = cellfun(@(x) regexp(x,sprintf('_%s(\\S{1,5})',stType),'tokens'),simHeaders,'un',0);
    intTypes = [intTypes{:}];intTypes = [intTypes{:}];
    intTypes = unique(intTypes,'stable');
    for jIntensity = 1:length(intTypes);
        intCols{jIntensity} = find(~cellfun(@isempty,cellfun(@(x) regexp(x,sprintf('_%s%s$',stType,intTypes{jIntensity})),simHeaders,'un',0)));
    end    
    
    distTypes = cellfun(@(x) regexp(x,'_dist(\S{1,4})_','tokens'),simHeaders,'un',0);
    distTypes = [distTypes{:}];distTypes = [distTypes{:}];
    distTypes = unique(distTypes,'stable');
    for kDist = 1:length(distTypes);
        distCols{kDist} = find(~cellfun(@isempty,cellfun(@(x) regexp(x,sprintf('_dist%s_',distTypes{kDist})),simHeaders,'un',0)));
    end    

    % NEXT: make subplots and start matching subsets for each subplot
    for jIntensity = 1:length(intTypes)
        
        for kDist = 1:length(distTypes)
            plotNo = kDist+4*(jIntensity-1);
            axh(plotNo) = subplot(3,4,plotNo);
        end
        
    end
    
    clear axh;
end