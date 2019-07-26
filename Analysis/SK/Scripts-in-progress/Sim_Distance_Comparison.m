% Channel-distance simulation
% 
% Load traces from individual channel distance simulation, and
% representative traces from recordings for comparison to macroscopic
% traces, then plot.

%% Load simulation data and split it up

simPath = 'C:\Users\Sammy\Dropbox\Goodman Lab\Posters Papers Proposals\2018-Katta-SpatiotemporalDynamics\source-data\sim-channel-distance\';
simFiles = {'channel-distance-displacement.xlsx','channel-distance-speed.xlsx','channel-distance-frequency.xlsx'};
simData = cell(0);
simDataMat = cell(0);
for iPanel = 1:3
    [~,~,simData{iPanel}] = xlsread(fullfile(simPath,simFiles{iPanel}));
    simDataMat{iPanel} = cell2mat(simData{iPanel}(2:end,:));
end

recData = cell(0);
recFile = 'representative_traces.xlsx';
[~,~,recInput] = xlsread(fullfile(simPath,recFile),'Traces');

recData{1} = recInput(:,1:7); % displacement (1um, 5, 10)
recData{2} = recInput(:,8:14); % speed (106um/s, 320, 3790)
recData{3} = recInput(:,15:18); % frequency (10Hz, 100, 500)

for iPanel = 1:3
    recDataMat{iPanel} = cell2mat(recData{iPanel}(2:end,:))*1e12; %pA
end

%% Prep figures for plots

% scrsz = get(groot,'ScreenSize');
% figure('Position',[200 200 800 600]);

% Plot rotated plots for individual channel currents to look at timing
isRotate = 0;

for iPanel = 1:3
    fh(iPanel,1) = figure('Position',[100 100 650 600]); % stim and macro currents
    fh(iPanel,2) = figure('Position',[500 100 800 600]); % channel currents
    if isRotate
        fh(iPanel,3) = figure('Position',[600 0 600 800]); % rotated channel currents
    end
    switch iPanel
        case 1 % displacement
            stType = 'disp';
        case 2 % speed
            stType = 'speed';
        case 3 % frequency
            stType = 'freq';
    end
    simTime = cell2mat(simData{iPanel}(2:end,1));
    simHeaders = simData{iPanel}(1,:);
    recTime = cell2mat(recData{iPanel}(2:end,1));
    recHeaders = recData{iPanel}(1,:);

    
    % find all subsets of columns now, use overlap in indices later to pull
    % out more specific subsets
    stimCols = ~cellfun(@isempty,cellfun(@(x) regexp(x,'^stim'),simHeaders,'un',0));
    currTotCols = ~cellfun(@isempty,cellfun(@(x) regexp(x,'^currTot'),simHeaders,'un',0));
    chCurrCols = ~cellfun(@isempty,cellfun(@(x) regexp(x,'^chcurr'),simHeaders,'un',0));
    onCols{1} = ~cellfun(@isempty,cellfun(@(x) regexp(x,'_on'),simHeaders,'un',0));
    onCols{2} = ~cellfun(@isempty,cellfun(@(x) regexp(x,'_off'),simHeaders,'un',0));

    intTypes = cellfun(@(x) regexp(x,sprintf('_%s(\\S{1,5})',stType),'tokens'),simHeaders,'un',0);
    intTypes = [intTypes{:}];intTypes = [intTypes{:}];
    intTypes = unique(intTypes,'stable');
    for jIntensity = 1:length(intTypes);
        intCols{jIntensity} = ~cellfun(@isempty,cellfun(@(x) regexp(x,sprintf('_%s%s$',stType,intTypes{jIntensity})),simHeaders,'un',0));
    end    
    
    distTypes = cellfun(@(x) regexp(x,'_dist(\S{1,4})_','tokens'),simHeaders,'un',0);
    distTypes = [distTypes{:}];distTypes = [distTypes{:}];
    distTypes = unique(distTypes,'stable');
    for kDist = 1:length(distTypes);
        distCols{kDist} = ~cellfun(@isempty,cellfun(@(x) regexp(x,sprintf('_dist%s_',distTypes{kDist})),simHeaders,'un',0));
    end
    
    % find subsets in experimental recordings too
    recIntTypes = cellfun(@(x) regexp(x,sprintf('_%s(\\S{1,5})',stType),'tokens'),recHeaders,'un',0);
    recIntTypes = [recIntTypes{:}];recIntTypes = [recIntTypes{:}];
    recIntTypes = unique(recIntTypes,'stable');
    for jIntensity = 1:length(recIntTypes);
        recIntCols{jIntensity} = ~cellfun(@isempty,cellfun(@(x) regexp(x,sprintf('_%s%s$',stType,recIntTypes{jIntensity})),recHeaders,'un',0));
    end    

    recOnCols{1} = ~cellfun(@isempty,cellfun(@(x) regexp(x,'_on'),recHeaders,'un',0));
    recOnCols{2} = ~cellfun(@isempty,cellfun(@(x) regexp(x,'_off'),recHeaders,'un',0));

    % make subplots, match column subsets for each subplot, and plot
    for jIntensity = 1:length(intTypes)
        figure(fh(iPanel,1));
        
        if strcmp(stType,'disp') || strcmp(stType,'speed') % plot on and off
            % plot stimulus traces
            axh1(jIntensity,1) = subplot(3,3,1+3*(jIntensity-1));
            thisCol = sum([stimCols; intCols{jIntensity}; onCols{1}],1)==3;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'b'); % on = blue
            hold on;
            thisCol = sum([stimCols; intCols{jIntensity}; onCols{2}],1)==3;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'m'); % off = magenta
            
            % plot experimental total current
            axh1(jIntensity,2) = subplot(3,3,2+3*(jIntensity-1));
            thisCol = sum([recIntCols{jIntensity}; recOnCols{1}],1)==2;
            plot(recTime, recDataMat{iPanel}(:,thisCol),'b'); % on = blue
            hold on;
            thisCol = sum([recIntCols{jIntensity}; recOnCols{2}],1)==2;
            plot(recTime, recDataMat{iPanel}(:,thisCol),'m'); % off = magenta
            
            % plot simulated total current
            axh1(jIntensity,3) = subplot(3,3,3+3*(jIntensity-1));
            thisCol = sum([currTotCols; intCols{jIntensity}; onCols{1}],1)==3;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'b');
            hold on;
            thisCol = sum([currTotCols; intCols{jIntensity}; onCols{2}],1)==3;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'m');
            
        elseif strcmp(stType,'freq') % only one current to plot
            % plot stimulus traces
            axh1(jIntensity,1) = subplot(3,3,1+3*(jIntensity-1));
            thisCol = sum([stimCols; intCols{jIntensity}],1)==2;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'b');
            
            % plot experimental total current
            axh1(jIntensity,2) = subplot(3,3,2+3*(jIntensity-1));
            thisCol = recIntCols{jIntensity};
            plot(recTime, recDataMat{iPanel}(:,thisCol),'b'); 
            
            % plot simulated total current
            axh1(jIntensity,3) = subplot(3,3,3+3*(jIntensity-1));
            thisCol = sum([currTotCols; intCols{jIntensity}],1)==2;
            plot(simTime, simDataMat{iPanel}(:,thisCol),'b');
        end

        % Plot individual channel currents in second window
        for kDist = 1:length(distTypes)
            figure(fh(iPanel,2));
            axh2(jIntensity,kDist) = subplot(3,4,kDist+4*(jIntensity-1));
            if strcmp(stType,'disp') || strcmp(stType,'speed')
                thisCol = sum([chCurrCols; intCols{jIntensity}; onCols{1}; distCols{kDist}],1)==4;
                plot(simTime, -simDataMat{iPanel}(:,thisCol),'b'); % on = blue
                hold on;
                thisCol = sum([chCurrCols; intCols{jIntensity}; onCols{2}; distCols{kDist}],1)==4;
                plot(simTime, -simDataMat{iPanel}(:,thisCol),'m'); % off = magenta
            elseif strcmp(stType,'freq')
                thisCol = sum([chCurrCols; intCols{jIntensity}; distCols{kDist}],1)==3;
                plot(simTime, -simDataMat{iPanel}(:,thisCol),'b');
            end
            
            % Alternative plot for looking at timing vs. distance
            if isRotate == 1
                figure(fh(iPanel,3));
                axh3(kDist,jIntensity) = subplot(4,3,kDist+4*(jIntensity-1));
                if strcmp(stType,'disp') || strcmp(stType,'speed')
                    thisCol = sum([chCurrCols; intCols{jIntensity}; onCols{1}; distCols{kDist}],1)==4;
                    plot(simTime, -simDataMat{iPanel}(:,thisCol),'b'); % on = blue
                    hold on;
                    thisCol = sum([chCurrCols; intCols{jIntensity}; onCols{2}; distCols{kDist}],1)==4;
                    plot(simTime, -simDataMat{iPanel}(:,thisCol),'m'); % off = magenta
                elseif strcmp(stType,'freq')
                    thisCol = sum([chCurrCols; intCols{jIntensity}; distCols{kDist}],1)==3;
                    plot(simTime, -simDataMat{iPanel}(:,thisCol),'b');
                end
            end
        end
        
    end
    
    % set x and y limits for each set of plots
    switch iPanel
        case 1
            set(axh1,'XLim',[-0.02 0.15])
            set(axh2,'XLim',[-0.02 0.15])
            set(axh1(:,1),'YLim',[0 10])
            set(axh1(:,2:3),'YLim',[-100 10])
%             set(axh1(1,2:3),'YLim',[-15 3])
%             set(axh1(2,2:3),'YLim',[-60 5])
%             set(axh1(3,2:3),'YLim',[-150 10])
            set(axh2,'YLim',[-2 0]);
        case 2
            set(axh1,'XLim',[-0.02 0.15])
            set(axh2,'XLim',[-0.02 0.15])
            set(axh1(:,1),'YLim',[0 10])
            set(axh1(:,2:3),'YLim',[-70 10])
%             set(axh1(1,2:3),'YLim',[-20 3])
%             set(axh1(2,2:3),'YLim',[-30 5])
%             set(axh1(3,2:3),'YLim',[-70 10])
            set(axh2,'YLim',[-2 0]);
            
        case 3
            set(axh1(1,:),'XLim',[-0.05 0.25])
            set(axh1(2,:),'XLim',[0.22 0.25])
            set(axh1(3,:),'XLim',[0.244 0.25])
            set(axh2(1,:),'XLim',[-0.05 0.25])
            set(axh2(2,:),'XLim',[0.22 0.25])
            set(axh2(3,:),'XLim',[0.244 0.25])
            
            set(axh1(:,1),'YLim',[0 10])
            set(axh1(:,2:3),'YLim',[-75 0])
%             set(axh1(1,2),'YLim',[-15 0])
%             set(axh1(1,3),'YLim',[-20 -10])
%             set(axh1(2,2),'YLim',[-60 -30])
%             set(axh1(2,3),'YLim',[-65 -60])
%             set(axh1(3,2),'YLim',[-75 -65])
%             set(axh1(3,3),'YLim',[-72.2 -71])
% 
            set(axh2,'YLim',[-2 0]);
%             set(axh3(:,1),'XLim',[-0.05 0.25])
%             set(axh3(:,2),'XLim',[0.22 0.25])
%             set(axh3(:,3),'XLim',[0.244 0.25])
%             set(axh3,'YLim',[-2 0]);

    end
    clear axh1 axh2 axh3;

end

plotfixer