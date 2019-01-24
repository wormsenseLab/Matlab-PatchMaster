% calculate

% theseCells = {'FAT186','FAT190','FAT241'};
% protList ={'WC_Probe';'NoPre'};
protList ={'NoPre'};
% protList ={'WC_Probe'};
matchType = 'first';
sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

% use allSteps(190124).xls
[postSweeps, postStim] = IdAnalysis(ephysData,protList,posteriorDistCells,'num','matchType',matchType, ...
'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
'recParameters', ephysMetaData,'sepByStimDistance',1,'saveSweeps',1);

clear protList sortSweeps matchType

%% plot
figure('Position',[500 100 800 600]);
stepSize = 10;

for i = 1:size(postSweeps,1)
    subplot(2,4,i)
    try
        plot(postSweeps{i,3}{abs(postSweeps{i,4}(:,1)-stepSize)<0.4}')
    catch
        continue
    end
    
end

suptitle('Posterior 250ms step');
