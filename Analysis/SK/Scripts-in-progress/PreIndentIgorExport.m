% PreIndentIgorExport.m

%% Load data
% load prepulseFatData(180108).mat

% or run this section to re-load the data given ephysData and
% ephysMetaData.
strainList = {'TU2769'};
internalList = {'IC2'};
resistCutoff = '<210'; % Rs < 210 MOhm

wtCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, 'RsM', resistCutoff);

clear strainList internalList resistCutoff ans;



protList ={'NoPrePulse'};
sortSweeps = {'position','position','position','position'};
matchType = 'full';
wtNoPreMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1,...
    'sortStimBy','time','recParameters',ephysMetaData);

clear protList sortSweeps matchType

protList ={'PrePulse'};
sortSweeps = {'position','position','position','position'};
matchType = 'full';
wtPreMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1,...
    'sortStimBy','time','recParameters',ephysMetaData);

clear protList sortSweeps matchType

%% Export to Igor

% Get traces and associated recording names
noPreTraces = vertcat(wtNoPreMRCs{:,2})';
preTraces = vertcat(wtPreMRCs{:,2})';

preCells = cellfun(@(x,y) repmat(x,size(y,1),1), wtPreMRCs(:,1), wtPreMRCs(:,2), 'un',0);
preCells = cellstr(vertcat(preCells{:}))';

noPreCells = cellfun(@(x,y) repmat(x,size(y,1),1), wtNoPreMRCs(:,1), wtNoPreMRCs(:,2), 'un',0);
noPreCells = cellstr(vertcat(noPreCells{:}))';


% Get sizes of on2 step (or just step for noPre), sort by size
preSteps = vertcat(wtPreMRCs{:,4});
[preStepsSort, preI] = sortrows(preSteps,1);

preRelSizes = round(preStepsSort(:,1)) - 5;
preAbsSizes = round(preStepsSort(:,1));
preN = preStepsSort(:,7);
preCellsSort = preCells(preI);
preTracesSort = preTraces(:,preI);



noPreSteps = vertcat(wtNoPreMRCs{:,3});
[noPreStepsSort, noPreI] = sortrows(noPreSteps,1);

noPreSizes = round(noPreStepsSort(:,1));
noPreN = noPreStepsSort(:,7);
noPreCellsSort = noPreCells(noPreI);
noPreTracesSort = noPreTraces(:,noPreI);


% Create name for Igor wave for each trace
preRelSizes = cellstr(num2str(preRelSizes));
preRelSizes = cellfun(@(x) regexprep(x,'-','neg'), preRelSizes,'un',0);
preRelSizes = cellfun(@(x) regexprep(x,'\s',''), preRelSizes,'un',0);

absPreCells = cellfun(@(x,y) sprintf('Abs_5umPre_%dum_%s',x,y), num2cell(preAbsSizes'), preCellsSort, 'un',0)';
relPreCells = cellfun(@(x,y) sprintf('Rel_5umPre_%sum_%s',x,y), preRelSizes', preCellsSort, 'un',0)';
absNoPreCells = cellfun(@(x,y) sprintf('Abs_NoPre_%dum_%s',x,y), num2cell(noPreSizes'), noPreCellsSort, 'un',0)';
relNoPreCells = cellfun(@(x,y) sprintf('Rel_NoPre_%dum_%s',x,y), num2cell(noPreSizes'), noPreCellsSort, 'un',0)';


% Write to table in preparation for writing out csv.
absPreTab = array2table(preTracesSort,'VariableNames',absPreCells');
relPreTab = array2table(preTracesSort,'VariableNames',relPreCells');
absNoPreTab = array2table(noPreTracesSort,'VariableNames',absNoPreCells');
relNoPreTab = array2table(noPreTracesSort,'VariableNames',relNoPreCells');

[fname,pname] = uiputfile('*.txt','Save Tables As');

writetable(absPreTab, [pname 'absPreTraces' '.txt']);
writetable(relPreTab, [pname 'relPreTraces' '.txt']);
writetable(absNoPreTab, [pname 'absNoPreTraces' '.txt']);
writetable(relNoPreTab, [pname 'relNoPreTraces' '.txt']);


% Write stimuli too
preStim = ephysData.FAT123.data{2,12}/0.408;
noPreStim = ephysData.FAT123.data{2,15}/0.408;

preRelStimSizes = -2:7;
preRelStimSizes = cellstr(num2str(preRelStimSizes'));
preRelStimSizes = cellfun(@(x) regexprep(x,'-','neg'), preRelStimSizes,'un',0);
preRelStimSizes = cellfun(@(x) regexprep(x,'\s',''), preRelStimSizes,'un',0)';

noPreStimSizes = 3:10;

preStimHeaders = cellfun(@(x) sprintf('Stim_5umPre_%sum',x), preRelStimSizes, 'un',0);
noPreStimHeaders = arrayfun(@(x) sprintf('Stim_NoPre_%dum',x), noPreStimSizes, 'un',0);

stimPreTab = array2table(preStim,'VariableNames',preStimHeaders);
stimNoPreTab = array2table(noPreStim,'VariableNames',noPreStimHeaders);

writetable(stimPreTab, [pname 'stimPreTraces' '.txt']);
writetable(stimNoPreTab, [pname 'stimNoPreTraces' '.txt']);