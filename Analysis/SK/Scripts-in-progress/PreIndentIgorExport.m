% PreIndentIgorExport.m

%% Load data
% Load prepulseWT_Recreated(181108).mat 
%(recreated with filtering/selection code at the current date).

% or run this section to re-load the data if ephysData and
% ephysMetaData are loaded (through FAT158).

strainList = {'TU2769'};
internalList = {'IC2'};
resistCutoff = '<210'; % Rs < 210 MOhm

wtCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, 'RsM', resistCutoff);

clear strainList internalList resistCutoff ans;

% Find MRCs for prepulse sweeps.
% Use TracePicks_Preto158_WT_IC2.xls
protList ={'PrePulse'};
sortSweeps = {'position','position','position','position'};
matchType = 'full';
wtPreMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1,...
    'sortStimBy','time','recParameters',ephysMetaData);

clear protList sortSweeps matchType

% Find MRCs for no-prepulse sweeps.
% Use TracePicks_NoPreto158_WT_IC2.xls
protList ={'NoPrePulse'};
sortSweeps = {'position','position','position','position'};
matchType = 'full';
wtNoPreMRCs = IdAnalysis(ephysData,protList,wtCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1,...
    'sortStimBy','time','recParameters',ephysMetaData);

clear protList sortSweeps matchType

% Exclude FAT116 because the pre-indentation step was 3um, rather than 5um
% like the rest of the data.
wtPreMRCs = wtPreMRCs(2:12,:); 

%% Collate data for export

% Get traces and associated recording names
noPreTraces = vertcat(wtNoPreMRCs{:,2})';
preTraces = vertcat(wtPreMRCs{:,2})';

preCells = cellfun(@(x,y) repmat(x,size(y,1),1), wtPreMRCs(:,1), wtPreMRCs(:,2), 'un',0);
preCells = cellstr(vertcat(preCells{:}))';

noPreCells = cellfun(@(x,y) repmat(x,size(y,1),1), wtNoPreMRCs(:,1), wtNoPreMRCs(:,2), 'un',0);
noPreCells = cellstr(vertcat(noPreCells{:}))';


% Get sizes of on2 step (or just step for noPre), sort by size
for i = 0:1
    preSteps = vertcat(wtPreMRCs{:,4+i});
    [preStepsSort, preI] = sortrows(preSteps,1+i);
    
    preRelSizes = round(preStepsSort(:,2));
    preAbsSizes = round(preStepsSort(:,1));
    preN = preStepsSort(:,7);
    preCellsSort = preCells(preI);
    preTracesSort = preTraces(:,preI);
    
    if i == 0
        on_PreSteps = preStepsSort(:,[1 2 6 13]);
        on_PreTraces = preTracesSort;
    else
        off_PreSteps = preStepsSort(:,[1 2 6 13]);
        off_PreTraces = preTracesSort;
    end
    
    
    % Column 4 for off
    noPreSteps = vertcat(wtNoPreMRCs{:,3+i});
    [noPreStepsSort, noPreI] = sortrows(noPreSteps,1+i);
    
    noPreSizes = round(noPreStepsSort(:,1));
    noPreN = noPreStepsSort(:,7);
    noPreCellsSort = noPreCells(noPreI);
    noPreTracesSort = noPreTraces(:,noPreI);
    
    if i == 0
        on_NoPreSteps = noPreStepsSort(:,[1 2 6 13]);
        on_NoPreTraces = noPreTracesSort;
    else
        off_NoPreSteps = noPreStepsSort(:,[1 2 6 13]);
        off_NoPreTraces = noPreTracesSort;
    end
end

% on/off_(no)PreSteps variable now contains: 
% absolute position, relative position (to 5um pre-step), current, nReps

%% Calculate summary statistics for each step size

% Relative displacement pre-step
whichSteps = on_PreSteps(:,[2,3]);

[~, sortIdx, eachDisp, dispStartIdx, dispEndIdx] = ...
    sortRowsTol(whichSteps, 0, 1);
thisMean = arrayfun(@(x,y) mean(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSD = arrayfun(@(x,y) std(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSEM = thisSD./sqrt(dispEndIdx-dispStartIdx+1); 
thisDisp = eachDisp(:,1);

relPreSummary = [thisDisp thisMean thisSD thisSEM];

% Absolute displacement pre-step
whichSteps = on_PreSteps(:,[1,3]);

[~, sortIdx, eachDisp, dispStartIdx, dispEndIdx] = ...
    sortRowsTol(whichSteps, 0, 1);
thisMean = arrayfun(@(x,y) mean(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSD = arrayfun(@(x,y) std(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSEM = thisSD./sqrt(dispEndIdx-dispStartIdx+1); 
thisDisp = eachDisp(:,1);

absPreSummary = [thisDisp thisMean thisSD thisSEM];

% No pre-step
whichSteps = off_PreSteps(:,[2,3]);
whichSteps(:,1) = -whichSteps(:,1);

[~, sortIdx, eachDisp, dispStartIdx, dispEndIdx] = ...
    sortRowsTol(whichSteps, 0, 1);
thisMean = arrayfun(@(x,y) mean(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSD = arrayfun(@(x,y) std(whichSteps(x:y,2)) ,dispStartIdx,dispEndIdx);
thisSEM = thisSD./sqrt(dispEndIdx-dispStartIdx+1); 
thisDisp = eachDisp(:,1);

noPreSummary = [thisDisp thisMean thisSD thisSEM];

relPreNorm = relPreSummary(:,1);
relPreNorm(:,2:4) = relPreSummary(:,2:4)/max(noPreSummary(:,2));
absPreNorm = absPreSummary(:,1);
absPreNorm(:,2:4) = absPreSummary(:,2:4)/max(noPreSummary(:,2));
noPreNorm = noPreSummary(:,1);
noPreNorm(:,2:4) = noPreSummary(:,2:4)/max(noPreSummary(:,2));


%% Format tables and export source data and summary stats as text files for 
% reading or for Igor import

% Turn arrays into tables and include variable names
out_Names = {'Abs_Position_um','Rel_Disp_um','Current_pA','nReps'};

out_PreSummRel = array2table(relPreSummary,'Var',{'Pre_Rel_Disp_um','Mean_Current_pA','SD','SEM'});
out_PreSummAbs = array2table(absPreSummary,'Var',{'Pre_Abs_Disp_um','Mean_Current_pA','SD','SEM'});
out_NoPreSumm = array2table(noPreSummary,'Var',{'NoPre_Disp_um','Mean_Current_pA','SD','SEM'});

out_PreNormRel = array2table(relPreNorm,'Var',{'Pre_Rel_Disp_um','Norm_Current','SD','SEM'});
out_PreNormAbs = array2table(absPreNorm,'Var',{'Pre_Abs_Disp_um','Norm_Current','SD','SEM'});
out_NoPreNorm = array2table(noPreNorm,'Var',{'NoPre_Disp_um','Norm_Current','SD','SEM'});

out_PreStepsOn = array2table(on_PreSteps,'VariableNames',out_Names);
out_PreStepsOff = array2table(off_PreSteps,'VariableNames',out_Names);
out_NoPreStepsOn = array2table(on_NoPreSteps,'VariableNames',out_Names);
out_NoPreStepsOff = array2table(off_NoPreSteps,'VariableNames',out_Names);

% Write tables to xls file
[fname, pname] = uiputfile({'*.xls','Excel file'},'Save source data as');

writetable(out_PreSummRel, fullfile(pname,fname),'Sheet','Summary','Range','A2');
writetable(out_PreSummAbs, fullfile(pname,fname),'Sheet','Summary','Range','F2');
writetable(out_NoPreSumm, fullfile(pname,fname),'Sheet','Summary','Range','K2');

writetable(out_PreNormRel, fullfile(pname,fname),'Sheet','Summary_Normalized','Range','A2');
writetable(out_PreNormAbs, fullfile(pname,fname),'Sheet','Summary_Normalized','Range','F2');
writetable(out_NoPreNorm, fullfile(pname,fname),'Sheet','Summary_Normalized','Range','K2');

writetable(out_PreStepsOn, fullfile(pname,fname),'Sheet','onPreSteps');
writetable(out_NoPreStepsOn, fullfile(pname,fname),'Sheet','onNoPreSteps');
writetable(out_PreStepsOff, fullfile(pname,fname),'Sheet','offPreSteps');
writetable(out_NoPreStepsOff, fullfile(pname,fname),'Sheet','offNoPreSteps');

%% Format wave names and export traces as text files for Igor
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
preStim = ephysData.FAT123.data{2,12}/0.408; % V to um factor
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