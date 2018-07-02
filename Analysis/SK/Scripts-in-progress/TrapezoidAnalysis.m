%TrapezoidAnalysis.m

protList = {'TrapRate','WC_ProbeLarge','WC_Probe8'};

matchType = 'first';
strainList = {'TU2769'};
internalList = {'IC2'};

stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 80];
resistCutoff = '<250';
extFilterFreq = [2.5 5];

antTrapCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
     'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
     'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
     'stimFilterFrequencykHz', extFilterFreq);
 
 clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

 %%

ExcludeSweeps(ephysData, protList, antTrapCells, 'matchType', matchType);

clear protList matchType;

%%

protList = {'TrapRate','WC_ProbeLarge','WC_Probe8'};
sortSweeps = {'velocity','velocity','velocity','magnitude'};
matchType = 'first';
% antTrapMRCs = IdAnalysis(ephysData,protList,antTrapCells,'num','matchType',matchType, ...
%     'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
%     'recParameters', ephysMetaData,'sepByStimDistance',0, 'limitStim',2);
postTrapMRCs = IdAnalysis(ephysData,protList,postTrapCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',0, 'limitStim',2);
clear protList sortSweeps matchType

%%

%NEXT: mean anterior dist, posterior dist comparison
%THEN: concatenate MRC stats for each, see if dropping rows with nReps <2
%helps with the glut of different velocities, else just pop each point in
%individually iwth its velocity 
%if that's too much, bin and include x error bars for velocity
%plot current vs velocity for ant and post
%plot tau1, tau2, tPeak vs. velocity for ant and post
