% ChargeExploration.m
% 
% 
% Exploring asymmetry in charge with multi-rep noise analysis traces.
% Variation in charge asymmetry with small steps, and whether charge
% actually returns to zero by the end of the on step even if there is no
% current detectable above noise.

strainList = {'TU2769'};
internalList = {'IC6'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200];
resistCutoff = '<250';
extFilterFreq = [2.5 5];

noiseTrapCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq, 'included', 1);
