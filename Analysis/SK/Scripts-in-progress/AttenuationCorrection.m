%% Testing stuff
a = cellfun(@(x) x(x(:,1)==9.9,[11 6]), anteriorMRCs(:,3), 'un', 0);
a=vertcat(a{:});

ant = a;

% 

%% Make list of approved traces (by selecting traces to exclude)

% protList = {'DispRate'};
% protList = {'PrePulse'};
% protList = {'WC_Probe';'WC_ProbeSmall';'WC_ProbeLarge'};
protList ={'WC_Probe';'NoPre'};

matchType = 'first';
strainList = {'TU2769'};
internalList = {'IC2'};
% cellTypeList = {'ALMR'};
stimPosition = {'anterior'};

wormPrep = {'dissected'};
cellDist = [40 200];
resistCutoff = '<250';
extFilterFreq = 2.5;

anteriorDistCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq);

stimPosition = {'posterior'};

posteriorDistCells = FilterRecordings(ephysData, ephysMetaData,...
    'strain', strainList, 'internal', internalList, ...
    'stimLocation', stimPosition, 'wormPrep', wormPrep, ...
    'cellStimDistUm',cellDist, 'RsM', resistCutoff, ...
    'stimFilterFrequencykHz', extFilterFreq);

clear cellDist strainList internalList cellTypeList stimPosition resistCutoff ans wormPrep;

%% Select sweeps

ExcludeSweeps(ephysData, protList, anteriorDistCells, 'matchType', matchType);
ExcludeSweeps(ephysData, protList, posteriorDistCells, 'matchType', matchType);

%% Find MRCs
sortSweeps = {'magnitude','magnitude','magnitude','magnitude'};

anteriorMRCs = IdAnalysis(ephysData,protList,anteriorDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1);

posteriorMRCs = IdAnalysis(ephysData,protList,posteriorDistCells,'num','matchType',matchType, ...
    'tauType','thalfmax', 'sortSweepsBy', sortSweeps, 'integrateCurrent',1 , ...
    'recParameters', ephysMetaData,'sepByStimDistance',1);

clear protList sortSweeps matchType

%% Pull in distance
whichMRCs = posteriorMRCs;
distVPeak = [];

thisAtt = attenuationData(:,[2 10]);

for iCell = 1:size(whichMRCs,1)
    thisCell = whichMRCs{iCell,3};
    whichStep = round(thisCell(:,1)) == 5;
    if any(whichStep)
        thisName = whichMRCs{iCell,1};
        hasAtt = strcmp(thisName,thisAtt(:,1));
        
        if any(hasAtt)
            distVPeak(iCell,:) = [thisCell(whichStep,[11 6]) thisAtt{hasAtt,2}];
        else
            distVPeak(iCell,:) = [thisCell(whichStep,[11 6]) nan];
        end
        
    end
    
end

clear a iCell thisCell whichMRCs whichStep thisName hasAtt

% Now calculate what the actual voltage was at the stimulus site based on
% the voltage attenuation factor. In this case, command voltage was -60mV
% for all steps. 
% 
% attenuation factor = Vm'/Vc so Vm' = Vc * attenuation factor
% 
% Then calculate what current would've been based on sodium reversal
% potential. 
% For IC2 solution used here, Ena = +94mV.
% 
% So Im' = Im * ((Vc-Ena)/(Vm'-Ena))
% 
Vc = -0.06; %in V
Ena = 0.094; % in V
Im = distVPeak(:,2);
Vatt = distVPeak(:,3);

Vm = Vc * Vatt;
distVPeak(:,4) = (Im * (Vc-Ena)) ./ (Vm-Ena);

%%

distVPeak_Post_5 = distVPeak;
