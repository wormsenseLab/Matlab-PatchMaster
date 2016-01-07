% IVRsCorrection.m
% 
% Document

function inIV = IVRsCorrection(ephysData, inIV, names)
% TODO: Currently last 15ms of step
steadyStateTime = 400:550;


% Based on IVq protocol as of 18-May-2015:
% IVq	(5kHz, 200?s per point, sweep interval = 400ms)
% time  (ms)	10          100         100
% V     (mV)	-60     -110 ++ 20      -60


if iscell(inIV)
    if length(inIV) == 1
        inIV = inIV{:};
    else
        error('TooManyStructs', 'inIV is a cell containing multiple structs. Please input a 1x1 struct');
    end
end


% Find actual voltage for each step at steady state by correcting from
% command voltage
for iCell = 1:length(names)
    cellName = names{iCell};
    
    for iStep = 1:12
        commandV = (-110 + 20*(iStep-1))/1000; % in V
        protRs = ephysData.(cellName).protRs;        
        rs = ephysData.(cellName).Rs(protRs)*1E6; % Rs for that WC_IVq in Ohms
        inIV.(cellName).actualV(iStep) = ...
            commandV - mean(inIV.(cellName).capCorrIV(steadyStateTime,iStep)*rs);
        inIV.(cellName).meanI(iStep) = mean(inIV.(cellName).capCorrIV(steadyStateTime,iStep));
    end
    
end

end