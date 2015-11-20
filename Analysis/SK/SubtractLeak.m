% SubtractLeak.m
%
%
%
%
% Created by Sammy Katta on 19th November, 2015.

function leakSubtract = SubtractLeak(nSteps, data, sf, varargin)

% INPUTS
% Validate inputs, allow optional inputs for more complicated usage of the
% function, and set defaults for optional inputs not passed in.

p = inputParser;

p.addRequired('nSteps', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('data', @(x) isnumeric(x));
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);

% default: include all traces
p.addOptional('approvedTraces', 1:nSteps, @(x) isnumeric(x)); 
% default: first 30ms used as baseline to find leak current
p.addParamValue('baseLength', 30, @(x) isnumeric(x) && isscalar(x) && x>0) 

p.parse(nSteps, stimData, sf, threshold, varargin{:});

approvedTraces = p.Results.approvedTraces;
baseLength = p.Results.baseLength * sf;

% ANALYSIS
% Subtract leak for approved traces if a list was given, or for all traces
% if no list was given. Skipped traces will be left as NaNs.
leakSubtract = NaN(nSteps,length(stimComI));
        
for iStep = 1:nSteps
    
    % Analyze only desired traces within this series
    if any(approvedTraces == iStep)
        
        % Subtract leak/baseline
        leak = mean(data(1:baseLength,iStep));
        leakSubtract(:,iStep) = data(:,iStep) - leak;
        
    end
end


end