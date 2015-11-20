% FindSteps.m
%
% OPTIONAL: approvedTraces, nStim, minStimInterval, roundedTo
% FUTURE: nStepsPerTrace? or minStepInterval?
%
% Created by Sammy Katta on 19th November, 2015.

function [stepSize, stepStarts, stepEnds] = FindSteps(nSteps, stimData, sf, threshold, varargin)

% INPUTS
% Validate inputs, allow optional inputs for more complicated usage of the
% function, and set defaults for optional inputs not passed in.
p = inputParser;

p.addRequired('nSteps', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('stimData', @(x) isnumeric(x));
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('threshold', @(x) isnumeric(x) && isscalar(x) && x>0)

p.addOptional('approvedTraces', 1:nSteps, @(x) isnumeric(x));

p.addParamValue('nStim', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParamValue('minStimInterval',300, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParamValue('roundedTo',0.1, @(x) isnumeric(x));
p.parse(nSteps, stimData, sf, threshold, varargin{:});

approvedTraces = p.Results.approvedTraces;
nStim = p.Results.nStim;
minStimInterval = p.Results.minStimInterval;
roundedTo = p.Results.roundedTo;


% ANALYSIS
% Find size of step (in um), as well as start and end indices.

% Initialize arrays as NaNs (can be placed in array). Non-included sweeps
% or sweeps that throw an error will be left as NaNs. You can later use
% nanmean to take the mean while ignoring NaN values.
stepSize = NaN(nSteps,1);
stepStarts = NaN(nSteps,1);
stepEnds = NaN(nSteps,1);

for iStep = 1:nSteps
    
    % Analyze only desired traces within this series
    if any(approvedTraces == iStep)
        
        % Figure out timepoints when stimulus command for when step starts
        % and ends, and get the size in microns.
        stepOn = stimData(:,iStep) - mean(stimData(1:10*sf,iStep));
        stepStart = find(stepOn > threshold);
        stepLength = length(stepStart);
        if stepLength == 0
            continue
        end
        stepStarts(iStep) = stepStart(1);
        stepEnds(iStep) = stepStart(1) + stepLength;
        stepSize(iStep) = mean(stepOn(stepStarts(iStep)+1:stepEnds(iStep)-1));
    end
end

% By default, round step size to nearest 0.1 to drop noise and allow 
% grouping of step sizes. roundedTo is an optional input that can be set
% larger or smaller depending on the range of step values used.
roundedTo = 1/roundedTo;
stepSize = round(stepSize*roundedTo)/roundedTo;

end