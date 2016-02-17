% FindSteps.m
% This function takes the stimulus command or other stimulus data for a
% Patchmaster series, ideally pre-scaled into um, with any bias removed,
% and locates the start, end, and size of isolated steps within each sweep.
% Staircases and more complicated step patterns are not supported.
%
% Data for skipped sweeps will be filled with NaNs.
%
% EXAMPLE USAGE:
%   stepSize = FindSteps(nSweeps, stimData, sf, threshold)
%   [stepSize, stepStarts, stepEnds] = ...
%         FindSteps(nSweeps, stimData, sf, threshold, approvedTraces)
%   [stepSize, stepStarts, stepEnds] = ...
%         FindSteps(nSweeps, stimData, sf, threshold, ')
%
% INPUTS:
%   nSweeps         scalar double      Total number of sweeps in the series.
%   stimData        double array       Stimulus data, in um, with
%                                      dimensions: [sweepLength nSweeps].
%                                      For absolute size, ensure that the
%                                      baseline is set at 0.
%   sf              scalar double      Sampling frequency, in kHz.
%   threshold
%
% OUTPUTS:
%   stepSize       double vector       Vector containing the magnitude of
%                                      the step for each sweep.
%   stepStarts     double vector       Vector containing the index of the
%                                      step's onset for each sweep.
%   stepEnds       double vector       Vector containing the index of the
%                                      step's offset for each sweep.
%
% OPTIONAL INPUTS:
%   approvedTraces double vector       List of sweep numbers in the given
%                                      series to be included in the
%                                      analysis. Output for sweeps not on
%                                      this list will be filled with NaNs.
%                                      If not provided, all sweeps will be
%                                      included.
%
% OPTIONAL PARAMETER/VALUE PAIRS:
%   'nStim'         double             Number of steps to be found per
%                                      sweep. Default is 1.
%   'minStimInterval double            Minimum time gap between steps if
%                                      there are multiple within the sweep.
%                                      Default is 300 ms.
%   'roundedTo'     double             Scale for rounding step size to.
%                                      Default rounds to the nearest 0.1.
%
% Created by Sammy Katta on 19th November, 2015.

function [stepSize, stepStarts, stepEnds] = FindSteps(nSweeps, stimData, sf, threshold, varargin)

% TODO: Make stepSize/Starts/Ends able to deal with multiple steps per
% sweep.
% TODO: Code ability to find multiple steps per sweep.
% NOTE: For staircases, you could find the step size, pass it back to the
% calling fxn, subtract step size from the stim trace, and run it through
% this fxn again to find the next step up.
% TODO: Think about finding off steps that are separated from on steps.

% INPUTS
% Validate inputs, allow optional inputs for more complicated usage of the
% function, and set defaults for optional inputs not passed in.
p = inputParser;

p.addRequired('nSweeps', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('stimData', @(x) isnumeric(x));
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);
p.addRequired('threshold', @(x) isnumeric(x) && isscalar(x) && x>0)

p.addOptional('approvedTraces', 1:nSweeps, @(x) isnumeric(x));

p.addParamValue('nStim', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParamValue('minStimInterval',300, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParamValue('roundedTo',0.1, @(x) isnumeric(x));
p.parse(nSweeps, stimData, sf, threshold, varargin{:});

approvedTraces = p.Results.approvedTraces;
nStim = p.Results.nStim;
minStimInterval = p.Results.minStimInterval;
roundedTo = p.Results.roundedTo;


% ANALYSIS
% Find size of step (in um), as well as start and end indices.

% Initialize arrays as NaNs (can be placed in array). Non-included sweeps
% or sweeps that throw an error will be left as NaNs. You can later use
% nanmean to take the mean while ignoring NaN values.
stepSize = NaN(nSweeps,1);
stepStarts = NaN(nSweeps,1);
stepEnds = NaN(nSweeps,1);

for iSweep = 1:nSweeps
    
    % Analyze only desired traces within this series
    if any(approvedTraces == iSweep)
        
        % Figure out timepoints when stimulus command for when step starts
        % and ends, and get the size in microns.
        stepOn = stimData(:,iSweep) - mean(stimData(1:10*sf,iSweep));
        stepStart = find(stepOn > threshold);
        stepLength = length(stepStart);
        if stepLength == 0
            continue
        end
        stepStarts(iSweep) = stepStart(1);
        stepEnds(iSweep) = stepStart(1) + stepLength;
        stepSize(iSweep) = mean(stepOn(stepStarts(iSweep)+1:stepEnds(iSweep)-1));
    end
end

% By default, round step size to nearest 0.1 to drop noise and allow
% grouping of step sizes. roundedTo is an optional input that can be set
% larger or smaller depending on the range of step values used.
roundedTo = 1/roundedTo;
stepSize = round(stepSize*roundedTo)/roundedTo;

end