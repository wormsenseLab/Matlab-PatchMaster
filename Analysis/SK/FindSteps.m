% findSteps.m
% This function takes the stimulus command or other stimulus data for a
% Patchmaster series, ideally pre-scaled into um, with any bias removed,
% and locates the start, end, and size of isolated steps within each sweep.
% Staircases and more complicated step patterns are not supported.
%
% Data for skipped sweeps will be filled with NaNs.
%
% EXAMPLE USAGE:
%   stepSize = findSteps(nSweeps, stimData, sf, threshold)
%   [stepSize, stepStarts, stepEnds] = ...
%         findSteps(nSweeps, stimData, sf, threshold, approvedTraces)
%   [stepSize, stepStarts, stepEnds] = ...
%         findSteps(nSweeps, stimData, sf, threshold, 'param1', 'value1')
%
% INPUTS:
%   nSweeps         scalar double      Total number of sweeps in the series.
%   stimData        double array       Stimulus data, in um, with
%                                      dimensions: [sweepLength nSweeps].
%                                      For absolute size, ensure that the
%                                      baseline is set at 0.
%   sf              scalar double      Sampling frequency, in kHz.
%   threshold       scalar double      Threshold for finding steps, in um.
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
%   'endTime'       double             Reverses step finding procedure by
%                                      first finding the end of the step
%                                      and then subtracting the length of 
%                                      time (in ms) set by this parameter
%                                      to find the start of the step.
%                                      Default is 0, for normal process.
%
% Created by Sammy Katta on 19th November, 2015.

function [stepSize, stepStarts, stepEnds] = findSteps(nSweeps, stimData, sf, threshold, varargin)

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

p.addParameter('nStim', 1, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('minStimInterval',300, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('roundedTo',0.1, @(x) isnumeric(x));
p.addParameter('endTime',0, @(x) isnumeric(x) && isscalar(x));
p.parse(nSweeps, stimData, sf, threshold, varargin{:});

approvedTraces = p.Results.approvedTraces;
nStim = p.Results.nStim;
minStimInterval = p.Results.minStimInterval;
roundedTo = p.Results.roundedTo;
endTime = p.Results.endTime;
%TODO: make endTime able to be a vector, for different step lengths in each
%sweep.

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
        
        if endTime == 0 % For normal step finding
            % Figure out timepoints when stimulus command for when step starts
            % and ends, and get the size in microns.
            stepZero = stimData(:,iSweep) - mean(stimData(1:10*sf,iSweep));
            stepOn = find(stepZero > threshold);
            stepLength = length(stepOn);
            if stepLength == 0
                continue
            end
            stepStarts(iSweep) = stepOn(1);
            stepEnds(iSweep) = stepOn(1) + stepLength;
            stepSize(iSweep) = mean(stepZero(stepStarts(iSweep)+1:stepEnds(iSweep)-1));
            
        else % Reverse step finding from end of step, given step length.
             % Useful for ramp on/step off stimuli.
             
            % Subtract baseline
            stepZero = stimData(:,iSweep) - mean(stimData(1:10*sf,iSweep));
            % Find where step signal is above threshold
            stepOn = find(stepZero > threshold); 
            % Calculate step length based on given time
            stepLength = endTime*sf;
            
            stepEnds(iSweep) = stepOn(end)+1; % first point below threshold
            stepStarts(iSweep) = stepOn(end) - stepLength;
            stepSize(iSweep) = mean(stepZero(stepStarts(iSweep)+1:stepEnds(iSweep)-1));
        end
                
    end
end

% By default, round step size to nearest 0.1 to drop noise and allow
% grouping of step sizes. roundedTo is an optional input that can be set
% larger or smaller depending on the range of step values used.
roundedTo = 1/roundedTo;
stepSize = round(stepSize*roundedTo)/roundedTo;

end