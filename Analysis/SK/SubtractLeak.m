% SubtractLeak.m
% This function takes the current/data trace from a given Patchmaster
% series and subtracts the leak during the baseline from the whole trace.
%
% EXAMPLE USAGE:
%   leakSubtract = SubtractLeak(data, sf)
%   leakSubtract = SubtractLeak(data, sf, approvedTraces)
%   leakSubtract = SubtractLeak(data, sf, 'baseLength', timeInMs)
%   [leakSubtract, leak] = SubtractLeak(data, sf, approvedTraces, ...
%         'baseLength', timeInMs)
%
% INPUTS:
%   data            double array       Traces on which to perform leak
%                                      subtraction with dimensions:
%                                      [sweepLength nSweeps].
%   sf              scalar double      Sampling frequency, in kHz.
%
% OUTPUT:
%   leakSubtract    double array       Leak-subtracted data, with
%                                      dimensions: [sweepLength nSweeps].
%                                      Skipped sweeps are filled with NaNs.
%
% OPTIONAL INPUTS:
%   approvedTraces  double vector      List of sweep numbers in the given
%                                      series to be included in the
%                                      analysis. Output for sweeps not on
%                                      this list will be filled with NaNs.
%                                      If not provided, all sweeps will be
%                                      included.
%
% OPTIONAL PARAMETER/VALUE PAIRS:
%   'BaseLength'    double             Length of time in ms at beginning of
%                                      sweep to use as baseline for
%                                      calculating leak for that sweep.
%                                      Default is 30ms.
%
% OPTIONAL OUTPUTS:
%   leak        double array           Magnitude of leak subtracted from
%                                      sweep. Skipped sweeps will be NaNs.
%
% Created by Sammy Katta on 19th November, 2015.

function [leakSubtract, leakSize] = SubtractLeak(data, sf, varargin)

% INPUTS
% Validate inputs, allow optional inputs for more complicated usage of the
% function, and set defaults for optional inputs not passed in.

p = inputParser;

p.addRequired('data', @(x) isnumeric(x));
p.addRequired('sf', @(x) isnumeric(x) && isscalar(x) && x>0);

% default: include all sweeps; check that sweep number is smaller than
% total sweep number
p.addOptional('approvedTraces', [], @(x) isnumeric(x) && max(x)<nSweeps);
% default: first 30ms used as baseline to find leak current
p.addParameter('BaseLength', 30, @(x) isnumeric(x) && isscalar(x) && x>0)

p.parse(data, sf, varargin{:});

nSweeps = size(data,2);

if p.Results.approvedTraces
    approvedTraces = p.Results.approvedTraces;
else
    approvedTraces = 1:nSweeps;
end
baseLength = p.Results.BaseLength * sf;

% ANALYSIS
% Subtract leak for approved traces if a list was given, or for all traces
% if no list was given. Skipped traces will be left as NaNs.
leakSubtract = NaN(length(data),nSweeps);
leakSize = NaN(nSweeps,1);

for iSweep = 1:nSweeps
    
    % Analyze only desired traces within this series
    if any(approvedTraces == iSweep)
        
        % Subtract leak/baseline
        leakSize(iSweep) = mean(data(1:baseLength,iSweep));
        leakSubtract(:,iSweep) = data(:,iSweep) - leakSize(iSweep);
        
    end
end


end