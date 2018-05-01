% matchProts.m
%
% Find the locations/indices of series matching the desired protocol name.
% matchProts is case-insensitive.
%
% USAGE:
%   protLoc = matchProts(ephysData, cellName, protName)
%   protLoc = matchProts(ephysData, cellName, protName, 'MatchType', 'first')
%
% INPUTS:
%   ephysData       struct      Imported data from ImportPatchData.
%
%   cellName        char        Name of the recording to look at (must be a
%                               field within ephysData).
%
%   protName        char/cell   String or cell array of strings with the
%                               names of protocols to be matched.
%
% OUTPUTS:
%   protLoc         double      Vector of the indices of all series that
%                               match the given set of protocol names, in
%                               numerical order.
%
% OPTIONAL PARAMETERS:
%   'MatchType;     char        String specifying how to match the protocol
%                               name. Defaults to 'full'.
%                                   'full' finds series where the name
%                                       exactly matches protName
%                                   'first' finds series where the start of
%                                       the series name matches protName
%                                   'last' finds series where the end of
%                                       the series name matches protName
% 
%   'IgnoreSpecial' logical     Ignore non-alphanumeric characters when
%                               matching protocol names. If 0/false, names 
%                               must match all special characters. Default
%                               is 1/true.
%
% Created by Sammy Katta on 23 February 2016.


function protLoc = matchProts(ephysData, cellName, protName, varargin)

% Parse inputs and allow proper input of parameter/value pairs.
p = inputParser;
p.addRequired('ephysData');
p.addRequired('cellName');
p.addRequired('protName');

% consider whether you want to specify default filtering parameters or use zero-length cells
% for no filtering
p.addParameter('MatchType', 'full', @(x) validateattributes(x,{'char'},{'nonempty'}));
p.addParameter('ignoreSpecial',true); % if 1/true, ignores anything other than alphanumeric characters

p.parse(ephysData, cellName, protName, varargin{:});

matchType = p.Results.MatchType;
stripFlag = p.Results.ignoreSpecial;
if ~islogical(stripFlag)
    stripFlag = logical(stripFlag);
end

% If only one protocol name was given, make it a cell.
if ischar(protName) && ~iscell(protName)
    protName = {protName};
end

% Find matches for the protocol name in the list of series names for the
% given recording, based on the user's specification.
%NOTE: In Matlab R2016b+, see if this can be replaced by fxn contains().
protLoc = cell(length(protName),1);
prots = ephysData.(cellName).protocols;

if stripFlag
    protName = regexprep(protName,'[^a-zA-Z0-9]','');
    prots = regexprep(prots,'[^a-zA-Z0-9]','');
end

for i = 1:length(protName)
    switch matchType
        case 'full' % match the string exactly
            protLoc{i} = find(strcmpi(protName{i},prots));
        case 'first' % match protName to the beginning of the series name
            protLoc{i} = find(strncmpi(protName{i},prots,length(protName{i})));
        case 'last' % match protName to the end of the series name
            flippedProts = cellfun(@fliplr, prots,'UniformOutput', false);
            protLoc{i} = find(strncmpi(fliplr(protName{i}),flippedProts,length(protName{i})));
    end
    
end

% Combine the list of all found series indices into one double vector,
% sorted numerically. Indices that matched with multiple protocol names are
% only included once.
protLoc = unique(sort([protLoc{:}]));

end