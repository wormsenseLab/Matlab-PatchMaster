% FilterRecordings.m
% 
% 
% Created by Sammy Katta on 4 Jan 2018.

function filteredCells = FilterRecordings(ephysData,ephysMetadata,varargin)

p = inputParser;
% metadata spreadsheet (optional or required? only necessary if using the filtering parameters)
p.addRequired('ephysData');
p.addRequired('ephysMetadata'); 

p.addOptional('allCells', cell(0));

p.addParameter('strain', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('externalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('internalSolution', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('cellType', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));
p.addParameter('wormPrep', cell(0), @(x) iscell(x) && ~isempty(x) && ischar(x{1}));

p.parse(ephysMetadata, varargin);

allCells = p.Results.allCells;

% doing it this way makes it easier to add new parameters straight into
% inputParser, rather than hardcoding parameter names.
filterParams = rmfield(p.Results,{'ephysData','ephysMetadata','allCells'});
filterHeaders = fieldnames(filterParams);

if isempty(allCells)
   allCells =  fieldnames(ephysData);
end

% find metadata rows for the cells specified by allCells and pull out the
% relevant parameters used for filtering.
allCellInd = cellfun(@(x) strcmp(ephysMetadata(:,1),x),allCells,'un',0)';
allCellInd = sum([allCellInd{:}],2);

% find columns in metadata spreadsheet whose headers match the filter
% parameter names (case- and whitespace-insensitive).
for i = 1:length(filterHeaders)
    paramInd(i) = find(strcmpi(regexprep(ephysMetadata(1,:), '\s+', ''),filterHeaders{i}));
end



end