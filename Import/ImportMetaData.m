% ImportMetaData.m
%
% This function asks user for an Excel (xlsx) file with metadata, and
% imports it as a cell array.
% Numbers are imported as doubles, strings as char cells.
% 
% OUTPUTS:
%   rawMetaData        cell array   Cell array in which numbers are
%                                   imported as doubles, and strings/text
%                                   as char cells.
% 
% OPTIONAL PARAMETERS:
%   StripHeader        0/1          If you want to strip the top (header)
%                                   row, set the flag to 1. Default = 0,
%                                   which keeps all rows.
%
% Created by Sammy Katta on 4 December 2015.

function rawMetaData = ImportMetaData(varargin)

p = inputParser;

p.addParamValue('StripHeader', 0, @(x) isnumeric(x) && x==any([0 1]));
p.parse(varargin{:});

stripHeader = p.Results.StripHeader; 



[fname, pname] = uigetfile('*.xls;*.xlsx','Pick Excel file with metadata');
[~, ~, rawMetaData] = xlsread(fullfile(pname,fname));

if stripHeader  % if file has unwanted header row, get rid of it.
    rawMetaData = rawMetaData(2:end,:);
end

end