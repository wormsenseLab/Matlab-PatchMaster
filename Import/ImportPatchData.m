% ImportPatchData.m
% 
% ImportPatchData runs a modified HEKA import function (through sigTOOL) 
% that imports Patchmaster .dat files variables and puts them into a Matlab
% variable. The traces inside are broken down by group (a biological cell) 
% and series (a pgf from a protocol).
% 
% EXAMPLE:
%   ephysData = ImportPatchData()
%   ephysData = ImportPatchData(previousData)
%   [ephysData, tree] = ImportPatchdata()
% 
% OPTIONAL INPUTS: 
%   previousData    struct        If you have previously created a data
%                                 struct with this function, input it
%                                 to append more fields (groups) to it.
%
% OUTPUT:
%   ephysData       struct        Data is output as a nested struct. Each
%                                 group is a struct inside this. Each group
%                                 struct contains the filename/date, the
%                                 list of pgfs run, and the data in a cell
%                                 array with dimensions (series, channel).
%                                 This allows referencing by group name
%                                 using dynamic field notation.
%
%   tree            cell          Optional output, if you want to be able
%                                 to see and reference the original
%                                 metadata tree to see what other
%                                 information is available.
% 
% IMPORTANT NOTES:
% You must run sigTOOL at the beginning of every MATLAB session, though you
% may close the sigTOOL window once it opens.
% sigTOOL;
% 
% You must also replace SigTOOL's ImportHEKA function with ImportHEKAtoMat,
% which actually outputs two Matlab variables (containing the data and
% metadata for your file) into the workspace.
% Place ImportHEKAtoMat.m in this folder:
% 'sigTOOL\sigTOOL Neuroscience Toolkit\File\menu_Import\group_NeuroScience File Formats'
% 
% Created by Sammy Katta, 28-May-2014

function [ephysData, tree] = ImportPatchData(varargin)

% Pick Patchmaster files to import
[filename, pathname] = uigetfile('*.dat', 'Pick Patchmaster files', 'MultiSelect', 'on');

% If only one file is picked, make it a 1x1 cell instead of char
if ~iscell(filename)
    filename = {filename};
end

% If user doesn't give an existing data struct, create a new one. If user 
% does input a struct, use it and add to it. If user inputs something that
% is not a struct, throw an error. 
switch nargin
    case 0
        ephysData = struct();
    case 1
        ephysData = varargin{1};
        if ~isstruct(ephysData)
            error('PatchDataAnalysis:inputNotStruct', ...
                'Input must be a struct to which you wish to append new data fields.')
        end
    otherwise
        error('PatchDataAnalysis:tooManyArguments',...
            'Input must be a single struct to which you wish to append new data fields.');
end

% Load files
for iFile = 1:length(filename)
    fName = fullfile(pathname, filename{iFile});
    
    % Run modified HEKA import function
    [~, tree, data] = ImportHEKAtoMat(fName);
    % data output contains weirdly padded cells, the following will fix that
    
    % Collapse the padded data into a single row and pull apart
    % multi-channel sweeps later, depending on which protocols/series are
    % recording multiple channels.
    for i = length(data):-1:1
        dCollapse(1:length(data{i}))= data{i};
    end
    
    % Save the original tree and data for now, for checking accuracy.
    [~,saveName] = fileparts(filename{iFile});
%     eval([['tree' saveName] ' = tree;']);
%     eval([['data' saveName] ' = data;']);
    
    % Split the data into series by recording name, etc. and assign into
    % the final data structure
    ephysData = SplitSeries(tree, dCollapse, ephysData, saveName);
    
end

end