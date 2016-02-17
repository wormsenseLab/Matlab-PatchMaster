% FilterProjectData.m
% 
% This function is used to filter out crap recordings or unidentified cells
% that will not be analyzed, and keep only data for the project(s) of
% interest. It will match data beginning with the given project prefix, 
% regardless of the rest of the name. It will also put the fields in
% alpha/numerical order after filtering.
% 
% USAGE:
%   ephysData = FilterProjectData(ephysData, projects)
% 
% INPUTS:
%   ephysData                       Struct array containing data imported
%                                   from Patchmaster by ImportPatchData.m.
% 
%   projects                        Cell array of strings containing
%                                   project name/prefix.
% 
% OUTPUTS:
%   ephysData                       Filtered struct array.

function ephysData = FilterProjectData(ephysData, projects)

dataFields = fieldnames(ephysData);
validProject = zeros(length(dataFields),1);

for iProj = 1:length(projects)
    isValid = strncmpi(dataFields,projects{iProj},length(projects{iProj}));
    validProject = validProject + isValid;
end

ephysData = rmfield(ephysData,dataFields(~logical(validProject)));
ephysData = orderfields(ephysData);

end