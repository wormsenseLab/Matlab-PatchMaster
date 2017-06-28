% ImportCalibImgTracks.m
% 
% This function reads in csv files containing tracking data from the
% trackCircle.py script created by Holger Fehlauer. That script finds the
% black bead in an image of a moving probe, fits a circle to it, and tracks
% the motion of the bead for matching to the photodiode data for
% calibration (at steady state).
%
% The data from the csv is saved in the field of the ephysData struct
% corresponding to the relevant recording.


function ephysData = ImportCalibImgTracks(ephysData, calibBase)

% keyboard;
% Pick Patchmaster files to import
[filename, pathname] = uigetfile('*.csv', 'Pick csv files output from trackCircle.py', 'MultiSelect', 'on');

% If only one file is picked, make it a 1x1 cell instead of char
if ~iscell(filename)
    filename = {filename};
end

for iFile = 1:length(filename)
   [~, recName] =  fileparts(filename{iFile});
   recName = strrep(recName,'result',''); 
   
   % Find name of existing recording matching tracked image and read in
   % image tracking data into that field of the ephysData struct.
   if ismember(lower(recName),lower(fieldnames(ephysData)))
       trackCol = find(~cellfun('isempty',(regexpi(calibBase(1,:),'Usable Image Tracking'))));
       try trackRow = find(~cellfun('isempty',(regexpi(calibBase(:,1),recName))));
       catch
           %If tracked image is not included in metadata sheet, let user know
           %and skip loading it.
           disp(sprintf('%s does not exist in metadata sheet.',recName));
           continue;
       end
       
       %If user has marked that the image tracking data is usable in the
       %metadata sheet, load and save the vector into the relevant recording
       %field.
       if strcmp(calibBase{trackRow,trackCol},'Y')
           ephysData.(calibBase{trackRow,1}).imgTrack = csvread(fullfile(pathname,filename{iFile}));
       else
           disp(sprintf('%s was excluded based on metadata',recName));
       end

   else

       % If tracked image doesn't match existing recording, let user know
       % and skip loading it.
       disp(sprintf('%s was not found in struct',recName));
       continue
   end
   
   
end

%read in ephysData, calibBase
%uiselect csv files to read in as separate variables
%match filename to calibData field names, case-insens, first or last
%if match exists and is marked as usable in calibBase
% then save var within that recording in subfield calibImgTrack
%if no match or not usable, disp message to that effect

end
