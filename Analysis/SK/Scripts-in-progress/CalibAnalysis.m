% Calibration Analysis Script
% For making loaded and unloaded calibration curves

%% Load data

calibData = ImportPatchData('incl',1);

% Keep only data with given project prefixes/names.
projects = {'PDC'};

%get rid of the crap/unrelated recordings
calibData = FilterProjectData(calibData, projects);
calibData = rmfield(calibData,{'PDC002_01_Testd','PDC002_05_Crapd'});

clear projects; fclose('all');
%% Get scale image and set bead size to find

% [scaleFile, scalePath] = uigetfile();
% scaleImage = imread(fullfile(scalePath,scaleFile));
% imtool(scaleImage); % use distance tool to select 50um distance
% 
% distConvert = 50/dist50um; % um per pixel
distConvert = 0.094;

beadSize = 22; % diameter in um
beadSizePx = beadSize/distConvert/2; % radius in pixels
beadSizeRange = [round(beadSizePx*0.9) round(beadSizePx*1.05)];

clear scaleFile scalePath

%%

pname = uigetdir();
calibFiles = dir(fullfile(pname,'*.tiff'));
calibFilesCell = {calibFiles(:).name}';
newStyle = 1; % 0 if stills from video, separated into a single set of steps per dat file
              % 1 if 45 stills for 3*15 steps

calibRecs = fieldnames(calibData);

allDists = [];
allXDists = [];

for iRec = 1:length(calibRecs)
   
    thisRec = calibRecs(iRec);  
    whichFiles = find(strncmpi(thisRec,calibFilesCell,9));
    
    if ~isempty(whichFiles)
        beadMotion = cell(length(whichFiles),1);
        
        for iImage = 1:length(whichFiles)
            thisFile = imread(fullfile(pname,calibFilesCell{whichFiles(iImage)}));
            [centers, radii] = imfindcircles(thisFile,beadSizeRange,'ObjectPolarity','dark','Sensitivity',0.9,'Method','twostage');
%             imshow(thisFile);
%             viscircles(centers,radii,'Enh',0,'LineSty','--','Color','k','LineW',1);
%             uiwait;
            centers = centers(1,:); radii = radii(1); % keep only the strongest circle
            beadMotion{iImage} = [centers radii];
            
        end
        
        beadMotion2 = [beadMotion(2:end); {nan(1,3)}];
        stimDist = cellfun(@(x,y) pdist([x(1:2);y(1:2)]),beadMotion,beadMotion2);
        
%         distTrace = cellfun(@(x,y) pdist([x(1:2);y(1:2)]),beadMotion,repmat(beadMotion(1),size(beadMotion)));
        stimLocs = find(stimDist>10);
        stimStart = [1;find(diff(stimLocs)>2)+1];
        
        if newStyle
            beadZeros = cellfun(@(x) repmat(x,10,1), beadMotion(stimLocs(stimStart)),'un',0);
            beadZeros = vertcat(beadZeros{:});
            beadZeros = mat2cell(beadZeros,ones(30,1),3);
        else
            beadZeros = repmat(beadMotion(1),length(stimLocs),1);
        end
        
        stimDist = cellfun(@(x,y) pdist([x(1:2);y(1:2)]),beadZeros,beadMotion(stimLocs));
%         stimDir = nanmean(cellfun(@(x,y) rad2deg(atan((y(2)-x(2))/(y(1)-x(1)))),beadZeros(2:end),beadMotion(stimLocs(2:end))));
% 
%         beadEndMotion = cell(0);
%         for iImage = 1:length(stimLocs)
%             thisFile = imread(fullfile(pname,calibFilesCell{whichFiles(stimLocs(iImage))}));
%             thisFileRot = ~im2bw(thisFile,0.2);
%             thisFileRot = imrotate(thisFileRot,stimDir);
%             stats = regionprops(thisFileRot,'Centroid','Area','Extrema');
%             [~, regionIdx] = sort(vertcat(stats(:).Area));
%             %largest region is white bg, second largest is black bead
%             beadEndMotion{iImage} = max(stats(regionIdx(end)).Extrema(end-1:end,:),[],1);   
%         end
%         
%         if newStyle
%             beadEndZeros = cellfun(@(x) repmat(x,10,1), beadEndMotion([1 11 21]),'un',0);
%             beadEndZeros = vertcat(beadEndZeros{:});
%             beadEndZeros = mat2cell(beadEndZeros,ones(30,1),2);
%         else
%             beadEndZeros = repmat(beadEndMotion(1),length(stimLocs),1);
%         end
%         stimDist = cellfun(@(x,y) pdist([x;y]),beadEndZeros,beadEndMotion');
        stimDist = reshape(stimDist,10,[]);
        stimDist = mean(stimDist,2);
        
        stimDist_um = stimDist*distConvert;
        stimDist_axial = stimDist_um/cosd(17);
        
        allDists = [allDists, stimDist_axial];
        allXDists = [allXDists, stimDist_um];
        
    end
end

%% Get PD measured stim size to compare with imaged displacement
% (without calibrating PD)

loadedMean = mean(allDists(:,1:5),2);
loadedFiles = fieldnames(calibData);
loadedFiles = loadedFiles(23:28);
sf = 5000; %Hz
chan = 3;
respTime = 0.2; % 200 ms
startTime = 0.05; %first 50ms

eachPos = [0 2.5 5 7.5 10 12.5 10 7.5 5 2.5 0];
stimStarts = (0:10)*sf;
stimEnds = (1:11).*sf;
stimMeans = [];
stimStartMeans = [];

for iRec = 1:length(loadedFiles)
   thisRec = [calibData.(loadedFiles{iRec}).data{chan, 1:3}];
   meanPD = mean(thisRec,2);
   meanPD = abs(meanPD-mean(meanPD(1:respTime*sf)));
   
   stimMeans(:,iRec) = arrayfun(@(x) mean(meanPD(x-respTime*sf:x)), stimEnds)';
   stimStartMeans(:,iRec) = arrayfun(@(x) mean(meanPD(x:x+startTime*sf)),stimStarts+1)';
end

stimMeansNorm = bsxfun(@rdivide,stimMeans,max(stimMeans));

%% Calculate PD vs. motion curve

eachPosPD = eachPos*cosd(18); %the angle is 17, but PD positions were calculated using 18 degrees
stimStarts = (0:10)*sf+1;
stimStarts(7:end) = stimStarts(7:end)+(0.3*sf); %middle step is 200ms longer, correct for that
respTime = 0.15; %150 ms, shorter bc PD takes time to move
startTime = 0.100; %first 50ms

for iRec = 1:length(loadedFiles)
    thisRec = [calibData.(loadedFiles{iRec}).data{chan,4:5}];
    wormSub = thisRec(:,1)-thisRec(:,2);
    
    zeroPD = bsxfun(@(x,y) x-y,thisRec, mean(thisRec(1:startTime*sf,:),1));
    zeroWormSub = zeroPD(:,1)-zeroPD(:,2);
    
    pdCalMeans(:,iRec) = arrayfun(@(x) mean(zeroPD(x:x+respTime*sf,1)), stimStarts)';
    pdCalSubMeans(:,iRec) = arrayfun(@(x) mean(wormSub(x:x+respTime*sf)), stimStarts)';
    pdCalZeroSubMeans(:,iRec) = arrayfun(@(x) mean(zeroWormSub(x:x+respTime*sf)), stimStarts)';
end


%% plot stuff
plot(eachPos(1:end-1),allDists(:,1:5),'r');
hold on;
plot(eachPos(1:end-1),allDists(:,6),'b');
xlabel('Commanded Position');
ylabel('Imaged Position');
plotfixer;
fplot(@(x) x,[0 15],'k:');

figure();
plot(distTrace*distConvert/cosd(17));
xlabel('Time (au)'); ylabel('Imaged Position (um)');

figure();
plot(allDists(:,1:5),stimMeansNorm(1:end-1,1:5));
xlabel('Imaged Displacement (um)');
ylabel('Normalized Photodiode Signal (V)')
hold on;
% plot(allDists(:,6),stimMeansNorm(1:end-1,6),'b');
plotfixer;
fplot(@(x) x/12.5,[0 14],'k:')