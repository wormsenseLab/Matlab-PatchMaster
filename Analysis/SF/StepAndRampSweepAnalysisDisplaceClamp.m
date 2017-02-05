%%% Sylvia Fechner
%%% Stanford University, 
%%% 20151115
%%% Update: 20160210
%%% Update: 20160901
%%% Update: 20170204
%%% Script to analyze data from FALCON in Displacement Clamp
%%% To commit to github
%%% go to my branch
%%% git add .
%%% git commit -m 'text'
%%% git push origin my branch


%%% ToDo 

%%% Currently, automatically saved data are saved in the same folder
%%% change automatically for Step and Ramp differently
%%% write into csv file: mac can't write to excel !! 
    % --> which values?
    % how to write each Indentation as col header??? in excel and csv
    % mean indentation from non-averaged signal?
    % normalized max current non-averaged?
    % write which simuli which series and which sweeps kept
    % n of averaged signal
%%% finding function in igor to load excel or csv files
%%% make it work for ramps as well
%%% running average
%%% ToDo maybe: plot (one) indentation over time; currently: Current over
%%% Filenumber

%%% ToDo maybe: make subplots for all Fivestep blocks?
%%% find number of nan values to find out the number of averages
%%% include legend again in current vs indentation
%%% TODO: change, when do take the baseline. current shifts now from -60 to
%%% -85 

%%  load dat.files 
clear all; close all; clc;
ephysData=ImportPatchData();
%load('ephysdata(20170130).mat')
%%
% load notes to get several values automatically needed for the conversion of the signals
loadFileMode = 1; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('Ephys-Meta-Sylvia.xlsx'); % be careful in which Folder saved.
end


%% Analysis Individual Recording 
close all; clc

%%% hardcoding part:
name = 'STF070'; % name of recording. placed into varaibel fiels names%
stimuli = 'FiveStep'; 
OnlyMechano = 0; % if = 0, then FALCON, if 1, then ForceClamp Only 
ReadFromSheet = 0; % if = 0, then command promp to delete block, if 1, then read from MetaDataSheet 
% protocol names:
% Single protocols: Step and Ramp-Hold; 
% Five sweeps per protocol:FiveStep, FiveRampHold; does not work with alternative names
%%%%%

Filenumber = 1; % wil be used to extract sampling freuqnecy; first file loaded, maybe change (ToDO: check if I did it)

    
Files = 1:length(ephysData.(name).protocols);% load all protocols  

% load all data from all protocols 
% load Current, Actuator Sensor, And Cantilver Signal of each step%
% if else statement included to also load protocols which don't have
% ForceClamp data; not necessary for Current
A=[];B=[]; C=[];D=[];   % to get an empty array, if I analyzed different files before
for i = Files(:,1):Files(:,end);
A{i} = ephysData.(name).data{1, i}; %Current
if isempty(ephysData.(name).data{3, i}) == 1
    continue
 else
B{i} = ephysData.(name).data{3, i}; % actuator Sensor
end
 if isempty(ephysData.(name).data{4, i}) ==1
    continue
 else
 C{i} = ephysData.(name).data{4, i}; % Cantilever Signal
 end
  if isempty(ephysData.(name).data{2, i}) == 1 % actuator input
     continue
 else
 D{i} = ephysData.(name).data{2, i}; % actuator SetPoint
  end
end


% find all files with certain protocol name: 
% ifelse statement: if not respective stimuli name (FiveStep or FiveRampHold), then empty array; 
for i = Files(:,1):Files(:,end);
   if find(strcmpi(ephysData.(name).protocols{1,i}, stimuli)) == 1;
        continue
   else 
         A{1, i} = []; B{1, i} = []; C{1, i} = []; D{1, i} = [];          
    end      
end


%compare Input Stimuli
isFiveStep = strcmp('FiveStep',stimuli); 
isFiveRamp = strcmp('FiveRampHold',stimuli);
isStep = strcmp('Step',stimuli);
isRamp = strcmp('Ramp-Hold',stimuli);
isFiveSine = strcmp('FiveSinus',stimuli);
isFifteenStep = strcmp('FifteenStep',stimuli); 
isIVStep = strcmp('IVStep',stimuli); 


% showing in command prompt: AllStimuli = patchmaster Filenumber;
% this helps to identify which "five Block to delete"

AllStimuliBlocks = (find(strcmpi(ephysData.(name).protocols, stimuli)))

if ReadFromSheet == 1;

    % finds meta datasheet, if deleting rec is preassigned
headers = raw(1,:);
FindRowRecording = strcmpi(raw,name); 
[FindRowRecording,col] = find(FindRowRecording,1); % 1 indicates to use only the first row with this name; in metadata sheet replicates of rec name
CellFiveBlock = find(strcmpi(headers, 'FiveBlockRec'));
AllFiveBlockUsed = raw(FindRowRecording,CellFiveBlock);
AllFiveBlockdeleted = [];
CellFiveBlockdeleted = find(strcmpi(headers, 'deletedFiveBlock'));
AllFiveBlockdeleted = raw(FindRowRecording,CellFiveBlockdeleted);
AllFiveBlockdeleted = [AllFiveBlockdeleted{:}] ;


SizeOfDeletedBlocks = size(AllFiveBlockdeleted,2);
if SizeOfDeletedBlocks == 1
display 'only one block deleted'
elseif SizeOfDeletedBlocks > 1
    AllFiveBlockdeleted = str2num(AllFiveBlockdeleted); % creates double, needed for for loop  
end
%end

%AllFiveBlockdeleted = cell2mat(AllFiveBlockdeleted);
% deleting whole blocks of FiveBlockStimuli; Whole block=Filenumber
if isnan(AllFiveBlockdeleted)==1
 display 'NaN value, no block deleted'
else
for i= 1:length(AllFiveBlockdeleted)
    A{1, AllFiveBlockdeleted(i)}  = []; 
    B{1, AllFiveBlockdeleted(i)}  = [];
    C{1, AllFiveBlockdeleted(i)}  = [];
    D{1, AllFiveBlockdeleted(i)}  = []; 
end
end

else
    
while 1
prompt = {'BlockNr (leave empty to quit)'};
dlg_title = 'Delete a block?';
num_lines = 1;
defaultans = {''};
IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);
FirstValue = str2num(IndValues{1});

if isempty(FirstValue) == 1 
     break
 else
    A{1, FirstValue}  = []; 
    B{1, FirstValue}  = [];
    C{1, FirstValue}  = [];
    D{1, FirstValue}  = []; 
 end
end

end

% removes all empty cells from the cell array
AShort = A(~cellfun('isempty',A)); BShort = B(~cellfun('isempty',B)); CShort = C(~cellfun('isempty',C)); DShort = D(~cellfun('isempty',D));

% concatenating all stimuli
Aall = []; Aall = cat(2,AShort{:}); 
Ball = []; Ball = cat(2,BShort{:}); 
Call = []; Call = cat(2,CShort{:});
Dall = []; Dall = cat(2,DShort{:});

% calculate sampling frequency
fs = ephysData.(name).samplingFreq{1, Files(:,AllStimuliBlocks(1))}; % sampling frequency from first Stimuli loaded; 
interval = 1/fs;   %%% get time interval to display x axis in seconds 
ENDTime = length(Aall)/fs; %%% don't remember why I complicated it
Time = (0:interval:ENDTime-interval)'; 



ActuSensor = [];
for i = 1:size(Ball,2),
ActuSensor(:,i) = Ball(:,i)*1.5; % 1.5 = sensitivity of P-841.10 from Physik Instrumente; travel distance 15 um; within 10 V; ToDo: measure real sensitivity
end



%%%%%% ForceClampSignals %%%%%%%

% to get Deflection of Cantilever: multiply with Sensitivity 
% get Sensitivity from Notes Day of Recording  
% FindRowStiff = strcmpi(raw,name); % name = recorded cell
% [Stiffrow,col] = find(FindRowStiff); % Siffrow: row correasponding to recorded cell
% 
% headers = raw(1,:);
% ind = find(strcmpi(headers, 'Sensitivity(um/V)')); % find col with Sensitivity
% Sensitivity = raw(Stiffrow,ind); 
% Sensitivity = cell2mat(Sensitivity);
% 
% indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
% Stiffness = raw(Stiffrow,indStiffness); 
% Stiffness = cell2mat(Stiffness);
%  
% 
% [Start,ActuSetPointZero,CantiDefl,Indentation,MeanIndentation,Force,MeanForce,normCantiDefl,allRiseTime,allOvershoot] = AnalyzeForceClampBOM(interval,ActuSensor,isFiveStep,isFifteenStep,Ball,Call,Dall,Sensitivity,Stiffness,fs);
% 
% %%%calculate Stiffness
% 
% MeanIndentationVer = MeanIndentation'
% MeanForceVer = MeanForce'
% 
% StiffnessWorm = MeanIndentationVer\MeanForceVer % is doing a leastsquarefit



[SlopeActu,MaxZeroSlopeActu,StdZeroSlopeActu,MaxZeroActu,StdZeroActu,MaxActuSensorPlateau,StdActuSensorPlateau,CellMaxActuFirst] = SlopeThreshold(ActuSensor);  


%calculate threshold (needed to determine the Onset (Start) of the Stimulus)
StartBase = [];
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1 || isIVStep == 1;
   StartBase = MaxZeroActu + 4*StdZeroActu;   %% play around and modifz 
else
   StartBase = MaxZeroSlopeActu + 4*StdZeroSlopeActu; %
end


% play around and modify
disp 'change threshold, if noise of' 
BelowPlateau = [];
BelowPlateau = MaxActuSensorPlateau - 10*StdActuSensorPlateau; 


%%%%%% CurrentSignals %%%%%%%
%[AvgMaxCurrent,AvgMaxCurrentMinus,AvgMaxCurrentOff,AvgMaxCurrentMinusOff,Start,StartOffBelow,Ende,EndeOff,ASubtractAvg,LengthRamp,LengthInms,EndeRamp] = AnalyzeCurrent(isFiveStep,isStep,isFifteenStep,ActuSensor,StartBase,Aall,ASubtract,fs,SlopeActu,BelowPlateau,CellMaxActuFirst,interval);
[LeakA, ASubtract, AvgMaxCurrent,AvgMaxCurrentMinus,AvgMaxCurrentOff,AvgMaxCurrentMinusOff,Start,StartOffBelow,Ende,EndeOff,ASubtractAvg,LengthRamp,LengthInms,EndeRamp,StartOffBelowShort,ActuSensorAvg] = AnalyzeCurrent(isFiveStep,isStep,isFifteenStep,ActuSensor,StartBase,Aall,fs,SlopeActu,BelowPlateau,CellMaxActuFirst,interval);


% modify for RampAndHold
LengthInms
LengthRamp

% calculate in pA
AverageMaxCurrentMinusppA = AvgMaxCurrentMinus*10^12;  ASubtractppA = ASubtract*10^12; % current in pA to visualize easier in subplot
AallppA=Aall*10^12; 


%%%%%% ForceClampSignals %%%%%%%

% to get Deflection of Cantilever: multiply with Sensitivity 
% get Sensitivity from Notes Day of Recording  
FindRowStiff = strcmpi(raw,name); % name = recorded cell
[Stiffrow,col] = find(FindRowStiff,1); % Siffrow: row correasponding to recorded cell

headers = raw(1,:);
ind = find(strcmpi(headers, 'Sensitivity(um/V)')); % find col with Sensitivity
Sensitivity = raw(Stiffrow,ind); 
Sensitivity = cell2mat(Sensitivity);

indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
Stiffness = raw(Stiffrow,indStiffness); 
Stiffness = cell2mat(Stiffness);


[ActuSetPoint,CantiDefl,MeanIndentation,Force,MeanForce,Indentation,normCantiDefl,allRiseTime,allOvershoot] = AnalyzeForceClamp(interval,Start,isFiveStep,isStep,isFifteenStep,Dall,Call,EndeOff,ActuSensor,Sensitivity,Stiffness,fs);


% Calculating Rise time and Overshoot on Cantilever Deflection signals
% shortened to the Onset of the step

% ToDo: needs to be modified for Ramp
% ToDo: include in ForceClamp function
%  CantiDeflShort = [];
%  MeanCantiDefl = [];
%  normCantiDefl = [];
%   for i = 1:size(CantiDefl,2);
%       EndeCanti(i) = Start(i)+1000;
%   CantiDeflShort(:,i) = CantiDefl(Start(i):EndeCanti(i),i); 
%   MeanCantiDefl(i) =  mean(CantiDefl(0.2*fs:0.4*fs,i));
%   normCantiDefl(:,i) = CantiDefl(:,i)/MeanCantiDefl(i);
%   end
 
% TimeShort = (0:interval:length(CantiDeflShort)/fs-interval)';  
% InfoSignal = stepinfo(CantiDeflShort, TimeShort, MeanCantiDefl, 'RiseTimeLimits', [0.0 0.63]); %%% over sorted data?? 
% allRiseTime = cat(1,InfoSignal.RiseTime);
% allOvershoot = cat(1,InfoSignal.Overshoot);


%% now figures
close all
xScatter = (1:length(MeanIndentation));
LengthInmsForPlot = LengthInms*1000;

%%%current with and without leak subtraction in a subplot %%%%
if OnlyMechano  == 0;
figure()
for i = 1:size(AallppA,2)
subplot(ceil(size(AallppA,2)/5),5,i)
plot(Time,AallppA(:,i))
%ylim([-5*10^-11 1*10^-11])
hold on
plot(Time,ASubtractppA(:,i))
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1;
title(round(MeanIndentation(i),1)) %% 
else
    title(round(LengthInmsForPlot(i),1)) 
end
end
suptitle({'Current (pA) with (red) and without (blue) leak subtraction';'Bold numbers: Indentation in µm'}) %('')
elseif OnlyMechano  == 1;
    'hello mechano'
end
%%%cantilever signals %%%%


figure()
for i = 1:size(Aall,2)
subplot(ceil(size(AallppA,2)/5),5,i)
plot(Time,CantiDefl(:,i))
%ylim([-5*10^-11 1*10^-11])
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1;
title(round(MeanIndentation(i),1)) %% 
else
    title(round(LengthInmsForPlot(i),1)) 
end
end
suptitle({'Cantilever Deflection'}) %('')





%%%normalized cantilever signals %%%%
figure()
for i = 1:size(Aall,2)
subplot(ceil(size(AallppA,2)/5),5,i)
plot(Time,ActuSensor(:,i))
ylim([-1 16])
%ylim([-5*10^-11 1*10^-11])
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
if isFiveStep == 1 || isStep == 1|| isFifteenStep == 1;
title(round(MeanIndentation(i),1)) %% 
else
    title(round(LengthInmsForPlot(i),1)) 
end
end
suptitle({'Actuator Sensor'}) %('')

%control plot

figure()
subplot(3,2,1)
scatter(xScatter, Start) 
%ylim([0 1000]) % ToDo: change it for Ramps
ylabel('Point')
xlabel('Filenumber')
title('control: Find OnSet of On-Stimulus')
hold on 
subplot(3,2,2)
scatter(xScatter, StartOffBelow) %change to start of below
%ylim([100 1000]) % ToDo: change it for Ramps
ylabel('Point')
xlabel('Filenumber')
title('control: Find OnSet of Off-Stimulus')
hold on 
subplot(3,2,3)
scatter(xScatter, LengthInmsForPlot) %change to start of below
%ylim([100 1000]) % ToDo: change it for Ramps
ylabel('Lengtg of OnSet Stimulus (ms)')
xlabel('Filenumber')
title('Length of Ramp')
subplot(3,2,4)
scatter(xScatter, AverageMaxCurrentMinusppA ) %change to start of below
%ylim([100 1000]) % ToDo: change it for Ramps
ylabel('Current (pA)')
xlabel('Filenumber')
title('control: compare Max current of each stimulus with traces')


%%
%%% Calculate Velocity For RampAndHoldStimuli
%%% TODo: Calculate Velocity from Indentation? Yes. Do it. Do it from fit

 if isRamp == 1 || isFiveRamp == 1;
 for i = 1:size(Indentation,2);
cftool(Time(Start(i):Start(i)+LengthRamp(i)),Indentation(Start(i):Start(i)+LengthRamp(i),i))
 end 
 else
    disp 'not opening fitTool for Steps currently'
 end
 
%% Refit a Recording
 if isRamp == 1 || isFiveRamp == 1;
  while 1
prompt = {'Enter number of recording, matches subplot (leave empty to quit, enter a number as long as you want to refit a ramp)','Start of Fit: add or sub some points','End of Fit: add or sub some points'};%,'SecondRec','ThirdRec','ForthRec'};
dlg_title = 'Refit a recording?';
num_lines = 1;
defaultans = {'','',''};%,'',''};
IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);
FirstRec = str2num(IndValues{1});
SecondRec = str2num(IndValues{2});
ThirdRec = str2num(IndValues{3});
%ForthRec = str2num(IndValues{3});

if isempty(FirstRec) == 1
    break
else
cftool(Time(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec),Indentation(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec,FirstRec));
end
  end
 else
    disp 'not opening fitTool for Steps currently'
 end
    
%% fit to calculate Stiffness

FittingModeOn = 0; 

 if FittingModeOn == 1;
     
 for i = 1:size(Indentation,2);
%cftool(Time(Start(i):Start(i)+LengthRamp(i)),Indentation(Start(i):Start(i)+LengthRamp(i),i))
cftool(Indentation(Start(i):Start(i)+LengthRamp(i),i),Force(Start(i):Start(i)+LengthRamp(i),i))
 end 
 else
    continue
 end


%% Refit single Recording Force vs Indentation

 if FittingModeOn == 1;
  while 1
prompt = {'Enter number of recording, matches subplot (leave empty to quit, enter a number as long as you want to refit a ramp)','Start of Fit: add or sub some points','End of Fit: add or sub some points'};%,'SecondRec','ThirdRec','ForthRec'};
dlg_title = 'Refit a recording?';
num_lines = 1;
defaultans = {'','',''};%,'',''};
IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);
FirstRec = str2num(IndValues{1});
SecondRec = str2num(IndValues{2});
ThirdRec = str2num(IndValues{3});
%ForthRec = str2num(IndValues{3});

if isempty(FirstRec) == 1
    break
else
 cftool(Indentation(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec,FirstRec),Force(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec,FirstRec));   
%cftool(Indentation(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec,i),Force(Start(FirstRec)+SecondRec:Start(FirstRec)+LengthRamp(FirstRec)+ThirdRec,FirstRec),i);
end
  end
 else
    disp 'not opening fitTool for Steps currently'
 end

%% delete single recordings 
%close all
%ToDo: has to be changed for ramps, because I want to average the current with
%same velocity 

ASubtractNew = ASubtractAvg;
AvgMaxCurrentMinusNew = AvgMaxCurrentMinus;
MeanIndentationNew = MeanIndentation;
LeakANew = LeakA;
%AverageMaxNormCurrentNew = 
%TO DO: Someting wrong with the order in command promt
% if I redo AverageMaxCurrentMinus= Nan, I have to reload it again or do it
% as for ASubtract new

while 1
prompt = {'Enter number of recording, matches subplot (leave empty to quit, enter a number as long as you want to delete a recording)'};%,'SecondRec','ThirdRec','ForthRec'};
dlg_title = 'Delete a recording?';
num_lines = 1;
defaultans = {''};%,'','',''};
IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);
FirstRec = str2num(IndValues{1});
%SecondRec = str2num(IndValues{2});
%ThirdRec = str2num(IndValues{3});
%ForthRec = str2num(IndValues{3});

if isempty(FirstRec) == 1
    break
else
  ASubtractNew(:,FirstRec) = NaN;
  AvgMaxCurrentMinusNew(:,FirstRec) = NaN;
  MeanIndentationNew(:,FirstRec) = NaN;
  LeakANew(:,FirstRec) = NaN;
  
end
end

LeakAppA = LeakANew*10^12;
%% Sort Data 
 
% change that mergeInd can be finally made by rounded MeanIndentation where
% files were deleted --> problems with NaN values --> maybe solution see
% end of the script

%%%%%% Sort Data

MeanSameIndCurrent=[]; MeanTraces=[]; SortCurrent=[]; SortCurrentOff =[];
NumberTracesPerInd=[];transAsub =[];SortASubtract = [];
MeanSameVelIndentation=[];
SortIndentation=[];
transIndentation =[];
MeanTracesIndentation =[];
% SortData & analyze to actual Force value, previously all Force values
% calculated from correspondent Indentation; SortForce at the end of the name means sorted by the
% actual Force value!!!

SortForce =[]; MeanSameIndForce=[]; SortForceTraces=[];transForce =[];MeanTracesForce =[];
  

if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1;
    disp 'Step Round and sort'
RoundMeanInd = round(MeanIndentation,1); % change to get it from MeanIndentation New with NaN values
[SortInd sorted_index] = sort(RoundMeanInd'); % get index of mean indentations
SortCurrent = AvgMaxCurrentMinusNew(sorted_index);
SortCurrentOff = AvgMaxCurrentMinusOff(sorted_index); %change that it works for deleting traces
%SortNormCurrent = AverageMaxCurrentMinusNew(sorted_index); % don't need it
%here, because normalize it afterwards ? maybe change it?
SortForce = MeanForce(sorted_index);
transAsub = ASubtractNew';
SortASubtract = transAsub(sorted_index,:);
transIndentation = Indentation'; 
SortIndentation = transIndentation(sorted_index,:);
transForce = Force'; 
SortForceTraces = transForce(sorted_index,:);
transCellMaxActuFirst = [];
transCellMaxActuFirst = CellMaxActuFirst';
MeanCellMaxActuFirst = [];
SortCellMaxActuFirst = [];
SortCellMaxActuFirst = transCellMaxActuFirst(sorted_index,:);

%signals sorted by Force
% application only on mean Averaged currents so far
SortedForceForForce = []; RoundMeanForce = [];tansIndSortForce=[];transMeanForceForce=[];
RoundMeanForce = round(MeanForce,1);
[SortedForceForForce sorted_indexForce] = sort(RoundMeanForce');
SortASubtractSortForce = []; SortIndentationSortForce =[];SortForceTracesSortForce=[];
SortASubtractSortForce = transAsub(sorted_indexForce,:);
tansIndSortForce = MeanIndentation';
 SortIndentationSortForce = tansIndSortForce(sorted_indexForce,:);%ToDo: update in delete series
transMeanForceForce = MeanForce';
SortForceTracesSortForce= transMeanForceForce(sorted_indexForce,:);

MergeForce= [];
MergeForce = builtin('_mergesimpts',SortedForceForForce,0.2,'average'); 
tolerance = 0.2;
L =[];
[~,FRowForce] = mode(SortedForceForForce); %gets the Frequency of the most frequent value
FindSameIndForce= NaN(FRowForce,length(MergeForce)); 
if size(Aall,2) > 5
FindSameIndInitialForce = {};
for L = 1:length(MergeForce);
FindSameIndInitialForce{L} = find([SortedForceForForce] >MergeForce(L)-tolerance & [SortedForceForForce]<MergeForce(L)+tolerance);
end
FindSameIndNaNForce = padcat(FindSameIndInitialForce{:});
FindSameIndForce = FindSameIndNaNForce;

for i = 1:length(MergeForce);
[r,c] = find(isnan(FindSameIndForce(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameIndForce(:,i)))>0
FindSameIndForce(r,i) =FindSameIndForce(r-1,i);
end
end

else %% 
    FindSameIndForce = [];
    for L = 1:length(MergeForce);
FindSameIndForce(:,L) = find([SortedForceForForce] >MergeForce(L)-tolerance & [SortedForceForForce]<MergeForce(L)+tolerance);
    end 
end

MeanTracesSortForce= [];MeanTracesIndentationSortForce=[];MeanTracesForceSortForce=[];
for k = 1:length(MergeForce);
MeanTracesSortForce(k,:) = nanmean(SortASubtractSortForce((FindSameIndForce(1,k)):(FindSameIndForce(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
MeanTracesIndentationSortForce(k,:) = nanmean(SortIndentationSortForce((FindSameIndForce(1,k)):(FindSameIndForce(end,k)),:),1);
MeanTracesForceSortForce(k,:) = nanmean(SortForceTracesSortForce((FindSameIndForce(1,k)):(FindSameIndForce(end,k)),:),1);
end

MeanTracesSortForce = MeanTracesSortForce';
MeanTracesIndentationSortForce = MeanTracesIndentationSortForce';
MeanTracesForceSortForce = MeanTracesForceSortForce';

%%%% calculate merged Indentations and sort signals for it %%%%
%SortASubtract = SortASubtract'; keep it as row, easier to calculate
MergeInd = [];
MergeInd = builtin('_mergesimpts',SortInd,0.2,'average'); %%% merge values with +/- 0.1 distance
tolerance = 0.2; % tolerance to group indentations
k =[];
[~,FRow] = mode(SortInd); %gets the Frequency of the most frequent value
FindSameInd= NaN(FRow,length(MergeInd)); % determine row length with highest probability of most frequent value
if size(Aall,2) > 5
FindSameIndInitial = {};
for k = 1:length(MergeInd);
FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
end

FindSameIndNaN = padcat(FindSameIndInitial{:});
FindSameInd = FindSameIndNaN;

for i = 1:length(MergeInd);
[r,c] = find(isnan(FindSameInd(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameInd(:,i)))>0
FindSameInd(r,i) =FindSameInd(r-1,i);
end
end

ind = find(isnan(FindSameIndNaN));
FindSameIndNaN(ind)=0;
LogicOfIndentations =  FindSameIndNaN > 0;
NumberOfAvergagesPerInd = sum(LogicOfIndentations);
NumberOfAvergagesPerInd = NumberOfAvergagesPerInd';

%%%%%%% if only one FiveStepProtcol was applied %%%%
else %% 
    FindSameInd = [];
    for k = 1:length(MergeInd);
FindSameInd(:,k) = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
    end    
end   
%FindLogicalNumberOfTraces = FindLogicalNumberOfTraces'
%NumberTracesPerInd = NumberTracesPerInd';
for k = 1:length(MergeInd);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
MeanSameIndCurrent(k) = nanmean(SortCurrent(FindSameInd(:,k))); %average MaxCurrent*-1 with same indentation
MeanSameIndCurrentOff(k) = nanmean(SortCurrentOff(FindSameInd(:,k))); %TODO: something wrong
MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
MeanTracesIndentation(k,:) = nanmean(SortIndentation((FindSameInd(1,k)):(FindSameInd(end,k)),:),1);
MeanTracesForce(k,:) = nanmean(SortForceTraces((FindSameInd(1,k)):(FindSameInd(end,k)),:),1);
MeanCellMaxActuFirst(k) = nanmean(SortCellMaxActuFirst(FindSameInd(:,k)));
end


else
    disp 'Ramp Round and Sort';
    
 Velocity=Velocity';
 StiffnessRamp = StiffnessRamp';
   % RoundMeanInd = round(Velocity,1);   
[SortVel sorted_index] = sort(Velocity); % get index of mean Velocity
SortCurrent = AvgMaxCurrentMinusNew(sorted_index);
SortIndentation = MeanIndentation(sorted_index);
SortStiffnessRamp = StiffnessRamp(sorted_index);
%SortNormCurrent = AverageMaxCurrentMinusNew(sorted_index); % Do I need this; yes, normalized to Off Response
SortForce = MeanForce(sorted_index);
transAsub = ASubtractNew';
SortASubtract = transAsub(sorted_index,:);
%SortASubtract = SortASubtract'; keep it as row, easier to calculate
MergeVel = [];
SortVel = SortVel';
MergeVel = builtin('_mergesimpts',SortVel,10,'average');%'average'); %%% TODO: does not work merged Velocity values with +/- 0.1 distance
tolerance = 10; % tolerance to group velocity
%TODo: find best value for velocity merge
k =[];
[~,FRow] = mode(SortVel); % TODO: does not work for Velocity gets the Frequency of the most frequent value
FindSameInd= NaN(10,length(MergeVel)); % determine row length with highest probability of most frequent value
%ToDo: changeValue: currently hardcoded.
if size(Aall,2) ~=5% > 5 %ToDo: change, it is not working, if single ramp is equal 5
FindSameIndInitial = {};
for k = 1:length(MergeVel);
FindSameIndInitial{k} = find([SortVel] >MergeVel(k)-tolerance & [SortVel]<MergeVel(k)+tolerance);
end

FindSameIndNaN = padcat(FindSameIndInitial{:});
FindSameInd = FindSameIndNaN;

for i = 1:length(MergeVel);
[r,c] = find(isnan(FindSameInd(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameInd(:,i)))>0
FindSameInd(r,i) =FindSameInd(r-1,i);
end
end


else %% if only one FiveStepProtcol was applied
    FindSameInd = [];
    for k = 1:length(MergeVel);
FindSameInd(:,k) = find([SortVel] >MergeVel(k)-tolerance & [SortVel]<MergeVel(k)+tolerance);
    end   
end

%FindLogicalNumberOfTraces = FindLogicalNumberOfTraces'
%NumberTracesPerInd = NumberTracesPerInd'
for k = 1:length(MergeVel);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
MeanSameIndCurrent(k) = nanmean(SortCurrent(FindSameInd(:,k))); %average MaxCurrent*-1 with same indentation
MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
MeanSameVelIndentation(k) = nanmean(SortIndentation(FindSameInd(:,k)));
MeanSameVelStiffnessRamp(k) = nanmean(SortStiffnessRamp(FindSameInd(:,k)));
end
end

%%%% maybe useful for later
% problem with not existing NaN values in FindSameInd
% FindSameIndValues = NaN(5,17); %size FindSameInd
% FindSameIndCurrent = NaN(5,17);
%     for k = 1:numel(FindSameInd);
% FindSameIndValues(k)= SortInd(FindSameInd(k)); %) >MergeInd(k)-tolerance & SortInd<MergeInd(k)+tolerance;
% FindSameIndCurrent(k) = AvgMaxCurrentMinusNew(FindSameInd(k));
%     end   
%     
%TODO: MeanTraces: first 30 values NAN; why
%ToDo: MeanSameIndCurrent for off current

MeanTraces = MeanTraces'; % transpose to column vector for export to igor
MeanSameIndCurrent = MeanSameIndCurrent';
MeanSameIndForce = MeanSameIndForce';
MeanTracesppA = MeanTraces*10^12; % to get current in pA
MeanTracesForce = MeanTracesForce';
MeanTracesIndentation = MeanTracesIndentation';




if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1;
    disp 'step has no velocity calulation yet - maybe make it'
else
MeanSameVelIndentation = MeanSameVelIndentation';
VelocityVer = Velocity';
%MeanSameVelStiffnessRamp = MeanSameVelStiffnessRamp'
MeanSameVelStiffnessRampVer = MeanSameVelStiffnessRamp'
StiffnessRampVer = StiffnessRamp';
%MergeVel!!!
end

% include here calculation of 
% FindLogicalNumberOfTraces = isnan(FindSameInd) == 0;
% TracesPerIndentation = sum(FindLogicalNumberOfTraces);
% TracesPerIndentation = TracesPerIndentation';
%missing calculation for numer of traces for only one block


%calculate max current of average traces

absMeanTraces = [];AvgMaxCurrentAVGMinus=[]; AvgMaxCurrentAVG=[];CellMinAVG=[];MinAVG=[];
 MinAOffAVG=[];AVGCellMinOff=[];CellMinOffAVGShort=[];
MaxActuSensorOnAVG =[];CellMaxActuFirstAVG = []; StartOffBelowShortAVG =[];StartOffBelowAVG =[];
AvgMaxCurrentOffAVG = [];

absMeanTraces = abs(MeanTraces);

for i = 1:size(absMeanTraces,2);
MinAVG(i) = max(absMeanTraces(Start(i):Ende(i),i));
CellMinAVG(i) = find([absMeanTraces(:,i)] == MinAVG(i),1,'first');
AvgMaxCurrentAVG(i) = mean(MeanTraces(CellMinAVG(i)-5:CellMinAVG(i)+5,i)); % Average from averaged traces and not average of the single peaks
StartOffBelowShortAVG(i) = find([ActuSensorAvg(MeanCellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first'); % %% find cell, where 1st value is lower than threshold; Onset Off-Stimulus
StartOffBelowAVG(i) = StartOffBelowShortAVG(i)+ MeanCellMaxActuFirst(i);
MinAOffAVG(i) = max(absMeanTraces(StartOffBelowAVG(i)-0.05*fs:StartOffBelowAVG(i)+0.01*fs,i)); %change not hardcode 
CellMinOffAVGShort(i) = find([absMeanTraces(StartOffBelowAVG(i)-0.05*fs:end,i)] == MinAOffAVG(i),1,'first');
CellMinOffAVGShort(i) = find([absMeanTraces(MeanCellMaxActuFirst(i):StartOffBelowAVG(i)+0.01*fs,i)] == MinAOffAVG(i),1,'first');  
CellMinOffAVG(i) = CellMinOffAVGShort(i)+MeanCellMaxActuFirst(i);
AvgMaxCurrentOffAVG(i) = mean(MeanTraces(CellMinOffAVG(i)-5:CellMinOffAVG(i)+5,i)); % ToDo Five or 10???
end

%calculate Mean On & Off Currents for Force dependence
absMeanTracesSortForce = [];
absMeanTracesSortForce = abs(MeanTracesSortForce);
MinAVGSortForce=[];CellMinAVGSortForce=[];AvgMaxCurrentAVGSortForce=[];
MinAOffAVGSortForce=[];CellMinOffAVGShortSortForce=[];AvgMaxCurrentOffAVGSortForce=[];CellMinOffAVGSortForce=[];
for i = 1:size(absMeanTracesSortForce,2);
MinAVGSortForce(i) = max(absMeanTracesSortForce(Start(i):Ende(i),i));
CellMinAVGSortForce(i) = find([absMeanTracesSortForce(:,i)] == MinAVGSortForce(i),1,'first');
AvgMaxCurrentAVGSortForce(i) = mean(MeanTracesSortForce(CellMinAVGSortForce(i)-5:CellMinAVGSortForce(i)+5,i)); % Average from averaged traces and not average of the single peaks
%StartOffBelowShortAVG(i) = find([ActuSensorAvg(MeanCellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first'); % %% find cell, where 1st value is lower than threshold; Onset Off-Stimulus
% StartOffBelowAVG(i) = StartOffBelowShortAVG(i)+ MeanCellMaxActuFirst(i);
MinAOffAVGSortForce(i) = max(absMeanTracesSortForce(StartOffBelowAVG(i)-0.05*fs:StartOffBelowAVG(i)+0.01*fs,i)); %change not hardcode 
CellMinOffAVGShortSortForce(i) = find([absMeanTracesSortForce(StartOffBelowAVG(i)-0.05*fs:end,i)] == MinAOffAVGSortForce(i),1,'first');
CellMinOffAVGShortSortForce(i) = find([absMeanTracesSortForce(MeanCellMaxActuFirst(i):StartOffBelowAVG(i)+0.01*fs,i)] == MinAOffAVGSortForce(i),1,'first');  
CellMinOffAVGSortForce(i) = CellMinOffAVGShortSortForce(i)+MeanCellMaxActuFirst(i);
AvgMaxCurrentOffAVGSortForce(i) = mean(MeanTracesSortForce(CellMinOffAVGSortForce(i)-5:CellMinOffAVGSortForce(i)+5,i));
end

AvgMaxCurrentAVG = AvgMaxCurrentAVG*-1
AvgMaxCurrentAVG = AvgMaxCurrentAVG';
AvgMaxCurrentOffAVG=AvgMaxCurrentOffAVG*-1
AvgMaxCurrentOffAVG = AvgMaxCurrentOffAVG'
MeanIndentationVer=MeanIndentation';
MeanForceVer=MeanForce';
AvgMaxCurrentMinus = AvgMaxCurrentMinus';
AvgMaxCurrentMinusOff =AvgMaxCurrentMinusOff';

AvgMaxCurrentAVGSortForce = AvgMaxCurrentAVGSortForce *-1;
AvgMaxCurrentAVGSortForce = AvgMaxCurrentAVGSortForce';
AvgMaxCurrentOffAVGSortForce =AvgMaxCurrentOffAVGSortForce*-1
AvgMaxCurrentOffAVGSortForce = AvgMaxCurrentOffAVGSortForce';

%figure()
%plot(CellMaxActuFirst)

figure()
plot(MeanTracesppA) % plot Current in pA
%xlim([0 0.6])
ylabel('Current (pA)')
xlabel('Time')
title('Mean Traces of certain cell')
legend(name)

MeanTracesSortForceppA = MeanTracesSortForce*10^12;
figure()
plot(MeanTracesSortForceppA) % plot Current in pA
%xlim([0 0.6])
ylabel('Current (pA)')
xlabel('Time')
title('Mean Traces of certain cell calcuated for Force')
legend(name)


StiffnessWorm = MeanIndentationVer\MeanForceVer % is doing a leastsquarefit
MeanLeak = nanmean(LeakANew)
SDLeak = nanstd(LeakANew)



%%
%%%% Figure Summary Analysis

xScatter = (1:length(MeanIndentation));
figure()
subplot(3,2,1)
scatter(xScatter, LeakAppA) 
title('control: Leak Current')
ylabel('Current (pA)')
xlabel('number of file (in recorded order)')


hold on 
subplot(3,2,2)
if isFiveStep == 1 || isFiveRamp == 1 || isFifteenStep == 1;
    i = 1;
while i <= length(MeanIndentation)
scatter(MeanIndentation(i:i+4), AverageMaxCurrentMinusppA(i:i+4),'LineWidth',2)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Cur vs Ind')
xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Indentation (um)')
%hold on 
% for j = 1:size(Aall,2)/5
% %legend(Files(j))  % include legend again
% end
end



else
scatter(MeanIndentation, AverageMaxCurrentMinusppA,'LineWidth',2)
title('Mean Cur vs Ind')
xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Indentation (um)')
end

hold on 
subplot(3,2,5)
scatter(MergeInd, AvgMaxCurrentAVG)
hold on
subplot(3,2,5)
scatter(MergeInd, AvgMaxCurrentOffAVG)
legend('ON-Current','OFF-Current')
hold on 
subplot(3,2,6)
scatter(MergeForce, AvgMaxCurrentAVGSortForce)
hold on 
subplot(3,2,6)
scatter(MergeForce, AvgMaxCurrentOffAVGSortForce)
legend('CurrentOn','CirrentOff')



hold on 
subplot(3,2,3)
scatter(xScatter, Start) 
%ylim([100 1000]) % ToDo: change it for Ramps
ylabel('Point')
xlabel('Filenumber')
title('control: Find OnSet of Stimulus')
hold on 
subplot(3,2,3)
scatter(xScatter, StartOffBelow) %change to start of below
%ylim([100 1000]) % ToDo: change it for Ramps
ylabel('Point')
xlabel('Filenumber')
title('control: Find OnSet of Stimulus')
hold on 
subplot(3,2,4)
if isFiveStep == 1 || isFiveRamp == 1 || isFifteenStep == 1;
    i = 1;
while i <= length(MeanIndentation)
scatter(Velocity(i:i+4), AverageMaxCurrentMinusppA(i:i+4),'LineWidth',2)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Cur vs Vel')
xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Velocity(um/s)')
%hold on 
% for j = 1:size(Aall,2)/5
% %legend(Files(j))  % include legend again
% end
end
else
scatter(Velocity, AverageMaxCurrentMinusppA,'LineWidth',2)
title('Mean Cur vs Ind')
%xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Velocity(um/s)')
end
%xlabel('Velocity')
%title('control: Current over time')
% hold on
% xscatterMergeInd = length(MergeInd)
% 
% set(p1,'markerfacecolor','r')
% set(p2,'markerfacecolor','g')
%ylim([100 1000]) % ToDo: change it for Ramps
% ylabel('Current (pA)')
hold on
xlabels{1} = 'Velocity (um/s)';
xlabels{2} = '';
ylabels{1} = 'Current (pA)';
ylabels{2} = 'Ind (um)';
subplot(3,2,5)
[ax,hl1,hl2] = plotxx(Velocity,AverageMaxCurrentMinusppA,Velocity,MeanIndentation,xlabels,ylabels)

%%%% maybe useful later
% figure()
% %subplot(3,2,5)
% i = 1
% while i <= length(MergeInd)
% scatter(FindSameIndValues(:,i), FindSameIndCurrent(:,i),'LineWidth',2)
% hold on
% i = i+1;
% end


%%
hold on 
subplot(3,2,4)
plot(Time,ASubtractppA)
%xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Time (s)')
title('Current')

%%% plotting ForceClamp signals in a subplot
allRiseTimeInms = allRiseTime*1000;

figure()
subplot(3,3,1)
plot(Time,CantiDefl)
%xlim([0 0.6])
xlabel('Time (s)')
title('Cantilever Deflection')
ylabel('Deflection (µm)')
hold on
subplot(3,3,2)
plot(Time,Indentation)
%xlim([0 0.6])
xlabel('Time (s)')
ylabel('Indentation (µm)')
title('Indentation')

hold on
subplot(3,3,3)
plot(Time,Force)
xlim([0 0.6])
xlabel('Time (s)')
ylabel('Force (µN)')
title('Force')
hold on
subplot(3,3,4)
plot(Time,ActuSetPoint)
%xlim([0 0.6])
xlabel('Time (s)')
ylabel('Displacement (µm)')
title('Displacement ActuSetPoint')
hold on
subplot(3,3,5)
plot(Time,ActuSensor)
%xlim([0 0.6])
xlabel('Time (s)')
ylabel('Displacement (µm)')
title('Displacement ActuSensor')
hold on
subplot(3,3,6)
plot(Time,normCantiDefl)
%xlim([0 0.6])
xlabel('Time (s)')
ylabel('normalized Deflection')
title('Cantilever Defl norm')
hold on
subplot(3,3,7)
scatter(MeanIndentation, allRiseTimeInms)  
%xlim([0 0.3])
title('RiseTime (CantiDefl)')
ylabel('Rise Time Tau (ms)')
xlabel('Indentation')
xlim([0 max(MeanIndentation)+1])
 hold on
 subplot(3,3,8)
 scatter(MeanIndentation, allOvershoot)  
 %xlim([0 max(MeanIndentation)+1])
 ylabel('% to steady state')
 xlabel('Indentation (µm)')
 title('Overshoot (CantiDefl)')
hold on 
subplot(3,3,9)
i = 1;
while i <= length(MeanIndentation)
scatter(MeanIndentation(i:i+4), MeanForce(i:i+4),'LineWidth',2)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Force vs Ind')
xlim([0 max(MeanIndentation)+1])
ylabel('Force (µN)')
xlabel('Indentation (µm)')
hold on 
for j = 1:size(Aall,2)/5
%legend(Files(j))  % include legend again
end
end


MeanTauIndentation = mean(allRiseTimeInms)
display 'include TauCalculation into delete singel series'

%ToDO: Nr of Avegrages wrong!
%rToDo: MergeIndRow = num2str(MergeIndRow);%how to write each Indentation as col header???



%%
%Export Data
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1; 
% ToDo change for Current ExportData = [MergeInd,MeanSameIndCurrent,NormMeanCurrent,MeanSameIndForce];
ExportMeanSameInd = [MergeInd,MeanSameIndForce,AvgMaxCurrentAVG,AvgMaxCurrentAVG,NumberOfAvergagesPerInd,AvgMaxCurrentOffAVG,AvgMaxCurrentOffAVG];
ExportMeanSingle = [MeanIndentationVer,MeanForceVer,AvgMaxCurrentMinus,AvgMaxCurrentMinusOff];
ExportMeanSameIndOFF = [NumberOfAvergagesPerInd,AvgMaxCurrentOffAVG,AvgMaxCurrentOffAVG];
ExportMeanSortForce = [MergeForce,AvgMaxCurrentAVGSortForce,AvgMaxCurrentOffAVGSortForce];
else
%XxportData = [MergeVel,MeanSameVelIndentation,MeanSameIndCurrent,NormMeanCurrent,MeanSameIndForce];  %MeanSameIndCurrent means here: current at same velocities
ExportMeanSameVelDataMechanics = [MergeVel,MeanSameVelStiffnessRampVer,MeanSameVelIndentation];
ExportMeanDataMechanics = [VelocityVer,StiffnessRampVer,MeanIndentationVer,MeanForceVer];
end

%%% write Matlabvariables
if isFiveStep == 1 || isStep == 1;
save(sprintf('FiveStep-%s.mat',name)); %save(sprintf('%sTEST.mat',name))
%if  isFifteenStep == 1;
    %save(sprintf('FifteenStep-%s.mat',name));
else 
save(sprintf('Ramp-%s.mat',name));
%end
end

%%%ToDO: include if RampHold, save RampSTF00X, if Step, StepXXX, otherwise,
%%%it overwrites the analysis

%%% write as csv, because cannot write with mac to excel

if isFiveStep == 1 || isStep == 1;
    %disp 'csv save file has to be written'
   % MergeIndRow = MergeInd';
%save Means of same indentation in csv file

filename = sprintf('StepSameInd-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'MergeInd-%s, MeanSameIndForce-%s, SameIndCurrent-%s,SameIndCurrentCOPY-%s, NrAVGperIndentation-%s,SameIndCurrentOFF-%s,SameIndCurrentOFF-COPY-%s \n',name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanSameInd, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.


filename = sprintf('StepSameIndOFF-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'NrAVGperIndentation-%s,SameIndCurOFF-%s, SameIndCurOFFCOPY-%s \n',name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanSameIndOFF, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('StepSortForce-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'MergeForce-%s,OnCurSortForce-%s, OffCurSortForce-%s \n',name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanSortForce, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

%save single Indentation values in separate csv file
filename = sprintf('StepSingleInd-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'Ind-%s, Force-%s, Cur-%s,CurOff-%s \n',name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanSingle, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

% how to include the Filenumber?
filename = sprintf('TracesAvgCurrent-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid,'CurInd1-%s, CurInd2-%s, CurInd3-%s, CurInd4-%s, CurInd5-%s, CurInd6-%s, CurInd7-%s, CurInd8-%s, CurInd9-%s, CurInd10-%s, CurInd11-%s, CurInd12-%s, CurInd13-%s, CurInd14-%s, CurInd15-%s, CurInd16-%s, CurInd17-%s, CurInd18-%s, CurInd19-%s \n',name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, MeanTraces, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('TracesAvgIndentation-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid,'Ind1-%s, Ind2-%s, Ind3-%s, Ind4-%s, Ind5-%s, Ind6-%s, Ind7-%s, Ind8-%s, Ind9-%s, Ind10-%s, Ind11-%s, Ind12-%s, Ind13-%s, Ind14-%s, Ind15-%s, Ind16-%s, Ind17-%s, Ind18-%s, Ind19-%s \n',name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, MeanTracesIndentation, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('TracesAvgForce-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid,'Force1-%s, Force2-%s, Force3-%s, Force4-%s, Force5-%s, Force6-%s, Force7-%s, Force8-%s, Force9-%s, Force10-%s, Force11-%s, Force12-%s, Force13-%s, Force14-%s, Force15-%s, Force16-%s, Force17-%s, Force18-%s, Force19-%s \n',name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, MeanTracesForce, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.


else %isFiveRamp == 1;
    
filename = sprintf('RampSameVel-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'MergeVel-%s, SameVelStiffness-%s, SameVelInd-%s \n',name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanSameVelDataMechanics, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('RampSingleVel-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'Vel-%s, Stiffness-%s, Ind-%s, Force-%s \n',name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportMeanDataMechanics, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

end

%%
% reload notes to update for FitMax values, Normalize traces

loadFileMode = 1; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('Ephys-Meta-Sylvia.xlsx'); % be careful in which Folder saved.
end


FitMax =[];FitMaxForce = [];
headers = raw(1,:);
indFitMax = find(strcmpi(headers, 'ONFitMax(A)'));
FitMax = raw(Stiffrow,indFitMax); %Stiffrow defined previously
FitMax = cell2mat(FitMax);
indFitMaxForce = find(strcmpi(headers, 'ONFitMaxForce(A)'));
FitMaxForce = raw(Stiffrow,indFitMaxForce); 
FitMaxForce = cell2mat(FitMaxForce);

AvgMaxCurrentNormInd = []; AvgMaxCurrentNormForce = [];
AvgMaxCurrentNormInd = AvgMaxCurrentAVG/FitMax;
AvgMaxCurrentNormForce = AvgMaxCurrentAVG/FitMaxForce;

ExportNormTraces = [];
%%% update .mat file
if isFiveStep == 1 || isStep == 1;
save(sprintf('FiveStep-%s.mat',name)); %save(sprintf('%sTEST.mat',name))
%if  isFifteenStep == 1;
    %save(sprintf('FifteenStep-%s.mat',name));
ExportNormTraces = [AvgMaxCurrentNormInd,AvgMaxCurrentNormForce];
filename = sprintf('NormAVGCurrent-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'NormCurFitInd-%s, NormCurFitForce-%s \n',name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportNormTraces, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

else 
save(sprintf('Ramp-%s.mat',name));
%end
end
%AvgMaxCurrentNorm = AvgMaxCurrentNorm'; 



%%

% int_cols = all(isnan(MeanIndentationNew)|round(MeanIndentationNew)==MeanIndentationNew,1);
% it = MeanIndentationNew(:,int_cols);
% 
% test = [1   NaN   2.2   3.2  4;
%      NaN 7.9   5.1   NaN  5;
%      3    5.5  NaN   4.1  NaN];
% int_cols = all(isnan(MeanIndentationNew)|round(MeanIndentationNew)==MeanIndentationNew,1);
% it = MeanIndentationNew(:,int_cols);
% flt = MeanIndentationNew(:,~int_cols);

% write excel sheet
% col_header={name,'MeanInd','MeanCurrent','NormMeanCurrent','MeanForce','TracesPerIndentation'};     %Row cell array (for column labels)
% %row_header(1:10,1)={'Time'};     %Column cell array (for row labels)
% xlswrite(name,MergeInd,'Sheet1','B2');     %Write data
% xlswrite(name,MeanSameIndCurrent,'Sheet1','C2');     %Write data
% xlswrite(name,NormMeanCurrent,'Sheet1','D2');     %Write data
% xlswrite(name,MeanSameIndForce,'Sheet1','E2');     %Write data
% xlswrite(name,TracesPerIndentation,'Sheet1','F2'); 
% xlswrite(name,col_header,'Sheet1','A1');     %Write column header
% %xlswrite('My_file.xls',row_header,'Sheet1','A2');      %Write row header
% col_header2={name,'Ind1','Ind2','Ind3','Ind4','Ind5','Ind6','Ind7','Ind8','Ind9','Ind10','Ind11','Ind12','Ind13','Ind14'}; %ToDO - get the values for the
% %Indentations
% xlswrite(name,MeanTraces,'Sheet2', 'B2');   
% xlswrite(name,col_header2,'Sheet2','B1'); 

%find number of nan values to find out the number of averages
%average traces 
%how to average two columns
%  for k = 1:length(MergeInd);
%      for i = 1:length(FindSameInd(:,k));
% % MeanTracesCurrent(:,k) = nanmean(SortASubtract(FindSameInd(:,k)))
% %B(:,nn+1) MeanTracesCurrent = nanmean(SortASubtract(:,1:2))
%      end
%  end




