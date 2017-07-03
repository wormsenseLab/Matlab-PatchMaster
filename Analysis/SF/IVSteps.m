%%  load dat.files 
clear all; close all; clc;
ephysData=ImportPatchData();
% load('ephysdata(20170130).mat')

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
name = 'STF080'; % name of recording. placed into varaibel fiels names%
stimuli = 'IVStep'; 
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
AllStimuliBlocks = (find(strcmpi(ephysData.(name).protocols, stimuli)))


%%%%
%deleting whole blocks of FiveBlockStimuli; Whole block=Filenumber
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

% removes all empty cells from the cell array
AShort = A(~cellfun('isempty',A)); BShort = B(~cellfun('isempty',B)); CShort = C(~cellfun('isempty',C)); DShort = D(~cellfun('isempty',D));

Aall = []; Aall = cat(2,AShort{:}); %concatenate all values from A to a double

%figure to check if stimulus is incomplete, if yes, delete inclomplete
%stmulusblock
for i = 1:size(Aall,2)
subplot(ceil(size(Aall,2)/10),10,i)
plot(Aall(:,i))
%ylim([-1 16])
end



%%
ActuSensor = [];
ActuSensorCellArray = cellfun(@(x) x*1.5,BShort,'un',0); %multiply all values in the cell with 1.5
%ActuSensor = mean(cat(3,ActuSensorCellArray{:}),3); %but I want to visualize it first
%TakemeanAcrossCellA = mean(cat(3,AShort{:}),3); 

ActuSensorB3D = cat(3, ActuSensorCellArray{:}); %1strows, 2nc col, 3rd different series
Aall3D = cat(3, AShort{:}); % current
Call3D = cat(3, CShort{:}); % cantilever signal
Dall3D = cat(3, DShort{:}); % Actuator Setpoint

ALeak = mean(Aall3D (2000:2200,:,:)); %leak current % toDo not hard coded;
ASubtract = bsxfun(@minus, Aall3D(:,:,:), ALeak);


fs = ephysData.(name).samplingFreq{1, Files(:,AllStimuliBlocks(1))}; % sampling frequency from first Stimuli loaded; 
interval = 1/fs;   %%% get time interval to display x axis in seconds 
ENDTime = length(Aall)/fs; %%% don't remember why I complicated it
Time = (0:interval:ENDTime-interval)'; 


%scatter = size(ALeak,2)*size(ALeak,3)
% figures to control leak
ASubtractppA = bsxfun(@times, ASubtract(:,:,:), 10^12);
figure()
 for i = 1:size(ASubtract,3)
   subplot(size(ASubtract,3),1,i)
     plot(ASubtract(:,:,i))
    xlim([2500 7000])
 end
 
 Aall3Dppa = bsxfun(@times, Aall3D(:,:,:), 10^12);
 figure()
 for i = 1:size(Aall3D,3)
   subplot(size(Aall3D,3),1,i)
     plot(Aall3Dppa(:,:,i))
   xlim([2500 7000])
 end
 
ALeakppA = bsxfun(@times, ALeak(:,:,:), 10^12); 
figure()
%plot(ALeak(1,:,:)) %scatter does not work ?
 for i = 1:size(Aall3D,3)
   subplot(size(Aall3D,3),1,i)
     plot(ALeakppA(:,:,i))
 end
 
% Avg over 3rd dimension (combine all series to one)
AvgIVq = mean(ASubtract,3);
ActuSensor = mean(ActuSensorB3D,3);

figure()
plot(Time, AvgIVq)
figure()
plot(Time, ActuSensor)

figure()
plot(AvgIVq)

% calculate the threshold for detecting the onset of the stmulus


[SlopeActu,MaxZeroSlopeActu,StdZeroSlopeActu,MaxZeroActu,StdZeroActu,MaxActuSensorPlateau,StdActuSensorPlateau,CellMaxActuFirst] = SlopeThreshold(ActuSensor); 

StartBase = [];
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1 || isIVStep == 1;
   StartBase = MaxZeroActu + 4*StdZeroActu;   %% play around and modifz 
else
   StartBase = MaxZeroSlopeActu + 4*StdZeroSlopeActu; %
end



      for i = 1:size(AvgIVq,2);
    Start(i) = find([ActuSensor(:,i)] > StartBase(i),1, 'first'); %% find cell, where 1st value is bigger than threshold; Onset On-Stimulus 
    Ende(i) = Start(i) + (fs/20);
   % EndeRamp(i)= StartOffBelow(i)-(0.300/interval); % Time= Points*interval
      end
  
       
absMeanTraces = [];AvgMaxCurrentAVGMinus=[]; AvgMaxCurrentAVG=[];CellMinAVG=[];MinAVG=[];
 MinAOffAVG=[];AVGCellMinOff=[];CellMinOffAVGShort=[];
MaxActuSensorOnAVG =[];CellMaxActuFirstAVG = []; StartOffBelowShortAVG =[];StartOffBelowAVG =[];
AvgMaxCurrentOffAVG = [];


MaxASub = []; StdASub =[]
MaxASub = max(AvgIVq(1500:1600,:,:));
StdASub = std(AvgIVq(1500:1600,:,:));
%ASubtract = bsxfun(@minus, Aall3D(:,:,:), ALeak);

ThresholdCur = MaxASub + 2*StdASub;

absMeanTraces = abs(AvgIVq);

for i = 1:size(absMeanTraces,2);
MinAVG(i) = max(absMeanTraces(Start(i):Ende(i),i));
CellMinAVG(i) = find([absMeanTraces(:,i)] == MinAVG(i),1,'first'); %error, if noise before
AvgMaxCurrentAVG(i) = mean(AvgIVq(CellMinAVG(i)-5:CellMinAVG(i)+5,i)); % Average from averaged traces and not average of the single peaks
%StartOffBelowShortAVG(i) = find([ActuSensorAvg(MeanCellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first'); % %% find cell, where 1st value is lower than threshold; Onset Off-Stimulus
%StartOffBelowAVG(i) = StartOffBelowShortAVG(i)+ MeanCellMaxActuFirst(i);
%MinAOffAVG(i) = max(absMeanTraces(StartOffBelowAVG(i)-0.05*fs:StartOffBelowAVG(i)+0.01*fs,i)); %change not hardcode 
%CellMinOffAVGShort(i) = find([absMeanTraces(StartOffBelowAVG(i)-0.05*fs:end,i)] == MinAOffAVG(i),1,'first');
%CellMinOffAVGShort(i) = find([absMeanTraces(MeanCellMaxActuFirst(i):StartOffBelowAVG(i)+0.01*fs,i)] == MinAOffAVG(i),1,'first');  
%CellMinOffAVG(i) = CellMinOffAVGShort(i)+MeanCellMaxActuFirst(i);
%AvgMaxCurrentOffAVG(i) = mean(absMeanTraces(CellMinOffAVG(i)-5:CellMinOffAVG(i)+5,i)); % ToDo Five or 10???
end


Voltage=[-100;-80;-60;-40;-20;0;20;40;60;80];
AvgMaxCurrentAVG = AvgMaxCurrentAVG'

headers = raw(1,:);
FindRowIndCellId = strcmpi(raw,name); % name = recorded cell
[RowCellId,col] = find(FindRowIndCellId,1);
Rs = []; TermRsI = []; MinusRs =[];
indRs = find(strcmpi(headers, 'Rs(MOhm)')); % find col with Sensitivity
Rs = raw(RowCellId,indRs); 
Rs = cell2mat(Rs)
RsinOHM = Rs*10E6;

%MinusRs = RsinOHM *-1;
TermRsI = RsinOHM * AvgMaxCurrentAVG;
%TermRsI = (TermRsI)';
VoltageInV=[-0.10;-0.080;-0.060;-0.040;-0.020;0;0.020;0.040;0.060;0.080];
Vcorrected = VoltageInV-TermRsI; %minus(Voltage,TermRsI); % Vcorrected =  Vcom -Rs*I;
VcorInMV = Vcorrected*1000;

figure()
plot(Voltage,AvgMaxCurrentAVG)
hold on
plot(VcorInMV,AvgMaxCurrentAVG)

figure()
plot(AvgIVq)

%% Export
save(sprintf('IVSteps-%s.mat',name)); %save(sprintf('%sTEST.mat',name))

ExportIVSteps = [];
ExportIVSteps = [Voltage,VcorInMV,AvgMaxCurrentAVG];

filename = sprintf('IVSteps-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'Vcom-%s, Vcorrected-%s, AvgMaxCurON-%s \n',name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportIVSteps , '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

