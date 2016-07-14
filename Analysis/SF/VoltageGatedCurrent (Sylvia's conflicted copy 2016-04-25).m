%%  load dat.files 
clear all; close all; clc;
ephysData=ImportPatchData();
%%
% load notes to get several values automatically needed for the conversion of the signals
loadFileMode = 0; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('Ephys-Meta-Sylvia.xlsx'); % be careful in which Folder saved.
end

%% Analysis Indivisual Recording 
close all; clc

%%% hardcoding part:
%%% change:  STF0XX, sampling fequency not yet fully automatic %%%%%% 

name = 'STF026'; % name of recording. placed into varaibel fiels names%
stimuli = 'OC_IVq'; % Single protocols: Step and Ramp-Hold; Five sweeps per protocol:FiveStep, FiveRampHold; does not work with alternative names
stimuli2 = 'WC_IVq';
Filenumber = 1; % wil be used to extract sampling freuqnecy; first file loaded, maybe change

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


A2 = A;
B2 = B;
C2 = C;
D2 = D;


% find all files with certain protocol name: 
% ifelse statement: if not respective stimuli name (FiveStep or FiveRampHold), then empty array; 
for i = Files(:,1):Files(:,end);
   if find(strcmpi(ephysData.(name).protocols{1,i}, stimuli)) == 1;
        continue
   else 
         A{1, i} = []; B{1, i} = []; C{1, i} = []; D{1, i} = [];          
    end      
end

% for stimuli 2
for i = Files(:,1):Files(:,end);
   if find(strcmpi(ephysData.(name).protocols{1,i}, stimuli2)) == 1;
        continue
   else 
         A2{1, i} = []; B2{1, i} = []; C2{1, i} = []; D2{1, i} = [];          
    end      
end

AllStimuliBlocks = (find(strcmpi(ephysData.(name).protocols, stimuli)))
AllStimuliBlocks2 = (find(strcmpi(ephysData.(name).protocols, stimuli2)))

% A{1,3}(:,1) % access column 1
AvgOCIVq = (A{1,AllStimuliBlocks(:,1)} + A{1,AllStimuliBlocks(:,2)} + A{1,AllStimuliBlocks(:,3)})/3;
AvgWCIVq = (A2{1,AllStimuliBlocks2(:,1)} + A2{1,AllStimuliBlocks2(:,2)} + A2{1,AllStimuliBlocks2(:,3)})/3;
AvgWCcorrected = AvgWCIVq - AvgOCIVq;

MeanAvgWC=[];
for k = 1:size(AvgWCcorrected,2);
   MeanAvgWC(k) = mean(AvgWCcorrected(150:500,k));
   end
MeanAvgWC = (MeanAvgWC)';

Voltage=[-80;-60;-40;-20;0;20;40;60;80];

figure()
subplot(2,3,1)
plot(AvgOCIVq)
title('control: OC IV')
hold on
subplot(2,3,2)
plot(AvgWCIVq)
title('control: WC')
ylabel('Current (pA)')
hold on
subplot(2,3,3)
plot(AvgWCcorrected)
title('control: corrected')
ylabel('Current (pA)')
%xlabel('number of file (in recorded order)')
hold on
subplot(2,3,4)
scatter(Voltage, MeanAvgWC)

