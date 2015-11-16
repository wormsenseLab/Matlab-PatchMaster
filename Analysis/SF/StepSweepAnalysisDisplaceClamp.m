%%% Sylvia Fechner
%%% Stanford University, 
%%% 20151115
%%% Script to analyze data from FALCON in Displacement Clamp
%%%% Five sweeps of steps within one series

%%  load dat.files 
clear all; close all; 
ephysData=ImportPatchData()

%% load notes to get several valus automatically needed for the coversion of the signals
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);

%% StepAnalysis 
%%% TO DO:
%%% check, if SetPoint and Sensor are pluged in correctly
%%% running average
%%% find beginning of step
%%% choose which traces of one recording are good
%%% average good recordings
close all % closes all figures
%%% hardcoding part:
%%% change: Files; STF0XX; rawStiff %%%%%% 
%%% ToDo -Analyze all steps -- look for step 
Files = 12:15; % Patchmaster File Numbers
name = 'STF012'; % name of recording. placed into varaibel fiels names
rawStiff = 2; %% To Do: change to find automatically %% change if two different cantilevers were used

%%% load Current, Actuator Sensor, And Cantilver Signal of each step%
A=[];B=[]; C=[];   %% to get an empty array, if I analyzed different files before
for i = Files(:,1):Files(:,end);
A{i} = ephysData.(name).data{1, i}; %Current
B{i} = ephysData.(name).data{3, i}; % actuator Sensor
C{i} = ephysData.(name).data{4, i}; % Cantilever Signal
D{i} = ephysData.(name).data{2, i}; % actuator SetPoint
end
Aall = cat(2,A{:});
Ball = cat(2,B{:});
Call = cat(2,C{:});
Dall = cat(2,D{:});

%%%%% after loading data
fs = ephysData.(name).samplingFreq{1, Files(:,1)}; %% sampling frequency from first file loaded; I currently assume it will be the same
interval = 1/fs;   %%%%% get time interval to display x axis in seconds
ENDTime = length(Aall)/fs; %% don't remember why I complicated it
Time = (0:interval:ENDTime-interval)'; 

%%%%%% Subtract leak current
LeakA = [];
ASubtract = [];
 for i = 1:size(Aall,2);
LeakA(i) = mean(Aall(1:100,i));  %%%% take the mean of the first 100 Points
ASubtract(:,i) = Aall(:,i) - LeakA(i); %%%% subtract leak  --> copy to Igor
 end

%%% getting max current (min value) for On response (needs to be done differently)

for i = 1:size(Ball,2);
   AverOnset(i)=  mean(Ball(1:100,i));
    MaxZeroActu(i)= max(Ball(1:100,i));
   StdZeroActu(i)= std(Ball(1:100,i));
end
StartBase = MaxZeroActu  + 2*StdZeroActu; % set a threshold to find the onset of the stimulus

MinA = [];
CellMin = [];
AverageMaxCurrent = [];
AverageMaxCurrentMinus = [];
for i = 1:size(Aall,2);
Start(i) = find([Ball(:,i)] > StartBase(i),1,'first'); %% find cells, where 1st values is bigger than threshold 
Ende(i) = Start(i) + 100; %% could change 100 to be dependent on fs 
MinA(i) = min(ASubtract(Start(i):Ende(i),i));
CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first'); % find cell with min value
Values = ASubtract(:,i);
%%% TODO: include taking the mean from an averaged signal, check peak
%%% finding function in igor
AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
 AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to falsify demontrating of increase in current
end
 
figure() % plot all Current leak Subtracted signals in one plot
plot(Time,ASubtract)
title('Current')


%%%% SetPoint
ActuSetPoint = []; SetPointDispl = [];
for i = 1:size(Dall,2),
ActuSetPoint(:,i) = Dall(:,i)*1.5;
SetPointDispl(i) = max(ActuSetPoint(:,i));
end

figure()  %%% plots all Cantilever deflections in one plot
plot(Time,ActuSetPoint)
title('Displacement SetPoint')

%%% to get Displacment of Actuator: multiply actuator sensor signal times 
%%% sensitivity of actuator: 1.5 %%%
ActuSensor = []; AmplitudeDispl = [];
for i = 1:size(Ball,2),
ActuSensor(:,i) = Ball(:,i)*1.5;
AmplitudeDispl(i) = max(ActuSensor(:,i));
%figure()
%figure()
% plot(Time,ActuSensor(:,i));
% title('Displacement Actuator')
end

figure()  %%% plots all Cantilever deflections in one plot
plot(Time,ActuSensor)
title('Displacement Actuator')

%%%%% to get Deflection of Cantilever: multiply with Sensitivity 
%%%% get Sensitivity from Notes of Recording day %%%%

headers = raw(1,:);
ind = find(strcmpi(headers, 'Sensitivity(um/V)'));
Sensitivity = raw(rawStiff,ind);  %%% rawstiff = if different cantilever were used change in hardcoding part
Sensitivity = cell2mat(Sensitivity);


Baseline =[];CantiZero = [];CantiDefl = [];
for i = 1:size(Call,2);
Baseline(i) = mean(Call(1:100,i));
CantiZero(:,i) = Call(:,i)-Baseline(i);
CantiDefl(:,i) = CantiZero(:,i)*Sensitivity;
%figure()
 %plot(Time,CantiDefl(:,i));
% title('Cantilever Deflection')
end

figure()  %%% plots all Cantilever deflections in one plot
plot(Time,CantiDefl)
title('Cantilever Deflection')

%%%% calculate Indentation = Actuator Sensor - Cantilever Deflection
for i = 1:size(Call,2);
Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(1000:2000,i));
end

figure() %%% plots all Indentations in one plot
plot(Time,Indentation)
title('Indentation')

%%% Calculating Force Signal 
indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
Stiffness = raw(rawStiff,indStiffness); 
Stiffness = cell2mat(Stiffness);

for i = 1:size(Aall,2);
Force(:,i) = CantiDefl(:,i) *Stiffness;
% figure()
% plot(Time,Force);%
% title('Force')
end
%%% 

figure() %%% plots all Force signals in one plot
plot(Time,Force)
title('Force')

Amplitude = [];
Amplitude =[AmplitudeDispl',MeanIndentation',AverageMaxCurrent', AverageMaxCurrentMinus'];%, (Files(:,1):Files(:,end))'];%; MeanIndentation(i)'];%, MeanIndentation(i), AverageMaxCurrent(i), AverageMaxCurrentMinus(i), i];

%%%%%% now make all plots
% %%%%%%% figures current with and without leak subtraction %%%%
figure()
for i = 1:size(Aall,2)
subplot(length(Files),5,i)
plot(Time,Aall(:,i))
ylim([-5*10^-11 1*10^-11])
hold on
plot(Time,ASubtract(:,i))
title(MeanIndentation(i))
end

xScatter = (1:length(MeanIndentation));
figure()
scatter(xScatter, LeakA) 
title('control: Leak Current')

xScatter = (1:length(MeanIndentation));
figure()
scatter(xScatter, Start) 
ylim([0 1000])
title('control: to see if thresholding is working')

figure()
i = 1;
while i <= length(MeanIndentation)
scatter(MeanIndentation(i:i+4), AverageMaxCurrentMinus(i:i+4),'LineWidth',10)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Current vs Indentation')
hold on 
for j = 1:length(Files)
legend(Files(j))
end
end

%% Average of variable numbers of columns
while 1
prompt = {'FirstValue','SecondValue','ThirdValue','InputReady'};
IndValues = inputdlg(prompt)
FirstValue = str2num(IndValues{1})
SecondValue = str2num(IndValues{2})
ThirdValue = str2num(IndValues{3})
a = IndValues{4}
b = ('done')

%xdatatemp = ASubtract(:,[1 6])
if isempty(ThirdValue) == 1
  ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue])
else
   ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue ThirdValue]) 
end
%elseif ~strcmp(a,b);
 %  msgbox('it is not done')
  if  strcmpi(a,b)
     msgbox('it is done');
     break;
  end
 end


% 
% prompt = {'FirstValue','SecondValue','ThirdValue'};
% IndValues = inputdlg(prompt)
% FirstValue = str2num(IndValues{1})
% SecondValue = str2num(IndValues{2})
% ThirdValue = str2num(IndValues{3})
% %xdatatemp = ASubtract(:,[1 6])
% if isempty(ThirdValue) == 1
%   ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue])
% else
%    ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue ThirdValue]) 
% end
% 

%%% Average current; but it is not leak subtracte
% E = []; F = []; AverInd = []
% for i = 1:size(A{FileAve1},2);
%     %if  length
%      %beak  
%     %else
% E = [A{FileAve1}(:,i) A{FileAve2}(:,i)];
% F = transpose(mean(E.'));
% StInd = transpose(std(E.'))
% AverInd(:,i) = F(:)
% end



%%% to delete a column of a cell array A{12}(:,5) = []
