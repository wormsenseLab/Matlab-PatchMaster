
%%  load dat.files 
clear all; close all; 
ephysData=ImportPatchData()

%% load notes from day of analysis
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);

%% StepAnalysis %%% devided into seperate section to avoid loading notes, which is the same for all recordings on one day
close all %% closes all figures
%%% hardcoding part:
%%%% change: Files; STF0XX; rawStiff %%%%%%
Files = 16:46;
name = 'STF010';
rawStiff = 2; %% change if two different cantilevers were used
%%%%% load Current, Actuator Sensor, And Cantilver Signal of each step
A=[];B=[]; C=[];   %% to get an empty array, if I analyzed different files before
for i = Files(:,1):Files(:,end);
A(:,i) = ephysData.(name).data{1, i}; %Current
B(:,i) = ephysData.STF010.data{3, i}; % actuator Sensor
C(:,i) = ephysData.STF010.data{4, i}; % Cantilever Signal
end
%%%%% after loading data
fs = 5000; %%currently always 5000 Hz %%% ephysData.STF009.samplingFreq{1, 15};  %%%% get sampling frequency
interval = 1/fs;   %%%%% get time interval to display x axis in seconds
ENDTime = length(A)/5000; %% don't remember why I complicated it
Time = (0:interval:ENDTime-interval)'; 

%%%%%% Subtract leak current
LeakA = [];
ASubtract = [];
for i = Files(:,1):Files(:,end);
LeakA(i) = mean(A(1:100,i));  %%%% take the mean of the first 100 Points
ASubtract(:,i) = A(:,i)-LeakA(i); %%%% subtract leak  --> copy to Igor
end

% %%%%%getting max current (min value) for On response (needs to be done differently)
MinA = [];
CellMin = [];
AverageMaxCurrent = [];
AverageMaxCurrentMinus = [];
for i = Files(:,1):Files(:,end);
 MinA(i) = min(ASubtract(760:860,i));
CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first');
Values = ASubtract(:,i);
AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5));
 AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1;
end
 









figure() % plot all Current leak Subtracted signals in one plot
plot(Time,ASubtract(:,Files))
title('Current')


%% to get Displacment of Actuator: multiply actuator sensor signal times 
%%%% sensitivity of actuator: 1.5 %%%
ActuSensor = []; AmplitudeDispl = [];
for i = Files(:,1):Files(:,end)
ActuSensor(:,i) = B(:,i)*1.5;
AmplitudeDispl(i) = max(ActuSensor(:,i));
%figure()
%figure()
%plot(Time,ActuSensor(:,i));
%itle('Displacement Actuator')
%title('Displacement Actuator')
end

%%%%% to get Deflection of Cantilever: multiply with Sensitivity 
%%%% get Sensitivity from Notes of Recording day %%%%

headers = raw(1,:);
ind = find(strcmpi(headers, 'Sensitivity(um/V)'));
Sensitivity = raw(rawStiff,ind);  %%% rawstiff = if different cantilever were used change in hardcoding part
Sensitivity = cell2mat(Sensitivity);

Baseline =[];CantiZero = [];CantiDefl = [];
for i = Files(:,1):Files(:,end);
Baseline(i) = mean(C(1:100,i));
CantiZero(:,i) = C(:,i)-Baseline(i);
CantiDefl(:,i) = CantiZero(:,i)*Sensitivity;
% figure()
% plot(Time,CantiDefl(:,i));
% title('Cantilever Deflection')
end

figure()  %%% plots all Cantilever deflections in one plot
plot(Time,CantiDefl(:,Files))
title('Cantilever Deflection')

%%%% calculate Indentation = Actuator Sensor - Cantilever Deflection
for i = Files(:,1):Files(:,end);
Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(1000:2000,i));
end

figure() %%% plots all Indentations in one plot
plot(Time,Indentation(:,Files))
title('Indentation')


%%% Calculating Force Signal 
indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
Stiffness = raw(rawStiff,indStiffness); 
Stiffness = cell2mat(Stiffness);

for i = Files(:,1):Files(:,end);
Force(:,i) = CantiDefl(:,i) *Stiffness;
% figure()
% plot(Time,Force);%
% title('Force')
end
%%% 

figure() %%% plots all Force signals in one plot
plot(Time,Force(:,Files))
title('Force')

Amplitude = [];
Amplitude =[AmplitudeDispl',MeanIndentation',AverageMaxCurrent', AverageMaxCurrentMinus', (1:i)'];%; MeanIndentation(i)'];%, MeanIndentation(i), AverageMaxCurrent(i), AverageMaxCurrentMinus(i), i];

%%%%%% now make all plots
% %%%%%%% figures current with and without leak subtraction %%%%
figure()
for i = Files(:,1):Files(:,end)
%plot(Time,ASubtract(:,i));
%title('Current - with leak substraction')
%hold on
%plot(Time,A(:,i));
%subplot(4,2,i)
%figure()
subplot(7,5,i-Files(:,1)+1)
plot(Time,A(:,i))
hold on
plot(Time,ASubtract(:,i))
title(MeanIndentation(i))
end

