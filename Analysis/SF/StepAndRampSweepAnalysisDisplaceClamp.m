%%% Sylvia Fechner
%%% Stanford University, 
%%% 20151115
%%% Update: 20160210
%%% Script to analyze data from FALCON in Displacement Clamp
%%% Five sweeps of steps within one series
%%% To commit to github
%%% go to my branch
%%% git add .
%%% git commit -m 'text'
%%% git push origin my branch
%%% ToDo 

%%% Currently, automatically saved data are saved in the same folder
%%% which AllStimuliBlocks used for analysis 
%%% write into excel sheet --> which values?
    %how to write each Indentation as col header??? in excel and csv
    % mac can't write to excel !!!
    % mean indentation from non-averaged signal?
    % normalized max current non-averaged?
    % write which simuli which series and which sweeps kept
    % n of averaged signal
%%% finding function in igor to load excel or csv files
%%% make it work for ramps as well
%%% running average
%%% ToDo: plot (one) indentation over time


%%% maybe: make subplots for all Fivestep blocks?
%%% find number of nan values to find out the number of averages
%%% include legend again in current vs indentation
%%% include Force vs indentation

%%  load dat.files 
clear all; close all; 
ephysData=ImportPatchData();

%% load notes to get several valus automatically needed for the coversion of the signals
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);

%% StepAnalysis 
close all % closes all figures
clc
%%% hardcoding part:
%%% change:  STF0XX, sampling fequency not yet fully automatic %%%%%% 

name = 'STF009'; % name of recording. placed into varaibel fiels names%
stimuli = 'Step'; %FiveStep or FiveRampHold
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

% find all files with certain protocol name: 
% ifelse statement: if not respective stimuli name (FiveStep or FiveRampHold), then empty array; 
for i = Files(:,1):Files(:,end);
   if find(strcmpi(ephysData.(name).protocols{1,i}, stimuli)) == 1;
        continue
   else 
         A{1, i} = []; B{1, i} = []; C{1, i} = []; D{1, i} = [];          
    end      
end


tf = strcmp('FiveStep',stimuli); %compare Input Stimuli
isStep = strcmp('Step',stimuli);
if tf == 1; %%% if single steps are used, avoid deleting Blocks with less than five
% replacing "broke" protocols with empty arrays (delete protocols with less then 5 stimuli) 
for i = 1:length(A);  
   if size(A{1, i},2) < 5 == 1  %% if less the five stimuli are within a protocol, the array is replaced by empty columns. this assumes that this happens only if I broke the protocol, because I forgot to download the wavetable in labview
       A{1, i} = [];  B{1, i} = []; C{1, i} = []; D{1, i} = [];
   else
       continue
   end
end
else
   disp 'step'
end

% showing in command prompt, which BlockStimuli will be analyzed and which
% are empty, because they have less than five stimuli
% no need to remove the empty array, because it was already deleted in the
% previous step
AllStimuliBlocks = (find(strcmpi(ephysData.(name).protocols, stimuli)))
LessThanFiveStimuli = [];
for i=1:length(AllStimuliBlocks)
LessThanFiveStimuli(i) = isempty(A{1,AllStimuliBlocks(i)});
end
%AllStimuliBlocks
if isStep == 1
    continue
else
LessThanFiveStimuli %dipslay in command window, if FiveSteps or FiveRamps
end


% deleting whole blocks of FiveBlockStimuli; Whole block=Filenumber
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
    A{1, FirstValue}  = [];  B{1, FirstValue}  = []; C{1, FirstValue}  = [];  D{1, FirstValue}  = []; 
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
fs = ephysData.(name).samplingFreq{1, Files(:,AllStimuliBlocks(1))}; % sampling frequency from first StimuliBlock loaded; 
interval = 1/fs;   %%% get time interval to display x axis in seconds TODO: Ramp is 2s long. change in graphs
ENDTime = length(Aall)/fs; %%% don't remember why I complicated it
Time = (0:interval:ENDTime-interval)'; 

%%%%%% Subtract leak current
LeakA = []; ASubtract = [];
 for i = 1:size(Aall,2);
LeakA(i) = mean(Aall(1:100,i));  %%% take the mean of the first 100 Points
ASubtract(:,i) = Aall(:,i) - LeakA(i); %%%
 end

 %%
%%%%%% CurrentSignals %%%%%%%
tf = strcmp('FiveStep',stimuli); %compare Input Stimuli
isStep = strcmp('Step',stimuli);
%%% getting max current (min value) for On response (needs to be done differently)
if tf == 1;
    disp 'FiveStepProtocol'
 AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
for i = 1:size(Ball,2);
   AverOnset(i)=  mean(Ball(1:100,i)); MaxZeroActu(i)= max(Ball(1:100,i)); StdZeroActu(i)= std(Ball(1:100,i));
end
StartBase = MaxZeroActu + 2*StdZeroActu; % set a threshold to find the onset of the stimulus

MinA = []; CellMin = [];AverageMaxCurrent = [];AverageMaxCurrentMinus = [];Start=[];   
 for i = 1:size(Aall,2);
 Start(i) = find([Ball(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
 Ende(i) = Start(i) + (fs/50); %% could change 100 to be dependent on fs 
 MinA(i) = min(ASubtract(Start(i):Ende(i),i));
 CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first'); % find cell with min value
 Values = ASubtract(:,i);
 AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
 AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
 end 
    else
    disp 'RampAndHold'  
    Steigung = [];
    BallAvg = tsmovavg(Ball,'s',5,1); %running avergage over 5 points overActuSensor SIgnal
     AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
   for j = 1:size(Ball,2)
    for i = 1:length(Ball)-1
        Steigung(i,j) = BallAvg(i+1,j) - BallAvg(i,j);
    end
   end  
   
SteigungAbs = [];
SteigungAbs = abs(Steigung);
for i = 1:size(Ball,2);
   AverOnset(i)=  mean(SteigungAbs(5:100,i)); % B = ActuSensor
    MaxZeroActu(i)= max(SteigungAbs(5:100,i));
   StdZeroActu(i)= std(SteigungAbs(5:100,i));
end

StartBase=[];StartBaseOff=[];
StartBase = MaxZeroActu + 4*StdZeroActu; %%% play around with; not perfekt for slower ramps; do simulation

% set a threshold to find the onset of the stimulus, calculated from Sensor
MinA = []; CellMin = [];AverageMaxCurrent = [];AverageMaxCurrentMinus = [];Start=[]; 
StartOffRes = [];EndeOffRes = [];CellMinOff = [];AverageMaxCurrentOff = [];AverageMaxCurrentMinusOff = [];
Velocity = []; VelocityOff = [];
ASubtractAvg = tsmovavg(ASubtract,'s',10,1);%average signal
for i = 1:size(Aall,2);
 Start(i) = find([SteigungAbs(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
 StartOffRes(i) = find([SteigungAbs(:,i)] > StartBase(i),1, 'last');
 Ende(i) = Start(i) + (fs/50); %% ToDO: could change 100 to be dependent on fs 
 EndeOffRes(i) = StartOffRes(i) + (fs/50);
 MinA(i) = min(ASubtractAvg(Start(i):Ende(i),i)); %%% which signal Do I want to take from the averaged one?
 MinAOff(i) = min(ASubtractAvg(StartOffRes(i):EndeOffRes(i),i));
 CellMin(i) = find([ASubtractAvg(:,i)] == MinA(i),1,'first'); % find cell with min value
 CellMinOff(i) = find([ASubtractAvg(:,i)] == MinAOff(i),1,'first');%%%ToDo change it to last? ToDo: calculate off current for Steps
 Values = ASubtractAvg(:,i);
 AverageMaxCurrent(i) = mean(Values(CellMin(i)-10:CellMin(i)+10)); % ToDo: dependent on fs average 11 cells 5+/-min value
 AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
AverageMaxCurrentOff(i) = mean(Values(CellMinOff(i)-10:CellMinOff(i)+10));
AverageMaxCurrentMinusOff(i) =  AverageMaxCurrentOff(i) *-1;
%Velocity(i) = max(SteigungAbs(Start(i):Ende(i),i));
%VelocityOff(i) = max(SteigungAbs(StartOffRes(i):EndeOffRes(i),i));
end 

end


%%

%%%%%% ForceClampSignals %%%%%%%
% SetPoint
ActuSetPoint = []; SetPointDispl = [];
for i = 1:size(Dall,2),
ActuSetPoint(:,i) = Dall(:,i)*1.5;
SetPointDispl(i) = max(ActuSetPoint(:,i));
end

% to get Displacment of Actuator: multiply actuator sensor signal times 
% sensitivity of actuator: 1.5 
ActuSensor = []; AmplitudeDispl = []; MeanDispl =[];
for i = 1:size(Ball,2),
ActuSensor(:,i) = Ball(:,i)*1.5;
AmplitudeDispl(i) = max(ActuSensor(:,i)); % ToDo: is it better to calculate the mean? Am I using this variabel?
%MeanDispl(i) = mean(ActuSensor(EndeOffRes(i)-1000:EndeOffRes(i)-500,i)); %
%ToDo Do I need MeanDispl for RampAndHold?
end


%%% Calculate Velocity For RampAndHoldStimuli
if tf == 1;
    disp 'FiveStepProtocol'
    continue
else
   SlopeSensor = [];
   SlopeSensorAvg = tsmovavg(ActuSensor,'s',5,1); %running avergage over 5 points overActuSensor SIgnal
   for j = 1:size(ActuSensor,2)
    for i = 1:length(ActuSensor)-1
        SlopeSensor(i,j) = (SlopeSensorAvg(i+1,j) - SlopeSensorAvg(i,j))/(Time(i+1) - Time(i));
    end
   end  
   
   SlopeSensorAbs = abs(SlopeSensor);
  for i = 1:size(Aall,2); 
  Velocity(i) = max(SlopeSensorAbs(Start(i):Ende(i),i));
VelocityOff(i) = max(SlopeSensorAbs(StartOffRes(i):EndeOffRes(i),i));
  end
end
%   figure()
%   plot(SlopeSensorAbs)

  
% to get Deflection of Cantilever: multiply with Sensitivity 
% get Sensitivity from Notes Day of Recording  
FindRowStiff = strcmpi(raw,name); % name = recorded cell
[Stiffrow,col] = find(FindRowStiff); % Siffrow: row correasponding to recorded cell

headers = raw(1,:);
ind = find(strcmpi(headers, 'Sensitivity(um/V)')); % find col with Sensitivity
Sensitivity = raw(Stiffrow,ind); 
Sensitivity = cell2mat(Sensitivity);

Baseline =[];CantiZero = [];CantiDefl = [];
for i = 1:size(Call,2);
Baseline(i) = mean(Call(1:100,i));
CantiZero(:,i) = Call(:,i)-Baseline(i);
CantiDefl(:,i) = CantiZero(:,i)*Sensitivity; % CantiDefl is in um
end

% calculate Indentation = Actuator Sensor - Cantilever Deflection (is in um)
if tf == 1;
    disp 'check, if it is correct'
    for i = 1:size(Call,2);
    Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(1000:2000,i));
end
else
MeanIndentation = [];
Indentation = [];
for i = 1:size(Call,2);
Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(EndeOffRes(i)-1000:EndeOffRes(i)-500,i));
end
end




%%% Calculating Force Signal: Cantilever Deflection * Stiffness 
indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
Stiffness = raw(Stiffrow,indStiffness); 
Stiffness = cell2mat(Stiffness);

Force = [];
MeanForce = [];
for i = 1:size(Aall,2);
Force(:,i) = CantiDefl(:,i)*Stiffness; % I kept it in uN, otherwise: CantiDefl(:,i)*10^-6 *Stiffness; 
MeanForce(i) = mean(Force(1000:2000,i));
end


%CantiDefl(:,1)*10^-6
% Calculating Rise time and Overshoot on Cantilever Deflection signals
% shortened to the Onset of the step
% ToDo: needs to be modified for Ramp

 CantiDeflShort = [];
 MeanCantiDefl = [];
 normCantiDefl = [];
  for i = 1:size(CantiDefl,2);
      EndeCanti(i) = Start(i)+1000;
  CantiDeflShort(:,i) = CantiDefl(Start(i):EndeCanti(i),i); 
  MeanCantiDefl(i) =  mean(CantiDefl(1000:2000,i));
  normCantiDefl(:,i) = CantiDefl(:,i)/MeanCantiDefl(i);
  end
 
TimeShort = (0:interval:length(CantiDeflShort)/fs-interval)';  
InfoSignal = stepinfo(CantiDeflShort, TimeShort, MeanCantiDefl, 'RiseTimeLimits', [0.0 0.63]); %%% over sorted data?? 
allRiseTime = cat(1,InfoSignal.RiseTime);
allOvershoot = cat(1,InfoSignal.Overshoot);



%% now figures
close all

% modify for RampAndHold
% calculate in pA
AverageMaxCurrentMinusppA = AverageMaxCurrentMinus*10^12; LeakAppA = LeakA*10^12; ASubtractppA = ASubtract*10^12; % current in pA to visualize easier in subplot
AallppA=Aall*10^12;

xScatter = (1:length(MeanIndentation));
figure()
subplot(2,2,1)
scatter(xScatter, LeakAppA) 
title('control: Leak Current')
ylabel('Current (pA)')
xlabel('number of file (in recorded order)')
% hold on
% subplot(1,2,2)
% scatter(xScatter, Start)  
% ylim([0 1000])
% ylabel('postion of cell')
% xlabel('number of file (in recorded order)')
% title('control: to see if thresholding is working')
%if stimuli = 'Six-Ramp-Hold' == 1
hold on 
subplot(2,2,2)
i = 1;
while i <= length(MeanIndentation)
scatter(MeanIndentation(i:i+4), AverageMaxCurrentMinusppA(i:i+4),'LineWidth',2)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Cur vs Ind')
xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Indentation')
hold on 
for j = 1:size(Aall,2)/5
%legend(Files(j))  % include legend again
end
end
hold on 
subplot(2,2,3)
scatter(xScatter, Start) 
ylim([400 1000]) % ToDo: change it for Ramps
ylabel('Point')
xlabel('Filenumber')
title('control: Find OnSet of Stimulus')
hold on 
subplot(2,2,4)
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
xlim([0 0.6])
xlabel('Time (s)')
title('Cantilever Deflection')
ylabel('Deflection (µm)')
hold on
subplot(3,3,2)
plot(Time,Indentation)
xlim([0 0.6])
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
xlim([0 0.6])
xlabel('Time (s)')
ylabel('Displacement (µm)')
title('Displacement ActuSetPoint')
hold on
subplot(3,3,5)
plot(Time,ActuSensor)
xlim([0 0.6])
xlabel('Time (s)')
ylabel('Displacement (µm)')
title('Displacement ActuSensor')
hold on
subplot(3,3,6)
plot(Time,normCantiDefl)
xlim([0 0.6])
xlabel('Time (s)')
ylabel('normalized Deflection')
title('Cantilever Defl norm')
hold on
subplot(3,3,7)
scatter(MeanIndentation, allRiseTimeInms)  
xlim([0 0.3])
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


%%%current with and without leak subtraction in a subplot %%%%


figure()
for i = 1:size(AallppA,2)
subplot(size(AallppA,2)/5,5,i)
plot(Time,AallppA(:,i))
%ylim([-5*10^-11 1*10^-11])
hold on
plot(Time,ASubtractppA(:,i))
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
title(round(MeanIndentation(i),1)) %% 
end
%hold on
suptitle({'Current (pA) with (red) and without (blue) leak subtraction';'Bold numbers: Indentation in µm'}) %('')

figure()
for i = 1:size(Aall,2)
subplot(size(Aall,2)/5,5,i)
plot(Time,normCantiDefl(:,i))
%ylim([-5*10^-11 1*10^-11])
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
title(round(MeanIndentation(i),1)) %% 
end


%msgbox('if you want to delete a whole block, run again');



%% delete single recordings 
%close all
%ToDo: has to be changed for ramps, because I want to average the current with
%same velocity 

ASubtractNew = ASubtract;
AverageMaxCurrentMinusNew = AverageMaxCurrent;
MeanIndentationNew = MeanIndentation;
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
  AverageMaxCurrentMinusNew(:,FirstRec) = NaN;
  MeanIndentationNew(:,FirstRec) = NaN;
  %MeanSameIndCurrent(:,FirstRec) = nan; % how does it work, when I am defining it later?
  %NormMeanCurrent(:,FirstRec) = nan; % has to be deleted; was saved from previous work in workspac
  %MeanSameIndForce(:,FirstRec) = nan; %ToDo: does not work
  %MeanIndentation(:,FirstRec) = nan; how to delete mergeInd?
end
end

%MeanNormSameIndCurrent
% Sort Data 

% delete mergeInd
% int_cols = all(isnan(MeanIndentationNew)|round(MeanIndentationNew)==MeanIndentationNew,1);
% it = MeanIndentationNew(:,int_cols);
% 
% test = [1   NaN   2.2   3.2  4;
%      NaN 7.9   5.1   NaN  5;
%      3    5.5  NaN   4.1  NaN];
% int_cols = all(isnan(MeanIndentationNew)|round(MeanIndentationNew)==MeanIndentationNew,1);
% it = MeanIndentationNew(:,int_cols);
% flt = MeanIndentationNew(:,~int_cols);

%%%%%% Sort Data


MeanTraces=[];
MeanSameIndCurrent=[];
MeanSameIndForce=[];
NumberTracesPerInd=[];


if tf == 1;
    disp 'FiveStepProtocol'
RoundMeanInd = round(MeanIndentation,1);
[SortInd sorted_index] = sort(RoundMeanInd'); % get index of mean indentations
SortCurrent = AverageMaxCurrentMinusNew(sorted_index);
%SortNormCurrent = AverageMaxCurrentMinusNew(sorted_index);
SortForce = MeanForce(sorted_index);
transAsub = ASubtractNew';
SortASubtract = transAsub(sorted_index,:);
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
else %% if only one FiveStepProtcol was applied
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
MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
%MeanNormSameIndCurrent
end

else
    disp 'FiveRampHold'
   % RoundMeanInd = round(Velocity,1);   
[SortVel sorted_index] = sort(Velocity); % get index of mean Velocity
SortCurrent = AverageMaxCurrentMinusNew(sorted_index);
%SortNormCurrent = AverageMaxCurrentMinusNew(sorted_index); % Do I need this; yes, normalized to Off Response
SortForce = MeanForce(sorted_index);
transAsub = ASubtractNew';
SortASubtract = transAsub(sorted_index,:);
%SortASubtract = SortASubtract'; keep it as row, easier to calculate
MergeVel = [];
SortVel = SortVel';
MergeVel = builtin('_mergesimpts',SortVel,20,'average');%'average'); %%% TODO: does not work merged Velocity values with +/- 0.1 distance
tolerance = 20; % tolerance to group velocity
%TODo: find best value for velocity merge

k =[];
[~,FRow] = mode(SortVel); % TODO: does not work for Velocity gets the Frequency of the most frequent value
FindSameInd= NaN(10,length(MergeVel)); % determine row length with highest probability of most frequent value
%ToDo: changeValue: currently hardcoded.
if size(Aall,2) > 5
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
%NumberTracesPerInd = NumberTracesPerInd';

for k = 1:length(MergeVel);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
MeanSameIndCurrent(k) = nanmean(SortCurrent(FindSameInd(:,k))); %average MaxCurrent*-1 with same indentation
MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
%MeanNormSameIndCurrent
end
end

%TODO: MeanTraces: first 30 values NAN; why


NormMeanCurrent=[];
MeanTraces = MeanTraces'; % transpose to column vector for export to igor
MeanSameIndCurrent = MeanSameIndCurrent';
MeanSameIndForce = MeanSameIndForce';
MeanTracesppA = MeanTraces*10^12; % to get current in pA
NormMeanCurrent = MeanSameIndCurrent/max(MeanSameIndCurrent); % normalize by fit values

% include here calculation of 
% FindLogicalNumberOfTraces = isnan(FindSameInd) == 0;
% TracesPerIndentation = sum(FindLogicalNumberOfTraces);
% TracesPerIndentation = TracesPerIndentation';
%missing calculation for numer of traces for only one block

figure()
plot(Time, MeanTracesppA) % plot Current in pA
%xlim([0 0.6])

ylabel('Current (pA)')
xlabel('Time')
title((name))

%%
%%% write to excel sheet
MergeIndRow = MergeInd';
%MergeIndRow = num2str(MergeIndRow);%how to write each Indentation as col header???

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

%%% write Matlabvariables
save(sprintf('%s.mat',name)); %save(sprintf('%sTEST.mat',name))

%%% write as csv, because cannot write with mac to excel
%testlabel = sprintf('%s-MergeInd',name)
filename = sprintf('%s.csv',name) ;
fid = fopen(filename, 'w');
% how to include the Filenumber?
fprintf(fid, 'MeanInd, MeanCurrent, NormMeanCurrent, MeanForce \n'); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
%ExportData = [MergeInd,MeanSameIndCurrent,NormMeanCurrent,MeanSameIndForce,TracesPerIndentation];
ExportData = [MergeInd,MeanSameIndCurrent,NormMeanCurrent,MeanSameIndForce];
dlmwrite(filename, ExportData, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('%sTraces.csv',name) ;
fid = fopen(filename, 'w');
%dlmwrite(filename,MergeIndRow,'-append', 'precision', '%.6f','\t')
% how to include the Filenumber?
fprintf(fid,'Ind1, Ind2, Ind3, Ind4, Ind5, Ind6, Ind7, Ind8, Ind9,Ind10,Ind11,Ind12,Ind13,Ind14 \n'); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, MeanTraces, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.



%%
%find number of nan values to find out the number of averages
%average traces 
%how to average two columns
%  for k = 1:length(MergeInd);
%      for i = 1:length(FindSameInd(:,k));
% % MeanTracesCurrent(:,k) = nanmean(SortASubtract(FindSameInd(:,k)))
% %B(:,nn+1) MeanTracesCurrent = nanmean(SortASubtract(:,1:2))
%      end
%  end



%%
%%%delete up to three row and replacing all values with nan
% 
% ASubtractNew = ASubtract;
% prompt = {'FirstValue','SecondValue','ThirdValue','ForthValue','FifthValue'};
% IndValues = inputdlg(prompt)
% FirstValue = str2num(IndValues{1})
% SecondValue = str2num(IndValues{2})
% ThirdValue = str2num(IndValues{3})
% ForthValue = str2num(IndValues{4})
% FifthValue = str2num(IndValues{5})
% 
% ASubtractNew(:,FirstValue) = nan;
% ASubtractNew(:,SecondValue) = nan;
% ASubtractNew(:,ThirdValue) = nan;
% ASubtractNew(:,ForthValue) = nan;
% ASubtractNew(:,FifthValue) = nan;


% %% Average of variable numbers of columns; not in a for loop yet
% while 1
% prompt = {'FirstValue','SecondValue','ThirdValue','InputReady'};
% dlg_title = 'Delete a recording?';
% num_lines = 1;
% defaultans = {'','','',''};
% IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);
% FirstValue = str2num(IndValues{1})
% SecondValue = str2num(IndValues{2})
% ThirdValue = str2num(IndValues{3})
% a = IndValues{4}
% b = ('done')
% 
% %xdatatemp = ASubtract(:,[1 6])
% if isempty(ThirdValue) == 1
%   ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue])
% else
%    ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue ThirdValue]) 
% end
% %elseif ~strcmp(a,b);
%  %  msgbox('it is not done')
%   if  strcmpi(a,b)
%      msgbox('it is done');
%      break;
%   end
%  end
% 
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
