%%% Sylvia Fechner
%%% Stanford University, 
%%% 20151115
%%% Update: 20160210
%%% Script to analyze data from FALCON in Displacement Clamp
%%% Five sweeps of steps within one series
%%% ToDo 

%%% does not work for single FiveStep currenlty; problem in replacing NaN
%%% values
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

%% StepAnaly13sis 
close all % closes all figures

%%% hardcoding part:
%%% change:  STF0XX, sampling fequency not yet fully automatic %%%%%% 

name = 'STF019'; % name of recording. placed into varaibel fiels names%
stimuli = 'Fivestep';
Filenumber = 1; % wil be used to extract sampling freuqnecy; first file loaded, maybe change

Files = 1:length(ephysData.(name).protocols);% load all protocols  

% load all data from all protocols 
% load Current, Actuator Sensor, And Cantilver Signal of each step%
% if else statement included to also load protocols which don't have ForceClamp data

A=[];B=[]; C=[];   %1% to get an empty array, if I analyzed different files before
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

% find all files with certain protocol name: currently: Fivestep
% ifelse statement: if not FiveStep, then empty array; 

for i = Files(:,1):Files(:,end);
   if find(strcmpi(ephysData.(name).protocols{1,i}, stimuli)) == 1;
        continue
   else 
         A{1, i} = [];
          B{1, i} = [];
           C{1, i} = [];
            D{1, i} = [];
    end      
end


%replacing "broke" protocols with empty arrays (delete protocols with less then 5 stimuli) 
for i = 1:length(A);  
   if size(A{1, i},2) < 5 == 1  %% if less the five stimuli are within a protocol, the array is replaced by empty columns. this assumes that this happens only if I broe the protocol, because I forgot to download the wavetable in labview
   % include message, which stimuli is broken
       A{1, i} = [];
       B{1, i} = [];
       C{1, i} = [];
       D{1, i} = [];
   else
       continue
   end
end


% deleting whole blocks of FiveStep
%
AllStimuliBlocks = (find(strcmpi(ephysData.(name).protocols, stimuli)))
% it also displays the broken stimuli. change that it is useful
while 1
%FirstValue = [];
prompt = {'BlockNr (leave empty to quit)'};
dlg_title = 'Delete a block?';
num_lines = 1;
defaultans = {''};
IndValues = inputdlg(prompt,dlg_title,num_lines,defaultans);

FirstValue = str2num(IndValues{1});%
% SecondValue = str2num(IndValues{2})
% ThirdValue = str2num(IndValues{3})
%a = IndValues{2};
%b = ('done');

 if isempty(FirstValue) == 1 % ColumnsOneIndentation = ASubtract(:,[FirstValue SecondValue ThirdValue]) 
     break
 else
    A{1, FirstValue}  = []; 
    B{1, FirstValue}  = []; 
    C{1, FirstValue}  = []; 
    D{1, FirstValue}  = []; 
 end
end


% removes all empty cells from the cell array
AShort = A(~cellfun('isempty',A));
BShort = B(~cellfun('isempty',B));
CShort = C(~cellfun('isempty',C));
DShort = D(~cellfun('isempty',D));


%concatenating all stimuli
Aall = [];
Aall = cat(2,AShort{:});
Ball = [];
Ball = cat(2,BShort{:});
Call = [];
Call = cat(2,CShort{:});
Dall = [];
Dall = cat(2,DShort{:});

%%%%% calculate sampling frequency
fs = ephysData.(name).samplingFreq{1, Files(:,Filenumber)}; % maybe change! %% sampling frequency from first file loaded; I currently assume it will be the same
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

ASubtractppA = ASubtract*10^12; % current in pA to visualize easier in subplot

%%% getting max current (min value) for On response (needs to be done differently)
% AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
% for i = 1:size(Ball,2);
%    AverOnset(i)=  mean(Ball(1:100,i));
%     MaxZeroActu(i)= max(Ball(1:100,i));
%    StdZeroActu(i)= std(Ball(1:100,i));
% end

AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
for i = 1:size(Ball,2);
   AverOnset(i)=  mean(Ball(1:100,i));
    MaxZeroActu(i)= max(Ball(1:100,i));
   StdZeroActu(i)= std(Ball(1:100,i));
end


StartBase = MaxZeroActu + 2*StdZeroActu; % set a threshold to find the onset of the stimulus

MinA = [];
CellMin = [];
AverageMaxCurrent = [];
AverageMaxCurrentMinus = [];
Start=[];
% %%
% for i = 1:size(Aall,2);
% Start(i) = find([Ball(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
% Ende(i) = Start(i) + 100; %% could change 100 to be dependent on fs 
% MinA(i) = min(ASubtract(Start(i):Ende(i),i));
% CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first'); % find cell with min value
% Values = ASubtract(:,i);
% AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
% AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % 
% end

 for i = 1:size(Aall,2);
 Start(i) = find([Ball(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
 Ende(i) = Start(i) + 100; %% could change 100 to be dependent on fs 
 MinA(i) = min(ASubtract(Start(i):Ende(i),i));
 CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first'); % find cell with min value
 Values = ASubtract(:,i);
 AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
 AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
 end

% or make subplots for all Fivestep blocks
 %figure() % plot all Current leak Subtracted signals in one plot
%plot(Time,ASubtract)
%plot(Time, AShort{:,1})
 %title('Current')


% SetPoint
ActuSetPoint = []; SetPointDispl = [];
for i = 1:size(Dall,2),
ActuSetPoint(:,i) = Dall(:,i)*1.5;
SetPointDispl(i) = max(ActuSetPoint(:,i));
end

% to get Displacment of Actuator: multiply actuator sensor signal times 
% sensitivity of actuator: 1.5 
ActuSensor = []; AmplitudeDispl = [];
for i = 1:size(Ball,2),
ActuSensor(:,i) = Ball(:,i)*1.5;
AmplitudeDispl(i) = max(ActuSensor(:,i));
end

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
CantiDefl(:,i) = CantiZero(:,i)*Sensitivity;
end

% calculate Indentation = Actuator Sensor - Cantilever Deflection
MeanIndentation = [];
Indentation = [];
for i = 1:size(Call,2);
Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(1000:2000,i));
end


%%% Calculating Force Signal: Cantilever Deflection * Stiffness 
indStiffness = find(strcmpi(headers, 'Stiffness (N/m)'));
Stiffness = raw(Stiffrow,indStiffness); 
Stiffness = cell2mat(Stiffness);

Force = [];
MeanForce = [];
for i = 1:size(Aall,2);
Force(:,i) = CantiDefl(:,i) *Stiffness;
MeanForce(i) = mean(Force(1000:2000,i));
end


%%% write to excel sheet% this can be deleted
%Amplitude = [];
%Amplitude =[AmplitudeDispl',MeanIndentation',AverageMaxCurrent', AverageMaxCurrentMinus'];%, (Files(:,1):Files(:,end))'];%; MeanIndentation(i)'];%, MeanIndentation(i), AverageMaxCurrent(i), AverageMaxCurrentMinus(i), i];


% Calculating Rise time and Overshoot on Cantilever Deflection signals
% shortened to the Onset of the step
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

xScatter = (1:length(MeanIndentation));
figure()
subplot(2,2,1)
scatter(xScatter, LeakA) 
title('control: Leak Current')
ylabel('Current')
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
scatter(MeanIndentation(i:i+4), AverageMaxCurrentMinus(i:i+4),'LineWidth',2)%,'filled') %% would be nice to see the change in leak
%set(h, 'SizeData', markerWidth^2)
hold on
i = i+5;
title('Mean Cur vs Ind')
xlim([0 max(MeanIndentation)+1])
ylabel('Current')
xlabel('Indentation')
hold on 
for j = 1:size(Aall,2)/5
%legend(Files(j))  % include legend again
end
end
hold on 
subplot(2,2,3)
scatter(xScatter, Start) 

hold on 
subplot(2,2,4)
plot(Time,ASubtractppA)
%xlim([0 max(MeanIndentation)+1])
ylabel('Current (pA)')
xlabel('Time (s)')
title('Current')

%%% plotting ForceClamp signals in a subplot
figure()
subplot(3,3,1)
plot(Time,CantiDefl)
xlim([0 0.6])
title('Cantilever Deflection')
ylabel('Deflection (µm)')
hold on
subplot(3,3,2)
plot(Time,Indentation)
xlim([0 0.6])
ylabel('Indentation (µm)')
title('Indentation')
hold on
subplot(3,3,3)
plot(Time,Force)
xlim([0 0.6])
ylabel('Force (µN)')
title('Force')
hold on
subplot(3,3,4)
plot(Time,ActuSetPoint)
xlim([0 0.6])
ylabel('Displacement (µm)')
title('Displacement ActuSetPoint')
hold on
subplot(3,3,5)
plot(Time,ActuSensor)
xlim([0 0.6])
ylabel('Displacement (µm)')
title('Displacement ActuSensor')
hold on
subplot(3,3,6)
plot(Time,normCantiDefl)
xlim([0 0.6])
%ylabel('Displacement (µm)')
title('Cantilever Defl norm')
hold on
subplot(3,3,7)
scatter(MeanIndentation, allRiseTime)  
xlim([0 0.3])
title('RiseTime (CantiDefl)')
ylabel('Rise Time Tau (s)')
xlabel('Indentation')
xlim([0 max(MeanIndentation)+1])
 hold on
 subplot(3,3,8)
 scatter(MeanIndentation, allOvershoot)  
 %xlim([0 max(MeanIndentation)+1])
 ylabel('% to steady state')
 xlabel('Indentation')
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
ylabel('Force')
xlabel('Indentation')
hold on 
for j = 1:size(Aall,2)/5
%legend(Files(j))  % include legend again
end
end

%%%current with and without leak subtraction in a subplot %%%%
figure()
for i = 1:size(Aall,2)
subplot(size(Aall,2)/5,5,i)
plot(Time,Aall(:,i))
ylim([-5*10^-11 1*10^-11])
hold on
plot(Time,ASubtract(:,i))
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
title(round(MeanIndentation(i),1)) %% 
end

figure()
for i = 1:size(Aall,2)
subplot(size(Aall,2)/5,5,i)
plot(Time,normCantiDefl(:,i))
%ylim([-5*10^-11 1*10^-11])
%RecNum = i; % include number of i within legend or title to easier
%determine the position of the plot
title(round(MeanIndentation(i),1)) %% 
end


msgbox('if you want to delete a whole block, run again');



%% delete single recordings %something is not working
close all
ASubtractNew = ASubtract;
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
  ASubtractNew(:,FirstRec) = nan;
  AverageMaxCurrentMinus(:,FirstRec) = nan;
end
end


% Sort Data 
RoundMeanInd = round(MeanIndentation,1);
[SortInd sorted_index] = sort(RoundMeanInd'); % get index of mean indentations
SortCurrent = AverageMaxCurrentMinus(sorted_index);
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
MeanTraces=[];
MeanSameIndCurrent=[];

%FindSameIndInitial = repmat({NaN},1,8)
FindSameIndInitial = {};

for k = 1:length(MergeInd);
FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
%MeanSameIndCurrent(k) = nanmean(SortCurrent(FindSameInd(:,k))); %average MaxCurrent*-1 with same indentation
%MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
%MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
end
FindSameIndNaN = padcat(FindSameIndInitial{:});
FindSameInd = FindSameIndNaN;

for i = 1:length(MergeInd);% FRow
[r,c] = find(isnan(FindSameInd(:,i))); % change code if only one set of parameters; maybe change length(MergeInd)
while sum(isnan(FindSameInd(:,i)))>0
FindSameInd(r,i) =FindSameInd(r-1,i);
end
end


for k = 1:length(MergeInd);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
MeanSameIndCurrent(k) = nanmean(SortCurrent(FindSameInd(:,k))); %average MaxCurrent*-1 with same indentation
MeanSameIndForce(k) = nanmean(SortForce(FindSameInd(:,k))); %average Force with same Indentation
MeanTraces(k,:) = nanmean(SortASubtract((FindSameInd(1,k)):(FindSameInd(end,k)),:),1); % mean traces in a row vector; problem with mean traces; problem, when inddentation oonly ones
end


MeanTraces = MeanTraces'; % transpose to column vector for export to igor
MeanSameIndCurrent = MeanSameIndCurrent';
MeanSameIndForce = MeanSameIndForce';
MeanTracesppA = MeanTraces*10^12; % to get current in pA
NormMeanCurrent = MeanSameIndCurrent/max(MeanSameIndCurrent); % normalize by fit values

figure()
plot(Time, MeanTracesppA) % plot Current in pA
xlim([0 0.6])

ylabel('Current (pA)')
xlabel('Time')
title((name))

%%
%%% write to excel sheet
MergeIndRow = MergeInd';
%MergeIndRow = num2str(MergeIndRow);%how to write each Indentation as col header???

col_header={name,'MeanInd','MeanCurrent','NormMeanCurrent','MeanForce'};     %Row cell array (for column labels)
%row_header(1:10,1)={'Time'};     %Column cell array (for row labels)
xlswrite(name,MergeInd,'Sheet1','B2');     %Write data
xlswrite(name,MeanSameIndCurrent,'Sheet1','C2');     %Write data
xlswrite(name,NormMeanCurrent,'Sheet1','D2');     %Write data
xlswrite(name,MeanSameIndForce,'Sheet1','E2');     %Write data
xlswrite(name,col_header,'Sheet1','A1');     %Write column header
%xlswrite('My_file.xls',row_header,'Sheet1','A2');      %Write row header
col_header2={name,'Ind1','Ind2','Ind3','Ind4','Ind5','Ind6','Ind7','Ind8','Ind9','Ind10'}; %ToDO - get the values for the
%Indentations
xlswrite(name,MeanTraces,'Sheet2', 'B2');   
xlswrite(name,col_header2,'Sheet2','B1'); 

%%% write Matlabvariables
save(sprintf('%s.mat',name)) %save(sprintf('%sTEST.mat',name))

%%% write as csv, because cannot write with mac to excel
%testlabel = sprintf('%s-MergeInd',name)
filename = sprintf('%s.csv',name) ;
fid = fopen(filename, 'w');
% how to include the Filenumber?
fprintf(fid, 'MeanInd, MeanCurrent, NormMeanCurrent, MeanForce \n') %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid)
ExportData = [MergeInd,MeanSameIndCurrent,NormMeanCurrent,MeanSameIndForce];
dlmwrite(filename, ExportData, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('%sTraces.csv',name) ;
fid = fopen(filename, 'w');
%dlmwrite(filename,MergeIndRow,'-append', 'precision', '%.6f','\t')
% how to include the Filenumber?
fprintf(fid,'Ind1, Ind2, Ind3, Ind4, Ind5, Ind6, Ind7, Ind8, Ind9,Ind10 \n') %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid)
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
