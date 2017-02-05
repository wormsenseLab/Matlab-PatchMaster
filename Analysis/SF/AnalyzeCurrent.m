function [LeakA, ASubtract, AvgMaxCurrent,AvgMaxCurrentMinus,AvgMaxCurrentOff,AvgMaxCurrentMinusOff,Start,StartOffBelow,Ende,EndeOff,ASubtractAvg,LengthRamp,LengthInms,EndeRamp,StartOffBelowShort,ActuSensorAvg] = AnalyzeCurrent(isFiveStep,isStep,isFifteenStep,ActuSensor,StartBase,Aall,fs,SlopeActu,BelowPlateau,CellMaxActuFirst,interval);%tf,isStep, SteigungAbs,

LeakA = []; ASubtract = []; MinA = []; CellMin = [];AvgMaxCurrent = [];AvgMaxCurrentMinus = [];Start=[];  Ende = []; EndeOff = [];  
CellMinOff = [];AvgMaxCurrentOff = [];AvgMaxCurrentMinusOff = []; 
StartOffBelow = [];ActuSensorAvg=[]; EndeRamp = [];LengthInms =[]; 


%%% this could go to the ForceClamp Section  
ActuSensorAvg = tsmovavg(ActuSensor,'s',10,1);
 
for i = 1:size(Aall,2);
StartOffBelowShort(i) = find([ActuSensorAvg(CellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first');% %% find cell, where 1st value is lower than threshold; Onset Off-Stimulus
StartOffBelow(i) = StartOffBelowShort(i)+ CellMaxActuFirst(i); % CellMaxActuFirst excludes the values before Max Value
end
    if isFiveStep == 1 || isStep ==1 || isFifteenStep ==1 ;
        disp 'Remember: only calculates length of Ramp of Step correctly, if Plateau phase is 300 ms'  
            for i = 1:size(Aall,2);
    Start(i) = find([ActuSensor(:,i)] > StartBase(i),1, 'first'); %% find cell, where 1st value is bigger than threshold; Onset On-Stimulus 
    EndeRamp(i)= StartOffBelow(i)-(0.300/interval); % Time= Points*interval
    
            end
else
    disp 'Remember: only calculates length of Ramp correctly (needed for Fitting), if Plateau phase is 500 ms'  
    for i = 1:size(Aall,2);
    Start(i) = find([SlopeActu(:,i)] > StartBase(i),1, 'first');% SlopeSensor
%     StartOffBelowShort(i) = find([ActuSensorAvg(CellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first');% this is the falling slope
%     StartOffBelow(i) = StartOffBelowShort(i)+ CellMaxActuFirst(i); %pretty good, checked values in Igor
    EndeRamp(i)= StartOffBelow(i)-(0.500/interval); % Time= Points*interval
    end
 end
%%%%%%%

%%%%%% Subtract leak current

  for i = 1:size(Aall,2);
 LeakA(i) = mean(Aall(Start(i)-0.02*fs:Start(i),i));  %%%
 ASubtract(:,i) = Aall(:,i) - LeakA(i); %%%
  end
  
ASubtractAvg = tsmovavg(ASubtract,'s',5,1);%average signal
absASubtractAvg = abs(ASubtractAvg);
   
for i = 1:size(Aall,2);
% StartOffBelowShort(i) = find([ActuSensorAvg(CellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first');% this is the falling slope
% StartOffBelow(i) = StartOffBelowShort(i)+ CellMaxActuFirst(i); %pretty good, checked values in Igor
Ende(i) = Start(i) + (fs/20);
EndeOff(i) = StartOffBelow(i) + (fs/20); % currently not in use
MinA(i) = max(absASubtractAvg(Start(i):Ende(i),i));
%min(ASubtractAvg(Start(i):Ende(i),i));
CellMin(i) = find([absASubtractAvg(:,i)] == MinA(i),1,'first');
%find([ASubtractAvg(:,i)] == MinA(i),1,'first');
MinAOff(i) = min(ASubtractAvg(StartOffBelow(i):StartOffBelow(i)+100,i)); 
CellMinOff(i) = find([ASubtractAvg(:,i)] == MinAOff(i),1,'first');  % change that it looks closer to seco
Values(:,i) = ASubtractAvg(:,i);% UNECESSARY? forget, why I did this%%%%
AvgMaxCurrentOff(i) = mean(ASubtractAvg(CellMinOff(i)-1:CellMinOff(i)+5,i)); % ToDo Five or 10???
AvgMaxCurrentMinusOff(i) =  AvgMaxCurrentOff(i) *-1;
AvgMaxCurrent(i) = mean(Values(CellMin(i)-1:CellMin(i)+5,i));%(i) = mean(Values(CellMin(i)-10:CellMin(i)+10));
AvgMaxCurrentMinus(i) = AvgMaxCurrent(i) *-1;
LengthRamp(i)=EndeRamp(i) - Start(i);
LengthInms(i) = LengthRamp(i)*interval;
end
end

%Velocity(i) = max(SteigungAbs(Start(i):Ende(i),i));
%VelocityOff(i) = max(SteigungAbs(StartOffRes(i):EndeOffRes(i),i));
%  AverageMaxCurrent(i) = mean(Values(CellMin(i)-10:CellMin(i)+10)); % ToDo: dependent on fs average 11 cells 5+/-min value
%  AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
% AverageMaxCurrentOff(i) = mean(Values(CellMinOff(i)-10:CellMinOff(i)+10));
% AverageMaxCurrentMinusOff(i) =  AverageMaxCurrentOff(i) *-1;
% %Velocity(i) = max(SteigungAbs(Start(i):Ende(i),i));
%VelocityOff(i) = max(SteigungAbs(StartOffRes(i):EndeOffRes(i),i));
 %  MinA = []; CellMin = [];AverageMaxCurrent = [];AverageMaxCurrentMinus = [];Start=[]; 
% StartOffRes = [];EndeOffRes = [];CellMinOff = [];AverageMaxCurrentOff = [];AverageMaxCurrentMinusOff = [];
% Velocity = []; VelocityOff = [];
 %Ende(i) = Start(i) + (fs/50); %% could change 100 to be dependent on fs 
% MinA(i) = min(ASubtractAvg(Start(i):Ende(i),i)); % find min current on avergaed signal (running average)
%CellMin(i) = find([ASubtractAvg(:,i)] == MinA(i),1,'first'); % find cell with min value
%Values = ASubtractAvg(:,i);
 %AvgMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
 %AveMaxCurrentMinus(i) = AvgMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
 
 
 
% %%
% %%%%%% CurrentSignals %%%%%%%
% tf = strcmp('FiveStep',stimuli); %compare Input Stimuli
% isStep = strcmp('Step',stimuli);
% %%% getting max current (min value) for On response (needs to be done differently)
% if tf == 1;
%     disp 'FiveStepProtocol'
%  AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
% for i = 1:size(Ball,2);
%    AverOnset(i)=  mean(Ball(1:100,i)); MaxZeroActu(i)= max(Ball(1:100,i)); StdZeroActu(i)= std(Ball(1:100,i));
% end
% StartBase = MaxZeroActu + 2*StdZeroActu; % set a threshold to find the onset of the stimulus
% 
% MinA = []; CellMin = [];AverageMaxCurrent = [];AverageMaxCurrentMinus = [];Start=[];   
%  for i = 1:size(Aall,2);
%  Start(i) = find([Ball(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
%  Ende(i) = Start(i) + (fs/50); %% could change 100 to be dependent on fs 
%  MinA(i) = min(ASubtract(Start(i):Ende(i),i));
%  CellMin(i) = find([ASubtract(:,i)] == MinA(i),1,'first'); % find cell with min value
%  Values = ASubtract(:,i);
%  AverageMaxCurrent(i) = mean(Values(CellMin(i)-5:CellMin(i)+5)); % average 11 cells 5+/-min value
%  AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
%  end 
%     else
%     disp 'RampAndHold'  
%     Steigung = [];
%     BallAvg = tsmovavg(Ball,'s',5,1); %running avergage over 5 points overActuSensor SIgnal
%      AverOnset = [];MaxZeroActu =[];StdZeroActu=[];
%    for j = 1:size(Ball,2)
%     for i = 1:length(Ball)-1
%         Steigung(i,j) = BallAvg(i+1,j) - BallAvg(i,j);
%     end
%    end  
%    
% SteigungAbs = [];
% SteigungAbs = abs(Steigung);
% for i = 1:size(Ball,2);
%    AverOnset(i)=  mean(SteigungAbs(5:100,i)); % B = ActuSensor
%     MaxZeroActu(i)= max(SteigungAbs(5:100,i));
%    StdZeroActu(i)= std(SteigungAbs(5:100,i));
% end
% 
% StartBase=[];StartBaseOff=[];
% StartBase = MaxZeroActu + 4*StdZeroActu; %%% play around with; not perfekt for slower ramps; do simulation
% 
% % set a threshold to find the onset of the stimulus, calculated from Sensor
% MinA = []; CellMin = [];AverageMaxCurrent = [];AverageMaxCurrentMinus = [];Start=[]; 
% StartOffRes = [];EndeOffRes = [];CellMinOff = [];AverageMaxCurrentOff = [];AverageMaxCurrentMinusOff = [];
% Velocity = []; VelocityOff = [];
% ASubtractAvg = tsmovavg(ASubtract,'s',10,1);%average signal
% for i = 1:size(Aall,2);
%  Start(i) = find([SteigungAbs(:,i)] > StartBase(i),1, 'first'); %% find cells, where 1st values is bigger than threshold 
%  StartOffRes(i) = find([SteigungAbs(:,i)] > StartBase(i),1, 'last');
%  Ende(i) = Start(i) + (fs/50); %% ToDO: could change 100 to be dependent on fs 
%  EndeOffRes(i) = StartOffRes(i) + (fs/50);
%  MinA(i) = min(ASubtractAvg(Start(i):Ende(i),i)); %%% which signal Do I want to take from the averaged one?
%  MinAOff(i) = min(ASubtractAvg(StartOffRes(i):EndeOffRes(i),i));
%  CellMin(i) = find([ASubtractAvg(:,i)] == MinA(i),1,'first'); % find cell with min value
%  CellMinOff(i) = find([ASubtractAvg(:,i)] == MinAOff(i),1,'first');%%%ToDo change it to last? ToDo: calculate off current for Steps
%  Values = ASubtractAvg(:,i);
%  AverageMaxCurrent(i) = mean(Values(CellMin(i)-10:CellMin(i)+10)); % ToDo: dependent on fs average 11 cells 5+/-min value
%  AverageMaxCurrentMinus(i) =  AverageMaxCurrent(i) *-1; % multiply with -1 to facilitate demontrating of increase in current
% AverageMaxCurrentOff(i) = mean(Values(CellMinOff(i)-10:CellMinOff(i)+10));
% AverageMaxCurrentMinusOff(i) =  AverageMaxCurrentOff(i) *-1;
% %Velocity(i) = max(SteigungAbs(Start(i):Ende(i),i));
% %VelocityOff(i) = max(SteigungAbs(StartOffRes(i):EndeOffRes(i),i));
% end 
% 
% end