function [Start,ActuSetPointZero,CantiDefl,Indentation,MeanIndentation,Force,MeanForce,normCantiDefl,allRiseTime,allOvershoot] = AnalyzeForceClampBOM(interval,ActuSensor,isFiveStep,isFifteenStep,Ball,Call,Dall,Sensitivity,Stiffness,fs);

 
%%% caclulate the threshold for the onset of the Stimulus    
 ActuSensorAvg = [];
  ActuSensorAvg = tsmovavg(ActuSensor,'s',5,1); 
  
    MaxZeroActu =[]; StdZeroActu=[]; %AvgOnsetActu = [];
for i = 1:size(ActuSensor,2);
   MaxZeroActu(i)= max(ActuSensorAvg(6:0.02*fs,i));
   StdZeroActu(i)= std(ActuSensorAvg(6:0.02*fs,i));
%    MaxActuSensorOn(i) = max(ActuSensor(:,i)); %ToDo: change to  value
%    where max is +/-; need it for OffResponse, is not working right now
%    CellMaxActuFirst(i) = find([ActuSensor(:,i)] == MaxActuSensorOn(i),1,'first'); 
end

StartBase = [];   
StartBase = MaxZeroActu + 3*StdZeroActu; 



% BelowPlateau = []; %gives the end of the stimulus % is not working right now
% BelowPlateau = MaxActuSensorPlateau - 8*StdActuSensorPlateau; 

%%% 
Start=[];  
%Ende = []; EndeOff = []; StartOffBelow = []; EndeStimulus = [];LengthInms =[]; LengthStimulus=[]; 
% StartOffBelowShort(i) = find([ActuSensorAvg(CellMaxActuFirst(i):end,i)] < BelowPlateau(i),1, 'first');% %% find cell, where 1st value is lower than threshold; Onset Off-Stimulus
% StartOffBelow(i) = StartOffBelowShort(i)+ CellMaxActuFirst(i); % CellMaxActuFirst excludes the values before Max Value

for i = 1:size(Ball,2);
    Start(i) = find([ActuSensorAvg(:,i)] > StartBase(i),1, 'first'); %% find cell, where 1st value is bigger than threshold; Onset On-Stimulus 
  %  EndeStimulus(i)= StartOffBelow(i)-(0.300/interval); % 
end

%%% SetPoint, Mean of max Displacement
ActuSetPoint = []; MeanSetPointDispl = []; BaselineSetPoint = [];ActuSetPointZero=[];StdSetPointDispl =[];
for i = 1:size(Dall,2),
    ActuSetPoint(:,i) = Dall(:,i); 
BaselineSetPoint(i) = mean(ActuSetPoint(1:0.02*fs,i));
ActuSetPointZero(:,i) = ActuSetPoint(:,i)-BaselineSetPoint(i);
MeanSetPointDispl(i) = mean(ActuSetPointZero(Start+0.02*fs:Start+0.28/interval,i)); %ToDo: vary StepSize, hardcoded
StdSetPointDispl(i)= std(ActuSetPointZero(Start+0.02*fs:Start+0.28/interval,i));
end


%%% Actuator Sensor, Mean of max Displacement
MeanActuSensorDispl = []; 
for i = 1:size(ActuSensor,2),
%MaxDisplSensor(i) = max(ActuSensor(:,i)); 
MeanActuSensorDispl(i) = mean(ActuSensor(Start+0.02*fs:Start+0.28/interval,i)); 
StdActuSensorDispl(i)= std(ActuSensor(Start+0.02*fs:Start+0.28/interval,i));
end

%%% Cantilever Signal, Mean of max Displacement
CantiAvg = []; Baseline =[]; CantiDeflCor = [];CantiDefl =[];
CantiAvg = tsmovavg(Call,'s',10,1);
for i = 1:size(Call,2);
 CantiDeflCor =  CantiAvg*Sensitivity; %Sensitivity from MetaData Sheet
Baseline(i) = mean(CantiDeflCor (11:0.02*fs,i));
CantiDefl(:,i) = CantiDeflCor(:,i)-Baseline(i); % CantiDefl is in um
end

%%% Indentation & Force
% calculate Indentation = Actuator Sensor - Cantilever Deflection (is in um)
% Calculating Force Signal: Cantilever Deflection * Stiffness
MeanIndentation = []; Indentation = []; Force = []; MeanForce = [];
if isFiveStep == 1 || isFifteenStep == 1;
    for i = 1:size(Call,2);
    Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i); % 
    MeanIndentation(i) = mean(Indentation(Start(i):Start(i)+0.28/interval,i));
    Force(:,i) = CantiDefl(:,i)*Stiffness;% I kept it in uN, otherwise: CantiDefl(:,i)*10^-6 *Stiffness; 
    MeanForce(i) = mean(Force(Start(i):Start(i)+0.28/interval,i));
    end
else
   disp 'has to be re written'
% for i = 1:size(Call,2);
% Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
% MeanIndentation(i) = mean(Indentation(EndeOff(i)-1000:EndeOff(i)-500,i));
%  Force(:,i) = CantiDefl(:,i)*Stiffness;% I kept it in uN, otherwise: CantiDefl(:,i)*10^-6 *Stiffness; 
%  MeanForce(i) = mean(Force(EndeOff(i)-1000:EndeOff(i)-500,i)); %cgange to make it independent of fs
% end
end

%%% normalize CantileverSignal
CantiDeflShort = []; MeanCantiDefl = [];  normCantiDefl = [];
 
  for i = 1:size(CantiDefl,2);
      EndeCanti(i) = Start(i)+0.28/interval;
  CantiDeflShort(:,i) = CantiDefl(Start(i):EndeCanti(i),i); %needed to use Stepinfo function 
  MeanCantiDefl(i) =  mean(CantiDefl(Start(i):Start(i)+0.28/interval,i));  
  normCantiDefl(:,i) = CantiDefl(:,i)/MeanCantiDefl(i);
  end
  
TimeShort = (0:interval:length(CantiDeflShort)/fs-interval)';  
InfoSignal = stepinfo(CantiDeflShort, TimeShort, MeanCantiDefl, 'RiseTimeLimits', [0.0 0.63]); %%% over sorted data?? 
allRiseTime = cat(1,InfoSignal.RiseTime);
allOvershoot = cat(1,InfoSignal.Overshoot);

end
 