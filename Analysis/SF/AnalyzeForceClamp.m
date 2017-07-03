function [ActuSetPoint,CantiDefl,MeanIndentation,Force,MeanForce,Indentation,normCantiDefl,allRiseTime,allOvershoot] = AnalyzeForceClamp(interval,Start,isFiveStep,isStep,isFifteenStep,Dall,Call,EndeOff,ActuSensor,Sensitivity,Stiffness,fs);
%Dall = concatenated SetPointsignal


ActuSetPoint = []; SetPointDispl = [];
for i = 1:size(Dall,2),
ActuSetPoint(:,i) = Dall(:,i)*1.5;%ToDo: extract from MeataDataSheet
SetPointDispl(i) = max(ActuSetPoint(:,i));
end

MaxDisplSensor = []; MeanDispl =[];
for i = 1:size(ActuSensor,2),
MaxDisplSensor(i) = max(ActuSensor(:,i)); % ToDo: is it better to calculate the mean? Am I using this variabel?
%MeanDispl(i) = mean(ActuSensor(EndeOffRes(i)-1000:EndeOffRes(i)-500,i)); %
%ToDo Do I need MeanDispl for RampAndHold?
end

% Baseline =[];CantiZero = [];CantiDefl = [];
% for i = 1:size(Call,2);
% Baseline(i) = mean(Call(1:100,i));
% CantiZero(:,i) = Call(:,i)-Baseline(i);
% CantiDefl(:,i) = CantiZero(:,i)*Sensitivity; % CantiDefl is in um
% end

%copied from BOM
CantiAvg = []; Baseline =[]; CantiDeflCor = [];CantiDefl =[];
CantiAvg = tsmovavg(Call,'s',10,1);
for i = 1:size(Call,2);
 CantiDeflCor =  CantiAvg*Sensitivity; %Sensitivity from MetaData Sheet
Baseline(i) = mean(CantiDeflCor (11:0.02*fs,i));
CantiDefl(:,i) = CantiDeflCor(:,i)-Baseline(i); % CantiDefl is in um
end

% calculate Indentation = Actuator Sensor - Cantilever Deflection (is in um)
% Calculating Force Signal: Cantilever Deflection * Stiffness 

MeanIndentation = []; Indentation = []; Force = []; MeanForce = [];
if isFiveStep == 1 || isStep == 1 || isFifteenStep == 1;
    for i = 1:size(Call,2);
    Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
    MeanIndentation(i) = mean(Indentation(Start(i):Start(i)+0.28/interval,i));
    Force(:,i) = CantiDefl(:,i)*Stiffness;% I kept it in uN, otherwise: CantiDefl(:,i)*10^-6 *Stiffness; 
    MeanForce(i) = mean(Force(Start(i):Start(i)+0.28/interval,i));
    end
else
for i = 1:size(Call,2);
Indentation(:,i) = ActuSensor(:,i) - CantiDefl(:,i);
MeanIndentation(i) = mean(Indentation(EndeOff(i)-1000:EndeOff(i)-500,i));
 Force(:,i) = CantiDefl(:,i)*Stiffness;% I kept it in uN, otherwise: CantiDefl(:,i)*10^-6 *Stiffness; 
 MeanForce(i) = mean(Force(EndeOff(i)-1000:EndeOff(i)-500,i)); %cgange to make it independent of fs
end
end
CantiDeflShort = [];
 MeanCantiDefl = [];
 normCantiDefl = [];
 
  for i = 1:size(CantiDefl,2);
      EndeCanti(i) = Start(i)+1000;
  CantiDeflShort(:,i) = CantiDefl(Start(i):EndeCanti(i),i); 
  MeanCantiDefl(i) =  mean(CantiDefl(Start(i):Start(i)+0.28/interval,i)); 
  %MeanCantiDefl(i) =  mean(CantiDefl(0.2*fs:0.4*fs,i)); 
  normCantiDefl(:,i) = CantiDefl(:,i)/MeanCantiDefl(i);
  end
  
 TimeShort = (0:interval:length(CantiDeflShort)/fs-interval)';  
InfoSignal = stepinfo(CantiDeflShort, TimeShort, MeanCantiDefl, 'RiseTimeLimits', [0.0 0.63]); %%% over sorted data?? 
allRiseTime = cat(1,InfoSignal.RiseTime);
allOvershoot = cat(1,InfoSignal.Overshoot);

end
 