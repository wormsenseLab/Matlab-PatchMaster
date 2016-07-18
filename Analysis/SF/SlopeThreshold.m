function [SlopeActu,MaxZeroSlopeActu,StdZeroSlopeActu,MaxZeroActu,StdZeroActu,MaxActuSensorPlateau,StdActuSensorPlateau,CellMaxActuFirst] = SlopeThreshold(ActuSensor) 
 SlopeActu = [];
 
    ActuSensorRunAvg = tsmovavg(ActuSensor,'s',5,1); %running avergage over 5 points overActuSensor SIgnal
    
   for j = 1:size(ActuSensor,2)
    for i = 1:length(ActuSensor)-1
        SlopeActu(i,j) = ActuSensorRunAvg(i+1,j) - ActuSensorRunAvg(i,j);
    end
   end
   
    AvgOnsetSlopeActu = [];MaxZeroSlopeActu =[];StdZeroSlopeActu=[];
    AvgOnsetActu = [];MaxZeroActu =[];StdZeroActu=[];    
        
for i = 1:size(ActuSensor,2);
   AvgOnsetSlopeActu(i)=  mean(SlopeActu(5:100,i)); % B = ActuSensor
   MaxZeroSlopeActu(i)= max(SlopeActu(5:100,i));
   StdZeroSlopeActu(i)= std(SlopeActu(5:100,i));
   AvgOnsetActu(i)=  mean(ActuSensor(1:100,i));
   MaxZeroActu(i)= max(ActuSensor(1:100,i)); 
   StdZeroActu(i)= std(ActuSensor(1:100,i));
   MaxActuSensorOn(i) = max(ActuSensor(:,i));
   CellMaxActuFirst(i) = find([ActuSensor(:,i)] == MaxActuSensorOn(i),1,'first'); 
end
ActuSensorPlateauAvg = []; AvgActuSensorPlateau = []; MaxActuSensorPlateau = [];StdActuSensorPlateau=[];

for i = 1:size(ActuSensor,2);
    ActuSensorPlateauAvg(:,i) = tsmovavg(ActuSensor(CellMaxActuFirst(i):CellMaxActuFirst(i)+500,i),'s',5,1);
 AvgActuSensorPlateau(i)=  mean(ActuSensorPlateauAvg(5:500,i)); % B = ActuSensor
  MaxActuSensorPlateau(i)= max(ActuSensorPlateauAvg(5:500,i));
 StdActuSensorPlateau(i)= std(ActuSensorPlateauAvg(5:500,i));
end
end
  

   
   