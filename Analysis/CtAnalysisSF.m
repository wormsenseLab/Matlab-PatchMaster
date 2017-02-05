%Analyze CT_protocols to calculate Cm,tau and Rs


%%  load dat.files 
clear all; close all; clc;
ephysData=ImportPatchData();

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

    name = 'STF031'; % name of recording. placed into varaibel fiels names
    
    % Look for pgf names ending in "WC_ct_neg and WC_ct_pos" and note their locations.
    protNameNeg = 'WC_ct_neg';
    protNamePos = 'WC_ct_pos';
  
    protLocNeg = []; protLocPos = [];
    protLocNeg = find(strncmp(protNameNeg,ephysData.(name).protocols,length(protNameNeg))); 
    protLocPos = find(strncmp(protNamePos,ephysData.(name).protocols,length(protNamePos))); 
    
    C = zeros(1,length(protLocNeg));
    tau = zeros(1,length(protLocNeg));
    Rs = zeros(1,length(protLocNeg));
    
    ctNeg = []; ctPos = [];ctNegLeak = []; ctPosLeak=[];meanCt = [];
   
    
        % Pull out capacity transient data, subtract leak at holding, multiply
        % the negative by -1 to overlay it on the positive, then plot, and
        % combine the two for the mean (visual check yourself if they really
        % are equivalent)
       for i= 1:length(protLocNeg);
           
        ctNeg = -1.*ephysData.(name).data{1,protLocNeg(i)};
        ctPos = ephysData.(name).data{1,protLocPos(i)};
       
        ctNegSubtract = bsxfun(@minus, ctNeg, mean(ctNeg(1:20,:)));
        ctPosSubtract = bsxfun(@minus, ctPos, mean(ctPos(1:20,:)));
        
      meanCt(:,i) = mean([ctNegSubtract ctPosSubtract],2); % gives the leak subtracted mean of each pair of POS-NEG CT Step combined
  
     end


  
 figure(); 
 plot(meanCt)      

% fs = ephysData.(name).samplingFreq{1,protLocNeg(1)}; %sampling frequency of first ProtLocNec File. Assumes that all WC_ct are recorded with the same sampling frequency
% interval = 1/fs;   %%% get time interval to display x axis in seconds 
% ENDTime = length(ctNeg)/fs; %%% don't remember why I complicated it
% Time = (0:interval:ENDTime-interval)';



   
deltaV = 10E-3; % V, 10 mV step

 fs = ephysData.(name).samplingFreq{1,protLocNeg(1)}; %sampling frequency of first ProtLocNec File. Assumes that all WC_ct are recorded with the same sampling frequency
 interval = 1/fs;   %%% get time interval to display x axis in seconds 
 ENDTime = length(ctNeg)/fs; %%% don't remember why I complicated it
 Time = (0:interval:ENDTime-interval)';
 TimeShort =  Time(53:150); 

 
IRsLeak = [];ICt = []; ICtShort = [];RmEstimate=[];ChargeQ=[]; Capacity=[]; 
for i=1:size(meanCt,2);

% current during last 10 points of voltage step, dependent on
        % Series resistance + leak; I offset comes from Membrane resistance
  IRsLeak(i) = mean(meanCt(140:149,i));  % Time= Points*interval; change to make it independent of sampling freuqency
          % current due to Rm + leak resistance for given deltaV, will be
        % subtracted from total current to get capacitance current
 RmEstimate(i) = deltaV/IRsLeak(i);   % estimate of membrane resistance, for comparison
 ICt(:,i) = meanCt(:,i)-IRsLeak(i) ; % capacitivit current Ic
 ICtShort(:,i) = ICt(53:150,i);
 ChargeQ(i) = trapz(ICtShort(:,i)/fs);
 Capacity(i) = ChargeQ(i)/deltaV;
end

ChargeQ = ChargeQ';
Capacity = Capacity';
protLocNeg = protLocNeg';
protLocPos = protLocPos';

MeanCapacity = mean(Capacity);

 figure(); 
 plot(ICtShort)      



    %correct until here!!!!!
%%

%%
FitStart=[];
FitEnde=[];
FitStart=2;
FitEnde=10;

%cftool(TimeShort(FitStart:FitEnde),ICtShort(FitStart:FitEnde,1))

 %for i = 1:size(meanCt,2);
       capFit = fit(TimeShort(FitStart:FitEnde),ICtShort(FitStart:FitEnde,1),'exp1');
            %     plot(capFit,t,ICt(intStart:intStart+minInd));
            
            % Calculate time constant in seconds for calculation
           tau(1) = -1/capFit.b;
 
           %end 
 figure()
plot(TimeShort(FitStart:FitEnde),ICtShort(FitStart:FitEnde,1))
hold on
plot(capFit)
%%
figure()
plot(TimeShort,ICtShort(:,7))


 for i = 1:size(Indentation,2);
cftool(Time(Start(i):Start(i)+LengthRamp(i)),Indentation(Start(i):Start(i)+LengthRamp(i),i))
 end 
 

    

      
  


        

