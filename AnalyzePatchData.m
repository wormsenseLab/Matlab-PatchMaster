% AnalyzePatchData.m

%% Import Data
ephysData = ImportPatchData();

%% Filter data by project
% Specify project name(s) to keep in cell array of strings. Toss out other 
% crap recordings or unidentified cells that will not be analyzed.
projects = {'FAT'};
dataFields = fieldnames(ephysData);
validProject = zeros(length(dataFields),1);

for iProj = 1:length(projects)
    isValid = strncmpi(dataFields,projects{iProj},length(projects{iProj}));
    validProject = validProject + isValid;
end

ephysData = rmfield(ephysData,dataFields(~logical(validProject)));

clearvars -except ephysData

%% Analyze capacity transient
allCells = fieldnames(ephysData);

for iCell = 1:length(allCells)
cellName = allCells{iCell}; %split into project name and cell numbers when feeding input

% TODO: Change pgf names in Patchmaster to reflect OC vs. WC.
protName = 'ct_neg';
% two alternatives for finding instances of the desired protocol
% find(~cellfun('isempty',strfind(ephysData.(cellName).protocols,'ct_neg')));
protLoc = find(strncmp('ct_neg',ephysData.(cellName).protocols,6));

Ct = zeros(1,length(protLoc));
tau = zeros(1,length(protLoc));
Rs = zeros(1,length(protLoc));

for i = 1:length(protLoc)
%     figure(); hold on;
    
    % Pull out capacity transient data, subtract leak at holding, multiply
    % the negative by -1 to overlay it on the positive, then plot, and 
    % combine the two for the mean (visual check yourself if they really 
    % are equivalent)
    ctNeg = -1.*ephysData.(cellName).data{1,protLoc(i)};
    ctPos = ephysData.(cellName).data{1,protLoc(i)+1};

    ctNeg = bsxfun(@minus, ctNeg, mean(ctNeg(1:20,:)));
    ctPos = bsxfun(@minus, ctPos, mean(ctPos(1:20,:)));

    meanCt = mean([ctNeg ctPos],2);
%     plot(mean(ctNeg,2),'b')
%     plot(mean(ctPos,2),'r')
%     plot(meanCt,'k')
    
    deltaV = 10E-3; % V, 10 mV step
    
    % current during last 10 points of voltage step, dependent on series
    % resistance + leak
    IRsLeak = mean(meanCt(140:149));
    
    % current due to Rs + leak resistance for given deltaV, will be 
    % subtracted from total current to get capacitance current
    rsLeakSub = deltaV/IRsLeak;   % estimate of series resistance, for comparison
    ICt = meanCt-IRsLeak;
%     plot(ICt,'g')

    % TODO: IMPORT SAMPLING FREQUENCY from metadata tree
    % Find the peak capacitative current, then find the point where it's
    % closest to zero, and calculate the area in between.
    % TODO: Change intStart to first zero crossing following t=50
    intStart = 52;
    intEnd = find(abs(ICt) == min(abs(ICt(intStart:150))));
    % trapz uses the trapezoidal method to integrate & calculate area under
    % the curve. But it assumes unit spacing, so divide by the sampling
    % frequency (10kHz in this case) to get units of seconds.
    intICt = trapz(ICt(intStart:intEnd))/10000;  
    % Calculate capacitance based on I = C*(dV/dt)
    Ct(i) = intICt/deltaV;

    % For fitting the curve to find the decay time constant, use peak cap 
    % current as the start point, and fit the next 2ms. This finds an
    % equation of the form Y=a*e^(bx), which is V(t) = Vmax*e^(-t/tau). So
    % the time constant tau is -1/capFit.b .
  
    % TODO: Try picking the end time as the time when it decays to 1/2e, to
    % make sure you have enough points to fit the fast component if it does
    % turn out to be a double exponential. Else, fit a shorter time than
    % 5ms. Either way, stick with exp1 because it's simpler (and because
    % you're focusing on the fast component). Compare the two.
    
    sampFreq = 1E4; % Hz
    fitStart = find(ICt == max(ICt(45:60)));  
    [~,fitInd] = min(abs(ICt(fitStart:150)-(ICt(fitStart)/(2*exp(1)))));
    
    fitTime = fitInd/sampFreq; % seconds
    t = 0:1/sampFreq:fitTime; % UPDATE with sampling frequency
    
    capFit = fit(t',ICt(fitStart:fitStart+fitInd),'exp1');
%     plot(capFit,t,ICt(intStart:intStart+minInd));
    
    % Calculate time constant in seconds.
    tau(i) = -1/capFit.b;
    % Calculate series resistance from tau = Rs*Cm, and divide by 1E6 for
    % units of megaohms.
    Rs(i) = tau(i)/Ct(i)/1E6;
end

ephysData.(cellName).Ct = Ct;
ephysData.(cellName).tau = tau;
ephysData.(cellName).Rs = Rs;

end

clearvars -except ephysData;

%%
for iCell = 1:length(allCells)
ephysData.(allCells{iCell}).Rs
end