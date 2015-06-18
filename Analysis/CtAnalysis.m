% CtAnalysis.m
%
% CtAnalysis uses the capacity transients in response to small voltage
% steps (ct_neg and ct_pos pgfs in Patchmaster) to calculate the
% capacitance, series resistance, and time constant of the recording.
% These values are output directly into the ephysData struct array, in the
% struct for the relevant recordings.
%
% EXAMPLE:
%   ephysData = CtAnalysis(ephysData)
% 
% 
% OUTPUT:
%   ephysData       struct        Data is output as a nested struct. Each
%                                 group/recording is a struct inside this.
%                                 Ct, Rs, and tau fields are added to each
%                                 group.
% 
% Created by Sammy Katta, 15-May-2015.

function ephysData = CtAnalysis(ephysData)

allCells = fieldnames(ephysData);

for iCell = 1:length(allCells)
    cellName = allCells{iCell}; %split into project name and cell numbers when feeding input
    
    % UPDATE: after June/July 2014 (FAT025), Patchmaster has separate pgfs for
    % 'OC_ct_neg' and 'WC_ct_neg' to make it easier to go back and figure out
    % which series resistance/capacitance measurements were important.
    
    % Look for pgf names ending in "ct_neg" and note their locations.
    protName = 'ct_neg';
    flippedProts = cellfun(@fliplr, ephysData.(cellName).protocols, ...
        'UniformOutput', false);
    protLoc = find(strncmp(fliplr(protName),flippedProts,length(protName)));
    
    
    C = zeros(1,length(protLoc));
    tau = zeros(1,length(protLoc));
    Rs = zeros(1,length(protLoc));
    
    for i = 1:length(protLoc)
        
        % Pull out capacity transient data, subtract leak at holding, multiply
        % the negative by -1 to overlay it on the positive, then plot, and
        % combine the two for the mean (visual check yourself if they really
        % are equivalent)
        ctNeg = -1.*ephysData.(cellName).data{1,protLoc(i)};
        ctPos = ephysData.(cellName).data{1,protLoc(i)+1};
        
        ctNeg = bsxfun(@minus, ctNeg, mean(ctNeg(1:20,:)));
        ctPos = bsxfun(@minus, ctPos, mean(ctPos(1:20,:)));
        
        meanCt = mean([ctNeg ctPos],2);
        %     figure(); hold on;
        %     plot(mean(ctNeg,2),'b');
        %     plot(mean(ctPos,2),'r');
        %     plot(meanCt,'k');
        
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
        % TODO: Softcode using V to find intStart = stepStart+x. Can't use
        % zero crossing or max b/c will contribute to C calculation error
        intStart = 52;
        intEnd = find(abs(ICt) == min(abs(ICt(intStart:150))));
        % trapz uses the trapezoidal method to integrate & calculate area under
        % the curve. But it assumes unit spacing, so divide by the sampling
        % frequency (10kHz in this case) to get units of seconds.
        intICt = trapz(ICt(intStart:intEnd))/10000;
        % Calculate capacitance based on I = C*(dV/dt)
        C(i) = intICt/deltaV;
        
        % For fitting the curve to find the decay time constant, use peak cap
        % current as the start point, and fit the next 2ms. This finds an
        % equation of the form Y=a*e^(bx), which is V(t) = Vmax*e^(-t/tau). So
        % the time constant tau is -1/capFit.b .
        
        % TODO: Try picking the end time as the time when it decays to 1/2e, to
        % make sure you have enough points to fit the fast component if it does
        % turn out to be a double exponential. Else, fit a shorter time than
        % 5ms. Either way, stick with exp1 because it's simpler (and because
        % you're focusing on the fast component) vs. exp2. Compare the two.
        
        sampFreq = 10000; % Hz
        fitStart = find(ICt == max(ICt(45:60)));
        if max(ICt(45:60)) ~= 0 % avoids crash when recording was lost and all values were 0
            [~,fitInd] = min(abs(ICt(fitStart:fitStart+30)-(ICt(fitStart)/(2*exp(1)))));
            
            fitTime = fitInd/sampFreq; % seconds
            t = 0:1/sampFreq:fitTime;
            
            capFit = fit(t',ICt(fitStart:fitStart+fitInd),'exp1');
            %     plot(capFit,t,ICt(intStart:intStart+minInd));
            
            % Calculate time constant in seconds for calculation
            tau(i) = -1/capFit.b;
            % Calculate series resistance from tau = Rs*Cm, and divide by 1E6 for
            % units of megaohms.
            Rs(i) = tau(i)/C(i)/1E6;
        end
        %         ephysData.(cellName).fitStart(i) = fitStart;
    end
    
    ephysData.(cellName).C = C;
    ephysData.(cellName).tau = tau/1E-3; % record tau in milliseconds
    ephysData.(cellName).Rs = Rs;
    
end

end