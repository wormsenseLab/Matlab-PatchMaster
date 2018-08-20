%roundVel.m
%
% Previously local fxn in IdAnalysis for rounding velocities  (for 
% combining same condition)across multiple orders of magnitude.

function trioRound = roundVel(velParam) %local function for rounding velocities
velSign = sign(velParam);

trioRound = round(abs(velParam),0); %round small velocities to nearest 1
trioRound(trioRound<500) = round(trioRound(trioRound<500)/2)*2; %round large to nearest 10
trioRound(trioRound>500 & trioRound<1500) = ...
    round(trioRound(trioRound>500 & trioRound<1500)/5)*5; %round medium to nearest 5
trioRound(trioRound>1500) = round(trioRound(trioRound>1500),-1); %round large to nearest 10

trioRound = trioRound .* velSign;
end