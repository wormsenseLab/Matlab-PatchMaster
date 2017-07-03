close all; clear all; clc
%Average voltage gated currents, which were corrected for Rs; bin voltages
%in xx mV

% load notes to get several values automatically needed for the conversion of the signals
loadFileMode = 0; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('AVERAGED-IVStepsSTF.xlsx'); % be careful in which Folder saved.
end

%%
name = 'GN701-'; %important '-' needs to be there
x=5; %averag value mV %check again and explain better

nameVcorrected= 'VoltageCor';
nameIVValues='IVAvgCurOn';

protNameVcor = strcat(name,nameVcorrected);
protNameIVValues = strcat(name,nameIVValues);

protLocVcor = find(strncmp(protNameVcor,text,length(protNameVcor)));
protLocVIV = find(strncmp(protNameIVValues,text,length(protNameIVValues)));

Vcorrected = []; headersVcorrected =[]; IVValues= []; headersIV=[];

for i=1:length(protLocVcor);
   Vcorrected(:,i) = numbers(:,protLocVcor(i));
   headersVcorrected{i} = text(:,protLocVcor(i));
   IVValues(:,i) = numbers(:,protLocVIV(i));
   headersIV{i} = text(:,protLocVIV(i));
end


 Vcorrected = vertcat(Vcorrected(:));
IVValues = vertcat(IVValues(:));

VcorrectedTrans =[];IVValuesTrans =[];
 VcorrectedTrans =  Vcorrected';
 IVValuesTrans = IVValues';
 
 
[SortVoltage sorted_indexV] = sort(VcorrectedTrans); 
SortIVValues = IVValuesTrans(sorted_indexV);
 
SortVoltage  = SortVoltage';
SortIVValues = SortIVValues';


MergeVoltage = builtin('_mergesimpts',SortVoltage,x,'average');
MergeVoltage(any(isnan(MergeVoltage),2),:)=[]; 


tolerance = 2;

k =[];
[~,FRow] = mode(SortVoltage); %gets the Frequency of the most frequent value
FindSameInd= NaN(FRow,length(MergeVoltage));

FindSameIndInitial = {};
for k = 1:length(MergeVoltage);
FindSameIndInitial{k} = find([SortVoltage] >MergeVoltage(k)-tolerance & [SortVoltage]<MergeVoltage(k)+tolerance);
end

FindSameIndNaN = padcat(FindSameIndInitial{:});
FindSameInd = FindSameIndNaN;


for i = 1:length(MergeVoltage);
[r,c] = find(isnan(FindSameInd(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameInd(:,i)))>0
FindSameInd(r,i) =FindSameInd(r-1,i); %r-1 %replaces the nan value with the previous number. maybe better: deleting nan. ichecked,  i does not matter, check again!
end
end

FinalMeanVoltage = []; FinalSTDVoltage = [];FinalMeanIVValues=[];FinalSTDIVValues=[];
for k = 1:length(MergeVoltage);
FinalMeanVoltage(k) = mean(SortVoltage(FindSameInd(:,k))); 
FinalSTDVoltage(k) = nanstd(SortVoltage(FindSameInd(:,k))); 
FinalMeanIVValues(k) = mean(SortIVValues(FindSameInd(:,k))); 
FinalSTDIVValues(k) = nanstd(SortIVValues(FindSameInd(:,k))); 
end

FinalMeanVoltage = FinalMeanVoltage';
FinalSTDVoltage = FinalSTDVoltage';
FinalMeanIVValues = FinalMeanIVValues';
FinalSTDIVValues = FinalSTDIVValues';

%Find number of Averages per Indentation or Force 
ind  = []; NumberOfAvergagesPerInd  = []; LogicOfIndentations =[];
ind = find(isnan(FindSameIndNaN));
FindSameIndNaN(ind)=0;   % replaces all nan values in matrix with 0
LogicOfIndentations =  FindSameIndNaN > 0; % find logials all greater 0
NumberOfAvergagesPerInd = sum(LogicOfIndentations); % sum up the logical to deterine the number of averages
NumberOfAvergagesPerInd = NumberOfAvergagesPerInd';

%Export signals to csv.

ExportAVGICSteps = [FinalMeanVoltage,FinalSTDVoltage,FinalMeanIVValues,FinalSTDIVValues, NumberOfAvergagesPerInd];

%%% write Matlabvariables
save(sprintf('AWG-IVSteps-%s.mat',name)); %save(sprintf('%sTEST.mat',name))

%%% write as csv, because cannot write with mac to excel
filename = sprintf('AWG-IVSteps-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'AWG-Vcor-%s, STD-Vcor-%s, AWG-IVSTeps-%s,STD-IVSteps-%s, NrAVGperIndentation-%s \n',name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportAVGICSteps, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

