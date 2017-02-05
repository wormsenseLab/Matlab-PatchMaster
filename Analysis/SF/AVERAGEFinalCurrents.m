%%% 20170104 
%%% Average Final Recordings 
%%% 

%% Load excel sheet with Averagde data
clear all; close all; clc;
loadFileMode = 1; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('AVERAGED-StepsSTF.xlsx'); % be careful in which Folder saved.
end
%%

name = 'GN653-';

nameInd = 'AVGInd';
nameNormOnCurInd = 'AVGNormCurInd';
nameForce = 'AVGForce';
nameNormOnForce= 'AVGNormCurForce';
nameIndOnCur = 'AVGCur';

protName = strcat(name,nameInd); 
protNameNormOnCur = strcat(name,nameNormOnCurInd);
protNameOnCur = strcat(name,nameIndOnCur);
protNameForce = strcat(name,nameForce);
protNameOnForceNorm = strcat(name,nameNormOnForce);

protLoc = []; protLocNormOnCur = []; protLocOnCur = []; protLocForce = [];
protLocNormOnForce = [];

protLoc = find(strncmp(protName,text,length(protName)));
protLocNormOnCur = find(strncmp(protNameNormOnCur,text,length(protNameNormOnCur))); 
protLocOnCur = find(strncmp(protNameOnCur,text,12)); %length(protName)
protLocForce = find(strncmp(protNameForce,text,12)); %length(protName)
protLocNormOnForce = find(strncmp(protNameOnForceNorm,text,length(protNameOnForceNorm))); 


Indentation = []; NormOnCurInd = []; OnCurrent = [];
Force = []; NormOnCurForce = [];

headersInd = {};
headersNormOnCur = {};
headersOnCur = {};
headersForce = {};
headersNormOnCurForce = {};

for i=1:length(protLoc);
   Indentation(:,i) = numbers(:,protLoc(i));
   headersInd{i} = text(:,protLoc(i));
   NormOnCurInd(:,i) = numbers(:,protLocNormOnCur(i));
   headersNormOnCur  {i} = text(:,protLocNormOnCur(i));
   OnCurrent(:,i) = numbers(:,protLocOnCur(i));
   headersOnCur{i} = text(:,protLocOnCur(i)); 
   Force(:,i) = numbers(:,protLocForce(i));
   headersForce{i} = text(:,protLocForce(i)); 
   NormOnCurForce(:,i) = numbers(:,protLocNormOnForce(i));
   headersNormOnCurForce{i} = text(:,protLocNormOnForce(i)); 
end

OnCurrentForce = [];
Indentation = vertcat(Indentation(:));
NormOnCurInd = vertcat(NormOnCurInd(:));
OnCurrent = vertcat(OnCurrent(:));
OnCurrentForce = OnCurrent;
Force = vertcat(Force(:));
NormOnCurForce = vertcat(NormOnCurForce(:));


SortTransNormOnCurInd = []; IndentationTrans = []; TransOnCurrent =[];TransForce = [];
SortForce = []; SortOnCurrent = []; TransNormOnCurForce = []; SortTransNormOnCurForce =[];SortOnCurrentForce=[];TransOnCurrentForce =[];

IndentationTrans = Indentation';
TransNormOnCurInd = NormOnCurInd';
TransOnCurrent = OnCurrent';
TransNormOnCurForce = NormOnCurForce';
TransForce = Force';
TransOnCurrentForce = OnCurrentForce'; 

[SortInd sorted_index] = sort(Indentation); 
[SortForce sorted_indForce] = sort(Force); 


SortOnCurrent = TransOnCurrent(sorted_index);
SortTransNormOnCurInd  = TransNormOnCurInd(sorted_index);

SortOnCurrentForce = TransOnCurrentForce(sorted_indForce);
SortTransNormOnCurForce = TransNormOnCurForce(sorted_indForce);

%SortInd = SortInd'
SmallIndLogic = SortInd <= 3.0 ;
SmallInd = find(SmallIndLogic);
SmallIndValues = SortInd(SmallInd);

BigIndLogic = SortInd > 3.0; 
BigInd = find(BigIndLogic);
BigIndValues = SortInd(BigInd);


%MergeSmallInd = builtin('_mergesimpts',SmallIndValues,0.2,'average');
%MergeBigInd = builtin('_mergesimpts',BigIndValues,0.5,'average');

MergeInd = builtin('_mergesimpts',SortInd,0.4,'average');
MergeInd(any(isnan(MergeInd),2),:)=[]; %needs to be removed, when big and small are used

display 'problem with merge force; findsameForce & MergeForce not same length'

MergeForce = builtin('_mergesimpts',SortForce,0.37,'average');
MergeForce(any(isnan(MergeForce),2),:)=[]; %needs to be removed, when big and small are used
%MergeInd = [MergeSmallInd;MergeBigInd];


tolerance = 0.2; % tolerance to group indentations
toleranceForce = 0.080; 
%MergeInd(any(isnan(MergeInd),2),:)=[]; %remove NaN values of the matrix.
%don't need it if i do logical indexing befor

k =[];
[~,FRow] = mode(SortInd); %gets the Frequency of the most frequent value
FindSameInd= NaN(FRow,length(MergeInd));

[~,FRowForce] = mode(SortForce);
FindSameForce=NaN(FRowForce,length(MergeForce)); 

FindSameIndInitial = {};
for k = 1:length(MergeInd);
FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
end

FindSameIndInitialForce = {};
for l = 1:length(MergeForce);
FindSameIndInitialForce{l} = find([SortForce] >MergeForce(l)-toleranceForce & [SortForce]<MergeForce(l)+toleranceForce);
end

FindSameIndNaN = padcat(FindSameIndInitial{:});
FindSameInd = FindSameIndNaN;

FindSameIndNaNForce = padcat(FindSameIndInitialForce{:});
FindSameForce = FindSameIndNaNForce;

%FindSameInd(isnan(FindSameInd)) = 1 ;


for i = 1:length(MergeInd);
[r,c] = find(isnan(FindSameInd(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameInd(:,i)))>0
FindSameInd(r,i) =FindSameInd(r-1,i); %r-1 %replaces the nan value with the previous number. maybe better: deleting nan. ichecked,  i does not matter, check again!
end
end

for i = 1:length(MergeForce);
[r,c] = find(isnan(FindSameForce(:,i))); % fails, if only one Block of recording; include it into if statement for this reason
while sum(isnan(FindSameForce(:,i)))>0
FindSameForce(r,i) =FindSameForce(r-1,i); %r-1 %replaces the nan value with the previous number. maybe better: deleting nan. ichecked,  i does not matter, check again!
end
end

LengthMergeInd = length(MergeInd)
LengthFindSameInd = length(FindSameInd)
LengthMergeForce = length(MergeForce)
LengthFindSameForce =length(FindSameForce)

FinalMeanNormCurrent = []; FinalMeanIndentation = []; FinalSTDNormCurrent = []; FinalSTDIndentation = []; FinalMeanOnCurrent = [];FinalSTDOnCurrent =[]; 
FinalMeanForce =[]; FinalSTDForce =[]; FinalMeanNormCurForce = [];FinalSTDNormCurForce =[];
for k = 1:length(MergeInd);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
FinalMeanNormCurrent(k) = nanmean(SortTransNormOnCurInd(FindSameInd(:,k))); %average with same indentation
FinalMeanIndentation(k) = nanmean(SortInd(FindSameInd(:,k))); %average  with same indentation
FinalSTDNormCurrent(k) = nanstd(SortTransNormOnCurInd(FindSameInd(:,k))); %average with same indentation
FinalSTDIndentation(k) = nanstd(SortInd(FindSameInd(:,k))); 
FinalMeanOnCurrent(k) = mean(SortOnCurrent(FindSameInd(:,k))); 
FinalSTDOnCurrent(k) = nanstd(SortOnCurrent(FindSameInd(:,k))); 
FinalMeanForce(k) = mean(SortForce(FindSameInd(:,k))); 
FinalSTDForce(k) = nanstd(SortForce(FindSameInd(:,k))); 
FinalMeanNormCurForce(k) = mean(SortTransNormOnCurForce(FindSameInd(:,k))); 
FinalSTDNormCurForce(k) = nanstd(SortTransNormOnCurForce(FindSameInd(:,k))); 
end


FinalMeanNormCurrent = FinalMeanNormCurrent';
FinalMeanIndentation =FinalMeanIndentation';
FinalSTDNormCurrent = FinalSTDNormCurrent';
FinalSTDIndentation = FinalSTDIndentation';
FinalMeanOnCurrent = FinalMeanOnCurrent';
FinalSTDOnCurrent = FinalSTDOnCurrent';
FinalMeanForce = FinalMeanForce';
FinalSTDForce = FinalSTDForce';
FinalMeanNormCurForce = FinalMeanNormCurForce';
FinalSTDNormCurForce =FinalSTDNormCurForce';


%calculated for sorted Force
FinalMeanNormCurSortForce = [];FinalSTDNormCurSortForce=[];FinalMeanOnCurrentSortForce=[];FinalSTDOnCurrentSortForce=[];FinalMeanForceSortForce=[];FinalSTDForceSortForce=[];
for k = 1:length(MergeForce);
%FindSameIndInitial{k} = find([SortInd] >MergeInd(k)-tolerance & [SortInd]<MergeInd(k)+tolerance);
FinalMeanNormCurSortForce(k) = nanmean(SortTransNormOnCurForce(FindSameForce(:,k))); %average with same indentation
FinalSTDNormCurSortForce(k) = nanstd(SortTransNormOnCurForce(FindSameForce(:,k)));
FinalMeanOnCurrentSortForce(k) = nanmean(SortOnCurrentForce(FindSameForce(:,k)));
FinalSTDOnCurrentSortForce(k) = nanstd(SortOnCurrentForce(FindSameForce(:,k)));
FinalMeanForceSortForce(k) = nanmean(SortForce(FindSameForce(:,k)));
FinalSTDForceSortForce(k) = nanstd(SortForce(FindSameForce(:,k)));
end


FinalMeanNormCurSortForce = FinalMeanNormCurSortForce';
FinalSTDNormCurSortForce = FinalSTDNormCurSortForce';
FinalMeanOnCurrentSortForce = FinalMeanOnCurrentSortForce';
FinalSTDOnCurrentSortForce=FinalSTDOnCurrentSortForce';
FinalMeanForceSortForce=FinalMeanForceSortForce';
FinalSTDForceSortForce=FinalSTDForceSortForce';


%Find number of Averages per Indentation or Force 
ind  = []; NumberOfAvergagesPerInd  = []; LogicOfIndentations =[];
ind = find(isnan(FindSameIndNaN));
FindSameIndNaN(ind)=0;   % replaces all nan values in matrix with 0
LogicOfIndentations =  FindSameIndNaN > 0; % find logials all greater 0
NumberOfAvergagesPerInd = sum(LogicOfIndentations); % sum up the logical to deterine the number of averages
NumberOfAvergagesPerInd = NumberOfAvergagesPerInd';

%Export signals to csv.

ExportAVGSignals = [FinalMeanIndentation,FinalMeanNormCurrent,FinalSTDIndentation,FinalSTDNormCurrent,FinalMeanOnCurrent,FinalMeanOnCurrent,FinalSTDOnCurrent,FinalMeanForce,FinalSTDForce,FinalMeanNormCurForce,FinalSTDNormCurForce, NumberOfAvergagesPerInd];

ExportAVGSignalsForce = [FinalMeanNormCurSortForce,FinalSTDNormCurSortForce,FinalMeanOnCurrentSortForce,FinalSTDOnCurrentSortForce,FinalMeanForceSortForce,FinalSTDForceSortForce];
%%% write Matlabvariables
save(sprintf('AWGsignalStep-%s.mat',name)); %save(sprintf('%sTEST.mat',name))

%%% write as csv, because cannot write with mac to excel

filename = sprintf('AWGsignalStep-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'AWGInd-%s, AWGNormOnCurInd-%s, AWGInd-STD-%s,AWGNormOnCurSTD-%s,ONCurrentMean-%s,ONCurrentMeanCOPY-%s,ONCurrentSTD-%s,AWGForce-%s,ForceSTD-%s,AWGNormOnCurForce-%s,NormOnCurForce-STD-%s, NrAVGperIndentation-%s \n',name,name,name,name,name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportAVGSignals, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('AWGsignalStepForce-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'AWGNormOnCurSortForce-%s, STDNormOnCurSortForce-%s,AWGOnCurSortForce-%s,STDOnCurSortForce-%s,AWGForceSortForce-%s, STDForceSortForce-%s, \n',name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportAVGSignalsForce, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.


