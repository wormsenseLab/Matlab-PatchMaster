% script to analzye voltage gated currents in On-cell and Whole cell mode.

% 1st series number for On and Whole cell IVqSteps need to be assigned in
% Meta data script
% Avereages three on cell recordings & three whole cell IVq
% the subtracts the averagee of the whole cell from On cell to reduce
% capacitance due to pipette 
% normalize current to input capacitance
% correct voltage for series resistance (maybe not necessary, but included if needed)
% if no On-cell current, because recording was immidiately broken in, NaN value
% or 0 @the IVq Meta-data sheet will skip the analysis for this particular
% cell
% same is true, if capaciance and Rs were not yet calculated, place NaN
% value
% at the beginning of for loop i=2:DetermineLengthCol, can be changed,
% if only the second or first half or only one particular cell should be
% analzyed, otherwise, script will analzye all recordings from the assigned
% metadata sheet


%%  load dat.files 
clear all; close all; clc;
ephysData=ImportPatchData();
%%
%clear all; close all; clc;
%load('ephysdata(20170130).mat')
%% 
% load notes to get several values automatically needed for the conversion of the signals
loadFileMode = 1; % change here, if you want to select a file or load always the same
if loadFileMode  == 0; % 
[filename,pathname] = uigetfile('*.*', 'Load file', 'MultiSelect', 'on'); 
[numbers, text, raw] = xlsread([pathname, filename]);
elseif loadFileMode == 1
[numbers, text, raw] = xlsread('Ephys-Meta-Sylvia-VGCNOTAnal.xlsx'); % be careful in which Folder saved.
end

%% Analysis Individual voltage gated current Recordings 
close all; clc

headers = raw(1,:);
CellID = find(strcmpi(headers, 'CellID'));
AllCol = raw(:,CellID);
DetermineLengthCol = length(AllCol);

%%% hardcoding part, if necessary:
%%% change here, if less analyzed
for i=2:DetermineLengthCol;% DetermineLengthCol %shorten here, if you want to analyze less
callCellID = raw(i,CellID);
name = callCellID{1}

%name = 'STF061'; if you want to analyze only one cell, get the col number of CellID and enter into for loop 

%now calculate mean for OC-curren7
FindRowIndCellId = strcmpi(raw,name); % name = recorded cell
[RowCellId,col] = find(FindRowIndCellId,1); % Siffrow: row correasponding to recorded cell
indOC = find(strcmpi(headers, 'OC_ivq Series')); % find col with Sensitivity
OCIVqInd = raw(RowCellId,indOC); 
OCIVqInd = cell2mat(OCIVqInd);

 Rs = []; TermRsI = []; MinusRs =[];
 indRs = find(strcmpi(headers, 'Rs(MOhm)')); % find col with Sensitivity
 Rs = raw(RowCellId,indRs);
 Rs = cell2mat(Rs)


if ischar(OCIVqInd) 
    
    display 'this cell has a NaN value'

elseif ischar(Rs)
    
     display 'no Rs and Cap cal'
    
elseif OCIVqInd ~= 0 


OCIVq_1 = [];   OCIVq_2 = []; OCIVq_3 = []; 
OCIVq_1 = ephysData.(name).data{1, OCIVqInd};
OCIVq_2 = ephysData.(name).data{1, OCIVqInd+1};
OCIVq_3 = ephysData.(name).data{1, OCIVqInd+2};

OCIVq_3D(:,:,1)=OCIVq_1;

% for i = OCIVqInd:OCIVqInd+2;
%    eval(['OCIVq_3D(:,:,index)=OCIVq_' num2str(index)';]);
% end

for index = 1:3;
   eval(['OCIVq_3D(:,:,index)=OCIVq_' num2str(index)';]);
end


%alternative
%test = cat(3,OCIVq_1,OCIVq_2,OCIVq_3)
%meantest = mean(test,3)

AvgOCIVq = mean(OCIVq_3D,3);

%now calculate mean for WC-current

indWC = find(strcmpi(headers, 'WC_ivq Series')); % find col with Sensitivity
WCIVqInd = raw(RowCellId,indWC); 
WCIVqInd = cell2mat(WCIVqInd)

WCIVq_1 = [];   WCIVq_2 = []; WCIVq_3 = []; 
WCIVq_1 = ephysData.(name).data{1, WCIVqInd}
WCIVq_2 = ephysData.(name).data{1, WCIVqInd+1}
WCIVq_3 = ephysData.(name).data{1, WCIVqInd+2}

WCIVq_3D(:,:,1)=WCIVq_1
%index = WCIVqInd:WCIVqInd+2;
for i = 1:3;
   eval(['WCIVq_3D(:,:,i)=WCIVq_' num2str(i)';]);
end


AvgWCIVq = mean(WCIVq_3D,3);

AvgWCcorrected = AvgWCIVq - AvgOCIVq;

IVValues = [];
for k = 1:size(AvgWCcorrected,2);
   IVValues(k) = mean(AvgWCcorrected(200:500,k)); %TODO: not harcoded
end

% figure()
% plot(AvgWCcorrected)
IVValuesTrans = IVValues';

% normalize current to input capacitance
IVValuesNormCapac = [];CapacInd =[];Capacitance =[];CapacitanceFarad=[];
CapacInd = find(strcmpi(headers, 'C_in (pF)'));
Capacitance = raw(RowCellId,CapacInd);
Capacitance = cell2mat(Capacitance);
CapacitanceFarad = Capacitance*10^-12;
IVValuesNormCapac = IVValuesTrans/CapacitanceFarad;


% Rs = []; TermRsI = []; MinusRs =[];
% indRs = find(strcmpi(headers, 'Rs(MOhm)')); % find col with Sensitivity
% Rs = raw(RowCellId,indRs); 
%Rs = cell2mat(Rs)
RsinOHM = Rs*10E6;
%MinusRs = RsinOHM *-1;
TermRsI = RsinOHM * IVValuesTrans;
%TermRsI = (TermRsI)';
Voltage=[-80;-60;-40;-20;0;20;40;60;80];
VoltageInV=[-0.080;-0.060;-0.040;-0.020;0;0.020;0.040;0.060;0.080];
Vcorrected = VoltageInV-TermRsI; %minus(Voltage,TermRsI); % Vcorrected =  Vcom -Rs*I;

ppAIVValues = IVValuesTrans*10^12;
ppAAvgOCIVq = AvgOCIVq*10^12;
ppAAvgWCIVq=AvgWCIVq*10^12;
ppAAvgWCcorrected=AvgWCcorrected*10^12;

Voltage=[-80;-60;-40;-20;0;20;40;60;80];
VcorrectedmV = Vcorrected*10E2;

% figure()
% plot(AvgWCcorrected)
% 
% figure()
% subplot(2,3,1)
% plot(ppAAvgOCIVq)
% title('control: OC IV')
% hold on
% subplot(2,3,2)
% plot(ppAAvgWCIVq)
% title('control: WC')
% ylabel('Current (pA)')
% hold on
% subplot(2,3,3)
% plot(ppAAvgWCcorrected)
% title('control: corrected')
% ylabel('Current (pA)')
% %xlabel('number of file (in recorded order)')
% hold on
% subplot(2,3,4)
% scatter(Voltage, ppAIVValues)
% hold on
% subplot(2,3,4)
% scatter(VcorrectedmV , ppAIVValues)


ExportTracesVGC = []; ExportVoltageGatedCurrent = [];
ExportVoltageGatedCurrent = [Voltage,VcorrectedmV,IVValuesTrans,IVValuesNormCapac];
ExportTracesVGC = [AvgOCIVq,AvgWCIVq,AvgWCcorrected];

  
save(sprintf('VoltageGC-%s.mat',name)); %save(sprintf('%sTEST.mat',name))

filename = sprintf('VolatgeGatedCurrents-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'Vcom-%s, Vcorrected-%s, IVValues-%s,IVValuesNormCapac-%s \n',name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, ExportVoltageGatedCurrent , '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.

filename = sprintf('TracesVGC-%s.csv',name) ;
fid = fopen(filename, 'w');
fprintf(fid, 'AVGCor1-%s, AVGCor2-%s, AVGCor3-%s,AVGCor4-%s, AVGCor5-%s, AVGCor6-%s,AVGCor7-%s, AVGCor8-%s, AVGCor9-%s \n',name,name,name,name,name,name,name,name,name); %, MergeInd,MeanSameIndCurrent, asdasd, ..\n); %\n means start a new line
fclose(fid);
dlmwrite(filename, AvgWCcorrected, '-append', 'delimiter', '\t'); %Use '\t' to produce tab-delimited files.


else
    display 'no OCIVq'
end

end