% Load touch test results

% [fname, pname] = uigetfile();
[pname, fname] = fileparts('C:\Users\Sammy\Dropbox\Goodman Lab\TouchTests\FAT_Habituation(161010).xlsx');
fname = [fname '.xlsx'];

resultFile = fullfile(pname,fname);
[~,genotypes] = xlsfinfo(resultFile);

touchTest = struct();
allResps = [];
allGenos = [];

for iSheet = 1:length(genotypes)
    thisResp = xlsread(resultFile,genotypes{iSheet});
    allResps = [allResps; thisResp];
    allGenos = [allGenos; repmat(genotypes(iSheet),size(thisResp,1),1)];
end

touchTest.response = allResps;
touchTest.genotype = allGenos;
touchTest.totResp = nansum(allResps,2);
touchTest.color = allGenos; touchTest.color(ismember(allGenos,{'GN212';'GN381'})) = {'fat'};

clear allResps iSheet allGenos fname pname thisFile thisResp;

%% Plot mean responses with gramm

% file:///C:/Users/Sammy/Dropbox/MATLAB/GenericFunctions/gramm/html/examples.html#32
% https://github.com/piermorel/gramm

fatColors = [140 191 216; 239 176 165; 197 221 129];

g = gramm('y',touchTest.totResp,'x',touchTest.genotype,'color',touchTest.color);
g.set_order_options('x',0);
g.set_names('x','Genotype','y','Touch Responses/10');
g.stat_violin('normalization','area','dodge',0,'fill','edge');
g.stat_boxplot('width',0.3);
g.set_title('Touch Test');
g.set_color_options('map','brewer_dark');
g.draw;

%% Plot habituation with gramm
g2 = gramm('y',touchTest.response,'x',touchTest.genotype,'color',touchTest.genotype);
g2.set_order_options('x',0);
g2.stat_summary();
g2.draw; 