% Plots.m

%% Plot current vs. cell-to-stimulus distance, colored by Rs

locDataFat(locDataFat==0)=NaN;
sortloc = sortrows(locDataFat,1); %Sort based on series resistance so you can make a colormap
clr=vals2colormap(sortloc(:,1),'jet',[50 300]); %Create colormap based on range of Rs values
scatter(sortloc(:,3),sortloc(:,4),25,clr,'filled'); %Create scatter and color by Rs values
caxis([50 300]) 
colormap(clr)
cbar=colorbar; %Add colorbar
cbar.Label.String = 'Series Resistance (MOhm)';
xlim([0 250]);
ylim([0 150]);
xlabel('Cell-Stimulus Distance (um)');
ylabel('Current at 9um Step (pA)');
plotfixer

%%
