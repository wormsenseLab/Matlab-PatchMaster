%  selectSweepsGUI.m
%
%
%

function keepSweep = selectSweepsGUI(data,leakSize)
nCols = 4;

%TODO: Set figure size so you can actually see the plots (maybe YLim on axes too)

% Create uipanel to hold subplots
%TODO: Add text at bottom w cellname, protname, instructions, finish button
handles.f = figure('Visible','off', ...
    'Units','normalized', 'Position', [0.05 0.3 0.9 0.6]);
handles.uip = uipanel('Position',[0 0.1 1 0.9]);
handles.hFinish = uicontrol('Style','pushbutton',...
    'String','Done', 'Units','normalized',...
    'Position',[0.9 0.025 0.1 0.05],...
    'Callback',@finishButton_Callback);
% Find how many sweeps must be plotted, divide up into rows of n plots each
[~, nSweeps, leak] = size(data);
nRows = fix(nSweeps/nCols)+1;

% Plot the raw traces and the leak subtracted traces for the given sweeps
%TODO: Plot against time
for iSweep = 1:nSweeps
    % Plot the sweep in its proper subplot
    handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
        'Parent', handles.uip);
    plot(data(:,iSweep,1)*1E12);
    hold on;
    % If the leak-subtracted traces are also there, plot them too
    if exist('leak','var')
        plot(data(:,iSweep,2)*1E12);
    end
    
    % Print leak size on subplot
    text(0.5,0.1, sprintf('%.1f pA',leakSize(iSweep)*1E12), ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
        'Units','normalized', 'FontSize',10);
    
    % Number the axis for later use, set the button down function (must be
    % done after plotting or at least after 'hold on', because plotting
    % resets the axes properties otherwise).
    set(gca, 'YLim', [-80 10], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
    
end

set(handles.f, 'Visible', 'on');

% Boolean of which sweeps to keep, will be updated as user clicks on plots
% to exclude those sweeps.
keepSweep = true(1,nSweeps);

uiwait;

    function toggleBGColor(src,~)
        currentSweep = str2double(get(gca,'Tag'));
        if keepSweep(currentSweep)
            set(gca,'Color','Red');
            keepSweep(currentSweep) = false;
        else
            set(gca,'Color','White');
            keepSweep(currentSweep) = true;
        end
    end

    function finishButton_Callback(src,~)
        uiresume;
        close(handles.f);
    end

end

