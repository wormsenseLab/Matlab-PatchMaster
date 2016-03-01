%  selectSweepsGUI.m
%
% Keyboard shortcuts: x = toggle Exclude All
%                     Enter = Done
%

function keepSweep = selectSweepsGUI(data,leakSize)
nCols = 4;

%TODO: Set figure size so you can actually see the plots (maybe YLim on axes too)

% Create uipanel to hold subplots
%TODO: Add text at bottom w cellname, protname, instructions, finish button
handles.f = figure('Visible','off', ...
    'Units','normalized', 'Position', [0.05 0.3 0.9 0.6],...
    'WindowKeyPressFcn',@shortcutKey_Press);
handles.uip = uipanel('Position',[0 0.1 1 0.9]);
handles.hFinish = uicontrol('Style','pushbutton',...
    'String','Done', 'Units','normalized',...
    'Position',[0.85 0.025 0.1 0.05],...
    'Callback',@finishButton_Callback);
handles.hExclude = uicontrol('Style','togglebutton',...
    'String','Exclude All','Units','normalized',...
    'Position',[0.75 0.025 0.1 0.05],...
    'Callback',@excludeButton_Callback);
% Find how many sweeps must be plotted, divide up into rows of n plots each
[~, nSweeps, leak] = size(data);
if mod(nSweeps,nCols) == 0
    nRows = fix(nSweeps/nCols);
else
    nRows = fix(nSweeps/nCols)+1;
end

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

    function excludeButton_Callback(hObject,eventdata)
        button_state = get(hObject,'Value');
        if button_state == get(hObject,'Max')
            for i=1:length(handles.plt)
                set(handles.plt(i),'Color','Red');
            end
            keepSweep(:)=false;
        elseif button_state == get(hObject,'Min')
            for i=1:length(handles.plt)
                set(handles.plt(i),'Color','White');
            end
            keepSweep(:) = true;
            
        end
        
    end


    function shortcutKey_Press(src,eventdata)
        switch eventdata.Key
            case 'return'
                finishButton_Callback()
            case 'x'
                set(handles.hExclude, 'Value', ~handles.hExclude.Value)
                excludeButton_Callback(handles.hExclude,[]);
%             case 'escape'
%                 uiresume;
%                 error('User terminated');
        end
        
    end


end


