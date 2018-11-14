%  selectSweepsGUI.m
%
% Keyboard shortcuts: x = toggle Exclude All
%                     Enter = Done
%

function [keepSweep, goBack] = selectSweepsGUI(data,dataType,channel,leakSize,sf,cellName,protName,nCols,pageNo,keepSweep)

goBack = 0; % pass this as 1 if user clicks "Previous" button
%TODO: Set figure size so you can actually see the plots (maybe YLim on axes too)

% Create uipanel to hold subplots
%TODO: Add text at bottom w cellname, protname, instructions, finish button
handles.f = figure('Visible','on', ...
    'Units','normalized', 'Position', [0.05 0.3 0.9 0.6],...
    'WindowKeyPressFcn',@shortcutKey_Press);
handles.uip = uipanel('Position',[0 0.1 1 0.9]);

% Print 
handles.cellTxt = uicontrol('Style','text',...
    'String',sprintf('Recording: %s',cellName), ...
    'Units','normalized', 'Position',[0.01 0.05 0.1 0.03], ...
    'HorizontalAlignment','left');
handles.protTxt = uicontrol('Style','text',...
    'String',sprintf('Protocol %s',protName), ...
    'Units','normalized', 'Position',[0.01 0.02 0.1 0.03],...
    'HorizontalAlignment','left');
handles.pageTxt = uicontrol('Style','text',...
    'String',sprintf('Sweeps %d-%d of %d',pageNo(1),pageNo(2), pageNo(3)), ...
    'Units','normalized', 'Position',[0.15 0.02 0.15 0.03],...
    'HorizontalAlignment','left');

handles.hPrevious = uicontrol('Style','pushbutton',...
    'String','Back to Previous (p)', 'Units','normalized',...
    'Position',[0.35 0.025 0.1 0.05],...
    'Callback',@previousButton_Callback);

handles.hInvert= uicontrol('Style','pushbutton',...
    'String','Invert Selection (n)', 'Units','normalized',...
    'Position',[0.6 0.025 0.1 0.05],...
    'Callback',@invertButton_Callback);

handles.hExclude = uicontrol('Style','togglebutton',...
    'String','Exclude All (x)','Units','normalized',...
    'Position',[0.725 0.025 0.1 0.05],...
    'Callback',@excludeButton_Callback);
handles.hFinish = uicontrol('Style','pushbutton',...
    'String','Done (Enter)', 'Units','normalized',...
    'Position',[0.85 0.025 0.1 0.05],...
    'Callback',@finishButton_Callback);

% Find how many sweeps must be plotted, divide up into rows of n plots each
[~, nSweeps, leak] = size(data);
if mod(nSweeps,nCols) == 0
    nRows = fix(nSweeps/nCols);
else
    nRows = fix(nSweeps/nCols)+1;
end

switch dataType
    
    % Plot the raw current traces and the leak subtracted traces 
    % for the given voltage clamp sweeps
    case 'A'
        for iSweep = 1:nSweeps
            % Plot the sweep in its proper subplot against time
            tVec = (0:length(data(:,iSweep,1))-1)/sf;
            handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
                'Parent', handles.uip);
            plot(tVec,data(:,iSweep,1)*1E12);
            hold on;
            % If the leak-subtracted traces are also there, plot them too
            if exist('leak','var')
                plot(tVec,data(:,iSweep,2)*1E12);
            end
            
            % Print leak size on subplot
            text(0.5,0.1, sprintf('%.1f pA',leakSize(iSweep)*1E12), ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
                'Units','normalized', 'FontSize',10);
            initializeBGColor();
            % Number the axis for later use, set the button down function (must be
            % done after plotting or at least after 'hold on', because plotting
            % resets the axes properties otherwise).
            set(gca, 'YLim', [-80 10], 'XLim',[0 tVec(end)], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
            
        end
     
    % Plot the voltage traces for current clamp sweeps
    case 'V'
        
        if isempty(strfind(protName,'Calib')) && channel == 1 % for non-PD calib traces, current clamp traces
            for iSweep = 1:nSweeps
                % Plot the sweep in its proper subplot against time
                tVec = (0:length(data(:,iSweep,1))-1)/sf;
                handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
                    'Parent', handles.uip);
                plot(tVec,data(:,iSweep,1)*1E3);
                hold on;
                
                % Print baseline voltage on subplot
                text(0.5,0.1, sprintf('%2.0f mV',leakSize(iSweep)*1E3), ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
                    'Units','normalized', 'FontSize',10);
                
                % Number the axis for later use, set the button down function (must be
                % done after plotting or at least after 'hold on', because plotting
                % resets the axes properties otherwise).
                set(gca, 'YLim', [-90 40], 'XLim',[0 tVec(end)], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
                
            end
        elseif isempty(strfind(protName,'Calib')) && channel == 3 % for PD signal for non-calib traces
            for iSweep = 1:nSweeps
                % Plot the sweep in its proper subplot against time
                tVec = (0:length(data(:,iSweep,1))-1)/sf;
                handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
                    'Parent', handles.uip);
                plot(tVec,data(:,iSweep,2)); %leak-subtracted trace
                hold on;
                
                % Print baseline voltage on subplot
                text(0.5,0.1, sprintf('%2.0f mV',leakSize(iSweep)*1E3), ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
                    'Units','normalized', 'FontSize',10);
                
                % Number the axis for later use, set the button down function (must be
                % done after plotting or at least after 'hold on', because plotting
                % resets the axes properties otherwise).
                set(gca, 'YLim', [-3 0.5], 'XLim',[0 tVec(end)], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
                
            end
        else % for PD calib traces
            for iSweep = 1:nSweeps
                % Plot the sweep in its proper subplot against time
                tVec = (0:length(data(:,iSweep,1))-1)/sf;
                handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
                    'Parent', handles.uip);
                plot(tVec,data(:,iSweep,1));
                hold on;
                              
                % Number the axis for later use, set the button down function (must be
                % done after plotting or at least after 'hold on', because plotting
                % resets the axes properties otherwise).
                set(gca, 'YLim',[-5 3], 'XLim',[0 tVec(end)], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
                
            end
        end
        
    case 'mV'
        
        if isempty(strfind(protName,'Calib')) && channel == 2 % for non-PD calib traces
            for iSweep = 1:nSweeps
                % Plot the sweep in its proper subplot against time
                tVec = 0:1/sf:length(data(:,iSweep,1))/sf-1/sf;
                handles.plt(iSweep) = subplot(nRows,nCols,iSweep,...
                    'Parent', handles.uip);
                plot(tVec,data(:,iSweep,1));
                hold on;
                
                % Print baseline voltage on subplot
                text(0.5,0.1, sprintf('%2.0f V',leakSize(iSweep)), ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','center', ...
                    'Units','normalized', 'FontSize',10);
                
                % Number the axis for later use, set the button down function (must be
                % done after plotting or at least after 'hold on', because plotting
                % resets the axes properties otherwise).
                set(gca, 'YLim', [-1 5], 'XLim',[0 tVec(end)], 'Tag', num2str(iSweep),'ButtonDownFcn', @toggleBGColor);
                
            end
        end
end


set(handles.f, 'Visible', 'on');

% Boolean of which sweeps to keep, will be updated as user clicks on plots
% to exclude those sweeps.
% keepSweep = true(1,nSweeps); 

uiwait;

    function initializeBGColor()
        currentSweep = str2double(get(gca,'Tag'));
        if ~keepSweep(iSweep)
            set(gca,'Color','Red');
        else
            set(gca,'Color','White');
        end
    end

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

    function previousButton_Callback(src,~)
        goBack = 1;
        uiresume;
        close(handles.f);
    end

% Invert keepSweep value and background color
    function invertButton_Callback(src,~) 
        keepSweep = ~keepSweep;
        for i=1:length(handles.plt)
            if keepSweep(i)
                set(handles.plt(i),'Color','White');
            else
                set(handles.plt(i),'Color','Red');
            end
        end
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
                finishButton_Callback();
            case 'x'
                set(handles.hExclude, 'Value', ~handles.hExclude.Value);
                excludeButton_Callback(handles.hExclude,[]);
            case 'n'
                invertButton_Callback();
            case 'p'
                previousButton_Callback();                
%             case 'escape'
%                 uiresume;
%                 error('User terminated');
        end
        
        
    end


end


