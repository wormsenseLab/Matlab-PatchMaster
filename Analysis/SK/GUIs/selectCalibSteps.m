function handles = selectCalibSteps(plotData,cellName)
% Creates a GUI for selecting ranges of data on a plot over which to
% average, in this case for the manual calibration steps.
% 
% USAGE:
% Click Begin Selection, then click on the graph to set a cursor at the
% beginning of the first range. Shift-click to set a cursor at the end of
% the range. To add more cursors, shift-click again, setting pairs of
% cursors down.
% 
% You can move the currently active cursor by simply clicking elsewhere, or
% move any previous cursor (making it the active one) by dragging it. Be
% sure to keep start/end cursor pairs together. Click Clear Cursors to
% remove all existing cursors at any time.
% 
% Click Finish Selection when all cursors are set.
% 
% INPUTS:
%   plotData            A column vector containing the data to be plotted,
%                       from which you will select ranges to average over.
% 
% OUTPUTS:
%   handles             A structure containing, among other things, the
%                       positions and values of the cursors.
% 
% BUGS:
%   Cursor datatips are supposed to keep you updated on which cursor is
%   which so you don't lose pairing, but they don't update properly.
% 
% Created by Sammy Katta on 23 February 2016.

%  Create and then hide the GUI as it is being constructed.
handles.f = figure('Visible','off','Position',[50,250,600,450]);

%  Construct the components.
handles.hText = uicontrol('Style','text','String','How many steps?',...
    'Position',[475,375,90,25]);
handles.hSteps = uicontrol('Style','edit','String','3',...
    'TooltipString','Type a number, then press enter.',...
    'Position',[475,350,90,25],...
    'Callback',@step_Callback);
handles.hBegin = uicontrol('Style','pushbutton',...
    'String','Begin Selection',...
    'Position',[475,300,90,25],...
    'Callback',@beginButton_Callback);
handles.hClear = uicontrol('Style','pushbutton',...
    'String','Clear Cursors',...
    'Position',[475,250,90,25],...
    'CallBack',@clearButton_Callback);
handles.hFinish = uicontrol('Style','pushbutton',...
    'String','Finish Selection',...
    'Position',[475,200,90,25],...
    'Callback',@finishButton_Callback);
handles.hNameText = uicontrol('Style','text',...
    'String',sprintf('Recording: \n %s',cellName), ...
    'Position',[475,100,90,30], ...
    'HorizontalAlignment','left');
handles.ha = axes('Units','Pixels','Position',[50,70,400,350]);
%    align([hText,hSteps,hClear,hFinish],'Center','None');

% Initialize the GUI.
% Change units to normalized so components resize
% automatically.
handles.f.Units = 'normalized';
handles.ha.Units = 'normalized';
handles.hText.Units = 'normalized';
handles.hSteps.Units = 'normalized';
handles.hBegin.Units = 'normalized';
handles.hClear.Units = 'normalized';
handles.hFinish.Units = 'normalized';
handles.hNameText.Units = 'normalized';

handles.f.Position = [0,1,0.8,0.8];

%Create a plot in the axes.
plot(handles.ha, plotData);
box off;
% Assign the GUI a name to appear in the window title.
handles.f.Name = 'Select calibration steps';
% Move the GUI to the center of the screen.
movegui(handles.f,'center')
% Make the GUI visible.
handles.f.Visible = 'on';

uiwait();

%  Callbacks for simple_gui. These callbacks automatically
%  have access to component handles and initialized data
%  because they are nested at a lower level.


    function step_Callback(hObject, eventdata)
        % hObject    handle to edit1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        % Return input text as double
        nSteps = str2double(get(hObject,'string'));
        if isnan(nSteps)
            errordlg('You must enter a numeric value','Invalid Input','modal')
            uicontrol(hObject)
            return
        else
            display(sprintf('%d steps.',nSteps));
        end
        
        
        % Allow data point selection, and create function to label cursor
        % datatips (to avoid mixups in which cursor belongs to which step)
        handles.dcm_obj = datacursormode(handles.f);
        set(handles.dcm_obj,'UpdateFcn',@cursorUpdateFxn)
        set(handles.dcm_obj, 'enable', 'on')
        
        
        %TODO: Fix this. Adding a new cursor updates datatip text for the
        %current cursor plus the previous one (or all the previous ones, if
        %reloading the figure). Or report the bug. Might not be solvable
        %unless you can find the specific datatip handles and write into
        %that string, or to only update the second line of the text.
        function txt = cursorUpdateFxn(~, event_obj)
            % Customizes text of data tips
            % read out data point
            pos = get(event_obj,'Position');
            
            nExistingCursors = size(getCursorInfo(handles.dcm_obj),2);
            
            if mod(nExistingCursors,2) % if this cursor is start of step
                txt = {[num2str((nExistingCursors-1)/2), ' Start'],...
                    ['V: ', num2str(pos(2))]};
                
            else % if this cursor is end of step
                txt = {[num2str(nExistingCursors/2-1), ' End'],...
                    ['V: ', num2str(pos(2))]};
            end
            
        end
        
    end

% Push button callbacks.

    function beginButton_Callback(source,eventdata) %delete all current datatips
        delete(findall(handles.ha,'Type','hggroup','HandleVisibility','off'));
        
                % Allow data point selection, and create function to label cursor
        % datatips (to avoid mixups in which cursor belongs to which step)
        handles.dcm_obj = datacursormode(handles.f);
        set(handles.dcm_obj,'UpdateFcn',@cursorUpdateFxn)
        set(handles.dcm_obj, 'enable', 'on')
        
        
        %TODO: Fix this. Adding a new cursor updates datatip text for the
        %current cursor plus the previous one (or all the previous ones, if
        %reloading the figure). Or report the bug. Might not be solvable
        %unless you can find the specific datatip handles and write into
        %that string, or to only update the second line of the text.
        function txt = cursorUpdateFxn(~, event_obj)
            % Customizes text of data tips
            % read out data point
            pos = get(event_obj,'Position');
            
            nExistingCursors = size(getCursorInfo(handles.dcm_obj),2);
            
            if mod(nExistingCursors,2) % if this cursor is start of step
                txt = {[num2str((nExistingCursors-1)/2), ' Start'],...
                    ['V: ', num2str(pos(2))]};
                
            else % if this cursor is end of step
                txt = {[num2str(nExistingCursors/2-1), ' End'],...
                    ['V: ', num2str(pos(2))]};
            end
            
        end
    end

    function clearButton_Callback(source,eventdata) %delete all current datatips
        delete(findall(handles.ha,'Type','hggroup','HandleVisibility','off'));
    end

    function finishButton_Callback(source,eventdata)
        set(handles.dcm_obj,'enable','Off');
        cursorInfo = fliplr(getCursorInfo(handles.dcm_obj));
        handles.cursorPoints = reshape([cursorInfo.Position],2,[])';
        uiresume()
    end



end