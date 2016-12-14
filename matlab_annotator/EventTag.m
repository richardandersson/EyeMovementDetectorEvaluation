function varargout = EventTag(varargin)
% EVENTTAG MATLAB code for EventTag.fig
%      EVENTTAG, by itself, creates a new EVENTTAG or raises the existing
%      singleton*.
%
%      H = EVENTTAG returns the handle to a new EVENTTAG or the handle to
%      the existing singleton*.
%
%      EVENTTAG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTTAG.M with the given input arguments.
%
%      EVENTTAG('Property','Value',...) creates a new EVENTTAG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EventTag_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EventTag_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EventTag

% Last Modified by GUIDE v2.5 03-Feb-2012 16:03:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EventTag_OpeningFcn, ...
                   'gui_OutputFcn',  @EventTag_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EventTag is made visible.
function EventTag_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EventTag (see VARARGIN)

% Choose default command line output for EventTag
clc
addpath([cd, filesep ,'HelpFunctions'])
handles.output = hObject;
% screenSize = get(0,'ScreenSize');
% handles.screenRes = [0 1680 0 1024];
% handles.screenDim = [0.38 0.30];
% handles.viewingDist = 0.67;
% handles.Fs = 250;       
handles.XSpan = 200;%*handles.Fs/1000; % 200 ms in number of samples
handles.XSpanIncrement = handles.XSpan/3;%str2double(handles.stepSize);%25;%*handles.Fs/1000;
handles.currentEvent = 'f'; % Fixation as default event
handles.currentEventSpan = handles.XSpan; % samples
handles.markersize = 15;
% handles.slider1 = 0;


% set(handles.axes2,'YLim',[0 40])
% set(handles.axes5,'YLim',[0 handles.screenRes(2)])

% Dialog box to get information about participants
prompt = {'Age:','Gender (m/f):','Experience of eye-tracking (number of years):' };
dlg_title = 'Questions';
num_lines = 1;
def = {'37','m','12'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ID = clock;
ID = [num2str(ID(1)),'_',num2str(ID(2)),'_',num2str(ID(3)),'_',num2str(ID(4))...
    ,'_',num2str(ID(5)),'_',num2str(ID(6))];
handles.answer = answer;
handles.ID = ID;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EventTag wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EventTag_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



% Update plot
[hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Open text-file (only mat-files allowed, contain a variable pos)
[FileName,PathName,FilterIndex] = uigetfile;
handles.PathName = PathName;
handles.FileName = FileName;

% Load a matrix 'pos', which contains eye-tracking data
load([PathName,FileName])

handles.screenRes = [0 ETdata.screenRes(1) 0 ETdata.screenRes(2)];
handles.screenDim = ETdata.screenDim;
handles.viewingDist = ETdata.viewDist;
handles.Fs = ETdata.sampFreq;
pos = ETdata.pos;

% Update Window title with file name
set(gcf, 'name', FileName)

% Add an empty column in the matrix only if data is not tagged
if size(pos,2) < 6
    handles.pos = [pos,zeros(size(pos,1),1)];
else
    handles.pos = pos;
end
handles.nSamples = size(handles.pos,1);

% Step size of slider press
set(handles.slider1,'Min',1)
set(handles.slider1,'Max',handles.nSamples*1000/handles.Fs); % ms
% stepSize = str2double(get(handles.stepSize,'String'))/(handles.nSamples*1000/handles.Fs);
stepSize = handles.XSpanIncrement/(handles.nSamples*1000/handles.Fs);
set(handles.slider1,'SliderStep',[stepSize,10*stepSize])

% Set slider position
set(handles.slider1,'Value',1)

% Start a time to log the time it takes to tag the trial
handles.startTime = clock;

% Set default YLim of velocity plot
set(handles.axes2,'YLim',[0 200]);

% Update plot
[hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);



% function stepSize_Callback(hObject, eventdata, handles)
% % hObject    handle to stepSize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of stepSize as text
% %        str2double(get(hObject,'String')) returns contents of stepSize as a double
% 
% % % Update step size of slider press
% % set(handles.slider1,'Min',1)
% % set(handles.slider1,'Max',handles.nSamples)
% % stepSize = str2double(get(handles.stepSize,'String'))/handles.nSamples;% in ms
% % stepSize = ceil(stepSize*handles.Fs/1000); % in samples
% % set(handles.slider1,'SliderStep',[stepSize,10*stepSize])
% 
% % Step size of slider press
% set(handles.slider1,'Min',1)
% set(handles.slider1,'Max',handles.nSamples*1000/handles.Fs); % ms
% stepSize = str2double(get(hObject,'String'))/(handles.nSamples*1000/handles.Fs);
% set(handles.slider1,'SliderStep',[stepSize,10*stepSize])

% INSTEAD USE 30% OF CURRENT XLIM

% % --- Executes during object creation, after setting all properties.
% function stepSize_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to stepSize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% set(hObject,'String','50')

% --- Executes on button press in zoomButton.
function zoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to zoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

get(hObject,'Value')
if get(hObject,'Value') == 1
    zoom on 
else
    zoom off
end
% Hint: get(hObject,'Value') returns toggle state of zoomButton


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% Turn off zoom
zoom off;

% GET KEY PRESSED
handles.currentEvent = eventdata.Key;
eventMembers = 'fsgpbu';
% Interception for non-event key presses
% For scrolling despite having clicked on data surface
if strcmp(eventdata.Key, 'leftarrow') % if left arrow - scroll left
    set(handles.slider1,'value', round(get(handles.slider1, 'value')-handles.XSpanIncrement)); % set slider
elseif strcmp(eventdata.Key, 'rightarrow') % if right arrow - scroll right
    set(handles.slider1,'value', round(get(handles.slider1,'value')+handles.XSpanIncrement)); % set slider
elseif strcmp(handles.currentEvent, 'home')
    set(handles.slider1, 'value', 1); % go to first value in ms
elseif strcmp(handles.currentEvent, 'end')
    set(handles.slider1, 'value', get(handles.slider1,'Max')) % go to last value in ms
elseif ismember(handles.currentEvent, eventMembers) % if event key press
    [x,y] = ginput(2); % get sample span via mouse clicks
    handles.currentEventSpan = round([x(1),x(2)]*handles.Fs/1000); % in samples
    if handles.currentEventSpan(2) > handles.nSamples % if stop outside limit
        handles.currentEventSpan(2) = handles.nSamples;
    elseif handles.currentEventSpan(1) < 1 % if start outside limit
        handles.currentEventSpan(1) = 1;
    end
    switch handles.currentEvent
        case 'f'
            eventType = 1; % Fixation
        case 's'
            eventType = 2; % Saccade
        case 'g'
            eventType = 3; % Glissade
        case 'p'
            eventType = 4; % pursuit
        case 'b'
            eventType = 5; % Blink
        case 'u'
            eventType = 6; % Undefined
    end
    set(handles.slider1,'value',round(x(2))-handles.XSpanIncrement); % update slider
    % Update pos
    handles.pos(handles.currentEventSpan(1):...
            handles.currentEventSpan(2),6) = eventType;
else % else unexpected key was pressed
    h = warndlg('Unknown event type');
%    disp(handles.currentEvent)
    return
end


% Update plot
[hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles);

% Update event plot
% [hObject, eventdata, handles] = updateEventFigure(hObject, eventdata, handles);


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in saveFile.
function saveFile_Callback(hObject, eventdata, handles)
% hObject    handle to saveFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Display a warning when data is unprocessed (0's)
if length(find(handles.pos(:,end)==0)) > 0
    % Construct a questdlg with three options
%     keyboard
    nonLabelledSamples = num2str(find(handles.pos(:,end)==0)'*1000/handles.Fs); % in ms
    choice = questdlg(['No all samples are labelled (see ',nonLabelledSamples,' ms). Would you really like to quit?'], ...
        'Question', ...
        'Yes','No, I want to continue','Yes');
    % Handle response
    switch choice
        case 'Yes'
            % Do nothing and continue
        case 'No, I want to continue'
            return        
    end
end
% Save data in structure ETdata 
ETdata.pos = handles.pos;
ETdata.screenDim = handles.screenDim;
ETdata.screenRes = [handles.screenRes(2),handles.screenRes(4)];
ETdata.viewDist = handles.viewingDist;
ETdata.sampFreq = handles.Fs;

% Dialog box to get information about participants
prompt = {'How confident are you about your results (1-7)? (1-very unconfident, 7- very confident ):',...
    'Comments/observations?'};
dlg_title = 'Confidence rating';
num_lines = [1 100; 10 100];
def = {'3','My comment'};
answer2 = inputdlg(prompt,dlg_title,num_lines,def);
results.confidence = answer2{1};
results.comments = answer2{2};

% Save info about tagger
fpData.info = handles.answer;
fpData.ID = handles.ID;

results.tagTime = etime(clock,handles.startTime);

uisave({'ETdata','results','fpData'},[handles.PathName,handles.FileName(1:end-4),'_labelled.mat'])

%pos = handles.pos;
%uisave('pos',[handles.PathName,handles.FileName(1:end-4),'_labelled.mat'])
% save([handles.PathName,handles.FileName(1:end-4),'_labelled.mat'],'pos')

%---------------------------------------------------
function [hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles)

% Clear axes
cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes4)
cla(handles.axes5)
cla(handles.axes7)

% Get current interval size (in ms)
offset = handles.XSpan;

% Get position of slider
sliderPos = ceil(get(handles.slider1,'Value'));
if sliderPos <= 0, 
    set(handles.slider1,'Value',1)
    sliderPos = ceil(get(handles.slider1,'Value'));
end
if sliderPos > handles.nSamples*1000/handles.Fs-offset, 
    set(handles.slider1,'Value',handles.nSamples*1000/handles.Fs-offset)
    sliderPos = ceil(get(handles.slider1,'Value'));
end

span = sliderPos:sliderPos+offset-1; % in ms
span = unique(ceil(span*handles.Fs/1000));% in samples

% Define a color for each event
colorTypes = {'r','g','y','m','k','c'};
handles.colorTypes = colorTypes;

%--------------------------------------------------------------------------
% XY-plot
%--------------------------------------------------------------------------
ms = 7;

% Plot XY-data in upper plot
hold(handles.axes1,'on')
plot(handles.axes1,handles.pos(span,4),handles.pos(span,5),'.-');
%plot(handles.axes1,handles.pos(span,4),handles.pos(span,5),'o','markerFaceColor','k','markerSize',ms,'markerEdgeColor','k');

% Mark beginning and end
plot(handles.axes1,handles.pos(span(1),4),handles.pos(span(1),5),'o','MarkerFaceColor','g');
plot(handles.axes1,handles.pos(span(end),4),handles.pos(span(end),5),'o','MarkerFaceColor','r');

% Plot all events in different color
for i = 1:length(colorTypes)
    idx =  find(handles.pos(span,6) == i);
    plot(handles.axes1,handles.pos(span(1)+idx-1,4),handles.pos(span(1)+idx-1,5),[colorTypes{i},'.'],'MarkerSize',handles.markersize);
   % plot(handles.axes1,handles.pos(span(1)+idx-1,4),handles.pos(span(1)+idx-1,5),[colorTypes{i},'o'],'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i});
end
axis(handles.axes1,'ij',handles.screenRes);
hold(handles.axes1,'off')

%--------------------------------------------------------------------------
% XY over time
%--------------------------------------------------------------------------
hold(handles.axes5,'on')
t_ms1 = span/handles.Fs*1000;

% Plot data
plot(handles.axes5,t_ms1,handles.pos(span,4),'-')
plot(handles.axes5,t_ms1,handles.pos(span,5),'--')

% Plot event in different colors
for i = 1:length(colorTypes)
    idx =  find(handles.pos(span,6) == i);
    t_ms = (span(1)+idx-1)/handles.Fs*1000;
    plot(handles.axes5,t_ms,handles.pos(span(1)+idx-1,4),[colorTypes{i},'.'],'MarkerSize',handles.markersize);
    plot(handles.axes5,t_ms,handles.pos(span(1)+idx-1,5),[colorTypes{i},'.'],'MarkerSize',handles.markersize);  
%         plot(handles.axes5,t_ms,handles.pos(span(1)+idx-1,4),[colorTypes{i},'o'],'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i});
%     plot(handles.axes5,t_ms,handles.pos(span(1)+idx-1,5),[colorTypes{i},'o'],'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i});
    
end

hold(handles.axes5,'off')
set(handles.axes5,'XLim',[t_ms1(1),t_ms1(end)])

% Cut axis if data outside of screen
axis(handles.axes5,'auto y')

if (max(handles.pos(span,4)) > handles.screenRes(2)) | ...
   (max(handles.pos(span,5)) > handles.screenRes(4))
    set(handles.axes5,'YLim',[0,handles.screenRes(2)])
end


%xlabel(handles.axes5,'Time (ms)')
ylabel(handles.axes5,'Coordinates (pixels)')
legend(handles.axes5,'x-coord','y-coord')

%--------------------------------------------------------------------------
% Velocity over time
%--------------------------------------------------------------------------

% Plot velocity-data in lower plot
% handles.vel = smooth(sqrt(diff(handles.pos(span,2)).^2+diff(handles.pos(span,3)).^2));
% keyboard

handles.vel = sqrt(diff(handles.pos(span,4)).^2+diff(handles.pos(span,5)).^2);
[alphaH,alphaV] = pixels2degrees(handles.vel,handles.viewingDist,[handles.screenRes(2),handles.screenRes(4)],handles.screenDim);
handles.vel = alphaH*handles.Fs;
handles.vel = [handles.vel(1);handles.vel]; % use to convert to degrees second./handles.sampFreq;
hold(handles.axes2,'on')

plot(handles.axes2,t_ms1,handles.vel,'.-')

for i = 1:length(colorTypes)
    idx =  find(handles.pos(span,6) == i);
    t_ms = (span(1)+idx-1)/handles.Fs*1000;
    
    %     idx = idx(1:end-1);
    plot(handles.axes2,t_ms,handles.vel(idx),[colorTypes{i},'.'],...
        'Linewidth',2,'MarkerSize',handles.markersize);
%       plot(handles.axes2,t_ms,handles.vel(idx),[colorTypes{i},'o'],...
%         'Linewidth',2,'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i});
    %     plot(handles.axes2,span(1)+idx-1,handles.vel(idx),[colorTypes{i},'.'],...

    
end


%If labelled by the sp-algorithm, the decision is plotted here
if size(handles.pos,2) == 7
        
        index = find(handles.pos(span,7));
        
        if ~isempty(index)
            
            for k = 1:length(index)
                
                decision = handles.pos(span(index(k)),7);
                textYpos = handles.vel(index(k)) + 10;
                textXpos = span(index(k))/handles.Fs*1000;
                axes(handles.axes2)
                text(textXpos,textYpos,num2str(decision));
            end
            
        end
end
    

hold(handles.axes2,'off')
set(handles.axes2,'XLim',[t_ms1(1),t_ms1(end)])
%xlabel(handles.axes2,'Time (ms)')
ylabel(handles.axes2,'Velocity (deg/s)')



%--------------------------------------------------------------------------
% Pupil size (y-coordinate) over time
%--------------------------------------------------------------------------

hold(handles.axes7,'on')
plot(handles.axes7,t_ms1,handles.pos(span,3),'.-')
% plot(handles.axes7,t_ms1,handles.pos(span,3),'o','markerFaceColor','k','markerSize',ms,'markerEdgeColor','k')


for i = 1:length(colorTypes)
    idx =  find(handles.pos(span,6) == i);
    t_ms = (span(1)+idx-1)/handles.Fs*1000;
    plot(handles.axes7,t_ms,handles.pos(span(1)+idx-1,3),[colorTypes{i},'.'],...
        'Linewidth',2,'MarkerSize',handles.markersize);
%     plot(handles.axes7,t_ms,handles.pos(span(1)+idx-1,3),[colorTypes{i},'o'],...
%         'Linewidth',2,'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i})
end
hold(handles.axes7,'off')
set(handles.axes7,'XLim',[t_ms1(1),t_ms1(end)])
xlabel(handles.axes7,'Time (ms)')
ylabel(handles.axes7,'Vertical pupil diameter (pixels)')


%--------------------------------------------------------------------------
% Zoomed in portion of data (axes 4)
%--------------------------------------------------------------------------

% Zoom in on saccade end after having labelled a saccade
samplesBeforeSaccade = 2;
samplesAfterSaccade =  30;
if strcmp(handles.currentEvent,'s')==1 
    span = (max(1,handles.currentEventSpan(2)-samplesBeforeSaccade)):...
           min(handles.nSamples,(handles.currentEventSpan(2)+samplesAfterSaccade));
end
samplesBeforeSaccade = 10;
if strcmp(handles.currentEvent,'g')==1
    span = (max(1,handles.currentEventSpan(2)-samplesBeforeSaccade)):...
           min(handles.nSamples,(handles.currentEventSpan(2)+samplesAfterSaccade));
end

% Plot XY-data 
hold(handles.axes4,'on')
plot(handles.axes4,handles.pos(span,4),handles.pos(span,5),'.-');
% plot(handles.axes4,handles.pos(span,4),handles.pos(span,5),'o','Linewidth',2,'markerFaceColor','k','markerSize',ms,'markerEdgeColor','k');

% Mark beginning and end
plot(handles.axes4,handles.pos(span(1),4),handles.pos(span(1),5),'o','MarkerFaceColor','g','markerSize',ms);
plot(handles.axes4,handles.pos(span(end),4),handles.pos(span(end),5),'o','MarkerFaceColor','r','markerSize',ms);

% Plot all events in different color
for i = 1:length(colorTypes)
    idx =  find(handles.pos(span,6) == i);
    plot(handles.axes4,handles.pos(span(1)+idx-1,4),handles.pos(span(1)+idx-1,5),[colorTypes{i},'.']);
%     plot(handles.axes4,handles.pos(span(1)+idx-1,4),handles.pos(span(1)+idx-1,5),[colorTypes{i},'o'],'markerFaceColor',colorTypes{i},'markerSize',ms,'markerEdgeColor',colorTypes{i});

end

% zoomSize = [min(handles.pos(eventSpan,4)),max(handles.pos(eventSpan,4)),...
%     min(handles.pos(eventSpan,5)),max(handles.pos(eventSpan,5))];
% if zoomSize(2)-zoomSize(1) <= 0 | zoomSize(4)-zoomSize(3) <= 0
%     zoomSize = [0 1 0 1];
% end


title(handles.axes4,'Zoomed in data')
axis(handles.axes4,'ij')%;,zoomSize);
axis(handles.axes4,'image')
% axis(handles.axes4,'ij',handles.screenRes);
% hold(handles.axes4,'off')


%---------------------------------------------------
% function [hObject, eventdata, handles] = updateEventFigure(hObject, eventdata, handles)
% 
% cla(handles.axes4)
% 
% 
% % Calculate event span and span with offset
% offset = 0.20; % Add some samples on each side of the event
% eventSpan = handles.currentEventSpan(1):handles.currentEventSpan(2);
% eventSpanOffset = (handles.currentEventSpan(1)-offset*length(eventSpan)):...
%     (handles.currentEventSpan(2)+offset*length(eventSpan));
% 
% % Take care of edge problems
% eventSpan = eventSpan(find(eventSpan > 0 & eventSpan <= handles.nSamples));
% eventSpanOffset = ceil(eventSpanOffset(find(eventSpanOffset > 0 & eventSpanOffset <= handles.nSamples)));
% 
% % Plot zoomed in portion in small plot
% hold(handles.axes4,'on')
% plot(handles.axes4,handles.pos(eventSpanOffset,4),handles.pos(eventSpanOffset,5),'.-b');
% plot(handles.axes4,handles.pos(eventSpan,4),handles.pos(eventSpan,5),[handles.colorTypes{handles.pos(eventSpan(1),6)},'.-']);
% hold(handles.axes4,'off')
% 
% zoomSize = [min(handles.pos(eventSpan,4)),max(handles.pos(eventSpan,4)),...
%     min(handles.pos(eventSpan,5)),max(handles.pos(eventSpan,5))];
% if zoomSize(2)-zoomSize(1) <= 0 | zoomSize(4)-zoomSize(3) <= 0
%     zoomSize = [0 1 0 1];
% end
% 
% 
% title(handles.axes4,'Most recently defined (sub)event')
% axis(handles.axes4,'ij',zoomSize);
% axis(handles.axes4,'image')


% --- Executes on button press in ylimp.
function ylimp_Callback(hObject, eventdata, handles)
% hObject    handle to ylimp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldY = get(handles.axes2,'YLim');
set(handles.axes2,'YLim',[0 oldY(2)+30])


% --- Executes on button press in ylimm.
function ylimm_Callback(hObject, eventdata, handles)
% hObject    handle to ylimm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldY = get(handles.axes2,'YLim');
set(handles.axes2,'YLim',[0 max(oldY(2)-30,1)])


% --- Executes on button press in XYm.
function XYm_Callback(hObject, eventdata, handles)
% hObject    handle to XYm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldY = get(handles.axes5,'YLim');
set(handles.axes5,'YLim',[min(oldY(1)+50,oldY(2)-50), max(oldY(1)+50,oldY(2)-50)])


% --- Executes on button press in XYp.
function XYp_Callback(hObject, eventdata, handles)
% hObject    handle to XYp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldY = get(handles.axes5,'YLim');
set(handles.axes5,'YLim',[oldY(1)-50 oldY(2)+50])


% --- Executes on button press in XM.
function XM_Callback(hObject, eventdata, handles)
% hObject    handle to XM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oldX = get(handles.axes5,'XLim');
newSpan = [oldX(1)+handles.XSpanIncrement oldX(2)-handles.XSpanIncrement];

% Make sure span does not exceed allowed interval
newSpan = [max(1,newSpan(1)),min(handles.nSamples*1000/handles.Fs,newSpan(2))];

% Make sure span is of positive length. If not, do nothing and return.
if  newSpan(2)-newSpan(1) <= 1000/handles.Fs
    return
end

% Set new limits
set(handles.axes5,'XLim',newSpan)
set(handles.axes2,'XLim',newSpan)

% Update spans
handles.XSpan = newSpan(2)-newSpan(1);
handles.XSpanIncrement = handles.XSpan/3;


% Update slider position
set(handles.slider1,'Value',newSpan(1))

% Step size of slider press
set(handles.slider1,'Min',1)
set(handles.slider1,'Max',handles.nSamples*1000/handles.Fs); % ms
stepSize = handles.XSpanIncrement/(handles.nSamples*1000/handles.Fs);
set(handles.slider1,'SliderStep',[stepSize,10*stepSize])

% Update handles structure
guidata(hObject, handles);

[hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles);

% --- Executes on button press in XP.
function XP_Callback(hObject, eventdata, handles)
% hObject    handle to XP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldX = get(handles.axes5,'XLim');
newSpan = [oldX(1)-handles.XSpanIncrement oldX(2)+handles.XSpanIncrement];
newSpan = [max(1,newSpan(1)),min(handles.nSamples*1000/handles.Fs,newSpan(2))];
set(handles.axes5,'XLim',newSpan)
set(handles.axes2,'XLim',newSpan)

handles.XSpan = newSpan(2)-newSpan(1);
handles.XSpanIncrement = handles.XSpan/3;

% Update slider position
set(handles.slider1,'Value',newSpan(1))

% Step size of slider press
set(handles.slider1,'Min',1)
set(handles.slider1,'Max',handles.nSamples*1000/handles.Fs); % ms
stepSize = handles.XSpanIncrement/(handles.nSamples*1000/handles.Fs);
set(handles.slider1,'SliderStep',[stepSize,10*stepSize])

% Update handles structure
guidata(hObject, handles);

[hObject, eventdata, handles] = updateFigure(hObject, eventdata, handles);

% --- Executes on button press in pupilM.
function pupilM_Callback(hObject, eventdata, handles)
% hObject    handle to pupilM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pupilP.
function pupilP_Callback(hObject, eventdata, handles)
% hObject    handle to pupilP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
