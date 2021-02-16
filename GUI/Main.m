function varargout = Main(varargin)
% MAIN MATLAB code for Main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main

% Last Modified by GUIDE v2.5 26-Nov-2013 09:51:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Main_OpeningFcn, ...
                   'gui_OutputFcn',  @Main_OutputFcn, ...
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


% --- Executes just before Main is made visible.
function Main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main (see VARARGIN)

% Choose default command line output for Main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

dir_main = dir(pwd);
cell_files = {};
for k = 1 : length(dir_main)
       cell_files{k,1} = dir_main(k,1).name;
end
set(handles.popup_file,'String',cell_files);
% UIWAIT makes Main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

enable_all_handles(handles,'Off')
set(handles.popup_file,'Enable','On')
set(handles.push_load,'Enable','On')


% --- Outputs from this function are returned to the command line.
function varargout = Main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in list_experiment.
function list_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to list_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_experiment


% --- Executes during object creation, after setting all properties.
function list_experiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_file.
function popup_file_Callback(hObject, eventdata, handles)
% hObject    handle to popup_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_file contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_file


% --- Executes during object creation, after setting all properties.
function popup_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_load.
function push_load_Callback(hObject, eventdata, handles)
% hObject    handle to push_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
filenames = get(handles.popup_file,'String');
n = get(handles.popup_file,'Value');

try
load(filenames{n,1});
disp([filenames{n,1} ' loaded.'])
enable_all_handles(handles,'On')
[path name ext] = fileparts(filenames{n,1});

set(handles.edit_save,'String',[name '_processed.mat'])

catch err
    disp('Error loading. Probably invalid file.')
    enable_all_handles(handles,'Off')
    set(handles.popup_file,'Enable','On')
    set(handles.push_load,'Enable','On')
end





% --- Executes on selection change in list_user.
function list_user_Callback(hObject, eventdata, handles)
% hObject    handle to list_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_user contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_user


% --- Executes during object creation, after setting all properties.
function list_user_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_signal.
function list_signal_Callback(hObject, eventdata, handles)
% hObject    handle to list_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_signal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_signal


% --- Executes during object creation, after setting all properties.
function list_signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in edit_ratio.
function edit_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns edit_ratio contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_ratio


% --- Executes during object creation, after setting all properties.
function edit_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
save(get(handles.edit_save,'String'),'Data')


% --- Executes on button press in push_load_dataset.
function push_load_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to push_load_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Data
global u_data
global y_data
global Fs
global I_beg
global I_end
global u_line
global y_line
global experiment
global user
global signal
global Bode_Mode

Bode_Mode = 'Amplitude';
experiments = get(handles.list_experiment,'String');
users = get(handles.list_user,'String');
signals = get(handles.list_signal,'String');

experiment = experiments{get(handles.list_experiment,'Value')};
user = users{get(handles.list_user,'Value')};
signal = signals{get(handles.list_signal,'Value')};

u_data = Data.(experiment).(user).(signal).Data(:,2);
y_data = Data.(experiment).(user).(signal).Data(:,1);
Fs = Data.(experiment).(user).(signal).Fs;


I_beg = 1;
I_end = length(u_data);
set(handles.edit_crop_beg,'String',num2str(I_beg))
set(handles.edit_crop_end,'String',num2str(I_end))

t = (0 : length(u_data)-1)/Fs;

t_beg = [t(I_beg) t(I_beg)];
t_end = [t(I_end) t(I_end)];
u_line = [min(u_data) max(u_data)];
y_line = [min(y_data) max(y_data)];

plot(handles.axes1,t,u_data, t_beg, u_line, t_end, u_line)
axis(handles.axes1,'tight')
set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
plot(handles.axes2,t,y_data, t_beg, y_line , t_end, y_line)
axis(handles.axes2,'tight')
set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')
%set(handles.text_datapoints,'String',num2str(length(u_data)))
%set(handles.text_datapoints_per_sample,'String',num2str(length(u_data)))
%set(handles.text_Fs,'String',num2str(Fs))


% --- Executes on slider movement.
function slider_beg_Callback(hObject, eventdata, handles)
% hObject    handle to slider_beg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global u_data
global y_data
global I_beg
global I_end
global u_line
global y_line

global Fs

pos_beg = get(hObject,'Value');
I_beg = round(pos_beg*length(u_data));
if(I_beg == 0)
   I_beg = 1; 
end

set(handles.edit_crop_beg,'String',num2str(I_beg))

t = (0:1:length(u_data)-1)/Fs;
t_beg = [t(I_beg),t(I_beg)];
t_end = [t(I_end),t(I_end)];


plot(handles.axes1,t,u_data, t_beg, u_line, t_end, u_line)
axis(handles.axes1,'tight')
set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
plot(handles.axes2,t,y_data, t_beg, y_line, t_end, y_line)
axis(handles.axes2,'tight')
set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')

% --- Executes during object creation, after setting all properties.
function slider_beg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_beg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_end_Callback(hObject, eventdata, handles)
% hObject    handle to slider_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global u_data
global y_data
global I_beg
global I_end
global u_line
global y_line
global Fs

pos_end = get(hObject,'Value');
I_end = round(pos_end*length(u_data));
if(I_end == 0)
   I_end = 1; 
end

set(handles.edit_crop_end,'String',num2str(I_end))

t = (0:1:length(u_data)-1)/Fs;
t_beg = [t(I_beg),t(I_beg)];
t_end = [t(I_end),t(I_end)];

plot(handles.axes1, t, u_data, t_beg, u_line, t_end, u_line)
axis(handles.axes1,'tight')
set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
plot(handles.axes2, t, y_data, t_beg, y_line, t_end, y_line)
axis(handles.axes2,'tight')
set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')


% --- Executes during object creation, after setting all properties.
function slider_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_delta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_delta as text
%        str2double(get(hObject,'String')) returns contents of edit_delta as a double


% --- Executes during object creation, after setting all properties.
function edit_delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_save_Callback(hObject, eventdata, handles)
% hObject    handle to edit_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_save as text
%        str2double(get(hObject,'String')) returns contents of edit_save as a double


% --- Executes during object creation, after setting all properties.
function edit_save_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_filter.
function popup_filter_Callback(hObject, eventdata, handles)
% hObject    handle to popup_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_filter


% --- Executes during object creation, after setting all properties.
function popup_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in edit_subdivide.
function edit_subdivide_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subdivide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns edit_subdivide contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edit_subdivide


% --- Executes during object creation, after setting all properties.
function edit_subdivide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subdivide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_export.
function push_export_Callback(hObject, eventdata, handles)
% hObject    handle to push_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global experiment
global user
global signal
global u_model
global y_model
global u_batch
global y_batch
global U_val
global Y_val
global u_val
global y_val
global u_val_batch
global y_val_batch
global Y_batch
global U_batch
global G_batch
global n_sub
global N_sample
global G_model
global Y_model
global U_model
global n_sub_val
global G_smooth
global Fs
global gamma
global delta
global filter
global method

Data = struct('u_batch',u_batch,'y_batch',y_batch,'U_val',U_val,'Y_val',Y_val,'U_model',U_model,'n_sub',n_sub);
save('test.mat','-struct','Data');

% 
% %plot generation
% n_sub = str2num(get(handles.edit_subdivide,'String'));
% 
% t = (0:1:length(u_data)-1)*(1/Fs);
% freq = (Fs/length(u_data))*(0:1:ceil(length(U_data)/2));
% 
% 
% % u = u_data(1:length(t));
% % y = y_data(1:length(t));
% % U = U_data(1:length(freq));
% % Y = Y_data(1:length(freq));
% % G = G_data(1:length(freq));
% % G_s = G_smooth(1:length(freq));
% % show_in_ws(u);
% % show_in_ws(y);
% % show_in_ws(U);
% % show_in_ws(Y);
% % show_in_ws(G);
% % show_in_ws(G_s);
% 
% if(strcmp(method,'Smooth - Average'))
%     mth = 'SA';
% else
%     mth = 'AS';
% end
% 
% file_name = ['Plots\' experiment '_' user '_' signal '_' num2str(round(Fs))...
%     '_' filter '_' num2str(gamma) '_' num2str(delta) '_' mth '_' num2str(n_sub)];
% 
% Data = struct( 'u',u_data,'y',y_data,'U',U_data,'Y',Y_data,'G',G_data,...
%     'G_sm',G_smooth,'Gamma',gamma,'Delta',delta,'Filter',filter,...
%     'Method',method,'Averaging',n_sub,'Fs',Fs,'experiment',...
%     experiment,'user', user, 'signal', signal);
% 
% save([file_name  '.mat'] ,'-struct','Data')
% clear Data
% 
% f1 = figure(1);
% set(f1,'Position', [0 0 800 1700])
% set(f1,'PaperPositionMode','auto')
% set(f1,'PaperType','A4')
% set(f1, 'visible','off');
% 
% %u(t)
% subplot(8,1,1)
% plot(t,u_data(1:length(t)))
% xlabel('Time [s]')
% ylabel('u(t)')
% grid on
% grid minor
% axis tight
% title_cell =     {  [experiment '\_' user '\_' signal],...
%                     ['Fs: ' num2str(Fs) ' Hz    Averaging: ' num2str(n_sub)], ...
%                     ['Filter: ' filter '   \gamma = ' num2str(gamma) '   \delta = ' num2str(delta)],...
%                     ['Method: ' method]}';
% title(title_cell)
% 
% disp(title_cell)
% 
% %y(t)
% subplot(8,1,2)
% plot(t,y_data(1:length(t)))
% xlabel('Time [s]')
% ylabel('y(t)')
% grid on
% grid minor
% axis tight
% 
% %|U(w)|
% subplot(8,1,3)
% loglog(freq,abs(U_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('|U(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(U(w))
% subplot(8,1,4)
% semilogx(freq,angle(U_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('\angle U(\omega)')
% grid on
% grid minor
% axis tight
% 
% %|Y(w)|
% subplot(8,1,5)
% loglog(freq,abs(Y_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel(' |Y(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(Y(w))
% subplot(8,1,6)
% semilogx(freq,angle(Y_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('\angle Y(\omega)')
% grid on
% grid minor
% axis tight
% 
% %|G(w)|
% subplot(8,1,7)
% h = loglog(freq,abs(G_data(1:length(freq))),freq,abs(G_smooth(1:length(freq))),'r');
% set(h(2),'LineWidth',1.1);
% xlabel('Frequency [Hz]')
% ylabel('|G(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(G(w))
% subplot(8,1,8)
% h = semilogx(freq,angle(G_data(1:length(freq))),freq,angle(G_smooth(1:length(freq))),'r');
% set(h(2),'LineWidth',1.1);
% xlabel('Frequency [Hz]')
% ylabel('\angle G(\omega)')
% grid on
% grid minor
% axis tight
% 
% 
% print(f1,file_name,'-dpdf')
% close(f1)
% disp('File Printed')


% --- Executes on selection change in popup_method.
function popup_method_Callback(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_method


% --- Executes during object creation, after setting all properties.
function popup_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_validation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_validation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_validation as text
%        str2double(get(hObject,'String')) returns contents of edit_validation as a double


% --- Executes during object creation, after setting all properties.
function edit_validation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_validation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_crop_beg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crop_beg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crop_beg as text
%        str2double(get(hObject,'String')) returns contents of edit_crop_beg as a double
global u_data
global y_data
global I_beg
global I_end
global u_line
global y_line
global Fs

I_beg = str2num(get(hObject,'String'));
set(handles.slider_beg,'Value',I_beg/length(u_data))

t = (0:1:length(u_data)-1)/Fs;
t_beg = [t(I_beg),t(I_beg)];
t_end = [t(I_end),t(I_end)];

plot(handles.axes1, t, u_data, t_beg, u_line, t_end, u_line)
axis(handles.axes1,'tight')
set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
plot(handles.axes2, t, y_data, t_beg, y_line, t_end, y_line)
axis(handles.axes2,'tight')
set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')



% --- Executes during object creation, after setting all properties.
function edit_crop_beg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_beg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_crop_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crop_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crop_end as text
%        str2double(get(hObject,'String')) returns contents of edit_crop_end as a double
global u_data
global y_data
global I_beg
global I_end
global u_line
global y_line
global Fs

I_end = str2num(get(hObject,'String'));
set(handles.slider_end,'Value',I_end/length(u_data))

t = (0:1:length(u_data)-1)/Fs;
t_beg = [t(I_beg),t(I_beg)];
t_end = [t(I_end),t(I_end)];

plot(handles.axes1, t, u_data, t_beg, u_line, t_end, u_line)
axis(handles.axes1,'tight')
set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
plot(handles.axes2, t, y_data, t_beg, y_line, t_end, y_line)
axis(handles.axes2,'tight')
set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')


% --- Executes during object creation, after setting all properties.
function edit_crop_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_analyse.
function push_analyse_Callback(hObject, eventdata, handles)
% hObject    handle to push_analyse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Global Variables
global u_data
global y_data
global I_beg
global I_end
global Fs
global Bode_Mode
global freq_plot
global x_min
global x_max

%
global U_model_av
global Y_model_av
global G_model_av
global U_val_av
global Y_val_av
global G_model_smooth_av
global Data

%Analysis Variables
thinning_ratio = str2num(get(handles.edit_ratio,'String'));
I_beg = str2num(get(handles.edit_crop_beg,'String'));
I_end = str2num(get(handles.edit_crop_end,'String'));
val = str2num(get(handles.edit_validation,'String'))/100;
n_sub = str2num(get(handles.edit_subdivide,'String'));
methods = get(handles.popup_method,'String');
method = methods{get(handles.popup_method,'Value')};
gamma = str2num(get(handles.edit_gamma,'String'));
delta = str2num(get(handles.edit_delta,'String'));
filters = get(handles.popup_filter,'String');
filter = filters{get(handles.popup_filter,'Value')};


%Cropping
disp('Cropping...')
u_model = u_data(I_beg:I_end);
y_model = y_data(I_beg:I_end);

%Thinning:
disp('Thinning...')
u_model = matrix_thinner(u_model,thinning_ratio);
y_model = matrix_thinner(y_model,thinning_ratio);
Fs = Fs/thinning_ratio;

%Validation Model
disp('Validation Model...')
N_model = floor((1-val)*length(u_model));

u_val = u_model(N_model + 1 : length(u_model));
y_val = y_model(N_model + 1: length(y_model));
u_model = u_model(1: N_model);
y_model = y_model(1: N_model);

%Batch Generation
disp('Batch Generation...')
N_sample = floor(N_model/n_sub);
N_val = length(u_val);
n_sub_val = floor(N_val/N_sample);

u_model_batch = reshape(u_model(1:N_sample*n_sub),N_sample,n_sub);
y_model_batch = reshape(y_model(1:N_sample*n_sub),N_sample,n_sub);
u_val_batch = reshape(u_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);
y_val_batch = reshape(y_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);

%Fourier Transform
disp('Fourier Transform...')
U_model_batch = zeros(size(u_model_batch));
Y_model_batch = zeros(size(y_model_batch));
for k = 1 : n_sub
    U_model_batch(:,k) = (1/N_sample)*fft(u_model_batch(:,k));
    Y_model_batch(:,k) = (1/N_sample)*fft(y_model_batch(:,k));
end

U_model_av = zeros(N_sample,1);
Y_model_av = zeros(N_sample,1);
for k = 1 : n_sub
    U_model_av = U_model_av + U_model_batch(:,k);
    Y_model_av = Y_model_av + Y_model_batch(:,k);
end
U_model_av = (1/n_sub)*U_model_av;
Y_model_av = (1/n_sub)*Y_model_av;

U_val_batch = zeros(size(u_val_batch));
Y_val_batch = zeros(size(y_val_batch));
for k = 1 : n_sub_val
    U_val_batch(:,k) = (1/N_sample)*fft(u_val_batch(:,k));
    Y_val_batch(:,k) = (1/N_sample)*fft(y_val_batch(:,k));
end

U_val_av = zeros(N_sample,1);
Y_val_av = zeros(N_sample,1);

for k = 1 : n_sub_val
    U_val_av = U_val_av + U_val_batch(:,k);
    Y_val_av = Y_val_av + Y_val_batch(:,k);
end
U_val_av = (1/n_sub_val)*U_val_av;
Y_val_av = (1/n_sub_val)*Y_val_av;


%TF Estimate
disp('Transfer Function Estimates...')
freq = (2*pi/N_sample)*(0:1:N_sample-1);


G_model_batch = zeros(size(U_model_batch));
G_model_smooth_batch = zeros(size(u_model_batch));

if(strcmp(method,'Smooth - Average'))
    % Generate G_batch
    for k = 1 : n_sub
        G_model_batch(:,k) = Y_model_batch(:,k)./U_model_batch(:,k);
    end
    
    % Average G_batch
    G_model_av = zeros(N_sample,1);
    for k = 1 : n_sub
        G_model_av = G_model_av + G_model_batch(:,k);
    end
    G_model_av = (1/n_sub)*G_model_av;
    
    
    % Smoothing of every G
    
    for k = 1 : n_sub
        [Psy_yu, Psy_u] = spect_filtered(u_model_batch(:,k),y_model_batch(:,k),freq,delta,gamma,filter,handles);
        G_model_smooth_batch(:,k) = Psy_yu./Psy_u;
        disp([num2str(k) '/' num2str(n_sub) ' Finished.'])
    end
    
    %Averaging of smoothed G
    G_model_smooth_av = zeros(N_sample,1);
    for k = 1 : n_sub
        G_model_smooth_av = G_model_smooth_av + G_model_smooth_batch(:,k);
    end
    G_model_smooth_av = (1/n_sub)*G_model_smooth_av;
            
        
else
    %Calculate G_model out of averaged Y_model and U_model
    G_model_av = Y_model_av./U_model_av;
    
    %Smooth averaged U and Y
    [Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq,gamma,delta,filter,handles);
    G_model_smooth_av = (Psy_yu./Psy_u)';
    
end


%Validation
disp('Validation...')
%Y_val_est = zeros(N_sample,1);
Y_val_est = G_model_smooth_av.*U_val_av;
% RMS_abs = sqrt((1/N_sample)*sum((abs(Y_val_est)-abs(Y_val_av)).^2));
% RMS_angle = sqrt((1/N_sample)*sum((angle(Y_val_est)-angle(Y_val_av)).^2));
RMS = sqrt((1/N_sample)*sum((abs(Y_val_est-Y_val_av)).^2));

%Data struct
disp('Creating structure...')
Parameters = struct('Thinning',thinning_ratio,...
                            'I_beg',I_beg,...
                            'I_end',I_end,...
                            'Validation',val,...
                            'Method',method,...
                            'Gamma',gamma,...
                            'Delta',delta,...
                            'Filter',filter);
        
        Raw = struct(   'u_data',u_data,...
                        'y_data',y_data);
        
        Model = struct( 'u_model',u_model,...
                        'y_model',y_model,...
                        'u_model_batch',u_model_batch,...
                        'y_model_batch',y_model_batch,...
                        'U_model_batch',U_model_batch,...
                        'Y_model_batch',Y_model_batch,...
                        'U_model_av',U_model_av,...
                        'Y_model_av',Y_model_av,...
                        'G_model_batch',G_model_batch,...
                        'G_model_smooth_batch',G_model_smooth_batch,...
                        'G_model_av',G_model_av,...
                        'G_model_smooth_av',G_model_smooth_av);
                    
        Validation = struct('u_val',u_val,...
                        'y_val',y_val,...
                        'u_val_batch',u_val_batch,...
                        'y_val_batch',y_val_batch,...
                        'U_val_batch',U_val_batch,...
                        'Y_val_batch',Y_val_batch,...
                        'U_val_av',U_val_av,...
                        'Y_val_av',Y_val_av,...
                        'Y_val_est',Y_val_est,...
                        'RMS',RMS);

                    Data = struct('Parameters',Parameters,'Raw',Raw,'Model',Model,'Validation',Validation);
                    show_in_ws(Data)
disp('Finished.')                   
%Plots
freq_plot = (Fs/N_sample)*(0:ceil(N_sample/2));
x_min = 20;
x_max = 20e3;
if(freq_plot(length(freq_plot))< x_max)
   x_max = freq_plot(length(freq_plot)); 
end


if(strcmp(Bode_Mode,'Amplitude'))
loglog(handles.axes3,freq_plot,abs(U_model_av(1:length(freq_plot))),freq_plot,abs(U_val_av(1:length(freq_plot))))
loglog(handles.axes4,freq_plot,abs(Y_model_av(1:length(freq_plot))),freq_plot,abs(Y_val_av(1:length(freq_plot))))
loglog(handles.axes5,freq_plot,abs(G_model_av(1:length(freq_plot))),freq_plot,abs(G_model_smooth_av(1:length(freq_plot))))
title_U = '|U(\omega)|';
title_Y = '|Y(\omega)|';
title_G = '|G(\omega)|';

else
    
semilogx(handles.axes3,freq_plot,angle(U_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(U_val_av(1:length(freq_plot)))*180/pi)
semilogx(handles.axes4,freq_plot,angle(Y_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(Y_val_av(1:length(freq_plot)))*180/pi)
semilogx(handles.axes5,freq_plot,angle(G_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(G_model_smooth_av(1:length(freq_plot)))*180/pi)
title_U = '\angle U(\omega)';
title_Y = '\angle Y(\omega)';
title_G = '\angle G(\omega)';

end

set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
axis(handles.axes1,'tight')
xlabel(handles.axes1,'Time [s]')
title(handles.axes1,'u(t)')

set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')
axis(handles.axes2,'tight')
xlabel(handles.axes2,'Time [s]')
title(handles.axes2,'y(t)')

set(handles.axes3,'XGrid','on')
set(handles.axes3,'YGrid','on')
set(handles.axes3,'XMinorGrid','on')
set(handles.axes3,'YMinorGrid','on')
axis(handles.axes3,'tight')
xlabel(handles.axes3,'Frequency [Hz]')
title(handles.axes3,title_U)
xlim(handles.axes3,[x_min x_max])
legend(handles.axes3,'System','Validation','Location','Best')

set(handles.axes4,'XGrid','on')
set(handles.axes4,'YGrid','on')
set(handles.axes4,'XMinorGrid','on')
set(handles.axes4,'YMinorGrid','on')
axis(handles.axes4,'tight')
xlabel(handles.axes4,'Frequency [Hz]')
title(handles.axes4,title_Y)
xlim(handles.axes4,[x_min x_max])
legend(handles.axes4,'System','Validation','Location','Best')

set(handles.axes5,'XGrid','on')
set(handles.axes5,'YGrid','on')
set(handles.axes5,'XMinorGrid','on')
set(handles.axes5,'YMinorGrid','on')
axis(handles.axes5,'tight')
xlabel(handles.axes5,'Frequency [Hz]')
title(handles.axes5,title_G)
xlim(handles.axes5,[x_min x_max])
legend(handles.axes5,'Non-smoothed','Smoothed','Location','Best')

 f1 = figure(1);
    plot(real(Y_val_av(1:length(freq_plot))),imag(Y_val_av(1:length(freq_plot))),real(Y_val_est(1:length(freq_plot))),imag(Y_val_est(1:length(freq_plot))))
    legend(['RMS = ' num2str(RMS)])
    axis tight
    grid on


% --- Executes on button press in push_bode.
function push_bode_Callback(hObject, eventdata, handles)
% hObject    handle to push_bode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Bode_Mode
global U_model_av
global Y_model_av
global G_model_av
global U_val_av
global Y_val_av
global G_model_smooth_av
global freq_plot
global x_max
global x_min


if(strcmp(Bode_Mode,'Amplitude'))
    
   Bode_Mode = 'Phase';
   set(handles.push_bode,'String','Show Amplitude')
    
elseif(strcmp(Bode_Mode,'Phase'))
    
    Bode_Mode = 'Amplitude';
    set(handles.push_bode,'String','Show Phase')
   
end

%Plots
if(freq_plot(length(freq_plot))< x_max)
   x_max = freq_plot(length(freq_plot)); 
end



if(strcmp(Bode_Mode,'Amplitude'))
loglog(handles.axes3,freq_plot,abs(U_model_av(1:length(freq_plot))),freq_plot,abs(U_val_av(1:length(freq_plot))))
loglog(handles.axes4,freq_plot,abs(Y_model_av(1:length(freq_plot))),freq_plot,abs(Y_val_av(1:length(freq_plot))))
loglog(handles.axes5,freq_plot,abs(G_model_av(1:length(freq_plot))),freq_plot,abs(G_model_smooth_av(1:length(freq_plot))))
title_U = '|U(\omega)|';
title_Y = '|Y(\omega)|';
title_G = '|G(\omega)|';

else
    
semilogx(handles.axes3,freq_plot,angle(U_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(U_val_av(1:length(freq_plot)))*180/pi)
semilogx(handles.axes4,freq_plot,angle(Y_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(Y_val_av(1:length(freq_plot)))*180/pi)
semilogx(handles.axes5,freq_plot,angle(G_model_av(1:length(freq_plot)))*180/pi,freq_plot,angle(G_model_smooth_av(1:length(freq_plot)))*180/pi)
title_U = '\angle U(\omega)';
title_Y = '\angle Y(\omega)';
title_G = '\angle G(\omega)';

end

set(handles.axes1,'XGrid','on')
set(handles.axes1,'YGrid','on')
set(handles.axes1,'XMinorGrid','on')
set(handles.axes1,'YMinorGrid','on')
axis(handles.axes1,'tight')
xlabel(handles.axes1,'Time [s]')
title(handles.axes1,'u(t)')

set(handles.axes2,'XGrid','on')
set(handles.axes2,'YGrid','on')
set(handles.axes2,'XMinorGrid','on')
set(handles.axes2,'YMinorGrid','on')
axis(handles.axes2,'tight')
xlabel(handles.axes2,'Time [s]')
title(handles.axes2,'y(t)')

set(handles.axes3,'XGrid','on')
set(handles.axes3,'YGrid','on')
set(handles.axes3,'XMinorGrid','on')
set(handles.axes3,'YMinorGrid','on')
axis(handles.axes3,'tight')
xlabel(handles.axes3,'Frequency [Hz]')
title(handles.axes3,title_U)
xlim(handles.axes3,[x_min x_max])
legend(handles.axes3,'System','Validation','Location','Best')

set(handles.axes4,'XGrid','on')
set(handles.axes4,'YGrid','on')
set(handles.axes4,'XMinorGrid','on')
set(handles.axes4,'YMinorGrid','on')
axis(handles.axes4,'tight')
xlabel(handles.axes4,'Frequency [Hz]')
title(handles.axes4,title_Y)
xlim(handles.axes4,[x_min x_max])
legend(handles.axes4,'System','Validation','Location','Best')

set(handles.axes5,'XGrid','on')
set(handles.axes5,'YGrid','on')
set(handles.axes5,'XMinorGrid','on')
set(handles.axes5,'YMinorGrid','on')
axis(handles.axes5,'tight')
xlabel(handles.axes5,'Frequency [Hz]')
title(handles.axes5,title_G)
xlim(handles.axes5,[x_min x_max])
legend(handles.axes5,'Non-smoothed','Smoothed','Location','Best')
