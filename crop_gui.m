function varargout = crop_gui(varargin)
% CROP_GUI MATLAB code for crop_gui.fig
%      CROP_GUI, by itself, creates a new CROP_GUI or raises the existing
%      singleton*.
%
%      H = CROP_GUI returns the handle to a new CROP_GUI or the handle to
%      the existing singleton*.
%
%      CROP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROP_GUI.M with the given input arguments.
%
%      CROP_GUI('Property','Value',...) creates a new CROP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crop_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crop_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crop_gui

% Last Modified by GUIDE v2.5 29-Oct-2013 14:51:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crop_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @crop_gui_OutputFcn, ...
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


% --- Executes just before crop_gui is made visible.
function crop_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crop_gui (see VARARGIN)

% Choose default command line output for crop_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



global Data
global Data_Length

Data_Matrix = load('Data_Matrix_thin.mat');
Data_Matrix_Length = load('Data_Matrix_Length_thin.mat');
Data = Data_Matrix.Data_Matrix;
Data_Length = Data_Matrix_Length.Data_Matrix_Length;


[size_m size_n] = size(Data);
Number_cell = {1 : round(size_n/2)};
set(handles.listbox1,'String',Number_cell{1,1}(:))

global k
global t
global Fs
global max1
global max2
global min1
global min2

k = get(handles.listbox1,'Value')
global current_data1
global current_data2
Fs = 44.1e3/4
t = (0:1:Data_Length(k)-1)/Fs;
current_data1 = Data(1:Data_Length(k),2*k-1);
current_data2 = Data(1:Data_Length(k),2*k);

    plot(handles.axes1,t,current_data1)
    plot(handles.axes2,t,current_data2)
    max1 = max(current_data1);
    min1 = min(current_data1);
    max2 = max(current_data2);
    min2 = min(current_data2);



% UIWAIT makes crop_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crop_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k
global t
global Fs
global max1
global max2
global min1
global min2

k = get(hObject,'Value')
global current_data1
global current_data2
global Data
global Data_Length

Fs = 44.1e3/4
t = (0:1:Data_Length(k)-1)/Fs;
current_data1 = Data(1:Data_Length(k),2*k-1);
current_data2 = Data(1:Data_Length(k),2*k);

    plot(handles.axes1,t,current_data1)
    plot(handles.axes2,t,current_data2)
    max1 = max(current_data1);
    min1 = min(current_data1);
    max2 = max(current_data2);
    min2 = min(current_data2);
    
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

global Data
global Data_Length
global k
global I_beg
global I_end
% 
% Data(I_end:length(Data),k*2-1) = 0;
% Data(1:I_beg,k*2-1) = 0;


% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k
global t
global Data_Length
global current_data1
global current_data2
global max1
global max2
global min1
global min2
global y1
global y2
global x_end
global x_beg
global I_end

position = get(hObject,'Value');
I_end = round(position*Data_Length(k));
if(I_end ==0)
   I_end = 1; 
end

x_beg  = [t(I_end),t(I_end)];
y1 = [min1,max1];
y2 = [min2,max2];

plot(handles.axes1,t,current_data1, x_beg, y1, 'g', x_end, y1, 'r')
plot(handles.axes2,t,current_data2, x_beg, y2, 'g', x_end, y2, 'r')



% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global k
global t
global Data_Length
global current_data1
global current_data2
global max1
global max2
global min1
global min2
global y1
global y2
global x_end
global x_beg
global I_end


position = get(hObject,'Value');
I_end = round(position*Data_Length(k));

if(I_end ==0)
   I_end = 1; 
end

x_end  = [t(I_end),t(I_end)];
y1 = [min1,max1];
y2 = [min2,max2];

plot(handles.axes1,t,current_data1, x_beg, y1, 'g', x_end, y1, 'r')
plot(handles.axes2,t,current_data2, x_beg, y2, 'g', x_end, y2, 'r')

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
