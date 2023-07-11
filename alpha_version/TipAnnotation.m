function varargout = TipAnnotation(varargin)
% TipAnnotation MATLAB code for TipAnnotation.fig
%      TipAnnotation, by itself, creates a new TipAnnotation or raises the existing
%      singleton*.
%
%      H = TipAnnotation returns the handle to a new TipAnnotation or the handle to
%      the existing singleton*.
%
%      TipAnnotation('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TipAnnotation.M with the given input arguments.
%
%      TipAnnotation('Property','Value',...) creates a new TipAnnotation or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TipAnnotation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TipAnnotation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TipAnnotation

% Last Modified by GUIDE v2.5 08-May-2021 15:49:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TipAnnotation_OpeningFcn, ...
                   'gui_OutputFcn',  @TipAnnotation_OutputFcn, ...
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


% --- Executes just before TipAnnotation is made visible.
function TipAnnotation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TipAnnotation (see VARARGIN)

% Choose default command line output for TipAnnotation
handles.output = hObject;
handles.ZoomMode = 0;
handles.PanMode = 0;
handles.traj = []; % x y z t
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TipAnnotation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TipAnnotation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\stabilization

% import button
    h = guidata(hObject);
    file_path = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\';
    data = load([file_path 'stabilization_c2\stabilization_data.mat']);
    data = data.data;
    data = data/255;
    
    trajs = dir([file_path 'gt']);
    gt_img = zeros(size(data));
    for ii = 3:length(trajs)
        traj = load([file_path 'gt\' trajs(ii).name]);
        traj = traj.traj;
        for jj = 1:size(traj,1)
            gt_img(traj(jj,1),traj(jj,2),traj(jj,3),traj(jj,4)) = 1;
        end
    end
    gt_img = imdilate(gt_img,strel('sphere',2));
    h.gt_img = gt_img;
    
    trajs = load([file_path 'tip_info.mat']);
    exp_img = zeros(size(data));
    for tt = 1:35
        for ii = 1:length(trajs.xCoord{tt})
            exp_img(trajs.xCoord{tt}(ii), trajs.yCoord{tt}(ii), trajs.zCoord{tt}(ii), tt) = 1;
        end
    end
    exp_img = imdilate(exp_img,strel('sphere',2));
    h.exp_img = exp_img;
    
    
    h.data = data;
    h.currentT = 1;
    h.currentZ = 1;
    h.currentX = 151;
    h.currentY = 151;
    h.traj = [];
    
    axes(handles.axes1);
    h.axes1Img = max(data(:,:,:,h.currentT),[],3);
    h.axes1Img = addTip_red(h.axes1Img,max(gt_img(:,:,:,h.currentT),[],3));
    h.axes1Img = addTip_green(h.axes1Img,max(exp_img(:,:,:,h.currentT),[],3));
    imshow(h.axes1Img);
    axes(handles.axes2);
    h.axes2Img = data(:,:,h.currentZ,h.currentT);
    imshow(h.axes2Img);
    linkaxes([handles.axes1 handles.axes2]);
    
    axes(handles.axes3);
    h.axes3Img = reshape(data(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = h.axes3Img';
    img_temp = reshape(gt_img(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = addTip_red(h.axes3Img,img_temp');
    img_temp = reshape(exp_img(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = addTip_green(h.axes3Img,img_temp');
    
    imshow(h.axes3Img);
    linkaxes([handles.axes1 handles.axes3],'x');
    
    axes(handles.axes4);
    h.axes4Img = reshape(data(:,h.currentY,:,h.currentT),301,41);
    h.axes4Img = addTip_red(h.axes4Img,reshape(gt_img(:,h.currentY,:,h.currentT),301,41));
    h.axes4Img = addTip_green(h.axes4Img,reshape(exp_img(:,h.currentY,:,h.currentT),301,41));
    imshow(h.axes4Img);
    linkaxes([handles.axes1 handles.axes4],'y');

    guidata(hObject,h);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    h = guidata(hObject);
    data = h.data;
    gt_img = h.gt_img;
    exp_img = h.exp_img;
    zValue = 42 - round(get(hObject,'Value'));  % z+1 = 42
    axes(handles.axes2);
    xLimits = h.axes2.XLim;
    yLimits = h.axes2.YLim;
    h.axes2Img = data(:,:,zValue,h.currentT);
    h.axes2Img = addTip_red(h.axes2Img,gt_img(:,:,zValue,h.currentT));
    h.axes2Img = addTip_green(h.axes2Img,exp_img(:,:,zValue,h.currentT));
    imshow(h.axes2Img);
    h.axes2.XLim = xLimits;
    h.axes2.YLim = yLimits;
    h.currentZ = zValue;
    guidata(hObject,h);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    h = guidata(hObject);
    data = h.data;
    gt_img = h.gt_img;
    exp_img = h.exp_img;
    frameValue = round(get(hObject,'Value'))
    
    axes(handles.axes1);
    xLimits = h.axes1.XLim;
    yLimits = h.axes1.YLim;
    h.axes1Img = max(data(:,:,:,frameValue),[],3);
    h.axes1Img = addTip_red(h.axes1Img,max(gt_img(:,:,:,frameValue),[],3));
    h.axes1Img = addTip_green(h.axes1Img,max(exp_img(:,:,:,frameValue),[],3));
    imshow(h.axes1Img);
    h.axes1.XLim = xLimits;
    h.axes1.YLim = yLimits;
    
    axes(handles.axes2);
    xLimits = h.axes2.XLim;
    yLimits = h.axes2.YLim;
    h.axes2Img = data(:,:,h.currentZ,frameValue);
    h.axes2Img = addTip_red(h.axes2Img,gt_img(:,:,h.currentZ,frameValue));
    h.axes2Img = addTip_green(h.axes2Img,exp_img(:,:,h.currentZ,frameValue));
    imshow(h.axes2Img);
    h.axes2.XLim = xLimits;
    h.axes2.YLim = yLimits;
    h.currentT = frameValue;
    
    axes(handles.axes3);
    xLimits = h.axes3.XLim;
    yLimits = h.axes3.YLim;
    h.axes3Img = reshape(data(h.currentX,:,:,frameValue),301,41);
    h.axes3Img = h.axes3Img';
    img_temp = reshape(gt_img(h.currentX,:,:,frameValue),301,41);
    h.axes3Img = addTip_red(h.axes3Img,img_temp'); 
    img_temp = reshape(exp_img(h.currentX,:,:,frameValue),301,41);
    h.axes3Img = addTip_green(h.axes3Img,img_temp');
    imshow(h.axes3Img);
    h.axes3.XLim = xLimits;
    h.axes3.YLim = yLimits;

    axes(handles.axes4);
    xLimits = h.axes4.XLim;
    yLimits = h.axes4.YLim;
    h.axes4Img = reshape(data(:,h.currentY,:,frameValue),301,41);
    h.axes4Img = addTip_red(h.axes4Img,reshape(gt_img(:,h.currentY,:,frameValue),301,41));
    h.axes4Img = addTip_green(h.axes4Img,reshape(exp_img(:,h.currentY,:,frameValue),301,41));
    imshow(h.axes4Img);
    h.axes4.XLim = xLimits;
    h.axes4.YLim = yLimits;
    guidata(hObject,h);
    


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    h = guidata(hObject);
    if h.ZoomMode == 0
        zoom on;
        h.ZoomMode = 1;
        if h.PanMode == 1
            pan off;
            h.PanMode = 0;
        end
    elseif h.ZoomMode == 1
        zoom off;
        h.ZoomMode = 0;
    end
    guidata(hObject,h);
        
    


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    h = guidata(hObject);
    if h.PanMode == 0
        pan on;
        h.PanMode = 1;
        if h.ZoomMode == 1
            zoom off;
            h.ZoomMode = 0;
        end
    elseif h.PanMode == 1
        pan off;
        h.PanMode = 0;
    end
    guidata(hObject,h);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% pick button
    h = guidata(hObject);
    data = h.data;
    gt_img = h.gt_img;
    exp_img = h.exp_img;
    [yValue,xValue,button] = ginput(1);
    if button == 1
        xValue = round(xValue);
        yValue = round(yValue);
        h.currentX = xValue;
        h.currentY = yValue;
        h.traj = [h.traj; xValue yValue h.currentZ h.currentT];
        h.traj

        axes(handles.axes1);
        xLimits = h.axes1.XLim;
        yLimits = h.axes1.YLim;
        h.axes1Img = max(data(:,:,:,h.currentT),[],3);
        h.axes1Img = addTip_red(h.axes1Img,max(gt_img(:,:,:,h.currentT),[],3));
        h.axes1Img = addTip_green(h.axes1Img,max(exp_img(:,:,:,h.currentT),[],3));
        h.axes1Img = plotLine(h.axes1Img,xValue,yValue);
        imshow(h.axes1Img);
        h.axes1.XLim = xLimits;
        h.axes1.YLim = yLimits;

        axes(handles.axes2);
        xLimits = h.axes2.XLim;
        yLimits = h.axes2.YLim;
        h.axes2Img = data(:,:,h.currentZ,h.currentT);
        h.axes2Img = addTip_red(h.axes2Img,gt_img(:,:,h.currentZ,h.currentT));
        h.axes2Img = addTip_green(h.axes2Img,exp_img(:,:,h.currentZ,h.currentT));
        h.axes2Img = plotLine(h.axes2Img,xValue,yValue);
        imshow(h.axes2Img);
        h.axes2.XLim = xLimits;
        h.axes2.YLim = yLimits;
        

        axes(handles.axes3);
        xLimits = h.axes3.XLim;
        yLimits = h.axes3.YLim;
        h.axes3Img = reshape(data(xValue,:,:,h.currentT),301,41);
        h.axes3Img = h.axes3Img';
        img_temp = reshape(gt_img(xValue,:,:,h.currentT),301,41);
        h.axes3Img = addTip_red(h.axes3Img,img_temp');
        img_temp = reshape(exp_img(xValue,:,:,h.currentT),301,41);
        h.axes3Img = addTip_green(h.axes3Img,img_temp');
        h.axes3Img = plotLine(h.axes3Img,h.currentZ,yValue);
        imshow(h.axes3Img);
        h.axes3.XLim = xLimits;
        h.axes3.YLim = yLimits;

        axes(handles.axes4);
        xLimits = h.axes4.XLim;
        yLimits = h.axes4.YLim;
        h.axes4Img = reshape(data(:,yValue,:,h.currentT),301,41);
        h.axes4Img = addTip_red(h.axes4Img,reshape(gt_img(:,yValue,:,h.currentT),301,41));
        h.axes4Img = addTip_green(h.axes4Img,reshape(exp_img(:,yValue,:,h.currentT),301,41));
        h.axes4Img = plotLine(h.axes4Img,xValue,h.currentZ);
        imshow(h.axes4Img);
        h.axes4.XLim = xLimits;
        h.axes4.YLim = yLimits;
        guidata(hObject,h);
    end
    
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% unpick button
    h = guidata(hObject);
    data = h.data;
    gt_img = h.gt_img;
    exp_img = h.exp_img;
    h.traj(end,:) = [];
    
    axes(handles.axes1);
    xLimits = h.axes1.XLim;
    yLimits = h.axes1.YLim;
    h.axes1Img = max(data(:,:,:,h.currentT),[],3);
    h.axes1Img = addTip_red(h.axes1Img,max(gt_img(:,:,:,h.currentT),[],3));
    h.axes1Img = addTip_green(h.axes1Img,max(exp_img(:,:,:,h.currentT),[],3));
    imshow(h.axes1Img);
    h.axes1.XLim = xLimits;
    h.axes1.YLim = yLimits;

    axes(handles.axes2);
    xLimits = h.axes2.XLim;
    yLimits = h.axes2.YLim;
    h.axes2Img = data(:,:,h.currentZ,h.currentT);
    h.axes2Img = addTip_red(h.axes2Img,gt_img(:,:,h.currentZ,h.currentT));
    h.axes2Img = addTip_green(h.axes2Img,exp_img(:,:,h.currentZ,h.currentT));
    imshow(h.axes2Img);
    h.axes2.XLim = xLimits;
    h.axes2.YLim = yLimits;


    axes(handles.axes3);
    xLimits = h.axes3.XLim;
    yLimits = h.axes3.YLim;
    h.axes3Img = reshape(data(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = h.axes3Img';
    img_temp = reshape(gt_img(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = addTip_red(h.axes3Img,img_temp');
    img_temp = reshape(exp_img(h.currentX,:,:,h.currentT),301,41);
    h.axes3Img = addTip_green(h.axes3Img,img_temp');
    
    imshow(h.axes3Img);
    h.axes3.XLim = xLimits;
    h.axes3.YLim = yLimits;

    axes(handles.axes4);
    xLimits = h.axes4.XLim;
    yLimits = h.axes4.YLim;
    h.axes4Img = reshape(data(:,h.currentY,:,h.currentT),301,41);
    h.axes4Img = addTip_red(h.axes4Img,reshape(gt_img(:,h.currentY,:,h.currentT),301,41));
    h.axes4Img = addTip_green(h.axes4Img,reshape(exp_img(:,h.currentY,:,h.currentT),301,41));
    imshow(h.axes4Img);
    h.axes4.XLim = xLimits;
    h.axes4.YLim = yLimits;
    guidata(hObject,h);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% save trajectory
    path = 'D:\dropbox\Modify Series Data\SL-092320-slice1-hippo-vessel-Modify Series\gt';
    ind = num2str(1000 + length(dir(path)) - 1);
    ind = ind(2:4);
    h = guidata(hObject);
    traj = h.traj;
    save([path '\' ind], 'traj');
% for show, no save
    h.traj = [];
    guidata(hObject,h);
    


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


function img = plotLine(img,x,y)
    if size(img,3) == 1
        img = repmat(img,1,1,3);
    end
    img(x,1:y-3,1) = 1;
    img(x,1:y-3,2:3) = 0;
    img(x,y+3:end,1) = 1;
    img(x,y+3:end,2:3) = 0;
    img(1:x-3,y,1) = 1;
    img(1:x-3,y,2:3) = 0;
    img(x+3:end,y,1) = 1;
    img(x+3:end,y,2:3) = 0;
    
function img = addTip_red(img,gt_img)
    if size(img,3) == 1
        img = repmat(img,1,1,3);
    end
    img(:,:,1) = img(:,:,1).*(1-gt_img) + gt_img;
    img(:,:,2) = img(:,:,2).*(1-gt_img);
    img(:,:,3) = img(:,:,3).*(1-gt_img);
    
    
function img = addTip_green(img,gt_img)
    if size(img,3) == 1
        img = repmat(img,1,1,3);
    end
    img(:,:,1) = img(:,:,1);
    img(:,:,2) = img(:,:,2) + gt_img/2;
    img(:,:,3) = img(:,:,3);
