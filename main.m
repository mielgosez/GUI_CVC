function varargout = main(varargin)
% MAIN MATLAB code for main.fig
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
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 10-Apr-2017 20:25:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)
%clc
scatter(0,0,100,0,'filled','b'), set(gca,'visible','off')
% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popCommunity.
function popCommunity_Callback(hObject, eventdata, handles)
% hObject    handle to popCommunity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popCommunity contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popCommunity


% --- Executes during object creation, after setting all properties.
function popCommunity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCommunity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in plotCommunity.
function plotCommunity_Callback(hObject, eventdata, handles)
% hObject    handle to plotCommunity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(handles.popCommunity,'String'); 
popChoice = contents{get(handles.popCommunity,'Value')};
%contents = cellstr( get(hObject,'String') )   ; % Collect all the values.
%popChoice = contents(get(hObject,'Value'))    ; % Detect the choice made by the client.
n = handles.numCommunities;
k = str2num(popChoice); % Choice made by the client.
indices = find(handles.communities(k,:) == 1);
coordinates = handles.positions(indices,:);
try
    convexHull = convhull(coordinates(:,1),coordinates(:,2));
    fill(coordinates(convexHull,1),coordinates(convexHull,2),[(k-1)/(n-1) 0 1-k/n],'facealpha',0.1),
    hold on;
    x_max = max(handles.positions(:,1)); %
    y_max = max(handles.positions(:,2));
    
    scatter(0.8*x_max,0.8*y_max-0.05*y_max*k,40,[(k-1)/(n-1) 0 1-k/n],'fill');
    hold on,
    text(0.85*x_max,0.8*y_max-0.05*y_max*k,['Community ' num2str(k)]);
    hold on;
catch ME % Sometimes there are collinear points.
    scatter(coordinates(:,1),coordinates(:,2),50,[(k-1)/(n-1) 0 1-k/n],'fill'), set(gca,'visible','off'), % Command deletes the axis
    hold on 
    x_max = max(handles.positions(:,1)); %
    y_max = max(handles.positions(:,2));
    scatter(0.8*x_max,0.8*y_max-0.05*y_max*k,40,[(k-1)/(n-1) 0 1-k/n],'fill');
    hold on,
    text(0.85*x_max,0.8*y_max-0.05*y_max*k,['Community ' num2str(k)]);
    hold on;
end
set(handles.tagCommunity,'String',['You have chosen the community ' num2str(popChoice)]);
%handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in clearGraph.
function clearGraph_Callback(hObject, eventdata, handles)
% hObject    handle to clearGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadGraph.
function loadGraph_Callback(hObject, eventdata, handles)
% hObject    handle to loadGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile();
load([path filesep file]); % Load adjacency list
load([path filesep 'iLCD' filesep 'VisCommunitiesTestIm_2.mat']); % Load communities
bg1   = biograph(AdjExt_knn)                                       ; % Biograph is instantiated.
h     = view(bg1)                                          ; % Show graph => set position of nodes.
nodes = get(h,'Nodes')                                     ; % Get all nodes' info.
handles.positions = zeros(length(nodes),2);
for i=1:length(nodes)
    handles.positions(i,:) = nodes(i).Position;
end
%% All this code is to shut down a figure window on biograph.
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))
%%
clear h nodes i bg1 names k child_handles; % Removes unnecessary variables
[handles.r,handles.c] = size(AdjExt_knn);
scatter(handles.positions(:,1),handles.positions(:,2),20,zeros(handles.r,1),'filled','b'), set(gca,'visible','off'), % Command deletes the axis
hold on 

handles.AdjExt_knn = AdjExt_knn; % Store adjacency list into handles.
handles.communities = VisCommunitiesExt; %Store communities into handles.
for i = 1:handles.r
    %text(handles.positions(i,1),handles.positions(i,2),num2str(i)),
    hold on
    for j = 1:(i-1)
        if AdjExt_knn(i,j) == 1
            plot(handles.positions([i j],1),handles.positions([i j],2),'black'),
            hold on
        end
    end
end
[handles.numCommunities,~] = size(VisCommunitiesExt);
set(handles.tagCommunity,'String',num2str(handles.numCommunities));
%str = '';
%for i = 1:numCommunities
%    if(i ~= numCommunities)
%        str = sprintf('%s\n%s',str,num2str(i));
%    else
%        str = sprintf('%s\n%s',str,num2str(numCommunities));
%    end
%end
n = handles.numCommunities;
%myarray = strcat('Data',cellfun(@num2str,num2cell((1:n)'),'UniformOutput',false));
myarray = strcat(cellfun(@num2str,num2cell((1:n)'),'UniformOutput',false));
set(handles.popCommunity,'String',myarray,'Value',1);
set(handles.popCommunity,'Enable','On');
%set(handles.popCommunity,'String',str);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
