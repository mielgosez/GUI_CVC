function varargout = main2(varargin)
% MAIN2 MATLAB code for main2.fig
%      MAIN2, by itself, creates a new MAIN2 or raises the existing
%      singleton*.
%
%      H = MAIN2 returns the handle to a new MAIN2 or the handle to
%      the existing singleton*.
%
%      MAIN2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN2.M with the given input arguments.
%
%      MAIN2('Property','Value',...) creates a new MAIN2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main2

% Last Modified by GUIDE v2.5 13-Apr-2017 17:32:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main2_OpeningFcn, ...
                   'gui_OutputFcn',  @main2_OutputFcn, ...
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


% --- Executes just before main2 is made visible.
function main2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main2 (see VARARGIN)
clc
scatter(0,0,100,0,'filled','b'), set(gca,'visible','off')
% Choose default command line output for main2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadGraphs.
function loadGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to loadGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.loadGraphs,'Visible','off');
[file,path] = uigetfile('*.mat','Select data file');
load([path filesep file]); % Load adjacency list
load([path filesep 'names.mat']);
load([path filesep 'iLCD' filesep 'VisCommunitiesTestIm_2.mat']); % Load communities
[~,handles.path_image] = uigetfile('*','Select the folder where images are stored');
[handles.r,handles.c] = size(AdjExt_knn);
[handles.numCommunities,~] = size(VisCommunitiesExt);
handles.AdjExt_knn = AdjExt_knn; % Store adjacency list into handles.
handles.communities = VisCommunitiesExt; %Store communities into handles.
handles.names = list_imgsTest; % Names of the images linked to the nodes.
fig = main2;
handles.texto = uicontrol(fig,'Style','text',...
    'String','Select a Community',...
    'Tag','titleText',...
    'Position',[900 580 200 15],'Visible','on');
handles.chkb0 = uicontrol(main2,'Style','checkbox',...
                  'String','No communities',...
                  'Tag','chkb0',...
                  'Value', 0,...
                  'Position',[900 (580-(0.03)*580) 120 15]);
for i = 1:handles.numCommunities
    handles.chkb{i} = uicontrol(main2,'Style','checkbox',...
                  'String',['Community ' num2str(i)],...
                  'Tag',['chkb' num2str(i)],...
                  'Value', 0,...
                  'Position',[900 (580-(i+1)*(0.03)*580) 120 15]);
end
plotComButton = uicontrol(main2,'Style','pushbutton',...
    'String','Plot',...
    'Tag','communityButton',...
    'Position',[900 100 100 15],...
    'Callback',@plotComButton_Callback);
%set(handles.popCommunity,'String',str);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotComButton_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
indexes = [];           % Listing communities.
indexes_nodes = [];     % Nodes belonging to a community.
indexes_nodes_null = [];
if(get(handles.chkb0,'Value'))
    indexes_nodes_null = find(sum(handles.communities)==0);% Nodes not belonging to any community.
end

for i = 1:handles.numCommunities
    if(get(handles.chkb{i},'Value'))
        indexes = [indexes i]; % Communities indices
        indexes_nodes = [indexes_nodes find(handles.communities(i,:) == 1)]; % Nodes' indices.
    end
end
indexes_nodes = unique(indexes_nodes);
bg1   = biograph(handles.AdjExt_knn([indexes_nodes_null indexes_nodes],[indexes_nodes_null indexes_nodes]))  ; % Biograph is instantiated.
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
%% Print nodes and edges
clear h nodes i bg1 names k child_handles; % Removes unnecessary variables
%%% Convert coordinates to canvas.
x_max = max(handles.positions(:,1));
x_min = min(handles.positions(:,1));
y_max = max(handles.positions(:,2));
y_min = min(handles.positions(:,2));
handles.positions(:,1) = 900*(handles.positions(:,1)-x_min)/(x_max-x_min)+10;
handles.positions(:,2) = 600*(handles.positions(:,2)-y_min)/(y_max-y_min)+10;
scatter(handles.positions(:,1),handles.positions(:,2),40,zeros(length([indexes_nodes_null indexes_nodes]),1),'filled','b'), set(gca,'visible','off'), % Command deletes the axis
hold on
scatter(handles.positions(1:length(indexes_nodes_null),1),handles.positions(1:length(indexes_nodes_null),2),40,zeros(length(indexes_nodes_null),1),'filled','black'), set(gca,'visible','off'), % Command deletes the axis
hold on
total_index = [indexes_nodes_null indexes_nodes];
%%
for i = 1:length(total_index)
    %text(handles.positions(i,1),handles.positions(i,2),num2str(i)),
    hold on
    %% This is a tag made to visualize the image linked to the node and provide with the necessary information (hovering).
    tagButton{i} = uicontrol(main2,'Style','pushbutton',...
                    'String','',...
                    'Tag',['nodeButton' num2str(i)],...
                    'Position',[handles.positions(i,1) handles.positions(i,2) 10 10],...
                    'TooltipString',['Community:' mat2str(find(handles.communities(:,total_index(i)) == 1)) sprintf('\n')],...
                    'Callback',{@print_info,i,handles.names,handles.path_image});
    %%
    for j = 1:(i-1)
        if handles.AdjExt_knn(total_index(i),total_index(j)) == 1
            plot(handles.positions([i j],1),handles.positions([i j],2),'black'),
            hold on
        end
    end
end
%% Print communities
adj_local = handles.AdjExt_knn([indexes_nodes_null indexes_nodes],[indexes_nodes_null indexes_nodes]); % Reduced adjacency matrix.
com_local = handles.communities(indexes,[indexes_nodes_null indexes_nodes])                          ; % Reduced community matrix.
n         = length(indexes)                                                                          ; % Number of communities.
for i = 1:length(indexes)
    indices     = find(com_local(i,:) == 1)   ; % Nodes in the community.
    coordinates = handles.positions(indices,:); % Coordinates for these nodes.
    try
        %% This covers nodes whoch do not belong to the highlighted community
        %convexHull = convhull(coordinates(:,1),coordinates(:,2));
        %p{i} = plot(coordinates(convexHull,1),coordinates(convexHull,2));
        %% This new approach just highlight the edges.
        
        for j=1:length(indices)
            local_index = find(adj_local(indices(j),indices) == 1); % Adjacent nodes in the community.
            for k = 1:length(local_index)
                p{j,k} = plot([handles.positions(indices(j),1) handles.positions(indices(local_index(k)),1)],...
                              [handles.positions(indices(j),2) handles.positions(indices(local_index(k)),2)]); % Print all the edges in a community.
                set(p{j,k},'Color',[(i-1)/(n-1) 0 1-i/n]), % Set a particular color for this community.
                set(p{j,k},'linewidth',2);
                hold on;
            end
        end
        x_max = max(handles.positions(:,1)); %
        y_max = max(handles.positions(:,2));
        scatter(0.8*x_max,0.8*y_max-0.05*y_max*i,40,[(i-1)/(n-1) 0 1-i/n],'fill'); %Legend
        hold on,
        text(0.85*x_max,0.8*y_max-0.05*y_max*i,['Community ' num2str(i)]); % Legend
        hold on;
    catch ME % Sometimes there are collinear points.
        disp(ME);
        disp([(i-1)/(n-1) 0 1-i/n]);
        scatter(coordinates(:,1),coordinates(:,2),50,[(i-1)/(n-1) 0 1-i/n],'fill'), set(gca,'visible','off'), % Command deletes the axis
        hold on
        x_max = max(handles.positions(:,1)); %
        y_max = max(handles.positions(:,2));
        scatter(0.8*x_max,0.8*y_max-0.05*y_max*i,40,[(i-1)/(n-1) 0 1-i/n],'fill');
        hold on,
        text(0.85*x_max,0.8*y_max-0.05*y_max*i,['Community ' num2str(i)]);
        hold on;
    end
end
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

function print_info(src,event,nodo,names,path_image)
    disp(names(nodo).nm);
    figure(1),imshow([path_image filesep names(nodo).nm '.jpg']),
    title(['Image ' handles.names(nodo).nm '.jpg'])
