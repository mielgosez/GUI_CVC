function varargout = main3(varargin)
% MAIN3 MATLAB code for main3.fig
%      MAIN3, by itself, creates a new MAIN3 or raises the existing
%      singleton*.
%
%      H = MAIN3 returns the handle to a new MAIN3 or the handle to
%      the existing singleton*.
%
%      MAIN3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN3.M with the given input arguments.
%
%      MAIN3('Property','Value',...) creates a new MAIN3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main3

% Last Modified by GUIDE v2.5 17-Apr-2017 18:49:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main3_OpeningFcn, ...
                   'gui_OutputFcn',  @main3_OutputFcn, ...
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


% --- Executes just before main3 is made visible.
function main3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main3 (see VARARGIN)
clc;
plot([],[]), set(gca,'visible','off');
% Choose default command line output for main3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plotGraph.
function plotGraph_Callback(hObject, eventdata, handles)
% hObject    handle to loadGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Switch off the previous button.
set(handles.plotGraph,'Visible','off');
%% Select all the files needed to create the graph.
[file,path] = uigetfile('*.mat','Select data file');
load([path filesep file]); % Load adjacency list
%load([path filesep 'names.mat']);
%load([path filesep 'iLCD' filesep 'VisCommunitiesTestIm_2.mat']); % Load communities
[~,handles.path_image] = uigetfile('*','Select the folder where images are stored');
%% Outliers deletion.
index_without_outliers = find( sum(AdjExt_knn) > 0);
AdjExt_knn = AdjExt_knn(index_without_outliers,index_without_outliers);
VisCommunitiesExt = VisCommunitiesExt(:,index_without_outliers);
list_imgsTrainExt = list_imgsTrainExt(index_without_outliers);
%% Feed the handles with information worth to keep later.
[handles.r,handles.c] = size(AdjExt_knn);
[handles.numCommunities,~] = size(VisCommunitiesExt);
handles.AdjExt_knn = AdjExt_knn        ; % Store adjacency list into handles.
handles.communities = VisCommunitiesExt; % Store communities into handles.
handles.names = list_imgsTrainExt      ; % Names of the images linked to the nodes.
%% Use the biograph to obtain the positions.
bg1   = biograph(handles.AdjExt_knn)      ; % Biograph is instantiated.
h     = view(bg1)                         ; % Show graph => set position of nodes.
nodes = get(h,'Nodes')                    ; % Get all nodes' info.
handles.positions = zeros(length(nodes),2);
for i=1:length(nodes)
    handles.positions(i,:) = nodes(i).Position;
end
%% All this code is to shut down a figure window on biograph.
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))
%% Plot all the nodes first.
clear h nodes i bg1 names k child_handles; % Removes unnecessary variables
p1 = scatter(handles.positions(:,1),handles.positions(:,2),40,zeros(handles.r,1),'filled','b'), set(gca,'visible','off'), % Command deletes the axis
hold on 
%% Plot all the edges then.
for i = 1:handles.r
    for j = 1:(i-1)
        if AdjExt_knn(i,j) == 1
            p3{i} = plot(handles.positions([i j],1),handles.positions([i j],2),'black'),
            hold on
        end
    end
end
%% Find nodes which do not belong to any community and print them in black.
indexes_nodes_null = find(sum(handles.communities)==0);% Nodes not belonging to any community.
p2 = scatter(handles.positions(indexes_nodes_null,1),...
        handles.positions(indexes_nodes_null,2),...
        50,...
        zeros(length(indexes_nodes_null),1),...
        'filled','black'), set(gca,'visible','off'), % Command deletes the axis
%% Create buttons such that a isolated node is selected and the communities connected to it can be inspected.
max_x = max(handles.positions(:,1));
min_x = min(handles.positions(:,1));
max_y = max(handles.positions(:,2));
min_y = min(handles.positions(:,2));
canvas_x = 951;
canvas_y = 621;
handles.positions2      = zeros(size(handles.positions));
handles.positions2(:,1) = canvas_x*(handles.positions(:,1)-min_x)/(max_x-min_x);
handles.positions2(:,2) = canvas_y*(handles.positions(:,2)-min_y)/(max_y-min_y);

for i = 1:length(indexes_nodes_null)
    tagButton{i} = uicontrol(main3,'Style','pushbutton' ,...
                    'String'  ,''                       ,...
                    'Tag'     ,['nodeButton' num2str(i)],...
                    'Position',[handles.positions2(indexes_nodes_null(i),1) handles.positions2(indexes_nodes_null(i),2) 10 10],...
                    'Callback',{@print_comm,...
                      indexes_nodes_null(i),... % Index of the non connected node.
                              handles.names,...
                         handles.path_image,...
                        handles.communities,...
                         handles.AdjExt_knn,...
                         indexes_nodes_null,...
                         });
end
%set(handles.popCommunity,'String',str);
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_comm(src,event, id,name,path_image,communities,AdjExt_knn,indexes_nodes_null)
%% Determine the nodes to be plotted
linked_nodes = find(AdjExt_knn(id,:) == 1); % Nodes associated to the selected not integrated node.
smaller_matrix = communities(:,linked_nodes); % This is a matrix smaller than the previous one, because it only takes into account the chosen nodes to verify communities.
smaller_matrix = smaller_matrix'; % We transpose the matrix so that when sum operation is performed, the results could be interpreted as the number of elements on each community.
% if the number of elements in a community of 'smaller_matrix' is bigger
% than one, then it means that the target node is somehow connected to it.
disp(smaller_matrix)
if(size(smaller_matrix,1) < 2)
    linked_communities = find(smaller_matrix>=1) ;
else
    linked_communities = find(sum(smaller_matrix)>=1)   ;           % Communities the target is connected to.
end
total_nodes_comm = []                                         ; %All the nodes linked to the  target and belonging to at least one communities.
total_nodes_void = find(sum(communities(:,linked_nodes)) == 0); %All the nodes linked to the target but not belonging to any community.
for j = 1:length(linked_communities)
    total_nodes_comm = [total_nodes_comm find(communities(linked_communities(j),:) == 1)];
end
total_nodes_comm = unique(total_nodes_comm);   
all_nodes = [id linked_nodes(total_nodes_void) total_nodes_comm];% Merging both classes of nodes. 
disp('All nodes');
disp(all_nodes);
%disp('AdjExt_knn');
%disp(AdjExt_knn(all_nodes,all_nodes));
disp('id');
disp(id);

%% Visualization with biograph
disp('Adjacency matrix')
disp(AdjExt_knn(all_nodes,all_nodes))
disp('Connectedness of node')
disp(AdjExt_knn(id,:))
bg1   = biograph(AdjExt_knn(all_nodes,all_nodes)); % Biograph is instantiated.
h     = view(bg1)                                ; % Show graph => set position of nodes.
nodes = get(h,'Nodes')                           ; % Get all nodes' info.
positions = zeros(length(nodes),2);
for j=1:length(nodes)
   positions(j,:) = nodes(j).Position;
end
%% All this code is to shut down a figure window on biograph.
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))
clear h nodes i bg1 names k child_handles; % Removes unnecessary variables
%% Plot the edges not belonging to any community.
for i = 1:length(all_nodes)
    for j = 1:(i-1)
        if AdjExt_knn(all_nodes(i),all_nodes(j)) == 1
            figure(1000),
            subplot(1,2,1),
            plot(positions([i j],1),positions([i j],2),'black'), % Remember that this is a row of x's and y's.
            hold on
        end
    end
end
%% Colouring the edges.
n         = length(linked_communities); 
for i = 1:length(linked_communities)
    indices     = find(communities(linked_communities(i),all_nodes) == 1)   ; % Nodes in the community. NOTE THAT THIS MATRIX IS REDUCED TO JUST CONSIDER RELEVANT NODES.
    coordinates = positions(indices,:); % Coordinates for these nodes.
    for j=1:length(indices)
        %disp(adj_local(all_nodes(indices(j)),all_nodes(indices)));
        local_index = find(AdjExt_knn(all_nodes(indices(j)),all_nodes(indices)) == 1); % Adjacent nodes in the community.
            for k = 1:length(local_index)
                figure(1000),
                subplot(1,2,1)
                p{j,k} = plot([positions(indices(j),1) positions(indices(local_index(k)),1)],...
                              [positions(indices(j),2) positions(indices(local_index(k)),2)]); % Print all the edges in a community.
                set(p{j,k},'Color',[(i-1)/(n) 0 1-i/n]), % Set a particular color for this community.
                set(p{j,k},'linewidth',2);
                hold on;
            end
    end
    %scatter(0.9*x_max,0.8*y_max-0.05*y_max*i,40,[(i-1)/(n-1) 0 1-i/n],'fill'); %Legend
    %hold on,
    %text(x_max,0.8*y_max-0.05*y_max*i,['Community ' num2str(i)]); % Legend
    %hold on;
    figure(1000),
    subplot(1,2,2),
    scatter(4,9-0.5*i,40,[(i-1)/(n) 0 1-i/n],'fill'), %Legend
    set(gca,'visible','off'), % Command deletes the axis
    xlim([0 10]), ylim([0 10]),
    hold on,
    figure(1000),
    subplot(1,2,2),
    text(5,9-0.5*i,['Community ' num2str(linked_communities(i))]); % Legend
    set(gca,'visible','off'), % Command deletes the axis
    xlim([0 10]), ylim([0 10]),
    hold on;
end
%% Plot all the nodes first.
figure(1000),
subplot(1,2,1),
p1 = scatter(positions(:,1),...
             positions(:,2),...
             60,...
             zeros(length(positions(:,1)),1),...
             'filled','b');
set(gca,'visible','off'), % Command deletes the axis
hold on 
figure(1000),
subplot(1,2,1),
p2 = scatter(positions(2:(length(total_nodes_void)+1),1),...
             positions(2:(length(total_nodes_void)+1),2),...
             60,...
             zeros(length(total_nodes_void),1),...
             'filled',...
             'black'); 
set(gca,'visible','off'), % Command deletes the axis
hold on 
figure(1000),
subplot(1,2,1),
p3 = scatter(positions(1,1),...
             positions(1,2),...
             100,...
             0,...
             'filled',...
             'green'); 
set(gca,'visible','off'), % Command deletes the axis
hold on 
%% Show the images linked to the selected nodes.   
for i = 1:length(all_nodes)    
    %tagButton{i} = uicontrol(figure(2000),'Style','pushbutton',...
    %                'String','',...
    %                'Tag',['nodeButton' num2str(i)],...
    %                'Units','normalized',...
    %                'Position',[0.9*(positions(i,1)-x_min)/(x_max-x_min) 0.9*(positions(i,2)-y_min)/(y_max-y_min) .05 .05],...
    %                'Callback',{@print_info,all_nodes(i),name,path_image});
    disp([path_image name(all_nodes(i)).nm '.jpg'])
    im = imread([path_image filesep name(all_nodes(i)).nm '.jpg']);
    figure(2000),
    plot3_images_in_coordinates(im,[0 positions(i,1) positions(i,2)]),
    set(gca,'visible','off'),
    hold on
    for j = 1:(i-1)
            if AdjExt_knn(all_nodes(i),all_nodes(j)) == 1
                figure(2000),
                plot3([0 0],positions([i j],1),positions([i j],2),'black'),
                set(gca,'visible','off'), % Command deletes the axis
                hold on
            end
    end    
end
%% Plotting the metrics for the communities and  selected node
% representativity
%id = 20; %Node
neigh_i = find(AdjExt_knn(id,:) == 1); %Neighborhood of i
Ci = find(sum(communities(:,neigh_i),2) > 0); %Communities
representativity = zeros(length(Ci),1);
for j = 1:length(Ci)
    for i = 1:length(neigh_i)
        representativity(j) = representativity(j) + communities(Ci(j),neigh_i(i));
    end
    representativity(j) = representativity(j)/length(neigh_i);
end
% Intrinsic cohesion
IC = zeros(length(Ci),1); % Intrinsic cohesion will be stored here per each community.
for i = 1:length(Ci)
    ciElements = find(communities(Ci(i),:) == 1); % Elements in community i.
    representativity_com = zeros(length(ciElements),1); % Representativity will be stored here
    for j=1:length(ciElements)
        neigh_j = find(AdjExt_knn(ciElements(j),:) == 1); %Neighbor
        for k = 1:length(neigh_j)
            representativity_com(j) = representativity_com(j) + communities(Ci(i),neigh_j(k));
        end
        representativity_com(j) = representativity_com(j)/length(neigh_j);
    end
    IC(i) = sum(representativity_com);
end
% Afilliation force
neighElements = find(AdjExt_knn(id,:) == 1); % Neighbors of the node id.
FA = zeros(length(Ci),1); % Affiliation force initialization.
for i = 1:length(Ci)
    ciElements = find(communities(Ci(i),:) == 1); % Elements in community i.
    representativity_com = zeros(length(ciElements),1); % Representativity will be stored here
    for j=1:length(ciElements)
        for k = 1:length(neighElements)
            representativity_com(j) = representativity_com(j) + communities(Ci(i),neighElements(k));
        end
        representativity_com(j) = representativity_com(i)/length(neigh_i);
    end
    FA(i) = sum(representativity_com);
end
figure(3000),
subplot(2,2,1),bar(Ci,representativity),title('Representativity vs Community'),
subplot(2,2,2),bar(Ci,IC),title('Intrinsec Cohesion vs Community'),
subplot(2,2,3),bar(Ci,FA),title('Affiliation Force vs Community'),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_info(src,event,nodo,names,path_image)
    figure(1),imshow([path_image filesep names(nodo).nm '.jpg']),
    disp('image');disp(nodo);
    title(['Image ' names(nodo).nm '.jpg'])