function varargout = GUI_kMeans_2b(varargin)
% GUI_KMEANS_2B MATLAB code for GUI_kMeans_2b.fig
%      GUI_KMEANS_2B, by itself, creates a new GUI_KMEANS_2B or raises the existing
%      singleton*.
%
%      H = GUI_KMEANS_2B returns the handle to a new GUI_KMEANS_2B or the handle to
%      the existing singleton*.
%
%      GUI_KMEANS_2B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_KMEANS_2B.M with the given input arguments.
%
%      GUI_KMEANS_2B('Property','Value',...) creates a new GUI_KMEANS_2B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_kMeans_2b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_kMeans_2b_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_kMeans_2b

% Last Modified by GUIDE v2.5 23-Mar-2020 11:06:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_kMeans_2b_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_kMeans_2b_OutputFcn, ...
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


% --- Executes just before GUI_kMeans_2b is made visible.
function GUI_kMeans_2b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_kMeans_2b (see VARARGIN)

% Choose default command line output for GUI_kMeans_2b
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_kMeans_2b wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_kMeans_2b_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Dataset.
function Dataset_Callback(hObject, eventdata, handles)
% hObject    handle to Dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Dataset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Dataset


% --- Executes during object creation, after setting all properties.
function Dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        feature_dir = [pwd '/lamdaDataset/features/'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all;
handles.targets = targets;
handles.feature_dir = feature_dir;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on selection change in QueryName1.
function QueryName1_Callback(hObject, eventdata, handles)
% hObject    handle to QueryName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns QueryName1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from QueryName1

imindex = get(hObject,'Value');
image_dir = handles.image_dir;
fname = [image_dir handles.filenames{imindex}];
axes(handles.axes1);
imshow(imread(fname)); axis image;
handles.q1Idx = imindex;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function QueryName1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QueryName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in QueryName2.
function QueryName2_Callback(hObject, eventdata, handles)
% hObject    handle to QueryName2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns QueryName2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from QueryName2
imindex = get(hObject,'Value');
image_dir = handles.image_dir;
fname = [image_dir handles.filenames{imindex}];
axes(handles.axes2);
imshow(imread(fname)); axis image;
handles.q2Idx = imindex;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function QueryName2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QueryName2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Create_pp_run_kMeans.
function Create_pp_run_kMeans_Callback(hObject, eventdata, handles)
% hObject    handle to Create_pp_run_kMeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes3,'reset');
num_Clusters = str2double(get(handles.number_of_Clusters, 'String'));

data = handles.data;
features = handles.features;
filenames = handles.filenames;   
targets = handles.targets;
    
queryIndex1 = handles.q1Idx;
queryIndex2 = handles.q2Idx; 
    
q_1 = data(queryIndex1,:); 
q_2 = data(queryIndex2,:); 
%q3 = data(queryIndex3,:); 

N = length(filenames); 
q1new = repmat(q_1,N,1);
q2new = repmat(q_2,N,1);
%q3new = repmat(q3,N,1);

dist_1 = xor(data, q1new);
dist_2 = xor(data, q2new);
%dist_3 = xor(data, q2new);

hamming_dist1 = sum(dist_1,2);
hamming_dist2 = sum(dist_2,2);
%hamming_dist3 = sum(dist_3,2);

n_hamming_dist1 = mat2gray(hamming_dist1);
n_hamming_dist2 = mat2gray(hamming_dist2);
%n_hamming_dist3 = mat2gray(hamming_dist3);
 
 
X = zeros(2,N);
X(1,:) = hamming_dist1;
X(2,:) = hamming_dist2;
%X(3,:) = n_hamming_dist3;
X = (X)';

 axes(handles.axes3);
 hold off; plot(X(:,1),X(:,2),'.');
 hold on; 
    


[K,L] = size(unique(X,'rows'));  %% Number of unique pareto points 
set(handles.num_pr_po,'String',num2str(K))

    
[idx,C] = kmeans(X,num_Clusters) ;  % Applay kmeans, C:Coordinates of the Clusters, D:Whinhin distance in a cluster 

% Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
DATA= cell(num_Clusters,1) ;

axes(handles.axes3);

for i = 1:num_Clusters
    DATA{i} = X(idx==i, :) ;    
    plot(X(idx==i, 1), X(idx==i, 2)  ,'.')  
        
end



    VAR = var(C,0,2);
    [VAR_sort , VAR_index] = sort(VAR);
    VAR_index = VAR_index(1:10, : );
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = [VAR_index, SUM ];
    LAST_sort = sortrows(LAST, 2);
    Best = LAST_sort(1,1);
      
    
    axes(handles.axes3);
    % Plot the best one 
    plot(DATA{Best,1}(:,1) , DATA{Best,1}(:,2), 's' )
    % Plot Best Cluster Center
    plot(C(Best,1), C(Best,2),  'kp','LineWidth',4); % Center of the best cluster
    
 
    %[C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS     
    Best_Points = DATA{Best,:};
  
    rtr_idx = ismember(X, Best_Points(:,:), 'rows');
    rtr_idx2 = find(rtr_idx); 
    
    Retrieved_Items =rtr_idx2;




%{
% Retrieval Indexes : DATA{C_index(1)}(:,:)  in the Pareto Points
rtr_idx = ismember(X, Best_Points(:,:), 'rows');
rtr_idx2 = find(rtr_idx); 

YY = X(rtr_idx2,:);

DD = pdist2(YY, Cmn, 'euclid');
DD_new = [rtr_idx2,DD];
DD_new_sort = sortrows(DD_new,2);

%}

predicted_labels = targets(Retrieved_Items,:);
union_of_query_labels = or(targets(queryIndex1, :), targets(queryIndex2, : ));
    
diff = ismember( predicted_labels, union_of_query_labels  , 'rows'); 
num_nz = nnz( diff(:,1) );
s = size(diff(:,1), 1);
    
for j=1:s;
        
    % Cummulative sum of the true-positive elements
    CUMM = cumsum(diff);          
    Precision_AT_K(j,1) = ( CUMM(j,1)  ) ./ j;              
    Recall_AT_K(j,1) = ( CUMM(j,1)  ) ./ (num_nz); % ???????????                    
       
end

acc = num_nz / s;   % accuracy of the best cluster 
%avg_Precision = sum(Precision_AT_K(:,1) ) / s;
avg_Precision = sum(Precision_AT_K(:,1) .* diff(:,1) ) / num_nz;
avg_Precision(isnan(avg_Precision))=0;


handles.Retrieved_Items = Retrieved_Items;
handles.X = X;
handles. Precision_AT_K =  Precision_AT_K;
handles.Recall_AT_K = Recall_AT_K;
handles.acc = acc;
handles.s = s;
handles.avg_Precision = avg_Precision;
    
guidata(hObject, handles);
    



% --- Executes during object creation, after setting all properties.
function FrontSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrontSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ImageSelector_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider






% --- Executes during object creation, after setting all properties.
function ImageSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FrontNum_Callback(hObject, eventdata, handles)
% hObject    handle to FrontNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrontNum as text
%        str2double(get(hObject,'String')) returns contents of FrontNum as a double


% --- Executes during object creation, after setting all properties.
function FrontNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrontNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImageNum_Callback(hObject, eventdata, handles)
% hObject    handle to ImageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageNum as text
%        str2double(get(hObject,'String')) returns contents of ImageNum as a double


% --- Executes during object creation, after setting all properties.
function ImageNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tictoc_Callback(hObject, eventdata, handles)
% hObject    handle to tictoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tictoc as text
%        str2double(get(hObject,'String')) returns contents of tictoc as a double


% --- Executes during object creation, after setting all properties.
function tictoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tictoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function main_Callback(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of main as text
%        str2double(get(hObject,'String')) returns contents of main as a double


% --- Executes during object creation, after setting all properties.
function main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR1_Callback(hObject, eventdata, handles)
% hObject    handle to mainR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR1 as text
%        str2double(get(hObject,'String')) returns contents of mainR1 as a double


% --- Executes during object creation, after setting all properties.
function mainR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR2_Callback(hObject, eventdata, handles)
% hObject    handle to mainR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR2 as text
%        str2double(get(hObject,'String')) returns contents of mainR2 as a double


% --- Executes during object creation, after setting all properties.
function mainR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR3_Callback(hObject, eventdata, handles)
% hObject    handle to mainR3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR3 as text
%        str2double(get(hObject,'String')) returns contents of mainR3 as a double


% --- Executes during object creation, after setting all properties.
function mainR3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL1_Callback(hObject, eventdata, handles)
% hObject    handle to mainL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL1 as text
%        str2double(get(hObject,'String')) returns contents of mainL1 as a double


% --- Executes during object creation, after setting all properties.
function mainL1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL2_Callback(hObject, eventdata, handles)
% hObject    handle to mainL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL2 as text
%        str2double(get(hObject,'String')) returns contents of mainL2 as a double


% --- Executes during object creation, after setting all properties.
function mainL2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL3_Callback(hObject, eventdata, handles)
% hObject    handle to mainL3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL3 as text
%        str2double(get(hObject,'String')) returns contents of mainL3 as a double


% --- Executes during object creation, after setting all properties.
function mainL3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nDCG_Value_Callback(hObject, eventdata, handles)
% hObject    handle to nDCG_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nDCG_Value as text
%        str2double(get(hObject,'String')) returns contents of nDCG_Value as a double


% --- Executes during object creation, after setting all properties.
function nDCG_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nDCG_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_48']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_48;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_64']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_64;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_96']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_96;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_256']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_256;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_Caffe_ssdh48_v3']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_Caffe_ssdh48_v3;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    maxFront = handles.maxFront;
    
    queryIndex1 = handles.q1Idx;
    queryIndex2 = handles.q2Idx;
    data = handles.data;
    q1 = data(queryIndex1,:);
    q2 = data(queryIndex2,:);
    
    
    %     Ranking with EMR
    N = length(handles.filenames);
tic    
    [H A landmarks Z] = EMRcomputeModel(handles.data);
    y1 = zeros(N,1);
    y1(queryIndex1) = 1;
    y2 = zeros(N,1);
    y2(queryIndex2) = 1;
    
    simEMR1 = EMRscore(H ,A, y1);
    simEMR2 = EMRscore(H ,A, y2);
    dist1 = 1-simEMR1;
    dist2 = 1-simEMR2;
       
t = toc;
set(handles.tictoc2,'String',num2str(t))

    X = zeros(2,N);
    X(1,:) = dist1;
    X(2,:) = dist2;
    
    X = (X)';
    
    [K,L] = size(unique(X,'rows'));  %% Number of unique pareto points 
    set(handles.num_pp_v2,'String',num2str(K))
    
    axes(handles.axes3);
    hold off; plot(X(:,1),X(:,2),'.');
    hold on; 
    
     
    
   
    [pf_idx] = pareto_fronts(X, maxFront);
    for k=1:maxFront
        plot(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , 'y-');
    end
    xlabel('c1');
    ylabel('c2'); 

   
    handles.pf_idx = pf_idx;
    handles.X = X;
       
    
    guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function tictoc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tictoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




 

function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double


% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double


% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit40_Callback(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40 as text
%        str2double(get(hObject,'String')) returns contents of edit40 as a double


% --- Executes during object creation, after setting all properties.
function edit40_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit41_Callback(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit41 as text
%        str2double(get(hObject,'String')) returns contents of edit41 as a double


% --- Executes during object creation, after setting all properties.
function edit41_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double


% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





    


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


filenames = handles.filenames;
features = handles.features;
targets = handles.targets;
X = handles.X;
Retrieved_Items = handles.Retrieved_Items;
acc = handles.acc;
s = handles.s; % Number of retrieved items
Precision_AT_K = handles.Precision_AT_K;
Recall_AT_K = handles.Recall_AT_K;
avg_Precision = handles.avg_Precision;

set(handles.acc_box,'String',num2str(acc));
set(handles.avg_prec,'String',num2str(avg_Precision));
set(handles.num_cluster_items,'String',num2str(s))

         
        
        
        
      
      
      axes(handles.axes11);
      hold off;
      %plot(Recall_AT_K, Precision_AT_K);
      x = linspace(0,s);
      plot(Recall_AT_K,Precision_AT_K );
      ylabel('Precision@k' ,'FontSize', 12)
      xlabel('Recall@k' ,'FontSize', 12) 
      hold on;
         
        
         cla(handles.axes13,'reset');
         cla(handles.axes14,'reset');
         cla(handles.axes15,'reset');
         cla(handles.axes16,'reset');
         cla(handles.axes17,'reset');
         cla(handles.axes18,'reset');         
         cla(handles.axes19,'reset');
         cla(handles.axes20,'reset');
         cla(handles.axes21,'reset');
         cla(handles.axes22,'reset');
         cla(handles.axes23,'reset');
         cla(handles.axes24,'reset');
         cla(handles.axes25,'reset');
         cla(handles.axes26,'reset');
         cla(handles.axes27,'reset');
         cla(handles.axes28,'reset');
         cla(handles.axes29,'reset');
         cla(handles.axes30,'reset');
         cla(handles.axes31,'reset');
         cla(handles.axes32,'reset');
         cla(handles.axes33,'reset');
         cla(handles.axes34,'reset');
         cla(handles.axes35,'reset');
         cla(handles.axes36,'reset');     
         
        axes(handles.axes13);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(1,1)}];
        imshow(imread(fname)); 
        set(handles.edit19,'string',num2str( handles.filenames{Retrieved_Items(1,1)}));
        axis image
      
        axes(handles.axes14);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(2,1)}];
        imshow(imread(fname)); 
        set(handles.edit20,'string',num2str( handles.filenames{Retrieved_Items(2,1)}));
        axis image
        
        axes(handles.axes15);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(3,1)}];
        imshow(imread(fname));
        set(handles.edit21,'string',num2str( handles.filenames{Retrieved_Items(3,1)}));
        axis image
        
        axes(handles.axes16);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(4,1)}];
        imshow(imread(fname)); 
        set(handles.edit22,'string',num2str( handles.filenames{Retrieved_Items(4,1)}));
        axis image
        
      
        axes(handles.axes17);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(5,1)}];
        imshow(imread(fname)); 
        set(handles.edit23,'string',num2str( handles.filenames{Retrieved_Items(5,1)}));
        axis image
        
       
        axes(handles.axes18);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(6,1)}];
        imshow(imread(fname)); 
        set(handles.edit24,'string',num2str( handles.filenames{Retrieved_Items(6,1)}));
        axis image
        
       
        axes(handles.axes19);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(7,1)}];
        imshow(imread(fname)); 
        set(handles.edit25,'string',num2str( handles.filenames{Retrieved_Items(7,1)}));
        axis image
        
      
        axes(handles.axes20);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(8,1)}];
        imshow(imread(fname)); 
        set(handles.edit26,'string',num2str( handles.filenames{Retrieved_Items(8,1)}));
        axis image
        
     
        axes(handles.axes21);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(9,1)}];
        imshow(imread(fname)); 
        set(handles.edit27,'string',num2str( handles.filenames{Retrieved_Items(9,1)}));
        axis image
        
     
        axes(handles.axes22);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(10,1)}];
        imshow(imread(fname)); 
        set(handles.edit28,'string',num2str( handles.filenames{Retrieved_Items(10,1)}));
        axis image
        
       
        axes(handles.axes23);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(11,1)}];
        imshow(imread(fname)); 
        set(handles.edit29,'string',num2str( handles.filenames{Retrieved_Items(11,1)}));
        axis image
        
    
        axes(handles.axes24);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(12,1)}];
        imshow(imread(fname)); 
        set(handles.edit30,'string',num2str( handles.filenames{Retrieved_Items(12,1)}));
        axis image
        
     
        axes(handles.axes25);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(13,1)}];
        imshow(imread(fname)); 
        set(handles.edit31,'string',num2str( handles.filenames{Retrieved_Items(13,1)}));
        axis image

        axes(handles.axes26);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(14,1)}];
        imshow(imread(fname)); 
        set(handles.edit32,'string',num2str( handles.filenames{Retrieved_Items(14,1)}));
        axis image
        
       
        axes(handles.axes27);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(15,1)}];
        imshow(imread(fname)); 
        set(handles.edit33,'string',num2str( handles.filenames{Retrieved_Items(15,1)}));
        axis image
       
        axes(handles.axes28);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(16,1)}];
        imshow(imread(fname)); 
        set(handles.edit34,'string',num2str( handles.filenames{Retrieved_Items(16,1)}));
        axis image
        
      
        axes(handles.axes29);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(17,1)}];
        imshow(imread(fname)); 
        set(handles.edit35,'string',num2str( handles.filenames{Retrieved_Items(17,1)}));
        axis image
      
        axes(handles.axes30);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(18,1)}];
        imshow(imread(fname)); 
        set(handles.edit36,'string',num2str( handles.filenames{Retrieved_Items(18,1)}));
        axis image
        
      
        axes(handles.axes31);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(19,1)}];
        imshow(imread(fname)); 
        set(handles.edit37,'string',num2str( handles.filenames{Retrieved_Items(19,1)}));
        axis image
        
        
        axes(handles.axes32);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(20,1)}];
        imshow(imread(fname)); 
        set(handles.edit38,'string',num2str( handles.filenames{Retrieved_Items(20,1)}));
        axis image
      
        axes(handles.axes33);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(21,1)}];
        imshow(imread(fname)); 
        set(handles.edit39,'string',num2str( handles.filenames{Retrieved_Items(21,1)}));
        axis image
        
     
        axes(handles.axes34);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(22,1)}];
        imshow(imread(fname)); 
        set(handles.edit40,'string',num2str( handles.filenames{Retrieved_Items(22,1)}));
        axis image
        
       
        axes(handles.axes35);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(23,1)}];
        imshow(imread(fname)); 
        set(handles.edit41,'string',num2str( handles.filenames{Retrieved_Items(23,1)}));
        axis image
        
        axes(handles.axes36);
        fname = [handles.image_dir handles.filenames{Retrieved_Items(24,1)}];
        imshow(imread(fname)); 
        set(handles.edit42,'string',num2str( handles.filenames{Retrieved_Items(24,1)}));
        axis image
            
      
 guidata(hObject, handles);
       


% --- Executes on selection change in hashCodeSelection_f.
function hashCodeSelection_f_Callback(hObject, eventdata, handles)
% hObject    handle to hashCodeSelection_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns hashCodeSelection_f contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hashCodeSelection_f


feature_dir = ['/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/features/'];
image_dir =['/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/imageFolder/']; 
data_dir = ['/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes'];
load([data_dir '/filenames']); % File names
load([data_dir '/targets']);   % Labels

hashCode_index = get(handles.hashCodeSelection_f, 'Value');

switch hashCode_index
           
    case 1
        load([data_dir '/hashCodes_64']); 
        data = hashCodes_64;
        load([feature_dir '/features_64']); 
        features = features_64;
        %data = features_64 > 0.5;
    case 2
       load([data_dir '/hashCodes_128']); 
       data = hashCodes_128;
       load([feature_dir '/features_128']); 
       features = features_128;
       %data = features_128 > 0.5;
    case 3
        load([data_dir '/hashCodes_256']); 
        data = hashCodes_256;
        load([feature_dir '/features_256']); 
        features = features_256;
        %data = features_256 > 0.5;
    case 4
        load([data_dir '/hashCodes_512']); 
        data = hashCodes_512;
        load([feature_dir '/features_512']); 
        features = features_512;
        %data = features_512 > 0.5;
    
        
       
end




set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
%set(handles.QueryName3,'String', filenames);


handles.image_dir = image_dir;
handles.data_dir = data_dir;
handles.filenames = filenames;
handles.targets = targets;
handles.data = data;
handles.features = features;



guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function hashCodeSelection_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hashCodeSelection_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_pp_v2_Callback(hObject, eventdata, handles)
% hObject    handle to num_pp_v2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_pp_v2 as text
%        str2double(get(hObject,'String')) returns contents of num_pp_v2 as a double


% --- Executes during object creation, after setting all properties.
function num_pp_v2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_pp_v2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function pushbutton11_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function number_of_Clusters_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_Clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_Clusters as text
%        str2double(get(hObject,'String')) returns contents of number_of_Clusters as a double

num_Clusters = str2double(get(handles.number_of_Clusters,'string'));
handles.num_Clusters = num_Clusters;

% --- Executes during object creation, after setting all properties.
function number_of_Clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_Clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_pr_po_Callback(hObject, eventdata, handles)
% hObject    handle to num_pr_po (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_pr_po as text
%        str2double(get(hObject,'String')) returns contents of num_pr_po as a double


% --- Executes during object creation, after setting all properties.
function num_pr_po_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_pr_po (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acc_box_Callback(hObject, eventdata, handles)
% hObject    handle to acc_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of acc_box as text
%        str2double(get(hObject,'String')) returns contents of acc_box as a double


% --- Executes during object creation, after setting all properties.
function acc_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acc_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avg_prec_Callback(hObject, eventdata, handles)
% hObject    handle to avg_prec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avg_prec as text
%        str2double(get(hObject,'String')) returns contents of avg_prec as a double


% --- Executes during object creation, after setting all properties.
function avg_prec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avg_prec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_cluster_items_Callback(hObject, eventdata, handles)
% hObject    handle to num_cluster_items (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_cluster_items as text
%        str2double(get(hObject,'String')) returns contents of num_cluster_items as a double


% --- Executes during object creation, after setting all properties.
function num_cluster_items_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_cluster_items (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
