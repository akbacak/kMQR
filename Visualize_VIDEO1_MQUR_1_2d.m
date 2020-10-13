
close all;
clear all;
clc;

load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/hashCodes_128.mat');
data = hashCodes_128;

load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/targets.mat');
targets = targets;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/filenames.mat');
filenames = filenames;
N = length(filenames);

queryIndex = xlsread('/home/ubuntu/keras/enver/dmlvh2/DMQVR/qLabels_1.xls'); 

queryIndex1 = 1;
queryIndex2 = 10;

q_1 = data(queryIndex1,:); 
q_2 = data(queryIndex2,:); 
%q3 = data(queryIndex3,:); 

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
input = unique(X, 'rows');
hold off; 
%plot(X(:,1),X(:,2),'o','MarkerFaceColor',[0 1 0]); 
plot(X(:,1),X(:,2),'k.');
ylabel('d_2 ', 'FontSize', 50);
xlabel('d_1 ', 'FontSize', 50);
hold on;
set(gca,'FontSize',35); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
union_of_query_labels = or(targets(queryIndex1, :), targets(queryIndex2, : ));
absolute_union_of_query_labels = nnz(union_of_query_labels);   
    
 for e = 1:N        
            MQUR_ALL(e,:) =  nnz( and(targets(e,:) , union_of_query_labels ) ) / absolute_union_of_query_labels ;
            
 end

MQUR_ONE = find(MQUR_ALL == 1);
%plot(X(MQUR_ONE,1),X(MQUR_ONE,2),'gs');
plot(X(MQUR_ONE,1),X(MQUR_ONE,2),'gd', 'LineWidth', 2);



%{
K= 20;% Number Of Clusters

[idx,C] = kmeans(X,K) ;  % Applay kmeans, C:Coordinates of the Clusters, D:Whinhin distance in a cluster 

% Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
DATA= cell(K,1) ;

for i = 1:K
    DATA{i} = X(idx==i, :) ;
    %VAR(i,:) = var(DATA{i,:},0 ,1);   %VARIANCES ALONG COLUMNS 
    %scatter3(X(idx==i, 1), X(idx==i, 2),  X(idx==i, 3)  ,'k.')    
end

    
    
    %{
    VAR = var(C,0,2);
    [VAR_sort , VAR_index] = sort(VAR);
    VAR_index = VAR_index(1:3, : );
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = [VAR_index, SUM ];
    LAST_sort = sortrows(LAST, 2);
    Best = LAST_sort(1,1);
    %}

    
    %{
    VAR = var(C,0,2);
   
    [VAR_sort , VAR_index] = sort(VAR);
    VAR = VAR_sort(1:3);
    VAR_index = VAR_index(1:3, : );
    %SUM = sum(sqrt(C.^2) ,2);  
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = VAR + SUM;
    [C_min, Best] = min(LAST);    
    %}

%{
[C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS 
C_var  = var(C,0,2); % Variances  along each column
C_mean = mean(C, 2); % Mean of columns along each row
C_new = [C, C_mean, C_var];
C_new_sorted = sortrows(C_new, [4 5 ]);
[C_min, Best] = min(var(C_new_sorted(:,[1 2 3]),0,2));
%}


[C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS 
Best_Points = DATA{Best}(:,:);

% Plot the best one 
plot(DATA{Best,1}(:,1) , DATA{Best,1}(:,2),    'o', 'MarkerSize',10 , 'MarkerFaceColor',[1 0 1]);

% Plot Best Cluster Center
plot(C(Best,1), C(Best,2), 'd', 'MarkerSize',20 , 'MarkerFaceColor',[1 0 0]);


% Retrieval Indexes : DATA{C_index(1)}(:,:)  in the Pareto Points

rtr_idx = ismember(X, Best_Points(:,:), 'rows');
Retrieved_Items = find(rtr_idx); 

plot(X(Retrieved_Items,1),X(Retrieved_Items,2),'d', 'MarkerFaceColor',[1 0 0]);


%{
%%%%%%%%%%%%%%%%%%%%  ReRanking %%%%%%%%%%% 

Retrieved_Items_features = features(Retrieved_Items(:,1), :);

f1 = features(queryIndex1,:); 
f2 = features(queryIndex2,:);
f3 = features(queryIndex3,:);


dist_f1 = pdist2(f1 , Retrieved_Items_features , 'euclid' );
dist_f2 = pdist2(f2 , Retrieved_Items_features , 'euclid' );
dist_f3 = pdist2(f3 , Retrieved_Items_features , 'euclid' );

[M2,B] = size(Retrieved_Items_features(:,1));
YY = zeros(3,M2);
YY(1,:) = dist_f1;
YY(2,:) = dist_f2;
YY(3,:) = dist_f3;
YY = (YY)';

for i = 1:M2
    
    plot3(YY(:, 1), YY(:,2), YY(:, 3)  ,'.')    
end

qf1 = YY(end-2,:); 
qf2 = YY(end-1,:); 
qf3 = YY(end,:);

Cmn_f = (qf1(:) + qf2(:) + qf3(:)).'/2;
plot3(Cmn_f(1), Cmn_f(2), Cmn_f(3), 's' , 'LineWidth', 2);

Dissim_f                = EuDist2(YY, Cmn_f, 'euclid'); % Rank by EMR then Rerank by Euclidian ! ?
Dissim_ff               = [Retrieved_Items_features, Dissim_f];
DissimSorted_f          = sortrows(Dissim_ff, 2);

Retrieved_Items_Ranked  = DissimSorted_f(:,1);


%}
%%%%%%%%%%%%%%%%%%%%%%% Metrics & MQUR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Retrieved_Items_Ranked = Retrieved_Items;
predicted_labels  = targets(Retrieved_Items,:);

[num_R, ~]  = size(Retrieved_Items);     
 for e = 1:num_R         
     MQUR_Ranked(e,:) =  nnz( and(predicted_labels(e,:) , union_of_query_labels ) ) / absolute_union_of_query_labels;       
 end
      
for ff = 1:num_R
    if  MQUR_Ranked(ff,:) ~= 1   %%%%%% MQUR !DEN FARKLI OLANLARI SEÇ 
        Retrieved_Items_Ranked(ff,:) = 0;   
        MQUR_Ranked(ff,:) = 0;
    end
end

Retrieved_Items_Ranked( all(~Retrieved_Items_Ranked,2), : ) = [];
MQUR_Ranked( all(~MQUR_Ranked,2), : ) = [];
plot(X(Retrieved_Items_Ranked,1),X(Retrieved_Items_Ranked,2),'bs','MarkerSize',10 , 'MarkerFaceColor',[1 1 0]);
 
 
predicted_labels_ranked  = targets(Retrieved_Items_Ranked,:);
%diff = ismember( predicted_labels_ranked, union_of_query_labels , 'rows' );
if isempty(MQUR_Ranked)
   MQUR_Ranked = 0;
end
  
num_nz = nnz( MQUR_Ranked(:,1) );
s = size(MQUR_Ranked(:,1), 1);
    
for j=1:s;        
    %Cummulative sum of the true-positive elements
    CUMM = cumsum(MQUR_Ranked);          
    Precision_AT_K(j,1) = ( CUMM(j,1)  ) ./ j;              
    Recall_AT_K(j,1) = ( CUMM(j,1)  ) ./ (num_nz); %              
end

avg_Precision = sum(Precision_AT_K(:,1) .* MQUR_Ranked(:,1) ) / num_nz;
avg_Precision(isnan(avg_Precision))=0;
% avg_Precision_OLD = sum(Precision_AT_K(:,1) ) / s;
 acc = num_nz / s;   % accuracy of the best cluster 


%{
plot(Recall_AT_K, Precision_AT_K);
hold off;
x = linspace(0,s);
plot( Precision_AT_K )
ylabel('Precision@k' ,'FontSize', 12)
xlabel('Number of Rterieved Items' ,'FontSize', 12) 
hold on;
%}


%}