
% https://www.mathworks.com/matlabcentral/answers/361878-getting-the-data-points-of-each-cluster-in-kmeans
close all;
clear all;
clc;

load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/dmlh2/hashCodes_64.mat');
data = hashCodes_64;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/features/features_64.mat');
features = features_64;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/targets.mat');
targets = targets;
load('lamdaDataset/hashCodes/filenames.mat');
filenames = filenames;


queryIndex1 = 562;
queryIndex2 = 1761;

q_1 = data(queryIndex1,:); 
q_2 = data(queryIndex2,:); 


N = length(filenames); 
q1new = repmat(q_1,N,1);
q2new = repmat(q_2,N,1);

dist_1 = xor(data, q1new);
dist_2 = xor(data, q2new);

hamming_dist1 = sum(dist_1,2);
hamming_dist2 = sum(dist_2,2);

n_hamming_dist1 = mat2gray(hamming_dist1);
n_hamming_dist2 = mat2gray(hamming_dist2);
 
 
X = zeros(2,N);
X(1,:) = hamming_dist1;
X(2,:) = hamming_dist2;

X = (X)';

K= 20; % Number Of Clusters

[idx,C] = kmeans(X,K) ;  % Applay kmeans, C:Coordinates of the Clusters, D:Whinhin distance in a cluster 

% Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
DATA= cell(K,1) ;
figure;
set(gca,'FontSize',30); 
hold on
for i = 1:K
    DATA{i} = X(idx==i, :) ;
    plot(X(idx==i, 1), X(idx==i, 2), 'o','MarkerFaceColor',[0 1 0],'MarkerSize',15)  
   
end
 ylabel('d_1 ', 'FontSize', 50)
 xlabel('d_2 ', 'FontSize', 50)
           

    VAR = var(C,0,2);
    [VAR_sort , VAR_index] = sort(VAR);
    VAR_index = VAR_index(1:10, : );
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = [VAR_index, SUM ];
    LAST_sort = sortrows(LAST, 2);
    Best = LAST_sort(1,1);
                                        
 
    %[C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS     
    Best_Points = DATA{Best,:};
  
    rtr_idx = ismember(X, Best_Points(:,:), 'rows');
    rtr_idx2 = find(rtr_idx); 
    
    Retrieved_Items =rtr_idx2;
  



%{ 
%%%%%%%%%%%%%%%%%%%%%%%%% Reranking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add queries to Feature Pareto space
rtr_idx2(end+1,:) = queryIndex1;
rtr_idx2(end+1,:) = queryIndex2;

rtr_idx2_features = features(rtr_idx2(:,1), :);

f1 = features(queryIndex1,:); 
f2 = features(queryIndex2,:);

dist_f1 = pdist2(f1 , rtr_idx2_features , 'euclid' );
dist_f2 = pdist2(f2 , rtr_idx2_features , 'euclid' );

[M2,B] = size(rtr_idx2_features(:,1));
YY = zeros(2,M2);
YY(1,:) = dist_f1;
YY(2,:) = dist_f2;
YY = (YY)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



avg_Precision = sum(Precision_AT_K(:,1) .* diff(:,1) ) / num_nz;
avg_Precision(isnan(avg_Precision))=0;
% avg_Precision_OLD = sum(Precision_AT_K(:,1) ) / s;
acc = num_nz / s;   % accuracy of the best cluster 


