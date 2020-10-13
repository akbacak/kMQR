
close all;
clear all;
clc;

load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/streetsDataset/hashCodes/hashCodes_16.mat');
data = hashCodes_16;
%load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/streetsDataset/features/features_64.mat');
%features = features_64;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/streetsDataset/hashCodes/targets.mat');
targets = targets;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/streetsDataset/hashCodes/filenames.mat');
filenames = filenames;
N = length(filenames);

queryIndex = xlsread('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/streetsDataset/streets_3d.xls');  % Reads randomly choosen query pairs from excell file
queryIndex = transpose( queryIndex ); 
queryIndex1 = queryIndex(1,:);        % First element of Query Triplets
queryIndex2 = queryIndex(2,:);        % Second element of Query Triplets
queryIndex3 = queryIndex(3,:); 
      
  
 K= 20; % Number Of Clusters

 for l = 1:240 % Number of Query Pairs , CAN TRY FOR DIFFERENT HAMMING RADIUS ALSO ?????? how?
              
       q_1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
       q_2 = data(queryIndex2,:);
       q_3 = data(queryIndex3,:);
       q1_rep{l,:} = repmat(q_1(l,:),N,1); % Make query matrix size to the same as data matrix size
       q2_rep{l,:} = repmat(q_2(l,:),N,1);  
       q3_rep{l,:} = repmat(q_3(l,:),N,1); 
       xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
       xor_data_q2new{l,:} = xor(data, q2_rep{l,:}); 
       xor_data_q3new{l,:} = xor(data, q3_rep{l,:});
       hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
       hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
       hamming_dist3{l,:} = sum(xor_data_q3new{l,:},2);
       %norm_hamming_dist1{l,:} =  hamming_dist1{l,:} / max( hamming_dist1{l,:}(:) ); % Normalize hamming  distances between 0&1
       %norm_hamming_dist2{l,:} =  hamming_dist2{l,:} / max( hamming_dist2{l,:}(:) );
       %dist1{l,:} = mat2gray(dist1{l,:}); % Normalize hamming  distances between 0&1
       %dist2{l,:} = mat2gray(dist2{l,:});     
        
       X = zeros(3,N);
       X(1,:) = hamming_dist1{l,:};
       X(2,:) = hamming_dist2{l,:};
       X(3,:) = hamming_dist3{l,:};
 
       X = (X)';
       
       
     %%%%%%%%%%%%%%%%%%%% Pareto Fronts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     maxFront =   3; 
     pf_idx(l,:) = pareto_fronts(X, maxFront);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %for e = 1:N        
        %           MQUR_ALL(e,:) =  nnz( and(targets(e,:) , union_of_query_labels ) ) / absolute_union_of_query_labels ;
            
        %end

        %MQUR_ONE = find(MQUR_ALL == 1);
        %scatter3(X(MQUR_ONE,1),X(MQUR_ONE,2),X(MQUR_ONE,3),'gs');

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [idx,C] = kmeans(X,K);  % Applay kmeans, C:Coordinates of the Clusters, 
                                        % D:Within distance in a cluster
                                        % idx, her bir nokatnın ait olduğu
                                        % cluster no'yu verir.                                        
 

    % Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
    DATA{l,:}= cell(K,1) ;
 
  
    for i = 1:K
        DATA{l,:}{i,:} = X(idx==i, :) ; % DİKKAT, BURDA 2 CELL ARRAY VAR!
       
    end
    
    %{
    VAR = var(C,0,2);
    [VAR_sort , VAR_index] = sort(VAR);
    VAR_index = VAR_index(1:1, : );
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = [VAR_index, SUM ];
    LAST_sort = sortrows(LAST, 2);
    Best = LAST_sort(1,1);
    %}
    
    [C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS   
    Best_Points(l,:) = DATA{l,:}(Best,:);
    
  
    rtr_idx{l,:} = ismember(X, Best_Points{l,:}(:,:),'rows'); 
    % Find indexes of 1 values, means we find Best_Points indexes in the X
    rtr_idx2{l,:} = find(rtr_idx{l,:}); 
    
    % Find indexes of 1 values, means we find Best_Points indexes in the X
    Retrieved_Items{l,:} = find(rtr_idx{l,:}); 
       
    
    %{
    % Add queries to Feature Pareto space
    rtr_idx2{l,:}(end+1,:) = queryIndex1(1,l);
    rtr_idx2{l,:}(end+1,:) = queryIndex2(1,l);
    rtr_idx2{l,:}(end+1,:) = queryIndex3(1,l);
 
    
    % RERANKING !!!!!!!!!!!!
    %  Use cnn features of retrieved items by hash codes
    rtr_idx2_features{l,:} = features(rtr_idx2{l,:}(:,1), :);
 
    % features of query pairs
    f1 = features(queryIndex1,:); 
    f2 = features(queryIndex2,:);
    f3 = features(queryIndex3,:);
     
    % Distance from each query pair features to retrieved items 
    dist_f1{l,:} = pdist2(f1(l,:) , rtr_idx2_features{l,:} , 'euclid' );
    dist_f2{l,:} = pdist2(f2(l,:) , rtr_idx2_features{l,:} , 'euclid' );
    dist_f3{l,:} = pdist2(f3(l,:) , rtr_idx2_features{l,:} , 'euclid' );
    
    % How many rows of trr_idx2
    [ M(l,:), ~] = size(rtr_idx2{l,:}(:,1));

    % Create  2xM zero vector for assigne each distance (dis_f1 and f2) to them
    % YY is Pareto space formed by cnn features of retrived itmes,
    % which are retrieved by has codes
    YY{l,:} = zeros(3,M(l,:)); 
    YY{l,:}(1,:) = dist_f1{l,:};
    YY{l,:}(2,:) = dist_f2{l,:};
    YY{l,:}(3,:) = dist_f3{l,:};
    YY{l,:} = (YY{l,:})';
    
    qf1(l,:) = YY{l,:}(end-2,:); 
    qf2(l,:) = YY{l,:}(end-1,:); 
    qf3(l,:) = YY{l,:}(end,:);
 
    Cmn_f(l,:) = (qf1(l,:) + qf2(l,:) + qf3(l,:)).'/2;
    %[Cmn_f] = center_3d(qf1(l,:), qf2(l,:), qf3(l,:), 'circumcenter');
    %Cmn_f = Cmn_f'; 
 
    % DD is the distance between each retrieved items to optimum point
    DD{l,:} = pdist2(YY{l,:}, Cmn_f, 'euclid');
 
    DD_new{l,:} = [rtr_idx2{l,:},DD{l,:}];
    DD_new_sort{l,:} = sortrows(DD_new{l,:},2);
    
    Retrieved_Items{l,:} = DD_new_sort{l,:}(:,1);
    % Now remove last two rows (qf1&qf2) from Retrieved items
    Retrieved_Items{l,:}(end) = [];
    Retrieved_Items{l,:}(end) = [];
    Retrieved_Items{l,:}(end) = [];

    %%%%%%%%%%%%% Choose First P Rterieved Items %%%%%%%%%%%%%%%%%%%%%%%%
    %P = 5;
    %Retrieved_Items{l,:} = Retrieved_Items{l,:}(1:P, :);
    %predicted_labels{l,:} = targets(Retrieved_Items(l, 1:P),:);

    %}
   
    Retrieved_Items_Ranked = Retrieved_Items;    
    predicted_labels{l,:} = targets(Retrieved_Items_Ranked{l,:} , :);
   
    
    
    
    predicted_labels_MQR{l,:} = targets(Retrieved_Items{l,:},:);
    union_of_query_labels_MQR{l,:} = or(targets(queryIndex(1,l), :), targets(queryIndex(2,l), : ));
    union_of_query_labels_MQR{l,:} = or(union_of_query_labels_MQR{l,:},  targets(queryIndex(3,l), : ));
    absolute_union_of_query_labels_MQR{l,:} = nnz(union_of_query_labels_MQR{l,:} );     
    
    diff_MQR{l,:} = ismember( predicted_labels_MQR{l,:}, union_of_query_labels_MQR{l,:}  , 'rows'); 
    if isempty(diff_MQR{l,:})
            diff_MQR{l,:} = 0;
    end
    num_nz_MQR(l,:) = nnz( diff_MQR{l,:}(:,1) );
    s_MQR{l,:} = size(diff_MQR{l,:}(:,1), 1);
    
    for j=1:s_MQR{l,:};
        
        % Cummulative sum of the true-positive elements
        CUMM_MQR{l,:} = cumsum(diff_MQR{l,:});          
        Precision_AT_K_MQR{l,:}(j,1) = ( CUMM_MQR{l,:}(j,1)  ) / j;              
        Recall_AT_K_MQR{l,:}{j,1} = ( CUMM_MQR{l,:}(j,1)  ) / (num_nz_MQR(l,:)); %?????????????                    
       
    end  
    
    acc_MQR(l,:) = num_nz_MQR(l,:) / s_MQR{l,:};   % accuracy of the best cluster 
    %avg_Precision_MQR(l,:) = sum(Precision_AT_K_MQR{l,:}(:,1) ) / s_MQR{l,:};
    avg_Precision_MQR(l,:) = sum(Precision_AT_K_MQR{l,:}(:,1) .* diff_MQR{l,:}) / num_nz_MQR(l,:);
    avg_Precision_MQR(isnan(avg_Precision_MQR))=0;

 
    %%%%%%%%%%%%%%%%%%%%% nDCG Scores of MQR %%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R_MQR(l,:)] = size(Best_Points{l,:}(:,1), 1);     
    for e = 1:  R_MQR(l,:)         
            MQUR_ALL_MQR{l,:}(e,:) =  nnz( and(predicted_labels_MQR{l,:}(e,:) , union_of_query_labels_MQR{l,:} ) ) / absolute_union_of_query_labels_MQR{l,:} ;
                                              
    end 
    
     DCG_MQR{l,:}(1,1) = MQUR_ALL_MQR{l,:}(1,1);
     i_DCG_MQR{l,:}(1,1) = 1;
     
    for e = 2:  R_MQR(l,:) 
        DCG_MQR{l,:}(e,:) = MQUR_ALL_MQR{l,:}(1,1) + cumsum(MQUR_ALL_MQR{l,:}(e,:) ./log2(e)) ;
        i_DCG_MQR{l,:}(e,:) = 1 + cumsum(1 /log2(e)) ;                                              
    end 
    
    nDCG_MQR{l,:} =   DCG_MQR{l,:} ./   i_DCG_MQR{l,:};  
    TnDCG_MQR{l,:} = transpose(nDCG_MQR{l,:});
    size_of_nDCG_MQR(l,:) = size( nDCG_MQR{l,:}, 1);
 
 
   %%%%%%%%%%%%%%%%%%%% Precision & Recall for DMQIR  %%%%%%%%%%%%%%%%%%%%
    
    % NOW FIND Precision values of first 10 front for each query pairs
    for j= 1:maxFront
        predicted_labels_DMQIR{l,j} =      targets(pf_idx{l,j}(:,4),:);
        union_of_query_labels_DMQIR{l,:} = or(targets(queryIndex(1,l), :), targets(queryIndex(2,l), : ));
        union_of_query_labels_DMQIR{l,:} = or(union_of_query_labels_DMQIR{l,:} , targets(queryIndex(3,l), : ));
        absolute_union_of_query_labels_DMQIR{l,:}(j,:) = nnz(union_of_query_labels_DMQIR{l,:} ); 
        
        diff_DMQIR{l,j} = ismember( predicted_labels_DMQIR{l,j}, union_of_query_labels_DMQIR{l,:}  , 'rows'); 
        if isempty(diff_DMQIR{l,j})
            diff_DMQIR{l,j} = 0;
        end
        num_nz_DMQIR(l,j) = nnz( diff_DMQIR{l,j}(:,1) );
        s_DMQIR{l,j} = size(diff_DMQIR{l,j}(:,1), 1);
        
        for jj=1:s_DMQIR{l,j};
        
        % Cummulative sum of the true-positive elements
        CUMM_DMQIR{l,j} = cumsum(diff_DMQIR{l,j});          
        Precision_AT_K_DMQIR{l,j}(jj,1) = ( CUMM_DMQIR{l,j}(jj,1)  ) / jj;              
        Recall_AT_K_DMQIR{l,j}{jj,1} = ( CUMM_DMQIR{l,j}(jj,1)  ) / (num_nz_DMQIR(l,j)); %?????????????                    
       
        end  
        
        acc_DMQIR(l,j) = num_nz_DMQIR(l,j) / s_DMQIR{l,j};   % accuracy of the best cluster 
        %avg_Precision_DMQIR(l,j) = sum(Precision_AT_K_DMQIR{l,j}(:,1) ) / s_DMQIR{l,j};
        avg_Precision_DMQIR(l,j) = sum(Precision_AT_K_DMQIR{l,j}(:,1) .* diff_DMQIR{l,j} ) / num_nz_DMQIR(l,j);
        avg_Precision_DMQIR(isnan(avg_Precision_DMQIR))=0;

        
        
    end  
    
        mean_avg_Precision_DMQIR(l,:) = sum(avg_Precision_DMQIR(l,:)) / maxFront;
        mean_acc_DMQIR(l,:) = sum(acc_DMQIR(l,:)) / maxFront;
       
 end

mAP_MQR = sum(avg_Precision_MQR(:,1)) / l;
avg_acc_MQR = sum(acc_MQR(:,1)) / l;

mAP_DMQIR = sum(mean_avg_Precision_DMQIR(:,1)) / l;
avg_acc_DMQIR = sum(mean_acc_DMQIR(:,1)) / l;
 


%%%%%%%%%%%%%%%%%%% MEAN AND AREA UNDER nDCG_MEAN_MQR CURVE %%%%%%%%%%%%%%%%%%%
 max_size_of_nDCG_MQR = max(size_of_nDCG_MQR);

for l=1:240 
    nDCG_zp_MQR{l,:}(1,:) = [TnDCG_MQR{l,:}(1,:)  , (zeros(1,  (max_size_of_nDCG_MQR -  size_of_nDCG_MQR(l,:) )  ) )  ];
  
end
nDCG_mean_kMQR = mean(cat(1, nDCG_zp_MQR{:}));

plot(nDCG_mean_kMQR);
     
% Area under nDCG Curve
x = 1:max_size_of_nDCG_MQR;   % Square brackets waste time here only
Area_Under_nDCG_mean_kMQR = trapz(x, nDCG_mean_kMQR); 
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% nDCG scores and Area_Under_n_DCG_mean_DMQIR %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for l = 1:240                 % Number of Query Pairs
              
        q1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
        q2 = data(queryIndex2,:);
        q1_rep{l,:} = repmat(q1(l,:),N,1); % Make query matrix size to the same as data matrix size
        q2_rep{l,:} = repmat(q2(l,:),N,1);      
        xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
        xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
        hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
        hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
        %norm_hamming_dist1{l,:} =  hamming_dist1{l,:} / max( hamming_dist1{l,:}(:) ); % Normalize hamming  distances between 0&1
        %norm_hamming_dist2{l,:} =  hamming_dist2{l,:} / max( hamming_dist2{l,:}(:) );
        %dist1{l,:} = mat2gray(dist1{l,:}); % Normalize hamming  distances between 0&1
        %dist2{l,:} = mat2gray(dist2{l,:});     
        
        X = zeros(2,N);
        X(1,:) = hamming_dist1{l,:};
        X(2,:) = hamming_dist2{l,:};
    
        X = (X)';
        [K,~] = size(unique(X,'rows')); % The number of PPs = K

                
             
        [pf_idx] = pareto_fronts(X, maxFront);   
        
        
        q1_label{l,:} = targets(queryIndex1(:,l), : ); % Label vector of Query 1
        q2_label{l,:} = targets(queryIndex2(:,l), : ); % Label vector of Query 2
        
        b{l,:} = or(q1_label{l,:} , q2_label{l,:});    % beta in the equation 7 
        absolute_b{l,:} = nnz(b{l,:});                 % Number of non-zero elements in the beta, nnz is a Matlab Func.
        
         
        for j  =1:maxFront
            
            Labels{l,j} = targets(pf_idx{j,1}(:,3),:);
                                     
            [R(l,j) , C(l,j)] = size(Labels{l,j});        
            
            
                      
            
            switch (mod(R(l,j) ,2) == 1)
                case 1
                     e_left(l,j) = (round(R(l,j) / 2) - 1) ;  
                     e_rigth(l,j)= (round(R(l,j) / 2)    ) ;
                case 0
                     e_left(l,j)  =  R(l,j) / 2 ;
                     e_rigth(l,j) =  R(l,j) / 2 ;
            end
            
            
           
            % ALL MQUR scores for each query pairs,  each front X Number of retrieved items               
            for e = 1:R(l,j)               
                MQUR_ALL{l,:}(j,e) =  nnz( and(Labels{l,j}(e,:) , b{l,:} ) ) /  absolute_b{l,:};
                                           
            end                   
            %MQUR score of the first retrieved items in each front
            MQUR_1(l,j) = MQUR_ALL{l,:}(j,1);
                        
            
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rtr_idx{j,1} = [];

    
            [v(j) , k(j)] = size(MQUR_ALL{l,1}(j,:)); % k(ll) is the size of column vector of ll th front
    
             for ff = 1:k(j)
                 if  MQUR_ALL{l,1}(j,ff) == 1    
                     
                     rtr_idx{1,1}(end+1,:) = pf_idx{j,1}( ff , 3);
                     %rtr_idx{j,1}(end+1,:) = pf_idx{j,1}( ff , 3);
                     %rtr_idx{ll,:}(end+1,:) = pf_idx{ll,1}( ff , :);
                 end
             end
             
             [kk(j) zz(j)] = size( rtr_idx{j,1}(:) );
             kk = kk';
             count_rtr = sum(kk(:)); 
          
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
           
            
           
            
             
            G{l,j} = [];
            rigth_G{l,j} = [];
            left_G{l,j}  = [];
            i_G{l,j} =[];
            i_rigth_G{l,j} = [];
            i_left_G{l,j}  = [];
            
            rigth_DG{l,j} = [];
            left_DG{l,j}  = [];
            i_DG{l,j} = [];
            i_rigth_DG{l,j} = [];
            i_left_DG{l,j}  = [];
            
            rigth_DCG{l,j} = [];
            left_DCG{l,j}  = [];
            i_rigth_DCG{l,j} = [];
            i_left_DCG{l,j}  = [];
            
            
             for s = 1:R(l,j)     
               G{l,j}(:,s) = MQUR_ALL{l,:}(j,s);
               i_G{l,j}(:,s) = 1;
             end
             
             switch (mod(R(l,j) ,2) == 1)
                case 1                    
                    rigth_G{l,j}   =   G{l,j}(: , (e_rigth(l,j)   )  : R(l,j));
                    i_rigth_G{l,j} = i_G{l,j}(: , (e_rigth(l,j)   )  : R(l,j));
                    
                    left_G{l,j}    =   G{l,j}(: , 1 : e_left(l,j) ); 
                    left_G{l,j}    =   fliplr(left_G{l,j});
                    i_left_G{l,j}   = i_G{l,j}(: , 1 : e_left(l,j) );
             
                    rigth_DG{l,j}(:,1)   =   rigth_G{l,j}(:,1);
                    i_rigth_DG{l,j}(:,1) = i_rigth_G{l,j}(:,1);
                    
             
             
                    %left_DG{l,j}(:,1)   =   left_G{l,j}(:,1);
                    left_DG   =   left_G;
                    i_left_DG{l,j}(:,1) = i_left_G{l,j}(:,1);
             
                    for d = 2 : e_rigth(l,j)
                         rigth_DG{l,j}(:,d)   =   rigth_G{l,j}(:,d) / log2(d); %%%%%%%%%%%%%%%
                         i_rigth_DG{l,j}(:,d) = i_rigth_G{l,j}(:,d) / log2(d); 
                    end                    
                    for  f = 2 : e_left(l,j)   
                         left_DG{l,j}(:,f)   =   left_G{l,j}(:,f) / log2(f);
                         i_left_DG{l,j}(:,f) = i_left_G{l,j}(:,f) / log2(f);
                    end
                    
                    
                 case 0
                    rigth_G{l,j}   =   G{l,j}(: , (e_rigth(l,j) +1)  : R(l,j));
                    i_rigth_G{l,j} = i_G{l,j}(: , (e_rigth(l,j) +1)  : R(l,j));
                    
                    left_G{l,j}    =   G{l,j}(: , 1 : e_left(l,j) ); 
                    left_G{l,j}    =   fliplr(left_G{l,j});
                    i_left_G{l,j}   = i_G{l,j}(: , 1 : e_left(l,j) );
             
                    rigth_DG{l,j}(:,1)   =   rigth_G{l,j}(:,1);
                    i_rigth_DG{l,j}(:,1) = i_rigth_G{l,j}(:,1);
                    %rigth_DG{l,j}(:,e_rigth(l,j)) =  MQUR_ALL{l,:}(j,e_rigth(l,j));
             
             
                    left_DG{l,j}(:,1)   =   left_G{l,j}(:,1);
                    i_left_DG{l,j}(:,1) = i_left_G{l,j}(:,1);
             
             
                    for h = 2 : e_rigth(l,j)
                         rigth_DG{l,j}(:,h)   =   rigth_G{l,j}(:,h) / log2(h);
                         i_rigth_DG{l,j}(:,h) = i_rigth_G{l,j}(:,h) / log2(h); 
                    end
                    for m = 2 : e_left(l,j)
                         left_DG{l,j}(:,m)   =   left_G{l,j}(:,m) / log2(m);
                         i_left_DG{l,j}(:,m) = i_left_G{l,j}(:,m) / log2(m);
                    end
             end
             
                                    
             rigth_DCG{l,j} = cumsum(rigth_DG{l,j}); 
             i_rigth_DCG{l,j} = cumsum( i_rigth_DG{l,j});
             n_rigth_DCG{l,j} = rigth_DCG{l,j} ./ i_rigth_DCG{l,j};
                                  
           
             left_DCG{l,j} = cumsum(left_DG{l,j}); 
             i_left_DCG{l,j} = cumsum( i_left_DG{l,j});
             n_left_DCG{l,j} = left_DCG{l,j} ./ i_left_DCG{l,j};
             
             flip_n_left_DCG{l,j} = fliplr(n_left_DCG{l,j});
             n_DCG{l,j} = horzcat(flip_n_left_DCG{l,j}, n_rigth_DCG{l,j});
                 
            
           
             
        end
 end
    
 for ll=1:l
    for jj=1:j
        [R_1(ll,jj) , C_1(ll,jj)] = size(n_rigth_DCG{ll,jj}); % Find size of each matrix in the n_rigth_DCG array
        [R_2(ll,jj) , C_2(ll,jj)] = size(n_left_DCG{ll,jj}); % Find size of each matrix in the n_left_DCG array
        
    end
end
max_1 = max(C_1(:)); % Find max size of n_rigth_DCG row
max_2 = max(C_2(:)); % Find max size of n_left_DCG row

for ll=1:l
    for jj=1:j
        n_rigth_DCG_zp{ll,jj}(1,:) =  [ n_rigth_DCG{ll,jj}(1,:) ,(zeros(1 ,max_1 - C_1(ll,jj) ))]; % Zero padding of all elements  in the n_rigth_DCG
        n_left_DCG_zp{ll,jj}(1,:)  =  [ n_left_DCG{ll,jj}(1,:)  ,(zeros(1 ,max_2 - C_2(ll,jj) ))]; % Zero padding of all elements  in the n_left_DCG
        
    end
end


n_rigth_DCG_zp_mean =[];
n_left_DCG_zp_mean =[];

for jj=1:j
    n_rigth_DCG_zp_sum  = 0;
    n_left_DCG_zp_sum  = 0;
    
    for ll=1:l
           
            n_rigth_DCG_zp_sum =  n_rigth_DCG_zp_sum  + n_rigth_DCG_zp{ll,jj}(1,:) ; % Sum of all rows for each righ front
            n_left_DCG_zp_sum =   n_left_DCG_zp_sum   + n_left_DCG_zp{ll,jj}(1,:) ;  % Sum of all rows for each left front
            
           
    end
    n_rigth_DCG_zp_mean(end+1,:) =  n_rigth_DCG_zp_sum/ll; % Mean value of each rigth front
    n_left_DCG_zp_mean(end+1,:) =  n_left_DCG_zp_sum/ll;   % Mean value of each left front
    
end

flip_n_left_DCG_zp_mean = fliplr(n_left_DCG_zp_mean);
n_DCG_mean_fronts = horzcat(flip_n_left_DCG_zp_mean, n_rigth_DCG_zp_mean);

nDCG_mean_DMQIR = mean(n_DCG_mean_fronts ,1); 


size_of_nDCG_mean_DMQIR = size(nDCG_mean_DMQIR,2);
% Area under nDCG Curve
x = 1:size_of_nDCG_mean_DMQIR;   % Square brackets waste time here only
Area_Under_nDCG_mean_DMQIR = trapz(x, nDCG_mean_DMQIR); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
