

close all;
clear all;
clc;

load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/hashCodes_64.mat');
data = hashCodes_64;
%load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/features/features_16.mat');
%features = features_16;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/targets.mat');
targets = targets;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/filenames.mat');
filenames = filenames;
N = length(filenames);

queryIndex = xlsread('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/qLabels_V2.xls'); 

queryIndex = transpose( queryIndex ); 
queryIndex1 = queryIndex(1,:);        % First element of Query Triplets
queryIndex2 = queryIndex(2,:);        % Second element of Query Triplets
       
maxFront = 3;  
K = 20; % Number Of Clusters

 for l = 1:240 % Number of Query Pairs 
       
     %%%%%%%%%%%%%% Creation of Pareto Points for MQR %%%%%%%%%%%%%%%%%%%%%
       q1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
       q2 = data(queryIndex2,:);
       q1_rep{l,:} = repmat(q1(l,:),N,1); % Make query matrix size to the same as data matrix size
       q2_rep{l,:} = repmat(q2(l,:),N,1);      
       xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
       xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
       hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
       hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
      
       X_MQR = zeros(2,N);
       X_MQR(1,:) = hamming_dist1{l,:};
       X_MQR(2,:) = hamming_dist2{l,:};
    
       X_MQR = (X_MQR)';
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%  Creation of Pareto Points for EMR %%%%%%%%%%%%%%%%%%%%%
        
           
        y1{l,:} = zeros(N,1);
        y1{l,:}(queryIndex1(l)) = 1; % Ranking  of query1 = 1
        y2{l,:} = zeros(N,1);
        y2{l,:}(queryIndex2(l)) = 1; % Ranking  of query2 = 1
    
    
     [H A landmarks Z] = EMRcomputeModel(data); % Compute EMR
        simEMR1{l,:} = EMRscore(H ,A, y1{l,:});  % EMR Ranking
        simEMR2{l,:} = EMRscore(H ,A, y2{l,:});  % EMR Ranking
        dist1{l,:}   = 1-simEMR1{l,:};      % Dissimilarity 
        dist2{l,:}   = 1-simEMR2{l,:};      % Dissimilarity 
      
           
        X_EMR = zeros(2,N);
        X_EMR(1,:) = dist1{l,:};
        X_EMR(2,:) = dist2{l,:};
    
        X_EMR = (X_EMR)';  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
 %%%%%%%%%%%%%%%%%%%%%%%%% Pareto Fronts for EMR %%%%%%%%%%%%%%%%%%%%%%%%%%
          pf_idx(l,:) = pareto_fronts(X_EMR, maxFront);
          
          q1_label{l,:} = targets(queryIndex1(l), : ); % Label vector of Query 1
          q2_label{l,:} = targets(queryIndex2(l), : ); % Label vector of Query 2
        
          b{l,:} = or(q1_label{l,:} , q2_label{l,:});    % beta in the equation 7 
          absolute_b{l,:} = nnz(b{l,:});                 % Number of non-zero elements in the beta, nnz is a Matlab Func.
       
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%% kMeans for MQR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [idx,C] = kmeans(X_MQR,K);  % Applay kmeans, C:Coordinates of the Clusters, 
                                        % D:Within distance in a cluster
                                        % idx, her bir nokatnın ait olduğu
                                        % cluster no'yu verir.                                        
 

    % Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
    DATA{l,:}= cell(K,1) ;
 
  
    for i = 1:K
        DATA{l,:}{i,:} = X_MQR(idx==i, :) ; % DİKKAT, BURDA 2 CELL ARRAY VAR!
       
    end
          
    
    VAR = var(C,0,2);
    [VAR_sort , VAR_index] = sort(VAR);
    VAR_index = VAR_index(1:10, : );
    SUM = sum(sqrt(C(VAR_index,:)).^2 ,2);
    LAST = [VAR_index, SUM ];
    LAST_sort = sortrows(LAST, 2);
    Best = LAST_sort(1,1);
                                        
     
    %[C_min, Best] = min(var(C,0,2)); %VARIANCES ALONG ROWS     
    Best_Points(l,:) = DATA{l,:}(Best,:);
    
  
    rtr_idx{l,:} = ismember(X_MQR, Best_Points{l,:}(:,:),'rows'); 
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
                                        
                                        
                                        
                                        
                                        
                                        
                                        
          
    
%%%%%%%%%%%%%%%%%%%% Precision & Recall for EMR  %%%%%%%%%%%%%%%%%%%%
    
    % NOW FIND Precision values of first 10 front for each query pairs
    for j= 1:maxFront
        predicted_labels_EMR{l,j} =      targets(pf_idx{l,j}(:,3),:);
        union_of_query_labels_EMR{l,:} = or(targets(queryIndex(1,l), :), targets(queryIndex(2,l), : ));
        absolute_union_of_query_labels_EMR{l,:}(j,:) = nnz(union_of_query_labels_EMR{l,:} ); 
        
        diff_EMR{l,j} = ismember( predicted_labels_EMR{l,j}, union_of_query_labels_EMR{l,:}  , 'rows'); 
        if isempty(diff_EMR{l,j})
            diff_EMR{l,j} = 0;
        end
        num_nz_EMR(l,j) = nnz( diff_EMR{l,j}(:,1) );
        s_EMR{l,j} = size(diff_EMR{l,j}(:,1), 1);
        
        for jj=1:s_EMR{l,j};
        
        % Cummulative sum of the true-positive elements
        CUMM_EMR{l,j} = cumsum(diff_EMR{l,j});          
        Precision_AT_K_EMR{l,j}(jj,1) = ( CUMM_EMR{l,j}(jj,1)  ) / jj;              
        Recall_AT_K_EMR{l,j}{jj,1} = ( CUMM_EMR{l,j}(jj,1)  ) / (num_nz_EMR(l,j)); %?????????????                    
       
        end  
        
        acc_EMR(l,j) = num_nz_EMR(l,j) / s_EMR{l,j};   % accuracy of the best cluster 
        %avg_Precision_EMR(l,j) = sum(Precision_AT_K_EMR{l,j}(:,1) ) / s_EMR{l,j};
        avg_Precision_EMR(l,j) = sum(Precision_AT_K_EMR{l,j}(:,1) .* diff_EMR{l,j} ) / num_nz_EMR(l,j);
        avg_Precision_EMR(isnan(avg_Precision_EMR))=0;

        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% nDCG scores for EMR           %%%%%%%%%%%%%%%%%%
        
         Labels{l,j} = targets(pf_idx{l,j}(:,3),:);
                                     
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
                     
                     rtr_idx{1,1}(end+1,:) = pf_idx{l,j}( ff , 3);
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
                    
             
             
                    left_DG{l,j}(:,1)   =   left_G{l,j}(:,1);
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
                 
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
             
        
    end  
    
        mean_avg_Precision_EMR(l,:) = sum(avg_Precision_EMR(l,:)) / maxFront;
        mean_acc_EMR(l,:) = sum(acc_EMR(l,:)) / maxFront;
       
 end
 
mAP_kMQR = sum(avg_Precision_MQR(:,1)) / l;
avg_acc_kMQR = sum(acc_MQR(:,1)) / l;

mAP_EMR = sum(mean_avg_Precision_EMR(:,1)) / l;
avg_acc_EMR = sum(mean_acc_EMR(:,1)) / l;


