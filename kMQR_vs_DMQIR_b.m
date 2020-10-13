% https://www.mathworks.com/matlabcentral/answers/361878-getting-the-data-points-of-each-cluster-in-kmeans
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
queryIndex = xlsread('//home/ubuntu/keras/enver/dmlvh2/DMQVR/qLabels_1.xls');   % Reads randomly choosen query pairs from excell file
queryIndex = transpose( queryIndex ); 
queryIndex1 = queryIndex(1,:);        % First element of Query Pair
queryIndex2 = queryIndex(2,:);        % Second element of Query Pair
       
maxFront = 1;  
K = 13; % Number Of Clusters

 for l = 1:250 % Number of Query Pairs 
       
     %%%%%%%%%%%%%% CREATION OF PARETO POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
       q_1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
       q_2 = data(queryIndex2,:);
       q1_rep{l,:} = repmat(q_1(l,:),N,1); % Make query matrix size to the same as data matrix size
       q2_rep{l,:} = repmat(q_2(l,:),N,1);      
       xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
       xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
       hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
       hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
      
       X = zeros(2,N);
       X(1,:) = hamming_dist1{l,:};
       X(2,:) = hamming_dist2{l,:};
    
       X = (X)';
       
     %%%%%%%%%%%%%%%%%%%% Pareto Fronts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          pf_idx(l,:) = pareto_fronts(X, maxFront);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



     %%%%%%%%%%%%%%%%%%%%%%%%%%%%% KMEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [idx, C] = kmeans(X,K) ;  % Applay kmeans, C:Coordinates of the Clusters, 
                                        % D:Within distance in a cluster
                                        % idx, her bir nokatnın ait olduğu
                                        % cluster no'yu verir.                                        
 

    % Get Ecah Clusters NOTE TRY TO GET IMAGE INDEXES IN C ALSO!
    DATA{l,:}= cell(K,1) ; 
    
    for i = 1:K
        DATA{l,:}{i,:} = X(idx==i, :) ; % DİKKAT, BURDA 2 CELL ARRAY VAR!       
    end  

    Cmn = [0,0];
    D(l,:) = pdist2(C, Cmn, 'euclid');
  
    % Find Closest index of Cluster, which is the first element of C_index, to the Cmn
    [out(l,:), C_index(l,:)] = sort(D(l,:)); 
    Best(l,:) = C_index(l,1);  % Best(l,:) is ACTUALLY Best{l,:} 
    % Points in the Best Cluster(Nearest cluster to Cmn)
    Best_Points(l,:) = DATA{l,:}(Best(l,:),:);
  
    rtr_idx{l,:} = ismember(X, Best_Points{l,:}(:,:),'rows'); 
    % Find indexes of 1 values, means we find Best_Points indexes in the X
    rtr_idx2{l,:} = find(rtr_idx{l,:}); 
    
    YY{l,:} = X(rtr_idx2{l,:},:); 
    % DD is the distance between each retrieved items to Cmn
    DD{l,:} = pdist2(YY{l,:}, Cmn, 'euclid'); 
    DD_new{l,:} = [rtr_idx2{l,:},DD{l,:}];
    DD_new_sort{l,:} = sortrows(DD_new{l,:},2);
    
    Retrieved_Items_MQR{l,:} = DD_new_sort{l,:}(:,1);    
    
    predicted_labels_MQR{l,:} = targets(Retrieved_Items_MQR{l,:},:);
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
        predicted_labels_DMQIR{l,j} =      targets(pf_idx{l,j}(:,3),:);
        union_of_query_labels_DMQIR{l,:} = or(targets(queryIndex(1,l), :), targets(queryIndex(2,l), : ));
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

for l=1:200 
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
 for l = 1:200                 % Number of Query Pairs
              
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
