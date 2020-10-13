clear all;
close all;
clc;

load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/hashCodes_128.mat');
data = hashCodes_128;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/features/features_128.mat');
features = features_128;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/targets.mat');
targets = targets;
load('/home/ubuntu/keras/enver/dmlvh2/DMQVR/myDataset/hashCodes/filenames.mat');
filenames = filenames;

queryIndex1 = 11;
queryIndex2 = 26;
N = length(filenames); 

    l=1;
    

    
            
        q1 = data(queryIndex1,:);         % q1 & q2 are query pairs in the loop
        q2 = data(queryIndex2,:);
        q1_rep{l,:} = repmat(q1(l,:),N,1); % Make query matrix size to the same as data matrix size
        q2_rep{l,:} = repmat(q2(l,:),N,1);      
        xor_data_q1new{l,:} = xor(data, q1_rep{l,:}); % xor of data and query matrices
        xor_data_q2new{l,:} = xor(data, q2_rep{l,:});       
        hamming_dist1{l,:} = sum(xor_data_q1new{l,:},2); % sum up rows to get hamming distances
        hamming_dist2{l,:} = sum(xor_data_q2new{l,:},2);
        norm_hamming_dist1{l,:} =  hamming_dist1{l,:} / max( hamming_dist1{l,:}(:) ); % Normalize hamming  distances between 0&1
        norm_hamming_dist2{l,:} =  hamming_dist2{l,:} / max( hamming_dist2{l,:}(:) );
       
        
        X = zeros(2,N);
        X(1,:) = hamming_dist1{l,:};
        X(2,:) = hamming_dist2{l,:};
    
        X = (X)';
        plot(X(:,1),X(:,2),  'o','MarkerFaceColor',[0 1 0],'MarkerSize',15)    
        set(gca,'FontSize',30); 
       
        maxFront = 1;
        
       [pf_idx] = pareto_fronts(X, maxFront);
        
          
       
      
       %{
       h(2)= plot(pf_idx{1,1}(:,1), pf_idx{1,1}(:,2), '-rs','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor',[1 0 0], 'LineWidth',3); 

       h(3)= plot(pf_idx{2,1}(:,1), pf_idx{2,1}(:,2), '-bo','MarkerSize',16,'MarkerEdgeColor','blue','MarkerFaceColor',[0 0 1], 'LineWidth',3); 

       h(4)= plot(pf_idx{3,1}(:,1), pf_idx{3,1}(:,2), '-gd','MarkerSize',18,'MarkerEdgeColor','green','MarkerFaceColor',[0 1 0], 'LineWidth',3);
       
       set(gca,'FontSize',46);  
        
       legend(h([2 3 4]),{'1st Pareto Front','2nd Pareto Front','3rd Pareto Front',  },'Location','southwest', 'FontSize', 34);
       %}
       
            
      
       ylabel('d_1 ', 'FontSize', 50)
       xlabel('d_2 ', 'FontSize', 50)
       
       
       
       
       
       
       
           
         
        