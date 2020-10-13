clear all;
close all;
clc;

load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/hashCodes_256.mat');
data = hashCodes_256;
load('/home/ubuntu/Desktop/Thesis_Follow_Up_2/dmqRetrieval/lamdaDataset/hashCodes/targets.mat');
targets = targets;
load('lamdaDataset/hashCodes/filenames.mat');
filenames = filenames;

    N = 2000;           % Number of samples in the Lamda Dataset
   
    queryIndex = xlsread('qLabels_V2.xls');  % Reads randomly choosen query pairs from excell file
    queryIndex = transpose( queryIndex ); 
    %queryIndex1 = queryIndex(1,:);        % First element of Query Pair
    %queryIndex2 = queryIndex(2,:);        % Second element of Query Pair
    queryIndex1 = 600;
    queryIndex2 = 1699;     
    
    
   l = 1;
              
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
        
         plot(X(:,1),X(:,2), 'k*','LineWidth',1); 
         set(gca,'FontSize',20);

         hold on;
