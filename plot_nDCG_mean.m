

close all;
clear all;
clc;

x=linspace(0,350,1);



%{
load('/home/ubuntu/Desktop/kMQR/Evaluations_Video1_2d/nDCG_mean_kMQR_16.mat');
plot(nDCG_mean_kMQR, '-.k', 'LineWidth',3); 
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Video2_3d/nDCG_mean_kMQR_32.mat');
plot(nDCG_mean_kMQR, '-g', 'LineWidth',3); 
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Video2_3d/nDCG_mean_kMQR_64.mat');
plot(nDCG_mean_kMQR, '-b',  'LineWidth',3);
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Video2_3d/nDCG_mean_kMQR_128.mat');
plot(nDCG_mean_kMQR, '-y',  'MarkerSize',6, 'MarkerFaceColor',[0 1 0],  'LineWidth',3');
hold on


load('/home/ubuntu/Desktop/kMQR/Evaluations_Video1_2d/nDCG_mean_kMQR_256.mat');
plot(nDCG_mean_kMQR, '-r',  'LineWidth',3);


load('//home/ubuntu/Desktop/kMQR/Evaluations_Video1_2d/nDCG_mean_kMQR_512.mat');
plot(nDCG_mean_kMQR, '-k',  'LineWidth',3');
hold on
%}

load('/home/ubuntu/Desktop/kMQR/Evaluations_Streets_2d/nDCG_mean_kMQR_512.mat');
plot(nDCG_mean_kMQR, '-k',  'LineWidth',3); 
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Streets_3d/nDCG_mean_kMQR_512.mat');
plot(nDCG_mean_kMQR, '-g',  'LineWidth',3); 
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Streets_4d/nDCG_mean_kMQR_512.mat');
plot(nDCG_mean_kMQR, '-b',  'LineWidth',3);
hold on

load('/home/ubuntu/Desktop/kMQR/Evaluations_Streets_5d/nDCG_mean_kMQR_512.mat');
plot(nDCG_mean_kMQR, '-r',  'LineWidth',3');
hold on






set(gca,'FontSize',46);


%title('Average nDCG scores along 10 Pareto fronts' ,'FontSize', 32)

ylabel('Average nDCG scores.' ,'FontSize', 38)
xlabel('Number of Retrieved Items' ,'FontSize', 38) 

%legend({ '16-bit','32-bit', '64-bit', '128-bit','256-bit', '512-bit',  },'Location','southwest' ,'FontSize', 38)

legend({ '2-queries', '3-queries', '4-queries',  '5-queries',},'Location','southwest' ,'FontSize', 38)
