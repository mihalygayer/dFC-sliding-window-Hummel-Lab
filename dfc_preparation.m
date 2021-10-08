% Mihály Gayer, semester student, Hummel Lab 2021.07.16, (Unige MSc Statistics). 
% Using parts of dFC code from Dr. Enrico Amico, cited in the code for 2 occasions

% 1. Preparation of the time series data from CONN Toolbox. Load data and transform time series into sliding window analysis form
% - commented, no need to run if TS is already saved
% 2. Load Time series data. 
% - Run1 (Ses1-Run1) 
% - Run2 (Ses1-Run2) and 
% - Run3 (Ses2-Run1)
% 3. Check stationarity inside the sliding window 
% 4. 

%% 1. Preparation of the time series data from CONN Toolbox
% addpath('/home/mgayer/matlab/conn')
% addpath('//opt/matlab/toolbox/spm12')
% addpath('//opt/matlab/toolbox/spm12/toolbox')
% addpath('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing')
% %% Use only 264 ROI names
% Index = find(contains(names,'5mm')); % Manon's comment, 2021.05.12
% Index=Index(1:264); % do not include rXYZ r: registered to MNI space, hence 264 ROIs. Manon's comment, 2021.05.12
% names=names(Index)
% save names.mat names
%%
% %conn
% 
% %% Load data
% % %% COND_Subject001_Session002
% % load('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/data/COND_Subject001_Session002.mat')
% % %% DATA_Subject001_Session002 - V structure
% % load('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/data/DATA_Subject001_Session002.mat')
% clearvars; clc
% 
% % ROI_Subject001_Session002
% %%
% cd '/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing'
% subjects=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37];
% 
% %% Example of one subject one run
%  load('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing/ROI_Subject001_Condition001.mat' )
% 
% %%
% save
%  %%
% % just use 'data'
% % d1data and d2data contain the first and the second order derivatives of your time series
% % https://www.nitrc.org/forum/message.php?msg_id=19444 
% 
% % use the condition weights only when it is a block design experiment
%     % https://www.nitrc.org/forum/message.php?msg_id=15732 % only for block design experiments
%     % https://www.nitrc.org/forum/message.php?msg_id=15735
% 
%  %%
%  load('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing/ROI_Subject001_Condition001.mat' )
%  data1=data;
%  load('/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing/ROI_Subject001_Condition002.mat' )
% data2=data;
% %% data contains the 3 time series from the three sessions, it is sufficient to load only one condition and use the corresponding time series for each runs
% sum(data1{10}==data2{10})
% sum(data1{100}==data2{100})
% 
% %% Condition 1 - meaning the 1st run out of the three in the preprocessed data
% TS_1=zeros(540,264,37); % number of time points, ROIs, subjects
% TS_2=TS_1; TS_3=TS_1;
% for subj_number=1:37
%     subj_id=subjects(subj_number);
%     
%     if subj_number <10
%        filename=sprintf('ROI_Subject00%d_Condition001.mat', subj_id);
%        load(['/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing/' filename])
%     else
%        filename=sprintf('ROI_Subject0%d_Condition001.mat', subj_id);
%        load(['/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/results/preprocessing/' filename])
%     end
% 
%     % https://www.nitrc.org/forum/message.php?msg_id=15732 % only for block design experiments
%     % https://www.nitrc.org/forum/message.php?msg_id=15735
% 
%     % RUN 1
%      idx = 1:540; 
%      data_reduced = cellfun(@(x)x(idx,:), data, 'uni',0);
%      Index = find(contains(names,'5mm')); % Manon's comment, 2021.05.12
%      Index=Index(1:264); % do not include rXYZ r: registered to MNI space, hence 264 ROIs. Manon's comment, 2021.05.12
%      data_reduced=data_reduced(Index);
%    for i=1:size(data_reduced,2)
%       TS_1(:, i, subj_number)=data_reduced{i};
%    end
%    
%        % RUN 2
%      idx = 541:1080; 
%      data_reduced = cellfun(@(x)x(idx,:), data, 'uni',0);
%      Index = find(contains(names,'5mm'));
%      Index=Index(1:264); % do not include rXYZ r: registered to MNI space, hence 264 ROIs
%      data_reduced=data_reduced(Index);
%    for i=1:size(data_reduced,2)
%       TS_2(:, i, subj_number)=data_reduced{i};
%    end
%    
%        % RUN 3
%      idx = 1081:1620; 
%      data_reduced = cellfun(@(x)x(idx,:), data, 'uni',0);
%      Index = find(contains(names,'5mm'));
%      Index=Index(1:264); % do not include rXYZ r: registered to MNI space, hence 264 ROIs
%      data_reduced=data_reduced(Index);
%    for i=1:size(data_reduced,2)
%       TS_3(:, i, subj_number)=data_reduced{i};
%    end
%    
% end
%% 2. Load Time series data
%% Save Load
% save 'TS_res_ses1_run1.mat' TS_1
% save 'TS_res_ses1_run2.mat' TS_2
% save 'TS_res_ses2_run1.mat' TS_3
% save 'ROI_Index_to_keep.mat' Index 
% save 'ROI_names.mat' names
% addpath '/home/mgayer/Sliding window analysis'
load TS_res_ses1_run1.mat
load TS_res_ses1_run2.mat
load TS_res_ses2_run1.mat

load ROI_Index_to_keep.mat
load ROI_names.mat

%% 3. Check stationarity
% clc;
% w_start=101;
% wSize = 60; 
% adftest_results=zeros(1,264); % 264 ROI
% for i=1:264
%     adftest_results(i)=adftest(TS_1(w_start:(w_start+wSize),i,1));
% end
% disp(1-sum(adftest_results)/264)
% %% Autocorrelation inside a random window - for one subject
% wSize = 60; 
% w_start = 1; 
% %autocorr(TS_1(w_start:w_start+wSize,1,1))

%% 3.
%% For all subjects and one Run:
% setup, based on one subject
wSize = 75; 
w_start = 1; 
wStep = 5;  % i.e. ~ jump in layer windows
clc;
subj_no=1;
TS=TS_2(:,:,subj_no); % Enrico's code: TS is 1200x360
N = size(TS,2); % 264: number of ROIs
numTP = size(TS,1);
% sliding window dFC
winit_last = numTP-wSize+1;
numLayers = length(1:wStep:winit_last); % number of sliding windows in total
mask_ut = triu(true(N),1); % ones on the upper triangular, zero otherwise

%% TS_1
% 34716 is the vector form of a FC matrix. iOne subject has numLayers times a FC vector
TS1_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*37); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
for subj_no=1:37 % repeat for 37 subjects
    TS=TS_1(:,:,subj_no); % Enrico's code: TS is 1200x360
    N = size(TS,2); % 264: number of ROIs
    numTP = size(TS,1);
    % sliding window dFC
    winit_last = numTP-wSize+1;
    numLayers = length(1:wStep:winit_last); % number of sliding windows in total

    mask_ut = triu(true(N),1); % ones on the upper triangular, zero otherwise
    dFCw_2D = nan(nnz(mask_ut),numLayers);
    j=1;

    %disp('Computing dFC on layer:')
    % For 1 subject compute each sliding window layered FC
    for l=1:numLayers 
        %disp(l)
        winit = ((l-1)*wStep)+1;
        wend = winit + wSize-1;
        ts_segment = TS(winit:wend,:);
        fc_segment = corr(ts_segment);                   
        fc_segment = (fc_segment + fc_segment')./2;
        fc_segment(logical(eye(size(fc_segment)))) = 0;
        fc_segment(isnan(fc_segment)) = 0;          
        dFCw_2D(:,l) = fc_segment(mask_ut); % transform into vector           
    end
    index_start=1+(subj_no-1)*numLayers; % storage for aggregating the information for each subject
    index_end=subj_no*numLayers;
    TS1_aggregate_dFCw_2D(:,index_start:index_end)=dFCw_2D;
end

%% TS_2
TS2_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*37); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
for subj_no=1:37 % repeat for 37 subjects
    TS=TS_2(:,:,subj_no); % Enrico's code: TS is 1200x360
    N = size(TS,2); % 264: number of ROIs
    numTP = size(TS,1);
    % sliding window dFC
    winit_last = numTP-wSize+1;
    numLayers = length(1:wStep:winit_last); % number of sliding windows in total

    mask_ut = triu(true(N),1); % ones on the upper triangular, zero otherwise
    dFCw_2D = nan(nnz(mask_ut),numLayers);
    j=1;

    %disp('Computing dFC on layer:')
    % For 1 subject compute each sliding window layered FC
    for l=1:numLayers 
        %disp(l)
        winit = ((l-1)*wStep)+1;
        wend = winit + wSize-1;
        ts_segment = TS(winit:wend,:);
        fc_segment = corr(ts_segment);                   
        fc_segment = (fc_segment + fc_segment')./2;
        fc_segment(logical(eye(size(fc_segment)))) = 0;
        fc_segment(isnan(fc_segment)) = 0;          
        dFCw_2D(:,l) = fc_segment(mask_ut); % transform into vector           
    end
    index_start=1+(subj_no-1)*numLayers; % storage for aggregating the information for each subject
    index_end=subj_no*numLayers;
    TS2_aggregate_dFCw_2D(:,index_start:index_end)=dFCw_2D;
end

%% TS_3
TS3_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*37); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
for subj_no=1:37 % repeat for 37 subjects
    TS=TS_3(:,:,subj_no); % Enrico's code: TS is 1200x360
    N = size(TS,2); % 264: number of ROIs
    numTP = size(TS,1);
    % sliding window dFC
    winit_last = numTP-wSize+1;
    numLayers = length(1:wStep:winit_last); % number of sliding windows in total

    mask_ut = triu(true(N),1); % ones on the upper triangular, zero otherwise
    dFCw_2D = nan(nnz(mask_ut),numLayers);
    j=1;

    %disp('Computing dFC on layer:')
    % For 1 subject compute each sliding window layered FC
    for l=1:numLayers 
        %disp(l)
        winit = ((l-1)*wStep)+1;
        wend = winit + wSize-1;
        ts_segment = TS(winit:wend,:);
        fc_segment = corr(ts_segment);                   
        fc_segment = (fc_segment + fc_segment')./2;
        fc_segment(logical(eye(size(fc_segment)))) = 0;
        fc_segment(isnan(fc_segment)) = 0;          
        dFCw_2D(:,l) = fc_segment(mask_ut); % transform into vector           
    end
    index_start=1+(subj_no-1)*numLayers; % storage for aggregating the information for each subject
    index_end=subj_no*numLayers;
    TS3_aggregate_dFCw_2D(:,index_start:index_end)=dFCw_2D;
end
%% Per subject connectivity changes in run1 run2 and run3
run_no=1;
subj_number=1;
n_clusters = 5;

if run_no==1
  subj_related_layers=zeros(size(TS1_aggregate_dFCw_2D,1), numLayers, 37);
  TS_aggregate_dFCw_2D=TS1_aggregate_dFCw_2D;
end

if run_no==2
    subj_related_layers=zeros(size(TS2_aggregate_dFCw_2D,1), numLayers, 37);
    TS_aggregate_dFCw_2D=TS2_aggregate_dFCw_2D;
end

if run_no==3
    subj_related_layers=zeros(size(TS3_aggregate_dFCw_2D,1), numLayers, 37);
    TS_aggregate_dFCw_2D=TS3_aggregate_dFCw_2D;
end
    
for subj_no=1:37 % access subject related layers and reformat
    index_start=1+(subj_no-1)*numLayers;
    index_end=subj_no*numLayers;
    subj_related_layers(:,:,subj_no)=TS_aggregate_dFCw_2D(:,index_start:index_end);
end
%% apply k means
aggregate_dFCw_2D=subj_related_layers(:,:,subj_number);
[idx, C, sumd] = kmeans(aggregate_dFCw_2D',n_clusters);
%% Visualise results
close all; t=tiledlayout(2, ceil(n_clusters/2+1));
for i=1:n_clusters
    FC_aux = zeros(N,N);
    FC_aux(mask_ut) = mean(aggregate_dFCw_2D(:,idx==i),2);
    FC_aux = FC_aux + FC_aux';
    
    nexttile
    imagesc(FC_aux); axis square; c = colorbar; %c.Label.String = 'Pearson''s correlation'; 
    title(['dFC Cluster ' int2str(i)]); colormap jet;
    xlabel('brain regions'); ylabel('brain regions'); caxis([-.7 .7]); 
    set(gcf,'color','w');
end
nexttile
scatter(1:n_clusters,sumd); xlabel('cluster number'); ylabel('Sum of k-means distance by cluster');
nexttile
hist(idx); xlabel('brain state'); ylabel('dwelling time (TR)'); title('Amount of time spent in each state')
nexttile
plot(1:numLayers, idx); ylabel('brain state number'); title('At what time points in which brain state');
title(t, sprintf('dFC sliding window states of subject %d with window length %d TRs (run %d)',subj_no, wSize, run_no))
%% By subject reformat the functional connectomes and calculate variability by brain region pairs
%% Choose which run
run_no=2;

numLayers = length(1:wStep:winit_last); % number of sliding windows in total
N = size(TS_1,2); % 264: number of ROIs
numTP = size(TS_1,1);
winit_last = numTP-wSize+1;
FCs=zeros(N,N,numLayers,37);
std_FC=zeros(N,N,37); % 264x264 ROIs

for subj_no=1:37
    for l=1:numLayers 

        if run_no==1; TS=TS_1(:,:,subj_no); end
        if run_no==2; TS=TS_2(:,:,subj_no); end
        if run_no==3; TS=TS_3(:,:,subj_no); end
        
        winit = ((l-1)*wStep)+1;
        wend = winit + wSize-1;
        
        ts_segment = TS(winit:wend,:);
        fc_segment = corr(ts_segment);                   
        fc_segment = (fc_segment + fc_segment')./2;
        fc_segment(logical(eye(size(fc_segment)))) = 0;
        fc_segment(isnan(fc_segment)) = 0; 
        
        FCs(:,:,l,subj_no)=fc_segment;
    end
    for i=1:N
        for j=1:N
            std_FC(i,j,subj_no)=std(FCs(i,j,:,subj_no));
        end
    end
end
%% uncomment the corresponding one
% overall_std_FC_run1=mean(std_FC,3);
% overall_std_FC_run2=mean(std_FC,3);
% overall_std_FC_run3=mean(std_FC,3);
% save overall_std_FC_run1.mat overall_std_FC_run1
% save overall_std_FC_run2.mat overall_std_FC_run2
% save overall_std_FC_run3.mat overall_std_FC_run3
%%
%histogram(overall_std_FC(:,:))
%% For one subject
subject_no=1
close all;
imagesc(std_FC(:,:,subj_no)); axis square; c = colorbar; c.Label.String = 'Std of the pearson''s correlation'; colormap jet;
    xlabel('brain regions'); ylabel('brain regions'); %caxis([0.2 0.307]); 
    title('Functional connectome variability over time (std.)');
    set(gcf,'color','w');
%% For all subjects
close all;
imagesc(overall_std_FC); axis square; c = colorbar; c.Label.String = 'Pearson''s correlation std.'; colormap jet;
    xlabel('brain regions'); ylabel('brain regions'); caxis([0.2 0.307]); 
    set(gcf,'color','w');
%% List those ROI pairs with highest connectivity variability
std_thr=0.3;
[roi_row,roi_column, values]=find(overall_std_FC>std_thr);
% Run 1
% high_var_roipairs_run1=[string(names(roi_row)); string(names(roi_column))]'; 
% save high_var_roipairs_run1.mat high_var_roipairs_run1

% Run2
%high_var_roipairs_run2=[string(names(roi_row)); string(names(roi_column))]';
%save high_var_roipairs_run2.mat high_var_roipairs_run2

% Run3
% high_var_roipairs_run3=[string(names(roi_row)); string(names(roi_column))]';
% save high_var_roipairs_run3.mat high_var_roipairs_run3

%% variability in the 

% Auditory =          1:13
% Cereb_cing =       14:32
% DMN =              33:89
% Dorsal =           90:100
% Frontop =         101:125
% Memory =          126:130
% Salience          131:148
% sensory_somato    149:183
% Subcortical       184:196
% Uncertain         197:224
% Ventral           225:233
% Visual            234:264
network=[repelem(1,13) repelem(2,19) repelem(3,57) repelem(4,11) repelem(5,25) repelem(6,5) repelem(7,18) repelem(8,35)...
         repelem(9,13) repelem(10,28) repelem(11,9) repelem(12, 31)];
between_FC_var_run1=zeros(12); % number of networks in the Power atlas
between_FC_var_run2=zeros(12);
between_FC_var_run3=zeros(12);
for i=1:12
    for k=1:12
    between_FC_var_run1(i,k)=mean(overall_std_FC_run1(network==i, network==k),'all');
    between_FC_var_run2(i,k)=mean(overall_std_FC_run2(network==i, network==k),'all');
    between_FC_var_run3(i,k)=mean(overall_std_FC_run3(network==i, network==k),'all');
    end
end
%% Visualise between network connectivity variability
xtick1=1:12;
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
ytick1=xtick1; yticklabels1=xticklabels1;
tiledlayout(1,3)
% Run1
nexttile
imagesc(between_FC_var_run1); axis square; %c = colorbar; c.Label.String = 'Std of the pearson''s correlation';
    colormap (flipud(parula));
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    xlabel('networks'); %ylabel('networks'); %caxis([0.2 0.307]); 
    title('Run 1 - Functional connectome variability over time (std.)');
    set(gcf,'color','w');

% Run2
nexttile
imagesc(between_FC_var_run2); axis square; %c = colorbar; c.Label.String = 'Std of the pearson''s correlation';
    colormap (flipud(parula));
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    xlabel('networks'); %ylabel('networks'); %caxis([0.2 0.307]); 
    title('Run 2 - Functional connectome variability over time (std.)');
    set(gcf,'color','w');


% Run3
nexttile
imagesc(between_FC_var_run3); axis square; c = colorbar; c.Label.String = 'Std of the pearson''s correlation';
    colormap (flipud(parula));
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    xlabel('networks'); %ylabel('networks'); %caxis([0.2 0.307]); 
    title('Run 3 - Functional connectome variability over time (std.)');
    set(gcf,'color','w');
%%
tiledlayout(1,2)
dif_2_1=((between_FC_var_run2-between_FC_var_run1)./between_FC_var_run1)*100;
nexttile
imagesc(dif_2_1); axis square; c = colorbar; c.Label.String = '% difference in std.';
    colormap (flipud(parula));
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    xlabel('networks'); %ylabel('networks'); %caxis([0.2 0.307]); 
    title('Difference of FC variability between Run 2 and Run 1','fontsize', 12);
    set(gcf,'color','w');

dif_3_1=((between_FC_var_run3-between_FC_var_run1)./between_FC_var_run1)*100;
nexttile
imagesc(dif_3_1); axis square; c = colorbar; c.Label.String = '% difference in std.';
    colormap (flipud(parula));
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    xlabel('networks'); %ylabel('networks'); %caxis([0.2 0.307]); 
    title('Difference of FC variability between Run 3 and Run 1', 'fontsize', 12);
    set(gcf,'color','w');
%%
clc;
between_FC_var_run1(8,8)
between_FC_var_run3(8,8)
%%
trial(8,8)/between_FC_var_run1(8,8)*100
%%
dif_3_1(8,8)
%% /between_FC_var_run1*100
pek=(between_FC_var_run3-between_FC_var_run1)./between_FC_var_run1*100
%% BY SAMPLE
%% aggregate information from First two Runs, TS_1 and TS_2
aggregate_dFCw_2D=[TS1_aggregate_dFCw_2D TS2_aggregate_dFCw_2D]; % put two tables next to each other horizontally

%% k-means clustering of the dFCs (as in Allen et al., 2014)   
n_clusters = 6;
[idx, C, sumd] = kmeans(aggregate_dFCw_2D',n_clusters, 'correlation'); % k means take input 91x34716 <- number of window layers x no_subjects X number of entries in the matrix
%% Visualise the brain states obtained from k means algorithm
close all;
figure, 
for i=1:n_clusters
    FC_aux = zeros(N,N);
    FC_aux(mask_ut) = mean(aggregate_dFCw_2D(:,idx==i),2);
    FC_aux = FC_aux + FC_aux';
    subplot(2,3,i),
    imagesc(FC_aux); axis square; c = colorbar; c.Label.String = 'Pearson''s correlation'; title(['dFC Cluster ' int2str(i)]); colormap jet;
    xlabel('brain regions'); ylabel('brain regions'); caxis([-.7 .7]); 
    set(gcf,'color','w');
end
suptitle(sprintf('dFC sliding window states with window length %d TRs',wSize))
%% number of windows per state
hist(idx)
[sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) sum(idx==5)] % 5 clusters
%% 6 clusters
[sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) sum(idx==5) sum(idx==6)] 
%% 
clc;
idx_TS1=idx(1:size(TS1_aggregate_dFCw_2D,2)); % number of windows per run
idx_TS2=idx(1+size(TS1_aggregate_dFCw_2D,2):end); % number of windows per run
% 5 clusters
distr_runs_per_state=[sum(idx_TS1==1) sum(idx_TS1==2) sum(idx_TS1==3) sum(idx_TS1==4) sum(idx_TS1==5);
                      sum(idx_TS2==1) sum(idx_TS2==2) sum(idx_TS2==3) sum(idx_TS2==4) sum(idx_TS2==5)]
%%
% 6 clusters
 distr_runs_per_state=[sum(idx_TS1==1) sum(idx_TS1==2) sum(idx_TS1==3) sum(idx_TS1==4) sum(idx_TS1==5) sum(idx_TS1==6);
                       sum(idx_TS2==1) sum(idx_TS2==2) sum(idx_TS2==3) sum(idx_TS2==4) sum(idx_TS2==5) sum(idx_TS2==6)]
%%
% 2 pm Friday
% Research question - fine tune
% Julie included page 13 6 networks, but 
% Clustering based on individual level or group level
% Write an email to Julia
% Motor learning literature
% Sensory motor within and between ghe others
% plots of inter 

%%%%%%%%%%%%%%%%%%%
%% MY CODE
%% Static FC overall time series
FC_overall=zeros(numel(Index));

for i=1:numel(Index)
    for j=1:numel(Index)
        series1=timeseries{i};
        series2=timeseries{j};
        FC_overall(i,j)=corr(series1,series2);
    end
end
%%

imagesc(FC_overall, [-1,1]); colormap jet; axis square; colorbar;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    
%%

%% Sliding window analysis - FIX window length

W=30;
step=2;
overall_length=540;
num_of_windows=floor((overall_length-W)/step);
SWA_FC=zeros(numel(Index),numel(Index),num_of_windows);
for k=0:num_of_windows-1

    
        for i=1:numel(Index)
                  for j=1:numel(Index)
                      series1=timeseries{i};
                      series2=timeseries{j};
        
                      series1=series1(1+k*step:W+k*step);
                      series2=series2(1+k*step:W+k*step);
                      FC(i,j)=corr(series1,series2);
                  end
        end
        
        SWA_FC(:,:,k+1)=FC;
end
%% Display FC in sliding window t
t=20;
imagesc(SWA_FC(:,:,t),[-1,1]); colormap jet; axis square; colorbar;
    set(gca,'xtick',[]); set(gca,'ytick',[]); 
    title(sprintf('Functional connectome of subject 1 in sliding window %d',t)); 
%%

W=400;
step=2;
overall_length=540;
num_of_windows=floor((overall_length-W)/step);
SWA_FC=zeros(numel(Index),numel(Index),num_of_windows);
window_position=[1,1+W];
c=1;
%%
while window_position(2)<overall_length
  
  window_position=[1+c*step, 1+c*step+W]; % [from, to]
  weights=zeros(1,overall_length); 
  weights(window_position(1):window_position(2))=1;
  
  idx = find( weights>0 );
  window_timeseries= cellfun(@(x)x(idx,:), timeseries, 'uni',0);
  
  test_results_this_window=zeros(1,numel(window_timeseries));
  
  for i=1:numel(window_timeseries)
  test_results_this_window(i)=adftest(window_timeseries{i});
  end
  disp(sum(test_results_this_window))
  FC=zeros(numel(timeseries),numel(timeseries));
  c=c+1;  
  
        for i=1:numel(timeseries)
                  for j=1:numel(timeseries)
                      FC(i,j)=corr(timeseries{i},timeseries{j});
                  end
        end
        FC=FC-diag(diag(FC)); % replace diagonal with zeros (subtraction)
        SWA_FC(:,:,k+1)=FC;
end
%% K means clustering
k=5;
for i=1:20 % update this 
    X=SWA_FC(:,:,i);
    out=kmeans(X,k);
end
%%
cluster_FC=SWA_FC(:,:,20);
trial=cluster_FC(:,:,out==1);

%%


%% https://www.nitrc.org/forum/message.php?msg_id=18633

% %% https://raphaelvallat.com/conn_ts.html
% 
% Info            = [];
% Info.wdir       = '/home/jbrancato/MasterProject/Code/conn_project_all_subjects_s2_s3_stimulation_effect/data';
% Info.session    = 1;        % 0 for all sessions, 1=session1, 2=session2, etc...
% Info.nsub       = 37;       % Total number of subjects
% 
% %%
% % ROI is the main structure containing ROI names, data
% ROI             = [];
% ROI.ROI1_data   = [];
% ROI.ROI2_data   = [];
% 
% % Select pair of ROIs
% ROI.ROI1_name   = 'NetDMN.MPFC';
% ROI.ROI2_name   = 'NetDMN.PC';
% 
% % Define outpath and outfilename
% Info.outdir     = pwd;
% Info.outfile    = [ Info.outdir '/TS_' ROI.ROI1_name  '_' ROI.ROI2_name '_RUN'
% 			num2str(Info.session) '.png' ];