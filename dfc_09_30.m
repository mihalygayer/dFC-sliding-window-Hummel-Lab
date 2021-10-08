% Mihály Gayer, semester student, Hummel Lab 2021.10.16, (Unige MSc Statistics). 
% Code prepared individually by me. Using small parts of dFC code from Dr. Enrico Amico, cited in the code for 2 occasions
% Main goal: - Choose window length and load time series. 
%            - Then create windowed layered correlations, for each of the 3 runs
%            - calculate FCV then carry out Sliding Window Analysis

% Folder and files: 'dFC' folder in Google Drive is sufficient to run every part of this code and obtain the analysis results: 
% https://drive.google.com/drive/folders/1FERWcWdkai4HKDl0Gnq5sR6CmcmMderM
%%%%%%
% 0. Use dfc_preparation.m script to obtain the time series data. Or load directly at step 1.a

% 1.a Load Time series data. Specify folder for Window Length. Each Time series for each window length is in different folder
% - Run1 (Ses1-Run1) 
% - Run2 (Ses1-Run2) and 
% - Run3 (Ses2-Run1)


% 1.b Load learning rate and intercept from behaviour data
% - log coefficients

% 1.c check stationarity of the timeseries in the windows. Not necessary,
% commented for now.

% 2.a Sliding window correlation layers. Obtain for each runs a 2D matrix
%    containing the sliding window correlation layers in vector form
% -  TS1_aggregate_dFCw_2D, with dimensions: vector length x (number of window layers * number of subjects)

% 2.b Aggregate information for clustering later. Choose which two runs.
%     ! Note: you must specify here which two runs you want to compare in the brain state analysis. 
%     Two runs are compared to one another when analysing brain states. This
%     only applies to the second section, for the brain state analysis.

% 3.a Per subject connectivity changes in run1 run2 and run3. Cluster for one subject. Visualise for one subject
% - choose a run, choose a subject, choose a n_clusters for k-means

% 3.b Obtain std_variability by subject over all the FC windows. For each run. overall_std_FC_run2.
% - first obtain 'FCs' 4 argument list with 264 x 264 x numLayers x numSubjects. With all layers of FC for all subjects per one run

% 3.c Over All FC entries: corr(FCV,Learning rate)

% 3.d How the FCV correlates with learning rate in the 12 networks? 
% - For each run obtain a 12x12 matrix with correlation of FCV with learning rate  

% 3.e Histogram of std FC (FCV) values for one subject
% - choose subject and run

% 3.f Fluctuations of FC values for one ROI pairs
% Choose subject and ROI pair (i,j)
% Change FCs_run1 to other run if want to change Run

% 3.g Within network FCV changes over the sliding windows
% - choose a network with i brain regions interval

% 3.h Static FC - For one subject

% 3.i One subject FCV matrix, 264x264 dimensions
% - choose subject
% - first one figure for one run
% - then 3 figures for each run (individual for one subject)

% 3.j Sample FCV for each runs
% - mean sample FCV in each run

% 3.k List those ROI pairs with highest connectivity variability
% - obtain: high_var_roipairs a table of ROI pairs of those with the highest FCV on the Sample level

% 3.l Sample mean FCV in the Power atlas networks (FCV uses std.)
%     Aggregate information in 12x12 matrix by taking the mean sample level variation
%     This is repeated for the 3 respective runs
%     Diagonal is within, non-diagonal represents within FC values

% 3.m Subject level FCV mean variation in the Power atlas networks
%     Obtain the same as in 3.l above, but available for each subject. In
%     the form of 12x12x35 (12x12 square matrix for the 35 subjects)


%%%%%%%%%%%%%
% - obtain 'overall_std_FC_run2'
% - histogram 1 vs all
% - Functional Connectome 1 vs all
% - Functional connectomes overall std for each run

% 3.c List those ROI pairs with highest connectivity variability
% - high_var_roipairs_run2.mat

% 3.d Visualise between network connectivity variability. Sample level, overall
% - change from one run to the next

% 4.a k-means clustering of the dFCs, aggregate level

% 4.b Visualise the number of brain states

%% 1.a Load Time series data - adapt the Path
clc; close all;

%cd 'C:\Users\user\Documents\Manon\Sliding Window Analysis\HPF_0.0222_50TR'; wSize = 50; 
addpath 'C:\Users\user\Documents\Manon\Sliding Window Analysis\HPF_0.01852_60TR'; wSize = 60; 
%cd 'C:\Users\user\Documents\Manon\Sliding Window Analysis\HPF_0.01587_70TR'; wSize = 70; 
%cd 'C:\Users\user\Documents\Manon\Sliding Window Analysis\HPF_0.01_111TR'; wSize = 111; 

load TS_res_ses1_run1.mat
load TS_res_ses1_run2.mat
load TS_res_ses2_run1.mat

load ROI_names.mat % 264 ROIs
% % check for missing data
% z=zeros(size(TS_1,3),3);
% for subj_no=1:size(TS_1,3)
% z(subj_no,1)=sum(sum(TS_1(:,:,subj_no)));
% z(subj_no,2)=sum(sum(TS_2(:,:,subj_no)));
% z(subj_no,3)=sum(sum(TS_3(:,:,subj_no)));
% end

% Exclude subj 15 and 33
TS_1=TS_1(:,:,[1:14 16:32 34:end]);
TS_2=TS_2(:,:,[1:14 16:32 34:end]);
TS_3=TS_3(:,:,[1:14 16:32 34:end]);
%
w_start=1;
wStep = 5;
%% 1.b Accuracy percent fits
%   addpath('L:\EconS_MG_2021\Data analysis\Fit - behavioural data\fits_withcoefficients')
addpath 'C:\Users\user\Documents\Manon\Sliding Window Analysis\fits_withcoefficients'
load('fits_SES12_accuracy_percent_subj.mat')
cfs=[];
for subj_no=1:37
 cf=fits_accuracy_percent_subj(subj_no).log_coeff; % Obtain log coefficients
 cfs=[cfs; cf];
end
cfs=cfs([1:14 16:32 34:37],:); % exclude subjects 15 and 33

% %% 1.c check stationarity
% % 1 subj,  1 window
% % 1 indicates that there is sufficient evidence to suggest that series is trend stationary. 
% % 0 indicates that there is not enough evidence
% 
% adftest_results=zeros(7,264); % 264 ROI
% kpss_results=zeros(7,264);
% for k=1:1%37 % for one subject, same results across all subjects
%     for j=1:7
%             if j>1; w_start=w_start+wSize; end
%         for i=1:264
%             adftest_results(j,i)=adftest(TS_2(w_start:(w_start+wSize-1),i,1));
%             kpss_results(j,i)=kpsstest(TS_2(w_start:(w_start+wSize-1),i,1));
%         end
%     end
% end
% clc; 
% disp(sum(adftest_results,2)/264)
% disp(sum(kpss_results,2)/264)

%% 2.a Sliding window correlation layers
% Parameter setup, based on one subject

clc;
subj_no=1;
TS=TS_1(:,:,subj_no); % Enrico's code: TS is 1200x360
N = size(TS,2); % 264: number of ROIs
numTP = size(TS,1);
% sliding window dFC
winit_last = numTP-wSize+1;
numLayers = length(1:wStep:winit_last); % number of sliding windows in total
mask_ut = triu(true(N),1); % ones on the upper triangular, zero otherwise

% TS_1
% 34716 is the vector form of a FC matrix. One subject has numLayers times a FC vector
TS1_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*size(TS_1,3)); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
TS1_dFCw_2D_bysubj = nan(nnz(mask_ut),numLayers,size(TS_1,3)); % same but store the 2D results into a subject by subject storage
for subj_no=1:size(TS_1,3) % repeat for 37 subjects
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
    TS1_dFCw_2D_bysubj(:,:,subj_no)=dFCw_2D;
end

% TS_2
TS2_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*size(TS_1,3)); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
TS2_dFCw_2D_bysubj = nan(nnz(mask_ut),numLayers,size(TS_1,3));
for subj_no=1:size(TS_1,3) % repeat for 37 subjects
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
    TS2_dFCw_2D_bysubj(:,:,subj_no)=dFCw_2D;
end

% TS_3
TS3_aggregate_dFCw_2D = nan(nnz(mask_ut),numLayers*size(TS_1,3)); % 37 subjects, put horizontally next to each other 37 times 34716*91 matrix. 
TS3_dFCw_2D_bysubj = nan(nnz(mask_ut),numLayers,size(TS_1,3));
for subj_no=1:size(TS_1,3) % repeat for 37 subjects
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
    TS3_dFCw_2D_bysubj(:,:,subj_no)=dFCw_2D;
end

%%%%%%%%
% 2.b Aggregate data for Clustering. Specify:
% Runs 1-2 aggregate information from First two Runs, TS_1 and TS_2
aggregate_dFCw_2D=[TS1_aggregate_dFCw_2D TS2_aggregate_dFCw_2D]; % put two tables next to each other horizontally
text_runs='Run 1 and Rund 2';
% % Runs 1-3 aggregate information from First two Runs, TS_1 and TS_2
% aggregate_dFCw_2D=[TS1_aggregate_dFCw_2D TS3_aggregate_dFCw_2D]; % put two tables next to each other horizontally
% text_runs='Run 1 and Rund 3';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.a FCV - Set subject, run and n# of clusters (5 by default)
% Per subject connectivity changes in run1 run2 and run3
run_no=2;
subj_number=1;
n_clusters = 5;

if run_no==1
    aggregate_dFCw_2D=TS1_dFCw_2D_bysubj(:,:,subj_number);
end

if run_no==2
    aggregate_dFCw_2D=TS2_dFCw_2D_bysubj(:,:,subj_number);
end

if run_no==3
    aggregate_dFCw_2D=TS3_dFCw_2D_bysubj(:,:,subj_number);
end
 
% apply k means
[idx, C, sumd] = kmeans(aggregate_dFCw_2D',n_clusters);
% Visualise results for one subject
clc; close all; t=tiledlayout(2, ceil(n_clusters/2+1));
for i=1:n_clusters
    FC_aux = zeros(N,N);
    FC_aux(mask_ut) = mean(aggregate_dFCw_2D(:,idx==i),2);
    FC_aux = FC_aux + FC_aux';
    
    nexttile
    imagesc(FC_aux); axis square; 
    title(['dFC Cluster ' int2str(i)]); colormap jet;
     caxis([-.9 .9]); 
    set(gcf,'color','w');
    if i==1; xlabel('brain regions'); ylabel('brain regions');
    end
    if i==4; c = colorbar; c.Label.String = 'Pearson''s correlation';
    end
end
nexttile % Sum of k-means distance by cluster
scatter(1:n_clusters,sumd); xlabel('cluster number'); ylabel('Sum of k-means distance by cluster');
nexttile
histogram(idx); xlabel('brain state'); ylabel('dwelling time (TR)'); title('Amount of time spent in each state')
nexttile
plot(1:numLayers, idx); ylabel('brain state number'); title('At what time points in which brain state');
title(t, sprintf('dFC sliding window states of subject %d with window length %d TRs (Run %d)',subj_number, wSize, run_no))

%% 3.b By subject reformat the functional connectomes and calculate variability by brain region pairs - Data preparation

numLayers = length(1:wStep:winit_last); % number of sliding windows in total
N = size(TS_1,2); % 264: number of ROIs
numTP = size(TS_1,1);
winit_last = numTP-wSize+1;
FCs=zeros(N,N,numLayers,size(TS_1,3));
std_FC=zeros(N,N,size(TS_1,3)); % 264x264 ROIs
std_FC_2D=zeros(length(fc_segment(mask_ut)),size(TS_1,3));
for r=1:3 % each runs
    run_no=r;
for subj_no=1:size(TS_1,3)
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
    
    % vector form for FC std, for each subject
    fc=std_FC(:,:,subj_no);
    std_FC_2D(:,subj_no)=fc(mask_ut);
end

if run_no==1; FCs_run1=FCs; std_FC_run1=std_FC; std_FC_2D_run1=std_FC_2D; overall_std_FC_run1=mean(std_FC,3); end
if run_no==2; FCs_run2=FCs; std_FC_run2=std_FC; std_FC_2D_run2=std_FC_2D; overall_std_FC_run2=mean(std_FC,3); end
if run_no==3; FCs_run3=FCs; std_FC_run3=std_FC; std_FC_2D_run3=std_FC_2D; overall_std_FC_run3=mean(std_FC,3); end
end
% uncomment if necessary to resave
% save overall_std_FC_run1.mat overall_std_FC_run1
% save overall_std_FC_run2.mat overall_std_FC_run2
% save overall_std_FC_run3.mat overall_std_FC_run3

%% 3.c Over All FC entries: corr(FCV,Learning rate)

% disp(mean(std_FC_2D_run1(:)))
% disp(mean(std_FC_2D_run2(:)))
% disp(mean(std_FC_2D_run3(:)))

for subj_no=1:size(TS_1,3)
FCV_run1(subj_no)=mean(std_FC_2D_run1(:,subj_no));
FCV_run2(subj_no)=mean(std_FC_2D_run2(:,subj_no));
FCV_run3(subj_no)=mean(std_FC_2D_run3(:,subj_no));
end

clc;
% corr(FCV_run1',cfs(:,1))
% corr(FCV_run2',cfs(:,1))
% corr(FCV_run3',cfs(:,1))
% 
% corr(FCV_run1',cfs(:,2))
% corr(FCV_run2',cfs(:,2))
% corr(FCV_run3',cfs(:,2))

% For each ROI pairs
learnrate_FCV_run1=zeros(N,N); learnrate_FCV_run2=zeros(N,N); learnrate_FCV_run3=zeros(N,N);
for i=1:N
    for j=1:N
        if j ~= i
            values=std_FC_run1(i,j,:);
            learnrate_FCV_run1(i,j)=corr(values(:),cfs(:,2));
            
            values=std_FC_run2(i,j,:);
            learnrate_FCV_run2(i,j)=corr(values(:),cfs(:,2));
            
            values=std_FC_run3(i,j,:);
            learnrate_FCV_run3(i,j)=corr(values(:),cfs(:,2));
            
        end
    end
end

clc; close all; tiledlayout(1,3); nexttile
imagesc(learnrate_FCV_run1); caxis([0.4 0.7]); colormap(flipud(gray)); xlabel('Brain regions'); ylabel('Brain regions');
title('Run 1, FCV correlation with learning rate');

nexttile
imagesc(learnrate_FCV_run2); caxis([0.4 0.7]); colormap(flipud(gray)); xlabel('Brain regions'); ylabel('Brain regions');
title('Run 2, FCV correlation with learning rate');

nexttile
imagesc(learnrate_FCV_run3); caxis([0.4 0.7]); colormap(flipud(gray)); xlabel('Brain regions'); ylabel('Brain regions');
title('Run 3, FCV correlation with learning rate');
colorbar 
set(gcf,'color','w');

% % Set Threshold for correlation values, only list those that are higher than thr
% thr=0.4; % List those 
% mrtx=learnrate_FCV_run1;
% indices = find(mrtx<0.5);
% mrtx(indices)=0;
% mrtx=max(mrtx, thr);
% mrtx(mrtx==thr)=0;
% close all;
% imagesc(mrtx)
%% 3.d How the FCV correlates with learning rate in the 12 networks? For each runs 
network=[repelem(1,13) repelem(2,19) repelem(3,57) repelem(4,11) repelem(5,25) repelem(6,5) repelem(7,18) repelem(8,35)...
         repelem(9,13) repelem(10,28) repelem(11,9) repelem(12, 31)];
between_FC_var_run1=zeros(12); % number of networks in the Power atlas
between_FC_var_run2=zeros(12);
between_FC_var_run3=zeros(12);
for i=1:12
    for k=1:12
    between_FC_var_run1(i,k)=mean(learnrate_FCV_run1(network==i, network==k),'all');
    between_FC_var_run2(i,k)=mean(learnrate_FCV_run2(network==i, network==k),'all');
    between_FC_var_run3(i,k)=mean(learnrate_FCV_run3(network==i, network==k),'all');
    end
end
%% 3.e Histogram of Functional Connectivity standard deviations for 1 subject and across all subjects, in one run
% Choose run and subject
run_no=2; 
subj_no=1;

if run_no==1; std_FC=std_FC_run1; end
if run_no==2; std_FC=std_FC_run2; end
if run_no==3; std_FC=std_FC_run3; end

close all;
tiledlayout(1,2)
% histogram of one subject's mean FC
nexttile
histogram(mean((FCs(:,:,:,subj_no)),3)); title(sprintf('Mean Functional Connectivity values of subject %d in Run %d', subj_no, run_no));
xlabel('Functional connectivity mean values for one subject'); ylabel('Occurence')

nexttile
mrtx=std_FC(:,:,subj_no);
histogram(mrtx(mask_ut)); title(sprintf('Functional Connectivity standard deviations (FCV) of subject %d in Run %d', subj_no, run_no));
xlabel('Functional connectivity standard deviation values for one subject'); ylabel('Occurence')
%% 3.f Fluctuations of FC values for one ROI pairs
% Choose subject and ROI pair (i,j)
% Change FCs_run1 to other run if want to change Run

subj_no=1;
i=167; % Sen-som. hand
j=240; % Visual

values=[FCs_run1(i,j,:,subj_no)]; values=values(:)';
close all;
plot(1:numLayers, values)
hold on
r_static=corr(TS_1(:,i,subj_no),TS_1(:,j,subj_no));
yline(r_static,'--'); text(0.7*numLayers,2.7*r_static, 'Static FC correlation', 'fontSize', 8);
title('Pairwise FC change over sliding windows between 2 ROIs (Sensory somatomotor and Visual)', 'fontSize',10)
xlabel('Sliding window number'); ylabel('Pearson correlation');
set(gcf,'color','w')

%% 3.g Within network FCV changes over the sliding windows
i=149:178; % Sen-som. hand
close all; 
t=tiledlayout(5,7)
for subj_no=1:size(TS_1,3)
%subj_no=12;
values=[FCs_run2(i,i,:,subj_no)]; %values=values(:)';
v=zeros(1,numLayers);
for l=1:numLayers
    v(l)=mean(triu(values(:,:,l)),'all');
end
nexttile
plot(1:numLayers, v)
hold on
r_static=mean(v);
yline(r_static,'--'); text(0.7*numLayers,1.3*r_static, 'Static FC', 'fontSize', 6);
%title('FC change over sliding windows, within network Sensory-somatomotor', 'fontSize',10)
ylim([0 0.4]);
title(sprintf('Subject %d',subj_no), 'fontSize', 5);
end
xlabel(t,'Sliding window number'); ylabel(t,'Pearson correlation (mean within network)');
set(gcf,'color','w')

%% 3.h Static FC - For one subject 
% (Keep in mind that this is calculated from the mean sliding window corr windows, not equal to the static FC calculated overall the time series (but close estimation) 
% static mean FC
% Change FCs_run1 to other run if you want to change Run
subject_no=1;

xtick1=[1 20 55 90 110 126 137 165 184 205 225 250];
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','M.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
ytick1=xtick1; yticklabels1=xticklabels1;

close all;
%tiledlayout(1,2); nexttile
imagesc(mean(FCs_run1(:,:,:,subj_no),3)); 
    c = colorbar; c.Label.String = 'Pearson correlation'; colormap jet; 
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('Static functional connectome of one subject','fontsize',12);
    set(gcf,'color','w'); 
    
%% 3.i One subject FCV matrix, 264x264
% Change FCs_run1 to other run if you want to change Run
subj_no=1;

close all; clc;
imagesc(mean(std_FC_run1(:,:,subj_no),3));% caxis([0.1 0.35]);
    c = colorbar; c.Label.String = 'SD of Pearson correlation'; colormap jet; 
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('FCV over all sliding windows for one subject (Std. of FC values)','fontsize',12);
    set(gcf,'color','w'); 
%% Same figure, for each of the 3 runs

mn=0; mx=0.55;

close all;
tiledlayout(1,3); nexttile
imagesc(std_FC_run1(:,:,subj_no)); axis square; %c = colorbar; c.Label.String = 'Standard deviation'; colormap jet;
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    %ylabel('brain regions'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('FCV in Run 1','fontsize',12);

nexttile
imagesc(std_FC_run2(:,:,subj_no)); axis square; %c = colorbar; c.Label.String = 'Standard deviation'; colormap jet;
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    %ylabel('brain regions'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('FCV in Run 1','fontsize',12);

nexttile
imagesc(std_FC_run3(:,:,subj_no)); axis square; 
    c = colorbar; c.Label.String = 'Standard Deviation of FC values'; colormap jet; 
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    %ylabel('brain regions'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('FCV in Run 1','fontsize',12);

%% 3.j Sample FCV for each runs - min max threshold
clc;
% min(overall_std_FC_run1(:))
% min(overall_std_FC_run2(:))
% min(overall_std_FC_run3(:))
% 
% max(overall_std_FC_run1(:))
% max(overall_std_FC_run2(:))
% max(overall_std_FC_run3(:))

% Set min and max values for visual plots
mn=0.22; mx=0.35

xtick1=[1 20 55 90 110 126 137 165 184 205 225 250];
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','M.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
ytick1=xtick1; yticklabels1=xticklabels1;

close all; 
tiledlayout(1,3);
nexttile
imagesc(overall_std_FC_run1); axis square; 
%c = colorbar; colormap jet; 
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    %ylabel('brain regions'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('Overall FCV in Run 1','fontsize',12);
    
nexttile
imagesc(overall_std_FC_run2); axis square; 
%c = colorbar;  colormap jet; 
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('Overall FCV in Run 2','fontsize',12);
    
nexttile
imagesc(overall_std_FC_run3); axis square; 
    c = colorbar; c.Label.String = 'Pearson''s correlation std.'; colormap jet;
    xlabel('ROIs, Power atlas networks'); ylabel('ROIs, Power atlas networks'); 
    caxis([mn mx]); 
    set(gcf,'color','w'); 
    set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
    title('Overall FCV in Run 3','fontsize',12);
    
%% 3.k List those ROI pairs with highest connectivity variability - sample level
% Change FCs_run1 to other run if you want to change Run
clc;

% Run2
% List those ROI pairs with FCV over 99 % of all. Most varying.
thr=prctile(overall_std_FC_run2(:),99); % 99th percentile value, set as threshold
[roi_row,roi_column]=find(triu(overall_std_FC_run2)>=thr); % Matrix is symmetric, use triu to avoid double counting
high_var_roipairs=[string(names(roi_row)); string(names(roi_column))]';
%save high_var_roipairs.mat high_var_roipairs_run2

%% STATISTICS take away: In General a Static FC Pairwise ROI value is less interpretable if has a higher variability over the sliding windows. High std. 
% For example these ROI pairs listed.


%% 3.l Sample level between and within network variability - average over all subjects (hence sample level)

% 1. Auditory =          1:13
% 2. Cereb_cing =       14:32
% 3. DMN =              33:89
% 4. Dorsal =           90:100
% 5. Frontop =         101:125
% 6. Memory =          126:130
% 7. Salience          131:148
% 8. sensory_somato    149:183
% 9. Subcortical       184:196
% 10.Uncertain         197:224
% 11.Ventral           225:233
% 12.Visual            234:264
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

%% 3.m Subj level between FCV - for all subjects
subj_no=1;
between_FC_var_run1_subj=zeros(12,12,35); % number of networks in the Power atlas
between_FC_var_run2_subj=zeros(12,12,35);
between_FC_var_run3_subj=zeros(12,12,35);
for subj_no=1:35
for i=1:12
    for k=1:12
    between_FC_var_run1_subj(i,k,subj_no)=mean(std_FC_run1(network==i, network==k, subj_no),'all');
    between_FC_var_run2_subj(i,k,subj_no)=mean(std_FC_run2(network==i, network==k, subj_no),'all');
    between_FC_var_run3_subj(i,k,subj_no)=mean(std_FC_run3(network==i, network==k, subj_no),'all');
    end
end
end
%% 3.n Learning rate ~ FCV between values
clc;
learning_FCV_between_run1=zeros(12,12); learning_FCV_between_run2=zeros(12,12); learning_FCV_between_run3=zeros(12,12);
for i=1:12
    for j=1:12
        if i~=j
            values=between_FC_var_run1_subj(i,j,:);
            learning_FCV_between_run1(i,j)=corr(values(:),cfs(:,2));
            
            values=between_FC_var_run2_subj(i,j,:);
            learning_FCV_between_run2(i,j)=corr(values(:),cfs(:,2));
            
            values=between_FC_var_run3_subj(i,j,:);
            learning_FCV_between_run3(i,j)=corr(values(:),cfs(:,2));
        end
    end
end
close all;
tiledlayout(1,3); nexttile
imagesc(learning_FCV_between_run1);  title('Run1'); caxis([-0.55 0.5]);
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);

nexttile
imagesc(learning_FCV_between_run2);  title('Run2'); caxis([-0.55 0.5]);
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);

nexttile
imagesc(learning_FCV_between_run3); colorbar; title('Run3'); caxis([-0.55 0.5]);
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
set(gcf,'color','w')
%%
close all;
imagesc(learning_FCV_between_run1);  title('Run1'); caxis([-0.55 0.5]); colorbar
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontsize',6); set(gca,'ytick',ytick1, 'yticklabel',yticklabels1);
set(gcf,'color','w')

%%
min(learning_FCV_between_run1(:))
min(learning_FCV_between_run2(:))
min(learning_FCV_between_run3(:))

max(learning_FCV_between_run1(:))
max(learning_FCV_between_run2(:))
max(learning_FCV_between_run3(:))

%% WRONG 3.d Visualise between network connectivity variability
close all;
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

    
    %% Difference between runs. Sample level connectivity variability differences
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
%% Memory retrieval within network - 6

transf=triu(true(5),1);
run1=overall_std_FC_run1(network==6, network==6); 
run1=run1(transf);

run2=overall_std_FC_run2(network==6, network==6); 
run2=run2(transf);

run3=overall_std_FC_run3(network==6, network==6); 
run3=run3(transf);

% Anova for three groups
values=[run1; run2; run3]
g=[repelem({'Run_1'},length(run1)), repelem({'Run_2'},length(run2)), repelem({'Run_3'},length(run3))]; close all;
anova1(values, g)

%% Anova for two groups Run1 and Run3
values=[run1; run3]
g=[repelem({'Run_1'},length(run1)), repelem({'Run_3'},length(run3))]; close all;
 anova1(values, g)
%% SAMPLE level Within, Between network Conn SD changes between runs - Bar graph like Elton and Gao 2015
% within:
close all

within=[diag(between_FC_var_run1) diag(between_FC_var_run2) diag(between_FC_var_run3)];
%Remove subcortical and uncertain
% xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
% xtick1=1:12;
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Vent.', 'Vis.'};
xtick1=1:11;
% xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Vent.', 'Vis.'};
% xtick1=1:10;
within=within([1:9 11 12],:)

within_std=zeros(11,3);
for i=1:11
    
    within2(i,1)=mean(between_FC_var_run1_subj(i,i,:));
    within2(i,2)=mean(between_FC_var_run2_subj(i,i,:));
    within2(i,3)=mean(between_FC_var_run3_subj(i,i,:));
    
    within_std(i,1)=std(between_FC_var_run1_subj(i,i,:));
    within_std(i,2)=std(between_FC_var_run2_subj(i,i,:));
    within_std(i,3)=std(between_FC_var_run3_subj(i,i,:));
end

b=bar(within2); set(b, {'DisplayName'}, {'Run 1','Run 2','Run 3'}'); %legend();
title('Within network connectivity variability changes across Runs');
ylabel('Connectivity variability'); xlabel('Power atlas connectivity networks')
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1);
set(gcf,'color','w');
% errorbar
[ngroups,nbars] = size(within);
% Get the x coordinate of the bars
hold on
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
errorbar(x',within,within_std,'k','linestyle','none');
hold off
%er = errorbar(x,within',errlow,errhigh); 

%% between
% make diag 0
between_FC_var_run1=between_FC_var_run1-diag(diag(between_FC_var_run1)); 
between_FC_var_run2=between_FC_var_run2-diag(diag(between_FC_var_run2)); 
between_FC_var_run3=between_FC_var_run3-diag(diag(between_FC_var_run3)); 
between=[sum(between_FC_var_run1,2) sum(between_FC_var_run2,2) sum(between_FC_var_run3,2)];

xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
xtick1=1:12;
%Remove subcortical and uncertain
% xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Vent.', 'Vis.'};
% xtick1=1:10;
% between=between([1:8 11 12],:)
close all; b=bar(between); set(b, {'DisplayName'}, {'Run 1','Run 2','Run 3'}'); legend('Location','NorthEastOutside');
title('Between network connectivity variability changes across Runs');
ylabel('Connectivity SD'); xlabel('Power atlas connectivity networks');
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1);
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behaviour ~ Subject level variability changes between runs

%%
network=[repelem(1,13) repelem(2,19) repelem(3,57) repelem(4,11) repelem(5,25) repelem(6,5) repelem(7,18) repelem(8,35)...
         repelem(9,13) repelem(10,28) repelem(11,9) repelem(12, 31)];
between_FC_var_run1=zeros(12,12,size(TS_1,3)); % number of networks in the Power atlas x number of subjects
between_FC_var_run2=zeros(12,12,size(TS_1,3));
between_FC_var_run3=zeros(12,12,size(TS_1,3));

for subj_no=1:size(TS_1,3)
    for i=1:12
        for k=1:12
            
        between_FC_var_run1(i,k, subj_no)=mean(std_FC_run1(network==i, network==k, subj_no),'all');
        between_FC_var_run2(i,k, subj_no)=mean(std_FC_run2(network==i, network==k, subj_no),'all');
        between_FC_var_run3(i,k, subj_no)=mean(std_FC_run3(network==i, network==k, subj_no),'all');
        end
    end
end
% Within
within=zeros(12,3,size(TS_1,3)); within_diff21=zeros(12,35); within_diff31=zeros(12,35);
for subj_no=1:size(TS_1,3)
within(:,:,subj_no)=[diag(between_FC_var_run1(:,:,subj_no)) diag(between_FC_var_run2(:,:,subj_no)) diag(between_FC_var_run3(:,:,subj_no))];

within_diff21(:,subj_no)=within(:,2,subj_no)-within(:,1,subj_no);

within_diff31(:,subj_no)=within(:,3,subj_no)-within(:,1,subj_no);
end

%% Between 
between=zeros(12,3,size(TS_1,3));
for subj_no=1:size(TS_1,3)
    for i=1:12
    between(i,1,subj_no)=sum(between_FC_var_run1_subj(i,:,subj_no))-between_FC_var_run1_subj(i,i,subj_no); % sum of between - within
    between(i,2,subj_no)=sum(between_FC_var_run2_subj(i,:,subj_no))-between_FC_var_run2_subj(i,i,subj_no);
    between(i,3,subj_no)=sum(between_FC_var_run3_subj(i,:,subj_no))-between_FC_var_run3_subj(i,i,subj_no);
    end 
end

between_diff21=zeros(12,35); between_diff31=zeros(12,35);
for subj_no=1:size(TS_1,3)
    between_diff21(:,subj_no)=between(:,2,subj_no)-between(:,1,subj_no);
    between_diff31(:,subj_no)=between(:,3,subj_no)-between(:,1,subj_no);
end

%% Within Run 2 Run 1
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
close all;
t=tiledlayout(2,6)
for i=1:12
    nexttile
    scatter(within_diff21(i,:), cfs(:,2), 4, 'filled')
    hold on
    lsline
    ylabel('Learning rate'); 
    [r, p] = corr(within_diff21(i,:)', cfs(:,2));
    
    if p<0.05; 
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(within_diff21(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','green','FontSize',6,'fontweight', 'bold' );
    else
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(within_diff21(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','black','FontSize',6,'fontweight', 'bold' );
    end
    %ylabel('Learning level (intercept)'); 
    xlabel('Change in FCV');
    ttle=xticklabels1(i); 
    title([ttle{1}]);
    
end
title(t, 'Subject level connectivity changes related to motor learning rate (Run 2 - Run 1), within networks');
%title(t, 'Subject level connectivity changes related to motor learning level (Run 2 - Run 1), within networks');
set(gcf,'color','w')

%% Within Run 3 Run 1
close all;
t=tiledlayout(2,6)
for i=1:12
    nexttile
    scatter(within_diff31(i,:), cfs(:,2), 4, 'filled')
    hold on
    lsline
    ylabel('Learning rate'); 
    [r, p] = corr(within_diff31(i,:)', cfs(:,2));
    
    
    if p<0.05; 
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(within_diff31(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','green','FontSize',6,'fontweight', 'bold' );
    else
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(within_diff31(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','black','FontSize',6,'fontweight', 'bold' );
    end
    %ylabel('Learning level (intercept)'); 
    xlabel('Change in FCV');
    ttle=xticklabels1(i); 
    title([ttle{1}]);
    
end
title(t, 'Subject level connectivity changes related to motor learning rate (Run 3 - Run 1), within networks');
%title(t, 'Subject level connectivity changes related to motor learning level (Run 3 - Run 1), within networks');
set(gcf,'color','w')
%% BETWEEN Delta FCV
% Run 2 Run 1
xticklabels1={'Aud.', 'Cer.', 'DMN.', 'Dor.','Fro.','Mem.', 'Sal.', 'Som.', 'Sub.', 'Unc.', 'Vent.', 'Vis.'};
close all;
t=tiledlayout(2,6)
for i=1:12
    nexttile
    scatter(between_diff21(i,:), cfs(:,2), 4, 'filled')
    hold on
    lsline
    ylabel('Learning rate'); 
    [r, p] = corr(between_diff21(i,:)', cfs(:,2));
    
    if p<0.05; 
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(between_diff21(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','green','FontSize',6,'fontweight', 'bold' );
    else
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(between_diff21(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','black','FontSize',6,'fontweight', 'bold' );
    end
    %ylabel('Learning level (intercept)'); 
    xlabel('Change in FCV');
    ttle=xticklabels1(i); 
    title([ttle{1}]);
    
end
title(t, 'Subject level connectivity changes related to motor learning rate (Run 2 - Run 1), between networks');
%title(t, 'Subject level connectivity changes related to motor learning level (Run 2 - Run 1), between networks');
set(gcf,'color','w')

%% Run 3 Run 1
close all;
t=tiledlayout(2,6)
for i=1:12
    nexttile
    scatter(between_diff31(i,:), cfs(:,2), 4, 'filled')
    hold on
    lsline
    ylabel('Learning rate'); 
    [r, p] = corr(between_diff31(i,:)', cfs(:,2));
    
    
    if p<0.05; 
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(between_diff31(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','green','FontSize',6,'fontweight', 'bold' );
    else
        r=num2str(round(r,2)); p=num2str(round(p,3));
        text(median(between_diff31(i,:)),0.9*max(cfs(:,2)),['r: ', r,',  p: ',p], 'Color','black','FontSize',6,'fontweight', 'bold' );
    end
    %ylabel('Learning level (intercept)'); 
    xlabel('Change in FCV');
    ttle=xticklabels1(i); 
    title([ttle{1}]);
    
end
title(t, 'Subject level connectivity changes related to motor learning rate (Run 3 - Run 1), between networks');
%title(t, 'Subject level connectivity changes related to motor learning level (Run 3 - Run 1), within networks');
set(gcf,'color','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. 
%% 4.b Run the clustering analysis with 
n_clusters = 5;
[idx, C, sumd] = kmeans(aggregate_dFCw_2D',n_clusters, 'replicates', 1); % k means take input 91x34716 <- number of window layers x no_subjects X number of entries in the matrix

%% 4.a k-means clustering of the dFCs, aggregate level(as in Allen et al., 2014)
% 
% close all;
% t=tiledlayout(2,5)
% for n_clusters=2:8%15
%     C=Cs(n_clusters).C;
%     sim_states=zeros(n_clusters,n_clusters);
%     for i=1:n_clusters
%         for k=1:n_clusters
%             if i~=k
%             sim_states(i,k)=corr(C(i,:)',C(k,:)');
%             end
%         end
%     end
% sim_states=tril(sim_states);
% nexttile
% imagesc(sim_states,[0, 0.9]); 
% myColorMap = jet(256); myColorMap(1,:) = 1; colormap(myColorMap); % make zero entries white
% title(sprintf('%d clusters',n_clusters))
% end
% c=colorbar; ylabel(c,'Correlation based similarity')
% title(t, 'Pairwise similarity between two brain states for each k-means cluster parameter')
% set(gcf,'color','w'); % make plot background white
%%
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


% Import the data 
% 5_6  3rd state 11 subj,  1st state corr 0.8 only 6 subjects
% 5_14 no run diff
% 5_15 no
% 5_23 concerns 10 subjects, corr 0.36; largest mean l.r.
% 5_36 no
% 5_38 YES

idx = readtable("C:\Users\user\Documents\Manon\Sliding Window Analysis\HPF_0.01852_60TR\idx13_5_43.csv", opts);
idx=idx.VarName1;
n_clusters = 5;
% Visualise the brain states obtained from k means algorithm
close all;
t=tiledlayout(2,ceil(n_clusters/2));
for i=1:n_clusters
    FC_aux = zeros(N,N);
    FC_aux(mask_ut) = mean(aggregate_dFCw_2D(:,idx==i),2);
    FC_aux = FC_aux + FC_aux';
    nexttile
    imagesc(FC_aux); axis square;  colormap(jet);
    title(sprintf('Brain state %d ' ,i)); colormap jet;
    xlabel('brain regions'); ylabel('brain regions'); caxis([-.7 .7]); 
    set(gcf,'color','w');
    if i==3; c = colorbar; c.Label.String = 'Pearson''s correlation'; end

end
title(t,['dFC sliding window states with window length ',num2str(wSize), ' TRs - comparing ',num2str(text_runs)])
% number of windows per state
window_idx=[idx [repelem(1,length(idx)/2) repelem(2,length(idx)/2)]'];

if length(histcounts(idx(1:length(idx)/2))) == length(histcounts(idx((1+length(idx)/2):end)))
    counts=[histcounts(idx(1:length(idx)/2)); 
           histcounts(idx((1+length(idx)/2):end))]; end
   
if length(histcounts(idx(1:length(idx)/2))) > length(histcounts(idx((1+length(idx)/2):end)))
    counts=[histcounts(idx(1:length(idx)/2));
        [histcounts(idx((1+length(idx)/2):end)) 0]]; % one of the runs have no windows in the last cluster, add 0 manually
end

if length(histcounts(idx(1:length(idx)/2))) < length(histcounts(idx((1+length(idx)/2):end)))
        counts=[[histcounts(idx(1:length(idx)/2)) 0 0];
        histcounts(idx((1+length(idx)/2):end))];
end
    
clc; nexttile,
h=bar(counts', 'stacked'); set(h, {'DisplayName'}, {'Run 1','Run 2'}'); legend('Location','northeastoutside');
   title('Dwelling time per brain state and across runs', 'fontsize', 8); xlabel('brain state'); ylabel('number of TRs spent overall');
   set(gcf,'color','w');

%% Subject level brain state, stacked bar
subjects_states_run1=zeros(size(TS_1,3),n_clusters);
subjects_states_run2=zeros(size(TS_1,3),n_clusters);
change_state_run1=zeros(size(TS_1,3),1);
change_state_run2=zeros(size(TS_1,3),1);

for subj_no=1:size(TS_1,3)
    start=1+(subj_no-1)*numLayers;
    fin=subj_no*numLayers;
    states=unique(idx(start:fin));
    change_state_run1(subj_no)=length(states);
    
    num=zeros(1,length(states));
    for i=1:length(num)
        num(i)=sum(states(i)==idx(start:fin))
    end
    
    for s=1:length(states)
    subjects_states_run1(subj_no,states(s))=num(s);
    end
    
    % Run 2
    start=length(idx)/2+1+(subj_no-1)*numLayers;
    fin=length(idx)/2+subj_no*numLayers;
    states=unique(idx(start:fin));
    change_state_run2(subj_no)=length(states);
    
    num=zeros(1,length(states));
    for i=1:length(num)
        num(i)=sum(states(i)==idx(start:fin))
    end
    
    for s=1:length(states)
    subjects_states_run2(subj_no,states(s))=num(s);
    end
end
%
close all; clc
xtick1=1:size(TS_1,3); 
xticklabels1=[1:14 16:32 34:37]; % 35 subjects

tiledlayout(2,1); nexttile
b=bar(subjects_states_run1, 'stacked'); 
set(b, {'DisplayName'}, {'State 1','State 2','State 3','State 4','State 5'}'); 
legend('location', 'northeastoutside'); 
set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontSize', 6);
title('Run1')
ylabel('Time spent in state (TR)')

nexttile
b=bar(subjects_states_run2, 'stacked'); 
set(b, {'DisplayName'}, {'State 1','State 2','State 3','State 4','State 5'}'); 
legend('location', 'northeastoutside');
title('Run3'); set(gca,'xtick',xtick1, 'xticklabel', xticklabels1, 'fontSize', 6);
set(gcf,'color','w')
ylabel('Time spent in state (TR)')
%% Dwelling time ~ Learning rate
clc;
subjects_states_run12=subjects_states_run1+subjects_states_run2;
dwell_learn=zeros(6,n_clusters);
for state=1:n_clusters
    %Run1
    indc=find(subjects_states_run1(:,state)); % only those subjects who are present in this state
    if ~isempty(indc)
    dwell_learn(1,state)=corr(subjects_states_run1(indc,state),cfs(indc,2));
    end

    % Or leave all subjects with 0 values
    dwell_learn(2,state)=corr(subjects_states_run1(:,state),cfs(:,2));
    
    %Run2
    indc=find(subjects_states_run2(:,state)); % only those subjects who are present in this state
    if ~isempty(indc)
      dwell_learn(3,state)=corr(subjects_states_run2(indc,state),cfs(indc,2));
    end

    % Or leave all subjects with 0 values
    dwell_learn(4,state)=corr(subjects_states_run2(:,state),cfs(:,2));
    
    %Run1 and 2 together
    indc=find(subjects_states_run12(:,state)); % only those subjects who are present in this state
    if ~isempty(indc)
    dwell_learn(5,state)=corr(subjects_states_run12(indc,state),cfs(indc,2));
    end

    % Or leave all subjects with 0 values
    dwell_learn(6,state)=corr(subjects_states_run12(:,state),cfs(:,2))
end
disp(sprintf('Max dwell~learn correlation: %d', max(abs(dwell_learn(:)))))
%%
clc; state_mean_learn=zeros(1,n_clusters);
for state=1:n_clusters
state_mean_learn(state)=sum(subjects_states_run12(:,state).*cfs(:,2))/sum(subjects_states_run12(:,state));
end
disp(state_mean_learn)
%% Set state and RUN
clc;
close all;
state=1
%tiledlayout(2,3)
%for state=1:n_clusters
%    nexttile
indc=find(subjects_states_run12(:,state)); 
scatter(subjects_states_run12(indc,state),cfs(indc,2)); 
% scatter(subjects_states_run1(:,state),cfs(:,2));
xlabel('Dwelling time in this state'); ylabel('Motor learning rate');
title(sprintf('Brain state %d',state));
hold on; lsline
[r, p]=corr(subjects_states_run12(indc,state),cfs(indc,2));
r=num2str(round(r,2));
text(0.9*max(subjects_states_run12(indc,state)),0.1*max(cfs(indc,2)),['r: ',r])
%end
%% Change state ~ learning rate
clc; 
corr(cfs(:,2),change_state_run1(:))
corr(cfs(:,2),change_state_run2(:))
x=change_state_run1(:)+change_state_run2(:);
corr(cfs(:,2),x)
%%
close all
scatter(cfs(:,2),change_state_run2(:))
%% % in one state
clc;
corr(max(subjects_states_run1')'/numLayers,cfs(:,2))
corr(max(subjects_states_run2')'/numLayers,cfs(:,2))
x=0.5*(max(subjects_states_run1')'/numLayers+max(subjects_states_run2')'/numLayers);
corr(x,cfs(:,2))
%% Run 1_2: Those that don't have state 4
disappeared=[1 3 5 11 12 13 14 16 17 19 21 25 31 34]; % reformat because of excluded subjects!
close all; clc;

scatter(1:14,cfs(disappeared,2))
mean(cfs(disappeared,2))

rest=[2 4 6 7 8 9 10 15 18 20 22 23 24 26 27 28 29 30 32 33 35];
mean(cfs(rest,2))
%% Run1_3
disappeared=[1 3 5 13 14 16 19 21 25 32 34]; % reformat because of excluded subjects!
close all; clc;

scatter(1:11,cfs(disappeared,2))
mean(cfs(disappeared,2))

rest=[2 4 6 7 8 9 10 11 12 15 18 20 22 23 24 26 27 28 29 30 31 33 35];
mean(cfs(rest,2))
