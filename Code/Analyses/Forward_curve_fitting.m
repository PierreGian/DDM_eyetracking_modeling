% Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis
% 12/2023
clear all

%%%%%%%%%%%% Initialize key parameters for batch script here
analysis_type = 'VECD'; %options are 'VECD' and 'VE'
text_extension = '/*.csv'; % asign txt or csv
dist2screen=1000; % 1000 for behave, 1450 for scan
%%%%%%%%%%%%%%%%%%%%%

%Path initialization
root_path = '/nfs/agency/';
master_path = [root_path,'code/DDM/DDM_Master_code/'];
addpath(genpath([root_path,'code/DDM/DDM_Master_code/']));
behave_stats_path = [root_path,'raw_data/DDM/']; %eye-tracking data from behavioral task
addpath('/nfs/agency/code/ForwardCurve/functions/')
addpath(genpath('/nfs/agency/code/ForwardCurve/psignifit/'))
% addpath('/nfs/agency/code/ForwardCurve/Weiwei/')
% addpath('/nfs/agency/code/DDM/')
path_crunch = [master_path,'eyetracking_crunch/behave/'];
path_edf = '/nfs/agency/raw_data/DDM/edf/';
path_csv = '/nfs/agency/raw_data/DDM/edats/';
path_output = [master_path 'outputs/'];
path_psychometric = [path_output 'Psychometric/'];
path_plots = [master_path,'plots/ForwardCurve/'];
behave_data_folder = {dir(path_edf).name}; %Folder with all subject names

%Outputting the paths corresponding to edf, csv, and crunch files
disp(['path_edf is ' path_edf])
disp(['path_csv is ' path_csv])
disp(['path_crunch is ' path_crunch])

%Vector of jump sizes - input required for modeling analyses
%X=[-2.5 -1.75 -1.25 -1 -0.75 -0.5 -0.25 -0.175 -0.1 0 0.1 0.175 0.25 0.5 0.75 1 1.25 1.75 2.5];
X=[-2.5 -1.75 -1 -0.5 -0.25 0 0.25 0.5 1 1.75 2.5];
Xsize = length(X);
Amp_fields = {'Amp4','Amp8'};
Dir_fields = {'Left','Right','Total'};
% find the index for 0 jump
zero_jump_index = find(X==0);

%Identifying list of subject names

Data_ToLoad = [path_output,'Subject_Data.mat'];
load(Data_ToLoad)
subject_list = fieldnames(Subject_Data)';

%Extract subject names from data struct

%Need to plot full case
N = length(subject_list);

%Select_idx = [1:length(subject_list)];
%Select_idx = [1,6];

biasandthresholds = cell(1,13);
biasandthresholds{1,1} = 'subject ID';
biasandthresholds{1,2} = 'bias4 right';
biasandthresholds{1,3} = 'threshold4 right';
biasandthresholds{1,4} = 'bias4 left';
biasandthresholds{1,5} = 'threshold4 left';
biasandthresholds{1,6} = 'bias4';
biasandthresholds{1,7} = 'threshold4';
biasandthresholds{1,8} = 'bias8 right';
biasandthresholds{1,9} = 'threshold8 right';
biasandthresholds{1,10} = 'bias8 left';
biasandthresholds{1,11} = 'threshold8 left';
biasandthresholds{1,12} = 'bias8';
biasandthresholds{1,13} = 'threshold8';


FCdata = struct();
FCdata.IDs = subject_list;
%FCdata.conds = conds_dat; %Add conditions later

BICpsych = struct();
BICpsych.Bayes = struct();
BICpsych.Bayes.Amp4 = zeros(1,length(X));
BICpsych.Bayes.Amp8 = zeros(1,length(X));

RTBICpsych = struct();
RTBICpsych.Amp4 = zeros(1,length(X));
RTBICpsych.Amp8 = zeros(1,length(X));

%for n=1:length(subject_list)
%for n = Select_idx
for n=1:length(subject_list)
%for n=[1]
    cur_name = subject_list{n};
    disp(['running ' cur_name])
    cur_subj = Subject_Data.(cur_name);
    totaltrials = cur_subj.counts.trials_count;
    
    if isfolder([path_plots analysis_type])
        disp(['saving jpgs in ' path_plots analysis_type])
    else
        mkdir([path_plots analysis_type])
        disp(['saving jpgs in ' path_plots analysis_type])
    end
    
    %single subject fitting
    [B4L(1,:),Y4L(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp4.Left.ForwPercent,X);
    [B4R(1,:),Y4R(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp4.Right.ForwPercent,X);
    [B4(1,:),Y4(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp4.Total.ForwPercent,X);
    [B8L(1,:),Y8L(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp8.Left.ForwPercent,X);
    [B8R(1,:),Y8R(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp8.Right.ForwPercent,X);
    [B8(1,:),Y8(1,:)] = S_curve_calculating_fitting(cur_subj.choice_data.Amp8.Total.ForwPercent,X);

    %[cur_bias,cur_thresh] = S_plot_w_bias_thresh(X,B8L,Y8L);

    %for amp 4
    gumbeldata4L(:,1) = X;
    gumbeldata4L(:,2) = Y4L*totaltrials;
    gumbeldata4L(:,3) = totaltrials;

    gumbeldata4R(:,1) = X;
    gumbeldata4R(:,2) = Y4R*totaltrials;
    gumbeldata4R(:,3) = totaltrials;

    gumbeldata4(:,1) = X;
    gumbeldata4(:,2) = Y4*totaltrials;
    gumbeldata4(:,3) = totaltrials;

    %amp 8
    gumbeldata8L(:,1) = X;
    gumbeldata8L(:,2) = Y8L*totaltrials;
    gumbeldata8L(:,3) = totaltrials;

    gumbeldata8R(:,1) = X;
    gumbeldata8R(:,2) = Y8R*totaltrials;
    gumbeldata8R(:,3) = totaltrials;

    gumbeldata8(:,1) = X;
    gumbeldata8(:,2) = Y8*totaltrials;
    gumbeldata8(:,3) = totaltrials;

    %call psignifit script - amp4
    [result4L] = psignifit_NCAP_Weiwei(gumbeldata4L);
    [result4R] = psignifit_NCAP_Weiwei(gumbeldata4R);
    [result4] = psignifit_NCAP_Weiwei(gumbeldata4);

    %call psignifit script - amp8
    [result8L] = psignifit_NCAP_Weiwei(gumbeldata8L);
    [result8R] = psignifit_NCAP_Weiwei(gumbeldata8R);
    [result8] = psignifit_NCAP_Weiwei(gumbeldata8);

    yfit4 = result4.psiHandle;
    yfit8 = result8.psiHandle;
    
    y_vec4 = Y4;
    y_est_vec4 = yfit4(X);
    y_vec8 = Y8;
    y_est_vec8 = yfit8(X);
    kk=4;
    BICpsych.Bayes.Amp4(n) = computeBIC(y_vec4,y_est_vec4,kk,length(X));
    BICpsych.Bayes.Amp8(n) = computeBIC(y_vec8,y_est_vec8,kk,length(X));
    

    %amp 4
%     biasandthresholds_savecsv_scan_left(X,result4L,subject_list{n},analysis_type,path2output);
%     biasandthresholds_savecsv_scan_right(X,result4R,subject_list{n},analysis_type,path2output);
%     biasandthresholds_savecsv_scan(X,result4,subject_list{n},analysis_type,path2output);
% 
%     %amp 8
%     biasandthresholds_savecsv_scan_left(X,result8L,subject_list{n},analysis_type,path2output);
%     biasandthresholds_savecsv_scan_right(X,result8R,subject_list{n},analysis_type,path2output);
%     biasandthresholds_savecsv_scan(X,result8,subject_list{n},analysis_type,path2output);

    biasandthresholds_savecsv_behave_Weiwei(X,B4(1,:),Y4(1,:),B8(1,:),Y8(1,:),subject_list{n},{analysis_type},path_output);
    
    % curve ploting
    h = figure;
    subplot(3,2,1)    %Y4L
    fitValues4L = plotPsych_NCAP_Weiwei(result4L);
    legend('4 degree Leftward','Location','northwest');

    subplot(3,2,2)      %Y8L
    fitValues8L = plotPsych_NCAP_Weiwei(result8L);
    legend('8 degree Leftward','Location','northwest');

    subplot(3,2,3)    %Y4R
    fitValues4R = plotPsych_NCAP_Weiwei(result4R);
    legend('4 degree Rightward','Location','northwest');

    subplot(3,2,4)      %Y8R
    fitValues8R = plotPsych_NCAP_Weiwei(result8R);
    legend('8 degree Rightward','Location','northwest');

    subplot(3,2,5)    %Y4
    fitValues4 = plotPsych_NCAP_Weiwei(result4);
    legend('4 degree Total','Location','northwest');

    subplot(3,2,6)      %Y8
    fitValues8 = plotPsych_NCAP_Weiwei(result8);
    legend('8 degree Total','Location','northwest');

    x_space = linspace(min(result4.data(:,1)),max(result4.data(:,1)),1000);

    sgtitle([subject_list{n} ' Behave ' analysis_type])

    set(h, 'Position', get(0, 'Screensize'));
    saveas(h, [path_plots analysis_type '/' subject_list{n} '_' analysis_type '.jpg']);
    close(h)
    clearvars results -except BICpsych
    clear -regexp B* Y* -except BICpsych

    psychometric_table = cell(1,Xsize);
    psychometric_table{1,1} = 'ts';
    x_idx = 2;
    for xx=X
        psychometric_table{1,x_idx} = num2str(xx);
        x_idx=x_idx+1;
    end
    psychometric_table{2,1} = '4L';
    curtab = result4L.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{2,p_idx} = pp;
        p_idx = p_idx+1;
    end
    psychometric_table{3,1} = '4R';
    curtab = result4R.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{3,p_idx} = pp;
        p_idx = p_idx+1;
    end
    psychometric_table{4,1} = '4';
    curtab = result4.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{4,p_idx} = pp;
        p_idx = p_idx+1;
    end
    psychometric_table{5,1} = '8L';
    curtab = result8L.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{5,p_idx} = pp;
        p_idx = p_idx+1;
    end
    psychometric_table{6,1} = '8R';
    curtab = result8R.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{6,p_idx} = pp;
        p_idx = p_idx+1;
    end
    psychometric_table{7,1} = '8';
    curtab = result8.data;
    probs = ((curtab(:,2)./curtab(:,3)).*100)';
    p_idx = 2;
    for pp=probs
        psychometric_table{7,p_idx} = pp;
        p_idx = p_idx+1;
    end

    output_fpath = [path_psychometric cur_name '_psychometric_ts.csv'];
    writecell(psychometric_table, output_fpath);

    fit_data = cell(6,length(x_space)+1);
    fit_data{1,1} = 'x';
    fit_data{2,1} = '4L';
    fit_data{3,1} = '4R';
    fit_data{4,1} = '4';
    fit_data{5,1} = '8L';
    fit_data{6,1} = '8R';
    fit_data{7,1} = '8';
    for xx=[1:length(x_space)]
        fit_data{1,xx+1} = x_space(xx);
        fit_data{2,xx+1} = fitValues4L(xx)*100;
        fit_data{3,xx+1} = fitValues4R(xx)*100;
        fit_data{4,xx+1} = fitValues4(xx)*100;
        fit_data{5,xx+1} = fitValues8L(xx)*100;
        fit_data{6,xx+1} = fitValues8R(xx)*100;
        fit_data{7,xx+1} = fitValues8(xx)*100;
    end

    output_fpath2 = [path_psychometric cur_name '_fitValues.csv'];
    writecell(fit_data, output_fpath2);

end

All_B4 = zeros(1,length(subject_list));
All_T4 = zeros(1,length(subject_list));
All_B8 = zeros(1,length(subject_list));
All_T8 = zeros(1,length(subject_list));

%Read cvs file and create *.mat object
csvFile = [path_output,'/biasandthresholds_',analysis_type,'.csv'];
Table_Params = readtable(csvFile);
TableS = Table_Params.subjectID;
for s=[1:length(TableS)]
    curID = TableS{s};
    row = Table_Params(s,:);
    if(~strcmp(curID,''))
        cIDX = find(strcmp(TableS,curID));
        All_B4(cIDX) = row.bias4;
        All_B8(cIDX) = row.bias8;
        All_T4(cIDX) = row.threshold4;
        All_T8(cIDX) = row.threshold8;
    end
end

FCdata.bias4 = All_B4;
FCdata.bias8 = All_B8;
FCdata.threshold4 = All_T4;
FCdata.threshold8 = All_T8;
nameMat = [path_output,'FCdata_behave_',analysis_type,'.mat'];
save(nameMat,'FCdata')

bic_struc = fullfile(path_psychometric , ['BICpsych.mat']);
save(bic_struc,'BICpsych')

function [BICret] = computeBIC(y,yhat,k,N)
    indices = find(isnan(y) == 1);
    y(indices) = [];
    yhat(indices) = [];
    RSS = sum((y-yhat).^2);
    N=N-length(indices);
    BICret = N*log((RSS/N))+k*log(N);
end