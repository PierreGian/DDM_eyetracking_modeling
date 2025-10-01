% Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis
% 12/2023
clear all

%%%%%%%%%%%% Initialize key parameters for batch script here
analysis_type = 'VE'; %options are 'VECD' and 'VE'
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

All_choice = struct();
All_choice.Amp4_L = zeros(length(subject_list),length(X));
All_choice.Amp4_R = zeros(length(subject_list),length(X));
All_choice.Amp4_Total = zeros(length(subject_list),length(X));
All_choice.Amp8_L = zeros(length(subject_list),length(X));
All_choice.Amp8_R = zeros(length(subject_list),length(X));
All_choice.Amp8_Total = zeros(length(subject_list),length(X));

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
    
    %%%%%%%%%%%% !! Need to recompute forward data in VE manner !! %%%%%%%%%%
    
    RT_range = [200 1200]; %RT bounds for this experiment
    VisualError=cur_subj.saccades.visualerror; %1x312 vector
    VisualError(isnan(VisualError))=0;
    ToPrint = ['Dealing with ',analysis_type,' case'];
    disp(ToPrint)
    saccade_idx = cur_subj.saccades.valid_idx;
    Responses = cur_subj.TaskData.mat_response;
    Accuracy = cur_subj.TaskData.mat_accuracy;
    Direction = cur_subj.TaskData.mat_dirs;
    Norm_RTs = zeros(1,length(Responses));
    RT = cur_subj.TaskData.mat_RTs;
    RT_toofast = RT<RT_range(1);
    RT_tooslow = RT>RT_range(2);
    RT(RT_toofast)=0; %filter RTs based on RT range
    RT(RT_tooslow)=0;
    ValidResponse = (Responses>0); %critical for choice_data
    x_relative = X;

    ForwardVE = VisualError>0;
    for nn=[1:length(VisualError)]
        if(strcmp(Direction{nn},'L'))
            ForwardVE(nn)=1-ForwardVE(nn);
        end
    end
    
    %Computation of valid indices for RT mean/SE extraction
    %all have to have true values (not zero)
    blocknb = 10;
    trialnb = 48;
    idx_valid = saccade_idx(1:blocknb*trialnb) & RT(1:blocknb*trialnb) & Accuracy(1:blocknb*trialnb) & VisualError(1:blocknb*trialnb); %only use data with detected saccades, accurate response times, & correct trials
    
    %Computation of valid indices for choice_data
    idx_include = saccade_idx(1:blocknb*trialnb) & ValidResponse(1:blocknb*trialnb) & VisualError(1:blocknb*trialnb); %choice_data indices to include %Potentially add RT as criterion as well??  

    %index variables to get all true values in vector
    Direction=cur_subj.TaskData.mat_dirs(1:blocknb*trialnb);
    Amplitude=cur_subj.TaskData.mat_amps(1:blocknb*trialnb);
    Jumpdir=cur_subj.TaskData.mat_jumpdir(1:blocknb*trialnb);
    Jumpsize=cur_subj.TaskData.mat_jumpsize(1:blocknb*trialnb);
    Resp=cur_subj.TaskData.mat_response(1:blocknb*trialnb);
    edat_forward=cur_subj.TaskData.mat_forwresp(1:blocknb*trialnb); %1 means forward, 0 means backward, 2 means excluded trials

    % index forward jumps and assign backwards jumps with (-) value
    forward_idx=strcmp(Jumpdir,'F'); %identify forward saccade jumps
    Jumpsize_relative = Jumpsize;
    Jumpsize_relative(~forward_idx)=Jumpsize_relative(~forward_idx)*-1;

    Amp4_idx=Amplitude==4; %how far the first target is
    Amp8_idx=Amplitude==8;
     
     edat_forward_idx = (edat_forward==1);
     fresp_zero_idx(1,:) = (Jumpsize(1,:) == 0) & edat_forward==1; %0 case, valid id, and edat_forward is 1
     bresp_zero_idx(1,:) = (Jumpsize(1,:) == 0) & edat_forward==0; %0 case, valid id, and edat_forward is 0
     %full_nonzero_idx = ~((Jumpsize(1,:) == 0) & idx_valid & (edat_forward==1 | edat_forward==0));

     forwmat_resp_zero_idx = ((Jumpsize(1,:) == 0) & Resp==1 & strcmp(Direction,'L')) | ((Jumpsize(1,:) == 0) & Resp==4 & strcmp(Direction,'R'));
     backmat_resp_zero_idx = ((Jumpsize(1,:) == 0) & Resp==1 & strcmp(Direction,'R')) | ((Jumpsize(1,:) == 0) & Resp==4 & strcmp(Direction,'L'));

     backward_idx = strcmp(Jumpdir,'B');
     bresp_rem_idx = backward_idx.*(bresp_zero_idx|fresp_zero_idx); 
     bresp_idx = backward_idx-bresp_rem_idx;%indices corresponding to non-zero backward responses
     bresp_idx = bresp_idx + backmat_resp_zero_idx; %add backward responses corresponding to 0 case based on actual responses
     bresp_idx = bresp_idx & idx_valid; %only keep valid indices
     fresp_rem_idx = forward_idx.*(bresp_zero_idx|fresp_zero_idx);
     fresp_idx = forward_idx-fresp_rem_idx;%indices corresponding to non-zero forward responses
     fresp_idx = fresp_idx + forwmat_resp_zero_idx; %add forward responses corresponding to 0 case based on actual responses
     %fresp_idx = forward_idx + bresp_rem_idx; %add 0 indices to forward indices vector
     fresp_idx = fresp_idx & idx_valid; %only keep valid indices

     RT_idx.forward = fresp_idx;
     RT_idx.backward = bresp_idx;
     RT_idx.toofast = RT_toofast;
     RT_idx.tooslow = RT_tooslow;

     left_idx = strcmp(Direction,'L');
     right_idx = strcmp(Direction,'R');

     choice_data = struct();
     Xsize = length(x_relative);
     Amp_fields = {'Amp4','Amp8'};
     Amp_vals = [4,8];
     Dir_fields = {'Left','Right','Total'};
     choice_data.Amp4 = struct();
     choice_data.Amp8 = struct();
     for f=1:length(Amp_fields)
        choice_data.(Amp_fields{f}).Left = struct();
        choice_data.(Amp_fields{f}).Right = struct();
        choice_data.(Amp_fields{f}).Total = struct();
        for ff=1:length(Dir_fields)
            choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwSumTrials = zeros(1,Xsize);
            choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwResp = zeros(1,Xsize);
            choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwPercent = zeros(1,Xsize);
        end
     end
     
     %Need to define/return two separate functions

     %target_dir is 1 for lefts, 4 for rights
     for f=1:length(Amp_fields)
         if(strcmp(Amp_fields{f},'Amp4'))
             c_amp=Amp4_idx;
             %disp('Amp4')
         else
             c_amp=Amp8_idx;
             %disp('Amp8')
         end
         for ff=1:2
            if(strcmp(Dir_fields{ff},'Left'))
                c_dir=left_idx;
                %disp('Left indices')
            else
                c_dir=right_idx;
                %disp('Right indices')
            end
            for xx=[1:length(x_relative)]
                 choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwSumTrials(1,xx) = sum(c_dir&c_amp&idx_include&(Jumpsize_relative==x_relative(xx)));
                 choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwResp(1,xx) = sum(c_dir&c_amp&idx_include&(Jumpsize_relative==x_relative(xx))&ForwardVE);
                 choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwPercent(1,xx) = choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwResp(xx)/choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwSumTrials(xx);
            end
        end
        last_dir = Dir_fields{3};
        choice_data.(Amp_fields{f}).(last_dir).ForwSumTrials = choice_data.(Amp_fields{f}).(Dir_fields{1}).ForwSumTrials+choice_data.(Amp_fields{f}).(Dir_fields{2}).ForwSumTrials;
        choice_data.(Amp_fields{f}).(last_dir).ForwResp = choice_data.(Amp_fields{f}).(Dir_fields{1}).ForwResp+choice_data.(Amp_fields{f}).(Dir_fields{2}).ForwResp;
        choice_data.(Amp_fields{f}).(last_dir).ForwPercent = choice_data.(Amp_fields{f}).(last_dir).ForwResp./choice_data.(Amp_fields{f}).(last_dir).ForwSumTrials;
     end

    %Saving data in big matrix
    All_choice.Amp4_L(n,:) = choice_data.Amp4.Left.ForwPercent;
    All_choice.Amp4_R(n,:) = choice_data.Amp4.Right.ForwPercent;
    All_choice.Amp4_Total(n,:) = choice_data.Amp4.Total.ForwPercent;
    All_choice.Amp8_L(n,:) = choice_data.Amp8.Left.ForwPercent;
    All_choice.Amp8_R(n,:) = choice_data.Amp8.Right.ForwPercent;
    All_choice.Amp8_Total(n,:) = choice_data.Amp8.Total.ForwPercent;

    %single subject fitting
    [B4L(1,:),Y4L(1,:)] = S_curve_calculating_fitting(choice_data.Amp4.Left.ForwPercent,X);
    [B4R(1,:),Y4R(1,:)] = S_curve_calculating_fitting(choice_data.Amp4.Right.ForwPercent,X);
    [B4(1,:),Y4(1,:)] = S_curve_calculating_fitting(choice_data.Amp4.Total.ForwPercent,X);
    [B8L(1,:),Y8L(1,:)] = S_curve_calculating_fitting(choice_data.Amp8.Left.ForwPercent,X);
    [B8R(1,:),Y8R(1,:)] = S_curve_calculating_fitting(choice_data.Amp8.Right.ForwPercent,X);
    [B8(1,:),Y8(1,:)] = S_curve_calculating_fitting(choice_data.Amp8.Total.ForwPercent,X);

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
    kk=5;
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
    %clearvars results -except BICpsych
    %clear -regexp B* Y* -except BICpsych

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

%% Compute averages here
Amp4_L_avg = nanmean(All_choice.Amp4_L,1);
Amp4_R_avg = nanmean(All_choice.Amp4_R,1);
Amp4_Total = nanmean(All_choice.Amp4_Total,1);
Amp8_L_avg = nanmean(All_choice.Amp8_L,1);
Amp8_R_avg = nanmean(All_choice.Amp8_R,1);
Amp8_Total = nanmean(All_choice.Amp8_Total,1);

%Get B4, Y4, B8, and Y8
[B4L(1,:),Y4L(1,:)] = S_curve_calculating_fitting(Amp4_L_avg,X);
[B4R(1,:),Y4R(1,:)] = S_curve_calculating_fitting(Amp4_R_avg,X);
[B4(1,:),Y4(1,:)] = S_curve_calculating_fitting(Amp4_Total,X);
[B8L(1,:),Y8L(1,:)] = S_curve_calculating_fitting(Amp8_L_avg,X);
[B8R(1,:),Y8R(1,:)] = S_curve_calculating_fitting(Amp8_R_avg,X);
[B8(1,:),Y8(1,:)] = S_curve_calculating_fitting(Amp8_Total,X);

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

biasandthresholds_savecsv_behave_Weiwei(X,B4(1,:),Y4(1,:),B8(1,:),Y8(1,:),'Group average',{analysis_type},path_output);

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

sgtitle(['Group Average Behave ' analysis_type])

set(h, 'Position', get(0, 'Screensize'));
saveas(h, [path_plots analysis_type '/Group_average_' analysis_type '.jpg']);
close(h)

%Psychometric results
fit_data2 = cell(6,length(x_space)+1);
fit_data2{1,1} = 'x';
fit_data2{2,1} = '4L';
fit_data2{3,1} = '4R';
fit_data2{4,1} = '4';
fit_data2{5,1} = '8L';
fit_data2{6,1} = '8R';
fit_data2{7,1} = '8';
fit_data2{1,2} = x_space(2);
fit_data2{2,2} = fitValues4L(2)*100;
fit_data2{3,2} = fitValues4R(2)*100;
fit_data2{4,2} = fitValues4(2)*100;
fit_data2{5,2} = fitValues8L(2)*100;
fit_data2{6,2} = fitValues8R(2)*100;
fit_data2{7,2} = fitValues8(2)*100;

output_fpath2 = [path_psychometric 'Average_fitValues.csv'];
writecell(fit_data2, output_fpath2);

% Full psychometric table
psychometric_table2 = cell(1,Xsize);
psychometric_table2{1,1} = 'ts';
x_idx = 2;
for xx=X
    psychometric_table2{1,x_idx} = num2str(xx);
    x_idx=x_idx+1;
end
psychometric_table2{2,1} = '4L';
curtab = result4L.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{2,p_idx} = pp;
    p_idx = p_idx+1;
end
psychometric_table2{3,1} = '4R';
curtab = result4R.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{3,p_idx} = pp;
    p_idx = p_idx+1;
end
psychometric_table2{4,1} = '4';
curtab = result4.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{4,p_idx} = pp;
    p_idx = p_idx+1;
end
psychometric_table2{5,1} = '8L';
curtab = result8L.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{5,p_idx} = pp;
    p_idx = p_idx+1;
end
psychometric_table2{6,1} = '8R';
curtab = result8R.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{6,p_idx} = pp;
    p_idx = p_idx+1;
end
psychometric_table2{7,1} = '8';
curtab = result8.data;
probs = ((curtab(:,2)./curtab(:,3)).*100)';
p_idx = 2;
for pp=probs
    psychometric_table2{7,p_idx} = pp;
    p_idx = p_idx+1;
end

output_fpath2 = [path_psychometric 'Average_psychometric_ts.csv'];
writecell(psychometric_table2, output_fpath2);

nameMat0 = [path_psychometric,'AllChoice_',analysis_type,'.mat'];
save(nameMat0,'All_choice')

%End of plotting group results
%% Final part of VE Forward curve analyses
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

% bic_struc = fullfile(path_psychometric , ['BICpsych.mat']);
% save(bic_struc,'BICpsych')

function [BICret] = computeBIC(y,yhat,k,N)
    indices = find(isnan(y) == 1);
    y(indices) = [];
    yhat(indices) = [];
    RSS = sum((y-yhat).^2);
    N=N-length(indices);
    BICret = N*log((RSS/N))+k*log(N);
end