% Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis
% 12/2023
clear all

%%%%%%%%%%%% Initialize key parameters for batch script here
analysis_type = 'VECD'; %options are 'VECD' and 'VE'
text_extension = '/*.csv'; % asign txt or csv
dist2screen=1000; % 1000 for behave, 1450 for scan
%%%%%%%%%%%%%%%%%%%%%

disp(['Running DDM modeling analysis for ' analysis_type ' case']);

%Path initialization
root_path = '/nfs/agency/';
master_path = [root_path,'code/DDM/DDM_Master_code/'];
addpath(genpath([root_path,'code/DDM/DDM_Master_code/']));
behave_stats_path = [root_path,'raw_data/behave_suite/']; %eye-tracking data from behavioral task
% addpath('/nfs/agency/code/ForwardCurve/Weiwei/')
% addpath('/nfs/agency/code/DDM/')
path_crunch = [master_path,'eyetracking_crunch/behave/'];
path_edf = '/nfs/agency/raw_data/DDM/edf/';
path_csv = '/nfs/agency/raw_data/DDM/edats/';
path_output = [master_path 'outputs/'];
path_DDM_output = [master_path 'outputs/DDM/'];
path_plots = [master_path,'plots/DDM/'];
path_RT_plots = [master_path,'plots/DDM/RT/'];
behave_data_folder = {dir(path_edf).name}; %Folder with all subject names

%Outputting the paths corresponding to edf, csv, and crunch files
disp(['path_edf is ' path_edf])
disp(['path_csv is ' path_csv])
disp(['path_crunch is ' path_crunch])

%Vector of jump sizes - input required for modeling analyses
%X=[-2.5 -1.75 -1.25 -1 -0.75 -0.5 -0.25 -0.175 -0.1 0 0.1 0.175 0.25 0.5 0.75 1 1.25 1.75 2.5];
X=[-2.5 -1.75 -1 -0.5 -0.25 0 0.25 0.5 1 1.75 2.5];
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

Select_idx = [1:length(subject_list)];

name_mat = [path_output,'FCdata_behave_',analysis_type,'.mat'];
load(name_mat)
FC_ID = FCdata.IDs;
All_B4 = FCdata.bias4;
All_T4 = FCdata.threshold4;
All_B8 = FCdata.bias8;
All_T8 = FCdata.threshold8;

All_bias = zeros(2,N);
All_threshold = zeros(2,N);
All_veW = nan(2,N);
All_DR = zeros(2,N);
All_DDMbias = zeros(2,N);
All_DDMbound = zeros(2,N);
All_DDMndt = zeros(2,N);
All_A_bound = zeros(2,N);
All_B_bound = zeros(2,N);
All_RT_s1 = zeros(2,N);
All_RT_smin1 = zeros(2,N);
All_RT_ratio_s1 = zeros(2,N);
All_RT_ratio_smin1 = zeros(2,N);

Txt_file = [path_output,'AllDDMparams_',analysis_type,'.csv'];
fileID = fopen(Txt_file,'w');

Str_format = '%s';
Nb_Params = 16;
Total_Cols = Nb_Params*2;
for inn = [1:Total_Cols]
    Str_format = [Str_format,'\t%s'];
end
Str_format = [Str_format,'\t%s\r\n'];

pref = {'Bias','Threshold','DDM_dr','DDM_bound','DDM_bias','DDM_ndt','DDM_ve_w','A_bound','B_bound','RT_s1','RT_smin1','RT_ratio_s1','RT_ratio_smin1','errFit','BIC','AIC'};
suffix = {'amp4_All','amp8_All'};
ToAppendStr = {'SubjectID','Accuracy'};
for ii = [1:Nb_Params]
    for jj = [1:2]
        ToAppendStr{1,end+1} = [pref{1,ii},'_',suffix{1,jj}];
    end
end
fprintf(fileID,Str_format,ToAppendStr{:});

%DDM inputs for DDM script

%fclose(fileID);

BICorig = struct();
BICorig.DDM = struct();
BICorig.DDM.Amp4 = zeros(1,length(X));
BICorig.DDM.Amp8 = zeros(1,length(X));

RTBICorig = struct();
RTBICorig.Amp4 = zeros(1,length(X));
RTBICorig.Amp8 = zeros(1,length(X));

for n=1:length(subject_list)
%for n = Select_idx
%for n=1:5
    cur_name = subject_list{n};
    disp(['running ' cur_name])
    mat_data = Subject_Data.(cur_name);
    rt_data = mat_data.RTs;
    choice_data = mat_data.choice_data;
    cur_accuracy = (sum(mat_data.TaskData.mat_accuracy)/length(mat_data.TaskData.mat_accuracy))*100;
    
    %% Format data for fit
    % data is a matrix with one row per stimulus stength
    % the 1st column is the stimulus strength
    % the 2nd column is the mean RT for "+"choices at that stimulus strength
    % the 3rd column is the SE of that mean RT
    % the 4th column is the mean RT for "-" choices at that stimulus strength
    % the 5th column is the SE of that mean RT
    % the 6th column is the number of "+" choices at that stimulus strength
    % the 7th column is the number of total trials at that stimulus strength
    
    data = struct();
    %Amplitude 4 data
    data.Amp4 = [X' rt_data.Amp4.Total.ForwardMean' rt_data.Amp4.Total.ForwardSEM' ...
        rt_data.Amp4.Total.BackwardMean' rt_data.Amp4.Total.BackwardSEM' ...
        choice_data.Amp4.Total.ForwResp' choice_data.Amp4.Total.ForwSumTrials'];
    %Amplitude 8 data
    data.Amp8 = [X' rt_data.Amp8.Total.ForwardMean' rt_data.Amp8.Total.ForwardSEM' ...
        rt_data.Amp8.Total.BackwardMean' rt_data.Amp8.Total.BackwardSEM' ...
        choice_data.Amp8.Total.ForwResp' choice_data.Amp8.Total.ForwSumTrials'];

    fields = fieldnames(data);
    labels = {'4 degrees','8 degrees'};
    for ff=1:length(fields)
        %turn 0s into Nans in 3rd and 5th column of the data matrix
        data_sz = size(data.(fields{ff}));
        for ii = [1:data_sz(1)]
           rowi = data.(fields{ff})(ii,:);
           if(rowi(3)==0)
               data.(fields{ff})(ii,3) = NaN;
           end
           if(rowi(5)==0)
               data.(fields{ff})(ii,5) = NaN;
           end
        end
    end
    
    for k = 1:length(fields)
        %% Run fit
        % The fitting code is flexible for difference scenarios.
        % These flags control which version of the DDM to use.
        fit_flags.conds = 1;
        fit_flags.type = 0;
        fit_flags.u0 = -1;
        % fit_flags.u0 = 1; % use for drift bias
        % fit_flags.barr = -1; uncomment this if you want to have no starting bias
        fit_flags.barr = 1; % uncomment for starting point if you want to have a starting bias
        fit_flags.resid = 1;
        fit_flags.nT1 = 1;
        fit_flags.subj_eq = 0;

        % the code tries to find the parameter values yielding maximum likelihood for the data
        % you need to start this function maximization somewhere, which you provide here
        % thetaGuess = [.1 20 400];   % no bias
        % thetaGuess = [.1 0.5 20 400];   % drift rate bias (2nd param)
        thetaGuess = [.1 20 1 400];   % starting point bias (3rd param)

        % you could change these values to control the function maximization algorithm
        opts = optimset('fminsearch');

        % for k = 1:length(fields)
        % this is where the work is done, returns best fit parameters
        [thetaFit,errFits(k),xflg,oput] = fminsearch('fitDiffMod', thetaGuess, opts, data.(fields{k}), fit_flags);
        thetafitvalues{k}=thetaFit;
        kkk = length(thetaFit);
        BICerr(k) = kkk*log(2*length(X))+2*errFits(k);
        AICerr(k) = 2*kkk+2*errFits(k);
        % end
        
        %% Plot results

        % typically, at this point, I would plot the fitted psychometric and chronometric functions from the best fit
        % parameters overlaid on the data. You can use the 'calcDiffMod' function to do this.
        % for k = 1:length(fields)
        % no bias
        % theta_cond = [thetaFit(1) 0 thetaFit(2) thetaFit(2) thetaFit(3) thetaFit(3)];
        % end
        % % drift rate bias
        % theta_cond = [thetaFit(1) thetaFit(2) thetaFit(3) thetaFit(3) thetaFit(4) thetaFit(4)];
        
        % % starting point bias
        theta_cond = [thetaFit(1) 0 thetaFit(2)-thetaFit(3) thetaFit(2)+thetaFit(3) thetaFit(4) thetaFit(4)];
        %NOTE: need to incorportate 
        
        % calculate the predictions
        % [t1pred, t2pred, ppred] = calcDiffMod(X, theta_cond, fit_flags.subj_eq);

        % calculate cmf/pmf fit
        x_fit = min(X):0.01:max(X);
        [t1pred{k}, t2pred{k}, ppred{k}] = calcDiffMod(x_fit, theta_cond, fit_flags.subj_eq);
        % calculate pmf
        pmf{k} = data.(fields{k})(:,6) ./ data.(fields{k})(:,7);
        % extract cmf
        cmf{k} = [data.(fields{k})(1:zero_jump_index-1,4); data.(fields{k})(zero_jump_index:end,2)];

    end
    %file_dir = [path2output_ddm subject_list{n} '_fits.txt'];
    %writecell(thetafitvalues, file_dir);
    
    for jj=[1:length(fields)]
        All_DR(jj,n) = thetafitvalues{jj}(1); %indexing the drift rate in vector
        All_DDMbias(jj,n) = thetafitvalues{jj}(3);
        All_DDMbound(jj,n) = thetafitvalues{jj}(2);
        All_DDMndt(jj,n) = thetafitvalues{jj}(4);
        All_A_bound(jj,n) = thetafitvalues{jj}(2)-thetafitvalues{jj}(3);
        All_B_bound(jj,n) = thetafitvalues{jj}(2)+thetafitvalues{jj}(3);
        u_s1 = All_DR(jj,n)*1;
        u_smin1 = All_DR(jj,n)*-1;
        All_RT_s1(jj,n) = (((All_A_bound(jj,n)+All_B_bound(jj,n))*coth((All_A_bound(jj,n)+All_B_bound(jj,n))*u_s1)) ./ u_s1) - (All_B_bound(jj,n)*coth((All_B_bound(jj,n)*u_s1)) ./ u_s1);
        All_RT_smin1(jj,n) = (((All_A_bound(jj,n)+All_B_bound(jj,n))*coth((All_A_bound(jj,n)+All_B_bound(jj,n))*u_smin1)) ./ u_smin1) - (All_A_bound(jj,n)*coth((All_A_bound(jj,n)*u_smin1)) ./ u_smin1);
        All_RT_ratio_s1(jj,n) = All_RT_s1(jj,n) / (All_RT_s1(jj,n) + All_DDMndt(jj,n));
        All_RT_ratio_smin1(jj,n) = All_RT_smin1(jj,n) / (All_RT_smin1(jj,n) + All_DDMndt(jj,n));
    end
    
    for k = 1:length(fields)
        temp = fields{k};
        tempname = strrep(temp,'Mean','');
        curvenames{k} = strrep(tempname,'_rt_','');
    end

    %get x fit index
    x_fit_idx = zeros(1,length(X));
    for xx = [1:length(X)]
        curX = X(xx);
        x_fit_idx(xx) = find(x_fit==curX);
    end
    %get y_vec & y_est_vec

    y_vec4 = pmf{1}';
    y_est_vec4 = ppred{1}(x_fit_idx);
    y_vec8 = pmf{2}';
    y_est_vec8 = ppred{2}(x_fit_idx);
    kk=4;
    BICorig.DDM.Amp4(n) = computeBIC(y_vec4,y_est_vec4,kk,length(X));
    BICorig.DDM.Amp8(n) = computeBIC(y_vec8,y_est_vec8,kk,length(X));


    h = figure;
    for m = 1:length(fields)
        % plot pmf and fit
        subplot(1,2,m)
        hold on;
        plot(X,pmf{m},'ro');
        plot(x_fit,ppred{m},'r');
        xlabel('Target Shift (degrees)');
        ylabel('Proportion forward responses');
        if max(ppred{m}) >= 0.75 && min(ppred{m}) <= 0.5
            threshold_DDM(m)= x_fit(find((ppred{m}>=0.75),1)) - x_fit(find((ppred{m}>=0.5),1));
            bias_DDM(m) = x_fit(find((ppred{m}>=0.5),1));
        else
            threshold_DDM(m) = x_fit(end) - x_fit(find((ppred{m}>=0.5),1));
            bias_DDM(m) = x_fit(find((ppred{m}>=0.5),1));
        end
        title(strrep(labels{m},'Mean_rt_',''),'interpreter','none');
        % legend(curvenames{m},'Location','northwest');
        % title([subject_list{n}  ' Psychometric Curve'])
    end
    sgtitle([cur_name  ' Psychometric Curve'])
    saveas(h, [path_plots cur_name '_DDM_forwardcurve.jpg']);
    close(h)

    p = figure;
    for m = 1:length(fields)
        % plot cmf and fit
        subplot(1,2,m)
        hold on;
        errorbar(X,cmf{m}',rt_data.(fields{m}).Total.FullSEM,'bo');
        errorbar(0,rt_data.(fields{m}).Total.ForwardMean(zero_jump_index),rt_data.(fields{m}).Total.ForwardSEM(zero_jump_index), 'ro');
        errorbar(0,rt_data.(fields{m}).Total.BackwardMean(zero_jump_index),rt_data.(fields{m}).Total.BackwardSEM(zero_jump_index), 'go');
        cmf_pred = [t2pred{m}(1:250) t1pred{m}(251:501)];
        plot(x_fit,cmf_pred,'b');
        xlabel('Target Shift (degrees)');
        ylabel('Reaction time (ms)');
        title(strrep(fields{m},'Mean_rt_',''),'interpreter','none');
        % legend(curvenames{m},'Location','northwest');
        % title([subject_list{n} ' ' curvenames{m} ' Manual Response Times'])
    end

    cmfpred1 = [t2pred{1}(1:250) t1pred{1}(251:501)];
    cmfpred2 = [t2pred{2}(1:250) t1pred{2}(251:501)];

    yrt_vec4 = cmf{1}';
    yrt_est_vec4 = cmfpred1(x_fit_idx);
    yrt_vec8 = cmf{2}';
    yrt_est_vec8 = cmfpred2(x_fit_idx);
    kk=4;
    RTBICorig.Amp4(n) = computeBIC(yrt_vec4,yrt_est_vec4,kk,length(X));
    RTBICorig.Amp8(n) = computeBIC(yrt_vec8,yrt_est_vec8,kk,length(X));

    sgtitle([cur_name ' Manual Response Times'])
    saveas(p, [path_RT_plots cur_name '_manual_rt.jpg']);
    close(p)
    writematrix ([bias_DDM threshold_DDM], [path_DDM_output subject_list{n} '_bias_threshold_DDM.csv']);
    All_bias(:,n) = bias_DDM';
    All_threshold(:,n) = threshold_DDM';
    
    %pref = {'Bias','Threshold','DDM_dr','DDM_bound','DDM_bias','DDM_ndt'};
    %suffix = {'amp4_All','amp8_All'};
    %write to file
    Str_format = '%s';
    Total_Cols = Nb_Params*2;
    for inn = [1:Total_Cols]
        Str_format = [Str_format,'\t%2.4f'];
    end
    Str_format = [Str_format,'\t%2.4f\r\n'];
    fprintf(fileID,Str_format,cur_name,cur_accuracy,bias_DDM,threshold_DDM,All_DR(:,n)',All_DDMbound(:,n)',All_DDMbias(:,n)',All_DDMndt(:,n)',All_veW(:,n)',All_A_bound(:,n)',All_B_bound(:,n)',All_RT_s1(:,n)',All_RT_smin1(:,n)',All_RT_ratio_s1(:,n)',All_RT_ratio_smin1(:,n)',errFits(1),errFits(2),BICerr(1),BICerr(2),AICerr(1),AICerr(2));
end

bic_struc = fullfile(path_DDM_output, ['BICorig.mat']);
save(bic_struc,'BICorig')
rtbic_struc = fullfile(path_DDM_output, ['RTBICorig.mat']);
save(rtbic_struc,'RTBICorig')

function [BICret] = computeBIC(y,yhat,k,N)
    indices = find(isnan(y) == 1);
    y(indices) = [];
    yhat(indices) = [];
    RSS = sum((y-yhat).^2);
    N=N-length(indices);
    BICret = N*log((RSS/N))+k*log(N);
end