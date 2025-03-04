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
path_plots = [master_path,'plots/DDM/'];
path_output_VE = [master_path 'outputs/VE/'];
path_RT_plots = [master_path,'plots/DDM/RT/'];
behave_data_folder = {dir(path_edf).name}; %Folder with all subject names

%Outputting the paths corresponding to edf, csv, and crunch files
disp(['path_edf is ' path_edf])
disp(['path_csv is ' path_csv])
disp(['path_crunch is ' path_crunch])

%Vector of jump sizes - input required for modeling analyses
%X=[-2.5 -1.75 -1.25 -1 -0.75 -0.5 -0.25 -0.175 -0.1 0 0.1 0.175 0.25 0.5 0.75 1 1.25 1.75 2.5];
X=[-2.5 -1.75 -1 -0.5 -0.25 0 0.25 0.5 1 1.75 2.5];
Amp_fields = {'amp4','amp8'};
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

amp4_filename = 'Amp4_data.txt';
amp8_filename = 'Amp8_data.txt';

amp4_dir = dir([path_output '/' amp4_filename]);
amp8_dir = dir([path_output '/' amp8_filename]);

amp4_data = readtable([path_output '/' amp4_dir(1).name],'Delimiter','\t');
amp8_data = readtable([path_output '/' amp8_dir(1).name],'Delimiter','\t');
allData.amp4_data = amp4_data;
allData.amp8_data = amp8_data;

accuracy = zeros(N,2);
bias = zeros(N,2);
threshold = zeros(N,2);
DDM_bias = zeros(N,2);
DDM_bound = zeros(N,2);
DDM_dr = zeros(N,2);
DDM_ndt = zeros(N,2);

for aa=[1:2]
    cur_amp = Amp_fields{aa};
    cur_table_name = [cur_amp '_data'];
    cur_table = allData.(cur_table_name);
    accuracy(:,aa) = table2array(cur_table(:,2));
    bias(:,aa) = table2array(cur_table(:,3));
    threshold(:,aa) = table2array(cur_table(:,4));
    DDM_dr(:,aa) = table2array(cur_table(:,5));
    DDM_bound(:,aa) = table2array(cur_table(:,6));
    DDM_bias(:,aa) = table2array(cur_table(:,7));
    DDM_ndt(:,aa) = table2array(cur_table(:,8));
end

ws = zeros(N,length(X)+3);

for n=1:length(subject_list)
    %for n = Select_idx
    %for n=1:5
    cur_name = subject_list{n};
    disp(['running ' cur_name])
    mat_data = Subject_Data.(cur_name).TaskData;
    rt_data = mat_data.mat_RTs;
    acc_data = mat_data.mat_accuracy;
    cur_amp = mat_data.mat_amps;
    ve_data = Subject_Data.(cur_name).saccades.visualerror;
    trial_len = length(ve_data);
    cur_bias_amp4 = bias(n,1);
    cur_bias_amp8 = bias(n,2);
    cur_ndt_amp4 = DDM_ndt(n,1);
    cur_ndt_amp8 = DDM_ndt(n,2);
    Jumpdir = mat_data.mat_jumpdir;
    Jumpsize = mat_data.mat_jumpsize;

    % index forward jumps and assign backwards jumps with (-) value
    forward_idx=strcmp(Jumpdir,'F'); %identify forward saccade jumps
    Jumpsize_relative = Jumpsize;
    Jumpsize_relative(~forward_idx)=Jumpsize_relative(~forward_idx)*-1;

    Xmat = zeros(trial_len,length(X)+3);
    Xmat(:,1) = ones(trial_len,1);
    ymat = acc_data';
    %ymat = rt_data';
    for tt=[1:trial_len]
        cAmp = cur_amp(tt);
        cShift = Jumpsize_relative(tt);
        cVE = ve_data(tt);
        if(cAmp==4)
            Xmat(tt,2) = cur_bias_amp4;
            %ymat(tt) = ymat(tt) - cur_ndt_amp4;
        else
            Xmat(tt,3) = cur_bias_amp8;
            %ymat(tt) = ymat(tt) - cur_ndt_amp8;
        end
        shift_idx = find(X==Jumpsize_relative(tt));
        col_shift_idx = shift_idx+3;
        if(~isnan(cVE))
            Xmat(tt,col_shift_idx) = cVE;
        end
    end

    w_fit = Xmat\ymat;
    ws(n,:) = w_fit;
end

w_acc = fullfile(path_output_VE, ['weights_acc.mat']);
save(w_acc,'ws')
%w_rt = fullfile(path_output_VE, ['weights_RT.mat']);
%save(w_rt,'ws')