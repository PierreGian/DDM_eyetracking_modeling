function [ToReturn] = extract_task_metrics(subj_name,subj_mat,X)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 12/2023
       % The purpose of this script is to load subject mat_data (struct
       % from "extract_task_information()") and extract metrics that may be
       % used for further analyses
       % INPUTS
       % "subj_name": subject name
       % "subj_mat": subject's task data
       % "X": List of all jumpsizes
       % OUTPUTs
       % ToReturn struct with count information pertaining to each jumpsize
       
       ToReturn = struct();
       ToReturn.counts = struct();
       ToReturn.saccades = struct();
       subj_data = subj_mat.TaskData;
       Accuracy = subj_data.mat_accuracy;
       Jumpsize = subj_data.mat_jumpsize;
       Jumpdir = subj_data.mat_jumpdir; %can be F or B
       Amplitude = subj_data.mat_amps;
       Direction = subj_data.mat_dirs; %can be L or R
       VE = subj_mat.saccades.visualerror;
       VE_altVE = subj_mat.saccades.visual_error_endsac_degree_filtered;
       Sac_amp = subj_mat.saccades.sac_amp;
       Sac_vel = subj_mat.saccades.sac_vel;
       Sac_acc = subj_mat.saccades.sac_acc;
       
       forward_idx=strcmp(Jumpdir,'F'); %identify forward saccade jumps
       Jumpsize_relative = Jumpsize;
       Jumpsize_relative(~forward_idx)=Jumpsize_relative(~forward_idx)*-1;
       
       Jumpsize_counts = zeros(1,length(X));
       correct_counts = zeros(1,length(X));
       percent_correct = zeros(1,length(X));
       Jumpsize_VE_all = zeros(8,length(X));
       Jumpsize_VE_all2 = zeros(8,length(X));
       Jumpsize_Samp_all = zeros(8,length(X));
       Jumpsize_Svel_all = zeros(8,length(X));
       Jumpsize_Sacc_all = zeros(8,length(X));
       VE_bools = [Amplitude==4 & strcmp(Direction,'L') & strcmp(Jumpdir,'F');Amplitude==4 & strcmp(Direction,'L') & strcmp(Jumpdir,'B');Amplitude==4 & strcmp(Direction,'R') & strcmp(Jumpdir,'F');Amplitude==4 & strcmp(Direction,'R') & strcmp(Jumpdir,'B')]; %1) amp4 left forward; 2) amp4 left backward; 3) amp4 right forward; 4) amp4 right backward
       VE_bools = [VE_bools;Amplitude==8 & strcmp(Direction,'L') & strcmp(Jumpdir,'F');Amplitude==8 & strcmp(Direction,'L') & strcmp(Jumpdir,'B');Amplitude==8 & strcmp(Direction,'R') & strcmp(Jumpdir,'F');Amplitude==8 & strcmp(Direction,'R') & strcmp(Jumpdir,'B')];%5) amp8 left forward; 6) amp8 left backward; 7) amp8 right forward; 8) amp8 right backward
       for xx=[1:length(X)]
          cur_X = X(xx);
          Jumpsize_counts(xx) = sum(Jumpsize_relative==cur_X);
          correct_counts(xx) = sum((Jumpsize_relative==cur_X) & Accuracy);
          percent_correct(xx) = correct_counts(xx)/Jumpsize_counts(xx);
          for bb=[1:size(VE_bools,1)]
              VE_cur = VE(VE_bools(bb,:) & Jumpsize_relative==cur_X);
              VE_cur(isnan(VE_cur))=[];
              Jumpsize_VE_all(bb,xx) = mean(VE_cur);
              VE_cur2 = VE_altVE(VE_bools(bb,:) & Jumpsize_relative==cur_X);
              VE_cur2(isnan(VE_cur2))=[];
              Jumpsize_VE_all2(bb,xx) = mean(VE_cur2);
              Samp_cur = Sac_amp(VE_bools(bb,:) & Jumpsize_relative==cur_X);
              Svel_cur = Sac_vel(VE_bools(bb,:) & Jumpsize_relative==cur_X);
              Sacc_cur = Sac_acc(VE_bools(bb,:) & Jumpsize_relative==cur_X);
              Samp_cur(isnan(Samp_cur))=[];
              Svel_cur(isnan(Svel_cur))=[];
              Sacc_cur(isnan(Sacc_cur))=[];
              Jumpsize_Samp_all(bb,xx) = mean(Samp_cur);
              Jumpsize_Svel_all(bb,xx) = mean(Svel_cur);
              Jumpsize_Sacc_all(bb,xx) = mean(Sacc_cur);
          end
       end
       
       ToReturn.counts.X_counts = Jumpsize_counts;
       ToReturn.counts.correct_resp_counts = correct_counts;
       ToReturn.counts.percent_correct = percent_correct;
       ToReturn.saccades.VE_bools = VE_bools;
       ToReturn.saccades.VE_all = Jumpsize_VE_all;
       ToReturn.saccades.VE_all_altVE = Jumpsize_VE_all2;
       ToReturn.saccades.Sac_amp = Jumpsize_Samp_all;
       ToReturn.saccades.Sac_vel = Jumpsize_Svel_all;
       ToReturn.saccades.Sac_acc = Jumpsize_Sacc_all;
end