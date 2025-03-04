function [] = write_Xcount_behave(subj,subj_data,X,path_output)
   % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 12/2023
   % The purpose of this script is to write the trial count outputs. 
   % A *.csv file is written to "output" and the proportions and counts are returned
   
   output_fpath = [path_output 'X_count/' subj '_X_count_behave.csv'];
   file = fullfile(output_fpath);
   
   %Overwriting file for each subject
   tabledata = cell(1,4);
   tabledata{1,1} = 'X';
   tabledata{1,2} = 'Correct_Resp_Count';
   tabledata{1,3} = 'X_counts';
   tabledata{1,4} = 'Percent_Correct';
   writecell(tabledata, output_fpath);
   
   correct_resp_ct = subj_data.counts.correct_resp_counts;
   X_count = subj_data.counts.X_counts;
   Percent_correct = subj_data.counts.percent_correct;
   
   for xx = [1:length(X)]
       NewIdx = xx+1; %index of new row
       tabledata{NewIdx,1} = num2str(X(xx));
       tabledata{NewIdx,2} = num2str(correct_resp_ct(xx));
       tabledata{NewIdx,3} = num2str(X_count(xx));
       tabledata{NewIdx,4} = num2str(Percent_correct(xx));
       writecell(tabledata, output_fpath);
   end
   
end