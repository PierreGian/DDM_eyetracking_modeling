function [no_saccade_counts,no_saccade_percent,good_saccade_percent,Total_trial_nb] = get_saccadecount_behave(subj,saccades,trials_count,path_output)
   % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 11/2023
   % The purpose of this script is to compute the proportions of saccades
   % that are considered valid + return the proportions of "no saccades"
   % based on the visual errors. A *.csv file is written to "output" and
   % the proportions and counts are returned
   
   output_fpath = [path_output 'saccade_count_behave.csv'];
   file = fullfile(output_fpath);
   
   %if file does not exist, then create file in path output
   if(~isfile(file))
       tabledata = cell(1,6);
       tabledata{1,1} = 'SubjectID';
       tabledata{1,2} = 'goodsac_nb';
       tabledata{1,3} = 'Perc_goodsac';
       tabledata{1,4} = 'no_VE';
       tabledata{1,5} = 'Perc_noVE';
       tabledata{1,6} = 'Total_trials';
       writecell(tabledata, output_fpath);
   end
   
   Total_trial_nb = size(saccades.visualerror,2);
   %computing numbers/percentages of "no saccades" and "good saccades"
   saccade_VE = saccades.visualerror;
   no_saccade_counts = sum(isnan(saccade_VE));
   no_saccade_percent = (no_saccade_counts/Total_trial_nb)*100;
   good_saccade_percent = (trials_count/Total_trial_nb)*100;
   
   %only add data to table if subject was not already referenced
   curTable = readcell(file);
   rowNb = size(curTable,1);
   oldSubject=0;
   for r=2:rowNb %skip header
      if(strcmp(curTable{r,1},subj))
          oldSubject=1; %identify old subjects
          break
      end
   end
   
   NewIdx = rowNb+1; %index of new row
   if(~oldSubject) %only add data corresponding to new subjects
       curTable{NewIdx,1} = subj;
       curTable{NewIdx,2} = num2str(trials_count);
       curTable{NewIdx,3} = num2str(good_saccade_percent);
       curTable{NewIdx,4} = num2str(no_saccade_counts);
       curTable{NewIdx,5} = num2str(no_saccade_percent);
       curTable{NewIdx,6} = num2str(Total_trial_nb);
       writecell(curTable, output_fpath);
   end
   
end

