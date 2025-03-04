function [ToReturn,saccades,totaltrialscount] = load_crunch_data_multiplefiles(data_input)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 11/2023
       % The purpose of this script is to load data from the eye-tracking
       % crunch *.mat objects that were preprocessed - multiple split files
       % case
       All_structs = load(data_input);
       data_labels = fieldnames(All_structs.all_data);
       info_labels = fieldnames(All_structs.all_info);
       struct_nb = size(data_labels,1);
       All_saccades = struct();
       ToReturn = struct();
       ToReturn.data = struct();
       All_counts = [];
       for dd=1:struct_nb
          data.data = All_structs.all_data.(data_labels{dd});
          data.info = All_structs.all_info.(info_labels{dd});
          saccades_fieldname = ['saccades' int2str(dd)];
          [cur_saccades,cur_trials_count] = load_crunch_data_singlefile(data,1);
          All_saccades.(saccades_fieldname) = cur_saccades;
          All_counts(end+1) = cur_trials_count;
          data_fields = fieldnames(data.data);
          for ff=[1:length(data_fields)]
             ToReturn.data.(data_fields{ff}) = data.data.(data_fields{ff});
          end
       end
       saccade_labels = fieldnames(All_saccades);
       saccades.visualerror_pixel_raw = All_saccades.(saccade_labels{1}).visualerror_pixel_raw;
       saccades.visualerror_pixel_filtered = All_saccades.(saccade_labels{1}).visualerror_pixel_filtered;
       saccades.visualerror_pixel_corrected = All_saccades.(saccade_labels{1}).visualerror_pixel_corrected;
       saccades.visualerror_degree_raw = All_saccades.(saccade_labels{1}).visualerror_degree_raw;
       saccades.visualerror_degree_filtered = All_saccades.(saccade_labels{1}).visualerror_degree_filtered;
       saccades.visualerror_degree_corrected = All_saccades.(saccade_labels{1}).visualerror_degree_corrected;
       saccades.visualerror = All_saccades.(saccade_labels{1}).visualerror;
       saccades.gaze_degree_pos = All_saccades.(saccade_labels{1}).gaze_degree_pos;
       saccades.exp_saccade2_degree = All_saccades.(saccade_labels{1}).exp_saccade2_degree;
       saccades.visual_error_endsac_degree_filtered = All_saccades.(saccade_labels{1}).visual_error_endsac_degree_filtered;
       saccades.valid_idx = All_saccades.(saccade_labels{1}).valid_idx;
       saccades.primary_amp_idx = All_saccades.(saccade_labels{1}).primary_amp_idx;
       saccades.is_saccade_event = All_saccades.(saccade_labels{1}).is_saccade_event;
       for dd=2:length(saccade_labels)
          cur_length = length(All_saccades.(saccade_labels{dd}).visualerror);
          for jj=1:cur_length
              saccades.visualerror_pixel_raw(end+1) = All_saccades.(saccade_labels{dd}).visualerror_pixel_raw(jj);
              saccades.visualerror_pixel_filtered(end+1) = All_saccades.(saccade_labels{dd}).visualerror_pixel_filtered(jj);
              saccades.visualerror_pixel_corrected(end+1) = All_saccades.(saccade_labels{dd}).visualerror_pixel_corrected(jj);
              saccades.visualerror_degree_raw(end+1) = All_saccades.(saccade_labels{dd}).visualerror_degree_raw(jj);
              saccades.visualerror_degree_filtered(end+1) = All_saccades.(saccade_labels{dd}).visualerror_degree_filtered(jj);
              saccades.visualerror_degree_corrected(end+1) = All_saccades.(saccade_labels{dd}).visualerror_degree_corrected(jj);
              saccades.visualerror(end+1) = All_saccades.(saccade_labels{dd}).visualerror(jj);
              saccades.gaze_degree_pos(end+1,:) = All_saccades.(saccade_labels{dd}).gaze_degree_pos(jj,:);
              saccades.exp_saccade2_degree(end+1) = All_saccades.(saccade_labels{dd}).exp_saccade2_degree(jj);
              saccades.visual_error_endsac_degree_filtered(end+1) = All_saccades.(saccade_labels{dd}).visual_error_endsac_degree_filtered(jj);
              saccades.valid_idx(end+1) = All_saccades.(saccade_labels{dd}).valid_idx(jj);
              saccades.primary_amp_idx(end+1) = All_saccades.(saccade_labels{dd}).primary_amp_idx(jj);
              saccades.is_saccade_event(end+1) = All_saccades.(saccade_labels{dd}).is_saccade_event(jj);
          end
%%%%%% valid_idx field used to be sorted by block and trial - not anymore          
%           validIDX_fields = fieldnames(All_saccades.(saccade_labels{dd}).valid_idx);
%           for ff=1:length(validIDX_fields)
%               saccades.valid_idx.(validIDX_fields{ff}) = All_saccades.(saccade_labels{dd}).valid_idx.(validIDX_fields{ff});
%           end
       end
       totaltrialscount = sum(All_counts);
end