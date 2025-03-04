function [saccades,trials_count] = load_crunch_data_singlefile(data,split)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 11/2023
       % The purpose of this script is to load data from the eye-tracking
       % crunch *.mat objects that were preprocessed - single file case
       % Note: this script does not take into account the 10 practice
       % trials
       if(split==0) %load file if not split file, otherwise, just pass data struct
            data = load(data);
       end
       saccades = struct();
       data_labels = fieldnames(data.data);
       check_block0 = @(x) contains(x,'_0'); %identify Block 0, Block 1-6 don't have a 0 in front
       Blocks_idx = ~(cellfun(check_block0,data_labels)); %returns experimental block idx
       Block_labels = data_labels(Blocks_idx)'; %names of experimental blocks
       Block_Nb = size(Block_labels,2);
       empty_saccades = cell(length(Block_labels),1);
       saccades.valid_idx = cell2struct(empty_saccades,Block_labels,1); %initialize valid saccade idx booleans
       Trial_Names = {};
       for bb=1:Block_Nb
            Block_data(bb) = data.data.(Block_labels{1,bb});
            Trial_Names{bb} = fieldnames(Block_data(bb))';
       end
       All_trial_nb = size(data.info.primary_saccade_amp,1);
       %Need to extract saccades with valid indices - need two pieces of
       %information: 1) Primary saccade amplitude, 2) visual error index
       %Rule is that BOTH need to not be NaNs
       Trial_count = 0;
       saccades.primary_amp_idx = [];
       for bb=1:Block_Nb
           Trial_Nb = size(Trial_Names{bb},2);
           saccades.valid_idx.(Block_labels{1,bb}) = zeros(1,Trial_Nb);
           for tt=1:Trial_Nb
                curTrial = Block_data(bb).(Trial_Names{bb}{tt}); %Extract current trial info
                Trial_count=Trial_count+1; %increment by 1
                Saccade_amp = curTrial.primary_saccade_amp; %Extract primary saccade amplitude
                Visual_error = curTrial.visual_error_corrected; %Extract Visual error
                if(~isnan(Saccade_amp) && ~isnan(Visual_error)) %both need to be valid
                    saccades.valid_idx.(Block_labels{1,bb})(1,tt) = 1;
                    saccades.primary_amp_idx(end+1) = 1;
                else
                    saccades.valid_idx.(Block_labels{1,bb})(1,tt) = 0;
                    saccades.primary_amp_idx(end+1) = 0;
                end
           end
       end
       Practice_count = All_trial_nb-Trial_count; %should be 10
       first_saccade_idx = Practice_count+1; %start counting at index 11
       
       % RULE: saccade must be detected by eyelink software -> no second
       % target for trial otherwise: extra layer of Quality Control
       
       % Message from SB: info.if_detect_saccade: this information is whether the crunch
       % code detects a valid saccade using the metrics that Weiwei put
       % into place - it's possible that the crunch code picks up a saccade
       % that the software did not pick up during the task - this is why
       % info.is_saccade_event is put into place as an extra QC measure to
       % make sure the software did detect the saccade during the task 
        
       % info.is_saccade_event: saccade detected by software during the
       % task ie. the task during the session detected the saccade and
       % displayed target 2 (this is seen when the software displays 'IStartSaccadeEvent'
       % - if the software does not detect a saccade during the task it will not display the second target 
       ss=first_saccade_idx;
       saccades.is_saccade_event = [];
       for bb=1:Block_Nb
           Trial_Nb = size(Trial_Names{bb},2);
           for tt=1:Trial_Nb
               saccades.is_saccade_event(end+1) = data.info.is_saccade_event(ss);
               if data.info.is_saccade_event(ss) == 0 
                    saccades.valid_idx.(Block_labels{1,bb})(1,tt) = 0; %switch to 0 if is_saccade_event() returns 0
               end
               ss=ss+1;
           end
       end
       
       % Message from SB: VisualError is the difference between the target2 location and the eye position at the time of target2 onset
       % 'visual error raw' is the raw data recorded by the eye tracker
       % a low-pass filter was used to filter out the low-frequency noise and generate 'filtered' data 
       % 'visual error filtered' is calculated based on the filtered data
       ss=first_saccade_idx;
       trials_count=0;
       saccades.visualerror_pixel_raw = [];
       saccades.visualerror_pixel_filtered = [];
       saccades.visualerror_pixel_corrected = [];
       saccades.visualerror_degree_raw = [];
       saccades.visualerror_degree_filtered = [];
       saccades.visualerror_degree_corrected = [];
       saccades.gaze_degree_pos = [];
       saccades.exp_saccade2_degree = [];
       saccades.visual_error_endsac_degree_filtered = [];
       saccades.sac_amp = [];
       saccades.sac_vel = [];
       saccades.sac_acc = [];
       y0=600;
       for bb=1:Block_Nb
           Trial_Nb = size(Trial_Names{bb},2);
           for tt=1:Trial_Nb
                if(saccades.valid_idx.(Block_labels{1,bb})(1,tt) == 0)
                    saccades.visualerror_pixel_raw(ss-Practice_count) = NaN;
                    saccades.visualerror_pixel_filtered(ss-Practice_count) = NaN;
                    saccades.visualerror_pixel_corrected(ss-Practice_count) = NaN;
                    saccades.visualerror_degree_raw(ss-Practice_count) = NaN;
                    saccades.visualerror_degree_filtered(ss-Practice_count) = NaN;
                    saccades.visualerror_degree_corrected(ss-Practice_count) = NaN;
                    saccades.gaze_degree_pos(ss-Practice_count,:) = [NaN,NaN];
                    saccades.exp_saccade2_degree(ss-Practice_count) = NaN;
                    saccades.visual_error_endsac_degree_filtered(ss-Practice_count) = NaN;
                    saccades.sac_amp(ss-Practice_count) = NaN;
                    saccades.sac_vel(ss-Practice_count) = NaN;
                    saccades.sac_acc(ss-Practice_count) = NaN;
                else
                    saccades.visualerror_pixel_raw(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_pixel_raw;
                    saccades.visualerror_pixel_filtered(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_pixel_filtered;
                    saccades.visualerror_pixel_corrected(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_corrected;
                    saccades.visualerror_degree_raw(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_degree_raw;
                    saccades.visualerror_degree_filtered(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_degree_filtered;
                    saccades.visualerror_degree_corrected(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_degree_corrected;
                    saccades.gaze_degree_pos(ss-Practice_count,:) = Block_data(bb).(Trial_Names{bb}{tt}).gaze_degree_pos(1,:);
                    saccades.exp_saccade2_degree(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).exp_saccade2_degree;
                    saccades.visual_error_endsac_degree_filtered(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).visual_error_endsac_degree_filtered;
                    saccades.sac_amp(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).primary_saccade_amp;
                    saccades.sac_vel(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).primary_saccade_peak_vel;
                    saccades.sac_acc(ss-Practice_count) = Block_data(bb).(Trial_Names{bb}{tt}).primary_saccade_peak_acc;
                    trials_count=trials_count+1;
                end
                ss=ss+1;
           end
       end
       saccades.visualerror = saccades.visualerror_degree_filtered;
       
       %Update format of saccades.valid_idx vector - 1x312 vector
       goodsaccade_idx = saccades.valid_idx;
       saccades.valid_idx = [];
       for bb=1:Block_Nb
           Trial_Nb = size(Trial_Names{bb},2);
           for tt=1:Trial_Nb
               saccades.valid_idx(end+1) = goodsaccade_idx.(Block_labels{1,bb})(1,tt);
           end
       end
       
end