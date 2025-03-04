function [ToReturn,gain] = extract_task_information(subj_name,data,saccades,block_nb)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 12/2023
       % The purpose of this script is to load subject data and extract all the task-relevant information
       cur_name = subj_name;
       ToReturn = struct();
       ToReturn.Block_names = {};
       ToReturn.mat_RTs = [];
       ToReturn.mat_response = [];
       ToReturn.mat_dirs = {};
       ToReturn.mat_amps = [];
       ToReturn.mat_jumpdir = {};
       ToReturn.mat_jumpsize = [];
       ToReturn.mat_detected = [];
       ToReturn.mat_accuracy = [];
       ToReturn.mat_forwresp = [];
       gain = [];
       
       count_trial = 0;
       for bb=[1:block_nb]
            Block_names{1,bb} = ['Block_',num2str(bb)];
            Trial_names = fieldnames(data.(Block_names{1,bb}));
            for tt=[1:length(Trial_names)]
                cur_trial = data.(Block_names{1,bb}).(Trial_names{tt});
                count_trial = count_trial+1;
                %Computation of good saccade indices
                ToReturn.mat_detected(end+1) = saccades.valid_idx(count_trial);

                if(isnan(cur_trial.manual_response))
                   ToReturn.mat_RTs(end+1) = 0;
                   ToReturn.mat_response(end+1) = 0;
                else
                   ToReturn.mat_RTs(end+1) = cur_trial.targ2_response_time-cur_trial.targ2_onset_time;
                   ToReturn.mat_response(end+1) = cur_trial.manual_response;
                end
                if(~isnan(cur_trial.targ1_pos_degree))
                    if(cur_trial.targ1_pos_degree>0)
                       ToReturn.mat_dirs{end+1} = 'R';
                       ToReturn.mat_amps(end+1) = fix(cur_trial.targ1_pos_degree);
                    else
                       ToReturn.mat_dirs{end+1} = 'L';
                       ToReturn.mat_amps(end+1) = fix(-1*cur_trial.targ1_pos_degree);
                    end
                else
                    ToReturn.mat_dirs{end+1} = NaN;
                end
                if(~isnan(cur_trial.targ2_jump_direction{1}))
                    ToReturn.mat_jumpdir{end+1} = cur_trial.targ2_jump_direction{1};
                    ToReturn.mat_jumpsize(end+1) = cur_trial.targ2_jump_size;
                else
                    ToReturn.mat_jumpdir{end+1} = NaN;
                    ToReturn.mat_jumpsize(end+1) = NaN;
                end
            end
       end
       %fix manual responses for NC_101 & NP 176- reversed
       if(strcmp(cur_name,'NC101') || strcmp(cur_name,'NP176')) %for some reason, responses need to be reversed here...
            rep1_idx = find(ToReturn.mat_response==1);
            rep4_idx = find(ToReturn.mat_response==4);
            ToReturn.mat_response(rep1_idx) = 4;
            ToReturn.mat_response(rep4_idx) = 1;
       end
       ToReturn.Block_names = Block_names;
       
       %Next extract accuracy info + forward responses
       target=-1;
       for cc=[1:length(ToReturn.mat_response)]
           if(~(ToReturn.mat_jumpsize(cc)==0)) 
               if((strcmp(ToReturn.mat_dirs{cc},'L') && strcmp(ToReturn.mat_jumpdir{cc},'F')) || (strcmp(ToReturn.mat_dirs{cc},'R') && strcmp(ToReturn.mat_jumpdir{cc},'B')))
                   target=1; %Left forward case & right backward case
               elseif((strcmp(ToReturn.mat_dirs{cc},'L') && strcmp(ToReturn.mat_jumpdir{cc},'B')) || (strcmp(ToReturn.mat_dirs{cc},'R') && strcmp(ToReturn.mat_jumpdir{cc},'F')))
                   target=4; %Left backward case & right forward case
               end
               ToReturn.mat_accuracy(cc) = (target==ToReturn.mat_response(cc));

               if(ToReturn.mat_detected(cc) && ~(ToReturn.mat_response(cc)==0)) %only consider good saccade indices
                   %Compute forward responses
                   if((strcmp(ToReturn.mat_dirs{cc},'L') && ToReturn.mat_response(cc)==1) || (strcmp(ToReturn.mat_dirs{cc},'R') && ToReturn.mat_response(cc)==4))
                       ToReturn.mat_forwresp(cc) = 1; %Two forward cases
                   elseif((strcmp(ToReturn.mat_dirs{cc},'L') && ToReturn.mat_response(cc)==4) || (strcmp(ToReturn.mat_dirs{cc},'R') && ToReturn.mat_response(cc)==1))
                       ToReturn.mat_forwresp(cc) = 0; %Two backward cases
                   end
               else
                   ToReturn.mat_forwresp(cc) = 2; 
               end
           else
               if(ToReturn.mat_response(cc)>0) %there needs to be a valid response for 0 jumpsize case
                    ToReturn.mat_accuracy(cc) = 1; %accuracy will be 1 regardless of response
                    %Compute forward responses
                    if(ToReturn.mat_detected(cc))
                        if((strcmp(ToReturn.mat_dirs{cc},'L') && ToReturn.mat_response(cc)==1) || (strcmp(ToReturn.mat_dirs{cc},'R') && ToReturn.mat_response(cc)==4))
                           ToReturn.mat_forwresp(cc) = 1; %Two forward cases
                        elseif((strcmp(ToReturn.mat_dirs{cc},'L') && ToReturn.mat_response(cc)==4) || (strcmp(ToReturn.mat_dirs{cc},'R') && ToReturn.mat_response(cc)==1))
                           ToReturn.mat_forwresp(cc) = 0; %Two backward cases
                        end
                    else
                        ToReturn.mat_forwresp(cc) = 2;
                    end
               else
                    ToReturn.mat_accuracy(cc) = 0; %if no response, then accuracy is 0
                    ToReturn.mat_forwresp(cc) = 2;
               end
           end

       end
       
       %Compute gain within subjects
       gaze_pos_x = saccades.gaze_degree_pos(:,1)';
       Trial_len = length(gaze_pos_x);
       target_pos = ToReturn.mat_amps;
       dirs_encode = strcmp(ToReturn.mat_dirs,'R')-strcmp(ToReturn.mat_dirs,'L');
       signed_target_pos = target_pos.*dirs_encode;
       for pp = [1:Trial_len]
           if(~isnan(gaze_pos_x(pp)) && ~isnan(target_pos(pp)))
               gain(1,pp) = (gaze_pos_x(pp)/signed_target_pos(pp))*100;
           else
               gain(1,pp) = NaN;
           end
       end
end