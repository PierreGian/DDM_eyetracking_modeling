function [RTs,choice_data,idx_include,RT_idx,Norm_RTs] = get_forwardcurve_data(x_relative,mat_subj,analysis_type,blocknb,trialnb)
    % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 11/2023
    % The purpose of this script is to generate the experimental data that
    % will be required for the forward curve modeling scripts
    % Two different possible cases: "VE" vs. "VECD" analysis types
    
    RT_range = [200 1200]; %RT bounds for this experiment
    VisualError=mat_subj.saccades.visualerror; %1x312 vector
    VisualError(isnan(VisualError))=0;
    ToPrint = ['Dealing with ',analysis_type,' case'];
    disp(ToPrint)
    saccade_idx = mat_subj.saccades.valid_idx;
    Responses = mat_subj.TaskData.mat_response;
    Accuracy = mat_subj.TaskData.mat_accuracy;
    Direction = mat_subj.TaskData.mat_dirs;
    Norm_RTs = zeros(1,length(Responses));
    RT = mat_subj.TaskData.mat_RTs;
    RT_toofast = RT<RT_range(1);
    RT_tooslow = RT>RT_range(2);
    RT(RT_toofast)=0; %filter RTs based on RT range
    RT(RT_tooslow)=0;
    ValidResponse = (Responses>0); %critical for choice_data
    
    %Computation of valid indices for RT mean/SE extraction
    %all have to have true values (not zero)
    idx_valid = saccade_idx(1:blocknb*trialnb) & RT(1:blocknb*trialnb) & Accuracy(1:blocknb*trialnb) & VisualError(1:blocknb*trialnb); %only use data with detected saccades, accurate response times, & correct trials
    
    %Computation of valid indices for choice_data
    idx_include = saccade_idx(1:blocknb*trialnb) & ValidResponse(1:blocknb*trialnb) & VisualError(1:blocknb*trialnb); %choice_data indices to include %Potentially add RT as criterion as well??  

    %index variables to get all true values in vector
    Direction=mat_subj.TaskData.mat_dirs(1:blocknb*trialnb);
    Amplitude=mat_subj.TaskData.mat_amps(1:blocknb*trialnb);
    Jumpdir=mat_subj.TaskData.mat_jumpdir(1:blocknb*trialnb);
    Jumpsize=mat_subj.TaskData.mat_jumpsize(1:blocknb*trialnb);
    Resp=mat_subj.TaskData.mat_response(1:blocknb*trialnb);
    edat_forward=mat_subj.TaskData.mat_forwresp(1:blocknb*trialnb); %1 means forward, 0 means backward, 2 means excluded trials

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
                 choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwResp(1,xx) = sum(c_dir&c_amp&idx_include&(Jumpsize_relative==x_relative(xx))&edat_forward_idx);
                 choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwPercent(1,xx) = choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwResp(xx)/choice_data.(Amp_fields{f}).(Dir_fields{ff}).ForwSumTrials(xx);
            end
        end
        last_dir = Dir_fields{3};
        choice_data.(Amp_fields{f}).(last_dir).ForwSumTrials = choice_data.(Amp_fields{f}).(Dir_fields{1}).ForwSumTrials+choice_data.(Amp_fields{f}).(Dir_fields{2}).ForwSumTrials;
        choice_data.(Amp_fields{f}).(last_dir).ForwResp = choice_data.(Amp_fields{f}).(Dir_fields{1}).ForwResp+choice_data.(Amp_fields{f}).(Dir_fields{2}).ForwResp;
        choice_data.(Amp_fields{f}).(last_dir).ForwPercent = choice_data.(Amp_fields{f}).(last_dir).ForwResp./choice_data.(Amp_fields{f}).(last_dir).ForwSumTrials;
     end
     
     RTs = struct(); %returning the struct RTs
     %Struct with 2 x 3 x 3 cases (4/8 degrees, Left/Right/Total, Forward/Backward,All) 
     Xsize = length(x_relative);
     RTs.Amp4 = struct();
     RTs.Amp8 = struct();
     for f=1:length(Amp_fields)
        RTs.(Amp_fields{f}).Left = struct();
        RTs.(Amp_fields{f}).Right = struct();
        RTs.(Amp_fields{f}).Total = struct();
        for ff=1:length(Dir_fields)
            RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardMean = zeros(1,Xsize);
            RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardSEM = zeros(1,Xsize);
            RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardMean = zeros(1,Xsize);
            RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardSEM = zeros(1,Xsize);
            RTs.(Amp_fields{f}).(Dir_fields{ff}).FullMean = zeros(1,Xsize);
            RTs.(Amp_fields{f}).(Dir_fields{ff}).FullSEM = zeros(1,Xsize);
        end
     end

     %Order is Amp4/Left, Amp4/Right, Amp4/Total, Amp8/Left, Amp8/Right,
     %Amp8/Total
     RTsize = length(RT);
     Allbooleans = zeros(6,RTsize);
     Allbooleans(1,:) = Amp4_idx & left_idx;
     Allbooleans(2,:) = Amp4_idx & right_idx;
     Allbooleans(3,:) = Amp4_idx;
     Allbooleans(4,:) = Amp8_idx & left_idx;
     Allbooleans(5,:) = Amp8_idx & right_idx;
     Allbooleans(6,:) = Amp8_idx;

     cBool = 1;
     for f=1:length(Amp_fields)
         for ff=1:length(Dir_fields)
            curBool = Allbooleans(cBool,:);
            for k=[1:length(x_relative)]
                %Encoding vector corresponding to boolean + with the right
                %encoding X jumpsize
                currentForward_idx = logical(curBool & fresp_idx & (Jumpsize_relative==x_relative(k)));
                currentBackward_idx = logical(curBool & bresp_idx & (Jumpsize_relative==x_relative(k)));
                currentAll_idx = logical(currentForward_idx | currentBackward_idx);
                %1) Dealing with Forward case - get RTs at each direction, then
                %mean & SE
                Forward_RT = RT(currentForward_idx);
                Forward_RT = Forward_RT(Forward_RT>0); %only consider non-zero RTs
                if(~isempty(Forward_RT) && length(Forward_RT)>2) %we need at least three trials for RT computation
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardMean(k) = nanmean(Forward_RT);
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardSEM(k) = nanstd(Forward_RT)/sqrt(length(Forward_RT));
                else
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardMean(k) = NaN;
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).ForwardSEM(k) = NaN;
                end
                
                %2) Dealing with Backward case
                Backward_RT = RT(currentBackward_idx);
                Backward_RT = Backward_RT(Backward_RT>0); %only consider non-zero RTs
                if(~isempty(Backward_RT) && length(Backward_RT)>2) %we need at least three trials for RT computation
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardMean(k) = nanmean(Backward_RT);
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardSEM(k) = nanstd(Backward_RT)/sqrt(length(Backward_RT));
                else
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardMean(k) = NaN;
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).BackwardSEM(k) = NaN;
                end
                
                %3) Dealing with all case
                All_RT = RT(currentAll_idx);
                All_RT = All_RT(All_RT>0); %only consider non-zero RTs
                if(~isempty(All_RT) && length(All_RT)>2) %we need at least three trials for RT computation
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).FullMean(k) = nanmean(All_RT);
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).FullSEM(k) = nanstd(All_RT)/sqrt(length(All_RT));
                else
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).FullMean(k) = NaN;
                    RTs.(Amp_fields{f}).(Dir_fields{ff}).FullSEM(k) = NaN;
                end
            end
         cBool=cBool+1; %Increment boolean by 1 at end of loop
     end
     
     cur_amp = Amp_vals(f);
     mid_idx = (length(x_relative)-1)/2+1;
     Ref_0_RT = RTs.(Amp_fields{f}).Total.FullMean(mid_idx);
     cAmp_idx = find(Amplitude==cur_amp);
     for a_idx = [1:length(cAmp_idx)]
        cur_idx = cAmp_idx(a_idx);
        if(~isnan(RT(cur_idx)))
            if(RT(cur_idx)<RT_range(1) || RT(cur_idx)>RT_range(2) || idx_valid(cur_idx)==0)
                Norm_RTs(cur_idx) = 0;
            else
                Norm_RTs(cur_idx) = (RT(cur_idx)/Ref_0_RT)*100;
            end
        else
            Norm_RTs(cur_idx) = NaN;
        end
     end
     end
     
end

