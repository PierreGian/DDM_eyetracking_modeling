function [data,info] = New_Eye_Tracking_Crunch(Cdr,dist,path2crunch,ascii)
% Original script from Sonia Bansal, edited by Weiwei Zhou & Pierre
% Gianferrara - UC Davis, Sensorimotor integration lab 12/2023
%
%input - Cdr: directory that contains all asc files (Cdr = cd for current folder)
%       dist: distance from eye to monitor (mm)
%output - data: sorted data for each trial in the block
%         info: contains data combined for all trials
% data and info contain raw and filtered data. Raw data is the one
% retrieved from the asc files. Filtered data is the one filtered by a
% low-pass filter with specific cut off frequency.

%%%%% All variables - information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Save = 1; % To save the converted data as .mat file, this should be 1.
PacketSize = 20000; % 20000 rows will be processed at once.
x0=960; % x0 - x position at the center of the screen
%y0=600; % y0 - y position at the center of the screen
y0 = 600;
% dist=1450; % mm % Subject distance from the screen
% Screen Resolution
scn_res_width=1920; % pixels
scn_res_height=1200; % pixels
% Screen Dimensions
scb_dim_width=518.48;  % mm
scb_dim_height=324.05; % mm
sampling_freq = 1000; % Eye tracker sampling frequency
cut_off_freq = 0.2*sampling_freq; % low pass filter cut off frequency
%Filtering good trials vs. bad trials
time_threshold_after_targ1_onset = 300;

disp('--- [MasterCrunch.m] ------------------------------');
%%%%% Three .asc files per subject (A, B, and C) %%%%%%%%%%%%%%%%%%%%%%
data=struct;
info=struct;
tic;
files = dir([Cdr,ascii]);

for j=1:length(files)
    FileName = files(j).name;

    %% Load the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Loading the data from ' FileName]);
    fid = fopen([Cdr '/' FileName],'rt');
    Data = fread(fid,'char');
    fclose(fid);
    disp(' ...done.');
    disp('Splitting the data into rows');
    Data = char(Data');
    [Data, RowNum] = Text2CellArray(Data,PacketSize);
    disp(['...done. (' num2str(RowNum) ' rows)']);

    %% Classify the rows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Classifying the rows');
    % Picking-up the rows needed with specific data-types
    DataTypes = [
        'MS'; % < 1> MSG
        'ST'; % < 2> START
        'EN'; % < 3> END
        'SS'; % < 4> SSACC
        'ES'; % < 5> ESACC
        'SF'; % < 6> SFIX
        'EF'; % < 7> EFIX
        'SB'; % < 8> SBLINK
        'EB'; % < 9> EBLINK
        ];
    Category = zeros(1,RowNum);
    char_12 = char(Data);
    char_12 = char_12(:,1:2); % Take first 2 characters of every row

    for c=1:9 % for each category
        F = find((char_12(:,1)==DataTypes(c,1)) & (char_12(:,2)==DataTypes(c,2)));
        if c==2
            TrialNum = length(F);
            FS = F; % 'START' rows
        elseif c==3
            TrialNum = min([length(F) TrialNum]);
            FE = F; % 'END' rows
        end
        Category(F) = ones(1,length(F))*c;
    end

    Trial = {'Trial_'};
    TrialName=cell(1,TrialNum);
    for b = 1:TrialNum
        TrialName(b) =strcat(Trial,num2str(b)) ;
    end

    % Classifying the rows
    MSG = Data(Category== 1);
    START = Data(Category== 2);
    END = Data(Category== 3);
    Data_SSACC = Data(Category== 4);
    Data_ESACC = Data(Category== 5);
    Data_SFIX = Data(Category== 6);
    Data_EFIX = Data(Category== 7);
    Data_SBLINK = Data(Category== 8);
    Data_EBLINK = Data(Category== 9);

    % Picking-up data-only rows that is comprised of only values.
    if TrialNum==0
        disp('No block is recorded in this file.');
        return;
    else
        Data = char(Data);
        F = [];
        for b=1:TrialNum
            F = [F (FS(b)+1):(FE(b)-1)];
        end
        Data = Data(F,:);
    end

    DataRowChar = [double('0123456789.- ') 9]; % numbers, dot, minus, space and tab
    DataRowNum = size(Data,1);

    % To save memory, divide the data into several packets and process each
    Fall = [];
    for p=1:ceil(DataRowNum/PacketSize)
        S = (p-1)*PacketSize+1;
        E = p*PacketSize;
        if E>DataRowNum
            E = DataRowNum;
        end
        Packet = double(Data(S:E,:));
        temp = zeros(size(Packet,1),1);
        for c=1:length(DataRowChar)
            temp = temp + sum(Packet==DataRowChar(c),2);
        end
        F = find(temp==size(Packet,2));
        F = F + S - 1;
        Fall = [Fall; F];
    end
    clear temp F;
    Data = [Data(Fall,:)];
    clear Fall;
    % report
    disp(['... done. (' num2str(TrialNum) ' trials found)']);

    %% Determining the raw data of each block %%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Interpreting the raw data');
    % To save memory, divide the data into several packets and process each
    DataRowNum = size(Data,1);
    RawData = [];

    for p=1:ceil(DataRowNum/PacketSize)
        S = (p-1)*PacketSize+1;
        E = p*PacketSize;
        if E>DataRowNum
            E = DataRowNum;
        end
        Packet = [double(Data(S:E,:))]';
        L = numel(Packet);
        % Finding the null-data, represented by '.'
        % (e.g., gaze-positions are lost during blinks)
        NullData = find((Packet(1:L-1)==46)&(Packet(2:L)==9)); % Finding "dot plus tab"
        NullData = [NullData find((Packet(1:L-1)==46)&(Packet(2:L)==32))]; % Finding "dot plus space"
        NullData = [NullData find((Packet(1:L-2)==46)&(Packet(2:L-1)==46)&(Packet(3:L)==46))]; % Finding "dot plus dot plus dot"
        NullData = [NullData find((Packet(1:L)==46)&(mod([1:L],size(Packet,1))==0))]; % Finding dot at the end of the rows
        Packet(NullData) = ones(1,length(NullData))*48; % Replace the dot with zero.
        %Data(S:E,:) = char(Packet');
        A=char(Packet');
        DATA_Chunck=zeros(size(A,1),5);
        for nn=1:size(A,1)%length(A)
            TEMP_string=str2num([A(nn,:)]);
            %             if length(TEMP_string)~=8
            %                 TEMP_string=[TEMP_string,[0 0 0]];
            %             end
            DATA_Chunck(nn,:)=TEMP_string(1:5);
        end
        %         RawData = [RawData; str2num(Data(S:E,:))];
        RawData = [RawData; DATA_Chunck];
    end

    clear Packet NullData;
    clear Data;
    % report
    disp('... done.');

    %% Determining START/END of each block %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Detecting start/end timestamps of each block');
    % Find the timestamps of START and END
    for b=1:TrialNum
        temp = double(char(START(b)));
        F = find(temp==9); % Find tab
        temp = temp((F(1)+1):F(2));
        STARTtime(b,:) = [b str2num(char(temp))];
        temp = double(char(END(b)));
        F = find(temp==9); % Find tab
        temp = temp((F(1)+1):F(2));
        ENDtime(b,:) = [b str2num(char(temp))];
    end

    % Add block numbers to the data
    RawData = [zeros(size(RawData,1),1) RawData];
    for b=1:TrialNum
        F = find( (RawData(:,2)>=STARTtime(b,2)) & (RawData(:,2)<=ENDtime(b,2)) );
        RawData(F,1) = ones(length(F),1)*b;
    end
    % report
    disp('... done.');

    %% Convert the saccade data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Interpreting the saccade data');
    % SSACC: Start of saccades
    SSACC = [];
    if isempty(Data_SSACC)==0
        SSACC = double(char(Data_SSACC));
        SSACC = SSACC(:,7:size(SSACC,2));
        F = find(SSACC(:,1)==82); % find 'R' (i.e., right eye)
        SSACC(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        SSACC = char(SSACC);
        SSACC = str2num(SSACC);
    end
    % ESACC: End of saccades
    ESACC = [];
    if isempty(Data_ESACC)==0
        ESACC = double(char(Data_ESACC));
        ESACC = ESACC(:,7:size(ESACC,2));
        F = find(ESACC(:,1)==82); % find 'R' (i.e., right eye)
        ESACC(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        ESACC = char(ESACC);
        ESACC_temp = [];
        for m =1:size(ESACC,1)
            temp5=str2num(ESACC(m,:));
            if not(isempty(temp5))
                ESACC_temp=[ESACC_temp;temp5];
            else
                continue
            end
        end
        ESACC = ESACC_temp;
    end

    % Appending block numbers to the data
    SSACC = [zeros(size(SSACC,1),1) SSACC];
    ESACC = [zeros(size(ESACC,1),1) ESACC];
    for b=1:TrialNum
        F = find( (SSACC(:,3)>=STARTtime(b,2)) & (SSACC(:,3)<=ENDtime(b,2)) );
        SSACC(F,1) = ones(length(F),1)*b;
        F = find( (ESACC(:,4)>=STARTtime(b,2)) & (ESACC(:,4)<=ENDtime(b,2)) );
        ESACC(F,1) = ones(length(F),1)*b;
    end

    % report
    temp = max([size(SSACC,1) size(ESACC,1)]);
    disp(['...done. (' num2str(temp) ' saccades)']);

    %% Convert the fixation data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Interpreting the fixation data');
    % SFIX: start of fixation
    SFIX = [];
    if isempty(Data_SFIX)==0
        SFIX = double(char(Data_SFIX));
        SFIX = SFIX(:,6:size(SFIX,2));
        F = find(SFIX(:,1)==76); % find 'L' (i.e., left eye)
        SFIX(F,1) = ones(length(F),1)*48; % Convert the 'L' into '0'
        F = find(SFIX(:,1)==82); % find 'R' (i.e., right eye)
        SFIX(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        SFIX = str2num(char(SFIX));
    end
    % EFIX: end of fixation
    EFIX = [];
    if isempty(Data_EFIX)==0
        EFIX = double(char(Data_EFIX));
        EFIX = EFIX(:,6:size(EFIX,2));
        F = find(EFIX(:,1)==76); % find 'L' (i.e., left eye)
        EFIX(F,1) = ones(length(F),1)*48; % Convert the 'L' into '0'
        F = find(EFIX(:,1)==82); % find 'R' (i.e., right eye)
        EFIX(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        EFIX=char(EFIX);
        EFIX = str2num(EFIX);
    end
    % Appending block numbers to the data
    SFIX = [zeros(size(SFIX,1),1) SFIX];
    EFIX = [zeros(size(EFIX,1),1) EFIX];
    for b=1:TrialNum
        F = find( (SFIX(:,3)>=STARTtime(b,2)) & (SFIX(:,3)<=ENDtime(b,2)) );
        SFIX(F,1) = ones(length(F),1)*b;
        F = find( (EFIX(:,4)>=STARTtime(b,2)) & (EFIX(:,4)<=ENDtime(b,2)) );
        EFIX(F,1) = ones(length(F),1)*b;
    end

    % report
    temp = max([size(SFIX,1) size(EFIX,1)]);
    disp(['...done. (' num2str(temp) ' fixations)']);

    %% Convert the blink data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Interpreting the blink data');
    % SBLINK: Start of blinks
    SBLINK = [];
    if isempty(Data_SBLINK)==0
        SBLINK = double(char(Data_SBLINK));
        SBLINK = SBLINK(:,8:size(SBLINK,2));
        F = find(SBLINK(:,1)==76); % find 'L' (i.e., left eye)
        SBLINK(F,1) = ones(length(F),1)*48; % Convert the 'L' into '0'
        F = find(SBLINK(:,1)==82); % find 'R' (i.e., right eye)
        SBLINK(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        SBLINK = str2num(char(SBLINK));
    end
    % EBLINK: End of blinks
    EBLINK = [];
    if isempty(Data_EBLINK)==0
        EBLINK = double(char(Data_EBLINK));
        EBLINK = EBLINK(:,8:size(EBLINK,2));
        F = find(EBLINK(:,1)==76); % find 'L' (i.e., left eye)
        EBLINK(F,1) = ones(length(F),1)*48; % Convert the 'L' into '0'
        F = find(EBLINK(:,1)==82); % find 'R' (i.e., right eye)
        EBLINK(F,1) = ones(length(F),1)*49; % Convert the 'R' into '1'
        EBLINK = str2num(char(EBLINK));
    end
    % Appending block numbers to the data
    SBLINK = [zeros(size(SBLINK,1),1) SBLINK];
    EBLINK = [zeros(size(EBLINK,1),1) EBLINK];
    if ~isempty(SBLINK)
        for b=1:TrialNum
            F = find( (SBLINK(:,3)>=STARTtime(b,2)) & (SBLINK(:,3)<=ENDtime(b,2)) );
            SBLINK(F,1) = ones(length(F),1)*b;
            F = find( (EBLINK(:,4)>=STARTtime(b,2)) & (EBLINK(:,4)<=ENDtime(b,2)) );
            EBLINK(F,1) = ones(length(F),1)*b;
        end
    end

    % report
    temp = max([size(SBLINK,1) size(EBLINK,1)]);
    disp(['...done. (' num2str(temp) ' blinks)']);

    %% Convert the MSG (message) data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Interpreting the MSG data']);
    MSGtemp = double(char(MSG));
    temp = size(MSGtemp,1);
    MSGtime = [];
    if isempty(MSG)==0 % If there are MSG rows
        L = size(MSGtemp,2);
        MSGtemp2 = (MSGtemp==32); % find spaces
        for c=1:L
            mask(:,c) = sum(MSGtemp2(:,1:c),2);
        end
        mask = (mask==0);
        MSGtime = (MSGtemp.*mask) + (1-mask)*32; % Timestamps and spaces
        MSGtime = MSGtime(:,5:L); % Omit the "MSG + tab" of each row
        MSGtime = str2num(char(MSGtime));
        clear MSGtemp2 mask L;
    end

    % Appending block numbers to the data
    MSGtime_trial = [zeros(temp,1) MSGtime];
    for b=1:TrialNum
        F = find( (MSGtime_trial(:,2)>=STARTtime(b,2)) & (MSGtime_trial(:,2)<=ENDtime(b,2)) );
        MSGtime_trial(F,1) = ones(length(F),1)*b;
    end
    % report
    disp(['...done. (' num2str(temp) ' MSG rows)']);

    % Retrive information from MSG
    UniqueTrial = unique(RawData(:,1));
    Targ1ON = zeros(length(UniqueTrial)+1,3);
    Targ2ON = zeros(length(UniqueTrial)+1,3);
    Targ2OFF = zeros(length(UniqueTrial)+1,2);
    Targ2RESP = zeros(length(UniqueTrial)+1,2);
    Targ2ACC = zeros(length(UniqueTrial)+1,2);
    Targ1POS = zeros(length(UniqueTrial)+1,2);
    Targ2POS = zeros(length(UniqueTrial)+1,2);
    jump_direction = cell(length(UniqueTrial)+1,1);
    jump_size = zeros(length(UniqueTrial)+1,1);
    trial_in_block = zeros(length(UniqueTrial)+2,2);
    block_id = zeros(length(UniqueTrial)+1,2);
    saccade_event = zeros(length(UniqueTrial)+1,2);
    trial_index = 0;
    for i=1:length(MSG)
        temp=regexpi(cell2mat(MSG(1,i)),'\s','split');
        GG1=find(ismember(temp,'Target1')==1, 1);
        GG2=find(ismember(temp,'Target2')==1, 1);
        GG3=find(ismember(temp,'Target2Cleared')==1, 1);
        GG4=find(ismember(temp,'Target2RESP')==1, 1);
        GG5=find(ismember(temp,'Target2ACC')==1, 1);
        GG6=find(ismember(temp,'X1')==1, 1);
        GG7=find(ismember(temp,'X2')==1, 1);
        GG8=find(ismember(temp,'TrialCode')==1, 1);
        GG9=find(ismember(temp,'TRIALID')==1, 1);
        GG10=find(ismember(temp,'BLOCKID')==1, 1);
        GG11=find(ismember(temp,'IStartSaccadeEvent')==1, 1);
        if ~isempty(GG1)
            Targ1ON(trial_index,:) = [str2double(temp{2}) str2double(temp{3}) str2double(temp{7})];
        elseif ~isempty(GG2)
            Targ2ON(trial_index,:) = [str2double(temp{2}) str2double(temp{3}) str2double(temp{6})];
        elseif ~isempty(GG3)
            Targ2OFF(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG4)
            Targ2RESP(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG5)
            Targ2ACC(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG6)
            Targ1POS(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG7)
            Targ2POS(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG8)
            jump_direction{trial_index} = temp{4}(4);
            jump_size(trial_index,:) = str2double(temp{4}(5:end));
        elseif ~isempty(GG9)
            trial_index = trial_index+1;
            trial_in_block(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG10)
            block_id(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        elseif ~isempty(GG11)
            saccade_event(trial_index,:) = [str2double(temp{2}) str2double(temp{4})];
        end
    end
    
    for row_n=[1:size(block_id,1)]
        cur_block = block_id(row_n,2);
        if(cur_block>6)
            %Figuring out direction
            X1pos = Targ1POS(row_n,2);
            X2pos = Targ2POS(row_n,2);
            if(X1pos<1000) %saccade to left
                if(X2pos>X1pos)
                    jump_direction{row_n} = 'B';
                else
                    jump_direction{row_n} = 'F';
                end
            else %saccade to right
                if(X2pos>X1pos)
                    jump_direction{row_n} = 'F';
                else
                    jump_direction{row_n} = 'B';
                end
            end
            %Figuring out target shift - need to compute error here
            jump_size(row_n) = rad2deg(2*(atan(((X2pos-x0)./2)/(dist*scn_res_width/scb_dim_width))))-rad2deg(2*(atan(((X1pos-x0)./2)/(dist*scn_res_width/scb_dim_width))));
            jump_size(row_n) = abs(round(jump_size(row_n) * 4)/4);
        end
    end

    invalid_index = find(Targ1ON(:,1)==0);
    if ~isempty(invalid_index)
        Targ1ON(invalid_index,:) = [];
        Targ2ON(invalid_index,:) = [];
        Targ2OFF(invalid_index,:) = [];
        Targ2RESP(invalid_index,:) = [];
        Targ2ACC(invalid_index,:) = [];
        Targ1POS(invalid_index,:) = [];
        Targ2POS(invalid_index,:) = [];
        jump_direction(invalid_index,:) = [];
        jump_size(invalid_index,:) = [];
        trial_in_block(invalid_index,:) = [];
        block_id(invalid_index,:) = [];
        saccade_event(invalid_index,:) = [];
        if trial_in_block(end,2) == 1
            trial_in_block(end,:) = [];
        end
    else
        Targ1ON(end,:) = [];
        Targ2ON(end,:) = [];
        Targ2OFF(end,:) = [];
        Targ2RESP(end,:) = [];
        Targ2ACC(end,:) = [];
        Targ1POS(end,:) = [];
        Targ2POS(end,:) = [];
        jump_direction(end,:) = [];
        jump_size(end,:) = [];
        trial_in_block(end-1:end,:) = [];
        block_id(end,:) = [];
        saccade_event(end,:) = [];
    end


    % Separate the raw data per trial number
    if_detect_saccade = zeros(length(UniqueTrial),1);
    targ1_onset_time = zeros(length(UniqueTrial),1);
    targ2_onset_time = zeros(length(UniqueTrial),1);
    targ2_offset_time = zeros(length(UniqueTrial),1);
    targ2_response_time = zeros(length(UniqueTrial),1);
    targ2_acc_time = zeros(length(UniqueTrial),1);
    targ1_pos_pixel = zeros(length(UniqueTrial),1);
    targ2_pos_pixel = zeros(length(UniqueTrial),1);
    targ1_pos_degree = zeros(length(UniqueTrial),1);
    targ2_pos_degree = zeros(length(UniqueTrial),1);
    visual_error_pixel_raw = zeros(length(UniqueTrial),1);
    visual_error_pixel_filtered = zeros(length(UniqueTrial),1);
    visual_error_degree_raw = zeros(length(UniqueTrial),1);
    visual_error_degree_filtered = zeros(length(UniqueTrial),1);
    visual_error_endsac_degree_filtered = zeros(length(UniqueTrial),1);
    gaze_degree_pos = zeros(length(UniqueTrial),2);
    exp_saccade2_degree = zeros(length(UniqueTrial),1);
    targ2_jump_direction = cell(length(UniqueTrial),1);
    targ2_jump_size = zeros(length(UniqueTrial),1);
    is_saccade_event = zeros(length(UniqueTrial),1);
    is_response_accurate = zeros(length(UniqueTrial),1);
    manual_response = zeros(length(UniqueTrial),1);

    primary_sac_onset = zeros(length(UniqueTrial),1);
    primary_sac_offset = zeros(length(UniqueTrial),1);
    primary_sac_amp = zeros(length(UniqueTrial),1);
    primary_sac_peak_vel = zeros(length(UniqueTrial),1);
    primary_sac_peak_acc = zeros(length(UniqueTrial),1);

    %% Create structure including all eye movement data under corresponding trial number
    disp('Creating Matrices for the eye movement data');
    for i=1:length(UniqueTrial)
        RawDatapertrial = RawData(RawData(:,1)==UniqueTrial(i),:);
        Saccadepertrial = ESACC(ESACC(:,1)==UniqueTrial(i),:);
        fixation_pertrial = EFIX(find(EFIX(:,1)==UniqueTrial(i)),:);
        blink_pertrial = EBLINK(EBLINK(:,1)==UniqueTrial(i),:);
        block_no = strcat('Block_', num2str(block_id(i,2)));
        trial_in_block_no = strcat('Trial_', num2str(trial_in_block(i,2)));
        trial_start_time = min(RawDatapertrial(:,2));
        % Removal of samples between 0.10 seconds before and 0.10 seconds after a blink artefact
        if ~isempty(blink_pertrial)
            blink_start_perTrial = blink_pertrial(:,3)-trial_start_time;
            blink_end_perTrial = blink_pertrial(:,4)-trial_start_time;
            blink_period = [blink_start_perTrial,blink_end_perTrial];
        else
            blink_start_perTrial = [];
            blink_end_perTrial = [];
            blink_period = [];
        end
        if ~isempty(fixation_pertrial)
            fixation_start_perTrial = fixation_pertrial(:,3)-trial_start_time;
            fixation_end_perTrial = fixation_pertrial(:,4)-trial_start_time;
            fixation_period = [fixation_start_perTrial,fixation_end_perTrial];
            %             fixation_end_pos_pixel = fixation_pertrial(:,6:7);
            %             fixation_end_pos_degree = [rad2deg(2*((atan(((fixation_end_pos_pixel(:,1)-x0)./2)/(dist*scn_res_width/scb_dim_width))))),rad2deg(2*((atan(((fixation_end_pos_pixel(:,2)-x0)./2)/(dist*scn_res_height/scb_dim_height)))))];
            %             fixation_start_pos_pixel=zeros(length(fixation_pertrial(:,3)),2);
            %             for kk = 1:length(fixation_pertrial(:,3))
            %                 fixation_start_pos_pixel(kk,:) = RawDatapertrial(RawDatapertrial(:,2)==fixation_pertrial(i,3),3:4);
            %             end
            %             fixation_start_pos_degree = [rad2deg(2*((atan(((fixation_start_pos_pixel(:,1)-x0)./2)/(dist*scn_res_width/scb_dim_width))))),rad2deg(2*((atan(((fixation_start_pos_pixel(:,2)-x0)./2)/(dist*scn_res_height/scb_dim_height)))))];
        else
            %             fixationstart_perTrial = [];
            %             fixation_end_perTrial = [];
            fixation_period = [];
            %             fixation_start_pos_pixel = [];
            %             fixation_start_pos_degree = [];
            %             fixation_end_pos_pixel = [];
            %             fixation_end_pos_degree = [];
        end

        data.(block_no).(trial_in_block_no).blink_start=blink_start_perTrial;
        data.(block_no).(trial_in_block_no).blink_end=blink_end_perTrial;
        data.(block_no).(trial_in_block_no).fixation_period = fixation_period;
        raw_gaze_pixel = RawDatapertrial(:,3:4);
        raw_gaze_degree = [rad2deg(2*((atan(((raw_gaze_pixel(:,1)-x0)./2)/(dist*scn_res_width/scb_dim_width))))),rad2deg(2*((atan(((raw_gaze_pixel(:,2)-x0)./2)/(dist*scn_res_height/scb_dim_height)))))];

        % Change timestamps - every block starts at 0ms
        time_perTrial = RawDatapertrial(:,2)-trial_start_time;
        data.(block_no).(trial_in_block_no).time_stamps = time_perTrial;
        % Change timestamps for saccade - every block starts at 0ms
        % these saccades were detected by the eye tracker
        saccade_start_perTrial = Saccadepertrial(:,3)-trial_start_time;
        data.(block_no).(trial_in_block_no).saccade_start=saccade_start_perTrial;
        saccade_end_perTrial = Saccadepertrial(:,4)-trial_start_time;
        data.(block_no).(trial_in_block_no).saccade_end=saccade_end_perTrial;
        sac_amplitude = Saccadepertrial(:,10);
        sac_vel = Saccadepertrial(:, 11);

        targ1_onset_time(i) = Targ1ON(i,1) - trial_start_time;
        targ2_onset_time(i) = Targ2ON(i,1) - trial_start_time;
        targ2_offset_time(i) = Targ2OFF(i,1) - trial_start_time;
        targ2_response_time(i) = Targ2RESP(i,1) - trial_start_time;
        targ2_acc_time(i) = Targ2ACC(i,1) - trial_start_time;
        manual_response(i) = Targ2RESP(i,2);
        is_response_accurate(i) = Targ2ACC(i,2);
        targ1_pos_pixel(i) = Targ1POS(i,2);
        targ2_pos_pixel(i) = Targ2POS(i,2);
        targ1_pos_degree(i) = rad2deg(2*((atan(((targ1_pos_pixel(i)-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
        targ2_pos_degree(i) = rad2deg(2*((atan(((targ2_pos_pixel(i)-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
        targ2_jump_direction(i) = jump_direction(i);
        targ2_jump_size(i) = jump_size(i);
        is_saccade_event(i) = saccade_event(i,1)~=0;

        data.(block_no).(trial_in_block_no).targ1_onset_time = targ1_onset_time(i);
        data.(block_no).(trial_in_block_no).targ2_onset_time = targ2_onset_time(i);
        data.(block_no).(trial_in_block_no).targ2_offset_time = targ2_offset_time(i);
        data.(block_no).(trial_in_block_no).targ2_response_time = targ2_response_time(i);
        data.(block_no).(trial_in_block_no).targ2_acc_time = targ2_acc_time(i);
        data.(block_no).(trial_in_block_no).manual_response = manual_response(i);
        data.(block_no).(trial_in_block_no).is_response_accurate = is_response_accurate(i);
        data.(block_no).(trial_in_block_no).targ2_jump_direction = targ2_jump_direction(i);
        data.(block_no).(trial_in_block_no).targ2_jump_size = targ2_jump_size(i);
        data.(block_no).(trial_in_block_no).targ1_pos_pixel = targ1_pos_pixel(i);
        data.(block_no).(trial_in_block_no).targ2_pos_pixel = targ2_pos_pixel(i);
        data.(block_no).(trial_in_block_no).targ1_pos_degree = targ1_pos_degree(i);
        data.(block_no).(trial_in_block_no).targ2_pos_degree = targ2_pos_degree(i);

        % filter out bad saccades and find primary saccade
        is_saccade_on_time = ones(size(saccade_start_perTrial));
        for s = 1:length(saccade_start_perTrial)
            sac_start = saccade_start_perTrial(s)+1;
            if (sac_start < targ1_onset_time(i) || sac_start > targ1_onset_time(i) + time_threshold_after_targ1_onset || sac_amplitude(s) == 0 || sac_vel(s) == 0)
                is_saccade_on_time(s) = 0;
            end
        end
        % saccades detection based on the gaze position
        [sac_onset,sac_offset,sac_amp,sac_peak_vel,sac_peak_acc,gaze_vel,gaze_acc,gaze_degree,gaze_pixel] = Saccade_Detection(raw_gaze_degree, raw_gaze_pixel,blink_period,sampling_freq, cut_off_freq, targ1_pos_degree(i));
        % find primary saccades in the range between target 1 onset and
        %  300 ms after target 1 onset
        if isempty(sac_onset)
            primary_sac_onset(i) = nan;
            primary_sac_offset(i) = nan;
            primary_sac_amp(i) = nan;
            primary_sac_peak_vel(i) = nan;
            primary_sac_peak_acc(i) = nan;
        else
            is_primary_sac = ones(size(sac_onset));
            for s = 1:length(sac_onset)
                sac_start = sac_onset(s)+1;
                if (sac_start < targ1_onset_time(i) || sac_start > targ1_onset_time(i) + time_threshold_after_targ1_onset)
                    is_primary_sac(s) = 0;
                end
            end
            if isempty(nonzeros(is_primary_sac))
                primary_sac_onset(i) = nan;
                primary_sac_offset(i) = nan;
                primary_sac_amp(i) = nan;
                primary_sac_peak_vel(i) = nan;
                primary_sac_peak_acc(i) = nan;
            else
                amp_temp = zeros(size(sac_onset));
                amp_temp(is_primary_sac==1) = sac_amp(is_primary_sac==1);
                [~,index] = max(amp_temp);
                is_primary_sac = zeros(size(sac_onset));
                is_primary_sac(index) = 1;
                primary_sac_onset(i) = sac_onset(is_primary_sac==1);
                primary_sac_offset(i) = sac_offset(is_primary_sac==1);
                primary_sac_amp(i) = sac_amp(is_primary_sac==1);
                primary_sac_peak_vel(i) = sac_peak_vel(is_primary_sac==1);
                primary_sac_peak_acc(i) = sac_peak_acc(is_primary_sac==1);
            end
        end
        % save to data structure
        data.(block_no).(trial_in_block_no).raw_gaze_pixel = raw_gaze_pixel;
        data.(block_no).(trial_in_block_no).raw_gaze_degree = raw_gaze_degree;
        data.(block_no).(trial_in_block_no).gaze_pixel_filtered = gaze_pixel;
        data.(block_no).(trial_in_block_no).gaze_degree_filtered = gaze_degree;
        data.(block_no).(trial_in_block_no).gaze_vel_filtered = gaze_vel;
        data.(block_no).(trial_in_block_no).gaze_acc_filtered = gaze_acc;
        % calculate the visual error from raw gaze data and filtered gaze
        % data
        if targ2_onset_time(i)< length(raw_gaze_pixel(:,1))
            visual_error_pixel_raw(i) = targ2_pos_pixel(i) - raw_gaze_pixel(targ2_onset_time(i)+1,1);
            visual_error_pixel_filtered(i) = targ2_pos_pixel(i) - gaze_pixel(targ2_onset_time(i)+1,1);
            visual_error_degree_raw(i) = targ2_pos_degree(i) - raw_gaze_degree(targ2_onset_time(i)+1,1);
            visual_error_degree_filtered(i) = targ2_pos_degree(i) - gaze_degree(targ2_onset_time(i)+1,1);
            if(~isnan(primary_sac_offset(i)))
                visual_error_endsac_degree_filtered(i) = targ2_pos_degree(i) - gaze_degree(primary_sac_offset(i)+1,1);
            else
                visual_error_endsac_degree_filtered(i) = nan;
            end
            gaze_degree_pos(i,:) = gaze_degree(targ2_onset_time(i)+1,:);
            gaze_degree_pos(i,2) = rad2deg(2*((atan((((gaze_pixel(targ2_onset_time(i)+1,2)-y0)/2)/(dist*scn_res_width/scb_dim_width))))));
            exp_saccade2_degree(i) = targ2_pos_degree(i) - targ1_pos_degree(i);
        else
            visual_error_pixel_raw(i) = nan;
            visual_error_pixel_filtered(i) = nan;
            visual_error_degree_raw(i) = nan;
            visual_error_degree_filtered(i) = nan;
            visual_error_endsac_degree_filtered(i) = nan;
            gaze_degree_pos(i,:) = [nan,nan];
            exp_saccade2_degree(i) = nan;
        end
        data.(block_no).(trial_in_block_no).primary_saccade_amp = primary_sac_amp(i);
        data.(block_no).(trial_in_block_no).primary_saccade_peak_vel = primary_sac_peak_vel(i);
        data.(block_no).(trial_in_block_no).primary_saccade_peak_acc = primary_sac_peak_acc(i);
        data.(block_no).(trial_in_block_no).primary_saccade_onset = primary_sac_onset(i);
        data.(block_no).(trial_in_block_no).primary_saccade_offset = primary_sac_offset(i);
        data.(block_no).(trial_in_block_no).visual_error_pixel_raw = visual_error_pixel_raw(i);
        data.(block_no).(trial_in_block_no).visual_error_pixel_filtered = visual_error_pixel_filtered(i);
        data.(block_no).(trial_in_block_no).visual_error_degree_raw = visual_error_degree_raw(i);
        data.(block_no).(trial_in_block_no).visual_error_degree_filtered = visual_error_degree_filtered(i);
        data.(block_no).(trial_in_block_no).visual_error_endsac_degree_filtered = visual_error_endsac_degree_filtered(i);
        data.(block_no).(trial_in_block_no).gaze_degree_pos = gaze_degree_pos(i,:);
        data.(block_no).(trial_in_block_no).exp_saccade2_degree = exp_saccade2_degree(i);

        if ~isnan(primary_sac_onset(i))
            if_detect_saccade(i) = 1;
        end
        data.(block_no).(trial_in_block_no).if_detect_saccade = if_detect_saccade(i);
    end

    %% find corrected visual error

    fieldstemp = fieldnames(data);
    fieldcount = 1;
    for i = 1:length(fieldstemp)
        if strcmp(fieldstemp{i},'Block_0')
            continue
        else
            fieldsdata{fieldcount,1} = fieldstemp{i};
            fieldcount = fieldcount+1;
        end
    end

    %create array to use for naming fields with 52 trials
    %soft code length of blocks bc DDM uses dif length 9/12/22
    for k = 1:length(fieldnames(data.(fieldsdata{1})))
        nums = num2str(k);
        tempfields{k} = ['Trial_' nums];
    end
    
    %safety net for data that cuts off middway through block SB 7/14/22
    for n = 1:length(fieldsdata)
        for k = 1:length(fieldnames(data.(fieldsdata{n}))) %fill the structure w/ temp field names
            blockidx(n).(tempfields{k}) = data.(fieldsdata{n}).(tempfields{k});
        end
        
        if k ~= length(fieldnames(data.(fieldsdata{1}))) %if k ~= 52 then that block doesn't have 52 trials 
            %length of first block 
            for p = (k+1):length(fieldnames(data.(fieldsdata{1}))) %for the remaining trials that weren't collected - fill tempfields w/ empty structure
                blockidx(n).(tempfields{p}) = struct(tempfields{p},{});
            end
        end
    end
    
    % find fixation error at center
    fixation_error = zeros(1,length(fieldsdata)); % fixation error for each block
    fixation_error_degree = zeros(1,length(fieldsdata));
    for i = 1:length(blockidx)
        fe = [];
        fe_degree = [];
        for k = 1:length(fieldnames(blockidx(i)))
            fieldsblock = fieldnames(blockidx);
            if isempty(blockidx(i).(fieldsblock{k}))
                continue
            else
          
            fixation_period = blockidx(i).(fieldsblock{k}).fixation_period; % fixation detected by the eyetracker
            if_saccade = blockidx(i).(fieldsblock{k}).if_detect_saccade; % if primary saccade was detected
            targ1_onset = blockidx(i).(fieldsblock{k}).targ1_onset_time;
            if (size(fixation_period,1)>=2)% at least 2 fixation detected
                temp_ind = find(fixation_period(:,1)<targ1_onset,1); % find fixation right before targ1 onset
                if ~isempty(temp_ind)
                    fixation_center_start = fixation_period(temp_ind,1); 
                    fixation_center_end = min(fixation_period(temp_ind,2),targ1_onset);
                    if if_saccade == 1
                        fe = [fe,mean(blockidx(i).(fieldsblock{k}).gaze_pixel_filtered(fixation_center_start:fixation_center_end)-960)]; % average the fixation location to compute the fixation error
                        fe_degree = [fe_degree,mean(blockidx(i).(fieldsblock{k}).gaze_degree_filtered(fixation_center_start:fixation_center_end)-960)];
                    end
                end
            end
            end
        end
        fe = fe(~isoutlier(fe));
        fe_degree = fe_degree(~isoutlier(fe_degree));
        fixation_error(i) = mean(fe,'omitnan');
        fixation_error_degree(i) = mean(fe_degree,'omitnan');
    end

    % compute offset gain
    offset_gain_left = zeros(1,length(fieldsdata)); % offset gain for left target for each block
    offset_gain_right = zeros(1,length(fieldsdata));% offset gain for right target for each block
    offset_gain_degree_left = zeros(1,length(fieldsdata));
    offset_gain_degree_right = zeros(1,length(fieldsdata));
    if dist == 1000
        target_pos_left = 439;
        target_pos_right = 1481;
        target_pos_degree_left = rad2deg(2*((atan(((target_pos_left-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
        target_pos_degree_right = rad2deg(2*((atan(((target_pos_right-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
    elseif dist == 1450
        target_pos_left = 584;
        target_pos_right = 1336;
        target_pos_degree_left = rad2deg(2*((atan(((target_pos_left-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
        target_pos_degree_right = rad2deg(2*((atan(((target_pos_right-x0)/2)/(dist*scn_res_width/scb_dim_width)))));
    end
    for i = 1:length(blockidx)
        og_left = [];
        og_right = [];
        og_degree_left = [];
        og_degree_right = [];
        for k = 1:length(fieldnames(blockidx(i)))
              fieldsblock = fieldnames(blockidx);
            if isempty(blockidx(i).(fieldsblock{k}))
                continue
            else
          
            targ1_pos = blockidx(i).(fieldsblock{k}).targ1_pos_pixel;
            jump_size = blockidx(i).(fieldsblock{k}).targ2_jump_size;
            fixation_period = blockidx(i).(fieldsblock{k}).fixation_period;
            targ2_onset = blockidx(i).(fieldsblock{k}).targ2_onset_time;
            targ2_offset = blockidx(i).(fieldsblock{k}).targ2_offset_time;
            if_saccade = blockidx(i).(fieldsblock{k}).if_detect_saccade;
            if (size(fixation_period,1)>=2) % at least 2 fixation detected
                fixation_num_during_targ2_onset = length(find(fixation_period(:,2)>targ2_onset & fixation_period(:,2)<targ2_offset)); % number of fixation detected during targ2 on the screen
                if if_saccade == 1 && jump_size ==0 && (targ1_pos == target_pos_left || targ1_pos == target_pos_right) && fixation_num_during_targ2_onset<=2 % discard trials if there are more than 2 fixation detected during tar2 on the screen
                    temp_ind = find(fixation_period(:,1)>targ2_onset,1); % find fixation after targ2 onset
                    if ~isempty(temp_ind)
                        fixation_targ2_start = fixation_period(temp_ind,1);
                        fixation_targ2_end = min(fixation_period(temp_ind,2),targ2_offset);
                        if targ1_pos == target_pos_right % right target
                            og_right = [og_right,mean(blockidx(i).(fieldsblock{k}).gaze_pixel_filtered(fixation_targ2_start:fixation_targ2_end),'omitnan')-fixation_error(i)-target_pos_right];
                            blockidx(i).(fieldsblock{k}).gaze_pixel_filtered(fixation_targ2_start:fixation_targ2_end)
                            og_degree_right = [og_degree_right,mean(blockidx(i).(fieldsblock{k}).gaze_degree_filtered(fixation_targ2_start:fixation_targ2_end),'omitnan')-fixation_error_degree(i)-target_pos_degree_right];
                            blockidx(i).(fieldsblock{k}).gaze_degree_filtered(fixation_targ2_start:fixation_targ2_end)
                        elseif targ1_pos == target_pos_left % left target
                            og_left = [og_left,mean(blockidx(i).(fieldsblock{k}).gaze_pixel_filtered(fixation_targ2_start:fixation_targ2_end),'omitnan')-fixation_error(i)-target_pos_left];
                            blockidx(i).(fieldsblock{k}).gaze_pixel_filtered(fixation_targ2_start:fixation_targ2_end)
                            og_degree_left = [og_degree_left,mean(blockidx(i).(fieldsblock{k}).gaze_degree_filtered(fixation_targ2_start:fixation_targ2_end),'omitnan')-fixation_error_degree(i)-target_pos_degree_left];
                            blockidx(i).(fieldsblock{k}).gaze_degree_filtered(fixation_targ2_start:fixation_targ2_end)
                        end
                    end
                end
            end
            end
        end
        offset_gain_left(i) = mean(og_left,'omitnan')/(960-target_pos_left);  % offset gain for left target
        offset_gain_right(i) = mean(og_right,'omitnan')/(target_pos_right-960); % offset gain for right target
        offset_gain_degree_left(i) = mean(og_degree_left,'omitnan')/(960-target_pos_degree_left);
        offset_gain_degree_right(i) = mean(og_degree_right,'omitnan')/(target_pos_degree_right-960); 
    end
    offset_gain_left(isnan(offset_gain_left))=0;
    offset_gain_right(isnan(offset_gain_right))=0;
    offset_gain_degree_left(isnan(offset_gain_degree_left))=0;
    offset_gain_degree_right(isnan(offset_gain_degree_right))=0;

    % correct the visual error
    for i = 1:length(blockidx)

        for k = 1:length(fieldnames(blockidx(i)))
            fieldsblock = fieldnames(blockidx);
            if isempty(blockidx(i).(fieldsblock{k}))
                continue
            else
            
            targ2_pos = blockidx(i).(fieldsblock{k}).targ2_pos_pixel;
            targ2_pos_degree = blockidx(i).(fieldsblock{k}).targ2_pos_degree;
            targ2_onset = blockidx(i).(fieldsblock{k}).targ2_onset_time;
            if_saccade = blockidx(i).(fieldsblock{k}).if_detect_saccade;
            gaze_raw = blockidx(i).(fieldsblock{k}).raw_gaze_pixel;
            gaze_filtered = blockidx(i).(fieldsblock{k}).gaze_pixel_filtered;
            gaze_filtered_degree = blockidx(i).(fieldsblock{k}).gaze_degree_filtered;
            if targ2_onset < length(gaze_raw(:,1))
                if targ2_pos < 960
                    visual_error_corrected = targ2_pos - (gaze_filtered(targ2_onset+1,1)-fixation_error(i)-offset_gain_left(i)*(960-targ2_pos)); % correct the visual error
                    visual_error_degree_corrected = targ2_pos_degree - (gaze_filtered_degree(targ2_onset+1,1)-fixation_error_degree(i)-offset_gain_degree_left(i)*(960-targ2_pos_degree));
                elseif targ2_pos > 960
                    visual_error_corrected = targ2_pos - (gaze_filtered(targ2_onset+1,1)-fixation_error(i)-offset_gain_right(i)*(targ2_pos-960));
                    visual_error_degree_corrected = targ2_pos_degree - (gaze_filtered_degree(targ2_onset+1,1)-fixation_error_degree(i)-offset_gain_degree_right(i)*(targ2_pos_degree-960));
                end
            else
                visual_error_corrected = nan;
                visual_error_degree_corrected = nan;
            end
                
            trial_in_block_no = strcat('Trial_', num2str(k));
            data.(fieldsdata{i,1}).(trial_in_block_no).visual_error_corrected = visual_error_corrected;
            data.(fieldsdata{i,1}).(trial_in_block_no).visual_error_degree_corrected = visual_error_degree_corrected;
            end
           
        end
    end

    % save data for all trials to info
    info.fixation_error = fixation_error;
    info.offset_gain_left = offset_gain_left;
    info.offset_gain_right = offset_gain_right;
    info.primary_saccade_amp = primary_sac_amp;
    info.primary_sac_peak_vel = primary_sac_peak_vel;
    info.primary_sac_peak_acc = primary_sac_peak_acc;
    info.primary_sac_onset = primary_sac_onset;
    info.primary_sac_offset = primary_sac_offset;
    info.visual_error_pixel_raw = visual_error_pixel_raw;
    info.visual_error_pixel_filtered = visual_error_pixel_filtered;
    info.visual_error_pixel_corrected = visual_error_corrected;
    info.visual_error_degree_raw = visual_error_degree_raw;
    info.visual_error_degree_filtered = visual_error_degree_filtered;
    info.visual_error_endsac_degree_filtered = visual_error_endsac_degree_filtered;
    info.gaze_degree_pos = gaze_degree_pos;
    info.visual_error_degree_corrected = visual_error_degree_corrected;
    info.exp_saccade2_degree = exp_saccade2_degree;
    info.targ1_onset_time = targ1_onset_time;
    info.targ2_onset_time = targ2_onset_time;
    info.targ2_offset_time = targ2_offset_time;
    info.targ2_response_time = targ2_response_time;
    info.targ2_acc_time = targ2_acc_time;
    info.targ1_pos_pixel = targ1_pos_pixel;
    info.targ2_pos_pixel = targ2_pos_pixel;
    info.targ1_pos_degree = targ1_pos_degree;
    info.targ2_pos_degree = targ2_pos_degree;
    info.targ2_jump_direction = targ2_jump_direction;
    info.targ2_jump_size = targ2_jump_size;
    info.is_saccade_event = is_saccade_event;
    info.is_response_accurate = is_response_accurate;
    info.if_detect_saccade = if_detect_saccade;
    info.manual_response = manual_response;

    if Save==1
        disp('Saving the data'); %per session
        start_of_file_name = 'Data_';
        File_name = strcat(start_of_file_name,FileName(1:end-4));
        save([path2crunch File_name],'data','info');
        disp(['...done. (' FileName ')']);
    end

end

