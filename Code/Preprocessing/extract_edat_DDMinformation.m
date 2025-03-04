function [EDAT_info] = extract_task_information(edat)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 12/2023
       % The purpose of this script is to load subject data and extract all the task-relevant information
       EDAT_info = struct();
       
        %first identify which cells are which
        for n = 1:length(edat)
            ToPrint = edat{n}.header;
            %disp(num2str(n))
            %disp(ToPrint)
            if(strcmp(edat{n}.header,'BlockNumBlock'))           block = n; end     % strcmp=string comparison
            if(strcmp(edat{n}.header,'TrialNumBlock'))              trial = n; end   % find the colomn number of each useful variable and give them a name
            if(strcmp(edat{n}.header,'Target2ACCBlock'))          t2acc = n; end
            if(strcmp(edat{n}.header,'Target2OnsetTimeBlock'))      t2onset = n; end
            if(strcmp(edat{n}.header,'FixationOffsetTimeBlock'))    offset= n;end %BlankScreen1OffsetTime??
            if(strcmp(edat{n}.header,'X1Block'))    x1_ind= n;end %X1 position
            if(strcmp(edat{n}.header,'X2Block'))    x2_ind= n;end %X2 position
        end
    EDAT_info.Blocks = edat{block}.col;
    EDAT_info.Trials = edat{trial}.col;
    EDAT_info.Target2ACC = edat{t2acc}.col;
    EDAT_info.Target2OnsetTime = edat{t2onset}.col;
    EDAT_info.Offset = edat{offset}.col;
    EDAT_info.X1 = edat{x1_ind}.col;
    EDAT_info.X2 = edat{x2_ind}.col;
end