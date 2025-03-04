function [edat,all_data,all_info] = Process_split_asc(subj,dist2screen,append_length,files,path_csv,path_edf,path_crunch)
       % Pierre Gianferrara, for NCAP project in the sensorimotor lab - UCDavis 11/2023
       % The purpose of this script is to preprocess *.asc split data into
       % separate files
       % Outputs: 1) edat after running Dennis Thompson's "readedat()" script
       % 2) all_data with all data struct objects in one place
       % 3) all_info with all info struct objects in one place
       path_crunch = [path_crunch,'split_data/'];
       all_data = struct();
       all_info = struct();
       edat = readedat([path_csv subj '/' files(1).name]); % Call the read edat file function, import all original data
       list_files = {[path_crunch 'Data_' subj '_trimmed.mat']};
       list_names = {[subj,'_trimmed.asc']};
       if(append_length==1)
           list_files{end+1} = [path_crunch 'Data_' subj '_appended.mat'];
           list_names{end+1} = [subj,'_appended.asc'];
       else
           for k=1:append_length
              list_files{end+1} = [path_crunch 'Data_' subj '_appended_',int2str(k),'.mat'];
              list_names{end+1} = [subj,'_appended_',int2str(k),'.asc'];
           end
       end
       if ~isfile([path_crunch 'Data_' subj '_trimmed.mat']) %crunching has not happened yet
           for ff = 1:length(list_files)
               data_nb = ['data',int2str(ff)];
               info_nb = ['info',int2str(ff)];
               f_name = list_names{ff};
               disp(f_name)
               [all_data.(data_nb),all_info.(info_nb)] = New_Eye_Tracking_Crunch([path_edf subj '/'],dist2screen,path_crunch,f_name);
           end
           disp('Saving the concatenated data'); %per session
           start_of_file_name = 'Data_';
           File_name = strcat(start_of_file_name,subj,'_ALL.mat');
           save([path_crunch File_name],'all_data','all_info');
           disp(['...done. (' File_name ')']);
       end
end