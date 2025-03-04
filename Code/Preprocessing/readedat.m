function [subject] = readedat(str)
% [subject] = readedat(str)
% str = type string - the name of the subject edat file to be read
% subject = type cell array - containing the data from the file
% Function will read in a ASCII text exported EDAT file
% Or a edat exported by my cygwin script
% Dennis Thompson, IRC - UCDavis 2/2016
% DMT UCDavis 9/2017 Rewrote from scratch. Cleaned up the logic 
%str = 'Hand_good.txt';
%str='Powershell_bad.csv';
%str = '/nfs/agency/raw_data/behave_suite/edats/NP101/TSDT_B_NP101_0.csv';
[fid,message] = fopen(str,'r');
if fid == -1
  error('Author:Function:OpenFile', 'Cannot open file: %s', str);
end
% Header has 3 posiable states
% Path, 3 lines of Junk, Header (output of using E-Prime in command Line mode)
% 3 lines of Junk, Header  (output of E-Prime GUI)
% Header (File has been cleaned up in some way)
% Assume delimiter is tab or comma


% read in first 5 lines
% first using commas then reread using tab delimiters
for n = 1:5
    comma = strsplit(fgetl(fid),',');
end
fseek(fid,0,'bof');
for n = 1:5
    tab = strsplit(fgetl(fid),'\t');
end
fseek(fid,0,'bof');

% the larger number of cells should determine if the file is tab or comma
% delimited
if size(comma,2) > size(tab,2) 
    delim = ','; 
else
    delim = '\t';
end

%return file to the beginning
fseek(fid,0,'bof');
% now look for the header (Assume first header is 'ExperimentName')
while ~feof(fid)
    rawline = strsplit(fgetl(fid),delim);
    if strcmp(rawline{1},'ExperimentName'),
        header = rawline;
        dataposition = ftell(fid);
        clear rawline;
        break
    end
end

% read all data in as text
n = 1;
while ~feof(fid)
    rawtext{n} = regexp(fgetl(fid),delim,'split');
    n = n + 1;
end
fclose(fid);
clear fid;
% and then convert numeric strings to doubles
% also transpose the array from n by k to k by n
for n = 1:size(rawtext,2),
    for k = 1:size(rawtext{n},2),
        % check for empty cells
        if isempty(rawtext{n}{k}), rawtext{n}{k} = NaN; end;
        % look for NULL
        if strcmp(rawtext{n}{k},'NULL'), rawtext{n}{k} = NaN; end;
         % look for string NaN
        if strcmp(rawtext{n}{k},'NaN'), rawtext{n}{k} = NaN; end;
        % look for string NAN
        if strcmp(rawtext{n}{k},'NAN'), rawtext{n}{k} = NaN; end;
        % look for numbers
        temp = str2double(rawtext{n}{k});       
        if isnan(temp)
            data{k}{n} = rawtext{n}{k}; % value is string
        else
            data{k}{n} = temp; % value is number
        end
    end
end

% remove [] () and '.' from headers
% place data in to fields
for n = 1:size(header,2),
    subject{n}.header = header{n}(isstrprop(header{n},'alphanum'));
    subject{n}.col = data{n};
end

% check col and make all data match (char or double
for n = 1:size(subject,2),
    test = cellfun(@isnumeric,subject{n}.col);
    % if all numbers - convert to array
    if sum(test)==length(test),
       col(test) = [subject{n}.col{test}];
       subject{n}.col = col;
    elseif sum(test) ~= 0 % is mixed number and strings - convert all numbers to string
        for k = find(test)
            subject{n}.col{k} = num2str(subject{n}.col{k});
        end
    end
end











