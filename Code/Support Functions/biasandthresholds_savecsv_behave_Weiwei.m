function [outputArg1,outputArg2] = biasandthresholds_savecsv_behave_Weiwei(X,B4,Y4,B8,Y8,subj,analysis_type,path2output)
%CURVE_BEHAVE_SAVECSV Summary of this function goes here
%   Detailed explanation goes here

outputpath = [path2output 'biasandthresholds_' analysis_type{1} '.csv'];

if ~isfile(outputpath)
    
    biasandthresholds = cell(2,5);
    
    biasandthresholds{1,1} = 'subject ID';
    biasandthresholds{1,2} = 'bias4';
    biasandthresholds{1,3} = 'threshold4';
    biasandthresholds{1,4} = 'bias8';
    biasandthresholds{1,5} = 'threshold8';
    writecell(biasandthresholds, outputpath);
    
end

csvpath = fullfile(outputpath);
biasandthresholds = readcell(csvpath);
heightdata = size(biasandthresholds);

[bias4, threshold4] = S_plot_w_bias_thresh(X,B4,Y4);

[bias8, threshold8] = S_plot_w_bias_thresh(X,B8,Y8);


for r = 1:heightdata(1)
    
    if r == 1
        continue
        
    elseif strcmp(biasandthresholds{r,1},subj) %if subj already in csv don't add again SB 12/13/21
        break
    elseif ~ismissing(biasandthresholds{r,1})
        continue
    else
        biasandthresholds{r,1} = subj;
        biasandthresholds{r,2} = bias4;
        biasandthresholds{r,3} = threshold4;
        biasandthresholds{r,4} = bias8;
        biasandthresholds{r,5} = threshold8;
        for c = 1:heightdata(2)
            biasandthresholds{r+1,c} = {};
        end
    end
end

biasandthresholds = cell2table(biasandthresholds(2:end,:),'VariableNames',biasandthresholds(1,:));

writetable(biasandthresholds, outputpath);

end

