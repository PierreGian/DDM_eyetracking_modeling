function Save_Behave_Bias_Threshold_pg(result,subj,outputpath)
%CURVE_BEHAVE_SAVECSV Summary of this function goes here
%   Detailed explanation goes here

csvpath = fullfile(outputpath);
biasandthresholds = readcell(csvpath);
row_index = size(biasandthresholds,1)+1;

% [bias4, threshold4] = S_plot_w_bias_thresh(X,B4,Y4);
[bias4_right, threshold4_right] = Psych_Threshold_Bias(result.result4R);
[bias4_left, threshold4_left] = Psych_Threshold_Bias(result.result4L);
[bias4, threshold4] = Psych_Threshold_Bias(result.result4);
[bias8_right, threshold8_right] = Psych_Threshold_Bias(result.result8R);
[bias8_left, threshold8_left] = Psych_Threshold_Bias(result.result8L);
[bias8, threshold8] = Psych_Threshold_Bias(result.result8);

if isempty(bias4_right)
    bias4_right = NaN;
end
if isempty(threshold4_right)
    threshold4_right = NaN;
end
if isempty(bias4_left)
    bias4_left = NaN;
end
if isempty(threshold4_left)
    threshold4_left = NaN;
end
if isempty(bias4)
    bias4 = NaN;
end
if isempty(threshold4)
    threshold4 = NaN;
end

if isempty(bias8_right)
    bias8_right = NaN;
end
if isempty(threshold8_right)
    threshold8_right = NaN;
end
if isempty(bias8_left)
    bias8_left = NaN;
end
if isempty(threshold8_left)
    threshold8_left = NaN;
end
if isempty(bias8)
    bias8 = NaN;
end
if isempty(threshold8)
    threshold8 = NaN;
end

biasandthresholds{row_index,1} = subj;
biasandthresholds{row_index,2} = num2str(bias4_right);
biasandthresholds{row_index,3} = num2str(threshold4_right);
biasandthresholds{row_index,4} = num2str(bias4_left);
biasandthresholds{row_index,5} = num2str(threshold4_left);
biasandthresholds{row_index,6} = num2str(bias4);
biasandthresholds{row_index,7} = num2str(threshold4);
biasandthresholds{row_index,8} = num2str(bias8_right);
biasandthresholds{row_index,9} = num2str(threshold8_right);
biasandthresholds{row_index,10} = num2str(bias8_left);
biasandthresholds{row_index,11} = num2str(threshold8_left);
biasandthresholds{row_index,12} = num2str(bias8);
biasandthresholds{row_index,13} = num2str(threshold8);

% biasandthresholds = cell2table(biasandthresholds(2:end,:),'VariableNames',biasandthresholds(1,:));

writecell(biasandthresholds, outputpath);
