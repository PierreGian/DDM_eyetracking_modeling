function [bias, threshold] = Psych_Threshold_Bias(result)

 x = linspace(min(result.data(:,1)),max(result.data(:,1)),1000);
 fitValues = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);
 if max(fitValues) >= 0.75 && min(fitValues) <= 0.5
    threshold = x(find((fitValues>=0.75),1)) - x(find((fitValues>=0.5),1));
    bias = x(find((fitValues>=0.5),1));
else
    threshold = x(end) - x(find((fitValues>=0.5),1));
    bias = x(find((fitValues>=0.5),1));
end