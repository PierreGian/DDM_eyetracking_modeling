function [B4,Y4,fcn] = S_curve_calculating_fitting(ForwPercent_amp4,X);
%% Calculate the S fitting curve
disp('Calculating CDF curve')
% X=[-2.5 -1.75 -1.25 -1 -0.75 -0.5 -0.25 -0.175 -0.1 0 0.1 0.175 0.25 0.5 0.75 1 1.25 1.75 2.5];

Y4=ForwPercent_amp4;
% Y4L=ForwPercent_amp4_L;
% Y4R=ForwPercent_amp4_R;
% 
% Y8=ForwPercent_amp8;
% Y8L=ForwPercent_amp8_L;
% Y8R=ForwPercent_amp8_R;

fcn = @(b,x) normcdf(x, b(1), b(2));                    % Objective Function

NRCF4 = @(b) norm(Y4 - fcn(b,X)); 
% NRCF4L = @(b) norm(Y4L - fcn(b,X));                     % Norm Residual Cost Function
% NRCF4R = @(b) norm(Y4R - fcn(b,X)); 

B4 = fminsearch(NRCF4, [0; 10])';  
clear fcn;
clear NRCF4;
% B4L = fminsearch(NRCF4L, [0; 10]);                          % Estimate Parameters
% B4R = fminsearch(NRCF4R, [0; 10]);
% 
% NRCF8 = @(b) norm(Y8 - fcn(b,X)); 
% NRCF8L = @(b) norm(Y8L - fcn(b,X));               
% NRCF8R = @(b) norm(Y8R - fcn(b,X)); 

% B8 = fminsearch(NRCF8, [0; 10]);  
% B8L = fminsearch(NbRCF8L, [0; 10]);                         
% B8R = fminsearch(NRCF8R, [0; 10]);
