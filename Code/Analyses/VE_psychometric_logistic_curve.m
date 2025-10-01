%Save % incorrect & % L/R VE similarity across target shifts
clear all

%%%%%%%%%%%% Initialize key parameters for batch script here
analysis_type = 'VE'; %options are 'VECD' and 'VE'
text_extension = '/*.csv'; % asign txt or csv
dist2screen=1000; % 1000 for behave, 1450 for scan
%%%%%%%%%%%%%%%%%%%%%

disp(['Running DDM modeling analysis for ' analysis_type ' case']);

%Path initialization
root_path = '/nfs/agency/';
master_path = [root_path,'code/DDM/DDM_Master_code/'];
addpath(genpath([root_path,'code/DDM/DDM_Master_code/']));
behave_stats_path = [root_path,'raw_data/behave_suite/']; %eye-tracking data from behavioral task
% addpath('/nfs/agency/code/ForwardCurve/Weiwei/')
% addpath('/nfs/agency/code/DDM/')
path_crunch = [master_path,'eyetracking_crunch/behave/'];
path_edf = '/nfs/agency/raw_data/DDM/edf/';
path_csv = '/nfs/agency/raw_data/DDM/edats/';
path_output = [master_path 'outputs/'];
path_plots = [master_path,'plots/DDM_VE/VE_compare2/'];
path_output_VE = [master_path 'outputs/Psychometric/'];
path_RT_plots = [master_path,'plots/DDM/RT/'];
behave_data_folder = {dir(path_edf).name}; %Folder with all subject names

%Outputting the paths corresponding to edf, csv, and crunch files
disp(['path_edf is ' path_edf])
disp(['path_csv is ' path_csv])
disp(['path_crunch is ' path_crunch])

%Vector of jump sizes - input required for modeling analyses
%X=[-2.5 -1.75 -1.25 -1 -0.75 -0.5 -0.25 -0.175 -0.1 0 0.1 0.175 0.25 0.5 0.75 1 1.25 1.75 2.5];
X=[-2.5 -1.75 -1 -0.5 -0.25 0 0.25 0.5 1 1.75 2.5];
Amp_fields = {'amp4','amp8'};
Dir_fields = {'Left','Right','Total'};
% find the index for 0 jump
zero_jump_index = find(X==0);

%Identifying list of subject names
Data_ToLoad = [path_output,'Subject_Data.mat'];
load(Data_ToLoad)
subject_list = fieldnames(Subject_Data)';

%Extract subject names from data struct

%Need to plot full case
N = length(subject_list);
X_n = length(X);

Select_idx = [1:length(subject_list)];

amp_inc = zeros(N,X_n);
amp_ve_all = zeros(N,X_n);

ToLoad = [path_output_VE,'AllChoice_VE.mat'];
load(ToLoad);
amp4_dat = All_choice.Amp4_Total;
amp8_dat = All_choice.Amp8_Total;

%amp_ve_all = ((amp4_dat+amp8_dat)./2).*2-1;
amp_ve_all = (amp4_dat+amp8_dat)./2;

%To plot on x-axis is X
% Group data and labels
groups = {
    nanmean(amp_ve_all,1)', 'Avg', 'b';
};

    SEs = {
        nanstd(amp_ve_all,0,1)./sqrt(N), 'Avg', 'b';
    };

step_fct = @(x) 1/(1+exp(-100000000000000000*(x-0)));

% Fit logistic models and plot
fitResults = cell(size(groups, 1), 1);
figure;
hold on;
for i = 1:size(groups, 1)
    [fitObj, yFit] = fitLogistic(X', groups{i, 1});
    plotFit(X, groups{i, 1}, fitObj, groups{i, 3}, groups{i, 2});
    fitResults{i} = struct('name', groups{i, 2}, 'fit', fitObj, 'color', groups{i, 3});
    disp(fitObj);
end
legend('show');
xlabel('Target Shift');
ylabel('VE response');
title('Logistic Function Fit');
hold off;

step_fit = 0.05;
h0=figure();
ax = gca;
[fitObj1, yFit] = fitLogistic(X', groups{1, 1});
xToFit = [-10:step_fit:10];
xToFit2 = [-10:0.001:10];
step_fs = zeros(1,length(xToFit2));
for xxx = [1:length(xToFit2)]
    step_fs(xxx) = step_fct(xToFit2(xxx));
end
yFit1 = fitObj1(xToFit);
yline(0.5,'--', 'Color','#999999','HandleVisibility', 'off')
hold on
%plot(X, nanmean(amp4_ve_LRdif,1), 'r.', 'HandleVisibility', 'off'); % raw data (not shown in legend)
errorbar(X,nanmean(amp_ve_all,1),SEs{1,1},'.','Color','#7703fc','MarkerSize',10,'MarkerFaceColor','blue', 'LineWidth', 1.1, 'HandleVisibility', 'off')
hold on
plot(xToFit, yFit1, '-','Color','#7703fc','LineWidth', 1.2, 'HandleVisibility', 'on');
hold on
plot(xToFit2, step_fs, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'on')
hold on
legend('Logistic fit','Theoretical switch','FontSize', 12,'Location','northwest');
xlabel('Target Shift (degrees)','FontSize', 16);
ylabel('% forward responses','FontSize', 16);
ylim([-0.25,1.5])
yticks([-1 0 1])
xlim([-3 3])
title('VE psychometric Group Fits','FontSize', 18);
box(ax, 'off');
hold off;
saveas(h0, [path_plots 'Group_VE_logistic_Avg.jpg']);
close(h0)

%toplot_idx = [14,1,8];
toplot_idx = [1:length(subject_list)];
% toplot_idx(17) = [];
% plotted = subject_list;
% plotted{17} = [];
inf_pt = zeros(length(toplot_idx),2);
growth_r = zeros(length(toplot_idx),2);

close all;

for ii = [1:length(toplot_idx)]
    if(ii==17)
        inf_pt(ii,1) = NaN;
        inf_pt(ii,2) = NaN;
        growth_r(ii,1) = NaN;
        growth_r(ii,2) = NaN;
        continue
    end
    c_idx = toplot_idx(ii);
    cur_name = subject_list{c_idx};
    h=figure();
    amp_y = amp_ve_all(c_idx,:)';
    ax = gca;
    ax.FontSize = 14;
    [fitObj1, yFit] = fitLogistic(X', amp_y);
    xToFit = [-10:step_fit:10];
    yFit1 = fitObj1(xToFit);
    inf_pt(ii,1) = fitObj1.x0;
    growth_r(ii,1) = fitObj1.k;
    yline(0.5,'--', 'Color','#999999','HandleVisibility', 'off')
    hold on
    plot(X, amp_y, '.','Color','#7703fc','MarkerSize',10 , 'HandleVisibility', 'off'); % raw data (not shown in legend)
    hold on
    plot(xToFit, yFit1, '-','Color','#7703fc','LineWidth', 1.5, 'HandleVisibility', 'on');
    hold on
    plot(xToFit2, step_fs, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'on')
    hold on
    legend('Logistic fit','Theoretical switch','FontSize', 12,'Location','northwest');

    %legend('Avg','FontSize', 14,'Location','southeast');
    xlabel('Target Shift (degrees)','FontSize', 16);
    ylabel('% forward responses','FontSize', 16);
    ylim([-0.25,1.5])
    yticks([-1 0 1])
    xlim([-3 3])
    c_title=['Individual VE Psychometric Fits'];
    title(c_title,'FontSize', 18);
    box(ax, 'off');
    hold off;
    saveas(h, [path_plots cur_name '_VE_logistic_Avg.jpg']);
    close(h)
end

% % Normalization plot
% xFine = linspace(-10, 0.01, 10);
% figure;
% subplot(2,1,1);
% hold on;
% for i = 1:length(fitResults)
%     f = createLogisticFunction(fitResults{i}.fit);
%     yOrig = f(xFine);
%     plot(xFine, yOrig, 'Color', fitResults{i}.color, 'LineWidth', 2);
%     plotAsymptoteLines(fitResults{i});
% end
% title('Original Fitted Logistic Function');
% xlabel('x'); ylabel('f(x)');
% grid on;


%% Helper Functions

% function [fitObj, yFit] = fitLogistic(x, y)
%     ftype = fittype('2/(1+exp(-k*(x-x0))) + b', 'independent', 'x', 'coefficients', {'k','x0','b'});
%     fitObj = fit(x, y, ftype, 'StartPoint', [0.5, 0, 0]);
%     yFit = fitObj(x);
% end

function [fitObj, yFit] = fitLogistic(x, y)
    ftype = fittype('1/(1+exp(-k*(x-x0)))', 'independent', 'x', 'coefficients', {'k','x0'});
    fitObj = fit(x, y, ftype, 'StartPoint', [0.5, 0]);
    yFit = fitObj(x);
end

function f = createLogisticFunction(fitObj)
    L = fitObj.L; k = fitObj.k; x0 = fitObj.x0; b = fitObj.b;
    f = @(x) L ./ (1 + exp(-k * (x - x0))) + b;
end

function g = normalizeLogisticFunction(f, fitObj)
    L = fitObj.L; b = fitObj.b;
    g = @(x) (f(x) - (b + L)) / (-L);  % Normalize to [0,1]
end

function plotFit(x, y, fitObj, color, groupName)
    plot(x, y, 'k-', 'HandleVisibility', 'off'); % raw data (not shown in legend)
    hold on;

    % Plot the fitted curve with DisplayName to show in legend
    hFit = plot(fitObj, color, x, y, [color '.']);
    hFit(1).DisplayName = groupName; % Only set DisplayName on curve
    hFit(2).HandleVisibility = 'off'; % disable marker legend entry
end

function plotAsymptotes(result, x)
    f = createLogisticFunction(result.fit);
    yFit = f(x);
    %%PLot log curve
    plot(x, yFit, '-', 'Color', result.color, 'LineWidth', 2, 'DisplayName', ['Logistic - ' result.name]);
    % Asymptote, inflection points
    yline(result.fit.b, '--', ['Upper Asymptote - ' result.name], 'Color', result.color, 'HandleVisibility', 'off');
    yline(result.fit.b + result.fit.L, '--', ['Lower Asymptote - ' result.name], 'Color', result.color, 'HandleVisibility', 'off');
    xline(result.fit.x0, '--', ['x0 - ' result.name], 'Color', result.color, 'HandleVisibility', 'off');
end

function plotAsymptoteLines(result)
    yline(result.fit.b, '--', ['b - ' result.name], 'Color', result.color);
    yline(result.fit.b + result.fit.L, '--', ['b+L - ' result.name], 'Color', result.color);
    xline(result.fit.x0, '--', ['x0 - ' result.name], 'Color', result.color);
end