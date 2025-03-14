load('output_td2.mat')
% load('output_td3.mat')
% load('output_td4.mat')
load('output_ekf.mat')

load('output_raps_bp50.mat')
load('output_raps_nbp50.mat')

load('output_raps_bcd.mat')
load('output_raps_bcdnb.mat')
load('output_raps_cvxbcd.mat')
% output_ekf = output;
% p_ekf = p;
% save('output_ekf.mat', 'output_ekf', 'p_ekf');
% output_td = output;
% p_td = p;
% save('output_td.mat', 'output_td', 'p_td');
% output_raps = output;
% p_raps = p;
% save('output_raps.mat', 'output_raps', 'p_raps');
%% Compare CVX and BCD
set(0,'defaultfigurecolor','w')
purple = [0.4940, 0.1840, 0.5560];
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
green = [0.4660, 0.6740, 0.1880];
orange = [0.9290, 0.6940, 0.1250];

disp('Compare Global and Local optimum')
% Local
hor_err = output_raps_bcd.hor_err(521:720);
ver_err = abs(output_raps_bcd.ned_err(3,521:720));
[f_bcd_hor,x_bcd_hor] = ecdf(hor_err);
[f_bcd_ver,x_bcd_ver] = ecdf(ver_err);
[fcost_bcd,xcost_bcd] = ecdf(output_raps_bcd.cost(521:720));
fprintf('RAPS Local Hor <= 1.5 m: %.2f%%\n', sum(hor_err <= 1.5) / length(hor_err) * 100);
fprintf('RAPS Local Ver <= 3 m: %.2f%%\n', sum(ver_err <= 3) / length(ver_err) * 100);
% Global
hor_err = output_raps_cvxbcd.hor_err(521:720);
ver_err = abs(output_raps_cvxbcd.ned_err(3,521:720));
[f_cvx_hor,x_cvx_hor] = ecdf(hor_err);
[f_cvx_ver,x_cvx_ver] = ecdf(ver_err);
[fcost_cvx,xcost_cvx] = ecdf(output_raps_cvxbcd.cost(521:720));
fprintf('RAPS Global Hor <= 1.5 m: %.2f%%\n', sum(hor_err <= 1.5) / length(hor_err) * 100);
fprintf('RAPS Global Ver <= 3 m: %.2f%%\n', sum(ver_err <= 3) / length(ver_err) * 100);
% Non-binary
hor_err = output_raps_bcdnb.hor_err(521:720);
ver_err = abs(output_raps_bcdnb.ned_err(3,521:720));
[f_bcdnb_hor,x_bcdnb_hor] = ecdf(hor_err);
[f_bcdnb_ver,x_bcdnb_ver] = ecdf(ver_err);
[fcost_bcdnb,xcost_bcdnb] = ecdf(output_raps_bcdnb.cost(521:720));
fprintf('RAPS Non-binary Hor <= 1.5 m: %.2f%%\n', sum(hor_err <= 1.5) / length(hor_err) * 100);
fprintf('RAPS Non-binary Ver <= 3 m: %.2f%%\n', sum(ver_err <= 3) / length(ver_err) * 100);

figure(1)
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(x_cvx_hor, f_cvx_hor, '-','Color',blue);
hold on
plot(x_bcd_hor, f_bcd_hor, '-','Color',red);
plot(x_bcdnb_hor, f_bcdnb_hor, '-','Color',green);
grid on
legend('Binary Global', 'Binary Local','Non-binary Local')
ylabel('Probability');
xlim([0,2.5])
title('a. Cumulative Probability of HE');
nexttile
plot(x_cvx_ver, f_cvx_ver, '-','Color',blue);
hold on
plot(x_bcd_ver, f_bcd_ver, '-','Color',red);
plot(x_bcdnb_ver, f_bcdnb_ver, '-','Color',green);
grid on
xlim([0,5])
legend('Binary Global', 'Binary Local','Non-binary Local')
xlabel('Error, m');
ylabel('Probability');
title('b. Cumulative Probability of VE');
nexttile
L = length(output_raps_cvxbcd.cost(521:720));
h_cvx = semilogx(xcost_cvx,fcost_cvx); hold on;
h_cvx.Color = blue;
h_cvx.Marker = '>';
h_cvx.LineWidth = 0.2;
h_cvx.MarkerIndices = 1:20:L;

h_bcd = semilogx(xcost_bcd,fcost_bcd); hold on;
h_bcd.Color = red;
h_bcd.Marker = 'o';
h_bcd.LineWidth = 0.2;
h_bcd.MarkerIndices = 1:20:L;

h_bcdnb = semilogx(xcost_bcdnb,fcost_bcdnb); hold on;
h_bcdnb.Color = green;
h_bcdnb.Marker = '*';
h_bcdnb.LineWidth = 0.2;
h_bcdnb.MarkerIndices = 1:20:L;
grid on
legend('Binary Global', 'Binary Local','Non-binary Local')
xlabel('Risk');
ylabel('Probability');
title('c. Cumulative Probability of Risk');
%------------------------------------%
figure(2)
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(output_raps_cvxbcd.cost(521:720),'.','Color',blue)
hold on
plot(output_raps_cvxbcd.cost_bcd(521:720),'.','Color',red)
grid on
legend('Global', 'Local')
ylabel('Risk')
title('a. Optimal objective value')
% set(gca, 'YScale', 'log') 
nexttile
plot(output_raps_cvxbcd.cost_bcd(521:720)-output_raps_cvxbcd.cost(521:720),'.','Color','k')
grid on
ylabel('\Delta Risk')
ylim([0,3])
title('b. Optimization error (C_\mu - C_g)')
nexttile
plot(output_raps_cvxbcd.comp_time(521:720),'.','Color',blue)
hold on
plot(output_raps_cvxbcd.comp_time_bcd(521:720),'.','Color',red)
hold off
legend('Global', 'Local')
set(gca, 'YScale', 'log')
set(gca, 'XGrid', 'on', 'YGrid', 'on')
xlabel('Epoch')
ylabel('Time, seconds')
title('c. Computation time')

% figure(3)
% plot(output_raps_cvxbcd.raps_flag(521:720))
%%
ind = ~isnan(output_raps_bp50.cost);
raps_b_hor_err = output_raps_bp50.hor_err(ind);
raps_b_ver_err = abs(output_raps_bp50.ned_err(3,ind));
raps_b_risk = output_raps_bp50.pos_risk(ind);
raps_b_meas = output_raps_bp50.num_meas_used(ind);
raps_b_penalty = output_raps_bp50.raps_penalty(ind);
raps_b_compt = output_raps_bp50.comp_time(ind);
raps_b_gdop = output_raps_bp50.GDOP(ind);
raps_b_horcov = zeros(1,length(ind));
raps_b_vercov = sqrt(output_raps_bp50.ned_cov(3,ind));
for i=1:length(ind)
    if (output_raps_bp50.ned_cov(1,i) == 400)
        raps_b_vercov(i) = NaN;
    else
        raps_b_horcov(i) = norm(output_raps_bp50.ned_cov(1:2,i));
    end
end
raps_b_horcov = raps_b_horcov(ind);

ind = ~isnan(output_raps_nbp50.cost);
raps_nb_hor_err = output_raps_nbp50.hor_err(ind);
raps_nb_ver_err = abs(output_raps_nbp50.ned_err(3,ind));
raps_nb_risk = output_raps_nbp50.pos_risk(ind);
raps_nb_meas = output_raps_nbp50.num_meas_used(ind);
raps_nb_penalty = output_raps_nbp50.raps_penalty(ind);
raps_nb_compt = output_raps_nbp50.comp_time(ind);
raps_nb_gdop = output_raps_nbp50.GDOP(ind);
raps_nb_horcov = zeros(1,length(ind));
raps_nb_vercov = sqrt(output_raps_nbp50.ned_cov(3,ind));
for i=1:length(ind)
    if (output_raps_nbp50.ned_cov(1,i) == 400)
        raps_nb_vercov(i) = NaN;
    else
        raps_nb_horcov(i) = norm(output_raps_nbp50.ned_cov(1:2,i));
    end
end
raps_nb_horcov = raps_nb_horcov(ind);

ind = ~isnan(output_td2.cost);
td2_hor_err = output_td2.hor_err(ind);
td2_ver_err = abs(output_td2.ned_err(3,ind));
td2_risk = output_td2.cost(ind);
td2_meas = output_td2.num_meas_used(ind);
td2_compt = output_td2.comp_time(ind);
td2_gdop = output_td2.GDOP(ind);
td2_horcov = NaN(1,length(ind));
td2_vercov = sqrt(output_td2.ned_cov(3,ind));
for i=1:length(ind)
    if (output_td2.ned_cov(1,i) == 400)
        td2_vercov(i) = NaN;
    else
        td2_horcov(i) = norm(output_td2.ned_cov(1:2,i));
    end
end
td2_horcov = td2_horcov(ind);
% ind = ~isnan(output_td2.hor_err);
% td3_hor_err = output_td3.hor_err(ind);
% td3_risk = output_td3.cost(ind);
% td3_meas = output_td3.num_meas_used(ind);
% ind = ~isnan(output_td4.hor_err);
% td4_hor_err = output_td4.hor_err(ind);
% td4_risk = output_td4.cost(ind);
% td4_meas = output_td4.num_meas_used(ind);

ind = ~isnan(output_ekf.cost);
ekf_hor_err = output_ekf.hor_err(ind);
ekf_ver_err = abs(output_ekf.ned_err(3,ind));
ekf_risk = output_ekf.cost(ind);
ekf_meas = output_ekf.num_meas_used(ind);
ekf_compt = output_ekf.comp_time(ind);
ekf_gdop = output_ekf.GDOP(ind);
ekf_horcov = zeros(1,length(ind));
ekf_vercov = sqrt(output_ekf.ned_cov(3,ind));
for i=1:length(ind)
    if (output_ekf.ned_cov(1,i) == 400)
        ekf_vercov(i) = NaN;
    else
        ekf_horcov(i) = norm(output_ekf.ned_cov(1:2,i));
    end
end
ekf_horcov = ekf_horcov(ind);

disp('Overall performance')
nonNaNCount = length(raps_b_hor_err);
fprintf('RAPS Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_b_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_b_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_b_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('RAPS Binary Hor Mean: %.2f\n', mean(raps_b_hor_err));
fprintf('RAPS Binary Hor RMS: %.2f\n', rms(raps_b_hor_err));
fprintf('RAPS Binary Hor Max: %.2f\n', max(raps_b_hor_err));
fprintf('RAPS Binary Ver Mean: %.2f\n', mean(raps_b_ver_err));
fprintf('RAPS Binary Ver RMS: %.2f\n', rms(raps_b_ver_err));
fprintf('RAPS Binary Ver Max: %.2f\n', max(raps_b_ver_err));

nonNaNCount = length(raps_nb_hor_err);
fprintf('RAPS Non-Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_nb_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_nb_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_nb_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Hor Mean: %.2f\n', mean(raps_nb_hor_err));
fprintf('RAPS Non-Binary Hor RMS: %.2f\n', rms(raps_nb_hor_err));
fprintf('RAPS Non-Binary Hor Max: %.2f\n', max(raps_nb_hor_err));
fprintf('RAPS Non-Binary Ver Mean: %.2f\n', mean(raps_nb_ver_err));
fprintf('RAPS Non-Binary Ver RMS: %.2f\n', rms(raps_nb_ver_err));
fprintf('RAPS Non-Binary Ver Max: %.2f\n', max(raps_nb_ver_err));

nonNaNCount = length(td2_hor_err);
fprintf('TD (lambda=2) Hor <= 1.0 m: %.2f%%\n', sum(td2_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('TD (lambda=2) Hor <= 1.5 m: %.2f%%\n', sum(td2_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('TD (lambda=2) Ver <= 3.0 m: %.2f%%\n', sum(td2_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('TD (lambda=2) Hor Mean: %.2f\n', mean(td2_hor_err));
fprintf('TD (lambda=2) Hor RMS: %.2f\n', rms(td2_hor_err));
fprintf('TD (lambda=2) Hor Max: %.2f\n', max(td2_hor_err));
fprintf('TD (lambda=2) Ver Mean: %.2f\n', mean(td2_ver_err));
fprintf('TD (lambda=2) Ver RMS: %.2f\n', rms(td2_ver_err));
fprintf('TD (lambda=2) Ver Max: %.2f\n', max(td2_ver_err));

nonNaNCount = length(ekf_hor_err);
fprintf('EKF Hor <= 1.0 m: %.2f%%\n', sum(ekf_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('EKF Hor <= 1.5 m: %.2f%%\n', sum(ekf_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('EKF Ver <= 3.0 m: %.2f%%\n', sum(ekf_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('EKF Hor Mean: %.2f\n', mean(ekf_hor_err));
fprintf('EKF Hor RMS: %.2f\n', rms(ekf_hor_err));
fprintf('EKF Hor Max: %.2f\n', max(ekf_hor_err));
fprintf('EKF Ver Mean: %.2f\n', mean(ekf_ver_err));
fprintf('EKF Ver RMS: %.2f\n', rms(ekf_ver_err));
fprintf('EKF Ver Max: %.2f\n', max(ekf_ver_err));
%%
purple = [0.4940, 0.1840, 0.5560];
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
green = [0.4660, 0.6740, 0.1880];
orange = [0.9290, 0.6940, 0.1250];

ekf_color = green;
td_color = orange;
raps_b_color = red;
raps_nb_color = blue;

%%
% Plot cost
figure(3)
tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(ekf_risk, '.', 'Color', ekf_color)
hold on
plot(td2_risk, '.', 'Color', td_color)
plot(raps_b_risk, '.', 'Color', raps_b_color)
plot(raps_nb_risk, '.', 'Color', raps_nb_color)
hold off
set(gca, 'XGrid', 'off', 'YGrid', 'on')
ylim([0.01,10^3])
axis tight
legend('KF', 'TD', 'RAPS-bi', 'RAPS-nb')
ylabel('Risk')
set(gca, 'YScale', 'log')
title('a. Risk over Time')
nexttile
plot(raps_b_penalty, '.', 'Color', raps_b_color)
hold on
plot(raps_nb_penalty, '.', 'Color', raps_nb_color)
axis tight
ylim([0,60])
grid on
legend('RAPS-bi', 'RAPS-nb')
ylabel('Value of the Penalty Term')
title('b. DiagRAPS Penalty Term')
nexttile
yyaxis left
plot(ekf_meas-td2_meas, '.', 'Color', td_color)
hold on
plot(ekf_meas-raps_b_meas, '.', 'Color', raps_b_color)
plot(ekf_meas-raps_nb_meas, '.', 'Color', raps_nb_color)
axis tight
ylim([0,60])
set(gca, 'YColor', 'k')
ylabel('No. of Meas.')
yyaxis right
plot(ekf_meas, '.', 'Color', purple)
set(gca, 'YColor', purple)
ylabel('Total No. of Meas.')
ylim([0,60])
grid on
xlabel('Epochs, sec')
title('c. No. of Meas. being Removed')
% nexttile
% plot(raps_nb_meas-raps_b_meas, '.', 'Color', 'k')


% Plot CDF
L = length(ekf_hor_err);
[f_ekf,x_ekf] = ecdf(ekf_hor_err);
[f_td,x_td] = ecdf(td2_hor_err);
[f_raps_b,x_raps_b] = ecdf(raps_b_hor_err);
[f_raps_nb,x_raps_nb] = ecdf(raps_nb_hor_err);

figure(4)
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
h_ekf = semilogx(x_ekf,f_ekf); hold on;
h_ekf.Color = ekf_color;
h_ekf.Marker = '*';
h_ekf.LineWidth = 0.2;
h_ekf.MarkerIndices = 1:200:L;

h_td = semilogx(x_td,f_td); hold on;
h_td.Color = td_color;
h_td.Marker = 'o';
h_td.LineWidth = 0.2;
h_td.MarkerIndices = 1:200:L;

h_raps_b = semilogx(x_raps_b,f_raps_b); hold on;
h_raps_b.Color = raps_b_color;
h_raps_b.Marker = '>';
h_raps_b.LineWidth = 0.2;
h_raps_b.MarkerIndices = 1:200:L;

h_raps_nb = semilogx(x_raps_nb,f_raps_nb); hold on;
h_raps_nb.Color = raps_nb_color;
h_raps_nb.Marker = '.';
h_raps_nb.LineWidth = 0.2;
h_raps_nb.MarkerIndices = 1:200:L;
xline(1.0)
xline(1.5)
grid on
xlim([0.1, 100])
legend([h_ekf,h_td,h_raps_b,h_raps_nb],{'KF', 'TD','RAPS-bi','RAPS-nb'})
legend('location','best');
currentXTicks = get(gca, 'XTick');
set(gca, 'XTick', unique([currentXTicks, 1.5]));
set(gca, 'XScale', 'log')
%title('CDF of Horizontal Error');
title('a. Cumulative Probability of HE');
ylabel('Probability');

nexttile
[f_ekf,x_ekf] = ecdf(ekf_ver_err);
[f_td,x_td] = ecdf(td2_ver_err);
[f_raps_b,x_raps_b] = ecdf(raps_b_ver_err);
[f_raps_nb,x_raps_nb] = ecdf(raps_nb_ver_err);
h_ekf = semilogx(x_ekf,f_ekf); hold on;
h_ekf.Color = ekf_color;
h_ekf.Marker = '*';
h_ekf.LineWidth = 0.2;
h_ekf.MarkerIndices = 1:200:L;

h_td = semilogx(x_td,f_td); hold on;
h_td.Color = td_color;
h_td.Marker = 'o';
h_td.LineWidth = 0.2;
h_td.MarkerIndices = 1:200:L;

h_raps_b = semilogx(x_raps_b,f_raps_b); hold on;
h_raps_b.Color = raps_b_color;
h_raps_b.Marker = '>';
h_raps_b.LineWidth = 0.2;
h_raps_b.MarkerIndices = 1:200:L;

h_raps_nb = semilogx(x_raps_nb,f_raps_nb); hold on;
h_raps_nb.Color = raps_nb_color;
h_raps_nb.Marker = '.';
h_raps_nb.LineWidth = 0.2;
h_raps_nb.MarkerIndices = 1:200:L;
xline(3)
grid on
legend([h_ekf,h_td,h_raps_b,h_raps_nb],{'KF', 'TD','RAPS-bi','RAPS-nb'})
legend('location','best');
xlim([0.1, 100])
currentXTicks = get(gca, 'XTick');
set(gca, 'XTick', unique([currentXTicks, 3.0]));
set(gca, 'XScale', 'log')
%title('CDF of Horizontal Error');
title('b. Cumulative Probability of VE');
xlabel('Error, meter')
ylabel('Probability');

% Define bin edges
binEdges = linspace(0, 0.1, 51);

% Create histograms to get bin counts and edges
[binCounts_b, binEdges_b] = histcounts(raps_b_compt, binEdges, 'Normalization', 'probability');
[binCounts_nb, binEdges_nb] = histcounts(raps_nb_compt, binEdges, 'Normalization', 'probability');
[binCounts_td, binEdges_td] = histcounts(td2_compt, binEdges, 'Normalization', 'probability');
[binCounts_ekf, binEdges_ekf] = histcounts(ekf_compt, binEdges, 'Normalization', 'probability');

% Calculate bin centers
binCenters_b = binEdges_b(1:end-1) + diff(binEdges_b) / 2;
binCenters_nb = binEdges_nb(1:end-1) + diff(binEdges_nb) / 2;
binCenters_td = binEdges_td(1:end-1) + diff(binEdges_td) / 2;
binCenters_ekf = binEdges_nb(1:end-1) + diff(binEdges_ekf) / 2;

% Plotting
% figure(5)
% plot(binCenters_b, binCounts_b, 'LineWidth', 2, 'Color', purple);
% hold on;
% plot(binCenters_nb, binCounts_nb, 'LineWidth', 2, 'Color', orange);
% grid on;
% xticks(0:0.02:0.15)
% xlabel('Time, second');
% ylabel('Probability');
% legend('RAPS-bi', 'RAPS-nb');

% Define bin edges
binEdges = linspace(0, 0.04, 50);
% Create histograms for each dataset
figure(5); % Open a new figure
% Histogram for RAPS-bi
histogram(raps_b_compt,binEdges,'Normalization','probability','FaceColor',raps_b_color,'FaceAlpha', 1);
hold on; % Retain the current plot when adding new plots
% Histogram for RAPS-nb
histogram(raps_nb_compt,binEdges,'Normalization','probability','FaceColor',raps_nb_color,'FaceAlpha', 0.7);
% Customize plot
grid on;
% xticks(0:0.02:0.1);
xlabel('Per epoch computation time, second');
ylabel('Probability');
legend('RAPS-bi', 'RAPS-nb');
hold off; % Release the plot hold

fprintf('KF Mean: %.5f, STD: %.5f, Max: %.5f, Min: %.5f.\n', mean(ekf_compt), std(ekf_compt), max(ekf_compt), min(ekf_compt))
fprintf('TD Mean: %.5f, STD: %.5f, Max: %.5f, Min: %.5f.\n', mean(td2_compt), std(td2_compt), max(td2_compt), min(td2_compt))
fprintf('RAPS-bi Mean: %.5f, STD: %.5f, Max: %.5f, Min: %.5f.\n', mean(raps_b_compt), std(raps_b_compt), max(raps_b_compt), min(raps_b_compt))
fprintf('RAPS-nb Mean: %.5f, STD: %.5f, Max: %.5f, Min: %.5f.\n', mean(raps_nb_compt), std(raps_nb_compt), max(raps_nb_compt), min(raps_nb_compt))

%% Split for feasible or non-feasible
% feasible
ind = output_raps_bp50.raps_flag == true;
L = sum(~isnan(output_raps_bp50.raps_flag));
fprintf('Percentage of feasible epoch (Binary): %.2f%%\n', sum(ind)/L*100);
raps_b_hor_err_f = output_raps_bp50.hor_err(ind);
raps_b_ver_err_f = abs(output_raps_bp50.ned_err(3,ind));
ind = output_raps_nbp50.raps_flag == true;
L = sum(~isnan(output_raps_nbp50.raps_flag));
fprintf('Percentage of feasible epoch (Non-Binary): %.2f%%\n', sum(ind)/L*100);
raps_nb_hor_err_f = output_raps_nbp50.hor_err(ind);
raps_nb_ver_err_f = abs(output_raps_nbp50.ned_err(3,ind));
% non-feasible
ind = output_raps_bp50.raps_flag == false;
raps_b_hor_err_nf = output_raps_bp50.hor_err(ind);
raps_b_ver_err_nf = abs(output_raps_bp50.ned_err(3,ind));
ind = output_raps_nbp50.raps_flag == false;
raps_nb_hor_err_nf = output_raps_nbp50.hor_err(ind);
raps_nb_ver_err_nf = abs(output_raps_bp50.ned_err(3,ind));
disp('Feasible epoch comparison')
nonNaNCount = length(raps_b_hor_err_f);
fprintf('RAPS Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_b_hor_err_f <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_b_hor_err_f <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_b_ver_err_f <= 3.0) / nonNaNCount * 100);
nonNaNCount = length(raps_nb_hor_err_f);
fprintf('RAPS Non-Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_nb_hor_err_f <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_nb_hor_err_f <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_nb_ver_err_f <= 3.0) / nonNaNCount * 100);

disp('Infeasible epoch comparison')
nonNaNCount = length(raps_b_hor_err_nf);
fprintf('RAPS Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_b_hor_err_nf <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_b_hor_err_nf <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_b_ver_err_nf <= 3.0) / nonNaNCount * 100);
nonNaNCount = length(raps_nb_hor_err_nf);
fprintf('RAPS Non-Binary Hor <= 1.0 m: %.2f%%\n', sum(raps_nb_hor_err_nf <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Hor <= 1.5 m: %.2f%%\n', sum(raps_nb_hor_err_nf <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Non-Binary Ver <= 3.0 m: %.2f%%\n', sum(raps_nb_ver_err_nf <= 3.0) / nonNaNCount * 100);

[f_raps_b_hor_f,x_raps_b_hor_f] = ecdf(raps_b_hor_err_f);
[f_raps_b_ver_f,x_raps_b_ver_f] = ecdf(raps_b_ver_err_f);
[f_raps_nb_hor_f,x_raps_nb_hor_f] = ecdf(raps_nb_hor_err_f);
[f_raps_nb_ver_f,x_raps_nb_ver_f] = ecdf(raps_nb_ver_err_f);
[f_raps_b_hor_nf,x_raps_b_hor_nf] = ecdf(raps_b_hor_err_nf);
[f_raps_b_ver_nf,x_raps_b_ver_nf] = ecdf(raps_b_ver_err_nf);
[f_raps_nb_hor_nf,x_raps_nb_hor_nf] = ecdf(raps_nb_hor_err_nf);
[f_raps_nb_ver_nf,x_raps_nb_ver_nf] = ecdf(raps_nb_ver_err_nf);
figure(6)
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
h_raps_b_hor_f = semilogx(x_raps_b_hor_f,f_raps_b_hor_f); hold on;
h_raps_b_hor_f.Color = raps_b_color;
h_raps_b_hor_f.Marker = '>';
h_raps_b_hor_f.LineWidth = 0.2;
h_raps_b_hor_f.MarkerIndices = 1:200:length(raps_b_hor_err_f);

h_raps_nb_hor_f = semilogx(x_raps_nb_hor_f,f_raps_nb_hor_f); hold on;
h_raps_nb_hor_f.Color = raps_nb_color;
h_raps_nb_hor_f.Marker = 'o';
h_raps_nb_hor_f.LineWidth = 0.2;
h_raps_nb_hor_f.MarkerIndices = 1:200:length(raps_nb_hor_err_f);
grid on
xlim([0.1, 100])
legend([h_raps_b_hor_f,h_raps_nb_hor_f],{'RAPS-bi','RAPS-nb'})
legend('location','best');
set(gca, 'XScale', 'log')
title('a. Hor. Feas.');
% xlabel('Horizontal error, meter');
ylabel('Probability');

nexttile
h_raps_b_ver_f = semilogx(x_raps_b_ver_f,f_raps_b_ver_f); hold on;
h_raps_b_ver_f.Color = raps_b_color;
h_raps_b_ver_f.Marker = '>';
h_raps_b_ver_f.LineWidth = 0.2;
h_raps_b_ver_f.MarkerIndices = 1:200:length(raps_b_hor_err_f);

h_raps_nb_ver_f = semilogx(x_raps_nb_ver_f,f_raps_nb_ver_f); hold on;
h_raps_nb_ver_f.Color = raps_nb_color;
h_raps_nb_ver_f.Marker = 'o';
h_raps_nb_ver_f.LineWidth = 0.2;
h_raps_nb_ver_f.MarkerIndices = 1:200:length(raps_nb_hor_err_f);
grid on
xlim([0.1, 100])
legend([h_raps_b_ver_f,h_raps_nb_ver_f],{'RAPS-bi','RAPS-nb'})
legend('location','best');
set(gca, 'XScale', 'log')
title('b. Ver. Feas.');
% xlabel('Horizontal error, meter');
ylabel('Probability');

nexttile
h_raps_b_hor_nf = semilogx(x_raps_b_hor_nf,f_raps_b_hor_nf); hold on;
h_raps_b_hor_nf.Color = raps_b_color;
h_raps_b_hor_nf.Marker = '>';
h_raps_b_hor_nf.LineWidth = 0.2;
h_raps_b_hor_nf.MarkerIndices = 1:200:length(raps_b_hor_err_nf);

h_raps_nb_hor_nf = semilogx(x_raps_nb_hor_nf,f_raps_nb_hor_nf); hold on;
h_raps_nb_hor_nf.Color = raps_nb_color;
h_raps_nb_hor_nf.Marker = 'o';
h_raps_nb_hor_nf.LineWidth = 0.2;
h_raps_nb_hor_nf.MarkerIndices = 1:200:length(raps_nb_hor_err_nf);
grid on
xlim([0.1, 100])
legend([h_raps_b_hor_nf,h_raps_nb_hor_nf],{'RAPS-bi','RAPS-nb'})
legend('location','best');
set(gca, 'XScale', 'log')
title('c. Hor. Infeas.');
xlabel('Error, meter');
ylabel('Probability');

nexttile
h_raps_b_ver_nf = semilogx(x_raps_b_ver_nf,f_raps_b_ver_nf); hold on;
h_raps_b_ver_nf.Color = raps_b_color;
h_raps_b_ver_nf.Marker = '>';
h_raps_b_ver_nf.LineWidth = 0.2;
h_raps_b_ver_nf.MarkerIndices = 1:200:length(raps_b_hor_err_nf);

h_raps_nb_ver_nf = semilogx(x_raps_nb_ver_nf,f_raps_nb_ver_nf); hold on;
h_raps_nb_ver_nf.Color = raps_nb_color;
h_raps_nb_ver_nf.Marker = 'o';
h_raps_nb_ver_nf.LineWidth = 0.2;
h_raps_nb_ver_nf.MarkerIndices = 1:200:length(raps_nb_hor_err_nf);
grid on
xlim([0.1, 100])
legend([h_raps_b_ver_nf,h_raps_nb_ver_nf],{'RAPS-bi','RAPS-nb'})
legend('location','best');
set(gca, 'XScale', 'log')
title('d. Ver. Infeas.');
xlabel('Error, meter');
ylabel('Probability');

%% GDOP vs Hor Acc
figure(7)
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
plot(ekf_horcov, ekf_hor_err, '.', 'Color', ekf_color)
hold on
plot(td2_horcov, td2_hor_err, '.', 'Color', td_color)
plot(raps_b_horcov, raps_b_hor_err, '.', 'Color', raps_b_color)
plot(raps_nb_horcov, raps_nb_hor_err, '.', 'Color', raps_nb_color)
plot([0.8,100],[0.8,100], 'Color', 'k')
xline(1.5)
legend('KF', 'TD', 'RAPS-bi', 'RAPS-nb', 'Consistent Line')
hold off
grid on
xlabel('HorSTD, meter');
ylabel('HE, meter');
xlim([0,100])
ylim([0.1,100])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

nexttile
plot(ekf_vercov, ekf_ver_err, '.', 'Color', ekf_color)
hold on
plot(td2_vercov, td2_ver_err, '.', 'Color', td_color)
plot(raps_b_vercov, raps_b_ver_err, '.', 'Color', raps_b_color)
plot(raps_nb_vercov, raps_nb_ver_err, '.', 'Color', raps_nb_color)
plot([0.8,100],[0.8,100], 'Color', 'k')
xline(3.0)
legend('KF', 'TD', 'RAPS-bi', 'RAPS-nb', 'Consistent Line')
hold off
grid on
xlabel('VerSTD, meter');
ylabel('VE, meter');
xlim([0,100])
ylim([0.1,100])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

disp('Horizontal')
fprintf('EKF Conservative Rate: %.2f%%\n', sum(ekf_hor_err <= ekf_horcov) / length(ekf_hor_err) * 100);
fprintf('TD Conservative Rate: %.2f%%\n', sum(td2_hor_err <= td2_horcov) / length(td2_hor_err) * 100);
fprintf('RAPS Binary Conservative Rate: %.2f%%\n', sum(raps_b_hor_err <= raps_b_horcov) / length(raps_b_hor_err) * 100);
fprintf('RAPS Non-binary Conservative Rate: %.2f%%\n', sum(raps_nb_hor_err <= raps_nb_horcov) / length(raps_nb_ver_err) * 100);
disp('Vertical')
fprintf('EKF Conservative Rate: %.2f%%\n', sum(ekf_ver_err <= ekf_vercov) / length(ekf_ver_err) * 100);
fprintf('TD Conservative Rate: %.2f%%\n', sum(td2_ver_err <= td2_vercov) / length(td2_ver_err) * 100);
fprintf('RAPS Binary Conservative Rate: %.2f%%\n', sum(raps_b_ver_err <= raps_b_vercov) / length(raps_b_ver_err) * 100);
fprintf('RAPS Non-binary Conservative Rate: %.2f%%\n', sum(raps_nb_ver_err <= raps_nb_vercov) / length(raps_nb_ver_err) * 100);