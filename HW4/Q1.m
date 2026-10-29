%init
close all;
clear all;
clc;

%function
function [IQ_ci, t_ci] = coherent_integration(IQ, T, point, time_length, range_length, channel_num)
    time_length = floor(time_length / point);

    IQ_ci = zeros(channel_num, range_length, time_length);
    t_ci = zeros(1, time_length);

    for i = 1:time_length
        idx = (i - 1) * point + 1 : i * point;
        t_ci(i) = T(idx(end));
        IQ_ci(:, :, i) = sum(IQ(:, :, idx), 3) / point;
    end
end

%config
MA_points = [1, 3, 5, 7];
CI_points = [1, 3, 5, 7];

%read
load Radar_sea.mat;

%data
IQ = I + 1i * Q;

%raw
fig1 = figure;

for channel = 1 : (nch / 2)
    RTI = 10 * log10(abs(squeeze(IQ(channel, :, :))) .^ 2);

    subplot(2, 2, channel)
    imagesc(T, R, RTI);
    axis xy;
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    shading flat;
    title("Raw RTI, Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
end
saveas(fig1, "Q1_1.png");

%ma
fig2 = figure;
channel = 1;

for i = 1:length(MA_points)
    IQ_ma = movmean(abs(IQ), MA_points(i), 3);
    RTI = 10 * log10(squeeze(IQ_ma(channel, :, :)) .^ 2);

    subplot(2, 2, i);
    imagesc(T, R,  RTI);
    axis xy;
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    shading flat;
    title("RTI, " + MA_points(i) + "-points MA ,Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
end
saveas(fig2, "Q1_2.png");

% Coherent Integration
fig3 = figure;
channel = 1;

for i = 1 : length(CI_points)
    [IQ_ci, T_ci] = coherent_integration(IQ, T, CI_points(i), ntrial, nrg1 - 1, nch / 2);
    RTI = 10 * log10(abs(squeeze(IQ_ci(channel, :, :))) .^ 2);

    subplot(2, 2, i);
    imagesc(T_ci, R, RTI);
    axis xy; 
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    title("RTI, " + CI_points(i) + "-points CI, Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
end
saveas(fig3, "Q1_3.png");

fig4 = figure;
[IQ_ma_ci, T_ci] = coherent_integration(IQ_ma, T, 5, ntrial, nrg1 - 1, nch / 2);
RTI = 10 * log10(abs(squeeze(IQ_ma_ci(channel, :, :))) .^ 2);
imagesc(T_ci, R, RTI);
axis xy; 
c = colorbar;
title(c, "dB");
colormap("jet");
title("RTI, 7-points MA, 7-points CI, Channel:" + channel);
xlabel("Time(s)");
ylabel("Range(km)");

saveas(fig4, "Q1_4.png");

save IQ_ma_ci IQ IQ_ma IQ_ci T_ci IQ_ma_ci;