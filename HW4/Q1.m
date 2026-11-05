%init
close all;
clear all;
clc;

%function
function IQ_ci = coherent_integration(IQ, point, time_length, range_length, fc, ts)
    IQ_ci = zeros(range_length, time_length - point);
    
    for r = 1 : range_length
        for n = 1 : time_length - point
            for k = 0 : point - 1
                IQ_ci(r, n) = IQ_ci(r, n) + (IQ(r, n + k) * cos(2 * pi * fc * k * ts));
            end
        end
    end

    IQ_ci = IQ_ci / point;
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
    RTI = 20 * log10(squeeze(IQ(channel, :, :)));

    subplot(2, 2, channel)
    imagesc(T, R, abs(RTI));
    axis xy;
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    title("RTI, Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
end

sgtitle("Raw RTI of 3 Channels, Date: " + year_date_time);
saveas(fig1, "Q1_1.png");

%ma
fig2 = figure;
channel = 1;

for i = 1:length(MA_points)
    IQ_ma = movmean(IQ, MA_points(i), 3);
    RTI = 20 * log10(squeeze(IQ_ma(channel, :, :)));

    subplot(2, 2, i);
    imagesc(T, R,  abs(RTI));
    axis xy;
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    shading flat;
    title("RTI, " + MA_points(i) + "-MA");
    xlabel("Time(s)");
    ylabel("Range(km)");
end

sgtitle("RTI Moving Average, Channels: " + channel + ", Date: " + year_date_time);
saveas(fig2, "Q1_2.png");

% Coherent Integration
fig3 = figure;
channel = 1;

for i = 1 : length(CI_points)
    IQ_ci = coherent_integration(squeeze(IQ(channel, :, :)), CI_points(i), ntrial, nrg1 - 1, 52e6, t_res);
    RTI = 20 * log10(IQ_ci);

    subplot(2, 2, i);
    imagesc(T(1:end - CI_points(i)), R, abs(RTI));
    axis xy;
    c = colorbar;
    %clim([-50 0]);
    title(c, "dB");
    colormap("jet");
    title("RTI, " + CI_points(i) + "-CI");
    xlabel("Time(s)");
    ylabel("Range(km)");
end

sgtitle("RTI Coherent Integration, Channels: " + channel + ", Date: " + year_date_time);
saveas(fig3, "Q1_3.png");

%overlay
fig4 = figure;
IQ_3 = squeeze(mean(IQ, 1));
RTI = 20 * log10(IQ_3);

imagesc(T, R, abs(RTI));
axis xy; 
c = colorbar;
title(c, "dB");
colormap("jet");
title("RTI, Channel: 1+2+3");
xlabel("Time(s)");
ylabel("Range(km)");

saveas(fig4, "Q1_4.png");

save IQs IQ_ma IQ_ci IQ_3;