%init
close all;
clear all;
clc;

%config
MA_points = [1, 3, 5, 7];
CI_points = [1, 2, 5, 10];

%read
load Radar_sea.mat;

%raw
fig1 = figure;
IQ = abs(I + 1i * Q);
RTI = 10 * log10(IQ .^ 2);

for channel = 1 : (nch / 2)
    subplot(2, 2, channel)
    h = pcolor(T, R,  squeeze(RTI(channel, :, :)));
    set(h, "LineStyle", "none");
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

i = 1;
for point = MA_points
    subplot(2, 2, i)
    RTI = 10 * log10(movmean(IQ, point, 3) .^ 2);
    h = pcolor(T, R,  squeeze(RTI(channel, :, :)));
    set(h, "LineStyle", "none");
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    shading flat;
    title(point + "-points MA RTI ,Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
    i = i + 1;
end
saveas(fig2, "Q1_2.png");

%Coherent Integration
%ma
fig3 = figure;
channel = 1;

i = 1;
for point = CI_points
    subplot(2, 2, i)
    RTI = pulsint(I + 1i * Q, "coherent");
    h = pcolor(T, R,  squeeze(RTI(channel, :, :)));
    set(h, "LineStyle", "none");
    c = colorbar;
    title(c, "dB");
    colormap("jet");
    shading flat;
    title(point + "-points Coherent Integration RTI, Channel:" + channel);
    xlabel("Time(s)");
    ylabel("Range(km)");
    i = i + 1;
end
saveas(fig3, "Q1_3.png");
