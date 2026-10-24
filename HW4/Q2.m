%init
close all;
clear all;
clc;

%load
load Radar_sea.mat;

%val
channel = 1;

fig1 = figure;
fig2 = figure;

i = 1;
for gate = [10, 50, 90, 130]
    figure(fig1);
    subplot(4, 2, i * 2 - 1)
    plot(T, squeeze(I(channel, gate, :)));
    grid on;
    xlim([min(T), max(T)]);
    title("I, Channel:" + channel + ", Range:" + range(gate) + "(km)");
    xlabel("Time(s)");
    ylabel("I(t)");

    subplot(4, 2, i * 2)
    plot(T, squeeze(Q(channel, gate, :)));
    grid on;
    xlim([min(T), max(T)]);
    title("Q, Channel:" + channel + ", Range:" + range(gate) + "(km)");
    xlabel("Time(s)");
    ylabel("Q(t)");

    figure(fig2);
    subplot(2, 2, i)
    scatter(squeeze(I(1, gate, :)), squeeze(Q(1, gate, :)), 1, "filled");
    axis equal;
    grid on;
    title("IQ Scatter Plot, Channel:" + 1 + ", Range:" + range(gate) + "(km)");
    xlabel("I");
    ylabel("Q");
    i = i + 1;
end

%save
saveas(fig1, "Q2_1.png");
saveas(fig2, "Q2_2.png");