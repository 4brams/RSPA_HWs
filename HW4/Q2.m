%init
close all;
clear all;
clc;

%load
load Radar_sea.mat;

%val
channel = 1;

fig1 = figure;

i = 1;
for gate = [10, 50, 90, 130]
    figure(fig1);
    subplot(4, 2, i * 2 - 1)
    hold on;
    plot(T, squeeze(I(channel, gate, :)));
    plot(T, squeeze(Q(channel, gate, :)));
    hold off;
    grid on;
    xlim([min(T), max(T)]);
    xlim([0, 20]);
    title("I, Channel:" + channel + ", Range:" + range(gate) + "(km)");
    xlabel("Time(s)");
    ylabel("Amp.(V)");

    subplot(4, 2, i * 2)
    scatter(squeeze(I(1, gate, :)), squeeze(Q(1, gate, :)), 1, "filled");
    axis equal;
    grid on;
    xlim([-800, 800]);
    ylim([-800, 800]);
    title("IQ Scatter Plot, Channel:" + 1 + ", Range:" + range(gate) + "(km)");
    xlabel("I");
    ylabel("Q");
    i = i + 1;
end

%save
saveas(fig1, "Q2_1.png");