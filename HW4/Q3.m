%init
close all;
clear all;
clc;

%config
N = 2 ^ 12;

%read
load Radar_sea.mat;

%val
fs = 1 / t_res;
F = ((-N / 2) : 1 : N / 2 - 1) / N * fs * ipp * cint;
S = movmean(I + 1i * Q, 5);
D = fftshift(fft(S - mean(S, 3), N, 3), 3) / N;

fig1 = figure;
%plot
for channel = 1 : (nch / 2)
    subplot(2, 2, channel)
    %surf(F,R,abs(squeeze(D(channel, :, :))) .^ 2);
    h = pcolor(F, R, abs(squeeze(D(channel, :, :))) .^ 2);
    set(h, "LineStyle", "none");
    c = colorbar;
    title(c, "Intensity");
    colormap("jet");
    shading flat;
    clim([0, 1]);
    title("RDI, Channel:" + channel + ", N:" + N);
    xlabel("Freq(Hz)");
    ylabel("Range(km)");
    %V = (D * 300000000) / (2 * fs);
    %v = (f * 300000000) / (2 * fs);
end

saveas(fig1, "Q3_1.png");