%init
close all;
clear all;
clc;

%config
N = 2 ^ 9;

%read
load Radar_sea.mat;
load IQ_ma_ci.mat;

%val
channel = 1;
fs = 1 / t_res;
F = ((-N / 2) : 1 : N / 2 - 1) / N * fs * ipp * cint;
sIQ_raw = squeeze(IQ(channel, :, :));
sIQ_ma = squeeze(IQ_ma(channel, :, :));
sIQ_ci = squeeze(IQ_ci(channel, :, :));
sIQ_ma_ci = squeeze(IQ_ma_ci(channel, :, :));

D_raw = abs(fftshift(fft(sIQ_raw - mean(sIQ_raw), N, 2), 2)) / N;
D_ma = abs(fftshift(fft(sIQ_ma - mean(sIQ_ma), N, 2), 2)) / N;
D_ci = abs(fftshift(fft(sIQ_ci - mean(sIQ_ci), N, 2), 2)) / N;
D_ma_ci = abs(fftshift(fft(sIQ_ma_ci - mean(sIQ_ma_ci), N, 2), 2)) / N;

fig1 = figure;

%raw
subplot(2, 2, 1)
imagesc(F, R, D_raw);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
shading flat;
clim([0, 10]);
title("Raw RDI, Channel:" + channel + ", N:" + N);
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 2)
imagesc(F, R, D_ma);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
shading flat;
clim([0, 10]);
title("RDI, 7-point MA, Channel:" + channel + ", N:" + N);
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 3)
imagesc(F, R, D_ci);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
shading flat;
clim([0, 10]);
title("RDI, 7-point CI, Channel:" + channel + ", N:" + N);
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 4)
imagesc(F, R, D_ma_ci);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
shading flat;
clim([0, 10]);
title("RDI, 7-point MA, 7-point CI, Channel:" + channel + ", N:" + N);
xlabel("Freq(Hz)");
ylabel("Range(km)");

saveas(fig1, "Q3_1.png");