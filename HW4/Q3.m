%init
close all;
clear all;
clc;

%config
N = 2 ^ 10;

%read
load Radar_sea.mat;
load IQs.mat;

%val
channel = 1;
fs = 1 / t_res;
F = ((-N / 2) : 1 : N / 2 - 1) / N * fs * ipp * cint;

sIQ_raw = squeeze(IQ(channel, :, :));
sIQ_ma = squeeze(IQ_ma(channel, :, :));
sIQ_ci = IQ_ci;
sIQ_3 = IQ_3;

D_raw = abs(fftshift(fft(sIQ_raw, N, 2), 2)) / N;
D_ma = abs(fftshift(fft(sIQ_ma, N, 2), 2)) / N;
D_ci = abs(fftshift(fft(sIQ_ci, N, 2), 2)) / N;
D_3 = abs(fftshift(fft(sIQ_3, N, 2), 2)) / N;

D_raw(:, N / 2 + 1) = (D_raw(:, N / 2 + 2) + D_raw(:, N / 2)) / 2;
D_ma(:, N / 2 + 1) = (D_ma(:, N / 2 + 2) + D_ma(:, N / 2)) / 2;
D_ci(:, N / 2 + 1) = (D_ci(:, N / 2 + 2) + D_ci(:, N / 2)) / 2;
D_3(:, N / 2 + 1) = (D_3(:, N / 2 + 2) + D_3(:, N / 2)) / 2;

fig1 = figure;

%raw
subplot(2, 2, 1)
imagesc(F, R, D_raw);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_raw, "all")]);
title("Raw RDI");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 2)
imagesc(F, R, D_ma);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_ma, "all")]);
title("RDI, 7-point MA");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 3)
imagesc(F, R, D_ci);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_ci, "all")]);
title("RDI, 7-point CI");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma ci
subplot(2, 2, 4)
imagesc(F, R, D_3);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_3, "all")]);
title("RDI, Channel:1+2+3");
xlabel("Freq(Hz)");
ylabel("Range(km)");

sgtitle("RDI, N:" + N + ", Date: " + year_date_time);
saveas(fig1, "Q3_1.png");

%ci
times = floor(ntrial / N);

D_raw = zeros(nrg1 - 1, N);
D_ma = zeros(nrg1 - 1, N);
D_ci = zeros(nrg1 - 1, N);
D_3 = zeros(nrg1 - 1, N);

for i = 1 : nrg1 - 1
    for j = 0 : times - 1
        D_raw(i, :) = D_raw(i, :) + abs(fftshift(fft(sIQ_raw(i, j * N + 1: (j + 1) * N), N, 2), 2)) / N;
        D_ma(i, :) = D_ma(i, :) + abs(fftshift(fft(sIQ_ma(i, j * N + 1: (j + 1) * N), N, 2), 2)) / N;
        D_ci(i, :) = D_ci(i, :) + abs(fftshift(fft(sIQ_ci(i, j * N + 1: (j + 1) * N), N, 2), 2)) / N;
        D_3(i, :) = D_3(i, :) + abs(fftshift(fft(sIQ_3(i, j * N + 1: (j + 1) * N), N, 2), 2)) / N;
    end
end

D_raw(:, N / 2 + 1) = (D_raw(:, N / 2 + 2) + D_raw(:, N / 2)) / 2;
D_ma(:, N / 2 + 1) = (D_ma(:, N / 2 + 2) + D_ma(:, N / 2)) / 2;
D_ci(:, N / 2 + 1) = (D_ci(:, N / 2 + 2) + D_ci(:, N / 2)) / 2;
D_3(:, N / 2 + 1) = (D_3(:, N / 2 + 2) + D_3(:, N / 2)) / 2;

fig2 = figure;

%raw
subplot(2, 2, 1)
imagesc(F, R, D_raw);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_raw, "all")]);
title("Raw RDI");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 2)
imagesc(F, R, D_ma);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_ma, "all")]);
title("RDI, 7-point MA");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma
subplot(2, 2, 3)
imagesc(F, R, D_ci);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_ci, "all")]);
title("RDI, 7-point CI");
xlabel("Freq(Hz)");
ylabel("Range(km)");

%ma ci
subplot(2, 2, 4)
imagesc(F, R, D_3);
axis xy; 
c = colorbar;
title(c, "Intensity");
colormap("jet");
clim([0, 5 * mean(D_3, "all")]);
title("RDI, Channel:1+2+3");
xlabel("Freq(Hz)");
ylabel("Range(km)");

sgtitle("RDI, " + times + "-times CI, N:" + N + ", Date: " + year_date_time);
saveas(fig2, "Q3_1.png");