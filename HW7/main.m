%init
close all;
clear all;
clc;

%func
function P_noise = noise(s, r, dim)
    s = sort(s, dim, "ascend");
    P_noise = mean(s(:, :, 1 : floor(r * length(s(1, 1, :)))), 3);
end

function p = mean_power(Sxx)
    p = sum(Sxx, 3) / length(Sxx(1, 1, :));
end

function v = mean_doppler_velocity(Sxx, P_noise, V)

    V = reshape(V, 1, 1, []);
    P_noise = reshape(P_noise, size(P_noise,1), size(P_noise,2), 1);

    S_clean = max(Sxx - P_noise, 0);

    numerator = sum(S_clean .* V, 3);
    denominator = sum(S_clean, 3);

    v = numerator ./ denominator;
    v(denominator == 0) = NaN;
end

function sigma = mean_spectral_width(Sxx, P_noise, V, mu)
    V = reshape(V, 1, 1, []);
    mu = reshape(mu, size(mu,1), size(mu,2), 1);
    P_noise = reshape(P_noise, size(P_noise, 1), size(P_noise, 2), 1);
    sigma = sqrt(sum((Sxx - P_noise) .* (V - mu) .^ 2, 3) ./ sum((Sxx - P_noise), 3));
end

%config
NFFT = 2 ^ 10;
range_idxs = [30, 20, 10, 5];
f_radar = 52e6;
c = 3e8;

%load
load Radar_Data.mat;

%val
IQ_size = size(IQ);
file_num = IQ_size(1);
CI_times = floor(ntrial / NFFT);
F = (-NFFT / 2 : NFFT / 2 - 1) / NFFT / t_res;
V = -(c * F) / (2 * f_radar);
dv = V(2) - V(1);

V_a_idxs = V < 3 & V > -3;
V_r_idxs = V < -3 & V > -15;
N_a = sum(V_a_idxs);   % number of velocity bins in atmospheric window
N_r = sum(V_r_idxs);   % number of velocity bins in rain window

%fft
IQ_FFT = fftshift(fft(IQ, NFFT, 3), 3) / NFFT;
IQ_FFT(:, :, NFFT / 2 + 1) = (IQ_FFT(:, :, NFFT / 2 + 2) + IQ_FFT(:, :, NFFT / 2)) / 2;

%ci
IQ_FFT_CI = zeros(file_num, nrg1 - 1, NFFT);
for time = 1 : CI_times
    S = fftshift(fft(IQ(:, :, (time - 1) * NFFT + 1 : time * NFFT), NFFT, 3), 3) / NFFT;
    IQ_FFT_CI = IQ_FFT_CI + abs(S) .^ 2;
end
IQ_FFT_CI(:, :, NFFT / 2 + 1) = (IQ_FFT_CI(:, :, NFFT / 2 + 2) + IQ_FFT_CI(:, :, NFFT / 2)) / 2;
IQ_FFT_CI = IQ_FFT_CI / CI_times;

%xorr
IQ_corr = zeros(file_num, nrg1 - 1, 2 * ntrial - 1);
lag = zeros(file_num, nrg1 - 1, 2 * ntrial - 1);
P_corr = zeros(file_num, nrg1 - 1);
P_signal_corr = zeros(file_num, nrg1 - 1);
for file_idx = 1 : file_num
    for r = 1 : nrg1 - 1
        [IQ_corr(file_idx, r, :), lag(file_idx, r, :)] = xcorr(squeeze(IQ(file_idx, r, :)), "unbiased");
        P_corr(file_idx, r) = real(IQ_corr(file_idx, r, ntrial));
        P_signal_corr(file_idx, r) = real((IQ_corr(file_idx, r, ntrial + 1) + IQ_corr(file_idx, r, ntrial - 1)) / 2);
    end
end
P_noise_corr = P_corr - P_signal_corr;
P_corr = P_corr / length(V);
P_noise_corr = P_noise_corr / length(V);
P_signal_corr = P_signal_corr / length(V);
lag = lag * t_res;

%total
% IQ_FFT_CI already stores averaged power spectra (NCI), no need to square again
Sxx = IQ_FFT_CI;
P_noise = noise(Sxx, 0.3, 3);
P = mean_power(Sxx);
P_signal = P - P_noise;
SNR = P_signal ./ P_noise;
v = mean_doppler_velocity(Sxx, P_noise, V);
w = mean_spectral_width(Sxx, P_noise, V, v);

%atm
P_a = mean_power(Sxx(:, :, V_a_idxs));
P_signal_a = P_a - P_noise;
SNR_a = P_signal_a ./ P_noise;
va = mean_doppler_velocity(Sxx(:, :, V_a_idxs), P_noise, V(V_a_idxs));
wa = mean_spectral_width(Sxx(:, :, V_a_idxs), P_noise, V(V_a_idxs), va);

%rain
P_r = mean_power(Sxx(:, :, V_r_idxs));
P_signal_r = P_r - P_noise;
SNR_r = P_signal_r ./ P_noise;
vr = mean_doppler_velocity(Sxx(:, :, V_r_idxs), P_noise, V(V_r_idxs));
wr = mean_spectral_width(Sxx(:, :, V_r_idxs), P_noise, V(V_r_idxs), vr);

%terminal speed
vt = vr - va;

%plot
figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = abs(squeeze(IQ_FFT(file_idx, :, :)));
    h = pcolor(V, R, D);
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    colorbar;
    colormap("jet");
    clim([0, 8 * mean(D, "all")]);
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("RDI\nNFFT: %d", NFFT));
saveas(gcf, "RDI.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = abs(squeeze(IQ_FFT_CI(file_idx, :, :)));
    h = pcolor(V, R, D);
    set(h,'EdgeColor', 'none');
    shading interp;
    colorbar;
    colormap("jet");
    clim([0, 8 * mean(D, "all")]);
    xlabel("v (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("RDI\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "RDI_CI.png");

figure;
i = 0;
for r = 1 : length(range_idxs)
    for file_idx = 1 : file_num
        i = i + 1;
        subplot(length(range_idxs), file_num, i);
        D = abs(squeeze(IQ_FFT(file_idx, range_idxs(r), :)));
        plot(V, D);
        grid on;
        xlabel("Velocity (m/s)");
        ylabel("Intensity");
        title(sprintf("Time: %s,\nRange: %.2f km", start_time(file_idx), R(range_idxs(r))));
    end
end
sgtitle(sprintf("Spectrum\nNFFT: %d", NFFT));
saveas(gcf, "Spectrum.png");

figure;
i = 0;
for r = 1 : length(range_idxs)
    for file_idx = 1 : file_num
        i = i + 1;
        subplot(length(range_idxs), file_num, i);
        re = abs(squeeze(IQ_corr(file_idx, range_idxs(r), :)));
        l = squeeze(lag(file_idx, range_idxs(r), :));
        plot(l, re);
        grid on;
        xlabel("lag (s)");
        ylabel("r");
        title(sprintf("Time: %s,\nRange: %.2f km", start_time(file_idx), R(range_idxs(r))));
    end
end
sgtitle(sprintf("Autocorrelation\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Autocorrelation.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    hold on;
    plot(10 * log10(P_corr(file_idx, :)), R);
    plot(10 * log10(P_signal_corr(file_idx, :)), R);
    plot(10 * log10(P_noise_corr(file_idx, :)), R);
    plot(10 * log10(SNR(file_idx, :)), R);
    hold off;
    grid on;
    legend(["Total", "Signal", "Noise", "SNR"]);
    xlabel("Power (db)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Power(Estimated by Autocorrelation Function)\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Power(Auto).png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    hold on;
    plot(10 * log10(P(file_idx, :)), R);
    plot(10 * log10(P_signal(file_idx, :)), R);
    plot(10 * log10(P_noise(file_idx, :)), R);
    plot(10 * log10(SNR(file_idx, :)), R);
    hold off;
    grid on;
    legend(["Total", "Signal", "Noise", "SNR"]);
    xlabel("Power (db)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Power(Estimated by Moment Method)\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Power(Moment).png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    hold on;
    plot(10 * log10(P_a(file_idx, :)), R);
    plot(10 * log10(P_signal_a(file_idx, :)), R);
    plot(10 * log10(P_noise(file_idx, :)), R);
    plot(10 * log10(SNR_a(file_idx, :)), R);
    hold off;
    grid on;
    legend(["Total", "Signal", "Noise", "SNR"]);
    xlabel("Power (db)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("ATM Power(Estimated by Moment Method)\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Power_a(Moment).png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    hold on;
    plot(10 * log10(P_r(file_idx, :)), R);
    plot(10 * log10(P_signal_r(file_idx, :)), R);
    plot(10 * log10(P_noise(file_idx, :)), R);
    plot(10 * log10(SNR_r(file_idx, :)), R);
    hold off;
    grid on;
    legend(["Total", "Signal", "Noise", "SNR"]);
    xlabel("Power (db)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Rain Power(Estimated by Moment Method)\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Power_r(Moment).png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = abs(squeeze(IQ_FFT_CI(file_idx, :, :)));
    hold on;
    h = pcolor(V, R, D);
    plot(v(file_idx, :), R, "w", "LineWidth", 1.5);
    hold off;
    grid on;
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    colorbar;
    colormap("jet");
    clim([0, 8 * mean(D, "all")]);
    xlim([-15, 15]);
    ylim([min(R), max(R)]);
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Mean Doppler Speed\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Mean_Doppler_Speed.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    plot(w(file_idx, :), R);
    grid on;
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Spectral Width\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Spectral_Width.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = abs(squeeze(IQ_FFT_CI(file_idx, :, :)));
    hold on;
    h = pcolor(V, R, D);
    plot(va(file_idx, :), R, "b", "LineWidth", 2);
    plot(vr(file_idx, :), R, "k", "LineWidth", 2);
    plot(vt(file_idx, :), R, "y", "LineWidth", 2);
    hold off;
    grid on;
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    colorbar;
    colormap("jet");
    clim([0, 8 * mean(D, "all")]);
    xlim([-15, 15]);
    ylim([min(R), 5.5]);
    legend(["RDI", "V_a", "V_r", "V_t"]);
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Terminal Speed\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Terminal_Speed.png");

figure;
i = 0;
for r = 1 : length(range_idxs)
    for file_idx = 1 : file_num
        i = i + 1;
        subplot(length(range_idxs), file_num, i);
        D = abs(squeeze(Sxx(file_idx, range_idxs(r), :)));
        fitted = P_noise(file_idx, range_idxs(r)) * ones(1, NFFT);
        fitted = fitted + (P_signal_a(file_idx, range_idxs(r)) * N_a * abs(dv) / (sqrt(2 * pi) * wa(file_idx, range_idxs(r)))) * exp(-0.5 * ((V - va(file_idx, range_idxs(r))) / wa(file_idx, range_idxs(r))).^2);
        fitted = fitted + (P_signal_r(file_idx, range_idxs(r)) * N_r * abs(dv) / (sqrt(2 * pi) * wr(file_idx, range_idxs(r)))) * exp(-0.5 * ((V - vr(file_idx, range_idxs(r))) / wr(file_idx, range_idxs(r))).^2);
        hold on;
        plot(V, D);
        plot(V, fitted);
        hold off;
        grid on;
        legend(["P.S.", "G.F."]);
        xlabel("Velocity (m/s)");
        ylabel("Intensity");
        title(sprintf("Time: %s,\nRange: %.2f km", start_time(file_idx), R(range_idxs(r))));
    end
end
sgtitle(sprintf("Gaussian Fit\nNFFT: %d, NCI: %d", NFFT, CI_times));
saveas(gcf, "GaussianFit.png");