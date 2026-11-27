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
range_idxs = [35, 25, 15, 5];
f_radar = 52e6;
c = 3e8;
rain_peak_noise_ratio = 2;

%load
load Radar_Data.mat;

%val
IQ_size = size(IQ);
file_num = IQ_size(1);
CI_times = floor(ntrial / NFFT);
F = (-NFFT / 2 : NFFT / 2 - 1) / NFFT / t_res;
V = -(c * F) / (2 * f_radar);
[V, velocity_order] = sort(V, "ascend");
dv = V(2) - V(1);

%fft
IQ_FFT = fftshift(fft(IQ, NFFT, 3), 3) / NFFT;
IQ_FFT(:, :, NFFT / 2 + 1) = (IQ_FFT(:, :, NFFT / 2 + 2) + IQ_FFT(:, :, NFFT / 2)) / 2;
IQ_FFT = IQ_FFT(:, :, velocity_order);

%ci
IQ_FFT_CI = zeros(file_num, nrg1 - 1, NFFT);
for time = 1 : CI_times
    S = fftshift(fft(IQ(:, :, (time - 1) * NFFT + 1 : time * NFFT), NFFT, 3), 3) / NFFT;
    IQ_FFT_CI = IQ_FFT_CI + abs(S) .^ 2;
end
IQ_FFT_CI(:, :, NFFT / 2 + 1) = (IQ_FFT_CI(:, :, NFFT / 2 + 2) + IQ_FFT_CI(:, :, NFFT / 2)) / 2;
IQ_FFT_CI = IQ_FFT_CI / CI_times;
IQ_FFT_CI = IQ_FFT_CI(:, :, velocity_order);

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

V_a_idxs = V < 1.5 & V > -1.5;
V_r_idxs = V < -4 & V > -15;

V_a_mask = false(file_num, nrg1 - 1, length(V));
V_r_mask = false(file_num, nrg1 - 1, length(V));

P_a = nan(file_num, nrg1 - 1);
P_signal_a = P_a;
SNR_a = P_a;
va = P_a;
wa = P_a;

P_r = nan(file_num, nrg1 - 1);
P_signal_r = P_r;
SNR_r = P_r;
vr = P_r;
wr = P_r;

for f = 1 : file_num
    for r = 1 : nrg1 - 1
        spectrum = reshape(Sxx(f, r, :), 1, []);
        noise_level = P_noise(f, r);
        spectrum_clean = max(spectrum - noise_level, 0);

        % Find atmospheric and rain peaks within their Doppler velocity bands.
        spectrum_a = spectrum;
        spectrum_a(~V_a_idxs) = -Inf;
        [a_max, a_idx] = max(spectrum_a);

        spectrum_r = spectrum;
        spectrum_r(~V_r_idxs) = -Inf;
        [r_max, r_idx] = max(spectrum_r);

        has_rain = false;
        if R(r) <= 5.8 && r_max >= rain_peak_noise_ratio * noise_level
            [val_min, min_loc] = min(spectrum(r_idx : a_idx));
            m_border = r_idx + min_loc - 1;
            has_rain = val_min < 0.9 * r_max && val_min < 0.9 * a_max;
        end
        has_rain = true;

        if has_rain
            l_border = find(V <= V(r_idx) & V >= -15 & spectrum <= noise_level * 1.1, 1, "last");
            if isempty(l_border)
                l_border = find(V >= -15, 1, "first");
            end

            r_border = find(V >= V(a_idx) & V <= 3.5 & spectrum <= noise_level * 1.1, 1, "first");
            if isempty(r_border)
                r_border = find(V <= 3.5, 1, "last");
            end

            V_r_mask(f, r, l_border : m_border) = true;
            V_a_mask(f, r, m_border + 1 : r_border) = true;
        else
            l_border = find(V <= V(a_idx) & V >= -3.5 & spectrum <= noise_level * 1.1, 1, "last");
            if isempty(l_border)
                l_border = find(V >= -3.5, 1, "first");
            end

            r_border = find(V >= V(a_idx) & V <= 3.5 & spectrum <= noise_level * 1.1, 1, "first");
            if isempty(r_border)
                r_border = find(V <= 3.5, 1, "last");
            end
            V_a_mask(f, r, l_border : r_border) = true;
        end

        mask_a = reshape(V_a_mask(f, r, :), 1, []);
        signal_a = spectrum_clean(mask_a);
        signal_power_a = mean(signal_a);
        power_a = signal_power_a + noise_level;
        if sum(signal_a) > 0
            P_a(f, r) = power_a;
            P_signal_a(f, r) = signal_power_a;
            SNR_a(f, r) = signal_power_a / noise_level;
            va(f, r) = sum(signal_a .* V(mask_a)) / sum(signal_a);
            wa(f, r) = sqrt(sum(signal_a .* (V(mask_a) - va(f, r)) .^ 2) / sum(signal_a));
        end

        if has_rain
            mask_r = reshape(V_r_mask(f, r, :), 1, []);
            signal_r = spectrum_clean(mask_r);
            signal_power_r = mean(signal_r);
            power_r = signal_power_r + noise_level;
            if sum(signal_r) > 0
                P_r(f, r) = power_r;
                P_signal_r(f, r) = signal_power_r;
                SNR_r(f, r) = signal_power_r / noise_level;
                vr(f, r) = sum(signal_r .* V(mask_r)) / sum(signal_r);
                wr(f, r) = sqrt(sum(signal_r .* (V(mask_r) - vr(f, r)) .^ 2) / sum(signal_r));
            end
        end
    end
end

N_a = sum(V_a_mask, 3);
N_r = sum(V_r_mask, 3);
rain_plot_mask = R >= 1.5 & R <= 5.5;

%terminal speed
vt = vr - va;

%plot
figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = 10 * log10(abs(squeeze(IQ_FFT(file_idx, :, :))));
    h = pcolor(V, R, D);
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    yline(5.1, 'm--', '0°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    c = colorbar;
    c.Label.String = 'dB';
    colormap("jet");
    %clim([0, 8 * mean(D, "all")]);
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("RDI\nNFFT: %d", NFFT));
saveas(gcf, "RDI.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = 10 * log10(abs(squeeze(IQ_FFT_CI(file_idx, :, :))));
    h = pcolor(V, R, D);
    set(h,'EdgeColor', 'none');
    shading interp;
    yline(5.1, 'm--', '0°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    c = colorbar;
    c.Label.String = 'dB';
    colormap("jet");
    %clim([0, 8 * mean(D, "all")]);
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
D = abs(squeeze(IQ_FFT(1, 30, :)));
plot(V, D);
grid on;
xlim([-20, 10]);
xlabel("Velocity (m/s)");
ylabel("Intensity");

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
    ylim([1.5, 5.5]);
    legend(["Total", "Signal", "Noise", "SNR"], "Location","best");
    xlabel("Power (db)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Rain Power(Estimated by Moment Method)\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Power_r(Moment).png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = 10 * log10(abs(squeeze(IQ_FFT_CI(file_idx, :, :))));
    hold on;
    h = pcolor(V, R, D);
    plot(v(file_idx, :), R, "w", "LineWidth", 1.5);
    hold off;
    grid on;
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    c = colorbar;
    c.Label.String = 'dB';
    colormap("jet");
    %clim([0, 8 * mean(D, "all")]);
    xlim([-15, 15]);
    ylim([min(R), max(R)]);
    yline(5.1, 'm--', '0°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
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
width_max = max([wa(:); wr(:)], [], "omitnan");
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    hold on;
    plot(wa(file_idx, :), R);
    rain_width = wr(file_idx, :);
    rain_width(~rain_plot_mask) = NaN;
    plot(rain_width, R);
    hold off;
    grid on;
    xlim([0, width_max * 1.05]);
    xlabel("Spectral Width (m/s)");
    ylabel("Range (km)");
    legend(["Atmosphere", "Rain"]);
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Atmospheric and Rain Spectral Width\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Spectral_Width_a_r.png");

figure;
for file_idx = 1 : file_num
    subplot(2, 3, file_idx);
    D = 10 * log10(abs(squeeze(IQ_FFT_CI(file_idx, :, :))));
    va_line = fillmissing(va(file_idx, :), "linear", "EndValues", NaN);
    vr_line = fillmissing(vr(file_idx, :), "linear", "EndValues", NaN);
    vr_line(~rain_plot_mask) = NaN;
    vt_line = vr_line - va_line;
    hold on;
    h = pcolor(V, R, D);
    plot(va_line, R, "b", "LineWidth", 2);
    plot(vr_line, R, "k", "LineWidth", 2);
    plot(vt_line, R, "y", "LineWidth", 2);
    hold off;
    grid on;
    set(h,'EdgeColor', 'none');
    shading interp;
    axis xy;
    c = colorbar;
    c.Label.String = 'dB';
    colormap("jet");
    %clim([0, 8 * mean(D, "all")]);
    xlim([-15, 15]);
    ylim([min(R), R(max(range_idxs))]);
    yline(5.1, 'm--', '0°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    legend(["RDI", "V_a", "V_r", "V_t"]);
    xlabel("Velocity (m/s)");
    ylabel("Range (km)");
    title("Time: " + start_time(file_idx));
end
sgtitle(sprintf("Terminal Speed\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Terminal_Speed.png");

figure;
hold on;
for file_idx = 1 : file_num
    va_line = fillmissing(va(file_idx, :), "linear", "EndValues", NaN);
    vr_line = fillmissing(vr(file_idx, :), "linear", "EndValues", NaN);
    vr_line(~rain_plot_mask) = NaN;
    vt_line = vr_line - va_line;
    plot(vt_line, R, "o-", "LineWidth", 2);
end
hold off;
grid on;
xlim([-15, 0]);
ylim([min(R), R(max(range_idxs))]);
yline(5.1, 'm--', '0°C', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
legend(start_time(:));
xlabel("Velocity (m/s)");
ylabel("Range (km)");
sgtitle(sprintf("Terminal Speed\nNFFT: %d , NCI: %d", NFFT, CI_times));
saveas(gcf, "Terminal_Speed2.png");

figure;
i = 0;
for r = 1 : length(range_idxs)
    for file_idx = 1 : file_num
        i = i + 1;
        subplot(length(range_idxs), file_num, i);
        D = abs(squeeze(Sxx(file_idx, range_idxs(r), :)));
        fitted = P_noise(file_idx, range_idxs(r)) * ones(1, NFFT);
        fitted = fitted + (P_signal_a(file_idx, range_idxs(r)) * N_a(file_idx, range_idxs(r)) * abs(dv) / (sqrt(2 * pi) * wa(file_idx, range_idxs(r)))) * exp(-0.5 * ((V - va(file_idx, range_idxs(r))) / wa(file_idx, range_idxs(r))).^2);
        if rain_plot_mask(range_idxs(r)) && N_r(file_idx, range_idxs(r)) > 0 && wr(file_idx, range_idxs(r)) > 0
            fitted = fitted + (P_signal_r(file_idx, range_idxs(r)) * N_r(file_idx, range_idxs(r)) * abs(dv) / (sqrt(2 * pi) * wr(file_idx, range_idxs(r)))) * exp(-0.5 * ((V - vr(file_idx, range_idxs(r))) / wr(file_idx, range_idxs(r))).^2);
        end
        hold on;
        plot(V, D);
        plot(V, fitted, "LineWidth", 1.2);
        hold off;
        grid on;
        xlim([-20, 10]);
        legend(["P.S.", "G.F."], "Location", "northwest");
        xlabel("Velocity (m/s)");
        ylabel("Intensity");
        title(sprintf("Time: %s,\nRange: %.2f km", start_time(file_idx), R(range_idxs(r))));
    end
end
sgtitle(sprintf("Gaussian Fit\nNFFT: %d, NCI: %d", NFFT, CI_times));
saveas(gcf, "GaussianFit.png");