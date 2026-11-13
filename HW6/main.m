%init
close all;
clear all;
clc;

%func
function P_noise = noise(Sxx, N)
    Sxx = sort(Sxx(:), "ascend");
    P_noise = mean(Sxx((1 : N * 0.33)')) * N;
    disp("Estimate Noise Power: " + P_noise);
end

function [P_signal, fittedSxx] = gaussianMomentSignalPower(f, Sxx, P_noise, peakCenters, halfWidthBins)
    df = abs(f(2) - f(1));
    noisePerBin = P_noise / numel(Sxx);
    P_signal = 0;
    fittedSxx = noisePerBin * ones(size(Sxx));

    for k = 1 : numel(peakCenters)
        idx = abs(f - peakCenters(k)) <= halfWidthBins * df;
        fp = f(idx);
        w = max(real(Sxx(idx)) - noisePerBin, 0);
        W = sum(w);

        f0 = sum(fp .* w) / W;
        sigma = sqrt(sum((fp - f0).^2 .* w) / W);

        P_signal = P_signal + W;
        fittedSxx = fittedSxx + (W * df / (sqrt(2 * pi) * sigma)) * exp(-0.5 * ((f - f0) / sigma).^2);
        
        disp("Gaussian peak " + k + ": f0=" + f0 + " Hz, sigma=" + sigma + " Hz, signal power=" + W);
    end

    disp("Estimate Signal Power: " + P_signal);
end

function [f, S, Sxx] = powerSpectrum(s, N, dt, name, peakCenters)
    f = (- N / 2 : N / 2 - 1) / N / dt;
    S = fftshift(fft(s, N)) / N;
    Sxx = S .* conj(S);
    P_noise = noise(Sxx, N);
    [P_signal, gaussianFit] = gaussianMomentSignalPower(f, Sxx, P_noise, peakCenters, 8);
    
    disp("Estimate SNR: " + P_signal / P_noise);

    figure;
    
    subplot(2, 1, 1);
    hold on;
    plot(f, abs(Sxx));
    plot(f, gaussianFit, "r--");
    hold off;
    grid on;
    xlim([-20, 20]);
    xlabel("f(Hz)");
    ylabel("Power(V^2)");
    legend(["S_{xx}", "Gaussian fit"]);
    
    subplot(2, 1, 2);
    hold on;
    plot(f, 10 * log10(abs(Sxx)));
    plot(f, 10 * log10(gaussianFit), "r--");
    hold off;
    grid on;
    xlim([-20, 20]);
    xlabel("f(Hz)");
    ylabel("Power(dB)");
    legend(["S_{xx}", "Gaussian fit"]);
    
    sgtitle(sprintf("Power Spectrum\nN-FFT: %d", N));
    saveas(gcf, name);
end

function corr(lag, r, name)
    figure;
    
    subplot(2, 1, 1);
    plot(lag, abs(r));
    grid on;
    xlim([0.25 * min(lag), 0.25 * max(lag)]);
    xlabel("Lag(s)");
    ylabel("|R_{ss}|");
    title("Amplitude");
    
    subplot(2, 1, 2);
    plot(lag, atan2d(imag(r), real(r)));
    grid on;
    xlim([0.25 * min(lag), 0.25 * max(lag)]);
    xlabel("Lag(s)");
    ylabel("Phase(deg)");
    title("Phase");
    
    sgtitle("Autocorrelation");
    saveas(gcf, name);
end

function print_arguments(lag, r, dt, n1, n2, name)
    idx0 = lag == 0;
    idx1 = lag == dt;
    R0 = r(idx0);
    fd = -angle(r(idx1)) / (2 * pi * dt);
    sigma1 = std(n1);
    sigma2 = std(n2);
    P_noise = sigma1 ^ 2 + sigma2 ^ 2;
    P_signal = real(R0) - P_noise;
    SNR_dB = 10 * log10(P_signal / P_noise);
    %disp(name);
    disp("Total Power: " + real(R0));
    disp("Signal Power: " + P_signal);
    disp("Noise Power: " + P_noise);
    disp("Frequency: " + fd + " Hz");
    disp("SNR: " + P_signal / P_noise + "");
end

%Q1
disp("Q1");
N = 2 ^ 9;
NFFT = 2 ^ 11;
A = 5;
B = 2;
f1 = 1;
f2 = -3;
fs = 50;
dt = 1 / fs;

t = 0 : dt : (N - 1) * dt;
n = 1 * randn(1, N);%mean:0 std:1
disp("Mean: " + mean(n));
disp("STD: " + std(n));
s = (A * cos(2 * pi * f2 * t) + B * cos(2 * pi * f1 * t) + n) + 1i * (A * sin(2 * pi * f2 * t) + B * cos(2 * pi * f1 * t) + n);

[r, lag] = xcorr(s, "unbiased");
%r = r / N;
lag = lag * dt;

powerSpectrum(s, NFFT, dt, "Q1_2.png", [f2, -f1, f1]);
corr(lag, r, "Q1_3.png");
print_arguments(lag, r, dt, n, n, "Q1");

%Q2
disp("Q2");
f3 = -0.5;

n = 2 * randn(1, N);%mean:0 std:2
disp("Mean: " + mean(n));
disp("STD: " + std(n));
s = (A * cos(2 * pi * f3 * t) + n) + 1i * (A * sin(2 * pi * f3 * t) + n);

[r, lag] = xcorr(s, "unbiased");
%r = r / N;
lag = lag * dt;

powerSpectrum(s, NFFT, dt, "Q2_2.png", f3);
corr(lag, r, "Q2_3.png");
print_arguments(lag, r, dt, n, n, "Q2");

%Q3
disp("Q3");
n1 = 2 * randn(1, N);%mean:0 std:2
n2 = 2 * randn(1, N);%mean:0 std:2
disp("n1 mean: " + mean(n1));
disp("n1 STD: " + std(n1));
disp("n1 mean: " + mean(n1));
disp("n1 STD: " + std(n1));
s = (A * cos(2 * pi * f3 * t) + n1) + 1i * (B * sin(2 * pi * f3 * t) + n2);

[r, lag] = xcorr(s, "unbiased");
%r = r / N;
lag = lag * dt;

powerSpectrum(s, NFFT, dt, "Q3_2.png", [f3, -f3]);
corr(lag, r, "Q3_3.png");
print_arguments(lag, r, dt, n1, n2, "Q3");