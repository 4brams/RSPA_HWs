%init
close all;
clear all;
clc;

%config
FILE_NAME = "SN_m_tot_V2.0.txt";
DELTA_T = 1;
LINE_WIDTH = 0.1;
FIGURE_WIDTH = 3000;
FIGURE_HEIGHT = 2000;

points = [1, 5, 9, 13];
%points = 1:10:2001;
%points = [1, 49, 299, 699];
T = [];

%func
function H = compute_H(h, w, m, dt)
    H = 0;
    for n = -m : m
        H = H + h(n + m + 1) * cos(w * n * dt);
    end
end

function H = systemResponse(f, h, m, dt)
    N = length(f);
    H = zeros(1, N);
    HIdx = 1;
    for w = 2*pi*f
        H(HIdx) = compute_H(h, w, m, dt);
        HIdx = HIdx + 1;
    end
end

function [t, f, sunspotMA, sunspotSpectrum, filterSystemResponse, sunPeriod] = processData(sunspotNum, period, point, delta_t)
    N = 2^11;
    margin = floor(point / 2);
    
    t = period((margin + 1) : (end - margin));
    f = ((-N / 2) : 1 : (N / 2 - 1)) / (delta_t * N);
    
    filter = zeros(1, point);
    for i = 1 : point
        filter(i) = 1 / i;
    end
    
    sunspotMA = conv(sunspotNum, filter, "valid");
    sunspotSpectrum = fftshift(fft(sunspotMA - mean(sunspotMA), N)) / N;
    filterSystemResponse = systemResponse(f, filter, margin, delta_t);
    
    [val, idx] = max(abs(sunspotSpectrum));
    sunPeriod = abs((1 / (f(idx))) / 12);
end

%main
fig1 = figure("Position", [0, 0, FIGURE_WIDTH, FIGURE_HEIGHT]);
hold on;
fig2 = figure;
hold on;

[sunspotNum, period, sizeofData] = readData(FILE_NAME);
legends = ["1-point", "5-points", "9-points", "13-points"];
idx = 0;
for point = points
    [t, f, sunspotMA, sunspotSpectrum, filterSystemResponse, T(end + 1)] = processData(sunspotNum, period, point, DELTA_T);

    figure(fig1);

    subplot(length(points), 3, idx * 3 + 1);
    plot(t, sunspotMA, LineWidth = LINE_WIDTH);
    grid on;
    xlim([min(period), max(period)]);
    xlabel("Time(Month)");
    ylabel("s(t)");
    title(point + "-MA Sumspots Number");

    subplot(length(points), 3, idx * 3 + 2);
    plot(f, abs(sunspotSpectrum), LineWidth = LINE_WIDTH);
    grid on;
    xlim([-0.035, 0.035]);
    xlabel("Frequency");
    ylabel("S(f)");
    title("Spectrum of " + point + "-MA Sumspots Number (T=" + T(end) + "years)");

    subplot(length(points), 3, idx * 3 + 3);
    plot(2 * pi * f, filterSystemResponse, LineWidth = LINE_WIDTH);
    grid on;
    xlim([min(2 * pi * f), max(2 * pi * f)]);
    xlabel("\omega");
    ylabel("H(\omega)");
    title(point + "-Point Filter System Response");

    figure(fig2);

    hold on;
    plot(t, sunspotMA, LineWidth = LINE_WIDTH);
    grid on;
    ylabel('s(t)');
    xlabel("Time(Month)");
    xlim([min(t), max(t)]);
    title("Sunspot Number of Different Point MA");

    idx = idx + 1;

    disp(point);
end

figure(fig2);
legend(legends);

%save
saveas(fig1, "Figure1.png");
saveas(fig2, "Figure2.png");