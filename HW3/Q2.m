%init
close all;
clear all;
clc;

%config
N1 = 2 ^ 10;
N2 = 2 ^ 9;
dt = 1e-2;
mu = 0;
sigma = 5;

%x
t1 = (- N1 / 2 : N1 / 2 - 1) * dt;
f1 = (- N1 / 2 : N1 / 2 - 1) / (N1 * dt);

x = -20 : 1 : 20;

t2 = (-N2 / 2 : N2 / 2 - 1) * dt;
f2 = (- N2 / 2 : N2 / 2 - 1) / (N2 * dt);

%y
y1 = sigma * randn(1, N1) + mu;
y2 = sigma * randn(1, N2) + mu;

g = (1 / (sigma * ((2 * pi) ^ (1 / 2)))) * exp(-((x - mu) .^ 2)/(2 * (sigma ^ 2)));

Y1 = abs(fftshift(fft(y1, N1))) / N1;
Y2 = abs(fftshift(fft(y2, N2))) / N2;

[r1, lags1] = xcorr(y1);
[r2, lags2] = xcorr(y2);

%plot
figure;

subplot(4, 2, 1);
plot(t1, y1);
grid on;
xlim([min(t1), max(t1)]);
xlabel("t");
ylabel("y1");
title("Random Numbers of Gauss Distrbution(N = " + N1 + ")");

subplot(4, 2, 2);
plot(t2, y2);
grid on;
xlim([min(t2), max(t2)]);
xlabel("t");
ylabel("y2");
title("Random Numbers of Gauss Distrbution(N = " + N2 + ")");

subplot(4, 2, 3);
hold on;
histogram(y1, x, 'Normalization', 'pdf');
xticks(x);
xticklabels(x);
plot(x, g);
hold off;
grid on;
xlim([min(x), max(x)]);
title("Histogram of y1");

subplot(4, 2, 4);
hold on;
histogram(y2, x, 'Normalization', 'pdf');
xticks(x);
xticklabels(x);
plot(x, g);
hold off;
grid on;
xlim([min(x), max(x)]);
title("Histogram of y2");

subplot(4, 2, 5);
plot(f1, Y1);
grid on;
xlabel("f");
ylabel("|Y1(w)|");
title("FFT of y1");

subplot(4, 2, 6);
plot(f2, Y2);
grid on;
xlabel("f");
ylabel("|Y2(w)|");
title("FFT of y1");

subplot(4, 2, 7);
plot(lags1, r1);
grid on;
xlim([min(lags1), max(lags1)]);
xlabel("f");
ylabel("|Y1(w)|");
title("Autocorrelation of y1");

subplot(4, 2, 8);
plot(lags2, r2);
grid on;
xlim([min(lags2), max(lags2)]);
xlabel("f");
ylabel("|Y2(w)|");
title("Autocorrelation of y2");

%save
saveas(gcf, "Q2_1.png");