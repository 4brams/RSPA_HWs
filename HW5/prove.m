%init
close all;
clear all;
clc;

%func
function y = sinc(x)
    y = sin(pi*x) ./ (pi*x);
    y(x == 0) = 1; 
end

%config
A = 1;
f0 = 50;
fr = 10;

%val
f_sampling = 10 * f0;
dt = 1 / f_sampling;
t = -0.25 : dt : 0.25 -dt;
x = 2 * pi * fr * t;
N = length(t);
f = (-N / 2 : N / 2 - 1) / N * f_sampling;

%signal
y_predict = 2 * sin(x);
r = A * cos(2 * pi * f0 * t + x);
os = 2 * sin(2 * pi * f0 * t);
e = r .* os;

fc = 30;
lpf = sinc(2.2 * fc * t) / N * f_sampling / 10;
y_real = conv(e, lpf, "same");

R = abs(fftshift(fft(r))) / N;
OS = abs(fftshift(fft(os))) / N;
E = abs(fftshift(fft(e))) / N;
LPF = abs(fftshift(fft(lpf))) / N;
Y_REAL = abs(fftshift(fft(y_real))) / N;
Y_PREDICT = abs(fftshift(fft(y_real))) / N;

%plot
fig1 = figure;

subplot(4, 2, 1);
plot(t, r);
grid on;
xlabel("t(s)");
ylabel("Amp.(V)");
title("s(t)");

subplot(4, 2, 2);
plot(f, R);
grid on;
xlabel("f(Hz)");
ylabel("|S(f)|");
title("s Spectrum");

subplot(4, 2, 3);
plot(t, e);
grid on;
xlim([-0.25, 0.25]);
xlabel("t(s)");
ylabel("Amp.(V)");
title("e(t)");

subplot(4, 2, 4);
hold on;
plot(f, E);
plot(f, LPF / max(LPF) * max(E) * 0.8);
hold off;
grid on;
xlabel("f(Hz)");
ylabel("|E(f)|");
title("e Spectrum");

subplot(4, 2, 5);
plot(t, y_real);
grid on;
xlim([-0.25, 0.25]);
xlabel("t(s)");
ylabel("Amp.(V)");
title("y_{real}(t)");

subplot(4, 2, 6);
plot(f, Y_REAL);
grid on;
xlabel("f(Hz)");
ylabel("|Y_{REAL}(f)|");
title("y_{real} Spectrum");

subplot(4, 2, 7);
plot(t, y_predict);
grid on;
xlim([-0.25, 0.25]);
xlabel("t(s)");
ylabel("Amp.(V)");
title("y_{wrong}(t)");

subplot(4, 2, 8);
plot(f, Y_PREDICT);
grid on;
xlabel("f(Hz)");
ylabel("|Y_{WRONG}(f)|");
title("y_{wrong} Spectrum");

saveas(fig1, "prove.png");