%init
close all;
clear all;
clc;

%config
A = 5;
T = 10;
N = 2 ^ 11;
dt = 1e-2;

%func
function y = f1(t, T, A)
    y = zeros(1, length(t));

    for i = 1:length(t)
        if((t(i) < (T / 2)) && (t(i) > (- T / 2)))
            y(i) = A;
        else
            y(i) = 0;
        end
    end
end

function y = f2(t, T, A)
    y = zeros(1, length(t));
    
    for i = 1:length(t)
        if((t(i) <= (T / 4)) && (t(i) >= (- T / 4)))
            y(i) = A;
        elseif((t(i) <= (T / 2)) && (t(i) > (T / 4)))
            y(i) = (- A / (T / 4)) * (t(i) - T / 2);
        elseif((t(i) < (- T / 4)) && (t(i) >= (- T / 2)))
            y(i) = (A / (T / 4)) * (t(i) + T / 2);
        else
            y(i) = 0;
        end
    end
end

function y = f3(t, T, A)
    for i = 1:length(t)
        if((t(i) <= (T / 2)) && (t(i) >= 0))
            y(i) = (-A / (T / 2)) * (t(i) - T / 2);
        elseif((t(i) < (0)) && (t(i) >= (- T / 2)))
            y(i) = (A / (T / 2)) * (t(i) + T / 2);
        else
            y(i) = 0;
        end
    end
end

%x
t = (- N / 2 : N / 2 - 1) * dt;
f = (- N / 2 : N / 2 - 1) / (N * dt);

%y
y1 = f1(t, T, A);
y2 = f2(t, T, A);
y3 = f3(t, T, A);

Y1 = abs(fftshift(fft(y1, N))) / N;
Y2 = abs(fftshift(fft(y2, N))) / N;
Y3 = abs(fftshift(fft(y3, N))) / N;

%plot
figure;

subplot(3, 2, 1);
plot(t, y1);
grid on;
xlim([min(t), max(t)]);
ylim([min(y1) - 1, max(y1) + 1]);
xlabel("t");
ylabel("y_1");
title("f_1(t)");

subplot(3, 2, 2);
plot(f, Y1);
grid on;
xlim([-2, 2]);
ylim([min(Y1) - 0.2, max(Y1)]);
xlabel("f");
ylabel("|Y_1|");
title("f_1(w)");

subplot(3, 2, 3);
plot(t, y2);
grid on;
xlim([min(t), max(t)]);
ylim([min(y2) - 1, max(y2) + 1]);
xlabel("t");
ylabel("y_2");
title("f_2(t)");

subplot(3, 2, 4);
plot(f, Y2);
grid on;
xlim([-2, 2]);
ylim([min(Y2) - 0.2, max(Y2)]);
xlabel("f");
ylabel("|Y_2|");
title("f_2(w)");

subplot(3, 2, 5);
plot(t, y3);
grid on;
xlim([min(t), max(t)]);
ylim([min(y3) - 1, max(y3) + 1]);
xlabel("t");
ylabel("y_3");
title("f_3(t)");

subplot(3, 2, 6);
plot(f, Y3);
grid on;
xlim([-2, 2]);
ylim([min(Y3) - 0.2, max(Y3)]);
xlabel("f");
ylabel("|Y_3|");
title("f_3(w)");

%save
saveas(gcf, "Q1_1.png");