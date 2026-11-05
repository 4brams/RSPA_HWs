%init
close all;
clear all;
clc;

%config
C = 1000;
N = floor(C ^ (1 / 2));

%even distribution
even_dis = rand(1, C);
mu = mean(even_dis);
sigma = std(even_dis);

fig1 = figure;
hold on;
h = histogram(even_dis, N, 'Normalization', 'pdf');
filter = [1/2, 1/2];
x = conv(h.BinEdges, filter, "valid");
e = ones(1, length(x));
plot(x, e);
hold off;
grid on;
xlabel("Number");
ylabel("Probability Density");
sgtitle("Even Distrbution, \mu = " + mu + ", \sigma = " + sigma + ", RMSE = " + rmse(h.Values, e));

%gauss distribution
gauss_dis = randn(1, C);
mu = mean(gauss_dis);
sigma = std(gauss_dis);

fig2 = figure;
hold on;
h = histogram(gauss_dis, N, 'Normalization', 'pdf');
filter = [1/2, 1/2];
x = conv(h.BinEdges, filter, "valid");
g = (1 / (sigma * ((2 * pi) ^ (1 / 2)))) * exp(-((x - mu) .^ 2)/(2 * (sigma ^ 2)));
plot(x, g);
hold off;
grid on;
xlabel("Number");
ylabel("Probability Density");
sgtitle("Gauss Distrbution, \mu = " + mu + ", \sigma = " + sigma + ", RMSE = " + rmse(h.Values, g));

%save
saveas(fig1, "even_distrbution.png");
saveas(fig2, "gauss_distrbution.png");
