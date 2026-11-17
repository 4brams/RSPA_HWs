%init
close all;
clear all;
clc;

%config
channel = 1;
FILE_ROUTES = ...
["./data/raw.chungli52.66a14a10",...
"./data/raw.chungli52.66a15f56",...
"./data/raw.chungli52.66a1583f",...
"./data/raw.chungli52.66a1666e",...
"./data/raw.chungli52.66a15127"];

%val
infile = FILE_ROUTES(1);
ghead_52_radar1;
prf = str2num(s19(5 : end));
ACQ_Time = ns2;
c = 3 * (10 ^ 8);
r_res = pw * c / 2 / 1000;
t_res = ipp * cint;
f_Ny = 1 / (2 * t_res);
max_range = min_range + r_res*(nrg1 - 1);
range = min_range : r_res : max_range;
lambda = c / (52 * (10 ^ 6));
max_Doppler = f_Ny * lambda / 2;
m = nch * nrg1;
n = ntrial;

T = 0 : t_res : t_res * (ntrial - 1);
R = range(2 : end);
IQ = zeros(length(FILE_ROUTES) ,nrg1 - 1, ntrial);
start_time = strings(1, length(FILE_ROUTES));

%read
for i = 1 : length(FILE_ROUTES)
    infile = FILE_ROUTES(i);
    ghead_52_radar1;
    
    IQx(1 : m, 1 : n) = fread(fid, [m, n], "short");

    I = zeros(nrg1 - 1, ntrial);
    Q = zeros(nrg1 - 1, ntrial);
    
    for r = 2 : nrg1
        I(r - 1, :) = IQx((nch * (r - 1)) + (2 * channel) - 1, :);
        Q(r - 1, :) = IQx((nch * (r - 1)) + (2 * channel), :);
    end
    
    IQ(i, :, :) = I + 1i * Q;
    start_time(i) = year_date_time;
end

save Radar_Data IQ T R nch nrg1 ntrial t_res ipp cint start_time;