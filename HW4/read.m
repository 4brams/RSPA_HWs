%init
close all;
clear all;
clc;

%config
FILE_ROUTE = "./data/raw.chungli52.66a1bb88";

%read
infile = FILE_ROUTE;
%outfile = "./swa_head.dat";
%outfile2 = "./sea_time.dat";
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
IQ(1 : m, 1 : n) = fread(fid, [m, n], "short");

I = zeros(nch / 2, nrg1 - 1, ntrial);
Q = zeros(nch / 2, nrg1 - 1, ntrial);

T = 0 : t_res : t_res * (ntrial - 1);
R = range(2 : end);

for channel = 1 : (nch / 2)
    for r = 2 : nrg1
        I(channel, r - 1, :) = IQ((nch * (r - 1)) + (2 * channel) - 1, :);
        Q(channel, r - 1, :) = IQ((nch * (r - 1)) + (2 * channel), :);
    end
end

save Radar_sea I Q T R nch nrg1 ntrial range t_res ipp cint;