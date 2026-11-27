% test_fix.m
close all; clear all; clc;

NFFT = 2 ^ 10;
range_idxs = [30, 20, 10, 5];
f_radar = 52e6;
c = 3e8;

load Radar_Data.mat;

IQ_size = size(IQ);
file_num = IQ_size(1);
CI_times = floor(ntrial / NFFT);
F = (-NFFT / 2 : NFFT / 2 - 1) / NFFT / t_res;
V = -(c * F) / (2 * f_radar);
[V, sort_idx] = sort(V, 'ascend');
dv = V(2) - V(1);

IQ_FFT = fftshift(fft(IQ, NFFT, 3), 3) / NFFT;
IQ_FFT(:, :, NFFT / 2 + 1) = (IQ_FFT(:, :, NFFT / 2 + 2) + IQ_FFT(:, :, NFFT / 2)) / 2;
IQ_FFT = IQ_FFT(:, :, sort_idx);

IQ_FFT_CI = zeros(file_num, nrg1 - 1, NFFT);
for time = 1 : CI_times
    S = fftshift(fft(IQ(:, :, (time - 1) * NFFT + 1 : time * NFFT), NFFT, 3), 3) / NFFT;
    IQ_FFT_CI = IQ_FFT_CI + abs(S) .^ 2;
end
IQ_FFT_CI(:, :, NFFT / 2 + 1) = (IQ_FFT_CI(:, :, NFFT / 2 + 2) + IQ_FFT_CI(:, :, NFFT / 2)) / 2;
IQ_FFT_CI = IQ_FFT_CI / CI_times;
IQ_FFT_CI = IQ_FFT_CI(:, :, sort_idx);

Sxx = IQ_FFT_CI;
s = sort(Sxx, 3, 'ascend');
P_noise = mean(s(:, :, 1 : floor(0.3 * NFFT)), 3);

V_a_idxs = V < 3 & V > -3;
V_r_idxs = V < -3 & V > -15;

V_a_mask = false(file_num, nrg1 - 1, length(V));
V_r_mask = false(file_num, nrg1 - 1, length(V));

P_a = zeros(file_num, nrg1 - 1);
P_signal_a = zeros(file_num, nrg1 - 1);
SNR_a = zeros(file_num, nrg1 - 1);
va = nan(file_num, nrg1 - 1);
wa = nan(file_num, nrg1 - 1);

P_r = nan(file_num, nrg1 - 1);
P_signal_r = zeros(file_num, nrg1 - 1);
SNR_r = nan(file_num, nrg1 - 1);
vr = nan(file_num, nrg1 - 1);
wr = nan(file_num, nrg1 - 1);

for f = 1 : file_num
    for r = 1 : nrg1 - 1
        spec = squeeze(Sxx(f, r, :))';
        pn = P_noise(f, r);
        sclean = max(spec - pn, 0);

        % 搜尋大氣峰 (在 -2.5 ~ 2.5 m/s 內)
        spec_a = spec;
        spec_a(~V_a_idxs) = -Inf;
        [a_max, a_idx] = max(spec_a);

        % 搜尋降雨峰 (在 -14 ~ -2 m/s 內)
        spec_r = spec;
        spec_r(~V_r_idxs) = -Inf;
        [r_max, r_idx] = max(spec_r);

        % 判斷降雨是否存在
        has_rain = false;
        if R(r) <= 5.8 && r_max > pn * 1.05
            [val_min, min_loc] = min(spec(r_idx : a_idx));
            m_border = r_idx + min_loc - 1;
            if val_min < 0.9 * r_max && val_min < 0.9 * a_max
                has_rain = true;
            end
        end

        if has_rain
            % 左邊界 l_border (向左延伸到雜訊水平或 -15 m/s)
            l_border = find(V <= V(r_idx) & V >= -15 & spec <= pn * 1.1, 1, 'last');
            if isempty(l_border)
                l_border = find(V >= -15, 1, 'first');
            end

            % 右邊界 r_border (向右延伸到雜訊水平或 +3.5 m/s)
            r_border = find(V >= V(a_idx) & V <= 3.5 & spec <= pn * 1.1, 1, 'first');
            if isempty(r_border)
                r_border = find(V <= 3.5, 1, 'last');
            end

            V_r_mask(f, r, l_border : m_border) = true;
            V_a_mask(f, r, m_border + 1 : r_border) = true;
        else
            l_border = find(V <= V(a_idx) & V >= -3.5 & spec <= pn * 1.1, 1, 'last');
            if isempty(l_border)
                l_border = find(V >= -3.5, 1, 'first');
            end
            r_border = find(V >= V(a_idx) & V <= 3.5 & spec <= pn * 1.1, 1, 'first');
            if isempty(r_border)
                r_border = find(V <= 3.5, 1, 'last');
            end
            V_a_mask(f, r, l_border : r_border) = true;
        end

        % atm
        mask_a = squeeze(V_a_mask(f, r, :))';
        sa = sclean .* mask_a;
        sum_sa = sum(sa);
        if sum_sa > 0
            P_signal_a(f, r) = sum_sa / NFFT;
            P_a(f, r)        = P_signal_a(f, r) + pn;
            SNR_a(f, r)      = P_signal_a(f, r) / pn;
            va(f, r)         = sum(sa .* V) / sum_sa;
            wa(f, r)         = sqrt(sum(sa .* (V - va(f, r)).^2) / sum_sa);
        end

        % rain
        if has_rain
            mask_r = squeeze(V_r_mask(f, r, :))';
            sr = sclean .* mask_r;
            sum_sr = sum(sr);
            if sum_sr > 0
                P_signal_r(f, r) = sum_sr / NFFT;
                P_r(f, r)        = P_signal_r(f, r) + pn;
                SNR_r(f, r)      = P_signal_r(f, r) / pn;
                vr(f, r)         = sum(sr .* V) / sum_sr;
                wr(f, r)         = sqrt(sum(sr .* (V - vr(f, r)).^2) / sum_sr);
            end
        end
    end
end

N_a = sum(V_a_mask, 3);
N_r = sum(V_r_mask, 3);
vt = vr - va;

disp('Test finished successfully without any error!');
exit;
