function H = compute_H(h, N, w, m, dt)
    H = ones(size(w));
    for n = -m : m
        H = H + h(n + m + 1) * cos(2 * pi * f * n * dt);
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
    filterSystemResponse = compute_H(f, point+1, delta_t);

    [val, idx] = max(abs(sunspotSpectrum));
    sunPeriod = abs((1 / (f(idx))) / 12);
end