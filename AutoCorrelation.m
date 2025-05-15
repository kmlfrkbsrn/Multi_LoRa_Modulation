% Parameters
N = 8;
f_m = 1e3;
t = linspace(-25e-3, 25e-3, 10000);
dt = t(2) - t(1);

% Original signal
gt = sinc(N*f_m*t) ./ sinc(f_m*t);

% Autocorrelation via convolution
%R = conv(gt, fliplr(gt), 'same') * dt;
tau = 1:10000;

[lags, R] = circ_autocorr(gt);
plot(lags, real(R));
xlabel('Lag'); ylabel('Autocorrelation');
title('Circular Autocorrelation of a[n]');


zero_crossings = [];
for i = 1:length(R)-1
    if R(i)*R(i+1) < 0
        tau_zero = tau(i) - R(i)*(tau(i+1)-tau(i)) / (R(i+1)-R(i));
        zero_crossings = [zero_crossings, tau_zero];
    end
end

%% --- Plot: Original Signal ---
%figure;
%plot(t*1e3, gt, 'b', 'LineWidth', 1.5);
%xlabel('Time (ms)');
%ylabel('g_T(t)');
%title('Original Signal g_T(t)');
%grid on;
%
% --- Plot: Autocorrelation with Zero Crossings ---
figure;
plot(lags, R, 'r', 'LineWidth', 1.5);
hold on;
yline(0, '--k');
plot(zero_crossings, zeros(size(zero_crossings)), 'ko', 'MarkerFaceColor', 'g');
xlabel('Sample Delay');
ylabel('Autocorrelation R(\tau)');
title('Autocorrelation of g_T(t) with Zero Crossings');
legend('R(\tau)', 'Zero Line', 'Orthogonal Points', 'Location', 'best');
grid on;
xticks(0:25:10000);

figure;
plot(gt(1:1000));
hold on;
plot(gt(26:1025));
plot(gt(51:1050));
title('R(\tau) = 0');
legend("Original g_T", "25 sample Delayed g_T", "50 sample Delayed g_T");



function [lags, R] = circ_autocorr(a)
    M = length(a);
    a = a / norm(a);  % normalize et (isteğe bağlı)
    R = zeros(1, M);
    for k = 0:M-1
        a_shifted = circshift(a, [0, k]); % circular shift
        R(k+1) = sum(a .* conj(a_shifted));  % dot product
    end
    lags = 0:M-1;
end


