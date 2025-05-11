% Parameters
N = 8;
f_m = 1e3;
t = linspace(-5e-3, 5e-3, 1000);
dt = t(2) - t(1);

% Original signal
gt = sinc(N*f_m*t) ./ sinc(f_m*t);

% Autocorrelation via convolution
R = conv(gt, fliplr(gt), 'same') * dt;
tau = linspace(-5e-3, 5e-3, length(R));

zero_crossings = [];
for i = 1:length(R)-1
    if R(i)*R(i+1) < 0
        tau_zero = tau(i) - R(i)*(tau(i+1)-tau(i)) / (R(i+1)-R(i));
        zero_crossings = [zero_crossings, tau_zero];
    end
end

% --- Plot: Original Signal ---
figure;
plot(t*1e3, gt, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('g_T(t)');
title('Original Signal g_T(t)');
grid on;

% --- Plot: Autocorrelation with Zero Crossings ---
figure;
plot(tau*1e3, R, 'r', 'LineWidth', 1.5);
hold on;
yline(0, '--k');
plot(zero_crossings*1e3, zeros(size(zero_crossings)), 'ko', 'MarkerFaceColor', 'g');
xlabel('Time Delay \tau (ms)');
ylabel('Autocorrelation R(\tau)');
title('Autocorrelation of g_T(t) with Zero Crossings');
legend('R(\tau)', 'Zero Line', 'Orthogonal Points', 'Location', 'best');
grid on;
