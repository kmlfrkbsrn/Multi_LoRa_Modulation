% Parameters
N = 8;
f_m = 1e3;
fs = 100e3;                % Sampling rate
T = 1 / f_m;               % Symbol duration (orthogonality spacing)
t = linspace(-5e-3, 5e-3, fs*10e-3); % 10 ms total duration
dt = t(2) - t(1);

% g_T(t) pulse
gt = sinc(N*f_m*t) ./ sinc(f_m*t);
gt(isnan(gt)) = N; % Fix divide-by-zero at t = 0

% Chirp parameters
BW = 10e3;
T_chirp = 2e-3;  % Chirp duration
k = BW / T_chirp;

% Time vector for chirps
tc = linspace(0, T_chirp, round(fs*T_chirp));

% Generate 2 chirps
chirp1 = cos(2*pi*(0.5*k*tc.^2));         % Upchirp
chirp2 = cos(2*pi*(0.5*k*tc.^2 + 2e3*tc)); % Upchirp + offset (user-specific)

% Embed chirps at different time positions using shifted g_T
signal = zeros(size(t));

% Find center indices
center1 = find(t >= -1.5e-3, 1);  % User 1 at -1.5ms
center2 = find(t >= 1.5e-3, 1);   % User 2 at +1.5ms

% Window pulses (g_T * chirp) and insert into main signal
g1 = sinc(N*f_m*(t - t(center1))) ./ sinc(f_m*(t - t(center1)));
g2 = sinc(N*f_m*(t - t(center2))) ./ sinc(f_m*(t - t(center2)));
g1(isnan(g1)) = N; g2(isnan(g2)) = N;

user1 = g1 .* interp1(tc + t(center1), chirp1, t, 'linear', 0);
user2 = g2 .* interp1(tc + t(center2), chirp2, t, 'linear', 0);

% Total signal
signal = user1 + user2;

% --- Plot Users and Sum ---
figure;
subplot(3,1,1); plot(t*1e3, user1); title('User 1 Signal'); xlabel('Time (ms)'); grid on;
subplot(3,1,2); plot(t*1e3, user2); title('User 2 Signal'); xlabel('Time (ms)'); grid on;
subplot(3,1,3); plot(t*1e3, signal); title('Combined Signal'); xlabel('Time (ms)'); grid on;

% --- Autocorrelation ---
R = conv(signal, fliplr(signal), 'same') * dt;
figure;
plot(t*1e3, R, 'r');
xlabel('Time Delay (ms)');
ylabel('Autocorrelation');
title('Autocorrelation of Combined Signal');
grid on;
