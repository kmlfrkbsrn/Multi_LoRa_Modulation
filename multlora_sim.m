%% 
%g_T(t) function
N = 8;
f_m = 1e3;
t = linspace(-5e-3, 5e-3, 1000);
gt = N * sinc(N*f_m*t) ./ sinc(f_m*t);
plot(t*1e3, gt);
xlabel('Time (ms)'); ylabel('g_T(t)');
title('Approximation of Pulse Shape g_T(t)');

hold on;

%D_N(t) function :Drichlect Kernel
N = 8;
f_m = 1e3;
t = linspace(-5e-3, 5e-3, 1000);
gt = sinc(N*f_m*t) ./ sinc(f_m*t);
plot(t*1e3, gt);
xlabel('Time (ms)'); ylabel('g_T(t)');
title('Approximation of Pulse Shape D_T(t)');
grid on;

%% TAYLOR SERIES APPROXIMATION

% Parametreler
N = 8;
f_m = 1e3; % Hz
t = linspace(-5e-3, 5e-3, 1000); % zaman (saniye)

% Orijinal g_T(t)
gt = N * sinc(N * f_m * t) ./ sinc(f_m * t);

% Taylor Açılımı (2. derece)
pi2 = pi^2;
num = 1 - (pi2 * (N * f_m * t).^2) / 6;
den = 1 - (pi2 * (f_m * t).^2) / 6;
gt_taylor = N * num ./ den;

% Grafik çizimi
figure;
plot(t * 1e3, gt, 'b', 'LineWidth', 2); hold on;
%plot(t * 1e3, gt_taylor, '--r', 'LineWidth', 2);
xlabel('Zaman (ms)');
ylabel('Genlik');
legend('g_T(t) (orijinal)', 'g_T(t) (Taylor approx)');
title('g_T(t) ve Taylor Yaklaşımı');
grid on;

%% TAYLOR SERIES APPROXIMATION -2 

% Parametreler
N = 8;
f_m = 1000; % Hz
t = linspace(-50e-3, 50e-3, 5000); % zaman (saniye)

% Orijinal g_T(t)
gt = N * sinc(N * f_m * t) ./ sinc(f_m * t);

% Taylor Açılımı (sadece küçük t için)
pi2 = pi^2;
num = 1 - (pi2 * (N * f_m * t).^2) / 6;
den = 1 - (pi2 * (f_m * t).^2) / 6;
gt_taylor = N * num ./ den;

% Grafik çizimi
figure;
plot(t * 1e3, gt, 'b', 'LineWidth', 2); hold on;

% Sadece küçük t aralığında Taylor çiz
idx = find(abs(t) <= 0.50e-3); % -0.5ms < t < 0.5ms
plot(t(idx) * 1e3, gt_taylor(idx), '--r', 'LineWidth', 2);

xlabel('Zaman (ms)');
ylabel('Genlik');
legend('g_T(t) (orijinal)', 'g_T(t) (Taylor approx)');
title('g_T(t) ve Taylor Yaklaşımı (Doğru Ölçek)');
grid on;

%% FOURIER APPROX. 

% Parametreler
N = 8;              % Kullanıcı sayısı
f_m = 1000;         % Hz
M = 20;             % Fourier terim sayısı 
t = linspace(-5e-3, 5e-3, 5000);  % zaman vektörü (saniye)

% Gerçek g_T(t) fonksiyonu
gt = N * sinc(N * f_m * t) ./ sinc(f_m * t);

% Fourier serisi yaklaşımı (geliştirilmiş)
fourier_approx = zeros(size(t));
for m = -M:M
    fourier_approx = fourier_approx + exp(1j * 2 * pi * m * f_m * t);
end
fourier_approx = real(fourier_approx) / (2*M + 1); % normalize et

% Grafik çizimi
figure;
plot(t*1e3, gt, 'b', 'LineWidth', 2); hold on;
plot(t*1e3, fourier_approx, '--r', 'LineWidth', 2);
xlabel('Zaman (ms)');
ylabel('Genlik');
title('g_T(t) ve Fourier Serisi Yaklaşımı ');
legend('Gerçek g_T(t)', 'Fourier Yaklaşımı');
grid on;

%% WINDOWING APPROX. 

% Parametreler
N = 8;              % Kullanıcı sayısı
f_m = 1000;         % Hz
t = linspace(-5e-3, 5e-3, 5000);  % zaman vektörü

% Gerçek g_T(t) fonksiyonu
gt = N * sinc(N * f_m * t) ./ sinc(f_m * t);

% Sinc fonksiyonu
sinc_func = sinc(f_m * t);

% Pencere fonksiyonları
hann_win = hann(length(t))';            % Hann penceresi
hamming_win = hamming(length(t))';       % Hamming penceresi
gaussian_win = gausswin(length(t), 2.5)';% Gaussian penceresi (spread ayarlı)

% Pencereli g_T(t) yaklaşımları
gt_hann = sinc_func .* hann_win;
gt_hamming = sinc_func .* hamming_win;
gt_gaussian = sinc_func .* gaussian_win;

% Grafik çizimi
figure;
plot(t*1e3, gt, 'k', 'LineWidth', 2); hold on;   % Gerçek g_T(t) siyah
plot(t*1e3, gt_hann, '--b', 'LineWidth', 2);      % Hann pencere mavi kesik
plot(t*1e3, gt_hamming, '--r', 'LineWidth', 2);   % Hamming pencere kırmızı kesik
plot(t*1e3, gt_gaussian, '--g', 'LineWidth', 2);  % Gaussian pencere yeşil kesik
xlabel('Zaman (ms)');
ylabel('Genlik');
title('Gerçek g_T(t) ve Pencere Yaklaşımları');
legend('Gerçek g_T(t)', 'Hann Pencereli', 'Hamming Pencereli', 'Gaussian Pencereli');
grid on;

