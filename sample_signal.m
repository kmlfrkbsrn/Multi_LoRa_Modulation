clear all;
close all;
clc;

% Zaman vektörü
t = linspace(-0.05, 0.05, 1000);  % saniye

% Parametreler
N = 4;
f_c = 1000 ;
f_m = 100; % Hz cinsinden, mesaj band genişliği
beta = 30 ;
theta = 2*pi*f_c*t+2*pi*f_m*t.^2+beta;           % Faz başlangıcı
phi = 2*pi*f_m*t;          % Açı farkı

% Sinc fonksiyonu (pi'li tanım)
sinc_func = @(x) sin(pi*x) ./ (pi*x);
%sinc_func(0) = 1; % tanımsızlık için düzeltme

% g_T(t): Pulse shaping zarfı
g_T = N * sinc_func(N * f_m * t) ./ sinc_func(f_m * t);

% Cos bileşeni
cos_component = cos(theta + (N-1)*phi/2) .* ones(size(t)); % sabit bir değer

% x(t): toplam sinyal
x_t = g_T .* cos_component;

% Grafikler
figure;

subplot(3,1,1);
plot(t, x_t, 'k', 'LineWidth', 1.5);
title('x(t) = g_T(t) * cos(\theta + (N-1)\phi/2)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, g_T, 'b', 'LineWidth', 1.5);
title('g_T(t) = N sinc(Nf_m t) / sinc(f_m t)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, cos_component, 'r', 'LineWidth', 1.5);
title('Cos(\theta + (N-1)\phi/2)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;


