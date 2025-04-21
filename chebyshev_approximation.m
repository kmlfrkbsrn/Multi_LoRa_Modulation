clear all;
close all;
clc;

% Chebyshev Polinomları ile g_T(t) Yaklaşımı

% Parametreler
N =4;  % Chebyshev polinomunun derecesi
t = linspace(0, 2*pi, 500);  % t aralığı
f = @(t) sinc(t);  % g_T(t) fonksiyonu, örneğin sinc(t)

% Fonksiyonu Chebyshev polinomları ile yakınsama (Yaklaşık olarak)
% İlk olarak Chebyshev polinomları hesaplanır
T = zeros(N, length(t));  % Polinomları tutmak için matris
for n = 1:N
    T(n, :) = cos(n * t);  % Birinci Tür Chebyshev polinomu
end

% Fonksiyonun yaklaşık katsayılarını hesapla
cheb_coeffs = zeros(1, N);  % Katsayılar
for n = 1:N
    cheb_coeffs(n) = trapz(t, f(t) .* T(n, :));  % Katsayıyı hesapla (inner product)
end

% Yaklaşık fonksiyonu hesapla
g_approx = zeros(1, length(t));  % Yaklaşık fonksiyon
for n = 1:N
    g_approx = g_approx + cheb_coeffs(n) * T(n, :);  % Yaklaşık fonksiyonu oluştur
end

% Grafik
figure;
hold on;
plot(t, f(t), 'k', 'DisplayName', 'Orijinal g_T(t)');  % Orijinal fonksiyon
plot(t, g_approx, 'r', 'DisplayName', 'Yaklaşık g_T(t)');  % Yaklaşık fonksiyon
legend show;
title('Chebyshev Polinomları ile g_T(t) Yaklaşımı');
xlabel('t');
ylabel('g_T(t) ve Yaklaşık Fonksiyon');
