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

%% 