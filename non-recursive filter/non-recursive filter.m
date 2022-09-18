clear all; 
clc; 
close all; 
format long;

N = 2048; 
t = 1:N; 
L = 256; 
tf = 257:512; 
Fs = 2048; 
f = [250 500 300 350 600 650]; 
A = [1.0, 0.2, 0.5, 0.3, 0.2, 0.8]; 
indx1 = 1; 
indx2 = 4;

fni = f ./ Fs; 

y = A(1) * cos(2 * pi * fni(1) .* t) + A(2) * sin(2 * pi * fni(2) .*t) + A( ...
    3) * -sin(2 * pi * fni(3) .*t) + A(4) * cos(2 * pi * fni(4) .*t) + A(5 ...
    ) * sin(2 * pi * fni(5) .*t) + A(6) * -sin(2 * pi * fni(6) .*t);

figure
plot(tf, y(tf), 'b');
title('График входного сигнала', FontSize=14);
axis tight;
grid on;

Xs = fft(y, N); 
MXs = (2 / N) * abs(Xs); 
figure
plot(t-1, MXs(t), 'r');
title('График модуля спектра входного сигнала', FontSize=14);
axis tight;
grid on;

f1s = fni(indx1); 
f2s = fni(indx2);
Wn = [1.8*f1s 2.2*f1s 1.8*f2s 2.2*f2s]; 
b = fir1(L, Wn, 'stop');

yf = filter(b, 1, y);
XF = fft(yf, N); 
MXF = (2 / N) * abs(XF); 
figure
plot(t-1, MXF(t), 'k');
title('График модуля спектра отфильтрованного сигнала', FontSize=14);
axis tight;
grid on;

figure
freqz(b, 1);
title('Графики АЧХ и ФЧХ', FontSize=14);
axis tight;
grid on;

fprintf(['Минимальный порядок фильтра L, при котором амплитуда удаленных ' ...
    'компонент сигнала не будет превышать значение 0,01: f(%d) = %d и f(' ...
    '%d) = %d равен: L = %d\n'], indx1, f(indx1), indx2, f(indx2), L);
fprintf('Вектор нормированных частот среза: [%f %f %f %f]\n', Wn);