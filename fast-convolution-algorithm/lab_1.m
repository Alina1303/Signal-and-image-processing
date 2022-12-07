close all; 
clc; 
clear all; 
format long;

t_beg_model_time = tic; %старт времени расчёта модели

%Исходные данные
Q_COUNTS = 4096; %число отсчётов сигнала
N_it = log2(Q_COUNTS); %число итераций БПФ
q = 1:Q_COUNTS; %счётчик данных
fprintf('Количество отсчётов БПФ = %d\n', Q_COUNTS);
fprintf('Количество итераций БПФ = %d\n', N_it);
    
%1 Генерация сигнала
te = 0:1/4095:1;
q = 1:4096;
E0 = 1600 * exp(70 * 1i * te)- 1458 * exp(-20 * te * 1i * pi) - 2102 * cos(3 ...
    * 2 * pi * te) + 2 * (750 * sin(2 * pi * te) - 500 * sin(3 * pi * te));


Et_Im = imag(E0); %выделение реальной и мнимой части
Et_Re = real(E0);

%2 Построение графика входного сигнала
figure    
subplot(2, 1, 1);        
plot(q, Et_Re(q), 'b');
title('График входного сигнала (Re)', FontSize=14);
axis tight;
grid on;                
subplot(2, 1, 2);        
plot(q, Et_Im(q), 'b');
title ('График входного сигнала (Im)', FontSize=14);
axis tight;
grid on; 

%3 Генерация шума для сигнала
Max_S = max(max(abs(Et_Re)), max(abs(Et_Im))); %расчёт максимального значения массива
rng(64801, 'twister');
N_por = round(Max_S); %порог шума
N0 = randi([-N_por N_por], 1, Q_COUNTS); %генерация данных шума

%4 Формирование зашумленного входного сигнала
S_Re(q) = Et_Re(q) + N0(q); %формирование сигнала (сигнал + шум)
S_Im(q) = Et_Im(q) + N0(q);

%5 Построение графика зашумленного входного сигнала
figure    
subplot(2, 1, 1);        
plot(q, S_Re(q), 'b');
title('График входного сигнала + шума (Re)', FontSize=14);
axis tight;
grid on;                
subplot(2, 1, 2);        
plot(q, S_Im(q), 'b');
title ('График входного сигнала + шума ( Im )', FontSize=14); 
axis tight;
grid on;             

%6 Алгоритм быстрой свёртки (БПФ)

[W_Re, W_Im] = func_gen_table_W(Q_COUNTS); %расчёт поворачивающих коэффициентов

%7 Запись значения поворачивающих множителей
write = sprintf('W_%d.txt', Q_COUNTS);
write_complex_to_dec(write, W_Re, W_Im);

W_Re = bitrevorder( W_Re ); %выстраивание комплексных коэффициентов W в 
% бит-реверсивном порядке
W_Im = bitrevorder( W_Im );

S_fft_Et(q) = 0; % инициализация массива произведений спектра сигнала на 
% спектр эталонного сигнала

flag_fft = -1; % флаг расчёта БПФ = 1 и ОБПФ = -1

%8 Процедура расчёта алгоритма прямого и обратного БПФ для входного сигнала
for step = 1:2 % step == 1 Расчёт БПФ, step == 2 Расчёт ОБПФ
   flaf_fft = -flag_fft; % флаг расчёта БПФ (flag_fft == 1) или ОБПФ (flag_fft == -1) для функции
   
   S_Re = bitrevorder(S_Re); % задание бит-реверсивного порядка при расчёте БПФ (ОБПФ) для входных данных
   S_Im = bitrevorder(S_Im);
  
   % РасчётБПФ - ОБПФ   
   for it = 1 : N_it   % число итераций БПФ
      count_W = 2^(it - 1);  % кол-во повторяющихся коэфф. W при расчёте итерации
      range_W = 2^(N_it - it); % диапазон используемых коэффициентов W
      for k = 1 : count_W     
         from_W = 1 + range_W * (k - 1);   % начало диапазона
         to_W = range_W * k;                 % конец диапазона

         if (it <= 3) %выведение значений переменных
            fprintf(['Iter %d, count_W = %d, range_W = %d, from_W = %d, ' ...
                'to_W = %d\n'], it, count_W, range_W, from_W, to_W);
         end
         % Подбор коэффициентов W для текущей итерации
         temp_W_Re(from_W : to_W) = W_Re(1 : range_W);
         temp_W_Im(from_W : to_W) = W_Im(1 : range_W);                           
      end       
 
%9 Расчёт базовых операций (бабочка) алгоритма БПФ - ОБПФ для текущей итерации 
      [S_Re, S_Im] = func_fft_float(S_Re, S_Im, temp_W_Re, temp_W_Im, flag_fft);        
   end
%10 Процедура согласованной фильтрации входного сигнала с помощью БПФ и
%ОБПФ
   if step == 1
      % Комплексное умножение результатов БПФ на значения спектра эталонного сигнала
      S_fft_Et(q) = (S_Re + 1i * S_Im(q)).*(Et_Re(q) + 1i * Et_Im(q));
      S_Re(q) = real(S_fft_Et(q));  % Реальная часть
      S_Im(q) = imag(S_fft_Et(q));  % Мнимая часть
   else % if step == 2
      % Масштабирование результатов ОБПФ - умножение на коэффициент 1/Q_COUNTS в соответствии с формулой
      S_Re = S_Re./Q_COUNTS;
      S_Im = S_Im./Q_COUNTS;
   end
end

%11 Построение графика абсолютных значений полученного сигнала
abs_S_ifft(q) = abs(S_Re(q) + 1i * S_Im(q));
figure    
plot(q, abs_S_ifft(q), 'b');
title ('Абсолютные значения ОБПФ');
axis tight;
grid on;   

%12 Результаты алгоритма быстрой свёртки
conv = sprintf('conv_%d.txt', Q_COUNTS);
write_complex_to_dec(conv, S_Re, S_Im );
          
t_end_model_time = toc(t_beg_model_time); %окончание времени расчёта
fprintf( ['Модель рассчитана!( %s ) Время расчёта: %f сек. (%f мин.)\n' ...
    ''], datestr( now ), t_end_model_time, t_end_model_time / 60 );

%Создать массив-вектор от 1 до 8 с шагом 1
n = (1:8)'; %массив-вектор
z = bitrevorder(n); %функция бит-реверс
fprintf('Число = %d, реверс = %d\n', n, z);
