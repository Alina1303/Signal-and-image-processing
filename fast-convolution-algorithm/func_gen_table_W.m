function[W_Re, W_Im] = func_gen_table_W(N)
    two_pi_N = 2 * pi / N;    
    NW = N / 2;  % NW - Количество весовых коэффициентов в таблице
    n = 1 : NW;    
    arg = two_pi_N * (n - 1);
    W_Re(n) = cos(arg);     % Реальная часть    
    W_Im(n) = -sin(arg);    % Мнимая часть    
end
