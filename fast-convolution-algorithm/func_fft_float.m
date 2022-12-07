function[T_Re, T_Im] = func_fft_float(S_Re, S_Im, W_Re, W_Im, flag_fft) 
    N = length(S_Re);
    i = 1 : 2 : N-1;
    k = 1 : N/2;
    T_Re(k)= S_Re(i)+ S_Re(i+1);
    T_Im(k)= S_Im(i)+ S_Im(i+1);
    Z_Re(k) = S_Re(i)- S_Re(i+1);
    Z_Im(k) = S_Im(i)- S_Im(i+1);
    T_Re(k+N/2)= Z_Re(k).*W_Re(k) - Z_Im(k).*W_Im(k)*flag_fft;
    T_Im(k+N/2)= Z_Re(k).*W_Im(k)*flag_fft + Z_Im(k).*W_Re(k);
end
