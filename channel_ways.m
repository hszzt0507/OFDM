clc
clear all
[error_bit_all(:,1),error_symbol_all(:,1)]=rs_reiligh_coding();
[error_bit_all(:,2),error_symbol_all(:,2)]=rs_awgn_coding();
[error_bit_all(:,3)]=qpsk_predict();
figure(2);
semilogy(0:20,error_bit_all(:,1),'r-^');
hold on
semilogy(0:20,error_bit_all(:,2),'b-*');
hold on
semilogy(0:20,error_bit_all(:,3),'g-s');
title('Comparison of religh channel and AWGN channel with (15.11)RS coding and qpsk performance');
legend('reiligh','awgn','predict');
xlabel('SNR');
ylabel('BER');