clc
clear all
[error_bit_all(:,1),error_symbol_all(:,1)]=rs_qpsk_coding();
[error_bit_all(:,2),error_symbol_all(:,2)]=rs_16qam_coding();
% [error_bit_all(:,3)]=qpsk_predict();
figure(2);
semilogy(0:20,error_bit_all(:,1),'r-^');
hold on
semilogy(0:20,error_bit_all(:,2),'r-s');
% hold on
% semilogy(0:20,error_bit_all(:,3),'g-s');
title('Comparison of qpsk coding and 16qam with Rs coding performance in rayligh channel');
legend('qpsk','16qam');
xlabel('SNR');
ylabel('BER');