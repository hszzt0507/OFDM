clc
clear all
[error_bit_all(:,1),error_symbol_all(:,1)]=con_coding();
[error_bit_all(:,2),error_symbol_all(:,2)]=con_rs_coding();
[error_bit_all(:,3),error_symbol_all(:,3)]=rs_coding();
[error_bit_all(:,4),error_symbol_all(:,4)]=no_coding();
figure(2);
semilogy(1:20,error_bit_all(:,4),'r-^');
hold on
semilogy(1:20,error_bit_all(:,3),'r-s');
hold on
semilogy(1:20,error_bit_all(:,1),'g-s');
hold on
semilogy(1:20,error_bit_all(:,2),'g-^');
title('Comparison of RS coding and RS with Convolutional coding performance');
legend('no coding','RS coding','Convolutional coding£¨2£¬1£¬7£©','Convolutional coding RS');
xlabel('SNR');
ylabel('BER');