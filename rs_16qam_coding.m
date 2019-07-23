function [error_bit_all,error_symbol_all]=rs_16qam_coding()
sonCarrierNum_temp=176;
symbols_Per_Carrier=1000;%每子载波含符号数/帧数
bits_Per_Symbol=1;%每符号含比特数
modulate_bit=4;%调制阶数(每个符号比特数)
IFFT_bin_length=2^ceil(log2(sonCarrierNum_temp));%FFT点数
PrefixRatio=1/4;%保护间隔与OFDM数据的比例 1/6~1/4
pilot_Inter=1;%插入导频间隔
% CP=PrefixRatio*IFFT_bin_length ;%每一个OFDM符号添加的循环前缀长度为1/4*IFFT_bin_length
CP=25;
SNR=0:1:20; %信噪比dB
nn=15;
kk=11;
%-------------------------------信源输入----------------------------------------------
inforSource=randi([0,1],1,sonCarrierNum_temp*symbols_Per_Carrier*bits_Per_Symbol);
%---------------------------信道编码--------------------------------------------------
msg4_temp=reshape(inforSource,4,[])';
msg4=bi2de(msg4_temp,'left-msb');%将原来的数据转换为4位16进制
msg4_togf=reshape(msg4,kk,[]).'; %带转换的矩阵，十一输入
msgGF=gf(msg4_togf,4);%转换为伽罗华域
msgrs=rsenc(msgGF,nn,kk); %(15,11）RS编码 11个输入 15个输出
msgrs1=reshape(msgrs.',1,length(msg4)/kk*nn);%将rs编码输出转成一行
msgrs2=de2bi(double(msgrs1.x),'left-msb');%十进制转二进制
source_coded_data_rs=reshape(msgrs2',1,length(msg4)/kk*nn*4);%待调制信号 输出一行信号（数据）
%----------------------------调制-----------------------------------------------------
data_temp1= reshape(source_coded_data_rs,modulate_bit,[])';   %以每组2比特进行分组，输出两列数据
modulate_data=qammod(bi2de(data_temp1),2^(modulate_bit));%输出一列数据
%-------------------------插入导频----------------------------------------------
modulate_data=reshape(modulate_data,60,[]);
[modulate_wide,modulate_length]=size(modulate_data);
modulate_data_temp=[modulate_data(1:30,:);zeros(1,modulate_length);modulate_data(31:60,:)];%在原来输出数据的中间插0
h1=commsrc.pn('GenPoly', [1 0 0 0 0 1 1],'NumBitsOut',61*modulate_length,'InitialConditions',[0 0 0 0 0 1]);
pn_code_temp=generate(h1);
pn_code=2*pn_code_temp-1;
pn_code=reshape(pn_code,61,[]);
modulate_data_pn=zeros(61,2*modulate_length);
for i=1:modulate_length
    modulate_data_pn(:,(2*i-1))=modulate_data_temp(:,i);
    modulate_data_pn(:,2*i)=pn_code(:,i);
end
modulate_data_pn(62:64,:)=0;
modulate_data_pn_out=[modulate_data_pn(31:64,:);modulate_data_pn(1:30,:)];
%-----------------------------------ifft-------------------------------------
time_signal_ifft=ifft(modulate_data_pn_out);
%-----------------------------------cp---------------------------------------
time_signal_cp=[time_signal_ifft(39:64,:);time_signal_ifft(1:64,:)];%把ifft的末尾CP个数补充到最前面
[time_signal_cp_wide,time_signal_cp_length]=size(time_signal_cp);
%-------------------------------------并串变换-------------------------------------
for ii=1:modulate_length
    time_signal_out(:,ii)=[time_signal_cp(:,2*ii-1);time_signal_cp(:,2*ii)];
end
time_signal_out_1=reshape(time_signal_out,[],1);
for a=1:21
% -----------------------------------信道----------------------------------------
% chan=comm.RayleighChannel('SampleRate',550000, ...
%     'PathDelays',[0 2e-6],'AveragePathGains',[0 -3],'MaximumDopplerShift',100,'RandomStream','mt19937ar with seed','Seed',8007);
chan=comm.RayleighChannel('SampleRate',550000, ...
    'PathDelays',[0 2e-6],'AveragePathGains',[0 -3],'MaximumDopplerShift',100);
Rayleigh_signal=chan(time_signal_out_1);
awgn_signal=awgn( Rayleigh_signal,a-1,'measured');%添加高斯白噪声
%   awgn_signal=time_signal_out_1;
%------------------------------------串并转换----------------------------------
   receive_signal_serial=awgn_signal;
   receive_signal_perallel=reshape(receive_signal_serial,time_signal_cp_wide,[]);
%------------------------------------去循环前缀---------------------------------
   receive_data=receive_signal_perallel(27:90,:);
%-----------------------------------fft---------------------------------------
frequency_data_no_cp=fft(receive_data);
frequency_data=[frequency_data_no_cp(35:64,:);frequency_data_no_cp(1:31,:)];
[frequency_data_wide,frequency_data_length]=size(frequency_data);
%---------------------------------信道估计----------------------------------------
channel_condition=zeros(frequency_data_wide,modulate_length);
estimate_data=zeros(frequency_data_wide,modulate_length);
for iii=1:modulate_length
    channel_condition(:,iii)=frequency_data(:,2*iii)./pn_code(:,iii);
    estimate_data(:,iii)=frequency_data(:,(2*iii-1))./channel_condition(:,iii);
end
real_data_temp=[estimate_data(1:30,:);estimate_data(32:61,:)];
%---------------------------------------解调-----------------------------------------
demodulate_data_temp=reshape(real_data_temp,1,[]);
demodulate_data=qamdemod(demodulate_data_temp,2^(modulate_bit));
demodulate_data_bits_temp=reshape(demodulate_data,[],1);
demodulate_data_bits=de2bi(demodulate_data_bits_temp);
real_data_temp1 = reshape(demodulate_data_bits',1,[]);
%--------------------------------------------译码-------------------------------------
[real_data_temp1_wide,real_data_length]=size(real_data_temp1);
yrsgs4=reshape(real_data_temp1 ,4,real_data_temp1_wide*real_data_length/4).';
yrsgs41=bi2de(yrsgs4,'left-msb');
yrsgs41=reshape(yrsgs41,nn,length(yrsgs41)/nn).';
ygsrsdecode=rsdec(gf(yrsgs41,4),nn,kk);
d1=reshape(ygsrsdecode.x',1,[]);
d2=de2bi(d1,'left-msb').';
rx_decode=reshape(d2,1,[]);
%-----------------------------------------------误码率----------------------------------
[error_num,error_ratio]=biterr(inforSource,rx_decode);
error_bit_all(a,1)=error_ratio;
%------------------------------------------------误符号率--------------------------------
error_symbol_data1= reshape(rx_decode,modulate_bit,[])';   %以每组2比特进行分组，输出两列数据
error_symbol_data_receive=bi2de(error_symbol_data1);
error_symbol_data2= reshape(inforSource,modulate_bit,[])';   %以每组2比特进行分组，输出两列数据
error_symbol_data_transmite=bi2de(error_symbol_data2);%输出一列数据
[error_symbol_num,error_symbol_ratio]=symerr(error_symbol_data_receive,error_symbol_data_transmite);
error_symbol_all(a,1)=error_symbol_ratio;
end
end

