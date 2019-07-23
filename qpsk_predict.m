function [error_ratio_all]=qpsk_predict()
for a=1:21
 error_ratio_all(a,1)=qfunc(sqrt(10^((a-1)/10)));
end
end

