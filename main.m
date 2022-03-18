%% 信号的参数
% 设置独立信源的个数
L = 2;
% 设置子阵的个数，这个参数在文献[6]中称为K
number_of_subarray = 3;
% 设置子阵中阵元的个数，这个参数在文献[1]中就是Mk
% 我们架设所有的子阵阵元个数都相同
antenna_in_subarray = 5;
% 总阵元个数
M = number_of_subarray * antenna_in_subarray;
% 设置真实 DOA
doa_vector = 10*(1:L);
% 设置信号功率和信噪比，均为分贝，并求得噪声功率
Ps = 0;
SNR = 20;
Pn = Ps - SNR;
% 设置信号时间采样的快拍数
N = 1000;

%% 幂方法的参数
% 设置幂方法的迭代次数Q
Q = 5;

%% 生成物理量
% 生成K个阵列输出xk(t)
xkg_command = "[";
for k = 1:number_of_subarray
    xkg_command = xkg_command + "x_"+num2str(k)+",";
end
xkg_command = extractBefore(xkg_command,strlength(xkg_command))+"] = xk_generator(Ps, L, N, SNR, number_of_subarray, doa_vector, antenna_in_subarray);";
eval(xkg_command);
