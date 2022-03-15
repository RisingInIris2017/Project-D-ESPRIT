% 设置独立信源的个数
L = 2;
% 设置子阵的个数
number_of_subarray = 3;
% 设置子阵中阵元的个数
antenna_in_subarray = 5;
% 总阵元个数
M = number_of_subarray * antenna_in_subarray;
% 设置真实 DOA
doa_vector = 10*(1:L);
% 设置信号功率和信噪比，均为分贝，并求得噪声功率
Ps = 0;
SNR = 20;
Pn = Ps - SNR;
% 生成各独立信号的时域采样构成的向量
signal_vector = signal_generator(L);