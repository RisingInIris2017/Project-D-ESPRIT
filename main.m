%% 信号的参数
% 设置独立信源的个数
L = 2;
% 设置子阵的个数
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
% 生成各独立信号的时域采样构成的向量
signal_vector = signal_generator(Ps, L, N);