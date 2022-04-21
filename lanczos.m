% 利用Lanczos算法，重新估计Us矩阵。
close all;clear;clc;
%% 参数设置部分
% 在这个脚本变为函数以后，这一部分应删去
% 设置独立信源的个数，在文献[2]中称为q
L = 2;
% 设置子阵的个数，这个参数在文献[1-6]中称为p，在文献[2]中称为K
number_of_subarray = 3;
% 设置子阵中阵元的个数，这个参数在文献[1]中就是Mk，在文献[1-6]中称为m
% 我们假设所有的子阵阵元个数都相同
antenna_in_subarray = 4;
% 总阵元个数
M = number_of_subarray * antenna_in_subarray;
% 设置真实 DOA
doa_vector = [10 30];
% 设置信号功率和信噪比，均为分贝，并求得噪声功率
Ps = 0;
SNR = 10;
Pn = Ps - SNR;
% 设置信号时间采样的快拍数
N = 500;

%% 构造J算子
% 这里由于各子阵都相同，因此Jk矩阵（不论上一横还是下一横）随着k取遍1:L，都相同
% Jk上一横
Jupperk = [eye(antenna_in_subarray-1,antenna_in_subarray-1),zeros(antenna_in_subarray-1,1)];
% Jk下一横
Jlowerk = [zeros(antenna_in_subarray-1,1),eye(antenna_in_subarray-1,antenna_in_subarray-1)];

% J上一横
Jupper = kron(eye(number_of_subarray),Jupperk);
% J下一横
Jlower = kron(eye(number_of_subarray),Jlowerk);

%% xk生成部分
xkg_command = "[";
for k = 1:number_of_subarray
    xkg_command = xkg_command + "x_"+num2str(k)+",";
end
xkg_command = extractBefore(xkg_command,strlength(xkg_command))+"] = xk_generator(Ps, L, N, SNR, number_of_subarray, doa_vector, antenna_in_subarray);";
eval(xkg_command);
clear xkg_command;

%% AC迭代的参数
% 在这个脚本变为函数后，这一段应删去
% 规定子阵拓扑图，生成文献[2-XB04]中的拉普拉斯矩阵
% 简单起见，我们设置三角形拓扑，1与2相连，2与3相连，1与3不相连
% 注意，这里不允许使用全连通图，否则W矩阵没有我们需要的特性。
laplacianMatrix = [1,-1,0;0,2,-1;0,0,1];
% 利用文献[2-XB04]中的公式W=I-aL生成W矩阵
% 拉普拉斯矩阵前面的标量系数alpha必须大于0
% 最佳的alpha取值为最大的特征值与第二小的特征值算术平均取倒数
eigLapMat = sort(eig(laplacianMatrix),'descend');
alpha = 2/(eigLapMat(1) + eigLapMat(number_of_subarray - 1));
W = eye(number_of_subarray) - alpha * laplacianMatrix;

%% Lanczos算法
% 暂时需要使用的阵列协方差矩阵，之后需用AC算法重写
x = [];
for k = 1:number_of_subarray
    x = [x;eval("x_"+num2str(k))];
end
Rxx = x*x';
% 迭代次数，按照文献[1]的说法，略大于L即可
num_of_lanczos = 2*L;
% 迭代初始值
beta = 1;
v = [zeros(M,1),rand(M,1)];
a_element = zeros(1,num_of_lanczos);
b_element = [1,zeros(1,num_of_lanczos - 1)];
for iter = 1:num_of_lanczos
    % 关于迭代下标的说明
    % 在这个定义下，v向量的下标从0开始取，alpha, beta, w的下标都是从1开始取。
    % 但是MATLAB中，v的下标也只能从1开始取，所以凡是涉及v的迭代，
    % 下标都要在文献中的下标基础上加一。
    w = Rxx*v(:,iter + 1);
    a_element(iter) = v(:,iter + 1)'*w;
    if(iter < num_of_lanczos)
        w = w - a_element(iter)*v(:,iter + 1)-b_element(iter)*v(:,iter);
        b_element(iter + 1) = vecnorm(w);
        v = [v,w/b_element(iter + 1)];
    end
end
a_part = diag(a_element);
b_part_upper = diag(b_element(2:length(b_element)));
size_of_b_part = length(b_part_upper);
b_part_upper = [[zeros(size_of_b_part,1),b_part_upper];zeros(1,size_of_b_part + 1)];
Tmat = a_part + b_part_upper + b_part_upper.';
% 尝试：经观察，Tmat的非对角元虚部都接近0，假如我们认为Tmat是一个Hermitian矩阵，我们对它作svd。
% 在 5 位有效数字内，它们的结果的确是一样的。
[~,Sigma,Vvec] = svd(Tmat);
Es = [];
for iter = 1:L
    Es = [Es, v(:,2:size(v,2))*Vvec(:,iter)];
end
%% ESRPIT
% 暂时使用中心化ESPRIT
Es_upper = Jupper*Es;
Es_lower = Jlower*Es;
psi = (Es_upper'*Es_upper)\(Es_upper'*Es_lower);
eig_psi = eig(psi).';
phase_psi = atan2(imag(eig_psi),real(eig_psi));
doa_estimate = sort(abs(asin(phase_psi/pi)*180/pi));
error = abs(doa_estimate - doa_vector)