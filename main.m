%% 信号的参数
% 设置独立信源的个数
L = 2;
% 设置子阵的个数，这个参数在文献[2]中称为K
number_of_subarray = 3;
% 设置子阵中阵元的个数，这个参数在文献[1]中就是Mk
% 我们假设所有的子阵阵元个数都相同
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

%% AC迭代的参数
% 规定子阵拓扑图，生成文献[2-XB04]中的拉普拉斯矩阵
% 简单起见，我们设置三角形拓扑，1与2相连，1与3相连，2与3不相连
% 注意，这里不允许使用全连通图，否则W矩阵没有我们需要的特性。
laplacianMatrix = [2,-1,-1;-1,1,0;-1,0,1];
% 利用文献[2-XB04]中的公式W=I-aL生成W矩阵
% 拉普拉斯矩阵前面的标量系数alpha必须大于0
% 最佳的alpha取值为最大的特征值与第二小的特征值算术平均取倒数
% 在上述拓扑图规定中，拉普拉斯矩阵必然是实对称的，满足svd与eig等效定理条件。
% 下面用svd()取代了eig()
eigLapMat = svd(laplacianMatrix);
alpha = 2/(eigLapMat(1) + eigLapMat(number_of_subarray - 1));
W = eye(number_of_subarray) - alpha * laplacianMatrix;

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
