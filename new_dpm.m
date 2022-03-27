% 利用文献2中的分布式幂方法，重新估计Us矩阵。

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
SNR = 20;
Pn = Ps - SNR;
% 设置信号时间采样的快拍数
N = 1000;

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
% 在这个脚本变为函数后，这一段应删去
xkg_command = "[";
for k = 1:number_of_subarray
    xkg_command = xkg_command + "x_"+num2str(k)+",";
end
xkg_command = extractBefore(xkg_command,strlength(xkg_command))+"] = xk_generator(Ps, L, N, SNR, number_of_subarray, doa_vector, antenna_in_subarray);";
eval(xkg_command);

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
% alpha = 2/(eigLapMat(1) + eigLapMat(number_of_subarray - 1));
alpha = 0.5;
W = eye(number_of_subarray) - alpha * laplacianMatrix;


%% 文献2中的DPM算法
% 迭代求Es的次数
PM = 10;
% 生成迭代初始值
% 迭代初始的Es（Us，信号子空间）矩阵，每antenna_in_subarray行归属一个子阵，每列归属一个信号
% Es = rand(antenna_in_subarray * number_of_subarray, L);
Es = 0.8*ones(antenna_in_subarray * number_of_subarray, L);
% 迭代初始的a矩阵，每列归属一个t时刻，每行均分给各个子阵，每层归属一个信号
a = zeros(N, number_of_subarray, L);
for pm = 1:PM
    for q = 1:L
        % 将属于第q个信号的向量提取出来
        eq = Es(:,q);
        % 重排成矩阵，每列是一个m维复向量eq,p
        eqp = reshape(eq,antenna_in_subarray,number_of_subarray);
        for p = 1:number_of_subarray
            eval("a(:,p,q)=x_"+num2str(p)+"'*eqp(:,p);");
        end
    end
    % 暂时用求算术平均代替AC
    at = number_of_subarray * mean(a,2);
    for q = 1:L
        for p = 1:number_of_subarray
            for t = 1:N
                Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L) = Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L) + eval("x_"+num2str(p)+"(:,t)*at(t,1,L);");
            end
            Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L) = Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L)/N;
        end  
    end
end
% 暂时用简单的方法归一化，之后需用AC替代
for q = 1:L
    for p = 1:number_of_subarray
        Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L) = Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L)/vecnorm(Es((p-1)*antenna_in_subarray+1:p*antenna_in_subarray,L));
    end  
end
% 暂时用J算子直接处理Es，检验Es的生成算法是否有效
% 之后需另外重写算法替代
Es_upper = Jupper*Es;
Es_lower = Jlower*Es;
eig_psi = eig((Es_upper'*Es_upper)\Es_upper'*Es_lower);
phase_psi=atan2(imag(eig_psi),real(eig_psi));
doa_estimate = sort(abs(asin(phase_psi/pi)'*180/pi));
error = abs(doa_estimate - doa_vector)