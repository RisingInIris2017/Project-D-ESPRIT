% 利用分布式幂方法估计Us矩阵。

%% 参数设置部分
% 在这个脚本变为函数以后，这一部分应删去
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

% %% 构造J算子
% % 这里由于各子阵都相同，因此Jk矩阵（不论上一横还是下一横）随着k取遍1:L，都相同
% % Jk上一横
% Jupperk = [eye(antenna_in_subarray-1,antenna_in_subarray-1),zeros(antenna_in_subarray-1,1)];
% % Jk下一横
% Jlowerk = [zeros(antenna_in_subarray-1,1),eye(antenna_in_subarray-1,antenna_in_subarray-1)];
% % 借助字符串拼接和eval来构造J算子
% command_upper = "blkdiag(";
% command_lower = "blkdiag(";
% for k=1:L
%     command_upper = command_upper + "Jupperk,";
%     command_lower = command_lower + "Jlowerk,";
% end
% command_upper = extractBefore(command_upper, strlength(command_upper)) + ");";
% command_lower = extractBefore(command_lower, strlength(command_lower)) + ");";
% % J上一横
% Jupper = eval(command_upper);
% % J下一横
% Jlower = eval(command_lower);

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

%% DPM算法
% 设置一个足够大的AC迭代次数Q
Q = 10;
% 设置一个足够大的AC迭代次数P
P = 10;
% 计算出迭代矩阵W^P
WP = W^P;
% 初始化U，这个矩阵在2.61式下方定义
U = [];
% 初始化ulk(q=0)
for l = 1:M
    for k = 1:number_of_subarray
        eval("ul_"+num2str(l)+"k_"+num2str(k)+"=rand(antenna_in_subarray,1);")
    end
end
% for q = 1:Q % 迭代
    % STEP 1：迭代求ul'
    for l = 1:M    % 这个循环是2.64式
        for k0 = 1:number_of_subarray
            eval("xtl_"+num2str(l)+"k_"+num2str(k0)+"=0;");
            for k = 1:number_of_subarray
                eval("xtl_"+num2str(l)+"k_"+num2str(k0)+"="...
                    +"xtl_"+num2str(l)+"k_"+num2str(k0)+"+number_of_subarray*WP(k0,k)*x_"+num2str(k)+"'*ul_"+num2str(l)+"k_"+num2str(k)+";");
            end
        end
    end
    for l = 1:M    % 这个循环是2.65式
        for k = 1:number_of_subarray
            eval("x_"+num2str(k)+"xtl_"+num2str(l)+"k_"+num2str(k)+"=zeros(antenna_in_subarray,N);");
            for t = 1:N
                eval("x_"+num2str(k)+"xtl_"+num2str(l)+"k_"+num2str(k)+"(:,t)=x_"+num2str(k)+"(:,t)*xtl_"+num2str(l)+"k_"+num2str(k)+"(t);");
            end
            eval("ul_"+num2str(l)+"k_"+num2str(k)+"="+"sum(x_"+num2str(k)+"xtl_"+num2str(l)+"k_"+num2str(k)+",2)/N;");
        end
    end
    for l = 1:M    % 2.59式，用ulk拼接得到ul矩阵
        eval("ul_"+num2str(l)+"=[];");
        for k = 1:number_of_subarray
            eval("ul_"+num2str(l)+"k_"+num2str(k)+"_tp="+"ul_"+num2str(l)+"k_"+num2str(k)+"';");
            eval("ul_"+num2str(l)+"=[ul_"+num2str(l)+",ul_"+num2str(l)+"k_"+num2str(k)+"_tp"+"];");
        end
    end
    % STEP 2: 迭代求ul
    % 以下将ul向量称为u_real_l，表示其不是带撇号的中间变量
    u_real_l_1=ul_1;
    for l = 2:M
        eval("U = [U u_real_l_"+num2str(l-1)+"]");
        eval("u_real_l_"+num2str(l)+"=ul_"+num2str(l)+"-U*U'*ul_"+num2str(l)+";");
    end
% end