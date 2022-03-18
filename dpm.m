% 利用分布式幂方法估计Us矩阵。
% 初始化Us
Us = [];
% 设置独立信源的个数，文献1-6称为K
L = 2;
% 初始化各子阵的阵元个数，文献1-6称为Mk
antenna_in_subarray = 3;
% 设置子阵的个数
number_of_subarray = 3;
% 计算求得总阵元个数M，用于l的迭代
M = number_of_subarray * antenna_in_subarray;
% 设置一个足够大的AC迭代次数Q
Q = 10;

%% 构造J算子
% 这里由于各子阵都相同，因此Jk矩阵（不论上一横还是下一横）随着k取遍1:L，都相同
% Jk上一横
Jupperk = [eye(antenna_in_subarray-1,antenna_in_subarray-1),zeros(antenna_in_subarray-1,1)];
% Jk下一横
Jlowerk = [zeros(antenna_in_subarray-1,1),eye(antenna_in_subarray-1,antenna_in_subarray-1)];
% 借助字符串拼接和eval来构造J算子
command_upper = "blkdiag(";
command_lower = "blkdiag(";
for k=1:L
    command_upper = command_upper + "Jupperk,";
    command_lower = command_lower + "Jlowerk,";
end
command_upper = extractBefore(command_upper, strlength(command_upper)) + ");";
command_lower = extractBefore(command_lower, strlength(command_lower)) + ");";
% J上一横
Jupper = eval(command_upper);
% J下一横
Jlower = eval(command_lower);

%% DPM算法
for l = 1:M
    for k = 1:L
        eval("ul_"+num2str(l)+"k_"+num2str(k)+"=rand(antenna_in_subarray,1);")
    end
end
