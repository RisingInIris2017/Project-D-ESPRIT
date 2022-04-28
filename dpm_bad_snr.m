% 没有模块化拆分的ESPRIT原理算法。
% 这里使用的是中心化ESPRIT，重构集中在信号生成部分。
close all;clear;clc;
%% 参数设置部分
% 真实的波达方向角
theta_real = [20 40];
% 信号参数
independent_signals_number = 2;
% 快拍数
numbers_of_samples = 1000;
t = 0:numbers_of_samples-1;
% 信号功率(dB)
Ps=0;

% 阵列参数
% 估计两个信号，所以只使用两个子阵，这两个子阵具有移不变特性
number_of_subarray = 2;
% 有几对阵元偶，那么每个子阵中的阵元就有几个，
% 所以“阵元偶个数”与“子阵中阵元个数”此处应相等
antenna_in_subarray = 8;
antennas_number = antenna_in_subarray * number_of_subarray;
% 采用线阵阵型布阵
% 阵元间距
distance_between_antenna=0.5;
% 子阵间位移
displacement_between_subarrays=0.25;
% 阵型向量，其第i个元素是第i个阵元对1号参考阵元的位移
d=0:distance_between_antenna:(antenna_in_subarray-1)*distance_between_antenna;
% 阵型确定后，该阵的阵列流形也就此确定
% Ax是子阵X的阵列流形，Ay是子阵Y的阵列流形
Ax=exp(-1i*2*pi*d.'*sin(deg2rad(theta_real)));
Ay=exp(-1i*2*pi*(d+displacement_between_subarrays).'*sin(deg2rad(theta_real)));

%% 信号生成
s = [sin(2*t);cos(2*t)];

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

%% 相同信噪比环境下不同ESPRIT的性能比较
% 三级迭代次数
iter2 = 20;
iter10 = 30;
iter20 = 40;
iter30 = 50;

% 记录均方误差
RMSE_2 = [];
RMSE_10 = [];
RMSE_20 = [];
RMSE_30 = [];
% 完整记录
test_SNR = -15:1:-14;
for SNR = test_SNR % 信噪比的步进
    disp(['当前信噪比：',num2str(SNR),'dB']);
    % 记录重复试验误差
    error_dpm_2 = [];
    error_dpm_10 = [];
    error_dpm_20 = [];
    error_dpm_30 = [];
    
    for repeat = 1:50 % 重复实验repeat次
        disp(['第',num2str(repeat),'次实验开始']);
        % 阵列响应生成
        % 子阵 X 经AWGN信道接收的信号
        X=awgn(Ax*s,SNR,'measured');
        % 子阵 Y 经AWGN信道接收的信号
        Y=awgn(Ay*s,SNR,'measured');

        % 迭代iter2次
        e_1_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_1_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        for iter = 1:iter2
            atn_1 = X'*e_1_x + Y'*e_1_y;
            atn_2 = X'*e_2_x + Y'*e_2_y;
            
            e_1_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_x = e_1_x + X(:,index)*atn_1(index);
            end
            e_1_x = e_1_x/numbers_of_samples;

            e_2_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_x = e_2_x + X(:,index)*atn_2(index);
            end
            e_2_x = e_2_x/numbers_of_samples;

            e_1_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_y = e_1_y + Y(:,index)*atn_1(index);
            end
            e_1_y = e_1_y/numbers_of_samples;

            e_2_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_y = e_2_y + Y(:,index)*atn_2(index);
            end
            e_2_y = e_2_y/numbers_of_samples;
        end
        e_1_x = e_1_x/vecnorm(e_1_x);
        e_1_y = e_1_y/vecnorm(e_1_y);
        e_2_x = e_2_x/vecnorm(e_2_x);
        e_2_y = e_2_y/vecnorm(e_2_y);
        e_1 = [e_1_x;e_1_y];
        e_2 = [e_2_x;e_2_y];
        Es_iter_2 = [e_1, e_2];
        
        % 迭代iter10次
        e_1_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_1_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        for iter = 1:iter10
            atn_1 = X'*e_1_x + Y'*e_1_y;
            atn_2 = X'*e_2_x + Y'*e_2_y;
            
            e_1_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_x = e_1_x + X(:,index)*atn_1(index);
            end
            e_1_x = e_1_x/numbers_of_samples;

            e_2_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_x = e_2_x + X(:,index)*atn_2(index);
            end
            e_2_x = e_2_x/numbers_of_samples;

            e_1_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_y = e_1_y + Y(:,index)*atn_1(index);
            end
            e_1_y = e_1_y/numbers_of_samples;

            e_2_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_y = e_2_y + Y(:,index)*atn_2(index);
            end
            e_2_y = e_2_y/numbers_of_samples;
        end
        e_1_x = e_1_x/vecnorm(e_1_x);
        e_1_y = e_1_y/vecnorm(e_1_y);
        e_2_x = e_2_x/vecnorm(e_2_x);
        e_2_y = e_2_y/vecnorm(e_2_y);
        e_1 = [e_1_x;e_1_y];
        e_2 = [e_2_x;e_2_y];
        Es_iter_10 = [e_1, e_2];
        
        % 迭代iter20次
        e_1_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_1_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        for iter = 1:iter20
            atn_1 = X'*e_1_x + Y'*e_1_y;
            atn_2 = X'*e_2_x + Y'*e_2_y;
            
            e_1_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_x = e_1_x + X(:,index)*atn_1(index);
            end
            e_1_x = e_1_x/numbers_of_samples;

            e_2_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_x = e_2_x + X(:,index)*atn_2(index);
            end
            e_2_x = e_2_x/numbers_of_samples;

            e_1_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_y = e_1_y + Y(:,index)*atn_1(index);
            end
            e_1_y = e_1_y/numbers_of_samples;

            e_2_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_y = e_2_y + Y(:,index)*atn_2(index);
            end
            e_2_y = e_2_y/numbers_of_samples;
        end
        e_1_x = e_1_x/vecnorm(e_1_x);
        e_1_y = e_1_y/vecnorm(e_1_y);
        e_2_x = e_2_x/vecnorm(e_2_x);
        e_2_y = e_2_y/vecnorm(e_2_y);
        e_1 = [e_1_x;e_1_y];
        e_2 = [e_2_x;e_2_y];
        Es_iter_20 = [e_1, e_2];

        % 迭代iter30次
        e_1_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_1_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        e_2_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        for iter = 1:iter30
            atn_1 = X'*e_1_x + Y'*e_1_y;
            atn_2 = X'*e_2_x + Y'*e_2_y;
            
            e_1_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_x = e_1_x + X(:,index)*atn_1(index);
            end
            e_1_x = e_1_x/numbers_of_samples;

            e_2_x = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_x = e_2_x + X(:,index)*atn_2(index);
            end
            e_2_x = e_2_x/numbers_of_samples;

            e_1_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_1_y = e_1_y + Y(:,index)*atn_1(index);
            end
            e_1_y = e_1_y/numbers_of_samples;

            e_2_y = zeros(antenna_in_subarray,1);
            for index = 1:numbers_of_samples
                e_2_y = e_2_y + Y(:,index)*atn_2(index);
            end
            e_2_y = e_2_y/numbers_of_samples;
        end
        e_1_x = e_1_x/vecnorm(e_1_x);
        e_1_y = e_1_y/vecnorm(e_1_y);
        e_2_x = e_2_x/vecnorm(e_2_x);
        e_2_y = e_2_y/vecnorm(e_2_y);
        e_1 = [e_1_x;e_1_y];
        e_2 = [e_2_x;e_2_y];
        Es_iter_30 = [e_1, e_2];
        
        % 求解psi矩阵
        Es_upper_2 = Jupper*Es_iter_2;
        Es_lower_2 = Jlower*Es_iter_2;
        psi_2 = (Es_upper_2'*Es_upper_2)\(Es_upper_2'*Es_lower_2);
        
        Es_upper_10 = Jupper*Es_iter_10;
        Es_lower_10 = Jlower*Es_iter_10;
        psi_10 = (Es_upper_10'*Es_upper_10)\(Es_upper_10'*Es_lower_10);
        
        Es_upper_20 = Jupper*Es_iter_20;
        Es_lower_20 = Jlower*Es_iter_20;
        psi_20 = (Es_upper_20'*Es_upper_20)\(Es_upper_20'*Es_lower_20);

        Es_upper_30 = Jupper*Es_iter_30;
        Es_lower_30 = Jlower*Es_iter_30;
        psi_30 = (Es_upper_30'*Es_upper_30)\(Es_upper_30'*Es_lower_30);
        
        % 求解DOA
        dpm_estimate_2 = sort( ... % 对应
            rad2deg( ... % 化为角度制便于比较
                abs( ... % 取反
                    asin(angle(eig(psi_2).')/pi ... % 求解公式
                  ))));
        error_dpm_2 = [error_dpm_2;abs(dpm_estimate_2 - theta_real)];
        
        dpm_estimate_10 = sort( ... % 对应
            rad2deg( ... % 化为角度制便于比较
                abs( ... % 取反
                    asin(angle(eig(psi_10).')/pi ... % 求解公式
                  ))));
        error_dpm_10 = [error_dpm_10;abs(dpm_estimate_10 - theta_real)];
        
        dpm_estimate_20 = sort( ... % 对应
            rad2deg( ... % 化为角度制便于比较
                abs( ... % 取反
                    asin(angle(eig(psi_20).')/pi ... % 求解公式
                  ))));
        error_dpm_20 = [error_dpm_20;abs(dpm_estimate_20 - theta_real)];
        
        dpm_estimate_30 = sort( ... % 对应
            rad2deg( ... % 化为角度制便于比较
                abs( ... % 取反
                    asin(angle(eig(psi_30).')/pi ... % 求解公式
                  ))));
        error_dpm_30 = [error_dpm_30;abs(dpm_estimate_30 - theta_real)];
        
    end
    RMSE_2 = [RMSE_2; [sqrt(mse(error_dpm_2(:,1))) sqrt(mse(error_dpm_2(:,2)))]];
    RMSE_10 = [RMSE_10; [sqrt(mse(error_dpm_10(:,1))) sqrt(mse(error_dpm_10(:,2)))]];
    RMSE_20 = [RMSE_20; [sqrt(mse(error_dpm_20(:,1))) sqrt(mse(error_dpm_20(:,2)))]];
    RMSE_30 = [RMSE_30; [sqrt(mse(error_dpm_30(:,1))) sqrt(mse(error_dpm_30(:,2)))]];
end
