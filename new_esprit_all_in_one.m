% 没有模块化拆分的ESPRIT原理算法。
% 这里使用的是中心化ESPRIT，重构集中在信号生成部分。
close all;clear;clc;
%% 参数设置部分
% 真实的波达方向角
theta_real = [0 40];
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

%% 相同信噪比环境下不同ESPRIT的性能比较
% 记录均方误差
RMSE_classical = [];
RMSE_dpm = [];
% 完整记录
totalerror_classical = [];
totalerror_dpm = [];
test_SNR = -15:5:20;
for SNR = test_SNR % 信噪比的步进
    disp(['当前信噪比：',num2str(SNR),'dB']);
    % 记录重复试验误差
    error_classical = [];
    error_dpm = [];
    
    for repeat = 1:50 % 重复实验repeat次
        disp(['第',num2str(repeat),'次实验开始']);
        % 阵列响应生成
        % 子阵 X 经AWGN信道接收的信号
        X=awgn(Ax*s,SNR,'measured');
        % 子阵 Y 经AWGN信道接收的信号
        Y=awgn(Ay*s,SNR,'measured');

        % 传统中心化ESPRIT
        % 合并两子阵的数据，得数据矩阵Z
        Z=[X;Y];
        % 求协方差矩阵R
        R=Z*Z'/numbers_of_samples;%协方差矩阵
        %特征值分解
        [EV1,D1]=eig(R);
        %取independent_signals_number个大特征值构成信号子空间
        Es=EV1(:,antennas_number-independent_signals_number+1:antennas_number);
        % 利用移不变特性，分解信号子空间
        Ex=Es(1:antenna_in_subarray,:);
        Ey=Es(antenna_in_subarray+1:antennas_number,:);
        %计算F矩阵，从F的特征值求解DOA
        F=pinv(Ex)*Ey;
        classical_estimate = sort(rad2deg(-asin(angle(eig(F).')/2/pi/displacement_between_subarrays)));%DOAestmate
        error_classical = [error_classical;abs(classical_estimate - theta_real)];
        
        % DPM-ESPRIT
        % 逐行更改中心化算法中“不够分布式”的地方
        % 首先，不允许直接获取Z矩阵，必须从X和Y（两个子阵的响应）入手

        % DPM求ei,k向量
        % 随机四个初始向量
        % 子阵X对1号信号的估计
        e_1_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        % 子阵X对2号信号的估计
        e_2_x = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        % 子阵Y对1号信号的估计
        e_1_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);
        % 子阵Y对2号信号的估计
        e_2_y = rand(antenna_in_subarray,1)+1i*rand(antenna_in_subarray,1);

        % 迭代
        for iter = 1:10
            % atn_1和atn_2存储着快拍数个标量，分别包含有整个阵列对1号信号和2号信号的估计的信息
            atn_1 = X'*e_1_x + Y'*e_1_y;
            atn_2 = X'*e_2_x + Y'*e_2_y;

            % 迭代子阵X对1号信号的估计
            e_1_x = zeros(antenna_in_subarray,1);
            % 累加求和取平均
            for index = 1:numbers_of_samples
                e_1_x = e_1_x + X(:,index)*atn_1(index);
            end
            e_1_x = e_1_x/numbers_of_samples;
            % 迭代子阵X对2号信号的估计
            e_2_x = zeros(antenna_in_subarray,1);
            % 累加求和取平均
            for index = 1:numbers_of_samples
                e_2_x = e_2_x + X(:,index)*atn_2(index);
            end
            e_2_x = e_2_x/numbers_of_samples;
            % 迭代子阵Y对1号信号的估计
            e_1_y = zeros(antenna_in_subarray,1);
            % 累加求和取平均
            for index = 1:numbers_of_samples
                e_1_y = e_1_y + Y(:,index)*atn_1(index);
            end
            e_1_y = e_1_y/numbers_of_samples;
            % 迭代子阵Y对2号信号的估计
            e_2_y = zeros(antenna_in_subarray,1);
            % 累加求和取平均
            for index = 1:numbers_of_samples
                e_2_y = e_2_y + Y(:,index)*atn_2(index);
            end
            e_2_y = e_2_y/numbers_of_samples;
        end
        % 归一化
        e_1_x = e_1_x/vecnorm(e_1_x);
        e_1_y = e_1_y/vecnorm(e_1_y);
        e_2_x = e_2_x/vecnorm(e_2_x);
        e_2_y = e_2_y/vecnorm(e_2_y);
        % 拼接
        e_1 = [e_1_x;e_1_y];
        e_2 = [e_2_x;e_2_y];
        Es = [e_1, e_2];

        % 构造J算子
        % 这里由于各子阵都相同，因此Jk矩阵（不论上一横还是下一横）随着k取遍1:L，都相同
        % Jk上一横
        Jupperk = [eye(antenna_in_subarray-1,antenna_in_subarray-1),zeros(antenna_in_subarray-1,1)];
        % Jk下一横
        Jlowerk = [zeros(antenna_in_subarray-1,1),eye(antenna_in_subarray-1,antenna_in_subarray-1)];

        % J上一横
        Jupper = kron(eye(number_of_subarray),Jupperk);
        % J下一横
        Jlower = kron(eye(number_of_subarray),Jlowerk);

        % 求解psi矩阵
        Es_upper = Jupper*Es;
        Es_lower = Jlower*Es;
        psi = (Es_upper'*Es_upper)\(Es_upper'*Es_lower);
        % 求解DOA
        dpm_estimate = sort( ... % 对应
            rad2deg( ... % 化为角度制便于比较
                abs( ... % 取反
                    asin(angle(eig(psi).')/pi ... % 求解公式
                  ))));
        error_dpm = [error_dpm;abs(dpm_estimate - theta_real)];
    end
    totalerror_classical = [totalerror_classical, error_classical];
    totalerror_dpm = [totalerror_dpm, error_dpm];
    RMSE_classical = [RMSE_classical; [sqrt(mse(error_classical(:,1))) sqrt(mse(error_classical(:,2)))]];
    RMSE_dpm = [RMSE_dpm; [sqrt(mse(error_dpm(:,1))) sqrt(mse(error_dpm(:,2)))]];
end
%% 画图
figure;
% 画两种ESPRIT对1号信号的估计均方误差随信噪比变化图
figure(1);
classical_1 = plot(test_SNR, RMSE_classical(:,1), 'bo-', 'LineWidth', 4);
hold on
dpm_1 = plot(test_SNR, RMSE_dpm(:,1), 'ro-', 'LineWidth', 4);
legend('经典ESPRIT','分布式幂方法ESPRIT','FontSize',25);
axis tight
set(gca,'linewidth',2,'fontsize',25);
xlabel('信噪比/dB','FontWeight','bold','FontSize',25);
ylabel('均方误差/角度','FontWeight','bold','FontSize',25);
title('两种ESPRIT对1号信号的估计均方误差随信噪比变化图','FontWeight','bold','FontSize',25);
% 画两种ESPRIT对2号信号的估计均方误差随信噪比变化图
figure(2);
classical_2 = plot(test_SNR, RMSE_classical(:,2), 'bo-', 'LineWidth', 4);
hold on
dpm_2 = plot(test_SNR, RMSE_dpm(:,2), 'ro-', 'LineWidth', 4);
legend('经典ESPRIT','分布式幂方法ESPRIT','FontSize',25);
axis tight
set(gca,'linewidth',2,'fontsize',25);
xlabel('信噪比/dB','FontWeight','bold','FontSize',25);
ylabel('均方误差/角度','FontWeight','bold','FontSize',25);
title('两种ESPRIT对2号信号的估计均方误差随信噪比变化图','FontWeight','bold','FontSize',25);
datatip(classical_1,0,RMSE_classical(4,1),'SnapToDataVertex','off');
datatip(classical_2,0,RMSE_classical(4,2),'SnapToDataVertex','off');
datatip(dpm_1,0,RMSE_dpm(4,1),'SnapToDataVertex','off');
datatip(dpm_2,0,RMSE_dpm(4,2),'SnapToDataVertex','off');
