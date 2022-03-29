function signal_vector = signal_generator(Ps, L, N)
    % 构造单音信号的时间序列
    % 为仿真简单起见，不再调制，直接检测单音正弦信号
    % 依照《阵列信号处理的理论与应用》p102，这些信号的频率应都相同，依据初相的不同彼此区分
    init_phase_vector = pi/2*(0:L-1);
        % 设置时间序列的采样点数
        t = 0:N-1;
        signal_vector = ones(L,N);
    figure
    for index = 1:L
        signal_vector(index,:)=10^(Ps/10)*sin(2e6*pi/(N-1)*t + init_phase_vector(index));
        plot(signal_vector(index,:));
        hold on
    end
end