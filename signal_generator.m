function signal_vector = signal_generator(L)
% 构造单音信号的时间序列
% 设置单频正弦信号的频率
omega_vector = 1:L;
    % 设置时间序列的采样点数
    N = 1000;
    t = 0:N-1;
    signal_vector = ones(L,N);
    % 设信号功率为0dB，用于规定信号幅值
    Ps = 0;
for index = 1:L
    signal_vector(index,:)=10^(Ps/10)*sin(2*pi*omega_vector(index)/(N-1)*t);
end
end