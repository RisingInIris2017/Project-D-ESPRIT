function [varargout] = xk_generator(Ps, L, N, SNR, number_of_subarray, doa_vector, antenna_in_subarray)
    % 生成各独立信号的时域采样构成的向量
    signal_vector = signal_generator(Ps, L, N);

    % 生成K个阵列流形
    amg_command = "[";
    for k = 1:number_of_subarray
        amg_command = amg_command + "A_"+num2str(k)+",";
    end
    amg_command = extractBefore(amg_command,strlength(amg_command))+"] = array_manifold_generator(number_of_subarray, antenna_in_subarray, doa_vector);";
    eval(amg_command);

    % 计算无噪情况下的输出xk(t)
    for k = 1:number_of_subarray
        eval("x_"+num2str(k)+"=A_"+num2str(k)+"*signal_vector;");
        % 简化问题，输出无噪的xk
        eval("varargout{k}=x_"+num2str(k)+";");
    end

    % 加噪
%     for k = 1:number_of_subarray
%         n = wgn(antenna_in_subarray,N,Ps-SNR);
%         eval("varargout{k}=x_"+num2str(k)+"+n;");
%     end
end