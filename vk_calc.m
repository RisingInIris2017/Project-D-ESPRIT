function vk = vk_calc(theta,zeta_k_i)
% ���ڼ���vk����
vk = ones(size(zeta_k_i,2),1);
    for index = 2:size(zeta_k_i,2)
        vk(index) = exp(1j*pi*zeta_k_i(:,index).'*[sin(theta);cos(theta)]);
    end
end

