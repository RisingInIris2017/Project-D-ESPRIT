function phi = phi_calc(theta,zeta_k)
 % ÓÃÓÚ¼ÆËãphi(theta, zeta_k)
    phi = exp(1j*pi*zeta_k.'*[sin(theta);cos(theta)]);
end