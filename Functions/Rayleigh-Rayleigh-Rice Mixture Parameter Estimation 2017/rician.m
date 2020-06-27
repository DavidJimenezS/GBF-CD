function P_ro_wc = rician(ro,sigma_c,nu)


I_0 = besseli(0,(ro.*nu)./(sigma_c^2));

P_ro_wc = (ro./(sigma_c^2)).*exp(-((ro.^2) + (nu^2))./(2*(sigma_c^2))).*I_0;

end