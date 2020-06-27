function P_ro_wn = rayleigh(ro,b_n)

P_ro_wn = (ro./(b_n^2)).*exp(-(ro.^2)./(2*(b_n^2)));