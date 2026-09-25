function d_eprime = de2_dh(ee,hh)

global  psi alpha_e

d_eprime = alpha_e*psi*ee*(hh^(psi-1));
end