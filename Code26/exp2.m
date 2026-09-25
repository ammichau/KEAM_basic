function eprime = exp2(ee,hh)

global delta_e psi alpha_e

eprime = (1-delta_e)*ee+alpha_e*ee*(hh^(psi));
end
