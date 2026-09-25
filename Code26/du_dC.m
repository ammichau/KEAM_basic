function [d_util] = du_dC(cc)

global  crra phi_c

d_util = phi_c*(cc^(-crra));
end
