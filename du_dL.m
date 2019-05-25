function [d_util] = du_dL(LL)

global  eta mu 

d_util = -mu*(LL^(eta));
end