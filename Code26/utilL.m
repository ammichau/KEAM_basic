function [util] = utilL(LL)

global  eta mu 

util = -mu*(LL^(1+eta))/(1+eta);
end