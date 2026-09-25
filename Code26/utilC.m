function [util] = utilC(cc)

global  crra phi_c

util = phi_c*(cc^(1-crra))/(1-crra);
end

