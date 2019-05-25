function ww = hprod(ftt)

global alpha_h z_h ybar_h 

ww = ybar_h+z_h*(ftt^alpha_h);
end

