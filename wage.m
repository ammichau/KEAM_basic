function ww = wage(ftt,ee)

global xi tau_wf gam_e 

ww = tau_wf*(ftt+gam_e*(ee^xi));
end

