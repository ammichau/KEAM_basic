%SimplerMod_May1 7_shell.m
clear all

%---------------------------------
%Solve the Baseline
%---------------------------------

%Set Main Path
    MAINdir = 'C:\Users\amichau9\Desktop\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines_betteronBC.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\amichau9\Desktop\Matlab\Solution\Baseline';   
    OUTPUTdir = 'C:\Users\amichau9\Desktop\Matlab\Output\Baseline';
    file = fullfile(OUTPUTdir,'SimulStats.xls')  ;      
 
    %Choose Experiment    
    kapscale=1;
    wagegapscale=1;
    rtoexpscale=1;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 0;
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)



%---------------------------------
%Experiment 1: Compstat on RoE
%---------------------------------
clear all
%Set Main Path
    MAINdir = 'C:\Users\ammic\Desktop\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines_betteronBC.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\ammic\Desktop\Matlab\Solution\RoE_incr';   
    OUTPUTdir = 'C:\Users\ammic\Desktop\Matlab\Output\RoE_incr';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale=1;
    wagegapscale=1;
    rtoexpscale=1.2;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)
    
%---------------------------------
%Experiment 2: Compstat on Wage Gap
%---------------------------------
clear all

%Set Main Path
    MAINdir = 'C:\Users\ammic\Desktop\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines_betteronBC.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\ammic\Desktop\Matlab\Solution\Wgap_dcr';   
    OUTPUTdir = 'C:\Users\ammic\Desktop\Matlab\Output\Wgap_dcr';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale=1;
    wagegapscale=0.8;
    rtoexpscale=1;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)    
    
%-----------------------------------------------
%Experiment 3: Compstat on Fixed Cost of Work
%-----------------------------------------------
clear all

%Set Main Path
    MAINdir = 'C:\Users\ammic\Desktop\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines_betteronBC.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\ammic\Desktop\Matlab\Solution\Kap_dcr';   
    OUTPUTdir = 'C:\Users\ammic\Desktop\Matlab\Output\Kap_dcr';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale=0.9;
    wagegapscale=1;
    rtoexpscale=1;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)     
    