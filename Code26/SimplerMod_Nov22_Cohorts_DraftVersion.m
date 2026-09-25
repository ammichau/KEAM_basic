%SimplerMod_May17_shell.m
clear all

%---------------------------------
%Uses solved Baseline as 1940s cohort calibration
%---------------------------------

%Set Main Path
    MAINdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');

%{
%---------------------------------
%Experiment 1: 1950's Cohort
    %---------------------------------
    clear all
    %Set Main Path
        MAINdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab';
        solver = fullfile(MAINdir,'SimplerMod_May17_splines.m');
        simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
    %Set Output Paths
        SOLUTIONdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Solution\Cohort1950';   
        OUTPUTdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Output\Cohort1950';
        file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
        
     %Choose Experiment    
        kapscale=  1.01;
        wagegapscale= 1.07; %1.07 gives decline as is
        rtoexpscale=  1.095;   %1.095;  
        
    %OPTIONS:
    %   %EXPORT DATA?
            sumstatout = 1;
            panelout = 1;    
            
    %Run
        cd(MAINdir);
        run(solver)
        run(simulator)

%}
%---------------------------------
%Experiment 1: 1960's Cohort
%---------------------------------
clear all
%Set Main Path
    MAINdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Solution\Cohort1960';   
    OUTPUTdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Output\Cohort1960';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale= 0.84;
    wagegapscale=  1.19;  %0.267/0.35;
    rtoexpscale=1.15;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    %run(solver)
    run(simulator)


%---------------------------------
%Experiment 1: 1970's Cohort
%---------------------------------
clear all
%Set Main Path
    MAINdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Solution\Cohort1970';   
    OUTPUTdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Output\Cohort1970';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale=1.0; %1.0
    wagegapscale= 1.4; %1.4
    rtoexpscale=1.35; %1.35
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)


%---------------------------------
%Experiment 1: 1980's Cohort
%---------------------------------
clear all
%Set Main Path
    MAINdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab';
    solver = fullfile(MAINdir,'SimplerMod_May17_splines.m');
    simulator = fullfile(MAINdir,'SimplerMod_May17_sim.m');
%Set Output Paths
    SOLUTIONdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Solution\Cohort1980';   
    OUTPUTdir = 'C:\Users\IRAMM03\Desktop\KEAM\Matlab\Output\Cohort1980';
    file = fullfile(OUTPUTdir,'SimulStats2.xls')  ;      
    
 %Choose Experiment    
    kapscale=1.037;
    wagegapscale=1.41;
    rtoexpscale=1.37;
    
%OPTIONS:
%   %EXPORT DATA?
        sumstatout = 1;
        panelout = 1;    
        
%Run
    cd(MAINdir);
    run(solver)
    run(simulator)    
