# KEAM_basic
Quick and dirty Matlab code to check KEAM mechanism

Stripped down model
%--------------------------------------------------------------------------------    
%HH problem with 
%   1) Nonemployed search intensity 
%   2) Employed intensive margin + quit choice
%   3) Three ages: young+ middle +old
%   4) No savings
%State (i,t,e,y,z)
%   -i= woman's type: 2x2x4
%       -omi = fixed wage type
%       -kapi =  fixed cost of work type
%       -kapy = extra fixed cost of work while young
%   -t = age
%   -e = experience
%   -y = husband's state (1=e, 2=recent u, 3=u)
%   -z = Macro state 1= high, 2=low

%--------------------------------------------------------------------------------

Goal: compstats on
  1) Aggregate Hours + Employment
  2) "Career" composition
  3) Business Cycle volatilities
  
Over potential mechanisms
  a) Returns to tenure
  b) Gender wage gap
  c) Fixed cost of participation
  
Includes:
    -SimplerMod_May17_splines.m  ---- Computes decision rules
    -SimplerMod_May17_sim.m   ----- Simulates model generated data
    -Bunch of functions.
