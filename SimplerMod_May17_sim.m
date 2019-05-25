% SimplerMod_vMay17.m
%clear all

%Set Paths
    MAINdir = 'C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab';
%Declare if this is the calibration script or the counterfactual script.
    calib= 1;
%--------------------------------------------------------------------------------    
%HH problem with 
%   1) Nonemployed search intensity 
%   2) Employed intensive margin
%   3) Two ages: young +old
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

%Simulates deterministic age life-cycle. 
%   Ex-ante types:
%       -Experience draw on first half of grid
%       -Log-normal ish interpolation on fixed type
%       -uniform on kappa
%--------------------------------------------------------------------------------

%load('paras')
%load('gH.mat', 'gH')
%load('gQ.mat', 'gQ')
%load('gS.mat', 'gS')

if calib==1
    
%Set panel size
Nsim=100; %Number of individuals PER SINGLE YEAR BIRTH COHORT to simulate
Ngen = 2019-1955; %Number of single year birth cohorts.
Tsim=Ngen*tsize;
Tsim_i= 39*tsize; %For the simple model: 25-39, 41-55, 56-64
Nsim_i=100*39;
ninter = 10; %how fine to make individ. type grids.

%Set initial distribution of states.
lamHbar=MarkovStationary(squeeze(lamH(1,:,:)));
Nw_init = 0.08; %Wife unemployment (Choosing high because some will probably take non-employment


%Set age paths for everybody
for t=1:14*tsize
    ageT(t)= 1;
end
for t=14*tsize+1:14*tsize+14*tsize
   ageT(t)=2; 
end
for t= 14*tsize+14*tsize+1:Tsim_i
    ageT(t)=3;
end

%Set individual types and interpolated policies.
    %First types across wages
    Wsimgrid=linspace(ftW(1,1),ftW(1,nI),ninter);
    muW = 0.7*ftW(1,1)+0.3*ftW(1,nI);
    sigW = (ftW(1,nI)-ftW(1,1))/2;
    pd0W = makedist('lognormal','mu',muW,'sigma',sigW);
    pdW = truncate(pd0W,ftW(1,1),ftW(1,nI)) ;
    cdfW(1) = cdf(pdW,Wsimgrid(1));
    ftW_i(1:floor(cdfW(1)*Nsim)) = Wsimgrid(1);
    for i=2:ninter-1
        cdfW(i)=cdf(pdW,Wsimgrid(i))-cdfW(i-1);
        ftW_i(floor(cdfW(i-1)*Nsim)+1:(floor(cdfW(i)*Nsim))) = Wsimgrid(i);
    end
    cdfW(ninter) = 1;
    %Individual fixed wages for simulation
    ftW_i(floor(cdfW(ninter-1)*Nsim)+1:Nsim) = Wsimgrid(ninter);
    clear muW sigW pd0W pdW cdfW Wsimgrid
    
    %Next, fixed types across kappa, lets do straight normal for now
    Ksimgrid=linspace(kapbarlow,kapbase,ninter);
    muK = 0.5*kapbarlow+0.5*kapbase;
    sigK = (kapbase-kapbarlow)/2;
    pd0K = makedist('normal','mu',muK,'sigma',sigK);
    pdK = truncate(pd0K,kapbarlow,kapbase) ;
    cdfK(1) = cdf(pdK,Ksimgrid(1));
    ftK_i(1:floor(cdfK(1)*Nsim)) = Ksimgrid(1);
    for i=2:ninter-1
        cdfK(i)=cdf(pdK,Ksimgrid(i))-cdfK(i-1);
        ftK_i(floor(cdfK(i-1)*Nsim)+1:(floor(cdfK(i)*Nsim))) = Ksimgrid(i);
    end
    cdfK(ninter) = 1;  
    %Individual fixed kappa types for simulation. Will be added/subtracted
        %from life-cycle kappas.
    ftK_i(floor(cdfK(ninter-1)*Nsim)+1:Nsim) = Ksimgrid(ninter);
    clear muK sigK pd0K pdK cdfK Ksimgrid
    
   %Last, lifecycle kappas. No interpolation here.
    Kysimgrid = kap(1,1:4);   
    Ky_i(1:floor(Nsim/nKy)) = Ksimgrid(1,1);
    Ky_inx(1:floor(Nsim/nKy)) = 1;
    for i=2:nKy-1
        Ky_i((i-1)*floor(Nsim/nKy)+1:i*floor(Nsim/nKy)) = Ksimgrid(1,i);
        Ky_inx((i-1)*floor(Nsim/nKy)+1:i*floor(Nsim/nKy)) = i;
    end
     Ky_i((nKy-1)*floor(Nsim/nKy)+1:nKy) = Ksimgrid(1,nKy);
     Ky_inx((nKy-1)*floor(Nsim/nKy)+1:nKy) = nKy;
     clear Kysimgrid
     
 %Now, interpolate policies
     %Add indices for types to make things easier
 gH0(1,1,:,:,:,:,:) = gH(1:4,:,:,:,:);
 gH0(1,2,:,:,:,:,:) = gH(5:8,:,:,:,:);
 gH0(2,1,:,:,:,:,:) = gH(9:12,:,:,:,:);
 gH0(2,2,:,:,:,:,:) = gH(13:16,:,:,:,:); 
 gS0(1,1,:,:,:,:,:) = gS(1:4,:,:,:,:);
 gS0(1,2,:,:,:,:,:) = gS(5:8,:,:,:,:);
 gS0(2,1,:,:,:,:,:) = gS(9:12,:,:,:,:);
 gS0(2,2,:,:,:,:,:) = gS(13:16,:,:,:,:); 
 gQ0(1,1,:,:,:,:,:) = gQ(1:4,:,:,:,:);
 gQ0(1,2,:,:,:,:,:) = gQ(5:8,:,:,:,:);
 gQ0(2,1,:,:,:,:,:) = gQ(9:12,:,:,:,:);
 gQ0(2,2,:,:,:,:,:) = gQ(13:16,:,:,:,:);  
 

 for i=1:Nsim
     for ie=1:nE
         for iy=1:nY
             for iz=1:nZ
                 for it=1:nT  %Room here to do cubic interp instead. 
                 ggH(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gH0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));
                 %gH1(1,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gH0(1,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %gH1(2,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gH0(2,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %ggH(i,it,ie,iy,iz) = interp1([ftW(1,1),ftW(1,nI)],squeeze(gH1(:,it,ie,iy,iz)),ftW_i(i));
                 ggS(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gS0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));
                 %gS1(1,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gS0(1,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %gS1(2,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gS0(2,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %ggS(i,it,ie,iy,iz) = interp1([ftW(1,1),ftW(1,nI)],squeeze(gS1(:,it,ie,iy,iz)),ftW_i(i));
                 ggQ(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gQ0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));
                 %gQ1(1,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gQ0(1,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %gQ1(2,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gQ0(2,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 %ggQ(i,it,ie,iy,iz) = interp1([ftW(1,1),ftW(1,nI)],squeeze(gQ1(:,it,ie,iy,iz)),ftW_i(i));                 
                 end
             end
         end
     end
 end
 
 clear gH0 gS0 gQ0 Ky_inx % gH1 gS1 gQ1
 
 %Set threshold for counting as NILF
    muS = mean(mean(mean(mean(mean(ggS)))));
    mnS = min(min(min(min(min(ggS)))));
    NiLFbar = 0.5*muS+0.5*mnS;
 
 %Set threshold for counting as FT
    muS = mean(mean(mean(mean(mean(ggH)))));
    mnS = min(min(min(min(min(ggH)))));
    FTbar = 0.5*muS+0.5*mnS;
    
  clear mnS muS
    
%-------------------------------------------------------------
%Types so far:
    %Age- deterministic and common to all
    %Ex-ante types + interpolated decision rules.
% Now need:
    %Business cycle shocks, common to all 
    %Spousal shocks- idiosyncratic
    %Initial experience    
%______________________________________________________________

%Draw business cycles
    %This is going to be for a hypothetical OLG structure
    %We have ages 25-64, 40 years
    %Need a "burn in" period of 38 years. 
    %I'll pick the business cycles to look something like data.
    %In data we start cohorts in 1930's. We start life @ 25. So, The first
    %macro states we need is 1955
    
    yqtime=(1955:0.25:2019);
    year=ones(size(yqtime));
    for iy=1955:2019
       for iq=1:4
          year((iy-1955)*4+iq) = iy;
       end
    end
    
zsim = ones(size(year,2),1);    
%Recession dates (1=expansion; 2=recession)
    %1957Q3-1958Q2
        strt = 1957.5;
        endd = 1958.25;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1960Q2-1961Q1
        strt = 1960.25;
        endd = 1961;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1969.Q4-1970Q4
        strt = 1969.75;
        endd = 1970.75;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1973Q4-1975Q1
        strt = 1973.75;
        endd = 1975;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1980Q1-1980Q3
        strt = 1980;
        endd = 1980.5;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1981Q3-1982Q4
        strt = 1981.5;
        endd = 1982.75;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %1990Q3-1991Q1
        strt = 1990.5;
        endd = 1991;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %2001Q1-2001Q4
        strt = 2001;
        endd = 2001.75;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    %2007Q4-2009Q2
        strt = 2007.75;
        endd = 2009.25;
        zsim((strt-1955)*4+1:(endd-1955)*4+1)=2;
    clear strt endd
    
%Initial experience type. Uniform for now.
exp_i= zeros(Nsim_i,Tsim);
rng(207);  %Same seed- area codes
esim0 = randi([1,floor(nE/3)],Nsim,1);
        for i=1:Ngen
            for ie=1:Nsim
               exp_i((i-1)*Nsim+ie,(i-1)*4+1) = egrid(esim0(ie,1));
            end
        end
 
 %Should also start with some in each state: unemployed + husbands across
 %states.
Hstat_i= zeros(Nsim_i,Tsim);    %Husband States
rng(585);  %Same seed
esim0 = rand(Nsim,1);
        for i=1:Ngen
            for ie=1:Nsim
               Hstat_i((i-1)*Nsim+ie,(i-1)*4+1) = 1; %Employed
               if (esim0(ie,1)<lamHbar(3)+lamHbar(2))
                Hstat_i((i-1)*Nsim+ie,(i-1)*4+1) = 2; %Recently Unemployed
               elseif esim0(ie,1)<lamHbar(3)
                Hstat_i((i-1)*Nsim+ie,(i-1)*4+1) = 3; %Currently Unemployed
               end
            end
        end
        
EmpI_i= zeros(Nsim_i, Tsim); UnempI_i= zeros(Nsim_i, Tsim);    %Wife State indicator
rng(612);  %Same seed
esim0 = rand(Nsim,1);
        for i=1:Ngen
            for ie=1:Nsim
                if (esim0(ie,1)<Nw_init)
                    EmpI_i((i-1)*Nsim+ie,(i-1)*4+1) = 0; %Unemployed
                    UnempI_i((i-1)*Nsim+ie,(i-1)*4+1) = 1; %Unemployed
                else
                    EmpI_i((i-1)*Nsim+ie,(i-1)*4+1) = 1; %Unemployed
                    UnempI_i((i-1)*Nsim+ie,(i-1)*4+1) = 0; %Unemployed                
               end
            end
        end        
 clear esim0;
 
 %Other draws: husband shock, female job loss/find shock for full life 
LHsim_i= zeros(Nsim_i,Tsim);
rng(812);  %Same seed
Lsim0 = rand(Nsim,Tsim_i);
        for i=1:Ngen
            for ie=1:Nsim
                for it=1:Tsim_i
                    LHsim_i((i-1)*Nsim+ie,(i-1)*4+it) = Lsim0(ie,it);
                end
            end
        end

LWsim_i= zeros(Nsim_i,Tsim);
rng(519);  %Same seed
Lsim0 = rand(Nsim,Tsim_i);
        for i=1:Ngen
            for ie=1:Nsim
                for i1=1:Tsim_i
                    LWsim_i((i-1)*Nsim+ie,(i-1)*4+it) = Lsim0(ie,it);
                end
            end
        end
    
clear Lsim0

Qtoss_i= zeros(Nsim_i,Tsim);
rng(911);  %Same seed
Qsim0 = rand(Nsim,Tsim_i);
        for i=1:Ngen
            for ie=1:Nsim
                for i1=1:Tsim_i
                    Qtoss_i((i-1)*Nsim+ie,(i-1)*4+it) = Qsim0(ie,it);
                end
            end
        end
clear Qsim0 

    save('Shocks_Types.mat','Nsim','Ngen','Tsim','Tsim_i','Nsim_t',...
         'ageT','ftW_i','ftK_i','Ky_i','NiLFbar','FTbar','yqtime','year',...
         'zsim','exp_i','Hstat_i','EmpI_i','UnempI_i','LHsim_i',...
         'LWsim_i','Qtoss_i') 

 else %If this is not the calibration loop, then we want everything the same except the policies.
    load('Shocks_Types.mat','Nsim','Ngen','Tsim','Tsim_i','Nsim_t',...
         'ageT','ftW_i','ftK_i','Ky_i','NiLFbar','FTbar','yqtime','year',...
         'zsim','exp_i','Hstat_i','EmpI_i','UnempI_i','LHsim_i',...
         'LWsim_i','Qtoss_i')    
     
     %Now, interpolate policies
     %Add indices for types to make things easier
 gH0(1,1,:,:,:,:,:) = gH(1:4,:,:,:,:);
 gH0(1,2,:,:,:,:,:) = gH(5:8,:,:,:,:);
 gH0(2,1,:,:,:,:,:) = gH(9:12,:,:,:,:);
 gH0(2,2,:,:,:,:,:) = gH(13:16,:,:,:,:); 
 gS0(1,1,:,:,:,:,:) = gS(1:4,:,:,:,:);
 gS0(1,2,:,:,:,:,:) = gS(5:8,:,:,:,:);
 gS0(2,1,:,:,:,:,:) = gS(9:12,:,:,:,:);
 gS0(2,2,:,:,:,:,:) = gS(13:16,:,:,:,:); 
 gQ0(1,1,:,:,:,:,:) = gQ(1:4,:,:,:,:);
 gQ0(1,2,:,:,:,:,:) = gQ(5:8,:,:,:,:);
 gQ0(2,1,:,:,:,:,:) = gQ(9:12,:,:,:,:);
 gQ0(2,2,:,:,:,:,:) = gQ(13:16,:,:,:,:);  
 

 for i=1:Nsim
     for ie=1:nE
         for iy=1:nY
             for iz=1:nZ
                 for it=1:nT  %Room here to do cubic interp instead. 
                 ggH(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gH0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));
                 ggS(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gS0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));
                 ggQ(i,it,ie,iy,iz) = interp2([ftW(1,1),ftW(1,nI)],[kapbarlow,kapbase],squeeze(gQ0(:,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i),ftW_i(i));              
                 end
             end
         end
     end
 end
 
 clear gH0 gS0 gQ0 % gH1 gS1 gQ1
     
end %ends the loop on whether it is calibration or experimential. 
 
 
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
% START SIMULATION
%___________________________________________________________________________________________

%Things to fill in. Zeros will mean unborn/dead
%States/outcomes
    Wage_i = zeros(Nsim_i, Tsim); Inc_w = zeros(Nsim_i, Tsim); Inc_hh(Nsim_i,Tsim); %Wife Wage, wife income, household income
    Hloss= zeros(Nsim_i, Tsim); % Husband job loss indicator
    Alive_i =zeros(Nsim_i, Tsim); Age_i= zeros(Nsim_i, Tsim); %Living indicator, Age group indicator
%Decisions
    Quit_i = zeros(Nsim_i, Tsim); Hours_i = zeros(Nsim_i, Tsim); Srch_i = zeros(Nsim_i, Tsim); %Quit, Hours, Search intensity
    FT_i = zeros(Nsim_i, Tsim); NiLF_i(Nsim_i,Tsim); NempI_i(Nsim_i,Tsim); %Indicator for Full Time, Indicator for Not in LF, not employed 
%Wife Flows
    EE_i = zeros(Nsim_i, Tsim); EU_i = zeros(Nsim_i, Tsim); UE = zeros(Nsim_i, Tsim);
    NE_i = zeros(Nsim_i, Tsim); EN_i = zeros(Nsim_i, Tsim); NU = zeros(Nsim_i, Tsim); UN = zeros(Nsim_i, Tsim);
    
    
%....and.... go!

for in=1:Ngen-1  %Loop over age year cohort generations. Toss last generation.
    for i=1:Nsim %Loop over individuals in that cohort- also is individual's fixed type
        %----------------------------------------------------------------------------------
        %First do period of birth
            %some things are set at stationary distribtution
            %All flows defined retrosp. as we see in data.
            ii = (in-1)*Nsim+i; %Individual's ID
            tt = (in-1)*4+1;    %Year-Quarter place in time
                Alive_i(ii,tt) = 1;
                Age_i(ii,tt) = ageT(1);
                e = exp_i(ii,tt);
                Wage_i(ii,tt) =  wage(ftW_i(i),e);
                iy = Hstat_i(ii,tt);
                iz = zsim(ii,tt);
                %Find quit policy
                Q= interp2(egrid,squeeze(ggQ(i,Age_i(ii,tt),:,iy,iz)),e);
                    if Q> Qtoss_i(ii,tt)
                        Quit_i(ii,tt)=1;
                        EmpI_i(ii,tt)=0;
                    end
            if (EmpI_i(ii,tt)==1 && Quit_i(ii,tt)<1) %Start employed, don't quit
                EmpI_i(ii,tt+1) = 1;
                %Find hours policy
                Hours_i(ii,tt) = interp(egrid,squeeze(ggH(i,Age_i(ii,tt),:,iy,iz)),e);
                if Hours_i(ii,tt) > FTbar
                    FT_i(ii,tt) = 1;
                end
                Inc_w(ii,tt) = Hours_i(ii,tt)*Wage_i(ii,tt);
                Inc_hh(ii,tt) = Inc_w(ii,tt)+wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                exp_i(ii,tt+1) = exp2(e,Hours_i(ii,tt));
            else  %Are Non-employed
                NempI_i(ii,tt) = 1;
                exp_i(ii,tt+1) = exp2(e,0);
                %Find search policy
                Srch_i(ii,tt) = interp(egrid,squeeze(ggS(i,Age_i(ii,tt),:,iy,iz)),e);
                if Srch_i(ii,tt)> NiLFbar
                    UnempI_i(ii,tt) = 1;
                else 
                    NiLF_i(ii,tt) = 1;
                end
                Inc_hh(ii,tt) = wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                %Draw job finding probability
                    if Srch_i(ii,tt)*findW(iz)>LWsim_i(ii,tt)
                        EmpI_i(ii,tt+1) = 1; %Found a job. Yay!  
                    else
                       NempI_i(ii,tt+1) = 1; %Did not find a job :(
                    end
            end
 %----------------------------------------------------------------------------------
 %  Now the rest of the life cycle with some retrospective variables.
        for it=2:Tsim_i
            tt = (in-1)*4+it;
            Alive_i(ii,tt) = 1;
            Age_i(ii,tt) = ageT(1);
            e = exp_i(ii,tt);
            Wage_i(ii,tt) =  wage(ftW_i(i),e);
            %Draw shocks
              %Wife job loss
                if LWsim_i(ii,tt)< lossW(iz)
                    EmpI_i(ii,tt) = 0;
                end                    
              %Husband Transition
                Hstat_i(ii,tt) = 1; %E
                if LHsim_i(ii,tt) < lamH(iz,iy,2)+ lamH(iz,iy,3) %Recently U
                    Hstat_i(ii,tt) = 2;
                end
                if LHsim_i(ii,tt) < lamH(iz,iy,3) %U
                    Hstat_i(ii,tt) = 3;
                    if Hstat_i(ii,tt-1) < 3
                       Hloss(ii,tt) = 1; 
                    end
                end                
                iy = Hstat_i(ii,tt);
                iz = zsim(ii,tt);
              %Find quit policy
                Q= interp2(egrid,squeeze(ggQ(i,Age_i(ii,tt),:,iy,iz)),e);
                    if Q> Qtoss_i(ii,tt)
                        Quit_i(ii,tt)=1;
                        EmpI_i(ii,tt)=0;
                    end                    
            if (EmpI_i(ii,tt)==1 && Quit_i(ii,tt)<1) %Start employed, don't quit
                EE_i(ii,tt) = EmpI(ii,tt-1);
                UE_i(ii,tt) = UnempI_i(ii,tt-1);
                NE_i(ii,tt) = NiLF_i(ii,tt-1);
                EmpI_i(ii,tt+1) = 1;
                %Find hours policy
                Hours_i(ii,tt) = interp(egrid,squeeze(ggH(i,Age_i(ii,tt),:,iy,iz)),e);
                if Hours_i(ii,tt) > FTbar
                    FT_i(ii,tt) = 1;
                end
                Inc_w(ii,tt) = Hours_i(ii,tt)*Wage_i(ii,tt);
                Inc_hh(ii,tt) = Inc_w(ii,tt)+wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                exp_i(ii,tt+1) = exp2(e,Hours_i(ii,tt));
            else  %Are Non-employed
                NempI_i(ii,tt) = 1;
                exp_i(ii,tt+1) = exp2(e,0);
                %Find search policy
                Srch_i(ii,tt) = interp(egrid,squeeze(ggS(i,Age_i(ii,tt),:,iy,iz)),e);
                if Srch_i(ii,tt)> NiLFbar
                    UnempI_i(ii,tt) = 1;
                    EU_i(ii,tt)=EmpI_i(ii,tt-1);
                    NU_i(ii,tt) = NiLF(ii,tt-1);
                else 
                    NiLF_i(ii,tt) = 1;
                    EN_i(ii,tt)=EmpI_i(ii,tt-1);
                    UN_i(ii,tt)=UnempI(ii,tt-1);
                end
                Inc_hh(ii,tt) = wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                %Draw job finding probability
                    if Srch_i(ii,tt)*findW(iz)>LWsim_i(ii,tt)
                        EmpI_i(ii,tt+1) = 1; %Found a job. Yay!                        
                    else
                       NempI_i(ii,tt+1) = 1; %Did not find a job :(
                    end
            end                      
                  
        end %it- lifecycle age
    end %i- indv. fixed type   
end %in- generation number
    







