% SimplerMod_vMay17.m
%clear all

%Set Paths
    MAINdir = 'C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab';
%Declare if this is the calibration script or the counterfactual script.
    calib= 1;
    %If counterfactual:
    %calib = 0;
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

load('paras')
load('gH.mat', 'gH')
load('gQ.mat', 'gQ')
load('gS.mat', 'gS')

if calib==1
    
%Set panel size
Nsim=100; %Number of individuals PER SINGLE YEAR BIRTH COHORT to simulate
Ngen = 2019-1955; %Number of single year birth cohorts.
ScratchT = (1973-1955)*tsize; %Where to start to begin timeseries in 1973 as in the data.
Tsim=Ngen*tsize;
Tsim_i= 39*tsize; %For the simple model: 25-39, 40-55, 56-64
Nsim_i=Nsim*Ngen;
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
    pd0W = makedist('normal','mu',muW,'sigma',sigW);
    pdW = truncate(pd0W,ftW(1,1),ftW(1,nI)) ;
    cdfW = cdf(pdW,Wsimgrid);
    rng(111);  %Same seed
    drw = rand(Nsim,1);
    ftW_i=ftW(1,nI).*ones(Nsim,1);
    for i=1:Nsim
        for j=ninter-1:-1:1
            if drw(i,1) < cdfW(1,j)
                ftW_i(i,1) = Wsimgrid(j);
            end
        end
    end

    clear muW sigW pd0W pdW cdfW Wsimgrid drw
    
    %Next, fixed types across kappa, lets do straight normal for now
        %can't use the cdf from now on or will be correlated with wage type
    kapbarlow = min(kapbar);
    kapbase = max(kapbar);
    Ksimgrid=linspace(kapbarlow,kapbase,ninter);
    muK = 0.5*kapbarlow+0.5*kapbase;
    sigK = (kapbase-kapbarlow)/2;
    pd0K = makedist('normal','mu',muK,'sigma',sigK);
    pdK = truncate(pd0K,kapbarlow,kapbase) ;
    cdfK = cdf(pdK,Ksimgrid);
    rng(222);  %Same seed
    drw = rand(Nsim,1);
    ftK_i=kapbase.*ones(Nsim,1);
    for i=1:Nsim
        for j=ninter-1:-1:1
            if drw(i,1) < cdfK(1,j)
                ftK_i(i,1) = Ksimgrid(j);
            end
        end
    end
    
    clear muK sigK pd0K pdK cdfK Ksimgrid drw
    
   %Last, lifecycle kappas. No interpolation here. Uniform
    Kysimgrid = kap(1,1:4);   
    cdfKy = zeros(1,4);
    cdfKy(1,1) = 1/size(Kysimgrid,2);
    for i=2:size(Kysimgrid,2)
        cdfKy(1,i) = cdfKy(1,i-1)+1/size(Kysimgrid,2);
    end
    rng(222);  %Same seed
    drw = rand(Nsim,1);
    Ky_i=kap(1,size(Kysimgrid,2)).*ones(Nsim,1);
    Ky_inx=size(Kysimgrid,2).*ones(Nsim,1);
    for i=1:Nsim
        for j=size(Kysimgrid,2)-1:-1:1
            if drw(i,1) < cdfKy(1,j)
                Ky_i(i,1) = Kysimgrid(j);
                Ky_inx(i,1) = j;
            end
        end
    end
    
     clear Kysimgrid drw cdfKy  
     
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
                 gH1(1,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gH0(1,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 gH1(2,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gH0(2,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 ggH(i,it,ie,iy,iz) = interp1([ftW(1,1),ftW(1,nI)],squeeze(gH1(:,it,ie,iy,iz)),ftW_i(i));

                 gS1(1,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gS0(1,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 gS1(2,it,ie,iy,iz) = interp1([kapbarlow,kapbase],squeeze(gS0(2,:,Ky_inx(i),it,ie,iy,iz)),ftK_i(i));
                 ggS(i,it,ie,iy,iz) = interp1([ftW(1,1),ftW(1,nI)],squeeze(gS1(:,it,ie,iy,iz)),ftW_i(i));

                 gQ1= (ftK_i(i)-kapbarlow)/(kapbase-kapbarlow).*squeeze(gQ0(:,1,Ky_inx(i),it,ie,iy,iz))+ (kapbase-ftK_i(i))/(kapbase-kapbarlow).*(squeeze(gQ0(:,2,Ky_inx(i),it,ie,iy,iz)));         
                 ggQ(i,it,ie,iy,iz) = gQ1(1)*(ftW_i(i)-ftW(1,1))/(ftW(1,nI)-ftW(1,1))+gQ1(2)*(ftW(1,nI)-ftW_i(i))/(ftW(1,nI)-ftW(1,1));         
                 %ggQ(i,it,ie,iy,iz) = max(min(ggQ(i,it,ie,iy,iz),1),0);                
                 end
             end
         end
     end
 end
 
 clear gH0 gS0 gQ0 gQ1 gH1 gS1 gQ1
 
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
esim0 = randi([floor(nE/3),nE],Nsim,1);
        for i=1:Ngen
            for ie=1:Nsim
               exp_i((i-1)*Nsim+ie,(i-1)*4+1) = egrid(esim0(ie,1));
               %exp_i((i-1)*Nsim+ie,(i-1)*4+1)=max(egrid);
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
               end
               if esim0(ie,1)<lamHbar(3)
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
                for it=1:Tsim_i
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
                for it=1:Tsim_i
                    Qtoss_i((i-1)*Nsim+ie,(i-1)*4+it) = Qsim0(ie,it);
                end
            end
        end
clear Qsim0 

    save('Shocks_Types.mat','Nsim','Ngen','Tsim','Tsim_i','Nsim_i',...
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
                 gQ1= squeeze(gQ0(:,1,Ky_inx(i),it,ie,iy,iz))+ (ftK_i(i)-kapbarlow)*(squeeze(gQ0(:,2,Ky_inx(i),it,ie,iy,iz)))-squeeze(gQ0(:,1,Ky_inx(i),it,ie,iy,iz))/(kapbase-kapbarlow);         
                 ggQ(i,it,ie,iy,iz) = gQ1(1)+ (ftW_i(i)-ftW(1,1))*(gQ1(2)-gQ1(1))/(ftW(1,nI)-ftW(1,1));         
                 ggQ(i,it,ie,iy,iz) = max(min(ggQ(i,it,ie,iy,iz),1),0);
                 end
             end
         end
     end
 end
 
 clear gH0 gS0 gQ0 kapbarlow kapbase % gH1 gS1 gQ1
     
end %ends the loop on whether it is calibration or experimential. 
 
 
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
% START SIMULATION
%___________________________________________________________________________________________

%Things to fill in. Zeros will mean unborn/dead
%States/outcomes
    Wage_i = zeros(Nsim_i, Tsim); Inc_w = zeros(Nsim_i, Tsim); Inc_hh=zeros(Nsim_i,Tsim); %Wife Wage, wife income, household income
    Hloss= zeros(Nsim_i, Tsim); Wloss= zeros(Nsim_i, Tsim);%  Husband job loss indicator, wife job loss
    Alive_i =zeros(Nsim_i, Tsim); Age_i= zeros(Nsim_i, Tsim); %Living indicator, Age group indicator
%Decisions
    Quit_i = zeros(Nsim_i, Tsim); Hours_i = zeros(Nsim_i, Tsim); Srch_i = zeros(Nsim_i, Tsim); %Quit, Hours, Search intensity
    FT_i = zeros(Nsim_i, Tsim); NiLF_i= zeros(Nsim_i,Tsim); NempI_i= zeros(Nsim_i,Tsim); %Indicator for Full Time, Indicator for Not in LF, not employed 
%Wife Flows
    EE_i = zeros(Nsim_i, Tsim); EU_i = zeros(Nsim_i, Tsim); UE_i = zeros(Nsim_i, Tsim);
    NE_i = zeros(Nsim_i, Tsim); EN_i = zeros(Nsim_i, Tsim); NU_i = zeros(Nsim_i, Tsim); UN_i = zeros(Nsim_i, Tsim);
    
    
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
                iz = zsim(tt);
                EmpI_i(ii,tt)=1;
               %Find quit policy
                ee=find(egrid<e,1,'last');                
                if ee<size(egrid,2)
                    Q = ggQ(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggQ(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                else
                    Q = ggQ(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                end   
                    if Q< Qtoss_i(ii,tt)
                        Quit_i(ii,tt)=1;
                        EmpI_i(ii,tt)=0;
                    end
            if (EmpI_i(ii,tt)>0 && Quit_i(ii,tt)<1) %Start employed, don't quit
                EmpI_i(ii,tt+1) = 1;
                %Find hours policy              
                if ee<size(egrid,2)
                    Hours_i(ii,tt) = ggH(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggH(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                else
                    Hours_i(ii,tt) = ggH(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                end                
                if Hours_i(ii,tt) > FTbar
                    FT_i(ii,tt) = 1;
                end
                Inc_w(ii,tt) = Hours_i(ii,tt)*Wage_i(ii,tt);
                Inc_hh(ii,tt) = Inc_w(ii,tt)+wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                exp_i(ii,tt+1) = exp2(e,Hours_i(ii,tt));
            else  %Are Non-employed
                NempI_i(ii,tt) = 1;
                Hours_i(ii,tt) = 0;
                exp_i(ii,tt+1) = exp2(e,0);
                %Find search policy                
                if ee<size(egrid,2)
                    Srch_i(ii,tt) = ggS(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggS(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                else
                    Srch_i(ii,tt) = ggS(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                end
                if Srch_i(ii,tt)> NiLFbar
                    UnempI_i(ii,tt) = 1;
                else 
                    NiLF_i(ii,tt) = 1;
                end
                Inc_hh(ii,tt) = wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                %Draw job finding probability
                    if (Srch_i(ii,tt)^(nu))*findW(iz)>LWsim_i(ii,tt)
                        EmpI_i(ii,tt+1) = 1; %Found a job. Yay!  
                    end
            end
            

            
 %----------------------------------------------------------------------------------
 %  Now the rest of the life cycle with some retrospective variables.
        for it=2:Tsim_i
            tt = (in-1)*4+it;
            if tt<Tsim+1
            Alive_i(ii,tt) = 1;
            Age_i(ii,tt) = ageT(it);
            e = exp_i(ii,tt);
            Wage_i(ii,tt) =  wage(ftW_i(i),e);
            %Draw shocks                
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
                iz = zsim(tt);
                
                if EmpI_i(ii,tt)>0
                    %Find quit policy
                     ee=find(egrid<e,1,'last');                
                     if ee<size(egrid,2)
                        Q = ggQ(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggQ(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                     else
                         Q = ggQ(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                     end                 
                    if Q< Qtoss_i(ii,tt)
                        Quit_i(ii,tt)=1;                                  
                    %Wife job loss
                    elseif LWsim_i(ii,tt)< lossW(iz)
                         Wloss(ii,tt) = 1;
                    end
                end
                           
            if (EmpI_i(ii,tt)>0 && Quit_i(ii,tt)<1 && Wloss(ii,tt)<1) %Start employed, don't quit or lose job
                EE_i(ii,tt) = EmpI_i(ii,tt-1);
                UE_i(ii,tt) = UnempI_i(ii,tt-1);
                NE_i(ii,tt) = NiLF_i(ii,tt-1);
                EmpI_i(ii,tt+1) = 1;
                %Find hours policy             
                if ee<size(egrid,2)
                    Hours_i(ii,tt) = ggH(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggH(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                else
                    Hours_i(ii,tt) = ggH(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                end  
                if Hours_i(ii,tt) > FTbar
                    FT_i(ii,tt) = 1;
                end
                Inc_w(ii,tt) = Hours_i(ii,tt)*Wage_i(ii,tt);
                Inc_hh(ii,tt) = Inc_w(ii,tt)+wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                exp_i(ii,tt+1) = exp2(e,Hours_i(ii,tt));
            else  %Are Non-employed
                EmpI_i(ii,tt+1) = 0;
                NempI_i(ii,tt) = 1;
                exp_i(ii,tt+1) = exp2(e,0);
                %Find search policy
                if ee<size(egrid,2)
                    Srch_i(ii,tt) = ggS(i,Age_i(ii,tt),ee,iy,iz)*(e-egrid(1,ee))/(egrid(1,ee+1)-egrid(1,ee))+ggS(i,Age_i(ii,tt),ee+1,iy,iz)*(egrid(1,ee+1)-e)/(egrid(1,ee+1)-egrid(1,ee));
                else
                    Srch_i(ii,tt) = ggS(i,Age_i(ii,tt),size(egrid,2),iy,iz);
                end
                if Srch_i(ii,tt)> NiLFbar
                    UnempI_i(ii,tt) = 1;
                    EU_i(ii,tt)=EmpI_i(ii,tt);
                    NU_i(ii,tt) = NiLF_i(ii,tt-1);
                else 
                    NiLF_i(ii,tt) = 1;
                    EN_i(ii,tt)=EmpI_i(ii,tt);
                    UN_i(ii,tt)=UnempI_i(ii,tt-1);
                end
                EmpI_i(ii,tt) = 0;
                Inc_hh(ii,tt) = wageH(Age_i(ii,tt))*BCwage(iz)*ym(iy);
                %Draw job finding probability
                    if (Srch_i(ii,tt)^(nu))*findW(iz)>LWsim_i(ii,tt)
                        EmpI_i(ii,tt+1) = 1; %Found a job. Yay!                        
                    end
            end                      
            end %don't go over grid
        end %it- lifecycle age
    end %i- indv. fixed type   
end %in- generation number

%********************************************************************************************************************************************************************************************
%CALCULATE AGGREGATE STATISTICS
%   1) By Age
%   2) By Career Type
%   3) Business Cycle Stats
%---------------------------------------------------------------------------------------------%---------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------%---------------------------------------------------------------------------------------------

%Cross Section things
Erate = zeros(1,3);
Nrate = zeros(1,3);
Urate = zeros(1,3);
EErate = zeros(1,3);
EUrate = zeros(1,3);
ENrate = zeros(1,3);
UErate = zeros(1,3);
NonEErate = zeros(1,3);
Qrate = zeros(1,3);
jLossrate = zeros(1,3);
mHours = zeros(1,3);
mWage_E = zeros(1,3);
mWage_U = zeros(1,3);
mWage_N = zeros(1,3);
mSearch_U = zeros(1,3);
mSearch_N = zeros(1,3);
mExp_U = zeros(1,3);
mExp_N = zeros(1,3);
mExp_E = zeros(1,3);
pop = zeros(1,3);
HstatE_E = zeros(1,3);
HstatE_U = zeros(1,3);
HstatE_N = zeros(1,3);
HstatS_E = zeros(1,3);
HstatS_U = zeros(1,3);
HstatS_N = zeros(1,3);
Ktype_E = zeros(1,3);
Ktype_U = zeros(1,3);
Ktype_N = zeros(1,3);

%Panel things


%Business Cycle Things
H_Ureccess = zeros(1,Tsim);
H_Uexp = zeros(1,Tsim);
H_Ereccess = zeros(1,Tsim);
H_Eexp = zeros(1,Tsim);


W_Ureccess = zeros(1,Tsim);
W_Uexp = zeros(1,Tsim);
W_Ereccess = zeros(1,Tsim);
W_Eexp = zeros(1,Tsim);
W_NEreccess = zeros(1,Tsim);
W_NEexp = zeros(1,Tsim);
W_NonEreccess = zeros(1,Tsim);
W_NonEexp = zeros(1,Tsim);
W_EEreccess = zeros(1,Tsim);
W_EEexp = zeros(1,Tsim);
W_ENreccess = zeros(1,Tsim);
W_ENexp = zeros(1,Tsim);
W_EUreccess = zeros(1,Tsim);
W_EUexp = zeros(1,Tsim);
W_NonEEreccess = zeros(1,Tsim);
W_NonEEexp = zeros(1,Tsim);

%------------------------------
%1) Summary statistics- by age
%-------------------------------
    %By age
    for in=1:Ngen-1  %Loop over age year cohort generations. Toss last generation.
        for i=1:Nsim %Loop over individuals in that cohort- also is individual's fixed type   
            ii = (in-1)*Nsim+i; %Individual's ID
        %Young
        a=1;    %age group
        for it=1:14*tsize
            tt = (in-1)*4+it;    %Year-Quarter place in time 
            if tt<Tsim+1 
            pop(a) = pop(a)+1;
            %Employed
                Erate(a) = Erate(a)+EmpI_i(ii,tt);
                mHours(a) = mHours(a)+Hours_i(ii,tt)*EmpI_i(ii,tt);
                mWage_E(a) = mWage_E(a) + Wage_i(ii,tt)*EmpI_i(ii,tt);
                mExp_E(a) = mExp_E(a)+exp_i(ii,tt)*EmpI_i(ii,tt);
                if Hstat_i(ii,tt)==1
                    HstatE_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                end
                Ktype_E(a)= Ktype_E(a)+ (Ky_i(i)+ftK_i(i))*EmpI_i(ii,tt);                
            %Unemp
                Urate(a) = Urate(a)+UnempI_i(ii,tt);
                mSearch_U(a) = mSearch_U(a)+Srch_i(ii,tt)*UnempI_i(ii,tt);
                mExp_U(a) = mExp_U(a)+exp_i(ii,tt)*UnempI_i(ii,tt);  
                mWage_U(a) = mWage_U(a) + Wage_i(ii,tt)*UnempI_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                end
                Ktype_U(a)= Ktype_U(a)+ (Ky_i(i)+ftK_i(i))*UnempI_i(ii,tt);             
            %NiLF
                Nrate(a) = Nrate(a)+NiLF_i(ii,tt);
                mSearch_N(a) = mSearch_N(a)+Srch_i(ii,tt)*NiLF_i(ii,tt);
                mExp_N(a) = mExp_N(a)+exp_i(ii,tt)*NiLF_i(ii,tt);   
                mWage_N(a) = mWage_N(a) + Wage_i(ii,tt)*NiLF_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                end
                Ktype_N(a)= Ktype_N(a)+ (Ky_i(i)+ftK_i(i))*NiLF_i(ii,tt);                  
             %All types
                EErate(a) = EErate(a)+EE_i(ii,tt);
                EUrate(a) = EUrate(a)+EU_i(ii,tt);
                ENrate(a) = ENrate(a)+EN_i(ii,tt);
                UErate(a) = UErate(a)+UE_i(ii,tt);
                NonEErate(a) = NonEErate(a)+NE_i(ii,tt)+UE_i(ii,tt);
                Qrate(a) = Qrate(a)+Quit_i(ii,tt);
                jLossrate(a) = jLossrate(a)+Wloss(ii,tt);                
            end
        end
    
        %Prime
        a=2;    %age group
        for it=14*tsize+1:14*tsize+14*tsize
            tt = (in-1)*4+it;    %Year-Quarter place in time 
            if tt<Tsim+1
            pop(a) = pop(a)+1;
            %Employed
                Erate(a) = Erate(a)+EmpI_i(ii,tt);
                mHours(a) = mHours(a)+Hours_i(ii,tt)*EmpI_i(ii,tt);
                mWage_E(a) = mWage_E(a) + Wage_i(ii,tt)*EmpI_i(ii,tt);
                mExp_E(a) = mExp_E(a)+exp_i(ii,tt)*EmpI_i(ii,tt);
                if Hstat_i(ii,tt)==1
                    HstatE_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                end
                Ktype_E(a)= Ktype_E(a)+ (Ky_i(i)+ftK_i(i))*EmpI_i(ii,tt);                
            %Unemp
                Urate(a) = Urate(a)+UnempI_i(ii,tt);
                mSearch_U(a) = mSearch_U(a)+Srch_i(ii,tt)*UnempI_i(ii,tt);
                mExp_U(a) = mExp_U(a)+exp_i(ii,tt)*UnempI_i(ii,tt);  
                mWage_U(a) = mWage_U(a) + Wage_i(ii,tt)*UnempI_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                end
                Ktype_U(a)= Ktype_U(a)+ (Ky_i(i)+ftK_i(i))*UnempI_i(ii,tt);             
            %NiLF
                Nrate(a) = Nrate(a)+NiLF_i(ii,tt);
                mSearch_N(a) = mSearch_N(a)+Srch_i(ii,tt)*NiLF_i(ii,tt);
                mExp_N(a) = mExp_N(a)+exp_i(ii,tt)*NiLF_i(ii,tt);   
                mWage_N(a) = mWage_N(a) + Wage_i(ii,tt)*NiLF_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                end
                Ktype_N(a)= Ktype_N(a)+ (Ky_i(i)+ftK_i(i))*NiLF_i(ii,tt);                  
             %All types
                EErate(a) = EErate(a)+EE_i(ii,tt);
                EUrate(a) = EUrate(a)+EU_i(ii,tt);
                ENrate(a) = ENrate(a)+EN_i(ii,tt);
                UErate(a) = UErate(a)+UE_i(ii,tt);
                NonEErate(a) = NonEErate(a)+NE_i(ii,tt)+UE_i(ii,tt);
                Qrate(a) = Qrate(a)+Quit_i(ii,tt);
                jLossrate(a) = jLossrate(a)+Wloss(ii,tt);                
            end
        end
        
        %Old
        a=3;    %age group
             for it=14*tsize+14*tsize+1:Tsim_i
            tt = (in-1)*4+it;    %Year-Quarter place in time 
            if tt<Tsim+1
            pop(a) = pop(a)+1;
            %Employed
                Erate(a) = Erate(a)+EmpI_i(ii,tt);
                mHours(a) = mHours(a)+Hours_i(ii,tt)*EmpI_i(ii,tt);
                mWage_E(a) = mWage_E(a) + Wage_i(ii,tt)*EmpI_i(ii,tt);
                mExp_E(a) = mExp_E(a)+exp_i(ii,tt)*EmpI_i(ii,tt);
                if Hstat_i(ii,tt)==1
                    HstatE_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_E(a)= HstatE_E(a)+ EmpI_i(ii,tt);
                end
                Ktype_E(a)= Ktype_E(a)+ (Ky_i(i)+ftK_i(i))*EmpI_i(ii,tt);                
            %Unemp
                Urate(a) = Urate(a)+UnempI_i(ii,tt);
                mSearch_U(a) = mSearch_U(a)+Srch_i(ii,tt)*UnempI_i(ii,tt);
                mExp_U(a) = mExp_U(a)+exp_i(ii,tt)*UnempI_i(ii,tt);  
                mWage_U(a) = mWage_U(a) + Wage_i(ii,tt)*UnempI_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_U(a)= HstatE_U(a)+ UnempI_i(ii,tt);
                end
                Ktype_U(a)= Ktype_U(a)+ (Ky_i(i)+ftK_i(i))*UnempI_i(ii,tt);             
            %NiLF
                Nrate(a) = Nrate(a)+NiLF_i(ii,tt);
                mSearch_N(a) = mSearch_N(a)+Srch_i(ii,tt)*NiLF_i(ii,tt);
                mExp_N(a) = mExp_N(a)+exp_i(ii,tt)*NiLF_i(ii,tt);   
                mWage_N(a) = mWage_N(a) + Wage_i(ii,tt)*NiLF_i(ii,tt);
               if Hstat_i(ii,tt)==1
                    HstatE_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                elseif Hstat_i(ii,tt)==2
                    HstatS_N(a)= HstatE_N(a)+ NiLF_i(ii,tt);
                end
                Ktype_N(a)= Ktype_N(a)+ (Ky_i(i)+ftK_i(i))*NiLF_i(ii,tt);                  
             %All types
                EErate(a) = EErate(a)+EE_i(ii,tt);
                EUrate(a) = EUrate(a)+EU_i(ii,tt);
                ENrate(a) = ENrate(a)+EN_i(ii,tt);
                UErate(a) = UErate(a)+UE_i(ii,tt);
                NonEErate(a) = NonEErate(a)+NE_i(ii,tt)+UE_i(ii,tt);
                Qrate(a) = Qrate(a)+Quit_i(ii,tt);
                jLossrate(a) = jLossrate(a)+Wloss(ii,tt);                
            end
             end
        
        end
    end

%Convert to rates and averages
mHours = mHours./Erate;
mWage_E = mWage_E./Erate;
mExp_E = mExp_E./Erate;
HstatE_E = Hstat_E./Erate;
Ktype_E = Ktype_E./Erate;
EErate = EErate./Erate;
ENrate = ENrate./Erate;
EUrate = EUrate./Erate;
Qrate = Qrate./Erate;
jLossrate = jLossrate./Erate;

mSearch_U = mSearch_U./Urate;
mExp_U = mExp_U./Urate;
mWage_U = mWage_U./Urate;
HstatE_U = Hstat_U./Urate;
Ktype_U = Ktype_U./Urate;
UErate = UErate./Urate;

mSearch_N = mSearch_N./Nrate;
mExp_N = mExp_N./Nrate;
mWage_N = mWage_N./Nrate;
HstatE_N = Hstat_N./Nrate;
Ktype_N = Ktype_N./Nrate;
NErate = NErate./Nrate;

Erate = Erate./pop;
Urate = Urate./pop;
Nrate = Nrate./pop;

header = {'' , 'Young', 'Middle', 'Old'}; 

diary = xlswrite('SimulStats.xls', header, 'CrossSection', 'A2'); 
diary = xlswrite('SimulStats.xls', header, 'CrossSection', 'E2'); 
diary = xlswrite('SimulStats.xls', header, 'CrossSection', 'I2'); 

t=1;

EmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mHours(1,t), mWage_E(1,t), mExp_E(1,t), HstatE_E(1,t), Ktype_E(1,t) };        
UnEmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_U(1,t), mWage_U(1,t), mExp_U(1,t), HstatE_U(1,t), Ktype_U(1,t) };        
NiLFXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_N(1,t), mWage_N(1,t), mExp_N(1,t), HstatE_N(1,t), Ktype_N(1,t) };       
FlowXstats =  {'EE rate', 'EN rate', 'EU rate', 'Quit rate', 'Job Loss Rate', 'UE rate', 'NE rate';...
            EErate(1,t), ENrate(1,t), EUrate(1,t), Qrate(1,t), jLossrate(1,t), UErate(1,t), NErate(1,t) };       
            
diary = xlswrite('SimulStats.xls', EmpXstats', 'CrossSection', 'A3');
diary = xlswrite('SimulStats.xls', UnEmpXstats(2,:)', 'CrossSection', 'F3');
diary = xlswrite('SimulStats.xls', NiLFXstats(2,:)', 'CrossSection', 'J3');
diary = xlswrite('SimulStats.xls', FlowXstats', 'CrossSection', 'A9');

t=2;

EmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mHours(1,t), mWage_E(1,t), mExp_E(1,t), HstatE_E(1,t), Ktype_E(1,t) };        
UnEmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_U(1,t), mWage_U(1,t), mExp_U(1,t), HstatE_U(1,t), Ktype_U(1,t) };        
NiLFXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_N(1,t), mWage_N(1,t), mExp_N(1,t), HstatE_N(1,t), Ktype_N(1,t) };       
FlowXstats =  {'EE rate', 'EN rate', 'EU rate', 'Quit rate', 'Job Loss Rate', 'UE rate', 'NE rate';...
            EErate(1,t), ENrate(1,t), EUrate(1,t), Qrate(1,t), jLossrate(1,t), UErate(1,t), NErate(1,t) };       
            
diary = xlswrite('SimulStats.xls', EmpXstats(2,:)', 'CrossSection', 'C3');
diary = xlswrite('SimulStats.xls', UnEmpXstats(2,:)', 'CrossSection', 'G3');
diary = xlswrite('SimulStats.xls', NiLFXstats(2,:)', 'CrossSection', 'K3');
diary = xlswrite('SimulStats.xls', FlowXstats(2,:)', 'CrossSection', 'C9');

t=3;

EmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mHours(1,t), mWage_E(1,t), mExp_E(1,t), HstatE_E(1,t), Ktype_E(1,t) };        
UnEmpXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_U(1,t), mWage_U(1,t), mExp_U(1,t), HstatE_U(1,t), Ktype_U(1,t) };        
NiLFXstats = {'Hours/Search', 'Wages', 'Experience', 'Husband Employed', 'Kappa';...
            mSearch_N(1,t), mWage_N(1,t), mExp_N(1,t), HstatE_N(1,t), Ktype_N(1,t) };       
FlowXstats =  {'EE rate', 'EN rate', 'EU rate', 'Quit rate', 'Job Loss Rate', 'UE rate', 'NE rate';...
            EErate(1,t), ENrate(1,t), EUrate(1,t), Qrate(1,t), jLossrate(1,t), UErate(1,t), NErate(1,t) };       
            
diary = xlswrite('SimulStats.xls', EmpXstats(2,:)', 'CrossSection', 'D3');
diary = xlswrite('SimulStats.xls', UnEmpXstats(2,:)', 'CrossSection', 'H3');
diary = xlswrite('SimulStats.xls', NiLFXstats(2,:)', 'CrossSection', 'L3');
diary = xlswrite('SimulStats.xls', FlowXstats(2,:)', 'CrossSection', 'D9');

%-------------------------------------
%2) Summary statistics- by career type
%-------------------------------------

%1) find career type
ageFT =zeros(Nsim_i,3);
ageEmp =zeros(Nsim_i,3);
ageAlive = zeros(Nsim_i,3);
cyclewoman = zeros(Nsim_i,1);
careerwoman = zeros(Nsim_i,1);
PTwoman = zeros(Nsim_i,1);
NiLFwoman = zeros(Nsim_i,1);

    for in=1:Ngen-1  %Loop over age year cohort generations. Toss last generation.
        for i=1:Nsim %Loop over individuals in that cohort- also is individual's fixed type   
            ii = (in-1)*Nsim+i; %Individual's ID
        for it=1:Tsim_i
            tt = (in-1)*4+it;    %Year-Quarter place in time 
            if tt<Tsim+1 && Age_i(ii,it)>0
                   ageEmp(ii,Age_i(ii,it)) = ageEmp(ii,Age_i(ii,it))+ EmpI_i(ii,tt); 
                   ageAlive(ii,Age_i(ii,it)) = ageAlive(ii,Age_i(ii,it))+ Alive_i(ii,tt);
                   ageFT(ii,Age_i(ii,it)) = ageFT(ii,Age_i(ii,it)) + FT_i(ii,tt);
            end
        end
        ageEmpR(ii,:) = ageEmp(ii,:)./ageAlive(ii,:);
        EmpR_i(ii) = mean(ageEmp(ii,:)./ageAlive(ii,:));
        if EmpR_i(ii)>0.7
            if sum(ageFT(ii,:))/sum(ageEmp(ii,:)) > 0.8
                careerwoman(ii) = 1;
            else
                PTwoman(ii) = 1;
            end
        elseif ageEmpR(ii,1)<0.5*ageEmpR(ii,2)
            cyclewoman(ii) = 1;
        elseif EmpR_i(ii)<0.5
            NiLFwoman(ii) =1;
        end
        end         
    end
 
ptWtot = sum(PTwoman);
cycleWtot = sum(cyclewoman);
careerWtot = sum(careerwoman);
nilfWtot = sum(NiLFwoman);

ptWrate = ptWtot/(ptWtot+cycleWtot+careerWtot+nilfWtot)
cycleWrate = cycleWtot/(ptWtot+cycleWtot+careerWtot+nilfWtot)
careerWrate = careerWtot/(ptWtot+cycleWtot+careerWtot+nilfWtot)
nilfWrate = nilfWtot/(ptWtot+cycleWtot+careerWtot+nilfWtot)

%-------------------------------------
%2) Summary statistics- cyclicality averages and stD
%-------------------------------------

%Generate indicators for recession and non to do easy matrix multiplication
%stats
    RecI = zeros(1,Tsim-ScratchT+1);
    ExpI = ones(1,Tsim-ScratchT+1);
    for t=1:Tsim-ScratchT+1
        if zsim(t+ScratchT-1,1)==2
            RecI(1,t) = 1 ; 
            ExpI(1,t) = 1 ; 
        end    
    end

%Need some more general indicators 
    H_EmpI=zeros(Nsim_i,Tsim);
    W_NonEmpI=zeros(Nsim_i,Tsim);
    W_NonEE=zeros(Nsim_i,Tsim);
    W_ENonE=zeros(Nsim_i,Tsim);

    for in=1:Ngen-1  %Loop over age year cohort generations. Toss last generation.
        for i=1:Nsim %Loop over individuals in that cohort- also is individual's fixed type   
            ii = (in-1)*Nsim+i; %Individual's ID
        for it=1:Tsim_i
            tt = (in-1)*4+it;    %Year-Quarter place in time 
            if tt<Tsim+1
                if (Hstat_i(ii,tt)==1 || Hstat_i(ii,tt)==2)
                    H_EmpI(ii,tt)=1; 
                end   
                if (NiLF_i(ii,tt)==1 || UnempI_i(ii,tt)==1)
                    W_NonEmpI(ii,tt)=1; 
                end  
                if (NE_i(ii,tt)==1 || UE_i(ii,tt)==1)
                    W_NonEE(ii,tt)=1; 
                end      
                if (EN_i(ii,tt)==1 || EU_i(ii,tt)==1)
                    W_ENonE(ii,tt)=1; 
                end                 
            end
        end
        end
    end
    
Pop_t = sum(Alive_i(:,ScratchT:Tsim),1);

%Husband Stuff
    H_Ereccess_t = (sum(H_EmpI(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    mH_Erecess = sum(H_Ereccess_t)/sum(RecI);
    H_Eexp_t = (sum(H_EmpI(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    mH_Eexp = sum(H_Eexp_t)/sum(ExpI);
    varH_Emp = var(sum(H_EmpI(:,ScratchT:Tsim),1)./Pop_t)*100;     

%Wife stuff
    W_Ureccess_t = (sum(UnempI_i(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    mW_Urecess = sum(W_Ureccess_t)/sum(RecI);
    W_Uexp_t = (sum(UnempI_i(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    mW_Uexp = sum(W_Uexp_t)/sum(ExpI);
    
    W_Erecess_t = (sum(EmpI_i(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    mW_Erecess = sum(W_Erecess_t)/sum(RecI);   
    W_Eexp_t = (sum(EmpI_i(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    mW_Eexp = sum(W_Eexp_t)/sum(ExpI);
    varW_Emp = var(sum(EmpI_i(:,ScratchT:Tsim),1)./Pop_t)*100; 

    W_Nrecess_t = (sum(NiLF_i(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    mW_Nrecess = sum(W_Nrecess_t)/sum(RecI);   
    W_Nexp_t = (sum(NiLF_i(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    mW_Nexp = sum(W_Nexp_t)/sum(ExpI);

    W_NEmprecess_t = (sum(W_NonEmpI(:,ScratchT:Tsim),1).*RecI(1,:))./sum(NiLF_i(:,ScratchT:Tsim),1);
    mW_NEmprecess = sum(W_NEmprecess_t)/sum(RecI);   
    W_NEmpexp_t = (sum(W_NonEmpI(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(NiLF_i(:,ScratchT:Tsim),1);
    mW_NEmpexp = sum(W_NEmpexp_t)/sum(ExpI);
    
    EErecess_t = (sum(EE_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_EErecess = sum(EErecess_t)/sum(RecI);   
    EEexp_t = (sum(EE_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_EEexp = sum(EEexp_t)/sum(ExpI);
    var_EE = var(sum(EE_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100; 
    
    ENrecess_t = (sum(EN_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_ENrecess = sum(ENrecess_t)/sum(RecI);   
    ENexp_t = (sum(EN_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_ENexp = sum(ENexp_t)/sum(ExpI);
    var_EN = var(sum(EN_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100; 
    
    EUrecess_t = (sum(EU_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_EUrecess = sum(EUrecess_t)/sum(RecI);   
    EUexp_t = (sum(EU_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_EUexp = sum(EUexp_t)/sum(ExpI);
    var_EU = var(sum(EU_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100; 
    
    ENonErecess_t = (sum(W_ENonE(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_ENonErecess = sum(ENonErecess_t)/sum(RecI);   
    ENonEexp_t = (sum(W_ENonE(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_ENonEexp = sum(ENonEexp_t)/sum(ExpI);
    var_ENonE = var(sum(W_ENonE(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100;     
    
    NonEErecess_t = (sum(W_NonEE(:,ScratchT:Tsim),1).*RecI(1,:))./sum(W_NonEmpI(:,ScratchT:Tsim),1);
    m_NonEErecess = sum(NonEErecess_t)/sum(RecI);   
    NonEEexp_t = (sum(W_NonEE(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(W_NonEmpI(:,ScratchT:Tsim),1);
    m_NonEEexp = sum(NonEEexp_t)/sum(ExpI);
    var_NonEE = var(sum(W_NonEE(:,ScratchT:Tsim),1)./sum(W_NonEmpI(:,ScratchT:Tsim),1))*100;  

    Hrs_recess_t = (sum(Hours_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Hrs_recess = sum(Hrs_recess_t)/sum(RecI);   
    Hrs_exp_t = (sum(Hours_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Hrs_exp = sum(Hrs_exp_t)/sum(ExpI);
    var_Hrs = var(sum(Hours_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100;    
    
    Quit_recess_t = (sum(Quit_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Quit_recess = sum(Quit_recess_t)/sum(RecI);   
    Quit_exp_t = (sum(Quit_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Quit_exp = sum(Quit_exp_t)/sum(ExpI);
    var_Quit = var(sum(Quit_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100;       
    
    Wage_recess_t = (sum(Wage_i(:,ScratchT:Tsim).*EmpI_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Wage_recess = sum(Wage_recess_t)/sum(RecI);   
    Wage_exp_t = (sum(Wage_i(:,ScratchT:Tsim).*EmpI_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(EmpI_i(:,ScratchT:Tsim),1);
    m_Wage_exp = sum(Wage_exp_t)/sum(ExpI);
    var_Wage = var(sum(Wage_i(:,ScratchT:Tsim).*EmpI_i(:,ScratchT:Tsim),1)./sum(EmpI_i(:,ScratchT:Tsim),1))*100;  
    
    Srch_recess_t = (sum(Srch_i(:,ScratchT:Tsim),1).*RecI(1,:))./sum(UnempI_i(:,ScratchT:Tsim),1);
    m_Srch_recess = sum(Srch_recess_t)/sum(RecI);   
    Srch_exp_t = (sum(Srch_i(:,ScratchT:Tsim),1).*ExpI(1,:))./sum(UnempI_i(:,ScratchT:Tsim),1);
    m_Srch_exp = sum(Srch_exp_t)/sum(ExpI);
    var_Srch = var(sum(Srch_i(:,ScratchT:Tsim),1)./sum(UnempI_i(:,ScratchT:Tsim),1))*100;   
    
    IncHH_recess_t = (sum(Inc_hh(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    m_IncHH_recess = sum(IncHH_recess_t)/sum(RecI);   
    IncHH_exp_t = (sum(Inc_hh(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    m_IncHH_exp = sum(IncHH_exp_t)/sum(ExpI);
    var_IncHH = var(sum(Inc_hh(:,ScratchT:Tsim),1)./Pop_t)*100;   
    
    IncW_recess_t = (sum(Inc_w(:,ScratchT:Tsim),1).*RecI(1,:))./Pop_t;
    m_IncW_recess = sum(IncW_recess_t)/sum(RecI);   
    IncW_exp_t = (sum(Inc_w(:,ScratchT:Tsim),1).*ExpI(1,:))./Pop_t;
    m_IncW_exp = sum(IncW_exp_t)/sum(ExpI);
    var_IncW = var(sum(Inc_w(:,ScratchT:Tsim),1)./Pop_t)*100;       
    
    Wshare_recess_t = ((sum(Inc_w(:,ScratchT:Tsim),1)./sum(Inc_hh(:,ScratchT:Tsim),1)).*RecI(1,:));
    m_Wshare_recess = m_IncW_recess/m_IncHH_recess  ;
    Wshare_exp_t = ((sum(Inc_w(:,ScratchT:Tsim),1)./sum(Inc_hh(:,ScratchT:Tsim),1)).*ExpI(1,:));
    m_Wshare_exp = m_IncW_exp/m_IncHH_exp;
    var_Wshare = var(((sum(Inc_w(:,ScratchT:Tsim),1)./sum(Inc_hh(:,ScratchT:Tsim),1))))*100;   
    
Rec_means = {'H_Emp', 'Emp', 'UnEmp', 'NiLF', 'NonEmp', 'EE', 'EU', 'EN', 'EnonE', 'NonEE', 'Hours', 'Quit', 'Search', 'Wage', 'Household Income', 'Wife Income', 'Wife Share of Inc' ;...
            mH_Erecess , mW_Erecess, mW_Urecess, mW_Nrecess, mW_NEmprecess, m_EErecess, m_ENrecess, m_EUrecess, m_ENonErecess, m_NonEErecess, m_Hrs_recess, m_Quit_recess, m_Srch_recess,...
            m_Wage_recess, m_IncHH_recess, m_IncW_recess, m_Wshare_recess  };
        
Exp_means = {'H_Emp', 'Emp', 'UnEmp', 'NiLF', 'NonEmp', 'EE', 'EU', 'EN', 'EnonE', 'NonEE', 'Hours', 'Quit', 'Search', 'Wage', 'Household Income', 'Wife Income', 'Wife Share of Inc' ;...
            mH_Eexp , mW_Eexp, mW_Uexp, mW_Nexp, mW_NEmpexp, m_EEexp, m_ENexp, m_EUexp, m_ENonEexp, m_NonEEexp, m_Hrs_exp, m_Quit_exp, m_Srch_exp,...
            m_Wage_exp, m_IncHH_exp, m_IncW_exp, m_Wshare_exp  };
        
Variance = {'H_Emp', 'Emp', 'UnEmp', 'NiLF', 'NonEmp', 'EE', 'EU', 'EN', 'EnonE', 'NonEE', 'Hours', 'Quit', 'Search', 'Wage', 'Household Income', 'Wife Income', 'Wife Share of Inc' ;...
            varH_Emp , varW_Emp, 0, 0, 0, var_EE, var_EN, var_EU, var_ENonE, var_NonEE, var_Hrs, var_Quit, var_Srch, var_Wage, var_IncHH, var_IncW, var_Wshare};
     
        
diary = xlswrite('SimulStats.xls', Rec_means', 'Cycle', 'A2');
diary = xlswrite('SimulStats.xls', Exp_means(2,:)', 'Cycle', 'C2');
diary = xlswrite('SimulStats.xls', Variance(2,:)', 'Cycle', 'D2');
 
%---------------------------------------------------------------------------------------------%---------------------------------------------------------------------------------------------
%EXPORT PANEL DATA
%---------------------------------------------------------------------------------------------%---------------------------------------------------------------------------------------------


%plot(egrid,squeeze(ggQ(1,1,:,1,1)))