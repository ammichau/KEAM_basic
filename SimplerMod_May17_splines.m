% SimplerMod_May17_splines.m
%clear all

%Set Paths
%    MAINdir = 'C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab';
%    SOLUTIONdir = 'C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab\Solution\RoE_incr';   

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

%Set Parameters
%--------------------------------------------------------------------------------
global beta crra eta mu r xi z_h tau_wf alpha_e delta_e psi phi_c gam_e ybar_h  alpha_h

 beta = 0.95;     % discount factor
 crra = 2 ;        % curvature of cons utility
 eta = 1.2;        % curvature of leisure
 mu = 0.4 ;        % utility weight leisure
 r = 0.0000;        % interest rate
 alpha_e = 0.03;    %Weight on hours in returns to experience
 delta_e = 0.01;   % depreciation human capital
 psi = 0.8;        % curvature on hours in building experience
 %Together should satisfy: (delta_h/alpha_h)^1/psi = 0.25
 lhome = 0.1;       %Time cost when at home.
 z_h = 0.45;        % homeproduction productivity
 ybar_h = 0.2;     %home production constant 0.2
 alpha_h = 0.5;     %curvature on fixed type productivity in home production
 nu = 0.5;         % parameter on search 
 nu_h = 0.5;         % curvature on hours in home production
 ym_pi = 0.8 ;     % prob. husband income
 xi = 0.8;         % curvature of experience in wage
 gam_e = 1.5;      %weight on experience in wage
 tau_wf =0.6;       %wage gap
 phi_c = 0.6;   %scale in utility fun
 nI=16; nT=4; nY=3; nE=30; nZ=2;
 bisectTol = 0.00001;
 VFtol = 0.03;
 maxViter=100;
 maxHSiter=100;
 
 %Compstat Parameters
 alpha_e = alpha_e*rtoexpscale;
 tau_wf = 1-wagegapscale*(1-tau_wf);
 
 %Recessions- z=1 expansion; z=2 recession
 piz= [0.9 , 0.1; 0.3 , 0.7];
 %Husband and wife job loss rate
 lossW = [0.03 ; 0.05 ];
 lamH(1,:,:,:) = [0.96, 0.0, 0.04 ; 0.3, 0.60, 0.1 ; 0.0, 0.95, 0.05]; %(iz,EE,ES,EU; 
 lamH(2,:,:,:) = [0.93, 0.0, 0.07; 0.25, 0.60, 0.15 ; 0.0, 0.88, 0.12];
 %Husband and wife job find rate
 findW = [0.92 ; 0.85 ].*1;
 %Husband Wage by age
 wageH = [0.8; 0.8; 0.6];
 BCwage = [1;0.85]; %wage penalty during recession?
 ym = [1;0.75;0]; %dude is either employed, recently unemployed or unemployed
 %Aging Probabilities- 18-24; 25-39; 40-54; 55-64; 65-80
    %Geometric duration: 1/p;
    %piT = [1/(9^12); 1/(14^12); 1/(14^12); 1/(9^12); 1/(14^12)];
    %For the simple model: 25-39, 41-55, 56-64
    tsize = 4; %quarterly
    piT = [1-1/(14*tsize), 1/(14*tsize), 0, 0 ; 0, 1-1/(14*tsize), 1/(14*tsize), 0; 0,0 1-1/(8*tsize) , 1/(8*tsize); ];
    tgrid=(1:1:nT-1);
%Individual types: 2 wage types, 2 fixed kappa types, 4 kappa draws when young 
    %(w_l, kapy_l, kap_1 ; w_l, kapy_l, kap_2, ... kap_4 ; w_l , kapy_h,
    %kap_1,... w_h, kapy_l kap_1...)
    igrid=[1:1:nI];
 %Wages   
    ftW=ones(1,nI);
    ftW(1,1:floor(nI/2)) = 0.8;
 %Disutility    
    kapbar(1) = -0.03;  %Fixed types
    kapbar(2) = 0.00;
    kap= ones(nT,nI).*(0.2); %Base fixed cost =0.2
    kap(2,:) = kap(2,:)*0.5; %Half the cost when middle
    kap(3,:) = kap(3,:)*0.5; %Half the cost when old
    kayscale=[3,2,1,0.9];
    nKy = 4; % # fixed kappa types when young
    for iw=1:2
    for j=1:4
        %four kappa types when young
        kap(:,(iw-1)*8+(j-1)+1) = kapscale.*kap(:,(iw-1)*8+(j-1)+1)./kayscale(j); %for now the types scale the cost linearly
        kap(:,(iw-1)*8+(j-1)+5) = kapscale.*kap(:,(iw-1)*8+(j-1)+5)./kayscale(j); %for now the types scale the cost linearly        
    end
        %two fixed kappa types:
        kap(:,(iw-1)*8+1:(iw-1)*8+5) = kap(:,(iw-1)*8+1:(iw-1)*8+5)+kapbar(iw); %       
    end
 
 
% ----- Generate a non-linear grid for experience 
        egrid = linspace(0, 1, nE);
        egrid= 2*(egrid.^1.5);

%---Check looks of comparative adv in home production--------------------------------  
%{
 for ie=1:nE
   hp(ie,1)=hprod(1) ;
   w(ie,1) = wage(1,egrid(ie));
    hp(ie,2)=hprod(nI) ;
   w(ie,2) = wage(nI,egrid(ie));
end

plot(egrid,w,'-',egrid,hp,'-.')
clear hp w
%}

INfile = fullfile(SOLUTIONdir,'paras.mat')  ;
save(INfile)

%-------------------------------------

% ----- Value functions
VE = zeros(nI,nT,nE,nY,nZ); VE0 = VE; VE1 = VE; V=VE; V0=V;
VU = VE; VU0 =VE; VU1 = VE;
%VR = VE;

%-- Policy Functions
gH = ones(nI,nT,nE,nY,nZ); gS = ones(nI,nT,nE,nY,nZ); gQ= zeros(nI,nT,nE,nY,nZ);

%numerical derivative pieces- first entry is aging vs staying the same
dVedh = ones(2,nZ,nY); dVUedh = ones(2,nZ,nY); 

% ----- Value at Retirement
    %Solve w/ method of undet. coefficients
    %{
    B = (1/(piT(nT)*beta))*(((piT(nT)*beta)^(-1/crra)) * ((1+r)^((crra-1)/crra))-1)^(-crra);
    for ia=1:nA
        VR(:,ia,:,:) =  B*(agrid(ia)^(1-crra))/(1-crra);
        VE(:,ia,:,:,nT) =  B*(agrid(ia)^(1-crra))/(1-crra);
        VE0(:,ia,:,:,nT) =  B*(agrid(ia)^(1-crra))/(1-crra);
    end
    %}
    
% ----- Value Function Iteration
        % Initial Guess
        for i=1:nI
            for it=1:nT-1
           for ie = 1:nE
                for iy = 1:nY
                    for iz = 1:nZ
                       VE0(i,it,ie,iy,iz) = utilC(wageH(it)*BCwage(iz)*ym(iy)+wage(ftW(i),egrid(ie))*0.5+hprod(ftW(i))*((1-0.5)^nu_h))+utilL(0.5)-kap(it,i); 
                        VU0(i,it,ie,iy,iz) = utilC(hprod(ftW(i))*0.2+ym(iy)*wageH(nT-1)*BCwage(iz)*lamH(iz,max(iy-1,1),1));
                       V(i,it,ie,iy,iz) = max(VE0(i,it,ie,iy,iz),VU0(i,it,ie,iy,iz));
                    end 
                end 
           end
           end 
        end
%plot(egrid,squeeze(VE0(1,2,:,1,1)),egrid,squeeze(VU0(1,2,:,1,1)))
%clear VE0          
%load('VEguess.mat','VE0','VE')
%clear VU0
%load('VUguess.mat','VU0')
% ----- Calculate the Value Function for t=T-1
   for it=nT-1:-1:1
 %try just one for now
 for i=1:nI
  
     maxViterAll = 3;
 for tVall = 1:maxViterAll
     
tVE=1;
while tVE<maxViter
    %Calculate policy functions: gA_Emp and gH for Employed
            for iy = 1:nY
                for iz = 1:nZ
                    for ie = 1:nE                        
                        %Bisection on h
                        hlow = 0;
                        hhigh= 1;
                        ith=1;
                        while ith<maxHSiter
                            h=0.5*(hlow+hhigh);
                        YY = wageH(it)*BCwage(iz)*ym(iy)+wage(ftW(i),egrid(ie))*h+hprod(ftW(i))*((1-h)^nu_h);
                        YYplot(ie,iy,iz) = YY;
                        dUCdh = du_dC(YY)*(wage(ftW(i),egrid(ie))-hprod(ftW(i))*nu_h*(1-h)^(nu_h-1));
                        dULdh = du_dL(h);
                        %Find gridpont for e' 
                        ee = exp2(egrid(ie),h); 
                        if egrid(nE)< ee
                            iee = nE-1;
                            %disp('experience grid is really not large enough');
                        elseif ee<egrid(1)
                            iee=1;
                        else
                            iee= min(nE-1,find(egrid>ee,1,'first'));
                        end
 
                        %Employed
                        for izz=1:nZ
                            for iyy=1:nY
                            %Numerical Derivative wrt e'
                            dVedh(1,izz,iyy) = (V(i,it,iee+1,iy,iz)-V(i,it,iee,iy,iz))/(egrid(iee+1)-egrid(iee)); %age
                            dVedh(1,izz,iyy) = dVedh(1,izz,iyy)*de2_dh(egrid(ie),h);
                            if it+1>nT-1
                              dVedh(2,izz,iyy)=0;  
                            else
                            dVedh(2,izz,iyy) = (V(i,it+1,iee+1,iy,iz)-V(i,it+1,iee,iy,iz))/(egrid(iee+1)-egrid(iee)); %age+1
                            dVedh(2,izz,iyy) = dVedh(2,izz,iyy)*de2_dh(egrid(ie),h);
                            end
                        %Repeat for Unemployment
                            %Numerical Derivative wrt e'
                            dVUedh(1,izz,iyy) = (VU0(i,it,iee+1,iy,iz)-VU0(i,it,iee,iy,iz))/(egrid(iee+1)-egrid(iee));
                            dVUedh(1,izz,iyy) = dVUedh(1,izz,iyy)*alpha_h*psi*(h^(psi-1))*egrid(ie);               
                            if it+1>nT-1
                              dVUedh(2,izz,iyy)=0;  
                            else
                            dVUedh(2,izz,iyy) = (VU0(i,it+1,iee+1,iy,iz)-VU0(i,it+1,iee,iy,iz))/(egrid(iee+1)-egrid(iee)); %age+1
                            dVUedh(2,izz,iyy) = dVUedh(2,izz,iyy)*de2_dh(egrid(ie),h);
                            end
                            end
                        end
 
                        %Total derivative:
                            TD = 0;
                            for izz=1:nZ
                                for iyy=1:nY
                                TD =TD+beta*piz(iz,izz)*lamH(izz,iy,iyy)*(1-piT(it))*((1-lossW(izz))*(dVedh(1,izz,iyy))...
                                    +lossW(izz)*(dVUedh(1,izz,iyy)))...
                                    +beta*piz(iz,izz)*lamH(izz,iy,iyy)*(piT(it))*((1-lossW(izz))*(dVedh(2,izz,iyy))...
                                    +lossW(izz)*(dVUedh(2,izz,iyy)));
                                end
                            end
                             TD = TD+dUCdh+dULdh;
                        if (abs(TD)<bisectTol || ith>maxHSiter-2 || h<0.001 )
                            gH(i,it,ie,iy,iz)=h;
                            V1 = 0;
                            for izz=1:nZ
                            for iyy=1:nY 
                                EVe = spline(egrid,squeeze(V(i,it,:,iyy,izz)),egrid(iee));
                                EVu = spline(egrid,squeeze(VU0(i,it,:,iyy,izz)),egrid(iee));
                                EVe2 = spline(egrid,squeeze(V(i,it+1,:,iyy,izz)),egrid(iee));
                                EVu2 = spline(egrid,squeeze(VU0(i,it+1,:,iyy,izz)),egrid(iee));                                
                                V1= V1+beta*piz(iz,izz)*lamH(izz,iy,iyy)*(1-piT(it))*...
                                    ((1-lossW(izz))*EVe + lossW(izz)*EVu)...
                                   +beta*piz(iz,izz)*lamH(izz,iy,iyy)*(piT(it))*...
                                   ((1-lossW(izz))*EVe2+lossW(izz)*EVu2);
                            end
                            end
                             ctdUtil(ie,iy,iz) = V1;
                            V1 = V1+ utilC(YY)+utilL(h)-kap(it,i);
                            flowUtil(ie,iy,iz) = utilC(YY)+utilL(h)-kap(it,i);
                            VE1(i,it,ie,iy,iz) = V1;
                            ith = 99999999;                               
                        elseif TD>0 %if TD>0 then there is more to gain by working harder, so move hlow up
                            hlow = 0.8*h+0.2*hlow;
                            ith = ith+1;
                        else
                            hhigh= 0.8*h+0.2*hhigh;
                            ith = ith+1;                            
                        end %if to check solution
                            
                        end %ends loop on hours
                        
                    end %ie
                        
                 end %iz
              end %iy


 V_err=abs(max(max(max(VE1(i,it,:,:,:)- VE0(i,it,:,:,:)))))
 V_err=V_err/(1+abs(max(max(max(VE0(i,it,:,:,:))))))
 plot(egrid,squeeze(gH(i,it,:,1,1)),'-',egrid,squeeze(gH(i,it,:,2,1)),'-.',egrid,squeeze(gH(i,it,:,1,2)),'--')
    VE0=VE1;
    tVE=tVE+1
    VE=VE1;
    tVE=tVE+1;
if V_err<VFtol

    tVE
    tVE=maxViter+1;
end

%save('VEguess.mat','VE0' ,'VE')
end %tVE

save('VEguess.mat','VE')

%load('VEguess.mat','VE0' ,'VE')
%---------------------UNEMPLOYED------------------------------
%clear VU0
%load('VUguess.mat','VU0')

 if tVall<maxViterAll   
tVU=1;
while tVU<maxViter
     %Calculate policy functions: gA_Un and gS for Unemployed
           for iy = 1:nY
                for iz = 1:nZ
                    for ie = 1:nE   
                        %Find gridpont for e' ; independent of choices.
                        ee = exp2(egrid(ie),h); 
                        if egrid(nE)< ee
                            iee = nE-1;
                            %disp('experience grid is really not large enough');
                        elseif ee<egrid(1)
                            iee=1;
                        else
                            iee= min(nE-1,find(egrid>ee,1,'first'));
                        end       
                        %Bisection on h
                        hlow = 0;
                        hhigh= 1;
                        ith=1;
                        while ith<maxHSiter
                            h=0.5*(hlow+hhigh);
                        YY = wageH(it)*BCwage(iz)*ym(iy)+hprod(ftW(i))*((1-h)^nu_h);
                        dUCdh = -du_dC(YY)*(hprod(ftW(i))*nu_h*(1-h)^(nu_h-1));
 
                        %Total derivative:
                            TD = 0;
                            for izz=1:nZ
                                for iyy=1:nY
                                TD =TD+beta*piz(iz,izz)*lamH(izz,iy,iyy)*((1-piT(it))*...
                                    (findW(izz)*nu*(h)^(nu-1))*max(V(i,it,iee,iyy,izz)-VU0(i,it,iee,iyy,izz),0)+...
                                    piT(it)*(findW(izz)*nu*(h)^(nu-1))*max(V(i,it+1,iee,iyy,izz)-VU0(i,it+1,iee,iyy,izz),0));
                                end
                            end
                             TD = TD+dUCdh;
                        if (abs(TD)<bisectTol || h<0.0001 || ith>maxHSiter-2)
                            gS(i,it,ie,iy,iz)=h;
                           %Interpolated values                             
                           V1 = 0;
                            for izz=1:nZ
                            for iyy=1:nY       
                                EVe = spline(egrid,squeeze(V(i,it,:,iyy,izz)),egrid(iee));
                                EVu = spline(egrid,squeeze(VU0(i,it,:,iyy,izz)),egrid(iee));
                                EVe2 = spline(egrid,squeeze(V(i,it+1,:,iyy,izz)),egrid(iee));
                                EVu2 = spline(egrid,squeeze(VU0(i,it+1,:,iyy,izz)),egrid(iee));                                   
                                V1= V1+beta*piz(iz,izz)*lamH(izz,iy,iyy)*(1-piT(it))*...
                                    ((1-findW(izz)*(h^(nu)))*EVu+(findW(izz)*(h^(nu)))*EVe)...
                                     +beta*piz(iz,izz)*lamH(izz,iy,iyy)*(piT(it))*...
                                     ((1-findW(izz)*(h^(nu)))*EVu2+(findW(izz)*(h^(nu)))*EVe2);
                            end
                            end
                            V1 = V1+ utilC(YY)+utilL(lhome);
                            VU1(i,it,ie,iy,iz) = V1;
                            ith = 99999999;                                             
                        elseif TD>0 %if TD>0 then there is more to gain by working harder, so move hlow up
                            hlow = 0.8*h+0.2*hlow;
                            ith = ith+1;
                        else
                            hhigh= 0.8*h+0.2*hhigh;
                            ith = ith+1;                            
                        end %if to check solution
                            
                        end %ends loop on hours
                        
                    end %ie
                        
                 end %iz
              end %iy

            
 VU_err=abs(max(max(max(VU1(i,it,:,:,:)- VU0(i,it,:,:,:)))))/abs(max(max(max(VU0(i,it,:,:,:)))))
 plot(egrid,squeeze(gS(i,it,:,1,1)),'-',egrid,squeeze(gS(i,it,:,2,1)),'-.',egrid,squeeze(gS(i,it,:,1,2)),'--')

    VU0=VU1;
    VU=VU1; 
    tVU=tVU+1
if VU_err<VFtol

    tVU
    tVU=maxViter+1;

end
end
end
 for iy = 1:nY
    for iz = 1:nZ
       for ie = 1:nE 
        V0(i,it,ie,iy,iz) = max(VU(i,it,ie,iy,iz),VE(i,it,ie,iy,iz));
        if VU(i,it,ie,iy,iz)>VE(i,it,ie,iy,iz)
           gQ(i,it,ie,iy,iz) =1;
        else
           gQ(i,it,ie,iy,iz) =0; 
        end
       end
    end
 end
 
         VallErr(tVall) = abs(max(max(max(V0(i,it,:,:,:)-V(i,it,:,:,:)))))
        V=V0;
 
end

 end
   end
   
   
save('VUguess.mat','VU')

INfile = fullfile(SOLUTIONdir,'policies.mat')  ;
save(INfile, 'gS', 'gH', 'gQ')
INfile = fullfile(SOLUTIONdir,'Vfuns.mat')  ;
save(INfile, 'VU', 'VE')

%**************************************************************************************************************************************
%---------------------------------------------------------------------------------------------------------------------------------------














plot(egrid,squeeze(gQ(1,2,:,1,1)),'-',egrid,squeeze(gQ(1,2,:,2,1)),'-.',egrid,squeeze(gQ(1,2,:,1,2)),'--')
plot(egrid,squeeze(gS(1,2,:,1,1)),'-',egrid,squeeze(gS(1,2,:,2,1)),'-.',egrid,squeeze(gS(1,2,:,1,2)),'--')
plot(egrid,squeeze(VU(1,1,:,1,1)),'-',egrid,squeeze(VE(1,1,:,1,1)),'-.',egrid,squeeze(V(1,1,:,1,1)),'--')   