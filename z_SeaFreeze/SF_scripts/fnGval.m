function  Results=fnGval(sp,input)
% function to return rho,vel, G, Cp, alpha S U A H K and Kp for G splines in either (P and T) or (m, P, and T)
%  ALL MKS with P in MPa  
%       exception:  partial molal and apparent volumes in cc/mol 
%
%    Usage: Results=fnGval(sp,input) 
%   Results=Result.rho,vel,G,Cp,alpha,S,U,A,H,Kt,Kp,Ks,mus,muw
%          where input is either a npts by (2 or 3) matrix of scatter points in [P  T m] or input is a cell of {P,T,m}  
%                 m in molality, P in MPa and T in K.  MW is molecular
%                 weight in kg/mol - needed for chemical potentials and
%                 partial molar quantities
%           rho in kg/m^3, vel in m/s, G in J/kg Cp in J/kg/K alpha in K^(-1)
%            mu is chemical potential in J/mol
% 
% JMB 2015-2019

nw=1000/18.01528;  % number of moles of water in a kilogram of water
nw=1000/18.015268; % TRF value
R=8.3144;  %kJ/mol/K


%Determine whether this is P-T or P-T-m spline for gridded or scattered points 
mu_flg=1;
if length(sp.knots)==2
    mu_flg=0;
end

flgTc=0;
%determine whether dimensionless T is used
if isfield(sp,'Tc')
    flgTc=1;
    Tc=sp.Tc;
end
    

% determine whether scattered or gridded 
if iscell(input) % gridded data
   nd=length(input);  
else % scattered data in columns
   [~,nd]=size(input);
end

%handle P T representations
if nd==1
    error('Check input PTm: fnGval finds one column of input. Input needs to be column centric: [P(:) T(:) optional m(:)]')
elseif nd==2   % spline in P and T only
    if iscell(input)  % gridded output
        P=input{1};T=input{2};    
        [Pm,Tm]=ndgrid(P,T);
        if flgTc
            tau=log(T/Tc);
            %taum=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=T;
        end

        G=sp_val(sp,{P,tau});
        d1T=sp_val(sp,[ 0 1],{P,tau});
        d2T=sp_val(sp,[ 0 2],{P,tau});
        d1P=sp_val(sp,[ 1 0],{P,tau});
        dPT=sp_val(sp,[ 1 1],{P,tau});
        d2P=sp_val(sp,[ 2 0],{P,tau});
        d3P=sp_val(sp,[ 3 0],{P,tau});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        if flgTc
            tau=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=Tm;
        end
        input=[input(:,1) tau];

        G=sp_val(sp,input);
        d1T=sp_val(sp,[ 0 1],input);
        d2T=sp_val(sp,[0 2],input);
        d1P=sp_val(sp,[1 0],input);
        dPT=sp_val(sp,[1 1],input);
        d2P=sp_val(sp,[2 0],input);
        d3P=sp_val(sp,[3 0],input);
    end

elseif nd==3  % spline in P, T and compositions  Dimensionless T not implemented here 11/2020
    if (isfield(sp,'MW'))
        M=sp.MW;
 %       default.  In the  future water may not be  the solvent and the
 %       first enty gives the MW for the solvent.
        if length(M)==2
           M=M(2);
        end
    else
        error('a compositional LBF needs the molecular weight set in the structure')
    end
    if (isfield(sp,'nu'))
        nu=sp.nu;
    else
        error('a compositional LBF needs the number of ions in solution (nu) set in the structure')
    end
    if (isfield(sp,'Go'))
        Goin=sp.Go;         
    else
        error('a compositional LBF needs the standard state of the solute at the reference P in the form of a univarient LBF')
    end
    if iscell(input) % gridded output
        P=input{1};T=input{2};  m=input{3}; 
        if(isstruct(Goin))
            Go=fnval(Goin,T);  % using a spline for the standard state vs T
        else
            Go=1;
        end
        m(m==0)=eps; % add eps to zero concentrations to avoid divide by zero
        mflg=0;
        if mu_flg
        if (m(1)>eps)
          m=[eps;m(:)];  % add in a zero value in order to calculate apparent quantities and remove it later
          mflg=1;
        end
        end   

        [Pm,Tm,mm]=ndgrid(P,T,m);
        if flgTc
            tau=log(T/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=T;
        end
        G=sp_val(sp,{P,tau,m});
        d1T=sp_val(sp,[ 0 1 0],{P,tau,m});
        d1Tm=sp_val(sp,[ 0 1 1],{P,tau,m});
        d2T=sp_val(sp,[ 0 2 0],{P,tau,m});
        d2Tm=sp_val(sp,[ 0 2 1],{P,tau,m});
        d1P=sp_val(sp,[ 1 0 0],{P,tau,m});
        dPm=sp_val(sp,[ 1 0 1],{P,tau,m});
        dPT=sp_val(sp,[ 1 1 0],{P,tau,m});
        d2P=sp_val(sp,[ 2 0 0],{P,tau,m});
        d3P=sp_val(sp,[ 3 0 0],{P,tau,m});
        dGdm=sp_val(sp,[ 0 0 1],{P,tau,m});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        mm=input(:,3); 
        if flgTc
            tau=log(Tm/Tc);
            tautm=Tm.^-1;
            tau2tm=-Tm.^-2;   
        else
            tau=Tm;
        end
        if(isstruct(Goin))
            Go=sp_val(Goin,Tm);  % using a spline for the standard state vs T
        else
            Go=1;
        end
        
        m=mm;
        if mu_flg
          mm0=zeros(size(mm))+eps;
          mflg=1;
        else
            mflg=0;
        end
        input=[input(:,1) tau input(:,3)];
        G=sp_val(sp,input);
        d1T=sp_val(sp,[ 0 1 0],input);
        d1Tm=sp_val(sp,[ 0 1 1],input);
        d2T=sp_val(sp,[0 2 0],input);
        d2Tm=sp_val(sp,[0 2 1],input);
        d1P=sp_val(sp,[1 0 0],input);
        dPm=sp_val(sp,[1 0 1],input);
        dPT=sp_val(sp,[1 1 0],input);
        d2P=sp_val(sp,[2 0 0],input);
        d3P=sp_val(sp,[3 0 0],input);
        dGdm=sp_val(sp,[ 0 0 1],input);
        
%        calculate zero concentration derivatives to determine apparent Volume and specific heat
        if mu_flg
           d2T0=sp_val(sp,[0 2 0],[Pm(:) tau(:) mm0(:)]);
           d1T0=sp_val(sp,[0 1 0],[Pm(:) tau(:) mm0(:)]);
           d1P0=sp_val(sp,[1 0 0],[Pm(:) tau(:) mm0(:)]);
           if flgTc
              d2T0=d2T0.*tautm.^2+d1T0.*tau2tm;
           end
           V0=1e-6*d1P0;  % 1e6 for MPa to Pa
           Cp0=-d2T0.*Tm(:);
        end
    end
end

if flgTc
    d2T=d2T.*tautm.^2+d1T.*tau2tm;
    d2Tm=d2Tm.*tautm.^2+d1Tm.*tau2tm;
    d1T=d1T.*tautm;
    d1Tm=d1Tm.*tautm;
    dPT=dPT.*tautm;
end

Cp=-d2T.*Tm;
Cv= Cp + Tm.*dPT.^2./d2P;
S=-d1T;
rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa

vel=real(sqrt(d1P.^2./(dPT.^2./d2T - d2P))); % MPa-Pa units conversion cancels. 
Ks=rho.*vel.^2/1e6;
% Problem with poor parameters that might produce an imaginary sound
% speed is overcome by reporting just the real component.  
% rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa
% Ks=rho.*vel.^2/1e6;
alpha=1e-6*dPT.*rho; % 1e6 for MPa to Pa

U=G-1e6*Pm./rho+Tm.*S;
A=U-Tm.*S;
H=G+Tm.*S;
Kt=-d1P./d2P;
Kp=d1P.*d2P.^(-2).*d3P -1;
Js=-dPT./d2T;

if iscell(input) % gridded output
  if mu_flg==1   
     f=1+M*mm;
     x=mm./(mm+nw);
     V=rho.^-1;
     mus=M*G + f.*dGdm;
     muw=G/nw-f.*mm.*dGdm/nw;   
     G_ss=Go+mus(:,:,1)-mus(1,:,1);       
     Vm=(M*d1P +f.*dPm);
     Cpm=Cp*M - f.* d2Tm.*Tm;
     Cpa=(Cp.*f -repmat(squeeze(Cp(:,:,1)),1,1,length(m)))./mm;
     Va=1e6*(V.*f - repmat(squeeze(V(:,:,1)),1,1,length(m)))./mm;
     Vex= Va-Vm(:,:,1); 
     gam=exp(1/R*1/nu*(mus-G_ss)./Tm - log(mm));
     phi=-nw*(muw-G(:,:,1)/nw)./mm./Tm/R/nu;
     aw=exp(-mm.*phi*2/nw);
     %Gex=nu*R*Tm.*mm.*(log(gam)-(phi-1)); % factor of m?
    Gex=G.*f-G(:,:,1)+mm.*(G_ss+nu*R*Tm.*(1-log(mm))); % not debugged

     if mflg  % remove the added zero concentration ppoint
       mm=mm(:,:,2:end);
       Va=Va(:,:,2:end);
       Cpa=Cpa(:,:,2:end);
       Vm=Vm(:,:,2:end);
       Cpm=Cpm(:,:,2:end);
       G=G(:,:,2:end);
       rho=rho(:,:,2:end);
       Cp=Cp(:,:,2:end);
       Cv=Cv(:,:,2:end);
       S=S(:,:,2:end);
       vel=vel(:,:,2:end);
       Ks=Ks(:,:,2:end);
       alpha=alpha(:,:,2:end);
        Js=Js(:,:,2:end);
       U=U(:,:,2:end);
       H=H(:,:,2:end);
        A=A(:,:,2:end);
       Kt=Kt(:,:,2:end);
       Kp=Kp(:,:,2:end);
       mus=mus(:,:,2:end);
       muw=muw(:,:,2:end);
       f=f(:,:,2:end);
       gam=gam(:,:,2:end);
       phi=phi(:,:,2:end);
       G_ss=G_ss(:,:,2:end);
       Vex=Vex(:,:,2:end);
       Gex=Gex(:,:,2:end);
       aw=aw(:,:,2:end);
      
     end
  else
    mus=[];
    muw=[];
    Cpa=[];
    Va=[];
    f=[];
    Vm=[];
    Cpm=[];
    gam=[];
    phi=[];
    G_ss=[];
    Vex=[];
    Gex=[];
    aw=[];
  end
else  % Scattered data  - broken at moment - need to calculate second set of points at zero concentration
   if mu_flg==1   
%        warning('The mixing quantities (Vex, Gex, gamma (gam),etc are currently broken - needs code revision')
      f=1+M*mm;
      V=rho.^-1;
      mus=M*G + f.*dGdm;
      muw=G/nw - 1/nw*f.*mm.*dGdm;
      Vm=M*d1P +f.*dPm;
      Cpm=Cp*M - f.* d2Tm.*Tm;
      Cpa=(Cp.*f - Cp0)./mm;
      Va=1e6*(V.*f - V0)./mm;
      Vex=Va-Vm(:,:,1);
      G_ss=Go+R*nu*Tm.*log(mm);
      gam=exp(1/R*1/nu*(mus-G_ss)./Tm - log(mm));
      phi=-nw*(muw-G(:,:,1)/nw)./mm./Tm/R/nu;
      Gex=nu*R*Tm.*(log(gam)-(phi-1)); 
      aw=exp(-mm.*phi*2/nw);

   else
     mus=[];
     muw=[];
     Cpa=[];
     Va=[];
     f=[];
     Vm=[];
     Cpm=[];
     gam=[];
     phi=[];
     G_ss=[];
     Vex=[];
     Gex=[];
     aw=[];
   end
end



Results.G=G;
Results.S=S;
Results.U=U;
Results.H=H;
Results.A=A;
Results.F=A;
Results.rho=rho;
Results.Cp=Cp;
Results.Cv=Cv;
Results.Kt=Kt;
Results.Ks=Ks;
Results.Kp=Kp;
Results.alpha=alpha;
Results.vel=vel;
Results.Js=Js;
Results.gamma=1e6*Results.alpha.*Results.Kt./Results.rho./Results.Cv;
if iscell(input)
   Results.P=P;
    Results.T=T; 
else
    Results.P=Pm;
    Results.T=Tm;
end
if mu_flg
    Results.Va=Va;
    Results.Cpa=Cpa;
    Results.mus=mus;
    Results.muw=muw;
    Results.f=f;
    Results.m=mm;
    Results.Vm=Vm;
    Results.Cpm=Cpm;
    Results.gam=gam;
    Results.phi=phi;
    Results.G_ss=G_ss;
    Results.Vex=Vex;
    Results.Gex=Gex;
    Results.aw=aw;
end

