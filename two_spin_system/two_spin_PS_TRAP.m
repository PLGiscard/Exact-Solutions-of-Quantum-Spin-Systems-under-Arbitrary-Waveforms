function [rho,TimeTot,elapsedtime]=two_spin_PS_TRAP(para,TMAX,NPrec,NInt)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% TMAX : en ms, maximum evolution time
% NPrec : number of numerical evalutation points per subinterval
% NInt: number of subintervals
% Path-sum trapezoidal 

%%%%%%%%%%%%%%%%%%%%%%% Setting beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for super-Gaussian pulse shape
n = para.n;
deltaf = para.deltaf/1000;
DeltaF = para.DeltaF/1000;
Phi0 = para.Phi0;
taup = para.taup*1000;
Omega1 = para.Omega(1)/1000;
Omega2 = para.Omega(2)/1000;
J = para.J/1000;
omega1 = para.omega1/1000;
deltat = para.deltat*1000;

%%%%%%%%%%%%%%%%%%%%%%% Building up the subintervals %%%%%%%%%%%%%%%%%%%%%%
TInt = linspace(0,TMAX,NInt+1);
TimeTot = zeros(1,NInt*NPrec);
U0 = eye(4);

Sigmaz = 0.5*[1,0;0,-1];
Lz = kron(Sigmaz,eye(2));
Sz = kron(eye(2),Sigmaz);
rho_0 = (Lz+Sz)/norm(Lz+Sz); % initial state
rho = zeros(4,4,NInt*NPrec);
rho(:,:,1) = rho_0;

tic;

for i=1:NInt
    T0 = TInt(i);
    T1 = TInt(i+1);
    
    [U11, U12, U13, U14, U21, U22, U23, U24, U31, U32, U33, U34, U41, ...
        U42, U43, U44, Time] = UPStwospins(T0,T1,NPrec,Omega1,Omega2,J,...
                            taup,DeltaF,Phi0,omega1,deltat,deltaf,n);

    for j=2:NPrec+1
        U = [U11(j), U12(j), U13(j), U14(j); U21(j), U22(j), U23(j), U24(j);...
             U31(j), U32(j), U33(j), U34(j); U41(j), U42(j), U43(j), U44(j)]*U0;
        rho(:,:,j+(i-1)*NPrec) = U*rho_0*U';
    end
    
    U0 = U;
    TimeTot(2+(i-1)*NPrec:i*NPrec+1) = Time(2:NPrec+1);

end

elapsedtime = toc;

end


function [U11, U12, U13, U14, U21, U22, U23, U24, U31, U32, U33, U34, U41, U42, U43, U44,Time] = UPStwospins(T0,T1,NT,Omega1,Omega2,J,taup,DeltaF,Phi0,omega1,deltat,deltaf,n)
% Computes the evolution operator over a time interval
NT=NT+1;
%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = linspace(T0,T1,NT); % Time points
dt = Time(2)-Time(1); % Define the true time step size

h11 = (1/2)*(pi*J+Omega1+Omega2) - (1/2)*(pi*J); 
h22 = (1/2)*(-pi*J+Omega1-Omega2)- (1/2)*(pi*J);
h33 = (1/2)*(-pi*J-Omega1+Omega2)- (1/2)*(pi*J); 
h44 = (1/2)*(pi*J-Omega1-Omega2) - (1/2)*(pi*J); 

beta = (1/2).*omega1.*exp(-2^(n+2).*((Time-deltat)/taup).^n ...
    + 1i.*(Phi0 + pi.*DeltaF.*(Time-deltat).^2./taup-2.*pi.*deltaf.*(Time-deltat)));
betaC = conj(beta);

H = triu(ones(NT,NT))-(1/2)*eye(NT);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of UII, that is U22, U33, U23 and U32 %%%%%%%%%%%%%%%%%%%%%%%%
eOmega12 = exp(-1i.*(Omega1+Omega2)/2.*Time);
fun = eOmega12.*beta;

BFourier(1,:) = cumtrapz(dt*fun);
[a,b] = meshgrid(BFourier(1,:));
BFourier = a-b;
BFourierC = conj(BFourier);

KerII = triu(-2*real(meshgrid(fun).*BFourierC));
dKerII = diag(diag(KerII));
KdK = KerII-dKerII./2;
PII = [KdK, KdK; KdK, KdK];
HII = -1i*[ h22.*H, (pi.*J).*H;  (pi.*J).*H, h33.*H];

KII = eye(2*NT) - dt*HII - dt*PII ;
opts.TRANSA=true;

v = zeros(2*NT,1);v(1)=1;
GII = linsolve(KII,v,opts);
g22 = GII(1:NT)'; 

g22(1) = g22(1)-1; 
g23 = GII(NT+1:2*NT)';

v(1)=0;v(NT+1)=1;
GII = linsolve(KII,v,opts);
g32 = GII(1:NT)';
g33 = GII(NT+1:2*NT)'; g33(1) = g33(1)-1;

U22 = cumtrapz(Time,g22./dt)+1;
U33 = cumtrapz(Time,g33./dt)+1;
U23 = cumtrapz(Time,g23./dt);
U32 = cumtrapz(Time,g32./dt);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U1II and U4II, that is U12, U13, U42 and U43 %%%%%%%%%%%%%%%%%%%%%%%%
U12 = -1i.*exp(-1i.*h11.*Time).*cumtrapz(Time,exp(-1i.*h11.*(-Time)).*betaC.*(U22+U32));
U13 = -1i.*exp(-1i.*h11.*Time).*cumtrapz(Time,exp(-1i.*h11.*(-Time)).*betaC.*(U23+U33));
U42 = -1i.*exp(-1i.*h44.*Time).*cumtrapz(Time,exp(-1i.*h44.*(-Time)).*beta.*(U22+U32));
U43 = -1i.*exp(-1i.*h44.*Time).*cumtrapz(Time,exp(-1i.*h44.*(-Time)).*beta.*(U23+U33));

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U11, U44, U14 and U41 %%%%%%%%%%%%%%%%%%%%%%%%
w = sqrt((pi*J)^2 + ((Omega1-Omega2)/2)^2);

coslist = cos(Time*w);
sinlist = sin(Time*w);
explistminus = exp(-1i*pi*J*Time);
explistplus = 1./explistminus;
explist2plus = explistplus.*explistplus;
ecos = explistplus.*coslist; eCcos = conj(ecos);
esin = explistplus.*sinlist; eCsin = conj(esin);

bcos = beta.*eCcos; 
mbcos = meshgrid(explist2plus.*bcos);
bCcos = betaC.*eCcos; 
mbCcos = meshgrid(explist2plus.*bCcos);

Int_bcos = cumtrapz(Time,bcos);
[a,b] = meshgrid(Int_bcos); Int_bcos = triu(a-b);
Int_bCcos = cumtrapz(Time,bCcos);
[a,b] = meshgrid(Int_bCcos); Int_bCcos = triu(a-b);

bsin = beta.*eCsin; 
mbsin = meshgrid(explist2plus.*bsin);
bCsin = betaC.*eCsin;
mbCsin = meshgrid(explist2plus.*bCsin);

Int_bsin = cumtrapz(Time,bsin);
[a,b] = meshgrid(Int_bsin); Int_bsin = triu(a-b);
Int_bCsin = cumtrapz(Time,bCsin);
[a,b] = meshgrid(Int_bCsin); Int_bCsin = triu(a-b);

a11 = Int_bcos + (1i/w) * (pi*J + (Omega1-Omega2)/2) * Int_bsin;
a12 = Int_bCcos + (1i/w) * (pi*J + (Omega1-Omega2)/2) * Int_bCsin;
a21 = Int_bcos + (1i/w) * (pi*J - (Omega1-Omega2)/2) * Int_bsin;
a22 = Int_bCcos + (1i/w) * (pi*J - (Omega1-Omega2)/2) * Int_bCsin;

A11 = a11.*(mbCcos -(1i/w) * (pi*J + (Omega1-Omega2)/2) * mbCsin) ...
    + a21.*(mbCcos -(1i/w) * (pi*J - (Omega1-Omega2)/2) * mbCsin);
A21 = a12.*(mbCcos -(1i/w) * (pi*J + (Omega1-Omega2)/2) * mbCsin) ...
    + a22.*(mbCcos -(1i/w) * (pi*J - (Omega1-Omega2)/2) * mbCsin);
A12 = a11.*(mbcos - (1i/w) * (pi*J + (Omega1-Omega2)/2) * mbsin) ...
    + a21.*(mbcos - (1i/w) * (pi*J - (Omega1-Omega2)/2) * mbsin);
A22 = a12.*(mbcos - (1i/w) * (pi*J + (Omega1-Omega2)/2) * mbsin) ...
    + a22.*(mbcos - (1i/w) * (pi*J - (Omega1-Omega2)/2) * mbsin);

dA11 = diag(diag(A11));
dA12 = diag(diag(A12));
dA21 = diag(diag(A21));
dA22 = diag(diag(A22));

P = (-1i)^2.*[A11-dA11./2, A12-dA12./2; A21-dA21./2, A22-dA22./2];
C = -1i.*[H.*h11, zeros(NT,NT); zeros(NT,NT), H*h44];

K14 = eye(2*NT) - dt*C - dt*P;
opts.TRANSA=true;

v = zeros(2*NT,1);v(1)=1;
G1144 = linsolve(K14,v,opts);
g11 = G1144(1:NT)'; g11(1) = g11(1)-1;
g14 = G1144(NT+1:2*NT)';

v(1)=0;v(NT+1)=1;
G1144 = linsolve(K14,v,opts);
g41 = G1144(1:NT)';
g44 = G1144(NT+1:2*NT)'; g44(1) = g44(1)-1;

U11 = cumtrapz(Time,g11/dt)+1;
U44 = cumtrapz(Time,g44/dt)+1;
U41 = cumtrapz(Time,g14/dt);
U14 = cumtrapz(Time,g41/dt);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U21, U24, U31 and U34 %%%%%%%%%%%%%%%%%%%%%%%%
bu11 = -1i*beta.*U11;
bu14 = -1i*beta.*U14;
bCu41 = -1i*betaC.*U41;
bCu44 = -1i*betaC.*U44;

a11 = bu11 + bCu41;
a12 = bu14 + bCu44;

A11 = cumtrapz(Time,a11.*(eCcos + (1i/w)*(pi*J + (Omega1-Omega2)/2)*eCsin));
A12 = cumtrapz(Time,a12.*(eCcos + (1i/w)*(pi*J + (Omega1-Omega2)/2)*eCsin));
A21 = cumtrapz(Time,a11.*(eCcos + (1i/w)*(pi*J - (Omega1-Omega2)/2)*eCsin));
A22 = cumtrapz(Time,a12.*(eCcos + (1i/w)*(pi*J - (Omega1-Omega2)/2)*eCsin));

U21 = A11.* (ecos - (1i/w)*((Omega1-Omega2)/2)*esin) ...
    + A21.* (-(1i/w)*pi*J*esin);
U24 = A12.* (ecos - (1i/w)*((Omega1-Omega2)/2)*esin) ...
    + A22.* (-(1i/w)*pi*J*esin);
U31 = A21.* (ecos + (1i/w)*((Omega1-Omega2)/2)*esin) ...
    + A11.* (-(1i/w)*pi*J*esin);
U34 = A22.* (ecos + (1i/w)*((Omega1-Omega2)/2)*esin) ...
    + A12.* (-(1i/w)*pi*J*esin);

end