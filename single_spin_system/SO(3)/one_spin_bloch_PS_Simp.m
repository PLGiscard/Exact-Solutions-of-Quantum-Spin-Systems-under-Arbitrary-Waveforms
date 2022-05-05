function [rho,TimeTot,elapsedTime] = one_spin_bloch_PS_Simp(para,TMAX,NPrec,NInt,rho0)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% para: parameters as generated by paragen
% TMAX : maximum evolution time
% NPrec : number of numerical evalutation points per path-sum interval
% NInt: number of time subintervals
% rho0: initial state

n = para.n;
deltaf = para.deltaf/1000;
DeltaF = para.DeltaF/1000;
Phi0 = para.Phi0;
taup = para.taup*1000;
Omega = para.Omega/1000;
omega1 = para.omega1/1000;
deltat = para.deltat*1000;

TInt = linspace(0,TMAX,NInt+1);
TimeTot = zeros(1,NInt*NPrec);

U0 = eye(3);
rho_0 = rho0;
rho = zeros(3,1,NInt*NPrec);
rho(:,:,1) = rho_0;

P = [1./sqrt(2),1i./sqrt(2),0;
    0,0,1;
    1./sqrt(2),-1i./sqrt(2),0];

TS = SimpsonW(NPrec+1);

tic;

for i=1:NInt
    T0 = TInt(i);
    T1 = TInt(i+1);
    
    [U11, U12, U13, U21, U22, U23, U31, U32, U33,Time] = BlochEvo(T0,T1,NPrec,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n,TS);

    for j=2:NPrec+1
        U = [U11(j), U12(j), U13(j); U21(j), U22(j), U23(j); U31(j), U32(j), U33(j)]*U0;
        rho(:,1,j+(i-1)*NPrec) = P'*U*P*rho_0;
    end
    
    U0 = U;
    TimeTot(2+(i-1)*NPrec:i*NPrec+1) = Time(2:NPrec+1);
end


elapsedTime = toc;


function  [U11, U12, U13, U21, U22, U23, U31, U32, U33, Time] = BlochEvo(T0,T1,NT,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n,TS)
%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = linspace(T0,T1,NT+1); % Time points
dt = Time(2)-Time(1); % Define the true time step size
Id = eye(NT+1);

beta = (1/2).*omega1.*exp(-2^(n+2).*((Time-deltat)/taup).^n+1i.*(Phi0+pi.*DeltaF.*(Time-deltat).^2./taup-2.*pi.*deltaf.*(Time-deltat)));

B = triu(meshgrid(beta)); BC = conj(B);
betaC = BC(1,:);

eOmega = exp(-1i.*Omega.*Time);
eOmegaC = conj(eOmega);

fun = eOmega.*beta; funC = conj(fun);
funmesh = meshgrid(fun);

BFourier = sum(dt*9/24.*TS.*funmesh.',1);
[a,b] = meshgrid(BFourier);
BFourier = a-b;
BFourierC = conj(BFourier);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U11 %%%%%%%%%%%%%%%%%%%%%%%%
% (1-Ker232)^*-1 yields the contribution of all walks from 2 to itself on
% the graph G\1. 
Ker232 = (-2).*funmesh.*BFourierC;
G232 = inv(Id - 9/24*dt* TS.*Ker232) - Id;

% Insertion of walks from 2 to itself in G\1 between edges 1->2 and 2->1
% and addition of self loop 1->1
Ker11 = Id - dt*9/24.*TS.*(1i.*Omega + (-2)*dt*StarProd(BC,StarProd(G232,B,TS),TS) + (-2).*dt.*StarProd(BC,B,TS));

v = zeros(NT+1,1);v(1)=1;
opts.UT = true;
opts.TRANSA = true;
G11 = linsolve(Ker11,v,opts)';
U11 = sum(TS.*meshgrid(G11).',1);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U33 %%%%%%%%%%%%%%%%%%%%%%%%
U33 = conj(U11); % By symmetry of the Hamiltonian

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U22 %%%%%%%%%%%%%%%%%%%%%%%%
% Walks from 2 to itself are walks from 2 to 1 and back or from 2 to 3 and
% back. Each contribution is conjugate to the other so the sum of the two
% yields twice the real part of one. Here we take the contribution of walks
% 2->3->2 from before, encapsulated in Ker232
Ker22 = Id-9/24*dt*TS.*(2*real(Ker232));

G22 = linsolve(Ker22,v,opts)'; % G22 is the *-resolvent of Ker22
U22 = sum(TS.*meshgrid(G22).',1);
U22 = U22(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U12 %%%%%%%%%%%%%%%%%%%%%%%%
% Walks 2->1 are walks from 2 to itself then 2->1 then 1 to itself on G\2.
U12 = -1i*sqrt(2)*(9/24*dt)*eOmegaC.*sum(TS.*meshgrid(fun.*U22).',1);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U32 %%%%%%%%%%%%%%%%%%%%%%%%
U32 = conj(U12); % By symmetry of the Hamiltonian

%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U21 %%%%%%%%%%%%%%%%%%%%%%%%
G21 = -1i*sqrt(2)*9/24*(betaC.*U11)*(TS.*G232);
G21 = G21(1,:)-1i.*sqrt(2).*betaC.*U11;

U21 = (9/24*dt)*sum(TS.*meshgrid(G21).',1);

%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U23 %%%%%%%%%%%%%%%%%%%%%%%%
U23 = conj(U21); % By symmetry of the Hamiltonian

%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U31 %%%%%%%%%%%%%%%%%%%%%%%%
U31 = 1i*sqrt(2)*(9/24*dt)*eOmega.*sum(TS.*meshgrid(funC.*U21).',1);

%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U13 %%%%%%%%%%%%%%%%%%%%%%%%
U13 = conj(U31); % By symmetry of the Hamiltonian

function TS = SimpsonW(NT)
% Creates the weight matrix for Simpson quadrature
TS = eye(NT) + diag(28/9*ones(NT-1,1),1) + diag(23/9*ones(NT-2,1),2);
TS(1,1:NT) = 1;

TS(2,1:NT)=28/9;
TS(3,1:NT)=23/9;

for i=7:NT
    TS(4:i-3,i)=24/9*ones(i-6,1);
end

TS(1:5,1:5) = [1,0,0,0,0; 1,1,0,0,0; 1,4,1,0,0; 1,28/9,23/9,1,0; 1,28/9,23/9,28/9,1].';

 function C=StarProd(A,B,TS)
 C = 9/24*(A*(TS.*B));