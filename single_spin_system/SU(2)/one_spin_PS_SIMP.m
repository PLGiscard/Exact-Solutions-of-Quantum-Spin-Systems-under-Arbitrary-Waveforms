function [rho,TimeTot,elapsedTime]=one_spin_PS_SIMP(para,TMAX,NPrec,NInt)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% para : parameters as set by paragen
% TMAX : in second, maximum evolution time
% NPrec : number of numerical evalutation points per path-sum interval
% NInt: number of time subintervals

%%%%%%%%%%%%%%%%%%%%%%% Setting beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for super-Gaussian pulse shape

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
U0 = eye(2);

rho_0 = [1,0;0,-1];
rho = zeros(2,2,NInt*NPrec);
rho(:,:,1) = rho_0;

tic;

TS = SimpsonW(NPrec+1); % Weight matrix for implementing Simpson quadrature rule

for i=1:NInt
    T0 = TInt(i);
    T1 = TInt(i+1);

    [U11,U12,U21,U22,Time] = UPSonespin(T0,T1,NPrec,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n,TS);
    
    for j=2:NPrec+1
        U = [U11(j), U12(j); U21(j), U22(j)]*U0;
        rho(:,:,j+(i-1)*NPrec) = U*rho_0*U';
    end
    
    U0 = U;
    TimeTot(2+(i-1)*NPrec:i*NPrec+1) = Time(2:NPrec+1);
    
end

elapsedTime = toc;

rho11 = reshape(rho(1,1,:),1,[]);
rho12 = reshape(rho(1,2,:),1,[]);
rho21 = reshape(rho(2,1,:),1,[]);
rho22 = reshape(rho(2,2,:),1,[]);

subplot(2,2,1)
plot(TimeTot,real(rho11),'-r');
ylim([-1 1])

subplot(2,2,2)
plot(TimeTot,real(rho12),'-r');
ylim([-1 1])

subplot(2,2,3)
plot(TimeTot,real(rho21),'-r');
ylim([-1 1])

subplot(2,2,4)
plot(TimeTot,real(rho22),'-r');
ylim([-1 1])

end
    
function [U11,U12,U21,U22,Time] = UPSonespin(T0,T1,NT,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n,TS)
%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = linspace(T0,T1,NT+1); % Time points
dt = Time(2)-Time(1); % Define the true time step size
Id = eye(NT+1);

beta = (1/2)*omega1*exp(-2^(n+2).*((Time-deltat)/taup).^n+1i.*(Phi0+pi.*DeltaF.*(Time-deltat).^2./taup-2.*pi.*deltaf.*(Time-deltat)));
eOmega = exp(-1i*Omega/2*Time);

fun = eOmega.*beta;
funmesh = meshgrid(fun);

BFourier = sum(dt*9/24.*TS.*funmesh.',1); % Effectively performs the FT of the windowed beta

[a,b] = meshgrid(BFourier);
BFourier = a-b;

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U11 %%%%%%%%%%%%%%%%%%%%%%%%
M = -1i*Omega/2 + (-1i)^2*conj(funmesh).*BFourier;
Ker11 = Id - dt*9/24.*TS.*M; % Simpson quadrature for the resolvent

v = zeros(NT+1,1);v(1)=1;
opts.UT = true; % Enforces upper triangularity 
opts.TRANSA = true;
G11 = linsolve(Ker11,v,opts)';
U11 =  sum(TS.*meshgrid(G11).',1);

%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U21 %%%%%%%%%%%%%%%%%%%%%%%%
U21 = -1i*dt*9/24*conj(eOmega).*sum(TS.*meshgrid(fun.*U11).',1);

%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U22 & U12 %%%%%%%%%%%%%%%%%%%%%%%%
U22 = conj(U11); % By symmetry of the Hamiltonian
U12 = -conj(U21);% By symmetry of the Hamiltonian

end

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
end
