function [rho,TimeTot,elapsedTime]=one_spin_PS_TRAP(para,TMAX,NPrec,NInt)
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

for i=1:NInt
    T0 = TInt(i);
    T1 = TInt(i+1);

    [U11,U12,U21,U22,Time] = UPSonespin(T0,T1,NPrec,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n);
    
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
    
function [U11,U12,U21,U22,Time] = UPSonespin(T0,T1,NT,Omega,taup,DeltaF,Phi0,omega1,deltat,deltaf,n)
%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = linspace(T0,T1,NT+1); % Time points
dt = (Time(2)-Time(1)); % Define the true time step size
Id = eye(NT+1);

beta = (1/2).*omega1.*exp(-2^(n+2).*((Time-deltat)/taup).^n+1i.*(Phi0+pi.*DeltaF.*(Time-deltat).^2./taup-2.*pi.*deltaf.*(Time-deltat)));
eOmega = exp(-1i.*Omega/2.*Time);

fun = eOmega.*beta;
BFourier = cumtrapz(dt*fun); % Effectively performs the FT of the windowed beta

% BFourier is the matrix 
% of B_{t,t'}(Omega/2pi) indexed by t' (columns)
% and t (lines). Here not upper triangular, this is enforced later.
[a,b] = meshgrid(BFourier);
BFourier = a-b;

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U11 %%%%%%%%%%%%%%%%%%%%%%%%
Ker121 = (-1i)^2*meshgrid(conj(fun)).*BFourier;
M = dt*(-1i*Omega/2 + Ker121);
dM = -1i*dt*Omega/2 * Id; % Diagonal matrix with the diagonal of M
Ker11 = triu(Id - M + dM/2); % Trapezoidal quadrature for the resolvent

v = zeros(NT+1,1);v(1)=1;
opts.UT = true; % Enforces upper triangularity 
opts.TRANSA = true;
G11 = linsolve(Ker11,v,opts)'; % G11 is the *-resolvent of Ker11 in (t,0)

G11(1) = 0; 
U11 =  1 + cumtrapz(G11);

%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U21 %%%%%%%%%%%%%%%%%%%%%%%%
U21 = -1i.*conj(eOmega).*dt.*cumtrapz(fun.*U11);

%%%%%%%%%%%%%%%%%%%%%%%% Calculation of U22 & U12 %%%%%%%%%%%%%%%%%%%%%%%%
U22 = conj(U11); % By symmetry of the Hamiltonian
U12 = -conj(U21);% By symmetry of the Hamiltonian