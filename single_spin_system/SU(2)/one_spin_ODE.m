function [rho,elapsedTime] = one_spin_ODE(para,TMAX,NT,tol)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% para : parameters as set by paragen
% TMAX : in second, maximum evolution time
% NT : number of numerical evalutation points 
% tol : relative accuracy required

options = odeset('RelTol',tol);

tic;

g0 = [1;0;0;-1];

[t, Sol]=ode45(@(t,g) myode(t,para,g),linspace(0,TMAX,NT),g0,options);

rho(1,1,:) = Sol(:,1);
rho(2,1,:) = Sol(:,2);
rho(1,2,:) = Sol(:,3);
rho(2,2,:) = Sol(:,4);

elapsedTime = toc;

Time = t;

subplot(2,2,1)
plot(Time*1000,real(Sol(:,1)));
ylim([-1 1])

subplot(2,2,2)
plot(Time*1000,real(Sol(:,2)));
ylim([-1 1])

subplot(2,2,3)
plot(Time*1000,real(Sol(:,3)));
ylim([-1 1])

subplot(2,2,4)
plot(Time*1000,real(Sol(:,4)));
ylim([-1 1])

end

function dg = myode(t,para,g)

[Cx,Cy] = chirp_fun(t,para);

H = [para.Omega/2,   Cx-1i*Cy;...
     Cx+1i*Cy,      -para.Omega/2];

H_hat = kron(eye(2),H)-kron(conj(H),eye(2));

gvec = [g(1); g(2); g(3); g(4)];

dg = -1i*H_hat*gvec;

end

function [Cx,Cy]=chirp_fun(t,para)

smfactor = para.n;
offs_f = para.deltaf;
bandwidth = para.DeltaF;
phi0 = para.Phi0;
tau_p = para.taup;
omega1 = para.omega1;
offs_t = para.deltat;

Cx = (1/2)*(exp(-(2^(smfactor+2))*((t-offs_t)/tau_p).^smfactor)).*(omega1*cos(phi0+(pi*bandwidth*((t-offs_t).^2)/tau_p)-2*pi*offs_f*(t-offs_t)));
Cy = (1/2)*(exp(-(2^(smfactor+2))*((t-offs_t)/tau_p).^smfactor)).*(omega1*sin(phi0+(pi*bandwidth*((t-offs_t).^2)/tau_p)-2*pi*offs_f*(t-offs_t)));

end
