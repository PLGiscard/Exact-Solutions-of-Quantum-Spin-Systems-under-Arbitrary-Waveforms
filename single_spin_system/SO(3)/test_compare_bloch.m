%% M. Foroozandeh, P.-L. Giscard, 04/2022
% Run test_compare_bloch to compute the one-spin system [in SO(3) 
% representation] solution by all methods available.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% para = paragen(Omega,alpha,DeltaF,taup,Phi0,deltat,deltaf,n);
para =  paragen(2*pi*1000,180,100000,0.001,0,0.0005,0,30);
TMAX = 0.001; % Max simulation time, in s

rho0 = [0;0;1]; % Initial state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NT_PCPA = 500; % Set the number of points to be used for the ode45 
               % reference solution and for PCPA (piecewise-constant propagator approximation)
tol = 1e-13;   % Relative and absolute tolerance for reference solution by ode45

% One additional time point must be added to PCPA and ode45 methods 
% to run comparisons at identical time points with PS Trap and PS Simp
[rho_ODE,Time_ODE,elapsedTime_ODE] = one_spin_bloch_cart_ODE(para,TMAX,NT_PCPA+1,tol,rho0);
fprintf('\n %s Finished in %5.3f seconds, Tolerance  = %e\n', 'ODE45', elapsedTime_ODE, tol);

figure(1);
sgtitle('Propagation via ODE45')
for i=1:3
    subplot(3,1,i)
    plot(Time_ODE*1000,real(rho_ODE(i,:)),'b');
    ylim([-1 1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho_PCPA,Time_PCPA,elapsedTime_PCPA] = one_spin_bloch_PCPA(para,0.001,NT_PCPA+1,rho0);
fprintf('\n %s Finished in %5.3f seconds, NT = %i', 'PCPA', elapsedTime_PCPA, NT_PCPA);

figure(2);
sgtitle('Propagation via PCPA')
for i=1:3
    subplot(3,1,i)
    plot(Time_PCPA*1000,real(rho_PCPA(i,:)),'b');
    ylim([-1 1])
end

% Compute the relative error on the Frobenius scalar product between
% density matrices
f = plotscalar(rho_PCPA,rho_ODE);
fprintf('\n Average accuracy of %s : %e\n\n', 'PCPA', sum(f)./length(f))
figure(3)
sgtitle('Accuracy of PCPA')
plot(1:length(f),real(f),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PATH-SUMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPrec = 200; % Number of points per path-sum interval
NInt = 1; % Number of time sub-intervals

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMPSON QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho_SIMP,Time_PS_SIMP,elapsedTime_PS] = one_spin_bloch_PS_Simp(para,1000*TMAX,NPrec,NInt,rho0);
fprintf('\n %s Finished in %5.3f seconds, NPrec = %i, NInt = %i', 'Path Sum Simpson', elapsedTime_PS, NPrec, NInt);
rho_PS_SIMP(1,:) = reshape(rho_SIMP(1,1,:),1,[]);
rho_PS_SIMP(2,:) = reshape(rho_SIMP(2,1,:),1,[]);
rho_PS_SIMP(3,:) = reshape(rho_SIMP(3,1,:),1,[]);

figure(4);
sgtitle('Propagation via Path-Sum Simpson')
for i=1:3
    subplot(3,1,i)
    plot(Time_PS_SIMP,real(rho_PS_SIMP(i,:)),'b');
    ylim([-1 1])
end

% Compute the relative error on the Frobenius scalar product between
% density matrices
if NT_PCPA==NInt*NPrec
    f = plotscalar(rho_PS_SIMP,rho_ODE);
    fprintf('\n Average accuracy of %s : %e\n\n', 'Path Sum Simpson', sum(f)./length(f))
    figure(5)
    sgtitle('Accuracy of PS Simpson')
    plot(1:length(f),real(f),'r')
else
    display('WARNING Number of points for ode45 must be equal to NInt*NPrec for CORRECT accuracy evaluation')
end


%%%%%%%%%%%%%%%%%%%%%%%%%% TRAPEZOIDAL QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%
[rho_TRAP,Time_PS_TRAP,elapsedTime_PS] = one_spin_bloch_PS_TRAP(para,1000*TMAX,NPrec,NInt,rho0);
fprintf('\n %s Finished in %5.3f seconds, NPrec = %i, NInt = %i', 'Path Sum Trapezoidal', elapsedTime_PS, NPrec, NInt);
rho_PS_TRAP(1,:) = reshape(rho_TRAP(1,1,:),1,[]);
rho_PS_TRAP(2,:) = reshape(rho_TRAP(2,1,:),1,[]);
rho_PS_TRAP(3,:) = reshape(rho_TRAP(3,1,:),1,[]);

figure(6);
sgtitle('Propagation via Path-Sum Trapezoidal')
for i=1:3
    subplot(3,1,i)
    plot(Time_PS_TRAP,real(rho_PS_TRAP(i,:)),'b');
    ylim([-1 1])
end

% Compute the relative error on the Frobenius scalar product between
% density matrices
if NT_PCPA==NInt*NPrec
    f = plotscalar(rho_PS_TRAP,rho_ODE);
    fprintf('\n Average accuracy of %s : %e\n\n', 'Path Sum Trapezoidal', sum(f)./length(f))
    figure(7)
    sgtitle('Accuracy of PS Trapezoidal')
    plot(1:length(f),real(f),'r')
else
    display(' ! WARNING Number of points for ode45 must be equal to NInt*NPrec for CORRECT accuracy evaluation !')
end


function E = plotscalar(a,b)
%E is the 1 minus the normalised scalar product between vectors 
%V = a(i) and W = b(i) for i running on all time points.
ll=size(a);ll=ll(2);
E=zeros(1,ll);

for i=1:ll
    V = a(:,i);
    W = b(:,i);
    E(i) = V'*W./(sqrt(V'*V)*sqrt(W'*W));
end
E = 1-E;

end

