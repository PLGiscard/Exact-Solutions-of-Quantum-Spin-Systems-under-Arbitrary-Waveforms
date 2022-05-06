%% M. Foroozandeh, P.-L. Giscard, 04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% para = paragen_2(Omega,J,alpha,DeltaF,taup,Phi0,deltat,deltaf,n)
para = paragen_2([-2*pi*600,2*pi*700],150,180,50000,0.001,0,0.0005,0,30);
TMAX = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE45 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NT_PCPA = 500; % Set the number of points to be used for the ode45 
               % reference solution and for PCPA (piecewise-constant propagator approximation)
tol = 1e-13;   % Relative and absolute tolerance for reference solution by ode45

% One additional time point must be added to PCPA and ode45 methods 
% to run comparisons at identical time points with PS Trap and PS Simp
[rho_ODE,Time_ODE,elapsedTime_ODE] = two_spin_ODE(para,TMAX,NT_PCPA+1,tol);
fprintf('\n %s Finished in %5.3f seconds, Tolerance  = %e\n\n', 'ODE45', elapsedTime_ODE, tol);

figure(1);
sgtitle('Propagation via ODE45')
for i=1:4
    for j=1:4
        subplot(4,4,j+4*(i-1))
        plot(Time_ODE*1000,real(squeeze(rho_ODE(i,j,:))),'b');
        ylim([-1 1])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho_PCPA,Time_PCPA,elapsedTime_PCPA] = two_spin_PCPA(para,TMAX,NT_PCPA+1);
fprintf('\n %s Finished in %5.3f seconds, NT = %i', 'PCPA', elapsedTime_PCPA, NT_PCPA);

figure(2);
sgtitle('Propagation via PCPA')
for i=1:4
    for j=1:4
        subplot(4,4,j+4*(i-1))
        plot(Time_PCPA*1000,real(squeeze(rho_PCPA(i,j,:))),'b');
        ylim([-1 1])
    end
end

f = plotfrob(rho_PCPA,rho_ODE);
fprintf('\n Average accuracy of %s : %e\n\n', 'PCPA', sum(f)./length(f))
figure(3)
sgtitle('Accuracy of PCPA')
plot(1:length(f),real(f),'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PATH-SUMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NPrec = 500; % Number of points per path-sum interval
NInt = 1; % Number of time sub-intervals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMPSON QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%
if NPrec<4
    display(' WARNING NPrec must be greater than or equal to 4 ')
end

[rho_PS_SIMP,Time_PS_SIMP,elapsedTime_PS]=two_spin_PS_SIMP(para,TMAX*1000,NPrec,NInt);
fprintf('\n %s Finished in %5.3f seconds, NPrec = %i, NInt = %i', 'Path Sum Simpson', elapsedTime_PS, NPrec, NInt);

figure(4);
sgtitle('Propagation via Path-Sum Simpson')
for i=1:4
    for j=1:4
        subplot(4,4,j+4*(i-1))
        plot(Time_PS_SIMP*1000,real(squeeze(rho_PS_SIMP(i,j,:))),'b');
        ylim([-1 1])
    end
end

if NT_PCPA==NInt*NPrec
    f = plotfrob(rho_PS_SIMP,rho_ODE);
    fprintf('\n Average accuracy of %s : %e\n\n', 'Path Sum Simpson', sum(f)./length(f))
    figure(5)
    sgtitle('Accuracy of PS Simpson')
    plot(1:length(f),real(f),'r')
else
    display(' WARNING Number of points for ode45 must be equal to NInt*NPrec for CORRECT accuracy evaluation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%% TRAPEZOIDAL QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%
[rho_PS_TRAP,Time_PS_TRAP,elapsedTime_PS]=two_spin_PS_TRAP(para,TMAX*1000,NPrec,NInt);
fprintf('\n %s Finished in %5.3f seconds, NPrec = %i, NInt = %i', 'Path Sum Trapezoidal', elapsedTime_PS, NPrec, NInt);


figure(6);
sgtitle('Propagation via Path-Sum Trapezoidal')
for i=1:4
    for j=1:4
        subplot(4,4,j+4*(i-1))
        plot(Time_PS_TRAP*1000,real(squeeze(rho_PS_TRAP(i,j,:))),'b');
        ylim([-1 1])
    end
end

if NT_PCPA==NInt*NPrec
    f = plotfrob(rho_PS_TRAP,rho_ODE);
    fprintf('\n Average accuracy of %s : %e\n\n', 'Path Sum Trapezoidal', sum(f)./length(f))
    figure(7)
    sgtitle('Accuracy of PS Trapezoidal')
    plot(1:length(f),real(f),'r')
else
    display('WARNING Number of points for ode45 must be equal to NInt*NPrec for CORRECT accuracy evaluation')
end

    

function f=frob(A,B)
f=trace(A'*B); % Frobenius scalar product between matrices A and B
end

function f=plotfrob(a,b)
%f is the 1 minus the normalised Frobenius scalar product between matrices
%  A = a(i) and B = b(i) for i running on all time points.
ll=size(a);ll=ll(3);
f=zeros(1,ll);

for i=1:ll
    A = a(:,:,i);
    B = b(:,:,i);
    f(i) = frob(A,B)/sqrt(frob(A,A)*frob(B,B));
end
f=1-f;

end

