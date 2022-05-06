function [rho,t_array,elapsedTime] = one_spin_bloch_PCPA(para,TMAX,NT,rho0)
%% M. Foroozandeh, P.-L. Giscard, 04/2022
% Propagates the initial state through PCPA (piecewise-constant propagator approximation)

t_array = linspace(0,TMAX,NT);
tres = t_array(2);

smfactor = para.n;
offs_f = para.deltaf;
bandwidth = para.DeltaF;
phi0 = para.Phi0;
tau_p = para.taup;
Omega = para.Omega;
omega1 = para.omega1;
offs_t = para.deltat;

tic,

waveform=chirp_fun(tau_p,bandwidth,phi0,omega1,offs_t,offs_f,smfactor,t_array);

rho_0 = rho0; % initial state

% This takes the offset and pulse information and run the numerical
% simulation and then plot the time evolution of the elements of the final density matrix

for i=1:length(waveform)
    
    total_rot_mat = Rrod(real(waveform(i)), imag(waveform(i)), Omega, tres);
    
    rho_tau_p = total_rot_mat*rho_0;
    rho(:,i)=rho_tau_p;
    
    rho_0=rho_tau_p; 
    
end

elapsedTime = toc;



end

function pulse=chirp_fun(tau_p,bandwidth,phi0,omega1,offs_t,offs_f,smfactor,t_array)

Cx = (exp(-(2^(smfactor+2))*((t_array-offs_t)/tau_p).^smfactor)).*(omega1*cos(phi0+(pi*bandwidth*((t_array-offs_t).^2)/tau_p)-2*pi*offs_f*(t_array-offs_t)));
Cy = (exp(-(2^(smfactor+2))*((t_array-offs_t)/tau_p).^smfactor)).*(omega1*sin(phi0+(pi*bandwidth*((t_array-offs_t).^2)/tau_p)-2*pi*offs_f*(t_array-offs_t)));

pulse = complex(Cx,Cy);

end

function total_rot_mat = Rrod(Cx_t, Cy_t, Omega, tres)
% Returns the rotation matrix associated with a time point in the pulse
% 
% Uses Rodrigues formula for calculation of the exponential of
% skew-symmetric matrix
%
% Input:
%     - Omega : spin resonance offset (rad/sec)
%     - tres : time resolution
%     - Cx_t and Cy_t : components of the pulse
%
% Output:
%     - total_rot_mat: rotation matrix associated with a time point in the pulse
%

Omega = Omega + eps; % avoid singularity at offs = 0, Cx_t = 0 and Cy_t = 0

R = [0 -Omega Cy_t; Omega 0 -Cx_t;-Cy_t Cx_t 0] * tres;

gamma = tres * norm([Cx_t Cy_t Omega]);

total_rot_mat = eye(3) + (sin(gamma)/gamma) * R + ((1-cos(gamma))/(gamma^2)) * R^2;

end
