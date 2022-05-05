function para = paragen(Omega,alpha,DeltaF,taup,Phi0,deltat,deltaf,n)
%% M. Foroozandeh, P.-L. Giscard, 04/2022

%%%%%%%%%%%%%%%%%%%%%%% Setting beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for a chirped pulse with super-Gaussian amplitude

if alpha >= 180 % flip angle

    Q = 5;

else

    Q = (2/pi)*log(2/(cosd(alpha)+1));

end

para.n = n; % Gaussian profile parameter. n=2 gives the true Gaussian, 
            % n>>1 yields a rectangle-like profile
para.deltat = deltat; % s
para.DeltaF = DeltaF; % Hz
para.deltaf = deltaf; % Hz
para.taup = taup; % s
para.Phi0 = Phi0; % rad
para.Omega = Omega; % rad/s
omega1=2*pi*sqrt(DeltaF*Q/(2*pi*taup));
para.omega1 = omega1; % rad/s

end