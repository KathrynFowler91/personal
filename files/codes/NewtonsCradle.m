function [t, theta, Lp] = NewtonsCradle()

% Newton's Cradle
%
% Solves the equation governing simple harmonic motion for a series of damped pendulums
%
% Created by Kathryn Fowler, PhD Student, The University of Manchester
% KathrynFowler91@outlook.com
% [29/5/2017]

%-------------------------------------------------------------------------%

% Variables
Np = 15;                % number of pendulums [#]
time_in_phase = 30;     % time for the pendulums to come back into phase [s]
num_in_phase = 30;      % number of swings for the pendulums to come back into phase [#]
dt = 0.05;               % time step [s]
time = time_in_phase; % total time [s]
t = (0:dt:time);        % time [s]
g = 9.81;               % gravity [m^2.s^-1]
b = 0.00001;                % damping constant []
m = 0.01;               % mass of bob [kg]
beta = b/m;             % dampening per unit mass [kg^-1]

% Initial Conditions
theta0 = pi/6;          % initial angle of pendulums [radians]
thetadot0 = 0;          % initial angular velocity of pendulums [radians]

% Model initiation
Nt = length(t);         % number of time steps
theta = zeros(Nt,Np);   % solution matrix

% Define lengths of pendulums
np = (num_in_phase:1:num_in_phase+Np-1); % number of swings for each pendulum
Lp = (g-beta.^2./4).*(time_in_phase./2./pi./np).^2; % length of each pendulum


% Solutions for Np pendulums

for i = 1:Np
   
   % dampening cases
   
   % no damping
   if beta==0;                  
       lamda = sqrt(g/Lp(i));
       C2 = atan(thetadot0/theta0/lamda);
       C1 = theta0/cos(C2);
       theta(:,i) = C1.*cos(lamda.*t-C2);
   
   % under damped
   elseif beta^2<4*g*Lp(i);     
       lamda0 = sqrt(g/Lp(i));
       gamma = -(beta/2/Lp(i));
       lamda = sqrt(lamda0^2-gamma^2);
       C2 = atan(thetadot0/theta0/lamda/gamma);
       C1 = theta0/cos(C2);
       theta(:,i) = C1.*exp(gamma.*t).*cos(lamda.*t-C2);
       
   % critically damped
   elseif beta^2==4*g*Lp(i);
       gamma = -beta/2/Lp(i);
       C1 = theta0;
       C2 = thetadot0;
       theta(:,i) = C1.*exp(gamma.*t)+C2.*t.*exp(gamma.*t);
       
   % over damped
   elseif beta^2>4*g*Lp(i);
       lamda_pos = -beta+sqrt(beta^2/4/Lp(i)^2-g/Lp(i));
       lamda_neg = -beta-sqrt(beta^2/4/Lp(i)^2-g/Lp(i));
       C2 = (thetadot0+theta0*lamda_neg)/(lamda_pos-lamda_neg);
       C1 = theta0-C2;
       theta(:,i) = C1.*exp(lamda_pos.*t)+C2.*t.*exp(lamda_neg.*t);
       
   else
       fprintf('\n no damping case \n')
   end
   
end

end

%-------------------------------------------------------------------------%

% Created by Kathryn Fowler, PhD Student, The University of Manchester
% KathrynFowler91@outlook.com
% [29/5/2017]

%-------------------------------------------------------------------------%

