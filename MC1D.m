%This Script is a simple 1D Electron Simulation tool
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')
global C
global x px Vx Vpx
global Efield
global nParticles
global Vvector

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²

% nParticles = 5;
% x(1, :) = zeros(1,nParticles);
% Vpx(1, :) = zeros(1, nParticles);
% Vx(1, :) = zeros(1, nParticles);
% px(1, :) = zeros(1,nParticles);
% 
% Efield = rand(1,nParticles)*0.4e07;

x = 0;
px = 0;
Vpx = 0;
Vx = 0;
Efield = 0.4e01;

tSim = 400;

%t = 0;
dt = 1e-6;
vdrift = 0;
vdrift_o = 0;
vavg = 0;

Vvector;

for i = 2:tSim
%     Vx(1:nParticles) = Vpx(1:nParticles) + ((C.q_0 .* Efield)./C.m_0 )*(dt);
%     ivp = ones(1,nParticles)*(i-1);
%     iv = ones(1,nParticles)*i;
%     x(1:nParticles) = px(1:nParticles) + (Vx(1:nParticles) .*dt);
    Vx = Vpx + ((C.q_0 * Efield)/C.m_0 )*(dt);
    x = px + Vpx*dt;
    ivp = i-1;
    iv = i;
  % update x
    Pscat = rand(1)*1;
    if(Pscat <= 0.05)
        %Vpx(1, :) = zeros(1, nParticles);
        Vx = 0;
    end
    Vvector = horzcat(Vvector,Vx);
    vavg = mean (Vvector);
    subplot(3,1,1) 
    plot ([ivp iv],[Vpx Vx],'b');
    plot ([ivp iv],[vavg-1 vavg], 'k');
    drawnow;
    xlabel('t');
    ylabel('v');
    legend(sprintf('Drift Velocity = %0.3f',vavg));
    hold all
    
    subplot(3,1,2)
    plot ([px x],[Vpx Vx],'r');
    plot ( [px x],[vavg-1 vavg], 'k');
    drawnow;
    xlabel('x');
    ylabel('v');
    legend(sprintf('Drift Velocity = %0.3f',vavg));
    hold on
    
    subplot(3,1,3)
    plot ([ivp iv],[px x],'m');
    drawnow;
    xlabel('t');
    ylabel('x');
    hold on
    
    
    %Pscat = 
    Vpx = Vx;
    px = x;
 
  
    %Vdrift0 = mean (Vpx);
    %Vdrift = (Vdrift + Vdrift0)/2;
end


