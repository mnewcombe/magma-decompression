function sum_res_squared = olivineMC(x0)

global H2O
global dist
global radius
global Solex_P
global Solex_H
global logDH
global Kd
global P0
global Pf

%{
load Fuegob_a_dist_H.txt
load Solex_Fuego_PH.txt
Solex_P = Solex_Fuego_PH(:,1);
Solex_H = Solex_Fuego_PH(:,2);
H2O_meas = Fuegob_a_dist_H(:,2);
H2O = H2O_meas;
dist = Fuegob_a_dist_H(:,1);
radius = 370/2;     % radius of xal in microns
T_C = 1030;
T_K = T_C+273;
T_hi = T_K+20;
T_lo = T_K-20;
DHhi = 9.6e-6*exp(-125000/(8.314*T_hi));
DHlo = 9.6e-6*exp(-125000/(8.314*T_lo));
logDHhi = log10(DHhi);
logDHlo = log10(DHlo);
%logDH = logDHhi - rand(1)*(logDHhi-logDHlo);  % select a logDH, accounting for temps between T_hi and T_lo
logDH = logDHhi;
%Kd = normrnd(7, 3);     % partition coefficient *10^4 taken from Le Voyer et al. 2014
Kd=1;
%x0(:,1) = rand(1)*4 + 2; % dP/dt in bar/s Lloyd 2014 gets 3 - 5 bar/s from embayments, so search 2 - 6 bar/s
%x0(:,2) = rand(1)*(3500-1250)+1250; %starting pressure in bar ranging from 1250 to 3500;
x0(1) = 0.1;
x0(2) = 1500;
%}


%H2O_meas = (1/1.3)*KilaueaIki2_dist_H2O(:,2);   %apply drift correction to SIMS data (from Herasil)
%plot(dist, H2O_meas, 'o')
%x0=3;
DH2O = 1e12*10^logDH;  % diffusivity in um2/s

% Specify initial conditions
%P0 = x0(2);  % initial pressure in bars from Ferguson
Pf = 20;     % final pressure in bars from Ferguson
dPbydt = x0(1);   % decompression rate in bars/s from Ferguson
t_tot = (P0-Pf)/dPbydt;

numx = 101;   %number of grid points in x
%numt = 1000000;  %number of time steps to be iterated
dx = radius/(numx - 1);
numt = ceil((2*DH2O*t_tot)/(0.9*dx^2)+1);
if numt>10^8
    sum_res_squared = 1000000;
elseif numt<10
    sum_res_squared = 1000000;
else
%dt = 0.01*radius*(dx^2)/DH2O;    % sets dt1/dx^2 <0.5 for stability
dt = t_tot/(numt - 1);
t = 0:dt:t_tot;
%dt = 0.01*radius*(dx^2)/DH2O;
%numt = (t_tot/dt)+1;
dP = -dt*dPbydt;
Pvec = P0:dP:Pf;
H2O_melt_vec = spline(Solex_P, Solex_H, Pvec);
H2O_ol_BC = H2O_melt_vec*Kd;     % calculates conc at edge of olivine at each timestep based on measured partition coefficient.
% Filter out values of Kd that produce an initial water concentration lower
% than the maximum measured water concentration


%Neumann = (2*DH2O*dt)/(dx^2);
    
%{
figure(1)
subplot(1, 3, 1)
plot(Pvec, H2O_melt_vec, '-b', 'linewidth', 2)
axis square
xlabel('Pressure (bars)')
ylabel('H_2O in melt (wt%)')
set(gca, 'FontSize', 12)
subplot(1, 3, 2)
plot(t, Pvec, '-b', 'linewidth', 2)
axis square
xlabel('Time (s)')
ylabel('Pressure (bars)')
set(gca, 'FontSize', 12)
subplot(1, 3, 3)
plot(t, H2O_ol_BC, '-b', 'linewidth', 2)
axis square
xlabel('Time (s)')
ylabel('Water at edge of olivine (ppm)')
set(gca, 'FontSize', 12)
%}
x = 0:dx:radius;   %vector of x values, to be used for plotting

C = zeros(numx,numt);   %initialize everything to zero

%specify initial conditions
C(:,1) = H2O_ol_BC(1);

%iterate difference equations
for j=1:numt-1
   t(j+1) = t(j) + dt;

   for i=2:numx-1
      C(i,j+1) = C(i,j) + (DH2O*dt/(dx^2))*(C(i+1,j) - 2*C(i,j) + C(i-1,j)); 
   end
   C(numx,j+1) = H2O_ol_BC(j+1);          % degassing boundary condition
   C(1,j+1) = C(2,j+1);  %C(1,j+1) found from no-flux condition
end

absdist = abs(dist);

H2Omodel = spline(x,C(:,numt),absdist);
residuals = H2O - H2Omodel;
res_squared = residuals.^2;
sum_res_squared = sum(res_squared)/length(H2Omodel);

%{
figure(1);
subplot(2, 2, 1)
hold on;
set(gca, 'fontsize', 12)
ph1=plot(x,C(:,1), '--k', 'linewidth', 2, 'displayname', 'Initial condition');
plot(-x,C(:,1), '--k', 'linewidth', 2);
ph2=plot(x,C(:,numt), '-r', 'linewidth', 2, 'displayname', 'D_H = 10^-^1^0 m^2/s');
plot(-x, C(:,numt),'-r', 'linewidth', 2);
ph5=plot(dist, H2O_meas, 'ok', 'markersize', 10, 'linewidth', 2, 'displayname', 'Measured water concentration (ppm)');
xlabel('Radial distance (\mum)');
ylabel('Concentration of water (ppm)');
%}
end
end
