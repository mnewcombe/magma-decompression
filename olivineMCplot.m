function sum_res_squared = olivineMCplot(x0)
global H2O
global dist
global Solex_P
global Solex_H
global radius
global logDH
global Kd
global P0
global Pf
global plottitle


DH2O = 1e12*10^logDH;  % diffusivity in um2/s

% Specify initial conditions
%P0 = x0(2);  % initial pressure in bars from Lloyd
Pf = 20;     % final pressure in bars from Lloyd
dPbydt = x0(1);   % decompression rate in bars/s from Ferguson
t_tot = (P0-Pf)/dPbydt;

numx = 101;   %number of grid points in x
%numt = 200000;  %number of time steps to be iterated
dx = radius/(numx - 1);
numt = ceil((2*DH2O*t_tot)/(0.9*dx^2)+1);
dt = t_tot/(numt - 1);
t = 0:dt:t_tot;
%dt = 0.01*radius*(dx^2)/DH2O;
%numt = (t_tot/dt)+1;
dP = -dt*dPbydt;
Pvec = P0:dP:Pf;
H2O_melt_vec = spline(Solex_P, Solex_H, Pvec);
H2O_ol_BC = H2O_melt_vec*Kd;     % calculates conc at edge of olivine at each timestep based on measured partition coefficient.

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
%figure(10)
%plot(absdist, H2Omodel, 'o')
residuals = H2O - H2Omodel;
res_squared = residuals.^2;
sum_res_squared = sum(res_squared)/length(H2Omodel);



figure(1)
title(plottitle)
subplot(2, 2, 1)
hold on;
ph1=plot(x,C(:,1), '--k', 'linewidth', 1, 'displayname', 'Initial condition');
plot(-x,C(:,1), '--k', 'linewidth', 1);
ph2=plot(x,C(:,numt), '-r', 'linewidth', 1, 'displayname', 'D_H = 10^-^1^0 m^2/s');
plot(-x, C(:,numt),'-r', 'linewidth', 1);
ph5=plot(dist, H2O, 'ok' ,'markersize', 3, 'linewidth', 1, 'displayname', 'Measured water concentration (ppm)');
xlabel('Radial distance (\mum)');
ylabel('Concentration of water (ppm)');
axis square
set(gca, 'fontsize', 14)

end