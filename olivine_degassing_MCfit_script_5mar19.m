
close all
clear all
global H2O_meas
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

load seguam_ol1.txt
load Solex_PHrelation_Seguam.txt
plottitle = 'Seguam olivine 1';
Solex_P = Solex_PHrelation_Seguam(:,1);
Solex_H = Solex_PHrelation_Seguam(:,2);
H2O_meas = seguam_ol1(:,2); 
dist = seguam_ol1(:,1);
radius = 600;     % radius of xal in microns
T_C = 1070;
P0 = 3250;          % pressure of water saturation, based on MI population
Pf = 20;            % final pressure based on groundmass measurements
logdPbydtmax = 1;
logdPbydtmin = -2;

T_K = T_C+273;
T_hi = T_K+20;
T_lo = T_K-20;
DHhi = 9.6e-6*exp(-125000/(8.314*T_hi));
DHlo = 9.6e-6*exp(-125000/(8.314*T_lo));

logDHhi = log10(DHhi);
logDHlo = log10(DHlo);

% Monte Carlo error analysis: Add random noise to the original H2O data to
% create p=# iterations synthetic H2O profiles. Fit each of these profiles.

p=100;

final_results = zeros(p,4+length(H2O)+1);
%results_ordered = zeros(n,4);
%sar = zeros(1,n);

for k = 1:p
k

one_sig = 1.09; %stdev of Saturday Herasil analyses

% Generate fake noisy data    
H2Onoise = addnoise(H2O_meas, one_sig);
H2O = H2Onoise;

crazy=1;
%do-while loop to get rands which are non-crazy
while(crazy==1) %if crazzzzzy, go back and get new vals

    % Select values of log DH and Kd by drawing from prior distributions
    logDH = logDHhi - rand(1)*(logDHhi-logDHlo);  % select a logDH, accounting for temps between T_hi and T_lo
    %Kd = normrnd(7, 1.5);     % partition coefficient *10^4 taken from Le Voyer et al. 2014 
    Kd = 5 + rand(1)*4;         % vary partition coefficient between 4 and 6
   %Kd=7;
    % Generate random starting value for dP/dt and P0
    %n=10;
    logdPbydttest = rand(1)*(logdPbydtmax - logdPbydtmin) + logdPbydtmin;
    x0(:,1) = 10^logdPbydttest;
    %x0(:,2) = rand(1)*(Pmax-Pmin)+Pmin; %starting pressure in bar
    %x0(:,2)=2250;
    H2O_melt_in = spline(Solex_P, Solex_H, P0);
    H2O_ol_test = H2O_melt_in*Kd;
    
    if (H2O_ol_test>max(H2O));
    crazy=0;
    end;   %%ooh i like these randos
    
end

%{
% Calculate misfit of model at each starting position
for i = 1:n
    
    sar(i) = olivineMC(x0(i,:));   %store values of sum of residuals, q1, Temp1, and q2.
end

% Sort results from lowest to highest misfit
results = [sar' x0];
[Y,I]=sort(results(:,1));
results_ordered=results(I,:);  %use the column indices from sort() to sort all columns of A.

clear x0
clear results

% Select the best result with the lowest misfit
x0 = results_ordered(1,2:3);
%}


x0orig = x0;
for j = 1
    
    options = optimset('TolFun',0.1, 'TolX', 0.1);
    [x, fval, exitflag] = fminsearch(@olivineMC, x0, options);
    x0 = x;
    olivineMCplot(x)
end

  final_results(k, 1:6) = [fval x0orig x0 P0 Kd logDH];
  final_results(k, 7:length(H2O)+6) = H2O';
  %final_results(k, end) = exitflag;
  
  clear x
  clear fval
  clear x0
   
end
%%
dlmwrite('final_results_Seguam_ol1_7Mar19.txt', final_results)
  clear x
  clear fval
  clear x0
  
P0 = final_results(:,4);
dPbydt = final_results(:,3);
fval = final_results(:,1);
pointsize = 100;

figure(1)
subplot(2, 2, 2)
plot(dPbydt/10, fval, 'o');
set(gca, 'fontsize', 14)
xlabel('dP/dt (MPa/s)')
ylabel('misfit')
%ylim([0 400])
%h = colorbar;
%ylabel(h, 'misfit', 'o')
axis square

%{
subplot(2, 2, 2)
scatter(final_results(:,6), final_results(:,5), pointsize, fval);
h = colorbar;
ylabel(h, 'misfit')
set(gca, 'fontsize', 14)
xlabel('log D_H (m^2/s)')
ylabel('Kd x 10000')
%}

subplot(2, 2, 3)
plot(dPbydt/10, final_results(:,5), 'o')
%scatter(dPbydt/10, final_results(:,5), pointsize, fval);
%h = colorbar;
%ylabel(h, 'misfit')
set(gca, 'fontsize', 14)
xlabel('dP/dt (MPa/s)')
ylabel('K_d x 10000')
axis square

subplot(2, 2, 4)
%scatter(dPbydt/10, final_results(:,6), pointsize, fval);
%h = colorbar;
%ylabel(h, 'misfit')
plot(dPbydt/10, final_results(:,6), 'o')
set(gca, 'fontsize', 14)
xlabel('dP/dt (MPa/s)')
ylabel('log D_H (m^2/s)')
axis square


set(gcf, 'PaperUnits', 'inches');
orient portrait
papersize = get(gcf, 'PaperSize');
width =7;         % Initialize a variable for width.
height = 7;          % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);
print -dpdf -r600 summary_fig
