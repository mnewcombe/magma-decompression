load final_results_Seguam_ol1_7Mar19.txt
final_results = final_results_Seguam_ol1_7Mar19;
P0 = final_results(:,4);
dPbydt = final_results(:,3);
fval = final_results(:,1);
pointsize = 100;

stats(1) = mean(dPbydt/10);
stats(2) = std(dPbydt/10);
stats(3) = mean(log10(dPbydt/10));
stats(4) = std(log10(dPbydt/10));

dlmwrite('stats.txt', stats)