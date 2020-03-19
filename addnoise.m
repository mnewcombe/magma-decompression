function H2Onoise = addnoise(H2O_meas, sig)
H2O = H2O_meas;
H2Onoise = zeros(length(H2O), 1);
H2Onoise(:,1) = H2O(:,1);

for i = 1:length(H2O)
      
    H2Onoise(i) = H2O(i) + normrnd(0, sig);

end

dlmwrite('H2Onoise.txt', H2Onoise)

end