% This function emulates trajectories for a given location (whose
% corresponding PC scores are provided in input) and returns the peaks of
% the first Ntot trajectories that show exactly NP peaks out of NS randomly
% selected points.



function Final_peaks = emulated_WPA_loc(Design_points, index_lr, Coeff, Std, Dnu, r, deg, PC_locat, Mn_locat, Ns, Np, error_sz, Ntot)

a=-3264; b=-3025; % boundaries for time interval at which emulation takes place

Final_peaks = zeros(Np, Ntot);
%rng(3214);

for i=1:Ntot
    p=-1;
    j=0;
    while p~=Np
        j=j+1;
        mytimes = a + (b-a)*sort(rand(Ns,1)); % vector of random times between a and b
        New_points=time2data(mytimes);
        [M_Pc, Cov_Pc, ~]= emulation_PCscores(Design_points, index_lr, Coeff, r, New_points, 'matern52', 'cov', Dnu);
        [mu, Cov] = emul_complete(M_Pc, Cov_Pc, PC_locat, Mn_locat, Std, r, 'cov');
        y=multiv_student(mu, Cov, deg, 1); % takes a sample from the emulated trajectory,
                                           % at the times randomly chosen in mytimes
        y = y + error_sz*(2*rand(Ns,1)-1);
        pks = extract_peaks(y);
        p=length(pks);
    end
    if mod(i,100)==0
        disp(['Number of stored trajectories: ' num2str(i)]);
    end
    Final_peaks(:,i)=pks;
end

end





%{
figure;
plot(mytimes, y, 'LineWidth', 1);
set(gca, 'FontSize', 20); 
xlim([-3264,-3025]);
ylim([16.6, 18.6]);
title(['Lat=' num2str(Coord(loc,1)) ', Lon=' num2str(Coord(loc,2))], 'FontSize', 18);
ylabel('Temp Delta (degrees)'); xlabel('time [ka]');
%}

