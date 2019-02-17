% This function emulates trajectories for a given location (whose
% corresponding PC scores are provided in input) and returns the peaks of
% the first Ntot trajectories that show exactly NP peaks out of NS randomly
% selected points.

function Final_peaks = emulated_WPA_loc(Design_points, index_lr, Coeff, Std, Dnu, r, deg, PC_locat, Mn_locat, Ns, Np, error_sz, Ntot)

a=-3264; b=-3025; % boundaries for time interval at which emulation takes place

Final_peaks = zeros(Np, Ntot);

for i=1:Ntot
    p=-1;
    j=0;
    while p~=Np
        j=j+1;
        mytimes = a + (b-a)*sort(rand(Ns,1)); % vector of random times between a and b
        New_points=time2data(mytimes);        % generate parameters corresponding to selected times
        % Carry out emulation for the first r PCs, with parameters above
        [M_Pc, Cov_Pc, ~]= emulation_PCscores(Design_points, index_lr, Coeff, r, New_points, 'matern52', 'cov', Dnu);
        % Compute linear combination of previously emulated PCs
        [mu, Cov] = emul_complete(M_Pc, Cov_Pc, PC_locat, Mn_locat, Std, r, 'cov');
        % Take a sample from the emulated trajectory, at the times randomly chosen in mytimes
        y=multiv_student(mu, Cov, deg, 1);
        % Perturb the sample of uniform errors of prescribed size provided in input
        y = y + error_sz*(2*rand(Ns,1)-1);
        pks = extract_peaks(y); % see custom function
        p=length(pks);
    end
    if mod(i,100)==0
        disp(['Number of stored trajectories: ' num2str(i)]);
    end
    Final_peaks(:,i)=pks;
end

end
