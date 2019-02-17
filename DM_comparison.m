% Nrep=100; required number of sets containing emulated warm peaks, for each location


function [Part_Mean, Part_Std] = DM_comparison(Design_points, index_lr, Coeff, Std, Dnu, Nrep)

PC_Val=dlmread('Data/Files in txt or csv/PC_Values_interpolated_for_data_locations.txt');
Mn_Val=dlmread('Data/Files in txt or csv/Mn_Values_interpolated_for data_locations.txt');
PI_Val=dlmread('Data/Files in txt or csv/Model_PI_Values_interpolated_for data_locations.txt');

[Coord, Plio_H, Std_H, Peak_size, Sample_Size, PI_H, ~, ~, ~] = read_Harry_data;
Nloc = size(Coord,1);
deg = size(Design_points,1) - 4;
r = 10;
Part_Mean = zeros(Nloc,1);
Part_Std  = zeros(Nloc,1);

rng(23478);

for loc=1:Nloc
    
    Pc = PC_Val(loc,:);
    Mn = Mn_Val(loc);
    Ns = Sample_Size(loc);
    Np = Peak_size(loc);
    
    if Np>1 && Ns/Np > 2 && Ns/Np <6 && loc~=45
        
        Peaks = emulated_WPA_loc(Design_points, index_lr, Coeff, Std, Dnu, r, deg, Pc, Mn, Ns, Np, 0, Nrep);
        Peaks = Peaks - PI_Val(loc);   % take the difference with pre-industrial
        
        Mn_Hr = Plio_H(loc) - PI_H(loc);       % Mean, Harry
        Std_Hr  = Std_H(loc);                % Variance, Harry    
        Means_Emul     = mean(Peaks,1);        % row vector 1 x Nrep
        StDevs_Emul = std(peaks, 0, 1);     % row vector 1 x Nrep
        Mn_Em  = mean(Means_Emul);
        Std_Em = mean(StDevs_Emul);
        
        Part_Mean(loc) = (Mn_Em - Mn_Hr)/sqrt(var(Means_Emul) + (0.75/Np) );
        Part_Std(loc) = (Std_Em - Std_Hr)/sqrt(2*var(StDevs_Emul));
        
    else
        Part_Mean(loc) = nan;
        Part_Std(loc) = nan;
    end
    
    disp(['Location number ' num2str(loc) ' completed.']);
end
    %Compat_Measure =  sqrt( (Part_Mean+Part_Std)/2 );

end


    
    
    