function [r_total,r_cell_contribution] = corrcoeff_contribution_percell(odor1,odor2)

Ncell = size(odor1,2);
Ntime = size(odor1,1);
r = corrcoef(odor1(:),odor2(:));
r_total= r(1,2);
r_cell_contribution = zeros(Ncell,1);
for i = 1:Ncell
    odor1_removed_cell = odor1;
    odor1_removed_cell(:,i) = [];
    odor2_removed_cell = odor2;
    odor2_removed_cell(:,i) = [];
    r = corrcoef(odor1_removed_cell(:),odor2_removed_cell(:));
    r_cell_contribution(i) = r_total - r(1,2);

end

    
    
