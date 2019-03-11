function [bac_x, bac_y] = shoving_loops(bac_x, bac_y, bac_r, bac_m, k, s_dist, overlap, dx, maxxSys)
% [sbac_y, I] = sort(bac_y);
% sbac_x = bac_x(I);
% sbac_r = bac_r(I);
% sbac_m = bac_m(I);

X = [bac_x, bac_y ]; % the coordinates of bacteria
Y = X;
r = max(bac_r); % the distance for overlap - as it is the minimum radius, there can be some overlap as is

[idx,D] = rangesearch(X,Y,2*r);

for j = 1:size(D)
    
    D_c = cell2mat(D(j));
    
    if length(D_c) > 1
        
        D_c(1) = [];
        
        index = cell2mat(idx(j));
        
        % checking to see if in the overlaping cases, the problem is the minimum radius used
        % this computes the distance between 2 particles
        % and compares it to the sum of their radii
        % if the distance is smaller than the  - there is overlap but then
        % if they are overlapping - by how much?
        
        sum_rad = k * (bac_r(index(1))*ones(length(index)-1,1) + bac_r(index(2:end)));
        
        A  =  (sum_rad - D_c')./D.c';
        
        for l = 1 : length(A)
            if A(l) > overlap
                overl_dist = A(l);
                % the mass ratio
                mass_1st = bac_m(index(1))*ones(length(index)-1,1);
                mass_neighbours = bac_m(index(2:end));
                aux_m1 = mass_1st./ (mass_1st + mass_neighbours);
                aux_m_neigh = 1 - aux_m1;
                % update the coordinates now
                bac_x(index(1)) = bac_x(index(1)) - sum(aux_m1) * overl_dist * sum(abs(bac_x(index(1)) - bac_x(index(2:end))));
                bac_y(index(1)) = bac_y(index(1)) - sum(aux_m1) * overl_dist * sum(abs(bac_y(index(1)) - bac_y(index(2:end))));
                % bac_x(index(1)) = bac_x(index(1)) - sum(aux_m1) * overl_dist * sum(abs(bac_x(index(1)) - bac_x(index(2:end))));
                bac_x(index(2:end)) = bac_x(index(2:end)) - aux_m_neigh .* overl_dist * (abs(bac_x(index(1)) - bac_x(index(2:end))));
                bac_y(index(2:end)) = bac_y(index(2:end)) - aux_m_neigh .* overl_dist * (abs(bac_y(index(1)) - bac_y(index(2:end))));
            end
        end
    end
end
end

