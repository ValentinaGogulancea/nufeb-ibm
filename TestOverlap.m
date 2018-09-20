function new_pos = TestOverlap(bac_x,bac_y,bac_r,bac_m,overlap,maxxSys,dx)

X = [bac_x bac_y]; % the coordinates of bacteria

Y = X;

r = max(bac_r); % the distance for overlap - as it is the maximum radius

sp = [];

[idx,D] = rangesearch(X,Y,2*r); % initial determining of the overlappings

% assembling a better matrix from the one outputted by rangesearch

for i = 1:length(D)
    D_num = cell2mat(D(i));
    index = cell2mat(idx(i));
    for j = 2: length(D_num)
        sum_radius = bac_r(index(1)) + bac_r(index(j));
        if sum_radius - D_num(j) > overlap
            sp = [sp;index(1) index(j)];
        end
    end
end

if ~isempty(sp)
    flag = sparse(sp(:,1),sp(:,2),1,length(bac_x),length(bac_x));
else
    flag = sparse(length(bac_x),length(bac_x));
end


while nnz(flag)~=0
    for k = 1:length(sp(:,1))
        % 
        if flag(sp(k,1),sp(k,2))
        flag(sp(k,1),sp(k,2)) = 0;
        overlap_x = bac_x(sp(k,1)) - bac_x(sp(k,2));
        overlap_y = bac_y(sp(k,1)) - bac_y(sp(k,2));
        % what the distance is
        dist_overlap = sqrt(overlap_x^2 + overlap_y^2);
        % the sum of the radii
        r_overlap =  (bac_r(sp(k,1)) +  bac_r(sp(k,2)));
        % normalized overlap
        
        r_norm = (r_overlap - dist_overlap)/dist_overlap;
        if r_norm < 1e2
        
        % computing how much to either side one must move
        a1 = bac_m(sp(k,1))/(bac_m(sp(k,1)) + bac_m(sp(k,2)));
        %
        a2 = bac_m(sp(k,2))/(bac_m(sp(k,1)) + bac_m(sp(k,2)));
        %
        bac_x(sp(k,1)) = bac_x(sp(k,1)) + r_norm * a1 * overlap_x;
        bac_y(sp(k,1)) = bac_y(sp(k,1)) + r_norm * a1 * overlap_y;
        %
         if  bac_x(sp(k,1)) < 0 
             bac_x(sp(k,1)) = bac_x(sp(k,1)) - r_norm * a1 * overlap_x;
         end
         if  bac_x(sp(k,2)) < 0 
             bac_x(sp(k,2)) = bac_x(sp(k,2)) + r_norm * a2 * overlap_x;
         end
        if  bac_y(sp(k,1)) < dx 
             bac_y(sp(k,1)) = bac_y(sp(k,1)) - r_norm * a1 * overlap_y;
        end
        bac_x(sp(k,2)) = bac_x(sp(k,2)) - r_norm * a2 * overlap_x;
        bac_y(sp(k,2)) = bac_y(sp(k,2)) - r_norm * a2 * overlap_y;
        
        if  bac_x(sp(k,2)) > maxxSys 
             bac_x(sp(k,2)) = bac_x(sp(k,2)) + r_norm * a2 * overlap_x;
        end
        if  bac_y(sp(k,2)) < dx 
             bac_y(sp(k,2)) = bac_y(sp(k,2)) + r_norm * a2 * overlap_y;
        end
        flag(sp(k,2),sp(k,1)) = 0;
        end  
        end
    end
   
    X = [bac_x bac_y];
    [sp, flag] = checkOVLP(X,r,bac_r,overlap);
   
end

new_pos = X;




