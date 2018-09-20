function [bac_x, bac_y] = bac_shovingloops(bac_x, bac_y, bac_r, bac_m, k, s_dist, overlap, dx, maxxSys)
[sbac_y, I] = sort(bac_y);
sbac_x = bac_x(I);
sbac_r = bac_r(I);
sbac_m = bac_m(I);

bac_n = length(sbac_x);
%%% Shoving until a stable system is achieved
shov = 1;       % boolean test for more shoving (1: shoving needed; 0: shoving ready)

xr = zeros(bac_n,1);
yr = zeros(bac_n,1);


while shov == 1 

    %%% One shoving step
    shov = 0;   % start assuming there is no need for shoving
    for i=1:bac_n-1  % browse n-1 cells
        for j=i+1:bac_n       % search overlapping cells among all other cells, "facing forward" ...
            dxx = sbac_x(i,1) - sbac_x(j,1);
            if abs(dxx) < s_dist
                dyy = sbac_y(i,1) - sbac_y(j,1);
                if abs(dyy) < s_dist
                    d = sqrt(dxx*dxx+dyy*dyy + 1e-20); % current distances between cell centers
                    % calculate overlap r0
                    r0 = k*(sbac_r(i,1) + sbac_r(j,1)) - d;
                    r0d = r0/d;
                    dix = dxx * r0d;
                    diy = dyy * r0d;

                    if r0 > overlap
                        shov = 1;
                        a1 = sbac_m(i)/(sbac_m(i)+sbac_m(j)); 
                        a2 = sbac_m(j)/(sbac_m(i)+sbac_m(j));
                        a1aux = a1; a2aux = a2;
                        if sbac_y(i) <= dx
                            a1aux = 1e-40/(1e-40 + sbac_m(j)); 
                        end
                        if sbac_y(j) <= dx
                            a2aux = 1e-40/(1e-40 + sbac_m(i)); 
                        end
                        xr(i) = xr(i) - dix*a1;
                        yr(i) = yr(i) - diy*a1aux;
                        xr(j) = xr(j) + dix*a2;
                        yr(j) = yr(j) + diy*a2aux;
                    end %if r0
                end %if
            end
        end %for j
    end %for i

    for i=1:bac_n  % browse n cells
        sbac_x(i) = sbac_x(i) - xr(i,1);
        sbac_y(i) = sbac_y(i) - yr(i,1);
        if sbac_x(i) < 0 
            sbac_x(i) = sbac_x(i) + xr(i,1);
        elseif sbac_x(i) > maxxSys 
            sbac_x(i) = sbac_x(i) + xr(i,1);
        end
        if sbac_y(i) < dx
            sbac_y(i) = sbac_y(i) + yr(i,1);
        end
        xr(i) = 0.0; % the coordinates resultant vector
        yr(i) = 0.0; % the coordinates resultant vector
    end %for i
end %while
bac_x(I) = sbac_x;
bac_y(I) = sbac_y;
end