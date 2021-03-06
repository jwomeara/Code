%KERNEL_DENSITY_ESTIMATE_WITH_REACHABILITY
%k number of nearest neighbours to use
% can formulate X_vec using X_vec = [ X1_vec', X2_vec', X3_vec' ];
% index is the row vector for which we will calculate the kernel density
function [f] = kernel_density_estimate_with_reachability( X_vec,index,radius,k_rd,n)

dim = size( X_vec,2 );
[id_y,dist_y] = neighborsWithinRadius( X_vec, X_vec(index,:), radius );
Y = X_vec(index,:);
X = X_vec(id_y,:);
y_star = cvMahaldist(Y',X');
f = 0;
for i=2:size(id_y,1)        
        [id_x,dist_x] = neighborsWithinRadius( X_vec,X_vec(id_y(i),:), radius );
        if( size(id_x,1) >= n )
            dist_x_sorted = sort( dist_x,'ascend');
            h = dist_x_sorted(n);
            X = X_vec(id_y(i),:);
            Y = X_vec(id_x,:);
            x_star = cvMahaldist( X',Y');
            x_star_sorted = sort(x_star,'ascend');
            d_r = x_star_sorted(k_rd);
            f = f + 1/((2*pi)^(dim/2)*(h)^dim) * exp( -1 * max(d_r^2,y_star(i)^2) / h * (2*d_r^2) );
        else
            f = -0.1;
        end    
end
% 1/m
%q_x = q_x / (size(id,1)-1);
f = f / (size(id_y,1)-1);
%f_fixed_h = f_fixed_h / (size(id,1)-1);

%else
%    q_x = -0.0238;
%    f = -0.0238;
%    f_fixed_h = -0.0238;
%end
end