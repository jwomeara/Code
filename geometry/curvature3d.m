function [curvature] = curvature3d( phi,Nx,dx,Ny,dy,Nz,dz,ghostcell_width )  

  % fill boundary cells
  phi(:,1:ghostcell_width,:) = ...
    phi(:,Nx+1:ghostcell_width+Nx,:);
  phi(:,Nx+ghostcell_width+1:end,:) = ...
    phi(:,ghostcell_width+1:2*ghostcell_width,:);
  phi(1:ghostcell_width,:,:) = ...
    phi(Ny+1:ghostcell_width+Ny,:,:);
  phi(Ny+ghostcell_width+1:end,:,:) = ...
    phi(ghostcell_width+1:2*ghostcell_width,:,:);
  phi(:,:,1:ghostcell_width) = ...
    phi(:,:,Nz+1:ghostcell_width+Nz);
  phi(:,:,Nz+ghostcell_width+1:end) = ...
    phi(:,:,ghostcell_width+1:2*ghostcell_width);

  % compute normal velocity
  phi_x  = (phi(3:end,:,:) - phi(1:end-2,:,:))/2/dx;
  phi_y  = (phi(:,3:end,:) - phi(:,1:end-2,:))/2/dy;
  phi_z  = (phi(:,:,3:end) - phi(:,:,1:end-2))/2/dz;
  phi_xy = (phi_y(3:end,:,:) - phi_y(1:end-2,:,:))/2/dx;
  phi_xz = (phi_z(3:end,:,:) - phi_z(1:end-2,:,:))/2/dx;
  phi_yz = (phi_z(:,3:end,:) - phi_z(:,1:end-2,:))/2/dy;
  
  phi_xx = (phi(3:end,:,:) - 2*phi(2:end-1,:,:) + phi(1:end-2,:,:))/dx/dx;
  phi_yy = (phi(:,3:end,:) - 2*phi(:,2:end-1,:) + phi(:,1:end-2,:))/dy/dy;
  phi_zz = (phi(:,:,3:end) - 2*phi(:,:,2:end-1) + phi(:,:,1:end-2))/dz/dz;
  
  phi_x  = phi_x(ghostcell_width:end-ghostcell_width+1, ...
                 1+ghostcell_width:end-ghostcell_width, ...
                 1+ghostcell_width:end-ghostcell_width);
  phi_y  = phi_y(1+ghostcell_width:end-ghostcell_width, ...
                 ghostcell_width:end-ghostcell_width+1, ...
                 1+ghostcell_width:end-ghostcell_width);
  phi_z = phi_z( 1+ghostcell_width:end-ghostcell_width, ...
                 1+ghostcell_width:end-ghostcell_width, ...
                 ghostcell_width:end-ghostcell_width+1 );
  
            
  phi_xy = phi_xy(ghostcell_width:end-ghostcell_width+1, ...
                  ghostcell_width:end-ghostcell_width+1, ...
                  1+ghostcell_width:end-ghostcell_width); %no deriv on z
  phi_xz = phi_xz(ghostcell_width:end-ghostcell_width+1, ...
                  1+ghostcell_width:end-ghostcell_width, ...%no deriv on y
                  ghostcell_width:end-ghostcell_width+1); 
  phi_yz = phi_yz(1+ghostcell_width:end-ghostcell_width, ...
                  ghostcell_width:end-ghostcell_width+1, ...
                  ghostcell_width:end-ghostcell_width+1);
              
  phi_xx = phi_xx(ghostcell_width:end-ghostcell_width+1, ...
                  1+ghostcell_width:end-ghostcell_width, ...
                  1+ghostcell_width:end-ghostcell_width);
              
  phi_yy = phi_yy(1+ghostcell_width:end-ghostcell_width, ...
                  ghostcell_width:end-ghostcell_width+1, ...
                  1+ghostcell_width:end-ghostcell_width );
  phi_zz = phi_zz(1+ghostcell_width:end-ghostcell_width, ...
                  1+ghostcell_width:end-ghostcell_width, ...
                  ghostcell_width:end-ghostcell_width+1 );
              
              
  abs_grad_phi_sq = (phi_x.*phi_x + phi_y.*phi_y +phi_z.*phi_z);
  abs_grad_phi_cube = abs_grad_phi_sq.^1.5;

  
  curvature = ( phi_x.*phi_x.*phi_yy ...
              - 2*phi_x.*phi_y.*phi_xy ...
              + phi_y.*phi_y.*phi_xx ...
              + phi_x.*phi_z.*phi_zz ...
              - 2*phi_x.*phi_z.*phi_xz ...
              + phi_z.*phi_z.*phi_xx ...
              + phi_y.*phi_y.*phi_zz ...
              - 2*phi_y.*phi_z.*phi_yz ) ...
                        ./ ...% grad_phi_sq.^1.5;
                    (abs_grad_phi_cube + (abs_grad_phi_cube == 0));
