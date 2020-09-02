function phi2 = zerogradientBC2(phi,nx,my)

phi(2:my-1,1) = phi(2:my-1,2);
phi(2:my-1,nx) = phi(2:my-1,nx-1);
phi(1,2:nx-1) = phi(2,2:nx-1);
phi(my,2:nx-1) = phi(my-1,2:nx-1);

phi(1,1) = (phi(2,1)+phi(1,2))/2; %phi(2,2);
phi(1,nx) = (phi(1,nx-1)+phi(2,nx))/2;%phi(2,nx-1);
phi(my,1) = (phi(my-1,1)+phi(my,2))/2;%phi(my-1,2);
phi(my,nx) = (phi(my,nx-1)+phi(my-1,nx))/2;%phi(my-1,nx-1);

phi2 = phi;
