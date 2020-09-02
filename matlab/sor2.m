	function [h] =sor2(h,h2,nx,ny,dx,dy,uu,vv,fxy,f,omega)
	% Using SOR to solve the Poisson equation: Laplace(h)= F(x,y)
	% for the balanced height h: h(nx,ny)
	% with the source term F(x,y): fxy(nx,ny) 
	% from the velocity u,v: uu(nx,ny) and vv(nx,ny)
	% and lateral b.c.: dh/dx= bc1= f*vv/g, dh/dy= bc2= -f*uu/g 
	%
	% omega: 0 < omega < 2 for stability of the iteration (under testing...)
	% stop iteration when diff is less than 0.001 or iterations > 10000
	%
	g= 9.8;
	n=0;
	diff=100;
	% note: for fastering the convergence, first guess of h needs not to be set to zero,
	% othereise, to zero at the first application...

	if(omega==0) 
	  omega=1;
	end
	%
	%uu = zeros(ny,nx);
	%vv = zeros(ny,nx);
	while( diff>0.001 && n<10000)
	
	  n=n+1;
	  h2(:,:)=h(:,:);

	% compute fxy

	  for j=1:ny
	  for i=1:nx

	    if (i==1)
	      ux= (uu(j,i+1)-uu(j,i))/dx;
	      vx= (vv(j,i+1)-vv(j,i))/dx;
	    elseif (i==nx)
	      ux= (uu(j,i)-uu(j,i-1))/dx;
	      vx= (vv(j,i)-vv(j,i-1))/dx;
	    else
	      ux= (uu(j,i+1)-uu(j,i-1))/(2*dx);
	      vx= (vv(j,i+1)-vv(j,i-1))/(2*dx);
	    end

	    if (j==1)
	      uy= (uu(j+1,i)-uu(j,i))/dy;
	      vy= (vv(j+1,i)-vv(j,i))/dy;
	    elseif (j==ny)
	      uy= (uu(j,i)-uu(j-1,i))/dy;
	      vy= (vv(j,i)-vv(j-1,i))/dy;
	    else
	      uy= (uu(j+1,i)-uu(j-1,i))/(2*dy);
	      vy= (vv(j+1,i)-vv(j-1,i))/(2*dy);
	    end

	    fvo= f*(vx-uy);
	    jco= ux*vy-uy*vx;
	    fxy(j,i)= (fvo+2*jco)/g;

	  end
	  end

	  for j=1:ny
	  for i=1:nx

	    bc1= f*vv(j,i)/g;
	    bc2= -f*uu(j,i)/g;

	    if(i==1) 
	      hL=h(j,i+1)-2*bc1*dx;   % using lateral b.c.
	    else
	      hL=h(j,i-1);
	    end

	    if(j==1) 
	      hD=h(j+1,i)-2*bc2*dy;   % using lateral b.c.
	    else
	      hD=h(j-1,i);
	    end

	    if(i==nx) 
	      hR=h(j,i-1)+2*bc1*dx;   % using lateral b.c.
	    else
	      hR=h(j,i+1);
	    end

	    if(j==ny) 
	      hU=h(j-1,i)+2*bc2*dy;   % using lateral b.c.
	    else
	      hU=h(j+1,i);
	    end

	  % for dx=dy 
	    if (dx==dy)
	      h2(j,i)=h(j,i)+omega*(hU+hD+hR+hL-4*h(j,i)-fxy(j,i)*dx*dx)/4;
	    else
	  % for dx ~= dy
	      ratio=(dx*dx)*(dy*dy)/(2*(dx*dx)+2*(dy*dy));
	      resij=(hR+hL-2*h(j,i))/(dx*dx) + (hU+hD-2*h(j,i))/(dy*dy) - fxy(j,i);
	      h2(j,i)=h(j,i)+omega*ratio*resij;
	    end

	  end
	  end

	  aveh=mean(abs(h2(:)));
	  diff2=abs((h2-h)/aveh);   % using fractional difference
	  diff=max(max(diff2(:,:)));
	  % adjust the threshold for diff when getting slow convergence...

	  %disp(max(max(h2)))
	  %disp(max(max(h)))
	  h=h2;    

	  % extract the constant value to reduce round-off error
	  h(:,:)=h(:,:) - h(1,1);

	  if (mod(n,100) == 0) 
	    fprintf('Iteration n= %d  diff=%f \n ',n,diff)
	  end

	end  

	end
