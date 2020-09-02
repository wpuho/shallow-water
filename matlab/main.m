	clear
	clc

%%%%%%%%%%%%%%%%%%% NAMELIST %%%%%%%%%%%%%%%%%%%%

	rad = 6378000;
	radx = rad/1000;
	d2r = pi/180;

	xc = 261;
	yc = 201;
	mxc = 151;
	myc = 201;

	nx = 401;
	my = 401;
	dx = 5*1000; % m
	dy = 5*1000; % m

	g = 9.8; % m/s
	UU = zeros(my,nx);
	VV = zeros(my,nx);
	UU(:,:) = -4; % m/s
	VV(:,:) = 0;

	H = 5*1000; % m
	hmax = 3*1000; % m
	sigmax = 80*1000; % m
	sigmay = 160*1000;% m

	rvm = 100*1000; 
	b = 2;
	Vmax = 20; % m/s
	dr = 500; % 

	omega = (360/86400)*d2r;
	f = 2*omega*sin(30*d2r);

	dt = 5;

	for i = 1:nx
	  lon(i) = i;
	end
	for j = 1:my
	  lat(j) = j;
	end

	mlon = lon(mxc);
	mlat = lat(myc);
	tylon = lon(xc);
	tylat = lat(yc);

%%%%%%%%%%%%%%%%%%%%%%% hM %%%%%%%%%%%%%%%%%%%%%%%

	hM = zeros(my,nx);

	for i = 1:nx
	  for j = 1:my

	    hM(j,i) = hmax*exp(-((i-mlon)*dx)^2/(sigmax^2)-((j-mlat)*dy)^2/(sigmay^2));

	  end
	end

%%%%%%%%%%%%%%%%%%%%%%% h1 %%%%%%%%%%%%%%%%%%%%%%%

	h1 = zeros(my,nx);
	for i = 1:nx

	h1(1,i) = 0;

	for j = 1:my

	  if (j==1)
	    h1(2,i) = h1(1,i) - f*UU(1,i)*dy/g;
	  elseif(j==my)
	    h1(my,i) = h1(my-1,i) - f*UU(my-1,i)*dy/g;
	  else
	    h1(j+1,i) = h1(j-1,i) - 2*f*UU(j,i)*dy/g;
	  end

	end
	end

%%%%%%%%%%%%%%%%%%%%%%% h2 %%%%%%%%%%%%%%%%%%%%%%

	r(1) = 0; % m
	Vt(1) = 1;
	n = 1;

	while (Vt(n)>=10^-7)
	  n = n + 1;
	  r(n) = r(n-1) + dr;
	  Vt(n) = Vmax*(r(n)/rvm)*exp((1-(r(n)/rvm)^b)/b);
	end

	Vt(1) = 0;
	mm = n;

	d = zeros(my,nx);
	Vtt = zeros(my,nx);
	theta = zeros(my,nx);
	u = zeros(my,nx);
	v = zeros(my,nx);

	for i = 1:nx
	for j = 1:my

	  rx = (i-xc)*dx;
	  ry = (j-yc)*dy;
	  d(j,i) = (rx^2 + ry^2)^0.5;

	  for nn = 1:mm-1

	    if (r(nn+1)>d(j,i) && r(nn)<d(j,i))
	      Vtt(j,i) = (Vt(nn+1)-Vt(nn))*((d(j,i)-r(nn))/(r(nn+1)-r(nn))) + Vt(nn);
	      theta = atan2(ry,rx);
	      u(j,i) = -Vtt(j,i)*sin(theta);
	      v(j,i) = Vtt(j,i)*cos(theta);
	    elseif (r(nn)==d(j,i))
	      Vtt(j,i) = Vt(nn);
	      theta = atan2((j-yc),(i-xc));
	      u(j,i) = -Vtt(j,i)*sin(theta);
	      v(j,i) = Vtt(j,i)*cos(theta);  
	    elseif (d(j,i)==r(mm)) 
	      Vtt(j,i) = Vt(mm);
	      theta = atan2((j-yc),(i-xc));
	      u(j,i) = -Vtt(j,i)*sin(theta);
	      v(j,i) = Vtt(j,i)*cos(theta); 
	    elseif (d(j,i)>r(mm))
	      Vtt(j,i) = Vmax*(d(j,i)/rvm)*exp((1-(d(j,i)/rvm)^b)/b);
	      theta = atan2((j-yc),(i-xc));
	      u(j,i) = -Vtt(j,i)*sin(theta);
	      v(j,i) = Vtt(j,i)*cos(theta); 
	    end

	  end

	end
	end

%%%%%%%%%%%%%%%%%%%% Initial %%%%%%%%%%%%%%%%%%%%%

	hh2 = [];
	Fxy = [];
	h2 = zeros(my,nx);
	h2 = sor2(h2,hh2,nx,my,dx,dy,u,v,Fxy,f,1);
	h1 = sor2(h1,hh2,nx,my,dx,dy,UU,VV,Fxy,f,1);
	h = h1 + h2;

	u = u + UU;
	v = v + VV;

	figure(1)
	h_H = h + H;

	contour(h_H-hM,300)
	colorbar
	hold on
	[j,i] = find(h_H==min(min(h_H)))
	plot(i,j,'kx')
	streamslice(u,v)
	%quiver(u,v);
	dd = 1;
	frames(dd) = getframe(gcf);

%%%%%%%%%%%%%%%%%%%% Integral %%%%%%%%%%%%%%%%%%%

	uphi1 = zeros(my,nx);
	vphi1 = zeros(my,nx);
	uphi2 = zeros(my,nx);
	vphi2 = zeros(my,nx);
	hphi1 = zeros(my,nx);
	hphi2 = zeros(my,nx);
	h_hM2 = zeros(my,nx);
	hux = zeros(my,nx);
	huy = zeros(my,nx);
	hxy = zeros(my,nx);
	u2 = zeros(my,nx);
	v2 = zeros(my,nx);

	for t = 1:3600*72/dt

	  uphi2 = timeintegral2(u,u,v,nx,my,dx,dy,dt);
	  vphi2 = timeintegral2(v,u,v,nx,my,dx,dy,dt);

	  for j = 2:my-1
	  for i = 2:nx-1

	    if (i==2)
	      u2(j,i) = uphi2(j,i) - dt*g*(h_H(j,i+1)-h_H(j,i))/(dx) + dt*f*vphi2(j,i);
	    elseif (i==nx-1)
	      u2(j,i) = uphi2(j,i) - dt*g*(h_H(j,i)-h_H(j,i-1))/(dx) + dt*f*vphi2(j,i);
	    else
	      u2(j,i) = uphi2(j,i) - dt*g*((h_H(j,i+1)-h_H(j,i-1))/(2*dx)) + dt*f*vphi2(j,i);
	    end

	    if (j==2)
	      v2(j,i) = vphi2(j,i) - dt*g*(h_H(j+1,i)-h_H(j,i))/(dy) - dt*f*uphi2(j,i);
	    elseif (j==my-1)
	      v2(j,i) = vphi2(j,i) - dt*g*(h_H(j,i)-h_H(j-1,i))/(dy) - dt*f*uphi2(j,i);
	    else
	      v2(j,i) = vphi2(j,i) - dt*g*((h_H(j+1,i)-h_H(j-1,i))/(2*dy)) - dt*f*uphi2(j,i);
	    end

	  end
	  end

	  u2 = zerogradientBC2(u2,nx,my);
	  v2 = zerogradientBC2(v2,nx,my);

	  hphi2 = timeintegral_h(h_H,hM,u2,v2,nx,my,dx,dy,dt);

	  for j = 1:my
	  for i = 1:nx

	    if (i==1)
	      hux(j,i) = (u2(j,i+1)-u2(j,i))/(dx);
	    elseif (i==nx)
	      hux(j,i) = (u2(j,i)-u2(j,i-1))/(dx);
	    else
	      hux(j,i) = (u2(j,i+1)-u2(j,i-1))/(2*dx);
	    end

	    if (j==1)
	      huy(j,i) = (v2(j+1,i)-v2(j,i))/(dy);
	    elseif (j==my)
	      huy(j,i) = (v2(j,i)-v2(j-1,i))/(dy);
	    else
	      huy(j,i) = (v2(j+1,i)-v2(j-1,i))/(2*dy);
	    end

	  end
	  end

	  for i = 1:nx
	  for j = 1:my

	    h2(j,i) = hphi2(j,i) - dt*(hphi2(j,i)-hM(j,i))*(hux(j,i)+huy(j,i));

	  end
	  end

	  u = u2;
	  v = v2;
	  h_H = h2;

	  if (mod(t,3600*6/dt)==0)

	    figure(t/(3600*6/dt)+1)
	    contour(h_H-hM,300)
	    colorbar
	    hold on
	    [j,i] = find(h_H==min(min(h_H)))
	    plot(i,j,'kx')
	    streamslice(u,v)

	    dd = dd + 1;
	    frames(dd) = getframe(gcf);

	  end

	end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
