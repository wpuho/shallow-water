	function calcuh2,h2,dr,Vmax,b,rvm,f,xc,yc
	  COMMON grid,nx,my,dx,dy
	  COMMON velocity,u,v
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	r = 0.
	Vt = 1.
	n = 1

	while Vt ge 10.^(-4) do begin

	  n = n + 1
	  r = r + dr
	  Vt = Vmax*(r/rvm)*exp((1.-(r/rvm)^b)/b)

	endwhile

	r = dindgen(n)
	Vt = dindgen(n)
	n = 0
	Vt(0) = 1.
	r(0) = 0.

	while Vt(n) ge 10.^(-4) do begin

	  n = n + 1
	  r(n) = r(n-1) + dr
	  Vt(n) = Vmax*(r(n)/rvm)*exp((1.-(r(n)/rvm)^b)/b)

	endwhile

	Vt(0) = 0.
	mm = n

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	d = dblarr(nx,my)
	Vtt = dblarr(nx,my)
	theta = dblarr(nx,my)
	u = dblarr(nx,my)
	v = dblarr(nx,my)

	for i = 0,nx-1 do begin
	;print,i
	for j = 0,my-1 do begin
	
	  rx = double((i+1)-xc)*dx
	  ry = double((j+1)-yc)*dy
	  d(i,j) = sqrt(rx^2 + ry^2)
	    for nn = 0,mm-1 do begin
	      
	      case nn of

	        ((r(nn+1) gt d(i,j)) && (r(nn) lt d(i,j))): Vtt(i,j) = (Vt(nn+1)-Vt(nn))*((d(i,j)-r(nn))/(r(nn+1)-r(nn))) + Vt(n)
	        
		(r(nn) eq d(i,j)): Vtt(i,j) = Vt(nn)
	
	      else: Vtt(i,j) = Vmax*(d(i,j)/rvm)*exp((1.-(d(i,j)/rvm)^b)/b)
	      
	      endcase

	     	
	    endfor

	    theta = atan2(rx,ry)
	    u(i,j) = -1.*Vtt(i,j)*sin(theta)
	    v(i,j) = Vtt(i,j)*cos(theta)

	endfor
	endfor

	h2 = sor(h2,u,v,f,1.)

	return,h2

	end	
