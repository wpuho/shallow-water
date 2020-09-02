	PRO main
	!P.BACKGROUND='FFFFFF'xl
	!P.COLOR='000000'xl
	COMMON grid,nx,my,dx,dy
	COMMON velocity,u,v

;;;;;;;;;;;;;;;;;;;; NAMELIST ;;;;;;;;;;;;;;;;;;;

	rad = 6378000.
	radx = rad/1000.
	pi = acos(-1.)
	d2r = pi/180.

	xc = 261
	yc = 201
	mxc = 151
	myc = 201

	nx = 401
	my = 401
	dx = 5.*1000.
	dy = 5.*1000.

	g = 9.8
	UU = dblarr(nx,my) + (-4.)
	VV = dblarr(nx,my)
	u = dblarr(nx,my)
	v = dblarr(nx,my)

	HH = 5.*1000.
	hmax = 3.*1000.
	sigmax = 80.*1000.
	sigmay = 160.*1000.

	rvm = 100.*1000.
	b = 2.
	Vmax = 20.
	dr = 1000.

	omega = (360./86400.)*d2r
	f = 2.*omega*sin(30.*d2r)

	dt = 5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	hM = dblarr(nx,my)
	hM = calcuhM(hM,hmax,sigmax,sigmay,mxc,myc)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	h1 = dblarr(nx,my)
	h1 = calcuh1(h1,UU,f,g)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	h2 = dblarr(nx,my)
	h2 = calcuh2(h2,dr,Vmax,b,rvm,f,xc,yc)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	hH  = h1 + h2 + HH
	hhM = dblarr(nx,my)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;openr,1,'height.txt'			;
	;readf,1,hhM				;
	;hH = hhM + hM				;
						;
	;openr,2,'u.txt'			;
	;readf,2,u				;
						;
	;openr,3,'v.txt'			;
	;readf,3,v				;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	u = u + UU
	v = v + VV

;	contour,h,nlevels=10,/isotropic

	hhM = hH - hM
	
	openw,11,'h.txt'
	printf,11,hhM
	
	openw,12,'u2.txt'
	printf,12,u

	openw,13,'v2.txt'
	printf,13,v

	print,'minu',min(u)
	print,'minv',min(v)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	timend = 51840 ; 3600*72/5
	dd = 1
	
	for t = 1,timend do begin
	 
	  hH = integral(hH,hM,dt,g,f)

	  if (t mod 4320) eq 0 then begin ; 3600*6/5

	    dd = dd + 1
	    hhM = hH - hM 

	    printf,11,hhM
	    print,min(hH)

	    printf,12,u
	    printf,13,v

	  endif

	endfor

	close,11
	close,12
	close,13

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	openw,12,'h1.txt'
;	printf,12,h1
;	close,12

;	openw,13,'h2.txt'
;	printf,13,h2
;	close,13

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	end
