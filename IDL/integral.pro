	function integral,hH,hM,dt,g,f
	  COMMON grid,nx,my,dx,dy
	  COMMON velocity,u,v
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	dd = 1
	dtt = double(dt)

	uphi2 = dblarr(nx,my)
	vphi2 = dblarr(nx,my)
	u2 = dblarr(nx,my)
	v2 = dblarr(nx,my)

	hphi2 = dblarr(nx,my)
	hhM = dblarr(nx,my)
	hux = dblarr(nx,my)
	huy = dblarr(nx,my)
	h2 = dblarr(nx,my)

	;for t = 1,3600*72/dt do begin

	  uphi2 = u
	  vphi2 = v

	  uphi2 = advection(uphi2,u,v,dtt)
	  vphi2 = advection(vphi2,u,v,dtt)

	  for i = 1,nx-2 do begin
	  for j = 1,my-2 do begin

	    case i of

	      1 : u2(i,j) = uphi2(i,j) - dtt*g*((hH(i+1,j)-hH(i,j))/(dx)) + dtt*f*vphi2(i,j)

	      nx-2 : u2(i,j) = uphi2(i,j) - dtt*g*((hH(i,j)-hH(i-1,j))/(dx)) + dtt*f*vphi2(i,j)

	      else : u2(i,j) = uphi2(i,j) - dtt*g*((hH(i+1,j)-hH(i-1,j))/(2.*dx)) + dtt*f*vphi2(i,j)

	    endcase

	    case j of

	      1 : v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j+1)-hH(i,j))/(dy)) - dtt*f*uphi2(i,j)

	      my-2 : v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j)-hH(i,j-1))/(dy)) - dtt*f*uphi2(i,j)

	      else : v2(i,j) = vphi2(i,j) - dtt*g*((hH(i,j+1)-hH(i,j-1))/(2.*dy)) - dtt*f*uphi2(i,j)

	    endcase

	  endfor
	  endfor

	  u2(0,1:my-2) = u2(1,1:my-2)
          u2(nx-1,1:my-2) = u2(nx-2,1:my-2)
          u2(1:nx-2,0) = u2(1:nx-2,1)
          u2(1:nx-2,my-1) = u2(1:nx-2,my-2)

          u2(0,0) = (u2(1,0)+u2(0,1))/2.
          u2(nx-1,0) = (u2(nx-2,0)+u2(nx-1,1))/2.
          u2(0,my-1) = (u2(0,my-2)+u2(1,my-1))/2.
          u2(nx-1,my-1) = (u2(nx-2,my-1)+u2(nx-1,my-2))/2.

          v2(0,1:my-2) = v2(1,1:my-2)
          v2(nx-1,1:my-2) = v2(nx-2,1:my-2)
          v2(1:nx-2,0) = v2(1:nx-2,1)
          v2(1:nx-2,my-1) = v2(1:nx-2,my-2)

          v2(0,0) = (v2(1,0)+v2(0,1))/2.
          v2(nx-1,0) = (v2(nx-2,0)+v2(nx-1,1))/2.
          v2(0,my-1) = (v2(0,my-2)+v2(1,my-1))/2.
          v2(nx-1,my-1) = (v2(nx-2,my-1)+v2(nx-1,my-2))/2.


	  hphi2 = hH

          hphi2 = advection(hphi2,u2,v2,dtt)

	  for i = 0,nx-1 do begin
	  for j = 0,my-1 do begin

	    case i of

	      0 : hux(i,j) = (u2(i+1,j)-u2(i,j))/(dx)

	      nx-1 : hux(i,j) = (u2(i,j)-u2(i-1,j))/(dx)

	    else : hux(i,j) = (u2(i+1,j)-u2(i-1,j))/(2.*dx)

	    endcase

	    case j of

	      0 : huy(i,j) = (v2(i,j+1)-v2(i,j))/(dy)

	      my-1 : huy(i,j) = (v2(i,j)-v2(i,j-1))/(dy)

	    else : huy(i,j) = (v2(i,j+1)-v2(i,j-1))/(2.*dy)

	    endcase

	  endfor
	  endfor

	  for i = 0,nx-1 do begin
	  for j = 0,my-1 do begin

	    h2(i,j) = hphi2(i,j) - dtt*(hphi2(i,j)-hM(i,j))*(hux(i,j)+huy(i,j))

	  endfor
	  endfor

	  u = u2
	  v = v2

	  hH = h2

	 ; if (t mod 3600*72/dt) eq 0 then begin

	 ;   dd = dd + 1
	 ;   hhM = hH - hM 

	 ;   printf,11,hhM
	 ;   print,min(hH)
	
	 ;   printf,12,u
	 ;   printf,13,v

	 ; endif

	;endfor

	return,hH

	end
