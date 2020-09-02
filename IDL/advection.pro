	function advection,phi,u,v,dtt
	  COMMON grid,nx,my,dx,dy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	uu = u
	vv = v

	phiu = dblarr(nx,my)

	for j = 0,my-1 do begin
	for i = 1,nx-2 do begin

	  if (uu(i,j) ge 0) then begin

	    if ((i eq 1) or (i eq nx-2)) then begin

	      phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i,j)-phi(i-1,j))/(dx))

	    endif else begin

	      phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i-1,j))/(2.*dx) -(phi(i+1,j)-3.*phi(i,j)+3.*phi(i-1,j)-phi(i-2,j))/(6.*dx))

	    endelse

	  endif else begin

	    if ((i eq 1) or (i eq nx-2)) then begin

	      phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i,j))/(dx))

	    endif else begin

	      phiu(i,j) = phi(i,j) - dtt*uu(i,j)*((phi(i+1,j)-phi(i-1,j))/(2*dx) -(phi(i+2,j)-3.*phi(i+1,j)+3.*phi(i,j)-phi(i-1,j))/(6.*dx))

	    endelse

	  endelse

	endfor
	endfor

	phiu(0,0:my-1) = phiu(1,0:my-1)
        phiu(nx-1,0:my-1) = phiu(nx-2,0:my-1)

	phiuv = dblarr(nx,my)


	for i = 0,nx-1 do begin
	for j = 1,my-2 do begin

	  if (vv(i,j) ge 0) then begin

	    if ((j eq 1) or (j eq my-2)) then begin

	     phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j)-phiu(i,j-1))/(dy))

	    endif else begin

	      phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j-1))/(2*dy) - (phiu(i,j+1)-3.*phiu(i,j)+3.*phiu(i,j-1)-phiu(i,j-2))/(6.*dy))

	    endelse

	  endif else begin

	    if ((j eq 1) or (j eq my-2)) then begin

	      phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j))/(dy))

	    endif else begin

	      phiuv(i,j) = phiu(i,j) - dtt*vv(i,j)*((phiu(i,j+1)-phiu(i,j-1))/(2.*dy) - (phiu(i,j+2)-3.*phiu(i,j+1)+3.*phiu(i,j)-phiu(i,j-1))/(6.*dy)) 

	    endelse

	  endelse

	endfor
	endfor

	phiuv(0:nx-1,0) = phiuv(0:nx-1,1)
        phiuv(0:nx-1,my-1) = phiuv(0:nx-1,my-2)

        phiuv(0,0) = (phiuv(1,0) + phiuv(0,1))/2.
        phiuv(nx-1,0) = (phiuv(nx-2,0) + phiuv(nx-1,1))/2.
        phiuv(0,my-1) = (phiuv(0,my-2) + phiuv(1,my-1))/2.
        phiuv(nx-1,my-1) = (phiuv(nx-2,my-1) + phiuv(nx-1,my-2))/2.

        phi = phiuv

	return,phi

	end
