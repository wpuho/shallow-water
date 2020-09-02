	function sor,h,u,v,f,omega
	  COMMON grid,nx,my,dx,dy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	  g = 9.8
	  n = 0
	  diff = 100.

	  h2 = dblarr(nx,my)

	  if (omega eq 0.) then omega = 1.

	  fxy = dblarr(nx,my)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	  while ((diff gt 0.001) && (n lt 10000)) do begin

	    n = n + 1
;	    h2 = h

	    for j = 0,my-1 do begin
	    for i = 0,nx-1 do begin

	      case i of

		0 : begin

		  ux = (u(i+1,j) - u(i,j))/dx
		  vx = (v(i+1,j) - v(i,j))/dx

		    end

		nx-1 : begin

		  ux = (u(i,j) - u(i-1,j))/dx
		  vx = (v(i,j) - v(i-1,j))/dx

		    end

		else : begin

		  ux = (u(i+1,j)-u(i-1,j))/(2*dx)
		  vx = (v(i+1,j)-v(i-1,j))/(2*dx)

		    end
		    
	      endcase

	      case j of

		0 : begin

		  uy = (u(i,j+1) - u(i,j))/dy
		  vy = (v(i,j+1) - v(i,j))/dy

		    end

		my-1 : begin

		  uy = (u(i,j) - u(i,j-1))/dy
		  vy = (v(i,j) - v(i,j-1))/dy

		    end

		else : begin

		  uy = (u(i,j+1)-u(i,j-1))/(2*dy)
		  vy = (v(i,j+1)-v(i,j-1))/(2*dy)

		    end

	      endcase

	      fvo = f*(vx-uy)
	      jco = ux*vy - uy*vx
	      fxy(i,j) = (fvo + 2.*jco)/g

	    endfor
	    endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	    for j = 0,my-1 do begin
	    for i = 0,nx-1 do begin

	      bc1 = f*v(i,j)/g
	      bc2 = -1.*f*u(i,j)/g

	      if (i eq 0) then begin

		hL = h(i+1,j) - 2.*bc1*dx

	      endif else begin

		hL = h(i-1,j)

	      endelse

	      if (j eq 0) then begin

		hD = h(i,j+1) - 2.*bc2*dy

	      endif else begin

		hD = h(i,j-1)

	      endelse

	      if (i eq nx-1) then begin

		hR = h(i-1,j) + 2.*bc1*dx

	      endif else begin

		hR = h(i+1,j)

	      endelse

	      if (j eq my-1) then begin

		hU = h(i,j-1) + 2.*bc2*dy

	      endif else begin

		hU = h(i,j+1)

	      endelse

	      if (dx eq dy) then begin

		h2(i,j) = h(i,j) + omega*(hU+hD+hR+hL-4.*h(i,j)-fxy(i,j)*dx*dy)/4.

	      endif else begin

		ratio = (dx*dx)*(dy*dy)/(2.*(dx*dx)+2.*(dy*dy))
		resij = (hR+hL-2.*h(i,j))/(dx*dx) + (hU+hD-2*h(i,j))/(dy*dy) - fxy(i,j)
		h2(i,j) = h(i,j) + omega*ratio*resij

	      endelse

	    endfor
	    endfor

	    aveh = mean(abs(h2))
	    diff2 = abs((h2-h)/aveh) 
	    diff = max(diff2)

	    h = h2

	    h = h - h(0,0)

	    if ((n mod 100) eq 0) then print,'n',n,'  diff',diff

	  endwhile

	  return,h

	end

