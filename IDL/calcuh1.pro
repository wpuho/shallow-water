	function calcuh1,h1,UU,f,g
	  COMMON grid,nx,my,dx,dy

	for i = 0,nx-1 do begin

	  h1(i,0) = 0.

	for j = 0,my-1 do begin

	  case j of

	    0: h1(i,j+1) = h1(i,j) - f*UU(i,j)*dy/g

	    my-1 : h1(i,my-1) = h1(i,my-2) - f*UU(i,my-2)*dy/g

	    else : h1(i,j+1) = h1(i,j-1) - 2*f*UU(i,j)*dy/g

	  endcase

	endfor
	endfor

	return,h1

	end
