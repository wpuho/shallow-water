	function calcuhM,hM,hmax,sigmax,sigmay,mxc,myc
	  COMMON grid,nx,my,dx,dy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	lon = dindgen(nx)
	lat = dindgen(my)

	for i = 0,nx-1 do begin

	  lon(i) = double(i+1)

	endfor

	for j = 0,my-1 do begin

	  lat(j) = double(j+1)

	endfor

	mlon = lon(mxc)
	mlat = lat(myc)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	for i = 0,nx-1 do begin
	for j = 0,my-1 do begin

	  hM(i,j) = hmax*exp(-(double((i+1)-mlon)*dx)^2/(sigmax^2)-(double((j+1)-mlat)*dy)^2/(sigmay^2))

	endfor
	endfor
	;contour,hM,nlevels=10,/isotropic
	;c_colors = [0,100,255]

	return,hM
	end
