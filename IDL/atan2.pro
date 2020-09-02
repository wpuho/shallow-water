	function atan2,x,y

	  pi = acos(-1.)

	  if (x gt 0) then val=atan(y/x)

	  if (x lt 0) && (y ge 0) then val = atan(y/x) + pi

	  if (x lt 0) && (y lt 0) then val = atan(y/x) - pi

	  if (x eq 0) && (y gt 0) then val = pi/2.

	  if (x eq 0) && (y lt 0) then val = -1.*pi/2.

	  if (x eq 0) && (y eq 0) then val = 0.

	  return,val

	end
