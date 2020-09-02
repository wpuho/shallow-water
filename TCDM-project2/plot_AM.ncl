begin
  fn = "SWM_AM_WT.txt"

  row = numAsciiRow(fn)
  col = numAsciiCol(fn)
 
  dr  = 5.
  XA = ispan(0,col-1,1)*dr

  indata = asciiread(fn,(/row,col/),"float")
  indata!0 = "Time"
  indata!1 = "R"
  indata&R = XA 

  mind = new(row,integer)
  do i=0,row-1
    mind(i) = dim_maxind(indata(i,:),0)
  end do

  wks = gsn_open_wks("x11","jhgfd")

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@cnFillOn = True
  res@trXMaxF  = 500.
  res@tiXAxisString = "Radius [km]"
  res@tiYAxisString = "Time [hour]"
  res@tmYLMode = "Manual"
  res@tmYLTickSpacingF = 12.

  gsres = True
  gsres@gsLineThicknessF = 5
  gsres@gsLineColor      = "white"

  plt = gsn_csm_contour(wks,indata,res)
  aa = gsn_add_polyline(wks,plt,mind*dr,ispan(0,row,1),gsres)
  draw(plt)
  frame(wks)
end
