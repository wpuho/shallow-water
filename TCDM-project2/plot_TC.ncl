begin
  fn = "SWM_TC.txt"

  row = numAsciiRow(fn)-1
  indata = asciiread(fn,(/row,15/),"float")
  xH = indata(:,4)
  yH = indata(:,5)
  mH = indata(:,6)
  xSPD = indata(:,7)
  ySPD = indata(:,8)
  mSPD = indata(:,9)
  xVOR = indata(:,10)
  yVOR = indata(:,11)
  mVOR = indata(:,12)
  UA   = indata(:,13)
  VA   = indata(:,14)

  fZ = "SWM_input_Z.bin"
  NX = 401
  NY = 401
  Z  = fbindirread(fZ,0,(/NY,NX/),"float")
  Z (0,0) = 0.1
  type = "x11"
  type@wkHeight = 1024
  type@wkWidth  = type@wkHeight

  wks = gsn_open_wks(type,"asda")
   
  res = True
  res@gsnFrame = False
  res@gsnDraw  = False
  zres  = res
  res@xyLineThicknessF = 5
  xyres = res
  xyres@tiXAxisString = "hours"

  xyres1 = xyres 
  xyres2 = xyres 
  xyres@tiYAxisString      = "Depth"
  xyres@tiYAxisFontColor   = "red"
  xyres@xyLineColor        = xyres@tiYAxisFontColor
  xyres@tmYLLabelFontColor = xyres@tiYAxisFontColor
  xyres1@tiYAxisString      = "SPD"
  xyres1@tiYAxisFontColor   = "darkgreen"
  xyres1@xyLineColor        = xyres1@tiYAxisFontColor
  xyres1@tmYRLabelFontColor = xyres1@tiYAxisFontColor
  xyres2@tiYAxisString      = "VOR"
  xyres2@tiYAxisFontColor   = "blue"
  xyres2@xyLineColor        = xyres2@tiYAxisFontColor
  xyres2@tmYRLabelFontColor = xyres2@tiYAxisFontColor

  zres@vpXF = 0.20
  zres@vpYF = 0.950
  zres@vpHeightF = 0.4
  zres@vpWidthF  = zres@vpHeightF 
  zres@cnFillOn = True
  zres@cnFillPalette = "MPL_Greys"
  zres@lbOrientation = "vertical"
  zres@cnLinesOn     = False
  zres@cnLineLabelsOn = False
  zres@cnInfoLabelOn  = False

  xyres@vpXF = zres@vpXF 
  xyres@vpYF = zres@vpYF - zres@vpHeightF - 0.05
  xyres@vpHeightF = zres@vpHeightF
  xyres@vpWidthF  = zres@vpHeightF
  xyres@trXMaxF   = row -1

  gsres = True
  gsres@gsLineThicknessF = 5.
  gsres@gsMarkerIndex    = 16

  plotz = gsn_csm_contour(wks,Z,zres)
  gsres@gsLineColor      = "red"
  gsres@gsMarkerColor    = gsres@gsLineColor
  gsres@gsLineColor      = "Blue"
  gsres@gsMarkerColor    = gsres@gsLineColor
  plotz@trackV = gsn_add_polyline(wks,plotz,xVOR,yVOR,gsres)
  plotz@trackVd = gsn_add_polymarker(wks,plotz,xVOR(0::24),yVOR(0::24),gsres)
  draw(plotz)

;  wmsetp("wdf", 1)
;  wmsetp("blw",0.025) ; thickness
  wmsetp("vcc",30)    ; Draw in blue.
  wmsetp("vcw",5.)   ; Increase the linewidth
  wmsetp("vrn",0.03)   ; length                    
  wmsetp("vch",0.01)   ; head size 

  do i=0,row-1,3
  wmvect(wks,xH(i),yH(i),UA(i),VA(i))
  end do
  wmvlbl(wks,1.,0.)

  plotTC = gsn_csm_xy3(wks,ispan(0,row-1,1),mH,mSPD,mVOR,xyres,xyres1,xyres2)

  draw(plotTC)
  frame(wks)

end 
