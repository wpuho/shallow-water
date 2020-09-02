load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; drawNDCGrid(wks) gsn_code.ncl/gsn_csm.ncl/contributed.ncl/shea_util.ncl
;=======================================================================
begin

fn  = "SWM_output_d00_07-00-00.bin"

fnZ = "SWM_input_Z.bin"
NX = 401
NY = 401

U   = fbindirread(fn,0,(/NY,NX/),"float")
V   = fbindirread(fn,1,(/NY,NX/),"float")
H   = fbindirread(fn,2,(/NY,NX/),"float")
SPD = fbindirread(fn,3,(/NY,NX/),"float")
VOR = fbindirread(fn,4,(/NY,NX/),"float") ; vorticity
STR = fbindirread(fn,5,(/NY,NX/),"float") ; streamfunction
PV  = fbindirread(fn,6,(/NY,NX/),"float") ; potential vorticity
WT  = fbindirread(fn,7,(/NY,NX/),"float") ; tangential wind
WR  = fbindirread(fn,8,(/NY,NX/),"float") ; radial wind


Z      = fbindirread(fnZ,0,(/NY,NX/),"float")
;-----------------------------------------------------------------------
  wks_type          = "x11"
  wks_type@wkWidth  = 1000
  wks_type@wkHeight = 1000


  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  zres = res 
  res@cnFillOn = True
  wres = res

  wres@cnLevelSelectionMode  = "ManualLevels"
  wres@gsnScalarContour = True
  wres@vcRefAnnoOn      = False
  wres@vcGlyphStyle     = "CurlyVector"
  wres@vcMinDistanceF   = 0.025
  wres@vcRefMagnitudeF  = 10.
  wres@vcRefLengthF     = 0.025
  wres@cnMaxLevelValF        = 70.
  wres@cnMinLevelValF        =  2.
  wres@cnLevelSpacingF       =  2.
 

  zres@cnLineColor          = "darkgreen"
  zres@cnLineThicknessF     = 3.
  zres@cnLineLabelsOn       = False
  zres@cnInfoLabelOn        = False
  zres@cnLevelSelectionMode = "ExplicitLevels"
  zres@cnLevels             = (/1.,1000,2000,3000/)
  zres@cnLineDrawOrder      = "PostDraw"

  out  = "hgfds"
  wks  = gsn_open_wks(wks_type,out)
  plot = new((/7/),"graphic")
  wres@gsnLeftString = "WindSpeed"
  plot(0) = gsn_csm_vector_scalar(wks,U,V,SPD,wres)
  plotz0  = gsn_csm_contour(wks,Z,zres)
  overlay(plot(0),plotz0)

  print("H")
  res@gsnLeftString = "Depth"
  plot(1) = gsn_csm_contour(wks,H,res)

  print("VOR")
  res@gsnLeftString = "VOR"
  plot(2) = gsn_csm_contour(wks,VOR,res)

  print("STR")
  res@gsnLeftString = "STR"
  plot(3) = gsn_csm_contour(wks,STR,res)

  print("PV")
  res@gsnLeftString = "PV"
  plot(4) = gsn_csm_contour(wks,PV,res)

  print("WT")
  res@gsnLeftString = "WT"
  plot(5) = gsn_csm_contour(wks,WT,res)

  print("WR")
  res@gsnLeftString = "WR"
  plot(6) = gsn_csm_contour(wks,WR,res)

  gsn_panel(wks,plot,(/3,3/),False)
  
end
