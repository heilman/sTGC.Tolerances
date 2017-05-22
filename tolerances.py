#!/usr/bin/env python

import math
import ROOT
from ROOT import *
from array import array
import numpy as np
from scipy import linalg

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetLineScalePS(0.1)

dimxtop = 743
dimxbot = 346.9
dimy    = 1325.6
scale   = 1.0
xcenter = dimxtop/2.
ycenter = dimy/2.
icenteralstrip = -1
# wstrip  = 2.7
# nstrips = 40
nstrips = 404
nstriptestpoints = 100

x0 = ((scale*dimxtop)-dimxbot)/2.
y0 = ((scale*dimy)-dimy)/2.
x1 = x0-(dimxtop-dimxbot)/2.
y1 = y0+dimy
x2 = x0+dimxbot+(dimxtop-dimxbot)/2.
y2 = y0+dimy
x3 = x0+dimxbot
y3 = y0
ystripmin = y0+15 ### 7.5*2 ??
ystripmax = ystripmin+1292.8
dy = (ystripmax-ystripmin)/nstrips
print "Strips dy =",dy
ncorners = 4
xcorners = array( 'f', [x0,x1,x2,x3] )
ycorners = array( 'f', [y0,y1,y2,y3] )
corners = TPolyLine(ncorners, xcorners, ycorners, "F")

def edgex(y,side):
   ### slope = (y1-y0)/(x1-x0)
   ### y-y1 = slope*(x-x1)
   ### y = slope*x + intercept
   if(side=="left"):
      slope = (y1-y0)/(x1-x0)
      intercept = -slope*x1+y1
   elif(side=="right"):
      slope = (y3-y2)/(x3-x2)
      intercept = -slope*x3+y3
   elif(side=="leftNP"):
      slope = (y1-y0)/(x1-x0)
      intercept = -slope*x1+y1
   elif(side=="rightNP"):
      slope = (y3np-y2np)/(x3np-x2np)
      intercept = -slope*x3np+y3np
   else:
      print "wrong choice:",side
      quit()
   x = (y-intercept)/slope
   return x

def ylinear(x,xa,ya,xb,yb):
   ### slope = (ya-yb)/(xa-xb)
   ### y-ya = slope*(x-xa)
   ### y = slope*x + intercept
   slope = (ya-yb)/(xa-xb)
   intercept = -slope*xa+ya
   y = slope*x + intercept
   return y

def yparabole(x,xL,yL,xC,yC,xR,yR):
   X = np.array([[xL*xL,xL,1],
                 [xC*xC,xC,1],
                 [xR*xR,xR,1]])
   Xinv = linalg.inv(X)
   Y = np.array([yL,yC,yR])
   A = Xinv.dot(Y)
   # d = (xL-xC)*(xL-xR)*(xC-xR)
   # a = (xR*(yC-yL)+xC*(yL-yR)+xL*(yR-yC))/d
   # b = (xR*xR*(yL-yC)+xC*xC*(yR-yL)+xL*xL*(yC-yR))/d
   # c = (xC*xR*(xC-xR)*yL+xR*xL*(xR-xL)*yC+xL*xC*(xL-xC)*yR)/d
   y = A[0]*x*x + A[1]*x + A[2]
   return y

def translation(y,dy):
   y += dy
   xL  = edgex(y,"left")
   xR  = edgex(y,"right")
   return (xL,xR,y)

def rotatepoint(y,axis,theta,side):
   x = edgex(y,side)-axis[0]
   y = y-axis[1]
   xr = x*math.cos(theta)-y*math.sin(theta)
   yr = x*math.sin(theta)+y*math.cos(theta)
   return (xr+axis[0],yr+axis[1])

def extend(xL,yL,xR,yR):
   slopeS = (yR-yL)/(xR-xL)
   interceptS = -slopeS*xR+yR
   slopeL = (y1-y0)/(x1-x0)
   interceptL = -slopeL*x1+y1
   slopeR = (y3-y2)/(x3-x2)
   interceptR = -slopeR*x3+y3
   xR = (interceptR-interceptS)/(slopeS-slopeR)
   yR = slopeR*xR+interceptR
   xL = (interceptL-interceptS)/(slopeS-slopeL)
   yL = slopeL*xL+interceptL
   return (xL,yL,xR,yR)

def pitchscale(y,istrip,pitch):
   yP = y + pitch*istrip
   xL  = edgex(yP,"left")
   xR  = edgex(yP,"right")
   return (xL,xR,yP)

def nonparallelscale(y,istrip,np):
   yL = y
   xL = edgex(yL,"left")
   yR = y + istrip*np
   xR = edgex(yR,"right")
   return (xL,yL,xR,yR)

def bowing(y,nxpoints,istrip,bowpar):
   yL = y
   xL = edgex(yL,"left")
   yR = y
   xR = edgex(yR,"right")
   xC = xL+(xR-xL)/2.
   yC = y + istrip*bowpar
   dxS = (xR-xL)/nxpoints
   xyvec = []
   for j in xrange(nxpoints):
      xj = xL+j*dxS
      yj = yparabole(xj,xL,yL,xC,yC,xR,yR)
      xyvec.append([xj,yj])
   xyvec.append([xR,yR])
   return xyvec

### find central strip y coordinate
y = ystripmin
xorigin = xcenter
yorigin = y
dycenter = 1.e20
for i in xrange(nstrips):
   if(abs(y-ycenter)<dycenter):
      dycenter       = abs(y-ycenter)
      yorigin        = y
      icenteralstrip = i
   y += dy


##################################################


### for the offset translation
offset = 30/1000.
# offset = 10000/1000.

### for the rotation
extendstrips = True
raxis = [180,100]   ### arbitrary chosen around the bottom left corner of the board 
theta = math.pi/45000 ### arbitrary chosen pi/500 of rotation
# theta = math.pi/60 ### arbitrary chosen pi/500 of rotation
piratio = "45k"
# piratio = "60"
pivot = TEllipse(raxis[0],raxis[1],2)

### for the pitch scaling (+-75um)
pitch = (75./nstrips)/1000.
# pitch = (20000/nstrips)/1000.
spitch = "75#mum"
# spitch = "20000#mum"

### for the non-parallelism (+-50um)
nonparallelism = (50./nstrips)/1000. ### (NP[um]/Nstrips)/um2mm
# nonparallelism = (30000./nstrips)/1000. ### (NP[um]/Nstrips)/um2mm
snonparallelism = "50#mum"
# snonparallelism = "30000#mum"

### for the bowing (+-50um)
nxpoints = 500
bowpar = (50./nstrips)/1000. ### (NP[um]/Nstrips)/um2mm
# bowpar = (20000./nstrips)/1000. ### (NP[um]/Nstrips)/um2mm

### add strips to the nominal pattern and the actual (distorted) patterns
strips      = []
translated  = []
rotated     = []
pitched     = []
nonparallel = []
bowings     = []

### histograms to contain the differences
hOffset      = TH1D("hOffset",      ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1)
hRotation    = TH1D("hRotation",    ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1)
hPitched     = TH1D("hPitched",     ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1)
hNonParallel = TH1D("hNonParallel", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1)
hBowing      = TH1D("hBowing",      ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1)

### the RMS... 
RMS     = {"offset":0, "rotation":0, "pitchscale":0, "nonparallelism":0, "bowing":0}
Npoints = {"offset":0, "rotation":0, "pitchscale":0, "nonparallelism":0, "bowing":0}
hRMS    = {"offset":0, "rotation":0, "pitchscale":0, "nonparallelism":0, "bowing":0}
hMean   = {"offset":0, "rotation":0, "pitchscale":0, "nonparallelism":0, "bowing":0}

### do the offset and rotation transformations
y = ystripmin
for i in xrange(nstrips):
   ### nominal pattern
   strip = TLine(edgex(y,"left"),y,edgex(y,"right"),y)
   strips.append(strip)

   ### translated pattern
   tcoord = translation(y,offset)
   xL=tcoord[0]
   xR=tcoord[1]
   yT=tcoord[2]
   tstrip = TLine(xL,yT,xR,yT)
   translated.append(tstrip)

   ### rotated pattern
   rcoordR = rotatepoint(y,raxis,theta,"right")
   rcoordL = rotatepoint(y,raxis,theta,"left")
   xR = rcoordR[0]
   yR = rcoordR[1]
   xL = rcoordL[0]
   yL = rcoordL[1]
   if(extendstrips):
      rcoord = extend(xL,yL,xR,yR)
      xL = rcoord[0]
      yL = rcoord[1]
      xR = rcoord[2]
      yR = rcoord[3]
   rstrip = TLine(xL,yL,xR,yR)
   rotated.append(rstrip)

   ### compare translated pattern to nominal one (trivial - all test poitns have the same hight so one is enough...)
   yNominal = y
   yActual  = yT ### assume that the x of the test point is the same as in for the nominal point
   hOffset.Fill(yActual-yNominal)
   RMS["offset"] += (yActual-yNominal)*(yActual-yNominal)
   Npoints["offset"] += 1

   ### compare rotated pattern to nominal one - need to loop over the test points
   xminstrip = edgex(y,"left")
   xmaxstrip = edgex(y,"right")
   dx = (xmaxstrip-xminstrip)/nstriptestpoints
   for p in xrange(nstriptestpoints):
      yNominal = y
      xR = rcoordR[0]
      yR = rcoordR[1]
      xL = rcoordL[0]
      yL = rcoordL[1]
      if(extendstrips):
         rcoord = extend(xL,yL,xR,yR)
         xL = rcoord[0]
         yL = rcoord[1]
         xR = rcoord[2]
         yR = rcoord[3]
      x = xminstrip + p*dx
      yActual  = ylinear(x,xL,yL,xR,yR)
      hRotation.Fill(yActual-yNominal)
      RMS["rotation"] += (yActual-yNominal)*(yActual-yNominal)
      Npoints["rotation"] += 1

   ### propogate to the next strip
   y += dy


### do the pitch scaling transformation
yup = yorigin
ydn = yorigin-dy
for i in xrange(nstrips/2):
   ### pitch-scaled pattern (going up)
   if(i==0): pcoordup = pitchscale(yup,i,0)
   else:     pcoordup = pitchscale(yup,i,+pitch)
   pcoorddn = pitchscale(ydn,i,-pitch)
   xLup = pcoordup[0]
   xRup = pcoordup[1]
   yPup = pcoordup[2]
   xLdn = pcoorddn[0]
   xRdn = pcoorddn[1]
   yPdn = pcoorddn[2]
   pstripup = TLine(xLup,yPup,xRup,yPup)
   pitched.append(pstripup)
   pstripdn = TLine(xLdn,yPdn,xRdn,yPdn)
   pitched.append(pstripdn)
   # print "yup=%g --> yPup=%g, ydn=%g --> yPdn=%g --> dyup=%g, dydn=%g" % (yup,yPup,ydn,yPdn,(yPup-yup),(ydn-yPdn))
   ### compare pitch-scaled pattern to nominal one
   yNominal = yup
   yActual  = yPup ### assume that the x of the test point is the same as in for the nominal point
   hPitched.Fill(yActual-yNominal)
   RMS["pitchscale"] += (yActual-yNominal)*(yActual-yNominal)
   Npoints["pitchscale"] += 1
   yNominal = ydn
   yActual  = yPdn
   hPitched.Fill(yActual-yNominal)
   RMS["pitchscale"] += (yActual-yNominal)*(yActual-yNominal)
   Npoints["pitchscale"] += 1
   yup += dy
   ydn -= dy


### do the non-parallel transformation
maxavpitch = 0
sumavpitch = 0
yup = yorigin
ydn = yorigin-dy
for i in xrange(nstrips/2):
   if(i==0): npcoordup = nonparallelscale(yup,i,0)
   else:     npcoordup = nonparallelscale(yup,i,+nonparallelism)
   npcoorddn = nonparallelscale(ydn,i,-nonparallelism)
   xLup = npcoordup[0]
   yLup = npcoordup[1]
   xRup = npcoordup[2]
   yRup = npcoordup[3]
   xLdn = npcoorddn[0]
   yLdn = npcoorddn[1]
   xRdn = npcoorddn[2]
   yRdn = npcoorddn[3]
   npstripup = TLine(xLup,yLup,xRup,yRup)
   nonparallel.append(npstripup)
   npstripdn = TLine(xLdn,yLdn,xRdn,yRdn)
   nonparallel.append(npstripdn)
   # print "yup=%g --> yRup=%g, ydn=%g --> yRdn=%g --> dyup=%g, dydn=%g" % (yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
   ### compare non-parallel pattern to nominal one
   xminstrip = edgex(yup,"left")
   xmaxstrip = edgex(yup,"right")
   dx = (xmaxstrip-xminstrip)/nstriptestpoints
   avpitchup = 0
   Nup = 0
   for p in xrange(nstriptestpoints):
      yNominal = yup
      x = xminstrip + p*dx
      yActual  = ylinear(x,xLup,yLup,xRup,yRup)
      hNonParallel.Fill(yActual-yNominal)
      RMS["nonparallelism"] += (yActual-yNominal)*(yActual-yNominal)
      Npoints["nonparallelism"] += 1
      Nup += 1
      avpitchup += abs(yActual-yNominal)
   avpitchup = avpitchup/Nup
   sumavpitch += avpitchup
   if(avpitchup>maxavpitch): maxavpitch = avpitchup
   xminstrip = edgex(ydn,"left")
   xmaxstrip = edgex(ydn,"right")
   dx = (xmaxstrip-xminstrip)/nstriptestpoints
   avpitchdn = 0
   Ndn = 0
   for p in xrange(nstriptestpoints):
      yNominal = ydn
      x = xminstrip + p*dx
      yActual  = ylinear(x,xLdn,yLdn,xRdn,yRdn)
      hNonParallel.Fill(yActual-yNominal)
      RMS["nonparallelism"] += (yActual-yNominal)*(yActual-yNominal)
      Npoints["nonparallelism"] += 1
      Ndn += 1
      avpitchdn += abs(yActual-yNominal)
   avpitchup = avpitchup/Nup
   sumavpitch += avpitchdn
   if(avpitchdn>maxavpitch): maxavpitch = avpitchdn
   yup += dy
   ydn -= dy


### do the bowing transformation
yup = yorigin
ydn = yorigin-dy
for i in xrange(nstrips/2):
   if(i==0): bcoordup = bowing(yup,nxpoints,i,0.)
   else:     bcoordup = bowing(yup,nxpoints,i,+bowpar)
   bcoorddn = bowing(ydn,nxpoints,i,-bowpar)
   stripU = TPolyLine(nxpoints+1)
   stripD = TPolyLine(nxpoints+1)
   for j in xrange(nxpoints+1):
      xJup = bcoordup[j][0]
      yJup = bcoordup[j][1]
      xJdn = bcoorddn[j][0]
      yJdn = bcoorddn[j][1]
      stripU.SetPoint(j,xJup,yJup)
      stripD.SetPoint(j,xJdn,yJdn)
   bowings.append(stripU)
   bowings.append(stripD)
   # print "[%d]  yup=%g-->yRup=%g, ydn=%g-->yRdn=%g --> dyup=%g, dydn=%g" % (i,yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
   ### compare bowing pattern to nominal one
   for j in xrange(nxpoints+1):
      yJup = bcoordup[j][1]
      yJdn = bcoorddn[j][1]
      yNominal = yup
      yActual = yJup
      delta = yActual-yNominal
      hBowing.Fill(delta)
      RMS["bowing"] += delta*delta
      Npoints["bowing"] += 1
      # print "up[%d] yup=%g --> yActual=%g --> dyup=%g" % (p,yup,yActual,delta)
      yNominal = ydn
      yActual = yJdn
      delta = yActual-yNominal
      hBowing.Fill(delta)
      RMS["bowing"] += delta*delta
      Npoints["bowing"] += 1
      # print "dn[%d] ydn=%g --> yActual=%g --> dydn=%g" % (p,ydn,yActual,delta)
   ##################
   ### propogate y...
   yup += dy
   ydn -= dy






############
### Draw ###
############
cnv = TCanvas("c","",1000,500)
cnv.Divide(2,1)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p1.SetPad(0.00, 0.00, 0.40, 1.00)
p2.SetPad(0.40, 0.00, 1.00, 1.00)
p1.cd()
ROOT.gStyle.SetLineScalePS(0.1)
frame = p1.DrawFrame(0, 0, scale*dimxtop, scale*dimy)
frame.SetTitle("QS1 Strip board with offset")
p1.SetTicks(1,1)
corners.SetLineColor(ROOT.kGreen+2)
corners.SetFillColor(ROOT.kGreen+2)
corners.Draw("f")
for strip in strips:
   strip.SetLineColor(ROOT.kRed-1)
   strip.SetLineWidth(1)
   strip.Draw("same")
for tstrip in translated:
   tstrip.SetLineColor(ROOT.kRed)
   tstrip.SetLineWidth(1)
   # tstrip.SetLineStyle(2)
   tstrip.Draw("same")
origin = TEllipse(xorigin,yorigin,2)
origin.SetLineColor(ROOT.kBlack)
origin.SetFillColor(ROOT.kBlack)
origin.Draw("same")
p1.RedrawAxis()
p1.Draw()
########
p2.cd()
ROOT.gStyle.SetLineScalePS(1)
p2.SetTicks(1,1)
hOffset.Draw()
ptxt = TPaveText(0.15,0.5,0.4,0.9,"NDC");
ptxt.SetFillStyle(4000); # will be transparent
ptxt.SetFillColor(0);
ptxt.SetTextFont(42);
ptxt.SetBorderSize(0);
ptxt.AddText("Offset = "+str(offset*1000.)+"#mum");
ptxt.Draw("same")
p2.RedrawAxis()
p2.Draw()
cnv.Update()
cnv.SaveAs("board_offset.pdf")
cnv.SaveAs("board_offset.png")
cnv.SaveAs("board_pdf(")
p1.SaveAs("board_only_offset.pdf")


cnv = TCanvas("c","",1000,500)
cnv.Divide(2,1)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p1.SetPad(0.00, 0.00, 0.40, 1.00)
p2.SetPad(0.40, 0.00, 1.00, 1.00)
p1.cd()
ROOT.gStyle.SetLineScalePS(0.1)
frame = p1.DrawFrame(0, 0, scale*dimxtop, scale*dimy)
frame.SetTitle("QS1 Strip board with rotation")
p1.SetTicks(1,1)
corners.SetLineColor(ROOT.kGreen+2)
corners.SetFillColor(ROOT.kGreen+2)
corners.Draw("f")
for strip in strips:
   strip.SetLineColor(ROOT.kRed-1)
   strip.SetLineWidth(1)
   strip.Draw("same")
for rstrip in rotated:
   rstrip.SetLineColor(ROOT.kRed)
   rstrip.SetLineWidth(1)
   # rstrip.SetLineStyle(2)
   rstrip.Draw("same")
pivot.SetLineColor(ROOT.kRed)
pivot.SetFillColor(ROOT.kRed)
pivot.SetLineWidth(1)
rstrip.Draw("same")
pivot.Draw("same")
pivotlabel = TText(50,90,"Rotation axis")
pivotlabel.SetTextColor(ROOT.kRed)
pivotlabel.SetTextFont(43)
pivotlabel.SetTextSize(10)
pivotlabel.Draw("same")
origin = TEllipse(xorigin,yorigin,2)
origin.SetLineColor(ROOT.kBlack)
origin.SetFillColor(ROOT.kBlack)
origin.Draw("same")
p1.RedrawAxis()
p1.Draw()
########
p2.cd()
ROOT.gStyle.SetLineScalePS(1)
p2.SetTicks(1,1)
hRotation.Draw()
ptxt = TPaveText(0.15,0.7,0.5,0.85,"NDC");
ptxt.SetFillStyle(4000); # will be transparent
ptxt.SetFillColor(0);
ptxt.SetTextFont(42);
ptxt.SetBorderSize(0);
ptxt.AddText("Rotation in #it{#theta}=#frac{#it{#pi}}{"+piratio+"}");
ptxt.AddText("arround ("+str(raxis[0])+"mm,"+str(raxis[1])+"mm)");
ptxt.Draw("same")
p2.RedrawAxis()
p2.Draw()
cnv.Update()
cnv.SaveAs("board_rotationxy.pdf")
cnv.SaveAs("board_rotationxy.png")
cnv.SaveAs("board_pdf")
p1.SaveAs("board_only_rotationxy.pdf")


cnv = TCanvas("c","",1000,500)
cnv.Divide(2,1)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p1.SetPad(0.00, 0.00, 0.40, 1.00)
p2.SetPad(0.40, 0.00, 1.00, 1.00)
p1.cd()
ROOT.gStyle.SetLineScalePS(0.1)
frame = p1.DrawFrame(0, 0, scale*dimxtop, scale*dimy)
frame.SetTitle("QS1 Strip board with pitch scale")
p1.SetTicks(1,1)
corners.SetLineColor(ROOT.kGreen+2)
corners.SetFillColor(ROOT.kGreen+2)
corners.Draw("f")
for strip in strips:
   strip.SetLineColor(ROOT.kRed-1)
   strip.SetLineWidth(1)
   strip.Draw("same")
for pstrip in pitched:
   pstrip.SetLineColor(ROOT.kRed)
   pstrip.SetLineWidth(1)
   # pstrip.SetLineStyle(2)
   pstrip.Draw("same")
origin = TEllipse(xorigin,yorigin,2)
origin.SetLineColor(ROOT.kBlack)
origin.SetFillColor(ROOT.kBlack)
origin.Draw("same")
p1.RedrawAxis()
p1.Draw()
########
p2.cd()
ROOT.gStyle.SetLineScalePS(1)
p2.SetTicks(1,1)
hPitched.Draw()
ptxt = TPaveText(0.12,0.65,0.4,0.85,"NDC");
ptxt.SetFillStyle(4000) # will be transparent
ptxt.SetFillColor(0)
ptxt.SetTextFont(42)
ptxt.SetBorderSize(0)
ptxt.AddText("Pitch scale = "+spitch)
ptxt.AddText("(across the board)")
ptxt.AddText("(origin at central strip)")
ptxt.Draw("same")
p2.RedrawAxis()
p2.Draw()
cnv.Update()
cnv.SaveAs("board_pitchscale.pdf")
cnv.SaveAs("board_pitchscale.png")
cnv.SaveAs("board_pdf")
p1.SaveAs("board_only_pitchscale.pdf")


cnv = TCanvas("c","",1000,500)
cnv.Divide(2,1)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p1.SetPad(0.00, 0.00, 0.40, 1.00)
p2.SetPad(0.40, 0.00, 1.00, 1.00)
p1.cd()
ROOT.gStyle.SetLineScalePS(0.1)
frame = p1.DrawFrame(0, 0, scale*dimxtop, scale*dimy)
frame.SetTitle("QS1 Strip board with non-parallelism")
p1.SetTicks(1,1)
corners.SetLineColor(ROOT.kGreen+2)
corners.SetFillColor(ROOT.kGreen+2)
corners.Draw("f")
for strip in strips:
   strip.SetLineColor(ROOT.kRed-1)
   strip.SetLineWidth(1)
   strip.Draw("same")
for npstrip in nonparallel:
   npstrip.SetLineColor(ROOT.kRed)
   npstrip.SetLineWidth(1)
   # npstrip.SetLineStyle(2)
   npstrip.Draw("same")
origin = TEllipse(xorigin,yorigin,2)
origin.SetLineColor(ROOT.kBlack)
origin.SetFillColor(ROOT.kBlack)
origin.Draw("same")
p1.RedrawAxis()
p1.Draw()
########
p2.cd()
ROOT.gStyle.SetLineScalePS(1)
p2.SetTicks(1,1)
hNonParallel.Draw()
ptxt = TPaveText(0.15,0.5,0.45,0.87,"NDC");
ptxt.SetFillStyle(4000) # will be transparent
ptxt.SetFillColor(0)
ptxt.SetTextFont(42)
ptxt.SetBorderSize(0)
ptxt.AddText("Non-parapllelism: "+snonparallelism)
ptxt.AddText("(across the board)")
ptxt.AddText("Changing the pitch")
ptxt.AddText("only at the right edge")
ptxt.AddText("(origin at central strip)")
ptxt.Draw("same")
p2.RedrawAxis()
p2.Draw()
cnv.Update()
cnv.SaveAs("board_parapllelism.pdf")
cnv.SaveAs("board_parapllelism.png")
cnv.SaveAs("board_pdf")
p1.SaveAs("board_only_parapllelism.pdf")


cnv = TCanvas("c","",1000,500)
cnv.Divide(2,1)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p1.SetPad(0.00, 0.00, 0.40, 1.00)
p2.SetPad(0.40, 0.00, 1.00, 1.00)
p1.cd()
ROOT.gStyle.SetLineScalePS(0.1)
frame = p1.DrawFrame(0, 0, scale*dimxtop, scale*dimy)
frame.SetTitle("QS1 Strip board with bowing")
p1.SetTicks(1,1)
corners.SetLineColor(ROOT.kGreen+2)
corners.SetFillColor(ROOT.kGreen+2)
corners.Draw("f")
for strip in strips:
   strip.SetLineColor(ROOT.kRed-1)
   strip.SetLineWidth(1)
   strip.Draw("same")
for bstrip in bowings:
   bstrip.SetLineColor(ROOT.kRed)
   bstrip.SetLineWidth(1)
   # bstrip.SetLineStyle(2)
   bstrip.Draw("same")
origin = TEllipse(xorigin,yorigin,2)
origin.SetLineColor(ROOT.kBlack)
origin.SetFillColor(ROOT.kBlack)
origin.Draw("same")
p1.RedrawAxis()
p1.Draw()
########
p2.cd()
ROOT.gStyle.SetLineScalePS(1)
p2.SetTicks(1,1)
hBowing.Draw()
ptxt = TPaveText(0.15,0.5,0.45,0.87,"NDC");
ptxt.SetFillStyle(4000) # will be transparent
ptxt.SetFillColor(0)
ptxt.SetTextFont(42)
ptxt.SetBorderSize(0)
ptxt.AddText("Bowing: "+snonparallelism)
ptxt.AddText("(across the board)")
ptxt.AddText("Maximum pitch")
ptxt.AddText("at the middle")
ptxt.AddText("(origin at central strip)")
ptxt.Draw("same")
p2.RedrawAxis()
p2.Draw()
cnv.Update()
cnv.SaveAs("board_bowing.pdf")
cnv.SaveAs("board_bowing.png")
cnv.SaveAs("board_pdf)")
p1.SaveAs("board_only_bowing.pdf")



### the RMS...
hMean["offset"] = hOffset.GetMean()*1000
hMean["rotation"] = hRotation.GetMean()*1000
hMean["pitchscale"] = hPitched.GetMean()*1000
hMean["nonparallelism"] = hNonParallel.GetMean()*1000
hMean["bowing"] = hBowing.GetMean()*1000
hRMS["offset"] = hOffset.GetRMS()*1000
hRMS["rotation"] = hRotation.GetRMS()*1000
hRMS["pitchscale"] = hPitched.GetRMS()*1000
hRMS["nonparallelism"] = hNonParallel.GetRMS()*1000
hRMS["bowing"] = hBowing.GetRMS()*1000
RMSquadrature = 0
hRMSquadrature = 0
for name,rms in RMS.iteritems():
   RMS[name] = math.sqrt(rms/Npoints[name])*1000
   RMSquadrature += RMS[name]*RMS[name]
   hRMSquadrature += hRMS[name]*hRMS[name]
RMSquadrature = math.sqrt(RMSquadrature)
hRMSquadrature = math.sqrt(hRMSquadrature)
print "RMS =",RMS
print "hRMS =",hRMS
print "hMean =",hMean
print "RMS quadrature =",RMSquadrature
print "hRMS quadrature =",hRMSquadrature
print "maxavpitch =",maxavpitch
print "sumavpitch =",sumavpitch