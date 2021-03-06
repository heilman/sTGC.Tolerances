#!/usr/bin/python

import math
import ROOT
from ROOT import *
from multiplet import Multiplet
import copy
import random
import string

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetLineScalePS(0.1)
# ROOT.gStyle.SetOptStat(0)

def rndstr():
   return ''.join(random.SystemRandom().choice(string.ascii_letters + string.ascii_uppercase + string.digits) for _ in range(10))

### make the quad
q0 = Multiplet("q0")
print ""


phi   = 10
theta = -98
psi   = 270
iret  = -1

def getH(q,name):
   h = q.slayers[0].histos[name].Clone()
   h.Reset()
   for i in xrange(len(q.zsboards)):
      h.Add(q.slayers[i].histos[name])
   return h

def plot(q,name,trans,text=[]):
   cname = name
   cname = cname.replace(".pdf","")
   cname = cname.replace("(","")
   cname = cname.replace(")","")
   cname = cname = cname+"."+trans
   # rootname = cname+".root"
   rootname = cname+".pdf"
   pdfonlyname = cname.replace(".","_")+"_only.pdf"
   pdfhonlyname = cname.replace(".","_")+"_honly.pdf"
   cnv = TCanvas(cname,"",1000,500)
   cnv.Divide(2,1)
   p1 = cnv.cd(1)
   p2 = cnv.cd(2)
   p1.SetPad(0.00, 0.00, 0.40, 1.00)
   p2.SetPad(0.40, 0.00, 1.00, 1.00)
   p1.cd()
   ROOT.gStyle.SetLineScalePS(0.1)
   # frame = p1.DrawFrame(0, 0, q.sboards[0].dimxtop, q.sboards[0].dimy)
   # frame.SetTitle("QS1 Strip board with offset")
   # p1.SetTicks(1,1)
   view = TView.CreateView(1)
   view.SetRange(0,0,0,q.sboards[0].dimxtop,q.sboards[0].dimy,100)
   # view.SetRange(0,0,0,q.sboards[0].dimy,q.sboards[0].dimy,q.sboards[0].dimy)
   view.SetView(phi,theta,psi,ROOT.Long(iret))

   x = TPolyLine3D(3)
   x.SetPoint(0,0,0,0)
   x.SetPoint(1,q.sboards[0].dimxtop*1.2,0,0)
   x.SetPoint(2,0,0,0)
   x.SetLineColor(ROOT.kBlack)
   x.Draw("same")
   xtxt = TPaveText(0.135,0.10,0.185,0.15,"NDC")
   xtxt.SetFillStyle(4000) # will be transparent
   xtxt.SetFillColor(0)
   xtxt.SetTextFont(42)
   xtxt.SetBorderSize(0)
   xtxt.AddText("x")
   xtxt.Draw("same")

   y = TPolyLine3D(3)
   y.SetPoint(0,0,0,0)
   y.SetPoint(1,0,q.sboards[0].dimy*1.1,0)
   y.SetPoint(2,0,0,0)
   y.SetLineColor(ROOT.kBlack)
   y.Draw("same")
   ytxt = TPaveText(0.22,0.89,0.27,0.94,"NDC")
   ytxt.SetFillStyle(4000) # will be transparent
   ytxt.SetFillColor(0)
   ytxt.SetTextFont(42)
   ytxt.SetBorderSize(0)
   ytxt.AddText("y")
   ytxt.Draw("same")

   z = TPolyLine3D(3)
   z.SetPoint(0,0,0,0)
   z.SetPoint(1,0,0,q.sboards[3].z*1.2)
   z.SetPoint(2,0,0,0)
   z.SetLineColor(ROOT.kBlack)
   z.Draw("same")
   ztxt = TPaveText(0.49,0.24,0.54,0.29,"NDC")
   ztxt.SetFillStyle(4000) # will be transparent
   ztxt.SetFillColor(0)
   ztxt.SetTextFont(42)
   ztxt.SetBorderSize(0)
   ztxt.AddText("z")
   ztxt.Draw("same")

   p1.Update()
   B = []
   S = []
   T = []
   for i in xrange(len(q.zsboards)):
      B.append(q.sboards[i].board.Clone(rndstr()))
      B[i].Draw("same")
      S.append([])
      T.append([])
      for s in xrange(len(q.slayers[i].strips)):
         S[i].append(q.slayers[i].strips[s].Clone(rndstr()))
         S[i][s].SetLineColor(S[i][s].GetLineColor())
         S[i][s].Draw("same")
         T[i].append(q.slayers[i].tstrips[trans][s].Clone(rndstr()))
         T[i][s].SetLineColor(ROOT.kGreen)
         T[i][s].Draw("same")

   p1.RedrawAxis()
   p1.Draw()
   p1.SaveAs( pdfonlyname.replace(")","").replace("(","") )

   ########
   p2.cd()
   # ROOT.gStyle.SetLineScalePS(1)
   p2.SetTicks(1,1)
   h = getH(q,trans)
   h.Draw()
   p2.RedrawAxis()
   ptxt = TPaveText(0.13,0.7,0.4,0.9,"NDC")
   ptxt.SetFillStyle(4000) # will be transparent
   ptxt.SetFillColor(0)
   ptxt.SetTextFont(42)
   ptxt.SetBorderSize(0)
   for txt in text: ptxt.AddText(txt)
   ptxt.Draw("same")
   # p2.SaveAs( pdfhonlyname.replace(")","").replace("(","") )

   cnvh = TCanvas(cname+"h","",500,500)
   h.Draw()
   ptxt.Draw("same")
   cnvh.RedrawAxis()
   cnvh.Update()
   cnvh.SaveAs( pdfhonlyname.replace(")","").replace("(","") );
   

   cnv.Update()
   cnv.SaveAs(rootname)
   cnv.SaveAs(name)



##### RMS calculation
# zweight = 0.5
# RMSlayer = {}
# RMSquadr2 = {}
# RMStotlayer2 = 0
# for i in xrange(len(q0.zsboards)):
#    for trans,sumd2 in q0.slayers[i].RMS.iteritems():
#       Npoints = q0.slayers[i].Npoints[trans]
#       rms = math.sqrt(sumd2*1.e6/Npoints) ## converted from mm to um
#       print "%s: RMS[%i]=%g" % (trans,i,rms)
#       if(trans=="YZrotationZ" or trans=="Zshift"): rms = zweight*rms
#       if(i==0):
#          RMSlayer.update({trans:rms})
#          RMSquadr2.update({trans:rms*rms})
#          RMStotlayer2 += rms*rms
#       else: RMSquadr2[trans] += rms*rms
# RMStotlayer = math.sqrt(RMStotlayer2)
# 
# RMStotquad2 = 0
# RMSquadr = {}
# for trans,sumrms2 in RMSquadr2.iteritems():
#    RMStotquad2 += sumrms2
#    RMSquadr.update({trans:math.sqrt(sumrms2)})
# RMStotquad = math.sqrt(RMStotquad2)

transformations = ["Translation", "TranslationCorr", "XYrotation", "Pitchscale", "Parallelism", "Bowing", "Zshift", "YZrotationY", "YZrotationZ"]
RMSy = {}
RMSz = {}
for trans in transformations:
   if(not (trans=="Zshift" or trans=="YZrotationY" or trans=="YZrotationZ")):
      for i in xrange(len(q0.zsboards)):
         sumd2   = q0.slayers[i].RMS[trans]
         Npoints = q0.slayers[i].Npoints[trans]
         # print "%s[%d]: sumd2=%gum^2, Npoints=%g" % (trans,i,sumd2,Npoints)
         if(i==0): RMSy.update({trans:{"sumd2":sumd2,"Npoints":Npoints}})
         else:
            RMSy[trans]["sumd2"]   += sumd2
            RMSy[trans]["Npoints"] += Npoints
   else:
      for i in xrange(len(q0.zsboards)):
         sumd2   = q0.slayers[i].RMS[trans]
         Npoints = q0.slayers[i].Npoints[trans]
         # print "%s[%d]: sumd2=%gum^2, Npoints=%g" % (trans,i,sumd2,Npoints)
         if(i==0): RMSz.update({trans:{"sumd2":sumd2,"Npoints":Npoints}})
         else:
            RMSz[trans]["sumd2"]   += sumd2
            RMSz[trans]["Npoints"] += Npoints

print ""
values = {}
for trans in transformations:
   if(not (trans=="Zshift" or trans=="YZrotationY" or trans=="YZrotationZ")):
      rms = math.sqrt(RMSy[trans]["sumd2"]*1.e6/RMSy[trans]["Npoints"])
      values.update({trans:rms})
      print "%s: RMS=%g (sumd2=%gum^2, Npoints=%g)" % (trans,rms,RMSy[trans]["sumd2"]*1.e6,RMSy[trans]["Npoints"])
   else:
      rms = math.sqrt(RMSz[trans]["sumd2"]*1.e6/RMSz[trans]["Npoints"])
      values.update({trans:rms})
      print "%s: RMS=%g (sumd2=%gum^2, Npoints=%g)" % (trans,rms,RMSz[trans]["sumd2"]*1.e6,RMSz[trans]["Npoints"])

totRMSy = 0
totRMSz = 0
totRMS  = 0
totRMSyNoOffset = 0
totRMSzNoOffset = 0
totRMSNoOffset  = 0
for trans,rms in values.iteritems():
   if(not (trans=="Zshift" or trans=="YZrotationY" or trans=="YZrotationZ")):
      if(not (trans=="TranslationCorr")): totRMSy += rms*rms
      if(not (trans=="Translation")):     totRMSyNoOffset += rms*rms
   else:
      if(not (trans=="TranslationCorr")): totRMSz += rms*rms
      if(not (trans=="Translation")):     totRMSzNoOffset += rms*rms
totRMSy = math.sqrt(totRMSy)
totRMSz = math.sqrt(totRMSz)
totRMS  = math.sqrt(totRMSy*totRMSy + (0.5*totRMSz)*(0.5*totRMSz))
totRMSyNoOffset = math.sqrt(totRMSyNoOffset)
totRMSzNoOffset = math.sqrt(totRMSzNoOffset)
totRMSNoOffset  = math.sqrt(totRMSyNoOffset*totRMSyNoOffset + (0.5*totRMSzNoOffset)*(0.5*totRMSzNoOffset))

print ""
plot(q0,"quadruplet.pdf(","Zshift",["Shift in Z","50#mum per layer"])
plot(q0,"quadruplet.pdf", "YZrotationZ",["Rotation in YZ","#it{#theta}=#frac{#it{#pi}}{1.8#times10^{4}}=0.01#circ"])
plot(q0,"quadruplet.pdf", "YZrotationY",["Rotation in YZ","#it{#theta}=#frac{#it{#pi}}{1.8#times10^{4}}=0.01#circ"])
plot(q0,"quadruplet.pdf", "Translation",["Offset in Y","75#mum per strip"])
plot(q0,"quadruplet.pdf", "TranslationCorr",["Offset correction in Y","Gaus(0,25#mum) per strip"])
plot(q0,"quadruplet.pdf", "XYrotation", ["Rotation in XY","#it{#theta}=#frac{#it{#pi}}{4.5#times10^{4}}=0.004#circ"])
plot(q0,"quadruplet.pdf", "Pitchscale", ["Pitch scale in Y","75#mum per layer"])
plot(q0,"quadruplet.pdf", "Parallelism",["Parallelism in XY","75#mum per layer"])
plot(q0,"quadruplet.pdf)","Bowing",     ["Bowing in XY","75#mum per layer"])


# print "RMS components per layer:",RMSlayer
# print "RMS total for per layer :",RMStotlayer
# print "RMS components per quad :",RMSquadr
# print "RMS total per quad      :",RMStotquad
print "\n-------------------SUMMARY-----------------------"
print "per layer:",values
print "-------------------------------------------------"
print "totRMSy =",totRMSy
print "totRMSz =",totRMSz
print "totRMS  =",totRMS
print "-------------------------------------------------"
print "totRMSyNoOffset =",totRMSyNoOffset
print "totRMSzNoOffset =",totRMSzNoOffset
print "totRMSNoOffset  =",totRMSNoOffset