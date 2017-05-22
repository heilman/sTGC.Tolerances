#!/usr/bin/env python
import math
import ROOT
from ROOT import *
from array import *
import numpy as np
from scipy import linalg

ROOT.gROOT.SetBatch(1)
# ROOT.gStyle.SetLineScalePS(0.1)
ROOT.gStyle.SetOptFit(ROOT.kTRUE);

'''
To compute the pointing resolution of a segment:
make a straight-line fit to (in the sTGC case) 8 points where
x(1)=ATLAS-Z position of first sTGC plane
...
x(8)=ATLAS-Z position of eighth sTGC plane
y(1)=y(2)=...=y(8)=0
dy(1)=dy(2)=...=dy(8)=res
where res is the single-hit resolution,
or if you want to take into account chamber non-conformities,
the sum in quadrature of the single-hit resolution
and the figure of merit from your study of non-conformities,
i.e. sqrt(mean^2+rms^2) of your final histogram.
Then, if you fit y=f(x)=p0+p1*x to these eight points,
the error of p1 is your pointing angle resolution.

in each quad it goes pad strip strip pad, pad strip strip pad.
And the gap is 2.8 mm, so the strips are 1.4 mm in the appropriate
direction +Z or -Z from the gap centres
'''

gGap = 2.8     # [mm]
dz   = gGap/2. # [mm]
zConfirm = [6993.5+dz, 7004.5-dz, 7015.5+dz, 7026.5-dz] # [mm]
zPivot   = [7327.5+dz, 7338.5-dz, 7349.5+dz, 7360.5-dz] # [mm]
z    = zConfirm+zPivot
y    = [0] * len(z)
yRes = [0.043] * len(z) # the non conformities in y [mm]
# zRes = [0.028] * len(z) # the non conformities in z [mm]
zRes = [0.082] * len(z) # the non conformities in z [mm]

zArr = array('f',z)
yArr = array('f',y)
zResArr = array('f',zRes)
yResArr = array('f',yRes)
points = { "z":zArr, "dz":zResArr, "y":yArr, "dy":yResArr}
# print points

graph = TGraphErrors(len(z),points["z"],points["y"],points["dz"],points["dy"])
graph.SetTitle("")
graph.GetXaxis().SetTitle("z [mm]")
graph.GetYaxis().SetTitle("y [mm]")
graph.SetMarkerColor(ROOT.kBlack);
graph.SetMarkerStyle(20);
graph.SetMarkerSize(0.3);
graph.Fit("pol1");
f1 = graph.GetFunction("pol1") # Access the fit resuts
f1.SetLineWidth(1)

cnv = TCanvas("cnv","",500,500)
graph.Draw("AP")
graph.SetMinimum(-1)
graph.SetMaximum(+1)
cnv.SaveAs("resolution.pdf")