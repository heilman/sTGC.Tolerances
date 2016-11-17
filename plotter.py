#!/usr/bin/env python

import ROOT
from ROOT import *

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetLineScalePS(0.1)

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
cnv.SaveAs("board.offset.pdf")
cnv.SaveAs("board.offset.png")
cnv.SaveAs("board.pdf(")