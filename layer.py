#!/usr/bin/python

import ROOT
from ROOT import *
import math
import copy
import numpy as np
from scipy import linalg
import progressbar
from time import sleep

############################################################
############################################################
############################################################
class StripLayer:
   'Class for strip layer'

   def __init__(self,name,board,color=ROOT.kGray+1,nstrips=404,npoints=100):
      # arguments
      self.name = name
      self.board = board
      self.nstrips = nstrips
      self.npoints = npoints
      self.color = color
      self.yorigin = -1
      self.xorigin = -1
      self.icenteralstrip = -1

      self.wstrip  = 2.7
      self.ystripmin = self.board.y0+15 ### 7.5*2 ??
      self.ystripmax = self.ystripmin+1292.8
      self.dy = (self.ystripmax-self.ystripmin)/self.nstrips

      ### for the transformations
      self.tstrips = {"Translation":[],"XYrotation":[],"Pitchscale":[],"Parallelism":[],"Bowing":[],"Zshift":[],"YZrotationY":[],"YZrotationZ":[]}
      self.Npoints = {"Translation":0, "XYrotation":0, "Pitchscale":0, "Parallelism":0,"Bowing":0, "Zshift":0, "YZrotationY":0, "YZrotationZ":0}
      self.RMS     = {"Translation":0, "XYrotation":0, "Pitchscale":0, "Parallelism":0,"Bowing":0, "Zshift":0, "YZrotationY":0, "YZrotationZ":0}
      self.histos  = {"Translation": TH1D(self.name+"_Strip_Translation", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "XYrotation":  TH1D(self.name+"_Strip_XYrotation",  ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Pitchscale":  TH1D(self.name+"_Strip_Pitchscale",  ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Parallelism": TH1D(self.name+"_Strip_Parallelism", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Bowing":      TH1D(self.name+"_Strip_Bowing",      ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Zshift":      TH1D(self.name+"_Strip_Zshift",      ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "YZrotationY": TH1D(self.name+"_Strip_YZrotationY", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "YZrotationZ": TH1D(self.name+"_Strip_YZrotationZ", ";z_{actual}-z_{nominal} [mm];Number of test points",50,-0.1,+0.1)}
      self.offset = -1
      self.extendstrips = True
      self.raxis = [-1,-1]
      self.thetaxy = -1
      self.piratioxy = ""
      self.pivot = TEllipse(self.raxis[0],self.raxis[1],2)
      self.pitch = -1
      self.spitch = ""
      self.nonparallelism = -1
      self.snonparallelism = ""
      self.zshift = -1
      self.szshift = ""
      self.thetayz = -1
      self.sthetayz = ""
      self.bowpar = -1
      self.sbowpar = ""

      #################
      self.strips = []
      self.origin()
      self.pattern()
      # print "in class self.strips=",self.strips
      #################

   def __del__(self):
      class_name = self.__class__.__name__
      # print class_name, "destroyed"

   def __cmp__(self, other):
      return cmp(self.name, other.name)
   def __copy__(self):
      # print '__copy__()'
      return StripLayer(self.name)
   def __deepcopy__(self, memo):
      # print '__deepcopy__(%s)' % str(memo)
      return StripLayer(copy.deepcopy(self.name, memo))

   def origin(self):
      ### find central strip y coordinate
      y = self.ystripmin
      self.xorigin = self.board.xcenter
      self.yorigin = y
      dycenter = 1.e20
      for i in xrange(self.nstrips):
         if(abs(y-self.board.ycenter)<dycenter):
            dycenter = abs(y-self.board.ycenter)
            self.yorigin = y
            self.icenteralstrip = i
         y += self.dy

   def pattern(self):
      y = self.ystripmin
      for i in xrange(self.nstrips):
         strip = TPolyLine3D(3)
         strip.SetPoint(0,self.board.edgex(y,"left"),y,self.board.z)
         strip.SetPoint(1,self.board.edgex(y,"right"),y,self.board.z)
         strip.SetLineColor(self.color)
         self.strips.append(strip)
         y += self.dy

   def ylinear(self,x,xa,ya,xb,yb):
      ### slope = (ya-yb)/(xa-xb)
      ### y-ya = slope*(x-xa)
      ### y = slope*x + intercept
      slope = (ya-yb)/(xa-xb)
      intercept = -slope*xa+ya
      y = slope*x + intercept
      return y

   def yparabole(self,x,xL,yL,xC,yC,xR,yR):
      X = np.array([[xL*xL,xL,1],
                    [xC*xC,xC,1],
                    [xR*xR,xR,1]])
      Xinv = linalg.inv(X)
      Y = np.array([yL,yC,yR])
      A = Xinv.dot(Y)
      y = A[0]*x*x + A[1]*x + A[2]
      return y

   def extend(self,xL,yL,xR,yR):
      slopeS = (yR-yL)/(xR-xL)
      interceptS = -slopeS*xR+yR
      slopeL = (self.board.y1-self.board.y0)/(self.board.x1-self.board.x0)
      interceptL = -slopeL*self.board.x1+self.board.y1
      slopeR = (self.board.y3-self.board.y2)/(self.board.x3-self.board.x2)
      interceptR = -slopeR*self.board.x3+self.board.y3
      xR = (interceptR-interceptS)/(slopeS-slopeR)
      yR = slopeR*xR+interceptR
      xL = (interceptL-interceptS)/(slopeS-slopeL)
      yL = slopeL*xL+interceptL
      return (xL,yL,xR,yR)

   def translation(self,y):
      y += self.offset
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      return (xL,xR,y)

   def xyrotattion(self,y,side):
      x = self.board.edgex(y,side)-self.raxis[0]
      y = y-self.raxis[1]
      xr = x*math.cos(self.thetaxy)-y*math.sin(self.thetaxy)
      yr = x*math.sin(self.thetaxy)+y*math.cos(self.thetaxy)
      return (xr+self.raxis[0],yr+self.raxis[1])

   def pitchscale(self,y,istrip,pitchsign):
      yP = y + pitchsign*istrip*self.pitch
      xL  = self.board.edgex(yP,"left")
      xR  = self.board.edgex(yP,"right")
      return (xL,xR,yP)

   def parallelism(self,y,istrip,npsign):
      yL = y
      xL = self.board.edgex(yL,"left")
      yR = y + npsign*istrip*self.nonparallelism
      xR = self.board.edgex(yR,"right")
      return (xL,yL,xR,yR)

   def bowing(self,y,nxpoints,istrip,bsign):
      yL = y
      xL = self.board.edgex(yL,"left")
      yR = y
      xR = self.board.edgex(yR,"right")
      xC = xL+(xR-xL)/2.
      yC = y + bsign*istrip*self.bowpar
      dx = (xR-xL)/nxpoints
      xyvec = []
      for i in xrange(nxpoints):
         xi = xL+i*dx
         yi = self.yparabole(xi,xL,yL,xC,yC,xR,yR)
         xyvec.append([xi,yi])
      xyvec.append([xR,yR])
      return xyvec

   def longidudinalshift(self,y,z,shiftsign):
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      z = z + shiftsign*self.zshift
      return (xL,xR,y,z)

   def yzrotation(self,y,z):
      ### do not change X !
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      rsigny = 1.
      rsignz = 1.
      if(y>self.yorigin):
         rsigny = -1.
         rsignz = +1.
      elif(y<self.yorigin):
         rsigny = +1.
         rsignz = -1.
      else:
         rsigny = 0.
         rsignz = 0.
      # yS = self.yorigin + rsign*abs(y)*math.cos(self.thetayz)
      # zS = self.board.z + rsign*abs(y)*math.sin(self.thetayz)
      dz = abs((y-self.yorigin)*math.sin(self.thetayz))
      zS = self.board.z + rsignz*dz
      yS = y            + rsigny*abs(dz*math.tan(self.thetayz))
      return (xL,xR,yS,zS)


   def transform(self,transname,color=ROOT.kRed-2,offset=30.,raxis=[180.,100.],thetaxy=math.pi/45000.,pitchum=75.,nonparallelismum=50.,zshift=10.,thetayz=math.pi/45000.,bowpar=50.):
      print "Transforming "+transname
	
      ### for the offset translation
      self.offset = offset/1000.

      ### for the rotation
      self.extendstrips = True
      self.raxis = raxis
      self.thetaxy = thetaxy
      self.piratioxy = str(thetaxy/math.pi)
      self.pivot = TEllipse(raxis[0],raxis[1],2)

      ### for the pitch scaling (+-75um)
      self.pitch = (pitchum/self.nstrips)/1000. ## per strip
      self.spitch = "75#mum"

      ### for the non-parallelism (+-50um)
      self.nonparallelism = (nonparallelismum/self.nstrips)/1000. ## per strip
      self.snonparallelism = "50#mum"

      ### for the offset translation
      self.zshift = zshift/1000. ## per layer
      self.szshift = "10#mum"

      ### for the y-z rotation
      self.thetayz = thetayz
      self.piratioyz = str(thetayz/math.pi)

      ### for the bowing
      self.bowpar = (bowpar/self.nstrips)/1000.
      self.sbowpar = "50#mum"
      self.nxpoints = 100

      #################################
      ### Start the transformations ###
      #################################
      bar = progressbar.ProgressBar(maxval=self.nstrips,widgets=[progressbar.Bar('=', 'Translation [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      y = self.ystripmin
      for i in xrange(self.nstrips):
         bar.update(i+1)
         sleep(0.001)
         ### translated pattern
         tcoord = self.translation(y)
         xL = tcoord[0]
         xR = tcoord[1]
         yT = tcoord[2]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yT,self.board.z)
         strip.SetPoint(1,xR,yT,self.board.z)
         self.tstrips["Translation"].append(strip)
         ### compare translated pattern to nominal one
         yNominal = y
         yActual  = yT ### assume that the x of the test point is the same as in for the nominal point
         delta = yActual-yNominal
         self.histos["Translation"].Fill(delta)
         self.RMS["Translation"] += delta*delta
         self.Npoints["Translation"] += 1
         ###################
         ### propogate y...
         y += self.dy
      bar.finish()

      bar = progressbar.ProgressBar(maxval=self.nstrips,widgets=[progressbar.Bar('=', 'XY Rotation [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      y = self.ystripmin
      for i in xrange(self.nstrips):
         bar.update(i+1)
         sleep(0.001)
         ### rotated pattern
         rcoordR = self.xyrotattion(y,"right")
         rcoordL = self.xyrotattion(y,"left")
         xR = rcoordR[0]
         yR = rcoordR[1]
         xL = rcoordL[0]
         yL = rcoordL[1]
         if(self.extendstrips):
            rcoord = self.extend(xL,yL,xR,yR)
            xL = rcoord[0]
            yL = rcoord[1]
            xR = rcoord[2]
            yR = rcoord[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yL,self.board.z)
         strip.SetPoint(1,xR,yR,self.board.z)
         self.tstrips["XYrotation"].append(strip)
         ### compare rotated pattern to nominal one
         xminstrip = self.board.edgex(y,"left")
         xmaxstrip = self.board.edgex(y,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            yNominal = y
            xR = rcoordR[0]
            yR = rcoordR[1]
            xL = rcoordL[0]
            yL = rcoordL[1]
            if(self.extendstrips):
               rcoord = self.extend(xL,yL,xR,yR)
               xL = rcoord[0]
               yL = rcoord[1]
               xR = rcoord[2]
               yR = rcoord[3]
            x = xminstrip + p*dx
            yActual = self.ylinear(x,xL,yL,xR,yR)
            delta = yActual-yNominal
            self.histos["XYrotation"].Fill(delta)
            self.RMS["XYrotation"] += delta*delta
            self.Npoints["XYrotation"] += 1
         ###################
         ### propogate y...
         y += self.dy
      bar.finish()


      ### do the pitch scaling transformation
      bar = progressbar.ProgressBar(maxval=self.nstrips/2,widgets=[progressbar.Bar('=', 'Pitch scaling [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         bar.update(i+1)
         sleep(0.001)
         ### pitch-scaled pattern (going up)
         if(i==0): pcoordup = self.pitchscale(yup,i,0.)
         else:     pcoordup = self.pitchscale(yup,i,+1.)
         pcoorddn = self.pitchscale(ydn,i,-1.)
         xLup = pcoordup[0]
         xRup = pcoordup[1]
         yPup = pcoordup[2]
         xLdn = pcoorddn[0]
         xRdn = pcoorddn[1]
         yPdn = pcoorddn[2]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLup,yPup,self.board.z)
         strip.SetPoint(1,xRup,yPup,self.board.z)
         self.tstrips["Pitchscale"].append(strip)
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLdn,yPdn,self.board.z)
         strip.SetPoint(1,xRdn,yPdn,self.board.z)
         self.tstrips["Pitchscale"].append(strip)
         ### compare pitch-scaled up pattern to nominal one
         yNominal = yup
         yActual  = yPup
         delta = yActual-yNominal
         self.histos["Pitchscale"].Fill(delta)
         self.RMS["Pitchscale"] += delta*delta
         self.Npoints["Pitchscale"] += 1
         ### compare pitch-scaled down pattern to nominal one
         yNominal = ydn
         yActual  = yPdn
         delta = yActual-yNominal
         self.histos["Pitchscale"].Fill(delta)
         self.RMS["Pitchscale"] += delta*delta
         self.Npoints["Pitchscale"] += 1
         # print "yup=%g --> yPup=%g, ydn=%g --> yPdn=%g --> dyup=%g, dydn=%g" % (yup,yPup,ydn,yPdn,(yPup-yup),(ydn-yPdn))
         ###################
         ### propogate y...
         yup += self.dy
         ydn -= self.dy
      bar.finish()

      ### do the non-parallel transformation
      bar = progressbar.ProgressBar(maxval=self.nstrips/2,widgets=[progressbar.Bar('=', 'Non-parallelism [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         bar.update(i+1)
         sleep(0.001)
         if(i==0): npcoordup = self.parallelism(yup,i,0.)
         else:     npcoordup = self.parallelism(yup,i,+1.)
         npcoorddn = self.parallelism(ydn,i,-1.)
         xLup = npcoordup[0]
         yLup = npcoordup[1]
         xRup = npcoordup[2]
         yRup = npcoordup[3]
         xLdn = npcoorddn[0]
         yLdn = npcoorddn[1]
         xRdn = npcoorddn[2]
         yRdn = npcoorddn[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLup,yLup,self.board.z)
         strip.SetPoint(1,xRup,yRup,self.board.z)
         self.tstrips["Parallelism"].append(strip)
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLdn,yLdn,self.board.z)
         strip.SetPoint(1,xRdn,yRdn,self.board.z)
         self.tstrips["Parallelism"].append(strip)
         # print "[%d]  yup=%g-->yRup=%g, ydn=%g-->yRdn=%g --> dyup=%g, dydn=%g" % (i,yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
         ### compare non-parallel pattern up to nominal one
         xminstrip = self.board.edgex(yup,"left")
         xmaxstrip = self.board.edgex(yup,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            x = xminstrip + p*dx
            yNominal = yup
            yActual = self.ylinear(x,xLup,yLup,xRup,yRup)
            delta = yActual-yNominal
            self.histos["Parallelism"].Fill(delta)
            self.RMS["Parallelism"] += delta*delta
            self.Npoints["Parallelism"] += 1
            # print "up[%d] yup=%g --> yActual=%g --> dyup=%g" % (p,yup,yActual,delta)
         ### compare non-parallel pattern down to nominal one
         xminstrip = self.board.edgex(ydn,"left")
         xmaxstrip = self.board.edgex(ydn,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            x = xminstrip + p*dx
            yNominal = ydn
            yActual = self.ylinear(x,xLdn,yLdn,xRdn,yRdn)
            delta = yActual-yNominal
            self.histos["Parallelism"].Fill(delta)
            self.RMS["Parallelism"] += delta*delta
            self.Npoints["Parallelism"] += 1
            # print "dn[%d] ydn=%g --> yActual=%g --> dydn=%g" % (p,ydn,yActual,delta)
         # print "RMS[%d]=%g" % (i,self.RMS["Parallelism"])
         ##################
	     ### propogate y...
         yup += self.dy
         ydn -= self.dy
      bar.finish()

      ### do the bowing transformation
      bar = progressbar.ProgressBar(maxval=self.nstrips/2,widgets=[progressbar.Bar('=', 'Bowing [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         bar.update(i+1)
         sleep(0.001)
         if(i==0): bcoordup = self.bowing(yup,self.nxpoints,i,0.)
         else:     bcoordup = self.bowing(yup,self.nxpoints,i,+1.)
         bcoorddn = self.bowing(ydn,self.nxpoints,i,-1.)
         stripU = TPolyLine3D(self.nxpoints+1)
         stripD = TPolyLine3D(self.nxpoints+1)
         for j in xrange(self.nxpoints):
            xJup = bcoordup[j][0]
            yJup = bcoordup[j][1]
            xJdn = bcoorddn[j][0]
            yJdn = bcoorddn[j][1]
            stripU.SetPoint(j,xJup,yJup,self.board.z)
            stripD.SetPoint(j,xJdn,yJdn,self.board.z)
         self.tstrips["Bowing"].append(stripU)
         self.tstrips["Bowing"].append(stripD)
         # print "[%d]  yup=%g-->yRup=%g, ydn=%g-->yRdn=%g --> dyup=%g, dydn=%g" % (i,yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
         ### compare bowing pattern up to nominal one
         for j in xrange(self.nxpoints):
            yJup = bcoordup[j][1]
            yJdn = bcoorddn[j][1]
            yNominal = yup
            yActual = yJup
            delta = yActual-yNominal
            self.histos["Bowing"].Fill(delta)
            self.RMS["Bowing"] += delta*delta
            self.Npoints["Bowing"] += 1
            # print "up[%d] yup=%g --> yActual=%g --> dyup=%g" % (p,yup,yActual,delta)
            yNominal = ydn
            yActual = yJdn
            delta = yActual-yNominal
            self.histos["Bowing"].Fill(delta)
            self.RMS["Bowing"] += delta*delta
            self.Npoints["Bowing"] += 1
            # print "dn[%d] ydn=%g --> yActual=%g --> dydn=%g" % (p,ydn,yActual,delta)
         ##################
	     ### propogate y...
         yup += self.dy
         ydn -= self.dy
      bar.finish()


      ### do the longitudinal shift transformation
      bar = progressbar.ProgressBar(maxval=self.nstrips,widgets=[progressbar.Bar('=', 'Longitudinal shift [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      y = self.ystripmin
      for i in xrange(self.nstrips):
         bar.update(i+1)
         sleep(0.001)
         scoords = self.longidudinalshift(y,self.board.z,+1)
         xL = scoords[0]
         xR = scoords[1]
         y  = scoords[2]
         zS = scoords[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,y,zS)
         strip.SetPoint(1,xR,y,zS)
         self.tstrips["Zshift"].append(strip)
         ### compare z-shifted pattern to nominal one
         zNominal = self.board.z
         zActual  = zS
         delta = zActual-zNominal
         self.histos["Zshift"].Fill(delta)
         self.RMS["Zshift"] += delta*delta
         self.Npoints["Zshift"] += 1
         ###################
         ### propogate y...
         y += self.dy
      bar.finish()

      ### do the yz-rotatied transformation
      bar = progressbar.ProgressBar(maxval=self.nstrips,widgets=[progressbar.Bar('=', 'YZ Rotation [', ']'), ' ', progressbar.Percentage()])
      bar.start()
      y = self.ystripmin
      for i in xrange(self.nstrips):
         bar.update(i+1)
         sleep(0.001)
         yzcoord = self.yzrotation(y,self.board.z)
         xL = yzcoord[0]
         xR = yzcoord[1]
         yT = yzcoord[2]
         zT = yzcoord[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yT,zT)
         strip.SetPoint(1,xR,yT,zT)
         self.tstrips["YZrotationY"].append(strip)
         self.tstrips["YZrotationZ"].append(strip)
         ### compare yz-rotated pattern to nominal one in y coordinate
         yNominal = y
         yActual  = yT ### assume that the x of the test point is the same as in for the nominal point
         delta = yActual-yNominal
         self.histos["YZrotationY"].Fill(delta)
         self.RMS["YZrotationY"] += delta*delta
         self.Npoints["YZrotationY"] += 1
         ### compare yz-rotated pattern to nominal one in z coordinate
         zNominal = self.board.z
         zActual  = zT ### assume that the x of the test point is the same as in for the nominal point
         delta = zActual-zNominal
         self.histos["YZrotationZ"].Fill(delta)
         self.RMS["YZrotationZ"] += delta*delta
         self.Npoints["YZrotationZ"] += 1
         ###################
         ### propogate y...
         y += self.dy
      bar.finish()
      

############################################################
############################################################
############################################################
class PadLayer:
   'Class for pad layer'

   def __init__(self,name,board,color=ROOT.kGray+1,nstrips=404,npoints=100):
      # arguments
      self.name = name
      self.board = board
      self.nstrips = nstrips
      self.npoints = npoints
      self.color = color
      self.yorigin = -1
      self.xorigin = -1
      self.icenteralstrip = -1

      self.wstrip  = 2.7
      self.ystripmin = self.board.y0+15 ### 7.5*2 ??
      self.ystripmax = self.ystripmin+1292.8
      self.dy = (self.ystripmax-self.ystripmin)/self.nstrips

      ### for the transformations
      self.tstrips = {"Translation":[],"XYrotation":[],"Pitchscale":[],"Parallelism":[],"Bowing":[],"Zshift":[],"YZrotationY":[],"YZrotationZ":[]}
      self.Npoints = {"Translation":0, "XYrotation":0, "Pitchscale":0, "Parallelism":0,"Bowing":0, "Zshift":0, "YZrotationY":0, "YZrotationZ":0}
      self.RMS     = {"Translation":0, "XYrotation":0, "Pitchscale":0, "Parallelism":0,"Bowing":0, "Zshift":0, "YZrotationY":0, "YZrotationZ":0}
      self.histos  = {"Translation": TH1D(self.name+"_Pad_Translation", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "XYrotation":  TH1D(self.name+"_Pad_XYrotation",  ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Pitchscale":  TH1D(self.name+"_Pad_Pitchscale",  ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Parallelism": TH1D(self.name+"_Pad_Parallelism", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Bowing":      TH1D(self.name+"_Pad_Bowing",       ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "Zshift":      TH1D(self.name+"_Pad_Zshift",      ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "YZrotationY": TH1D(self.name+"_Pad_YZrotationY", ";y_{actual}-y_{nominal} [mm];Number of test points",50,-0.1,+0.1),
                      "YZrotationZ": TH1D(self.name+"_Pad_YZrotationZ", ";z_{actual}-z_{nominal} [mm];Number of test points",50,-0.1,+0.1)}
      self.offset = -1
      self.extendstrips = True
      self.raxis = [-1,-1]
      self.thetaxy = -1
      self.piratioxy = ""
      self.pivot = TEllipse(self.raxis[0],self.raxis[1],2)
      self.pitch = -1
      self.spitch = ""
      self.nonparallelism = -1
      self.snonparallelism = ""
      self.zshift = -1
      self.szshift = ""
      self.thetayz = -1
      self.sthetayz = ""
      self.bowpar = -1
      self.sbowpar = ""

      #################
      self.strips = []
      self.origin()
      self.pattern()
      #################

   def __del__(self):
      class_name = self.__class__.__name__
      # print class_name, "destroyed"

   def __cmp__(self, other):
      return cmp(self.name, other.name)
   def __copy__(self):
      # print '__copy__()'
      return PadLayer(self.name)
   def __deepcopy__(self, memo):
      # print '__deepcopy__(%s)' % str(memo)
      return PadLayer(copy.deepcopy(self.name, memo))

   def origin(self):
      ### find central strip y coordinate
      y = self.ystripmin
      self.xorigin = self.board.xcenter
      self.yorigin = y
      dycenter = 1.e20
      for i in xrange(self.nstrips):
         if(abs(y-self.board.ycenter)<dycenter):
            dycenter = abs(y-self.board.ycenter)
            self.yorigin = y
            self.icenteralstrip = i
         y += self.dy

   def pattern(self):
      y = self.ystripmin
      for i in xrange(self.nstrips):
         strip = TPolyLine3D(3)
         strip.SetPoint(0,self.board.edgex(y,"left"),y,self.board.z)
         strip.SetPoint(1,self.board.edgex(y,"right"),y,self.board.z)
         strip.SetLineColor(self.color)
         self.strips.append(strip)
         y += self.dy

   def ylinear(self,x,xa,ya,xb,yb):
      ### slope = (ya-yb)/(xa-xb)
      ### y-ya = slope*(x-xa)
      ### y = slope*x + intercept
      slope = (ya-yb)/(xa-xb)
      intercept = -slope*xa+ya
      y = slope*x + intercept
      return y

   def yparabole(self,x,xL,yL,xC,yC,xR,yR):
      X = np.array([[xL*xL,xL,1],
                    [xC*xC,xC,1],
                    [xR*xR,xR,1]])
      Xinv = linalg.inv(X)
      Y = np.array([yL,yC,yR])
      A = Xinv.dot(Y)
      y = A[0]*x*x + A[1]*x + A[2]
      return y

   def extend(self,xL,yL,xR,yR):
      slopeS = (yR-yL)/(xR-xL)
      interceptS = -slopeS*xR+yR
      slopeL = (self.board.y1-self.board.y0)/(self.board.x1-self.board.x0)
      interceptL = -slopeL*self.board.x1+self.board.y1
      slopeR = (self.board.y3-self.board.y2)/(self.board.x3-self.board.x2)
      interceptR = -slopeR*self.board.x3+self.board.y3
      xR = (interceptR-interceptS)/(slopeS-slopeR)
      yR = slopeR*xR+interceptR
      xL = (interceptL-interceptS)/(slopeS-slopeL)
      yL = slopeL*xL+interceptL
      return (xL,yL,xR,yR)

   def translation(self,y):
      y += self.offset
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      return (xL,xR,y)

   def xyrotattion(self,y,side):
      x = self.board.edgex(y,side)-self.raxis[0]
      y = y-self.raxis[1]
      xr = x*math.cos(self.thetaxy)-y*math.sin(self.thetaxy)
      yr = x*math.sin(self.thetaxy)+y*math.cos(self.thetaxy)
      return (xr+self.raxis[0],yr+self.raxis[1])

   def pitchscale(self,y,istrip,pitchsign):
      yP = y + pitchsign*istrip*self.pitch
      xL  = self.board.edgex(yP,"left")
      xR  = self.board.edgex(yP,"right")
      return (xL,xR,yP)

   def parallelism(self,y,istrip,npsign):
      yL = y
      xL = self.board.edgex(yL,"left")
      yR = y + npsign*istrip*self.nonparallelism
      xR = self.board.edgex(yR,"right")
      return (xL,yL,xR,yR)

   def bowing(self,y,nxpoints,istrip,bsign):
      yL = y
      xL = self.board.edgex(yL,"left")
      yR = y
      xR = self.board.edgex(yR,"right")
      xC = xL+(xR-xL)/2.
      yC = y + bsign*istrip*self.bowpar
      dx = (xR-xL)/nxpoints
      xyvec = []
      for i in xrange(nxpoints):
         xi = xL+i*dx
         yi = self.yparabole(xi,xL,yL,xC,yC,xR,yR)
         xyvec.append([xi,yi])
      xyvec.append([xR,yR])
      return xyvec

   def longidudinalshift(self,y,z,shiftsign):
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      z = z + shiftsign*self.zshift
      return (xL,xR,y,z)

   def yzrotation(self,y,z):
      ### do not change X !
      xL = self.board.edgex(y,"left")
      xR = self.board.edgex(y,"right")
      rsigny = 1.
      rsignz = 1.
      if(y>self.yorigin):
         rsigny = -1.
         rsignz = +1.
      elif(y<self.yorigin):
         rsigny = +1.
         rsignz = -1.
      else:
         rsigny = 0.
         rsignz = 0.
      # yS = self.yorigin + rsign*abs(y)*math.cos(self.thetayz)
      # zS = self.board.z + rsign*abs(y)*math.sin(self.thetayz)
      dz = abs((y-self.yorigin)*math.sin(self.thetayz))
      zS = self.board.z + rsignz*dz
      yS = y            + rsigny*abs(dz*math.tan(self.thetayz))
      return (xL,xR,yS,zS)


   def transform(self,color=ROOT.kRed-2,offset=30.,raxis=[180.,100.],thetaxy=math.pi/45000.,pitchum=75.,nonparallelismum=50.,zshift=10.,thetayz=math.pi/45000.,bowpar=50.):
      # print "transform() called !"
	
      ### for the offset translation
      self.offset = offset/1000.

      ### for the rotation
      self.extendstrips = True
      self.raxis = raxis
      self.thetaxy = thetaxy
      self.piratioxy = str(thetaxy/math.pi)
      self.pivot = TEllipse(raxis[0],raxis[1],2)

      ### for the pitch scaling (+-75um)
      self.pitch = (pitchum/self.nstrips)/1000. ## per strip
      self.spitch = "75#mum"

      ### for the non-parallelism (+-50um)
      self.nonparallelism = (nonparallelismum/self.nstrips)/1000. ## per strip
      self.snonparallelism = "50#mum"

      ### for the offset translation
      self.zshift = zshift/1000. ## per layer
      self.szshift = "10#mum"

      ### for the y-z rotation
      self.thetayz = thetayz
      self.piratioyz = str(thetayz/math.pi)

      ### for the bowing
      self.bowpar = (bowpar/self.nstrips)/1000.
      self.sbowpar = "50#mum"
      self.nxpoints = 100

      #################################
      ### Start the transformations ###
      #################################

      y = self.ystripmin
      for i in xrange(self.nstrips):
         ### translated pattern
         tcoord = self.translation(y)
         xL = tcoord[0]
         xR = tcoord[1]
         yT = tcoord[2]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yT,self.board.z)
         strip.SetPoint(1,xR,yT,self.board.z)
         self.tstrips["Translation"].append(strip)
         ### compare translated pattern to nominal one
         yNominal = y
         yActual  = yT ### assume that the x of the test point is the same as in for the nominal point
         delta = yActual-yNominal
         self.histos["Translation"].Fill(delta)
         self.RMS["Translation"] += delta*delta
         self.Npoints["Translation"] += 1
         ###################
         ### propogate y...
         y += self.dy

      y = self.ystripmin
      for i in xrange(self.nstrips):
         ### rotated pattern
         rcoordR = self.xyrotattion(y,"right")
         rcoordL = self.xyrotattion(y,"left")
         xR = rcoordR[0]
         yR = rcoordR[1]
         xL = rcoordL[0]
         yL = rcoordL[1]
         if(self.extendstrips):
            rcoord = self.extend(xL,yL,xR,yR)
            xL = rcoord[0]
            yL = rcoord[1]
            xR = rcoord[2]
            yR = rcoord[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yL,self.board.z)
         strip.SetPoint(1,xR,yR,self.board.z)
         self.tstrips["XYrotation"].append(strip)
         ### compare rotated pattern to nominal one
         xminstrip = self.board.edgex(y,"left")
         xmaxstrip = self.board.edgex(y,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            yNominal = y
            xR = rcoordR[0]
            yR = rcoordR[1]
            xL = rcoordL[0]
            yL = rcoordL[1]
            if(self.extendstrips):
               rcoord = self.extend(xL,yL,xR,yR)
               xL = rcoord[0]
               yL = rcoord[1]
               xR = rcoord[2]
               yR = rcoord[3]
            x = xminstrip + p*dx
            yActual = self.ylinear(x,xL,yL,xR,yR)
            delta = yActual-yNominal
            self.histos["XYrotation"].Fill(delta)
            self.RMS["XYrotation"] += delta*delta
            self.Npoints["XYrotation"] += 1
         ###################
         ### propogate y...
         y += self.dy


      ### do the pitch scaling transformation
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         ### pitch-scaled pattern (going up)
         if(i==0): pcoordup = self.pitchscale(yup,i,0.)
         else:     pcoordup = self.pitchscale(yup,i,+1.)
         pcoorddn = self.pitchscale(ydn,i,-1.)
         xLup = pcoordup[0]
         xRup = pcoordup[1]
         yPup = pcoordup[2]
         xLdn = pcoorddn[0]
         xRdn = pcoorddn[1]
         yPdn = pcoorddn[2]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLup,yPup,self.board.z)
         strip.SetPoint(1,xRup,yPup,self.board.z)
         self.tstrips["Pitchscale"].append(strip)
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLdn,yPdn,self.board.z)
         strip.SetPoint(1,xRdn,yPdn,self.board.z)
         self.tstrips["Pitchscale"].append(strip)
         ### compare pitch-scaled up pattern to nominal one
         yNominal = yup
         yActual  = yPup
         delta = yActual-yNominal
         self.histos["Pitchscale"].Fill(delta)
         self.RMS["Pitchscale"] += delta*delta
         self.Npoints["Pitchscale"] += 1
         ### compare pitch-scaled down pattern to nominal one
         yNominal = ydn
         yActual  = yPdn
         delta = yActual-yNominal
         self.histos["Pitchscale"].Fill(delta)
         self.RMS["Pitchscale"] += delta*delta
         self.Npoints["Pitchscale"] += 1
         # print "yup=%g --> yPup=%g, ydn=%g --> yPdn=%g --> dyup=%g, dydn=%g" % (yup,yPup,ydn,yPdn,(yPup-yup),(ydn-yPdn))
         ###################
         ### propogate y...
         yup += self.dy
         ydn -= self.dy


      ### do the non-parallel transformation
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         if(i==0): npcoordup = self.parallelism(yup,i,0.)
         else:     npcoordup = self.parallelism(yup,i,+1.)
         npcoorddn = self.parallelism(ydn,i,-1.)
         xLup = npcoordup[0]
         yLup = npcoordup[1]
         xRup = npcoordup[2]
         yRup = npcoordup[3]
         xLdn = npcoorddn[0]
         yLdn = npcoorddn[1]
         xRdn = npcoorddn[2]
         yRdn = npcoorddn[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLup,yLup,self.board.z)
         strip.SetPoint(1,xRup,yRup,self.board.z)
         self.tstrips["Parallelism"].append(strip)
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xLdn,yLdn,self.board.z)
         strip.SetPoint(1,xRdn,yRdn,self.board.z)
         self.tstrips["Parallelism"].append(strip)
         # print "[%d]  yup=%g-->yRup=%g, ydn=%g-->yRdn=%g --> dyup=%g, dydn=%g" % (i,yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
         ### compare non-parallel pattern up to nominal one
         xminstrip = self.board.edgex(yup,"left")
         xmaxstrip = self.board.edgex(yup,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            x = xminstrip + p*dx
            yNominal = yup
            yActual = self.ylinear(x,xLup,yLup,xRup,yRup)
            delta = yActual-yNominal
            self.histos["Parallelism"].Fill(delta)
            self.RMS["Parallelism"] += delta*delta
            self.Npoints["Parallelism"] += 1
            # print "up[%d] yup=%g --> yActual=%g --> dyup=%g" % (p,yup,yActual,delta)
         ### compare non-parallel pattern down to nominal one
         xminstrip = self.board.edgex(ydn,"left")
         xmaxstrip = self.board.edgex(ydn,"right")
         dx = (xmaxstrip-xminstrip)/self.npoints
         for p in xrange(self.npoints):
            x = xminstrip + p*dx
            yNominal = ydn
            yActual = self.ylinear(x,xLdn,yLdn,xRdn,yRdn)
            delta = yActual-yNominal
            self.histos["Parallelism"].Fill(delta)
            self.RMS["Parallelism"] += delta*delta
            self.Npoints["Parallelism"] += 1
            # print "dn[%d] ydn=%g --> yActual=%g --> dydn=%g" % (p,ydn,yActual,delta)
         # print "RMS[%d]=%g" % (i,self.RMS["Parallelism"])
         ##################
	     ### propogate y...
         yup += self.dy
         ydn -= self.dy


      ### do the bowing transformation
      yup = self.yorigin
      ydn = self.yorigin-self.dy
      for i in xrange(self.nstrips/2):
         if(i==0): bcoordup = self.bowing(yup,self.nxpoints,i,0.)
         else:     bcoordup = self.bowing(yup,self.nxpoints,i,+1.)
         bcoorddn = self.bowing(ydn,self.nxpoints,i,-1.)
         stripU = TPolyLine3D(self.nxpoints+1)
         stripD = TPolyLine3D(self.nxpoints+1)
         for j in xrange(self.nxpoints):
            xJup = bcoordup[j][0]
            yJup = bcoordup[j][1]
            xJdn = bcoorddn[j][0]
            yJdn = bcoorddn[j][1]
            stripU.SetPoint(j,xJup,yJup,self.board.z)
            stripD.SetPoint(j,xJdn,yJdn,self.board.z)
         self.tstrips["Bowing"].append(stripU)
         self.tstrips["Bowing"].append(stripD)
         # print "[%d]  yup=%g-->yRup=%g, ydn=%g-->yRdn=%g --> dyup=%g, dydn=%g" % (i,yup,yRup,ydn,yRdn,(yRup-yup),(ydn-yRdn))
         ### compare bowing pattern up to nominal one
         for j in xrange(self.nxpoints):
            yJup = bcoordup[j][1]
            yJdn = bcoorddn[j][1]
            yNominal = yup
            yActual = yJup
            delta = yActual-yNominal
            self.histos["Bowing"].Fill(delta)
            self.RMS["Bowing"] += delta*delta
            self.Npoints["Bowing"] += 1
            # print "up[%d] yup=%g --> yActual=%g --> dyup=%g" % (p,yup,yActual,delta)
            yNominal = ydn
            yActual = yJdn
            delta = yActual-yNominal
            self.histos["Bowing"].Fill(delta)
            self.RMS["Bowing"] += delta*delta
            self.Npoints["Bowing"] += 1
            # print "dn[%d] ydn=%g --> yActual=%g --> dydn=%g" % (p,ydn,yActual,delta)
         ##################
	     ### propogate y...
         yup += self.dy
         ydn -= self.dy


      ### do the longitudinal shift transformation
      y = self.ystripmin
      for i in xrange(self.nstrips):
         scoords = self.longidudinalshift(y,self.board.z,+1)
         xL = scoords[0]
         xR = scoords[1]
         y  = scoords[2]
         zS = scoords[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,y,zS)
         strip.SetPoint(1,xR,y,zS)
         self.tstrips["Zshift"].append(strip)
         ### compare z-shifted pattern to nominal one
         zNominal = self.board.z
         zActual  = zS
         delta = zActual-zNominal
         self.histos["Zshift"].Fill(delta)
         self.RMS["Zshift"] += delta*delta
         self.Npoints["Zshift"] += 1
         ###################
         ### propogate y...
         y += self.dy

      ### do the yz-rotatied transformation
      y = self.ystripmin
      for i in xrange(self.nstrips):
         yzcoord = self.yzrotation(y,self.board.z)
         xL = yzcoord[0]
         xR = yzcoord[1]
         yT = yzcoord[2]
         zT = yzcoord[3]
         strip = TPolyLine3D(3)
         strip.SetPoint(0,xL,yT,zT)
         strip.SetPoint(1,xR,yT,zT)
         self.tstrips["YZrotationY"].append(strip)
         self.tstrips["YZrotationZ"].append(strip)
         ### compare yz-rotated pattern to nominal one in y coordinate
         yNominal = y
         yActual  = yT ### assume that the x of the test point is the same as in for the nominal point
         delta = yActual-yNominal
         self.histos["YZrotationY"].Fill(delta)
         self.RMS["YZrotationY"] += delta*delta
         self.Npoints["YZrotationY"] += 1
         ### compare yz-rotated pattern to nominal one in z coordinate
         zNominal = self.board.z
         zActual  = zT ### assume that the x of the test point is the same as in for the nominal point
         delta = zActual-zNominal
         self.histos["YZrotationZ"].Fill(delta)
         self.RMS["YZrotationZ"] += delta*delta
         self.Npoints["YZrotationZ"] += 1
         ###################
         ### propogate y...
         y += self.dy
