#!/usr/bin/python

import ROOT
from ROOT import *
import math
from array import array
import copy

class StripBoard:
   'Class for strip board'

   def __init__(self,name,z,center=[0,0],color=ROOT.kGreen+2):
      # arguments
      self.name = name
      self.z = z
      self.center = center
      self.color = color

      self.dimxtop = 743
      self.dimxbot = 346.9
      self.dimy    = 1325.6
      self.scale   = 1.0
      self.xcenter = self.center[0] + self.dimxtop/2.
      self.ycenter = self.center[1] + self.dimy/2.
      
      self.x0 = self.center[0] + ((self.scale*self.dimxtop)-self.dimxbot)/2.
      self.y0 = self.center[1] + ((self.scale*self.dimy)-self.dimy)/2.
      self.x1 = self.center[0] + self.x0-(self.dimxtop-self.dimxbot)/2.
      self.y1 = self.center[1] + self.y0+self.dimy
      self.x2 = self.center[0] + self.x0+self.dimxbot+(self.dimxtop-self.dimxbot)/2.
      self.y2 = self.center[1] + self.y0+self.dimy
      self.x3 = self.center[0] + self.x0+self.dimxbot
      self.y3 = self.center[1] + self.y0

      self.ncorners = 5
      self.xcorners = array( 'f', [self.x0,self.x1,self.x2,self.x3,self.x0] )
      self.ycorners = array( 'f', [self.y0,self.y1,self.y2,self.y3,self.y0] )
      self.zcorners = array( 'f', [self.z, self.z, self.z, self.z, self.z ] )
      self.board = TPolyLine3D(self.ncorners, self.xcorners, self.ycorners, self.zcorners, "F")
      self.board.SetLineColor(color)

   def __del__(self):
      class_name = self.__class__.__name__
      # print class_name, "destroyed"

   def __cmp__(self, other):
      return cmp(self.name, other.name)
   def __copy__(self):
      # print '__copy__()'
      return StripBoard(self.name)
   def __deepcopy__(self, memo):
      # print '__deepcopy__(%s)' % str(memo)
      return StripBoard(copy.deepcopy(self.name, memo))

   def edgex(self,y,side):
      ### slope = (y1-y0)/(x1-x0)
      ### y-y1 = slope*(x-x1)
      ### y = slope*x + intercept
      if(side=="left"):
         slope = (self.y1-self.y0)/(self.x1-self.x0)
         intercept = -slope*self.x1+self.y1
      elif(side=="right"):
         slope = (self.y3-self.y2)/(self.x3-self.x2)
         intercept = -slope*self.x3+self.y3
      else:
         print "wrong choice:",side
         quit()
      x = (y-intercept)/slope
      return x