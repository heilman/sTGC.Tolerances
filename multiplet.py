#!/usr/bin/python

import ROOT
from ROOT import *
import math
from array import array
from board import StripBoard
from layer import StripLayer
import copy

class Multiplet:
   'Class for strip/pad multiplet'

   def __init__(self,name,finflate=1,dzStripLayers=10.97,center=[0,0,0],angles=[0,0]):
      # arguments
      self.name = name
      self.center = center
      self.centerxy = [self.center[0],self.center[1]]
      self.angles = angles
      self.dzStripLayers = dzStripLayers
      self.zsboards = [self.center[2] + 0.*self.dzStripLayers,
                       self.center[2] + 1.*self.dzStripLayers,
                       self.center[2] + 2.*self.dzStripLayers,
                       self.center[2] + 3.*self.dzStripLayers]
      self.sboards = []
      self.slayers = []
      self.nslayers = -1
      self.finflate = finflate

      #################
      ### add all items
      self.add()
      #################

   def __del__(self):
      class_name = self.__class__.__name__
      # print class_name, "destroyed"

   def __cmp__(self, other):
      return cmp(self.name, other.name)
   def __copy__(self):
      # print '__copy__()'
      return Multiplet(self.name)
   def __deepcopy__(self, memo):
      # print '__deepcopy__(%s)' % str(memo)
      return Multiplet(copy.deepcopy(self.name, memo))

   def add(self):
      ### declare the multiplet objects
      for i in xrange(len(self.zsboards)):
         name = self.name+"."+str(i)
         sboard = StripBoard(name,self.zsboards[i],self.centerxy)
         self.sboards.append(sboard)
         self.slayers.append(StripLayer(name,self.sboards[i]))
         ### call all transformations
         self.slayers[i].transform( ### OPTION A
                                    name+".strip", 
                                    ROOT.kRed-2,
                                    75.*self.finflate,            ## offset
                                    25.*self.finflate,            ## offset correction 
                                    [180.,100.],                  ## raxis(xy)
                                    math.pi/45000.*self.finflate, ## theta(xy)
                                    75.*self.finflate,            ## pitchum
                                    75.*self.finflate,            ## nonparallelismum
                                    50.*self.finflate,            ## zshift
                                    math.pi/18000.*self.finflate, ## theta(yz)
                                    75.*self.finflate             ## bowpar
                                  )
         # self.slayers[i].transform( ### OPTION B
         #                            name+".strip", 
         #                            ROOT.kRed-2,
         #                            150.*self.finflate,            ## offset
         #                            25.*self.finflate,            ## offset correction 
         #                            [180.,100.],                  ## raxis(xy)
         #                            math.pi/45000.*self.finflate, ## theta(xy)
         #                            75.*self.finflate,            ## pitchum
         #                            75.*self.finflate,            ## nonparallelismum
         #                            50.*self.finflate,            ## zshift
         #                            math.pi/18000.*self.finflate, ## theta(yz)
         #                            75.*self.finflate             ## bowpar
         #                          )
         # self.slayers[i].transform( ### OPTION C
         #                            name+".strip", 
         #                            ROOT.kRed-2,
         #                            75.*self.finflate,            ## offset
         #                            25.*self.finflate,            ## offset correction 
         #                            [180.,100.],                  ## raxis(xy)
         #                            math.pi/45000.*self.finflate, ## theta(xy)
         #                            150.*self.finflate,            ## pitchum
         #                            150.*self.finflate,            ## nonparallelismum
         #                            50.*self.finflate,            ## zshift
         #                            math.pi/18000.*self.finflate, ## theta(yz)
         #                            150.*self.finflate             ## bowpar
         #                          )
         # self.slayers[i].transform( ### OPTION D
         #                            name+".strip", 
         #                            ROOT.kRed-2,
         #                            75.*self.finflate,            ## offset
         #                            25.*self.finflate,            ## offset correction 
         #                            [180.,100.],                  ## raxis(xy)
         #                            math.pi/22500.*self.finflate, ## theta(xy)
         #                            75.*self.finflate,            ## pitchum
         #                            75.*self.finflate,            ## nonparallelismum
         #                            50.*self.finflate,            ## zshift
         #                            math.pi/18000.*self.finflate, ## theta(yz)
         #                            75.*self.finflate             ## bowpar
         #                          )
      ### layers counter
      self.nslayers = len(self.slayers)