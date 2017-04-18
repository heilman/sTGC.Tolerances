#!/usr/bin/python

import ROOT
from ROOT import *
import math
from array import array
from board import Board
from layer import Layer
import copy

class Multiplet:
   'Class for strip/pad multiplet'

   def __init__(self,name,dzStripLayers=10.97,dzPadLayers=10.97,center=[0,0,0],angles=[0,0]):
      # arguments
      self.name = name
      self.center = center
      self.centerxy = [self.center[0],self.center[1]]
      self.angles = angles
      self.dzStripLayers = dzStripLayers
      self.dzPadLayers = dzPadLayers
      self.zsboards = [self.center[2] + 0.*self.dzStripLayers,
                       self.center[2] + 1.*self.dzStripLayers,
                       self.center[2] + 2.*self.dzStripLayers,
                       self.center[2] + 3.*self.dzStripLayers]
      self.zpboards = [self.center[2] + 0.*self.dzPadLayers,
                       self.center[2] + 1.*self.dzPadLayers,
                       self.center[2] + 2.*self.dzPadLayers,
                       self.center[2] + 3.*self.dzPadLayers]
      self.sboards = []
      self.pboards = []
      self.slayers = []
      self.players = []
      self.nlayers = -1

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
      for i in xrange(len(self.zboards)):
         name = self.name+"."+str(i)
         self.sboards.append(StripBoard(name,self.zsboards[i],self.centerxy))
         self.pboards.append(PadBoard(name,self.zpboards[i],self.centerxy))
         self.slayers.append(StripLayer(name,self.sboards[i]))
         self.players.append(PadLayer(name,self.pboards[i]))
         self.slayers[i].transform()
         self.players[i].transform()
      ### layers counter
      self.nlayers = len(self.layers)