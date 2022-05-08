from __future__ import absolute_import, division, print_function
import math
import copy
import csv
import mmtbx.model
import numpy as np
import iotbx
from numpy import nan
from libtbx.utils import null_out
from collections import OrderedDict
from scitbx.array_family import flex
from mmtbx.conformation_dependent_library import generate_protein_fragments


class AEV(object):
  """
  Smith J S, Isayev O, Roitberg A E. ANI-1: an extensible neural network potential with DFT
  accuracy at force field computational cost[J]. Chemical science, 2017, 8(4): 3192-3203.
  """

  def __init__(self,
               model,
               rs_values=(2.0, 3.8, 5.2, 5.5, 6.2, 7.0, 8.6, 10.0),
               # probe distances (A) for radial
               radial_eta=4,
               cutoff=8.1,
               # radial cutoff distance
               ts_values=(0.392699, 1.178097, 1.963495, 2.748894),
               # probe angles (rad) for angular
               angular_rs_values=(3.8, 5.2, 5.5, 6.2),
               # probe distances (A) for angular
               angular_eta=4,
               angular_zeta=8,
               # parameters for probe angles
               ):
    self.hierarchy = model.get_hierarchy()
    self.rs_values = rs_values
    self.crystal_symmetry = model.crystal_symmetry()
    self.radial_eta = radial_eta
    self.cutoff = cutoff
    self.ts_values = ts_values
    self.angular_rs_values = angular_rs_values
    self.angular_eta = angular_eta
    self.angular_zeta = angular_zeta
    self.EAEVs = format_class(length_of_radial=len(self.rs_values))
    self.MAEVs = format_class(length_of_radial=len(self.rs_values))
    self.BAEVs = format_class(length_of_radial=len(self.rs_values))
    self.center_atom = None
    self.chain_hierarchy = None
    self.generate_AEV()

  def get_values(self):
    result = OrderedDict()
    result['B'] = list(self.BAEVs.values())[0]
    result['M'] = list(self.MAEVs.values())[5]
    result['E'] = list(self.EAEVs.values())[-1]
    return result

  def empty_dict(self):
    results = OrderedDict()
    for atom in self.hierarchy.atoms():
      if atom.name == ' CA ':
        res_name = atom.format_atom_record()[22:26]
        results.setdefault(res_name, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
      self.BAEVs.update(results)
      self.MAEVs.update(results)
      self.EAEVs.update(results)
    return 0

  def generate_ca(self, length=5):
    hierarchy = self.chain_hierarchy
    for chain in hierarchy.chains():
      con = chain.conformers()[0]
      b = int(length)
      for a in range(len(con.atoms())):
        rc = []
        if b <= len(con.atoms()):
          b = a + int(length)
          for atom in con.atoms()[a:b]:
            if atom.name == ' CA ':
              rc.append(atom)
        if len(rc) == int(b - a):
          yield rc

  def generate_AEV(self):
    """
    Traversal the C-alpha atoms list and generate AEV values of every chain's
    every atom. Every C-alpha atom has 3 direction AEVs: forward(BAEVs),
    backward(EAEVs) and all direction(MAEVS).
    """

    chain_list = []
    self.empty_dict()
    for chain in self.hierarchy.chains():
      chain_list.append(chain.id)
    chain_list = list(set(chain_list))
    for chain in chain_list:
      chain_hierarchy = self.hierarchy.deep_copy()
      asc = chain_hierarchy.atom_selection_cache()
      sel = asc.selection("protein and name CA and chain " + chain)
      self.chain_hierarchy = chain_hierarchy.select(sel)
      for atomlist in self.generate_ca():
        self.center_atom = atomlist[0]
        self.BAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[-1]
        self.EAEVs.update(self.calculate(atomlist))
        self.center_atom = atomlist[2]
        self.MAEVs.update(self.calculate(atomlist))

  def cutf(self, distance):
    """
    Formula (2), page 3194
    """
    if distance <= self.cutoff:
      Fc = 0.5 * math.cos(math.pi * distance / self.cutoff) + 0.5
    else:
      Fc = 0
    return Fc

  def calculate(self, atom_list):
    """
    Formula (3) and (4), page 3194
    """
    n = self.radial_eta
    assert self.radial_eta == self.angular_eta
    l = self.angular_zeta
    AEVs = format_class()
    res_name = self.center_atom.format_atom_record()[22:26]
    AEVs.setdefault(res_name, [])
    if atom_list != []:
      atom1 = self.center_atom
      atomlist = copy.copy(atom_list)
      if atom1 in atomlist:
        atomlist.remove(atom1)
      for Rs in self.rs_values:
        GmR = 0
        for atom2 in atomlist:
          R = atom1.distance(atom2)
          f = self.cutf(R)
          if f != 0:
            mR = math.exp(- n * ((R - Rs) ** 2)) * f
            GmR += mR
        AEVs[res_name].append(GmR)
      for Rs in self.angular_rs_values:
        for zetas in self.ts_values:
          i = 0
          GmA = 0
          zeta_list = []
          atomlist = copy.copy(atom_list)
          if atom1 in atomlist:
            atomlist.remove(atom1)
          for atom2 in atomlist[::-1]:
            if atom2 in atomlist:
              atomlist.remove(atom2)
            for atom3 in atomlist:
              Rij = atom1.distance(atom2)
              Rik = atom1.distance(atom3)
              ZETAijk = atom1.angle(atom2, atom3)
              if ZETAijk != 0:
                i += 1
                fk = self.cutf(Rik)
                fj = self.cutf(Rij)
                if fk != 0 and fj != 0:
                  zeta_list.append(ZETAijk)
                  mA = (((1 + math.cos(ZETAijk - zetas))) ** l) * \
                       math.exp(- n * ((((Rij + Rik) / 2) - Rs) ** 2)) * fj * fk
                  GmA += mA
          GmA = GmA * (2 ** (1 - l))
          AEVs[res_name].append(GmA)
    return AEVs
class format_class(OrderedDict):
  def __init__(self, length_of_radial=None):
    OrderedDict.__init__(self)
    self.length_of_radial=length_of_radial


  def __repr__(self):
    outl = '...\n'
    for key, item in sorted(self.items()):
      outl += '  %s :' % (key)
      for i, v in enumerate(item):
        if i==self.length_of_radial: outl+='|'
        outl += '%0.4f, ' % v
      outl += '\n'
    return outl
