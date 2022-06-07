#!/usr/bin/env python
# coding: utf-8


from __future__ import absolute_import, division, print_function
import os, sys
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



def compare_helices(file1, file2):
  pdb_inp_1 = iotbx.pdb.input(file_name = file1)
  helices_1 = pdb_inp_1.extract_secondary_structure().helices
  pdb_inp_2 = iotbx.pdb.input(file_name = file2)
  helices_2 = pdb_inp_2.extract_secondary_structure().helices
  output = []
  find_in_h2 = []
  format = "%-5s %-6s %-10s %-18s %-18s"
  for h1 in helices_1:
    status = 'unget'
    helix_id = h1.helix_id
    chain_1 = h1.start_chain_id
    resseq_1 = int(h1.start_resseq)
    for h2 in helices_2:
      chain_2 = h2.start_chain_id
      resseq_2 = int(h2.start_resseq)
      if chain_1==chain_2:
        if (abs(resseq_1 - resseq_2) < 4):
          status = 'get'
          if h1.start_resname == h2.start_resname:
            if h1.start_resseq == h2.start_resseq:
              length = 'identical'
              start = h1.start_resname + h1.start_resseq
              end = h1.end_resname + h1.end_resseq
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
            elif h1.start_resseq > h2.start_resseq:
              length = 'longer'
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
            else:
              length = 'shorter'
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
          else:
            if h1.start_resseq == h2.start_resseq:
              length = 'same'
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
            elif h1.start_resseq > h2.start_resseq:
              length = 'longer'
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
            else:
              length = 'shorter'
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
    if status == 'unget':
      status = 'less'
      length = 'None'
      start = h1.start_resname + h1.start_resseq
      end = h1.end_resname + h1.end_resseq
      out = format % (helix_id, status, length, start, end)
      output.append(out)
      find_in_h2.append(h2)
  for h3 in helices_2:
    if h3 not in find_in_h2:
      status = 'more'
      length = 'None'
      helix_id = 'None'
      start = h3.start_resname + h3.start_resseq
      end = h3.end_resname + h3.end_resseq
      out = format % (helix_id, status, length, start, end)
      output.append(out)
  title = format % ('ID', 'Status', 'length', 'start_res', 'end_res')
  file_name = file1[-12:-8]
  print('file name:', file_name)
  print(title)
  for record in output:
    print(record)
  print('\n')
  return output

def run(args):
  for root, dirs, files in os.walk(args[0]):
    for file in files:
      file_name_1 = root +  file
      file_name_2 = args[1] + file
      compare_helices(file_name_1, file_name_2)

if __name__ == '__main__':
  run(args=sys.argv[1:])

