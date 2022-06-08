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



def compare_helices(file1, file2, count_dict):
  pdb_inp_1 = iotbx.pdb.input(file_name = file1)
  helices_1 = pdb_inp_1.extract_secondary_structure().helices
  pdb_inp_2 = iotbx.pdb.input(file_name = file2)
  helices_2 = pdb_inp_2.extract_secondary_structure().helices
  output = []
  find_in_h2 = []
  format = "%-2s %-9s %6s %-16s %-18s"
  for h1 in helices_1:
    status = 'unget'
    helix_id = h1.helix_id
    chain_1 = h1.start_chain_id
    resseq_1 = int(h1.start_resseq)
    for h2 in helices_2:
      abs_length = h1.length - h2.length
      chain_2 = h2.start_chain_id
      resseq_2 = int(h2.start_resseq)
      if chain_1==chain_2:
        if (abs(resseq_1 - resseq_2) <= 4):
          status = 'get'
          if h1.start_resseq == h2.start_resseq:
            if h1.end_resseq == h2.end_resseq:
              status = 'identical'
              length = abs_length
              start = h1.start_resname + h1.start_resseq
              end = h1.end_resname + h1.end_resseq
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['identical'] += 1
            elif h1.end_resseq > h2.end_resseq:
              status = 'longer'
              length = abs_length
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['longer'] += 1
              if abs_length > count_dict['max_long']:
                count_dict['max_long'] = abs_length
            else:
              status = 'shoter'
              length = abs_length
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['shorter'] += 1
              if abs_length < count_dict['max_short'] and abs_length!=-5:
                count_dict['max_short'] = abs_length
          else:
            if h1.length == h2.length:
              status = 'same'
              length = abs_length
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['same'] += 1
            elif h1.length > h2.length:
              status = 'longer'
              length = abs_length
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['longer'] += 1
              if abs_length > count_dict['max_long']:
                count_dict['max_long'] = abs_length
            else:
              status = 'shorter'
              length = abs_length
              start = h1.start_resname + h1.start_resseq + '(' + h2.start_resname + h2.start_resseq  + ')'
              end = h1.end_resname + h1.end_resseq + '(' + h2.end_resname + h2.end_resseq + ')'
              out = format % (helix_id, status, length, start, end)
              output.append(out)
              find_in_h2.append(h2)
              count_dict['shorter'] += 1
              if abs_length < count_dict['max_short'] and abs_length!=-5:
                count_dict['max_short'] = abs_length
    if status == 'unget':
      status = 'lost'
      length = 'None'
      start = h1.start_resname + h1.start_resseq
      end = h1.end_resname + h1.end_resseq
      out = format % (helix_id, status, length, start, end)
      output.append(out)
      count_dict['lost'] += 1
  for h3 in helices_2:
    if h3 not in find_in_h2:
      status = 'more'
      length = 'None'
      helix_id = 'None'
      start = h3.start_resname + h3.start_resseq
      end = h3.end_resname + h3.end_resseq
      out = format % (helix_id, status, length, start, end)
      output.append(out)
      count_dict['more'] += 1
  title = format % ('ID', 'Status', 'length', 'start_res', 'end_res')
  file_name = file1[-12:-8]
  print('file name:', file_name)
  print(title)
  for record in output:
    print(record)
  print('\n')
  return output

def run(args):
  count_dict = {'longer':0, 'shorter':0, 'same':0, 'identical':0, 'lost':0, 'more':0, 'max_long':0, 'max_short':0}
  for root, dirs, files in os.walk(args[0]):
    for file in files:
      file_name_1 = root +  file
      file_name_2 = args[1] + file
      compare_helices(file_name_1, file_name_2, count_dict)
  print(count_dict)
if __name__ == '__main__':
  run(args=sys.argv[1:])
