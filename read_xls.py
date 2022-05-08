from __future__ import absolute_import, division, print_function
import openpyxl
import sys, os
import mmtbx.secondary_structure
import csv
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import null_out
from matplotlib import pyplot as plt
from matplotlib import patches as patch
from libtbx import group_args
from mmtbx.validation import ramalyze
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.conformation_dependent_library import generate_protein_threes

fmt_helix = "HELIX%5d%4d %3s %s%5d  %3s %s%5s%3s                                 %3d"

def get_ids(atoms, i):
  ai = atoms[i]
  return group_args(
    resname = ai.parent().resname,
    resseq  = ai.parent().parent().resseq_as_int(),
    cid     = ai.parent().parent().parent().id)


# read xls file ( modified helices records)
def read_xls(file_name):
  helix_dict = {}
  workbook = openpyxl.load_workbook(file_name)
  for sheet in workbook.worksheets:
    data = list(sheet.columns)
    j = 1
    for datas in data[1][1:]:
      residue = data[3][j].value
      if datas.value != None:
        key = str(datas.value)
        helix_dict[key] = []
      if residue != None and str(residue) != 'delete':
        residue = [int(i) for i in residue.split('-')]
        helix_dict[key].append(residue)
      j += 1
    print(helix_dict)
  return helix_dict

# output pdb files according to xls file
def output_helix_records(helix_dic):
  input_path = './PDB/standard_set/'
  save_path = './PDB/manual_set/'
  for name in helix_dic.keys():
    cntr = 0
    helix_recs = []
    file_name = input_path + name + '.pdb'
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    hierarchy = pdb_inp.construct_hierarchy()
    hierarchy.remove_alt_confs(always_keep_one_conformer=True)
    save_file = save_path + name +'_new.pdb'
    for value in helix_dic[name]:
      cntr += 1
      resseq = (' resseq ' + str(int(value[0])) + ':' + str(int(value[1])) + ' ')
      ase = hierarchy.atom_selection_cache().selection(resseq)
      pdb_hierarchy = hierarchy.select(ase)
      atoms = pdb_hierarchy.atoms()
      start = get_ids(atoms=atoms, i=0)
      stop = get_ids(atoms=atoms, i=-1)
      length = stop.resseq - start.resseq + 1
      h_rec = fmt_helix % (
        cntr, cntr,
        start.resname, start.cid, start.resseq,
        stop.resname, stop.cid, stop.resseq,
        1, length)
      helix_recs.append(h_rec)
    helix_recs_str = "\n".join(helix_recs)
    print(helix_recs_str)
    with open(save_file, "w") as fo:
      fo.write(helix_recs_str)
      fo.write("\n")
      fo.write(hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry()))

def run(args):
  helix_dic = read_xls(args[0])
  output_helix_records(helix_dic)

if __name__ == '__main__':
  run(args=sys.argv[1:])
