from __future__ import absolute_import, division, print_function
import sys, os
import mmtbx.secondary_structure
from scitbx.array_family import flex
from libtbx.utils import null_out
import iotbx.pdb
from scitbx.matrix import col
from libtbx import group_args
from mmtbx.validation import ramalyze
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.conformation_dependent_library import generate_protein_threes

fmt_helix="HELIX%5d%4d %3s %s%5d  %3s %s%5s%3s                                 %3d"

def get_indexed_residue(hierarchy,index=0):
  if not hierarchy:
    return None
  count=0
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        for residue in conformer.residues():
          if count==index:
            return residue
          count+=1
def get_atom_from_residue(residue=None,
  atom_name=None):
  if not residue:
    return None,None
  # just make sure that atom_name is there
  for atom in residue.atoms():
    if atom.name.replace(" ","")==atom_name.replace(" ",""):
      return atom.name,atom.xyz

def get_helices(pdb_inp, hierarchy):
  ss = pdb_inp.extract_secondary_structure().filter_annotation()
  helices = ss.helices
  helices = split_helices_with_hbond(helices, hierarchy)
  alpha_helices = []
  for helix in helices:
    helix_type = helix.helix_class
    h_bond = helix.get_n_defined_hbonds()
    length = helix.length
    if helix_type == 'alpha' and length > 4:
      alpha_helices.append(helix)
  for i, h in enumerate(alpha_helices):
    h.set_new_serial(serial=i+1, adopt_as_id=True)
  return alpha_helices

def bad_h_bond_index(helix, helix_hierarchy):
  max_h_bond_length = 3.7
  next_i_dict={
    1:4,   # alpha:  O of residue i H-bonds to N of residue i+4
    3:5,   # pi:     O of residue i H-bonds to N of residue i+5
    5:3,   # 3_10:   O of residue i H-bonds to N of residue i+3
   }
  helix_class=helix.get_class_as_int()
  next_i=next_i_dict[helix_class]
  bad_h_bond_index = []
  number_of_poor_h_bonds=0
  asc = helix_hierarchy.atom_selection_cache()
  for i in range(helix.length-next_i):
    cur_residue=get_indexed_residue(
      helix_hierarchy,index=i)
    cur_atom,cur_xyz=get_atom_from_residue(
      residue=cur_residue,
      atom_name=' O  ')
    next_residue=get_indexed_residue(
      helix_hierarchy,index=i+next_i)
    next_atom,next_xyz=get_atom_from_residue(
      residue=next_residue,
      atom_name=' N  ')
    if cur_xyz and next_xyz:  # both present at least
      dd=col(cur_xyz)-col(next_xyz)
      dist=dd.length()
      if dist > max_h_bond_length:
        bad_h_bond_index.append(i)
        print(cur_residue.resname, cur_residue.resseq, dist)
  return bad_h_bond_index

def split_helices_with_hbond(helices, hierarchy):
  asc = hierarchy.atom_selection_cache()
  new_helices = []
  for i,h in enumerate(helices):
    all_length = int(h.length) - 1
    selected_h = hierarchy.select(asc.selection(h.as_atom_selections()[0]))
    i_split_res = bad_h_bond_index(h, selected_h)
    rgs = selected_h.only_model().only_chain().residue_groups()
    for j in i_split_res:
      if j < all_length:
        h1 = h.deep_copy()
        h1.end_resname = rgs[j].atom_groups()[0].resname.strip()
        h1.end_resseq = rgs[j].resseq
        h1.end_icode = rgs[j].icode
        h1.length = int(h1.end_resseq) - int(h1.start_resseq) + 1
        h1.erase_hbond_list()
        h.start_resname = rgs[j+1].atom_groups()[0].resname.strip()
        h.start_resseq = rgs[j+1].resseq
        h.start_icode = rgs[j+1].icode
        h.length = int(h.end_resseq) - int(h.start_resseq) + 1
        h.erase_hbond_list()
        if h1.length > 4:
          new_helices.append(h1)
    if h.length > 4:
      new_helices.append(h)
  if new_helices != []:
    for i, h in enumerate(new_helices):
      h.set_new_serial(serial=i+1, adopt_as_id=True)
    helices = new_helices
  return helices

def write_pdb(file_name):
  save_path = './PDB/training_set/'
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  hierarchy = pdb_inp.construct_hierarchy()
  helix_list = get_helices(pdb_inp, hierarchy)
  alpha_helices = split_helices_with_hbond(helix_list, hierarchy)
  helix_recs = []
  for helix in alpha_helices:
    helix_recs.append(helix.as_pdb_str())
  helix_recs_str = '\n'.join(helix_recs)
  name = file_name[-12:-8]
  save_file = save_path + name + '_new.pdb'
  with open(save_file, "w") as fo:
    fo.write(helix_recs_str)
    fo.write("\n")
    fo.write(hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry()))

def run(args):
  output_path = args[0]
  for root, dirs, files in os.walk(output_path):
    for file in files:
      print(file)
      file_name = root + file
      assert os.path.isfile(file_name)
      write_pdb(file_name)


if __name__ == '__main__':
  run(args=sys.argv[1:])