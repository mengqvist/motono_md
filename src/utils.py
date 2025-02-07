
from Bio.PDB import PDBParser, Select, PDBIO
import re
import matplotlib as mpl
import numpy as np



class CleanSelect(Select):
    """
    Class for selecting the parts of the structure that I want to save.
    Subclassing the Select class from Bio.PDB,
    as described about midway down at this page:
    https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    """
    def accept_model(self, model):
        """
        Keep only first model (important for NMR structures)
        """
        if model.get_id() == 0:
            return 1
        return 0

    def accept_chain(self, chain):
        """
        Accept all chains
        """
        return 1

    def accept_residue(self, residue):
        """
        Keep only amino acid residues, this cleans away DNA and RNA
        Remove all HETATOMS
        """
        # the 22 natural aa
        aa = ['ALA', 'ARG',
             'ASN', 'ASP',
             'CYS', 'GLN',
             'GLU', 'GLY',
             'HIS', 'ILE',
             'LEU', 'LYS',
             'MET', 'PHE',
             'PRO', 'PYL',
             'SEC', 'SER',
             'THR', 'TRP',
             'TYR', 'VAL']

        # keep amino acids
        if residue.get_resname() in aa:

            # skip HETATOMS
            hetatm_flag, resseq, icode = residue.get_id()
            if hetatm_flag == " ":
                return 1
        return 0

    def accept_atom(self, atom):
        """
        Remove hydrogens
        Keep only the first atom position for disordered atoms.
        """
        # first check whether it is a hydrogen or not
        hydrogen = re.compile("[123 ]*H.*")
        name = atom.get_id()
        if not hydrogen.match(name):

            # now skip all alternate locations (keep only first instance)
            if (not atom.is_disordered()) or atom.get_altloc() == 'A':
                atom.set_altloc(' ')  # Eliminate alt location ID before output.
                return 1
        return 0
        
        
        
def get_cmap(data, vmin=None, vmax=None, cmap='viridis'):
    """Convert data values to colors
    """
    if vmin is None:
        vmax = np.nanmin(list(data.values()))
    if vmax is None:
        vmax = np.nanmax(list(data.values()))
        
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = mpl.cm.get_cmap(cmap)
    
    return cmap, norm
        
        
def continuous_colors_for_pymol(data, cmap, norm):
    """Convert data values to colors
    """  
    color_dict = {}
    for pos in data.keys():
        val = data[pos]
        color = mpl.colors.rgb2hex(cmap(norm(val))).replace('#', '0x')
        color_dict[pos] = color
        
    return color_dict        
        
        
def make_pymol_coloring_script(color_dict, outfile, bg_color='grey80', display_as='cartoon'):
    """Make coloring script for pymol
    """
    outlines = []
    outlines.append(f'bg_color {bg_color}')
    outlines.append(f'as {display_as}')
    outlines.append('set cartoon_discrete_colors, 1')
    outlines.append('set cartoon_oval_length, 1.0')
    outlines.append('set cartoon_oval_witdth, 0.2')
    outlines.append('hide (hydro)')
    
    for pos in color_dict.keys():
        outlines.append('color {}, resi {}'.format(color_dict[pos], pos))
        
    with open(outfile, 'w') as f:
        f.write('\n'.join(outlines))
    
