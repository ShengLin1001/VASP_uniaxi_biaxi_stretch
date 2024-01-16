from numpy import sqrt
from numpy import array, dot
from os    import makedirs, path

from ase import Atoms
from ase.build import bulk, fcc111, surface
from ase.io import write
from ase.visualize import view

###################################################################################
################################# My function #####################################
###################################################################################

# move atoms
def move_atoms(self, translate_matrix = array([0.1, 0.1, 0.0])):
    my_cell = array(self.get_cell())
    self.positions[:,:] += dot(translate_matrix, my_cell)
    self.wrap()
    print('complete move atoms')
    
# stretch unit cell
def uni_stretch(self, stretch_factor = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003], base_filename = 'POSCAR_uniaxial', special_filename = '', my_path = './dump/'):
    self_copy = self.copy()
    my_cell = array(self_copy.get_cell())
    for i in stretch_factor:
        temp = my_cell.copy()
        temp[:2,0] = my_cell[:2,0] * i
        self_copy.set_cell(temp, scale_atoms=True)
        formatted_i = f'{i:.3f}'
        filename = f'{my_path}{base_filename}_{special_filename}_{formatted_i}'
        print(path.dirname(filename))
        makedirs(path.dirname(filename), exist_ok=True)
        write(filename, self_copy, format='vasp', direct=True, vasp5=True)
    print('complete uniaxial stretch')

def bi_stretch(self, stretch_factor = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003], base_filename = 'POSCAR_biaxial', special_filename = '', my_path = './dump/'):
    self_copy = self.copy()
    my_cell = array(self_copy.get_cell())
    for i in stretch_factor:
        temp = my_cell.copy()
        temp[:2,:2] = my_cell[:2,:2] * i
        self_copy.set_cell(temp, scale_atoms=True)
        formatted_i = f'{i:.3f}'
        filename = f'{my_path}{base_filename}_{special_filename}_{formatted_i}'
        print(path.dirname(filename))
        makedirs(path.dirname(filename), exist_ok=True)
        write(filename, self_copy, format='vasp', direct=True, vasp5=True)
    print('complete biaxial stretch')


# generate the POSCAR
def generate_hcp_film(symbol = 'Au', structure = 'hcp', my_a = 2.95 , my_covera = sqrt(8.0/3.0), num_layers = 12, my_vacuum = 20, slice_plane = (0,0,1), my_periodic = True ):
    # Create an bulk hcp structure for gold
    bulk_hcp = bulk('Au', 'hcp', a = my_a, covera = my_covera)
    # slice, Set the number of layers to 12
    num_rep_z = int(num_layers/len(bulk_hcp.numbers))
    slab_hcp = surface(bulk_hcp, slice_plane , num_rep_z, vacuum = my_vacuum, tol=1e-6, periodic=my_periodic)
    slab_hcp.wrap()
    print('complete generate hcp film')
    return slab_hcp

def generate_fcc_film(symbol = 'Au', structure = 'fcc', a_fcc = 2.95*sqrt(2.0), num_layers = 12, my_vacuum = 20, slice_plane = (1,1,1), my_periodic = True ):
    bulk_fcc = bulk('Au', 'fcc', a=a_fcc)
    # slice, Set the number of layers to 12
    num_rep_z = int(num_layers/len(bulk_fcc.numbers))
    slab_fcc = surface(bulk_fcc, slice_plane, num_rep_z, vacuum=my_vacuum, tol=1e-6, periodic=my_periodic)
    slab_fcc.wrap()
    print('complete generate fcc film')
    return slab_fcc

###################################################################################
#################################### find path     ################################
###################################################################################

# 获取当前文件的完整路径
current_file_path = __file__

# 获取当前文件所在的目录
current_directory = path.dirname(current_file_path)

print("当前文件的完整路径:", current_file_path)
print("当前文件所在的目录:", current_directory)

###################################################################################
#################################### HCP structure ################################
###################################################################################

# Define the lattice parameters for hcp and fcc structures
a_hcp = 2.95  # lattice constant
num_layers = 12

gold_hcp = generate_hcp_film('Au', 'hcp', a_hcp, num_layers = num_layers, my_vacuum = 20, slice_plane = (0,0,1), my_periodic = True)

# move atoms along any special directions in scale units
move_atoms(gold_hcp)

# uniaxial and biaxial stretch
uni_stretch(gold_hcp, special_filename='gold_hcp', my_path = './dump/hcp/')
bi_stretch(gold_hcp, special_filename='gold_hcp', my_path = './dump/hcp/') 

###################################################################################
#################################### FCC structure ################################
###################################################################################

a_fcc = 2.95*sqrt(2.0)

gold_fcc = generate_fcc_film('Au', 'fcc', a_fcc=a_fcc, num_layers = num_layers, my_vacuum = 20, slice_plane = (1,1,1), my_periodic = True )

# move atoms along any special directions in scale units
move_atoms(gold_fcc)

# uniaxial and biaxial stretch
uni_stretch(gold_fcc, special_filename='gold_fcc', my_path = './dump/fcc/')
bi_stretch(gold_fcc, special_filename='gold_fcc', my_path = './dump/fcc/') 

