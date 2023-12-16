# Functions used in Colab Notebook 2 (https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB2_PAH_HF.ipynb)


def FF_optimization(molecule_num:int):
  '''Performs a MMFF94 optimization of a PAH and displays 3D rendering of result. 

    :param int molecule_num: Number corresponding to a PAH from the menu.

    PAH Options:
  
    .. list-table:: 
       :widths: 10 15 5 20 15 15
       :header-rows: 1
    
       * - **Table Index**
         - **PAH**
         - **Rings**
         - **Dimensions (Bohr)**
         - **pi electrons**
         - **Fuse-type**
       * - 1
         - benzene
         - 1
         - 8 x 8 x 3
         - 6
         - linear
       * - 2
         - naphthalene
         - 2
         - 12 x 8 x 3
         - 10
         - linear
       * - 3
         - anthracene
         - 3
         - 16 x 8 x 3
         - 14
         - linear
       * - 4
         - tetracene
         - 4
         - 20 x 8 x 3
         - 18
         - linear
       * - 5
         - pentacene
         - 5
         - 24 x 8 x 3
         - 22
         - linear
       * - 6
         - hexacene
         - 6
         - 28 x 8 x 3
         - 26
         - linear
       * - 7
         - heptacene
         - 7
         - 32 x 8 x 3
         - 30
         - linear
       * - ---
         - ---
         - ---
         - ---
         - ---
         - ---
       * - 8
         - pyrene
         - 4
         - 16 x 12 x 3
         - 16
         - non-linear
       * - 9
         - perylene
         - 5
         - 20 x 12 x 3
         - 20
         - non-linear
       * - 10
         - coronene
         - 7
         - 20 x 16 x 3
         - 24
         - non-linear

    Example:

    >>> # Perform a geometry optimization of heptacene
    >>> FF_optimization(7)

    .. image:: figures/heptacene.png
       :scale: 70 %
       :align: center

  '''
  global molecule_smiles, molecule_smiles_copy # global for the subsequent code block
  molecule_smiles = [] # list of tuples of rdkit object (from Smiles) and molecule name
  cyclobut        = Chem.MolFromSmiles("c1ccc1"), 'cylobut'
  ring1           = Chem.MolFromSmiles("c1ccccc1"), 'Benzene'
  ring2           = Chem.MolFromSmiles("c1c2ccccc2ccc1"), 'Naphthalene'
  ring3           = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1"), 'Anthracene'
  ring4           = Chem.MolFromSmiles("c1c2cc3cc4ccccc4cc3cc2ccc1"), 'Tetracene'
  ring5           = Chem.MolFromSmiles("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1"), 'Pentacene'
  ring6           = Chem.MolFromSmiles("c1c2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2ccc1"), 'Hexacene'
  ring7           = Chem.MolFromSmiles("C1=CC=C2C=C3C=C4C=C5C=C6C=C7C=CC=CC7=CC6=CC5=CC4=CC3=CC2=C1"), 'Heptacene'
  ring4a          = Chem.MolFromSmiles("c1cc2cccc3c2c4c1cccc4cc3"), 'Pyrene'
  ring5a          = Chem.MolFromSmiles("c1ccc5cccc4c5c1c2cccc3cccc4c23"), 'Perylene'
  ring7a          = Chem.MolFromSmiles("c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67"), 'Coronene'
  molecule_smiles.append(cyclobut)
  molecule_smiles.append(ring1)
  molecule_smiles.append(ring2)
  molecule_smiles.append(ring3)
  molecule_smiles.append(ring4)
  molecule_smiles.append(ring5)
  molecule_smiles.append(ring6)
  molecule_smiles.append(ring7)
  molecule_smiles.append(ring4a)
  molecule_smiles.append(ring5a)
  molecule_smiles.append(ring7a)
  molecule_smiles_copy = molecule_smiles.copy()

  molecule = molecule_smiles[molecule_num][0] # user-selected molecule from the list
  global mol, xyz_geom
  mol = Chem.Mol(molecule)
  mol = AllChem.AddHs(mol, addCoords=True)
  AllChem.EmbedMolecule(mol)
  AllChem.MMFFOptimizeMolecule(mol) # uses MMFF94 Force Field to optimize structure

  # format structure into readable xyz coordinates
  xyz_geom = []
  for lines in Chem.MolToXYZBlock(mol).split("\n")[2:]:
      strip = lines.strip()
      if strip:
          xyz_geom.append(strip)
  xyz_geom = "\n".join(xyz_geom)

  # render 3D ball-and-stick model of structure
  display(Markdown('<br>'))
  display(Markdown('               MMFF94 Optimized Geometry'))
  view = py3Dmol.view(
      data=Chem.MolToMolBlock(mol),
      style={"stick": {}, "sphere": {"scale": 0.3}},
      width=400,
      height=300,
  )
  view.zoomTo()
  view.show()

def run_calculation(xyz_struc, basis="STO-3G"):
  '''Performs a HF/STO-3G single-point calculation on a PAH structure. Energy results are displayed in a Markdown table.

    :param str xyz_struc: geometry of PAH in xyz format
    :param str basis: basis set; see options: https://pyscf.org/_modules/pyscf/gto/basis.html
    
    Heptacene Example:
    
    .. image:: figures/heptacene_scf.png

  '''
  clear_output(wait=True)
  display(widgets.VBox([top_panel, output]))
  global mf, mol_quantum, molecule_index, molecule_name
  # create pyscf object from input structure
  mol_quantum    = gto.M(atom=xyz_struc, basis=basis, unit="ANG", symmetry=False)
  mol_quantum.build();
  molecule_index = PAH_d_options[PAH_dropdown.value-1][1]
  molecule_name  = PAH_d_options[PAH_dropdown.value-1][0].split('-')[1]
  print(f"Peforming HF/STO-3G single-point caculation on {molecule_name}...")
  mf = scf.RHF(mol_quantum).run(verbose=0)
  print("SCF Converged!")
  mos        = mf.mo_coeff
  elec_ener  = mf.energy_elec()[0] # electronic energy
  nuc_ener   = mf.energy_nuc()     # nuclear     "
  total_ener = mf.energy_tot()     # total       "
  # Markdown table display:
  display(Markdown('<br>'))
  display(Markdown(f'''Component|Energy (Ha)
  :---:|:---:
  electronic  | {elec_ener:.6f}
  nuclear     | {nuc_ener:.6f}
  **total**   | {total_ener:.6f}
  '''))
  display(Markdown('<br>'))

def get_HL_gap(mf):
  '''Finds HOMO-LUMO gap of molecule from HF/STO-3G calculation result.

    :param pyscf.scf.RHF mf: pyscf HF calculation object
    :return: Float representing HL gap in electron-volts.
  '''
  global homo_num, lumo_num
  homo_num  = np.count_nonzero(mf.mo_occ == 2) - 1 # HOMO index
  lumo_num  = np.count_nonzero(mf.mo_occ == 2)     # LUMO index
  homo_en   = mf.mo_energy[homo_num]               # HOMO energy
  mos_occ   = mf.mo_coeff[:,:lumo_num]
  lumo_en   = mf.mo_energy[lumo_num]               # LUMO energy
  mos_unocc = mf.mo_coeff[:,lumo_num:]
  hl_gap    = abs(lumo_en - homo_en)
  hl_gap_ev = hl_gap * 27.2114                   # 1 Ha = 27.2114 eV conversion
  return hl_gap_ev

def elec_properties(b):
  '''Prints dipole moment, HL gap, and Mulliken charges from HF calculation result.

    :param b: Button click
  '''
  HL_gap = get_HL_gap(mf)
  capture_output = mf.dip_moment(verbose=3)
  print(f"HOMO-LUMO Gap (eV): {HL_gap:.4f}")
  print("===============================")
  capture = mf.mulliken_pop(verbose=3)
  for atom in mol.GetAtoms():
    atom.SetProp('molAtomMapNumber', str(atom.GetIdx()+1))
  display(mol)
  properties_button.disabled = True

def click_scf(b):
  run_calculation(xyz_geom)
  display(Markdown('---'))
  display(properties_button)
  display(Markdown('<br>'))
  scf_button.disabled = True
  properties_button.disabled = False

def mo_diagram(mf, show_diagram, show_table, savefig, filename):
  '''Displays MO energy level diagram and dataframe.

    :param pyscf.scf.RHF mf: pyscf HF calculation object.
    :param bool show_diagram: Display matplotlib energy diagram.
    :param bool show_table: Display same data as Pandas dataframe.
    :param bool savefig: Save figure as .png file.
    :param str filename: filename to save .png file to (w/o file extension).

    Example:

    .. image:: figures/heptacene_diagram.png
       :scale: 70 %

  '''
  # Find HOMO-LUMO gap
  homo_num  = np.count_nonzero(mf.mo_occ == 2) - 1
  lumo_num  = np.count_nonzero(mf.mo_occ == 2)
  homo_en   = mf.mo_energy[homo_num]
  mos_occ   = mf.mo_coeff[:,:lumo_num]
  lumo_en   = mf.mo_energy[lumo_num]
  mos_unocc = mf.mo_coeff[:,lumo_num:]
  hl_gap    = abs(lumo_en - homo_en)
  hl_gap_ev = hl_gap * 27.2114

  # Create a list MO labels (i.e. HOMO, HOMO-1, etc.)
  MO_labels = []
  for i in range(len(mf.mo_occ)-lumo_num):
    if i == 0:
      MO_labels.append('LUMO')
    else:
      MO_labels.append(f'LUMO+{i}')
  MO_labels.reverse()
  for i in range(homo_num+1):
    if i == 0 :
      MO_labels.append('HOMO')
    else:
      MO_labels.append(f'HOMO-{i}')
  MO_labels.reverse()

  table = pd.DataFrame({"En (Ha)": mf.mo_energy, "Occ.": mf.mo_occ.astype('int'), "Type": MO_labels})
  table.index.name = 'MO #'

  global frontier_nums, frontier_labels
  frontier_nums   = np.arange(homo_num-4,lumo_num+5,1)
  frontier_labels = table.iloc[frontier_nums]['Type'].tolist()
  energies        = table.iloc[homo_num-4:lumo_num+5][::-1]['En (Ha)']

  MO_labels2      = []
  # Create a list MO labels (i.e. HOMO, HOMO-1, etc.)
  for i in range(len(mf.mo_occ)-lumo_num):
        MO_labels2.append(f'LUMO+{i}')
  MO_labels2.reverse()
  for i in range(homo_num+1):
        MO_labels2.append(f'HOMO-{i}')
  MO_labels2.reverse()

  # pandas dataframe for the matplotlib diagram
  table2 = pd.DataFrame({"En (Ha)": mf.mo_energy, "Occ.": mf.mo_occ.astype('int'), "Type": MO_labels2})
  datat  = table2.values.tolist()

  y = [] # energy values
  orb = 5
  for i in datat:
    if int(i[2][5:]) < orb:
        y.append(round(i[0],4))
  x = [1]*len(y)
  for i in range(0,len(y)):
        if y[i] in y[:i]:
            count = list(y[:i]).count(y[i])
            x[i] += count*0.3 # x value offset for degenerate orbitals (avoids overlap)

  # Building up the diagram
  fig = plt.figure(figsize=(7, 6))
  ax  = fig.add_subplot(1, 1, 1)
  plt.ylabel("Energy (Ha)", labelpad=7)
  plt.scatter(x[orb:], y[orb:], marker=0, s=1200, linewidths=4, color='#F97306', label='virtual')
  plt.scatter(x[:orb], y[:orb], marker=0, s=1200, linewidths=4, color='green', label='occupied')
  annotations  = []
  annotations2 = []
  for i in range(len(x)):
        if i < orb:
            if i != 0:
                annotations.append('HOMO-' + str(i))
                annotations2.append('LUMO+' + str(i))
            else:
                annotations.append('HOMO')
                annotations2.append('LUMO')
  all_annotations = annotations[::-1]+annotations2
  for i, label in enumerate(all_annotations):
    if list(y).count(y[i])==1:
        plt.annotate(label, (x[i] + 0.005, y[i]-0.01), size=8)
        plt.text(x[i]-0.42, y[i]-0.01, "{:.3f}".format(y[i]), size=8)
    else:
        plt.annotate(label, (x[i] - 0.20, y[i] + 0.015), size=8)
        if y[i] not in y[:i]:
            plt.text(x[i]-0.42, y[i]-0.01, "{:.3f}".format(y[i]), size=8)
  plt.rcParams["legend.markerscale"] = 0.50
  plt.legend(loc='upper right', handletextpad=-0.3)
  plt.title(f'{molecule_name} Frontier Orbitals')
  plt.xticks([])
  plt.xlim([0.3, 3])
  plt.ylim([min(y) - 0.1, max(y) + 0.1])

  # display HOMO-LUMO gap value
  plt.text(0.80, 0.730, f'H-L Gap: {hl_gap:.3f} Ha\n            ({hl_gap_ev:.3f} ev)',
        size=15, horizontalalignment='center', verticalalignment='center',
        transform=ax.transAxes)

  # adding the stick molecule
  mol_copy     = molecule_smiles_copy[molecule_index][0]
  mol_copy_pil = MolToImage(mol_copy)
  imagebox     = OffsetImage(mol_copy_pil, zoom=0.5, interpolation='bicubic')
  imagebox.image.axes = ax
  ab = AnnotationBbox(imagebox, (0.5,0.53),
                    xybox=(120., -80.),
                    xycoords='figure fraction',
                    boxcoords="offset points",
                    frameon=False)
  ax.add_artist(ab)
  plt.tight_layout()

  # save/display control events
  if savefig == True:
    plt.savefig(f'{filename}.png', dpi=800)
  if show_diagram == True:
    plt.show()
  elif show_diagram == False:
    plt.close()
    display(Markdown('          '))
  if show_table == True:
    display(table.iloc[homo_num-4:lumo_num+5][::-1])

def on_click_save(b):
  mo_diagram(mf_mo, show_diagram=False, show_table=False, savefig=True, filename=filename_text_mpl.value)

def get_orbital(i, label):
  '''Saves a .cube file of orbital and prob. density from the above calculation.

    :param int i: index to iterate over in for loop
    :param str label: MO label, i.e. HOMO, HOMO-1, etc
  '''
  tools.cubegen.orbital(mol_quantum, f'{label}.cube', mf.mo_coeff[:,i],nx=80,ny=80,nz=80)
  square_cube = cube_tools.cube(f'{label}.cube')
  pos = np.ma.array(square_cube.data, mask = square_cube.data < 0)
  good_isovalues.append([f'{label}',pos.mean()])
  square_cube.square_cube(power=2)
  pos = np.ma.array(square_cube.data, mask = square_cube.data < 0)
  good_isovalues.append([f'{label}_square',pos.mean()])
  square_cube.write_cube(f'{label}_square.cube')

def draw_orbital(label, psi_square=False, show_noninteractive_png=False):
  '''Renders isosurface of selected orbital from HF/STO-3G calculation.

    :param str label: orbital designation, e.g. HOMO, HOMO-1, etc
    :param bool psi_square: show orbital (False) or probability density (True)
    :param bool show_noninteractive_png: render a static image (True) or interactive 3D display (False)   

    Example:

    >>> # Display probability density of HOMO-3 orbital of heptacene (calculated in a previous cell)
    >>> draw_orbital('HOMO-3', psi_square=True, show_noninteractive_png=False)

    .. image:: figures/heptacene_HOMO3.png
       :scale: 55 %
       :align: center

  '''
  view = py3Dmol.view(width=500,height=500)
  if psi_square:
    index  = int(np.where(good_isovalues == f'{label}_square')[0])
    isoval = float(good_isovalues[index][1])
    print('-----------------------------------------')
    print(f'isovalue: {isoval:.4f}')
    print('rendering...')
    view.addVolumetricData(grid_data[f'{label}_square'], "cube", {'isoval': isoval, 'color': "red", 'opacity': 0.90})
  else:
    with open(f'{label}.cube') as f:
      cube_data = f.read()
    index  = int(np.where(good_isovalues == f'{label}')[0])
    isoval = float(good_isovalues[index][1])*2
    print('----------------------------------------')
    print(f'isovalue: {isoval:.4f}')
    print('rendering...')
    view.addVolumetricData(grid_data[f'{label}'], "cube", {'isoval': isoval, 'color': "red", 'opacity': 0.90})
    view.addVolumetricData(grid_data[f'{label}'], "cube", {'isoval': -isoval, 'color': "blue", 'opacity': 0.90})
  view.addModel(Chem.MolToMolBlock(mol), 'mol')
  view.setStyle({'stick':{}, 'sphere': {"scale":0.3}})
  if show_noninteractive_png:
    view.zoomTo()
    view.show()
    view.png()
  else:
    view.zoomTo()
    view.show()
