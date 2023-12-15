# Functions used in Colab Notebook 3 (https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB3_DFT_PIB.ipynb)

def edge_cleaner(func_3d, nx, ny, nz, num_edges=1):
  '''Sets outermost layer (or two) of grid values of a 3D function to zero.
  
    This is needed for the density because of the discontinuities in the numeric 2nd derivatives at the box edges.

    :param np.array func_3d: Flattened 3D grid of points.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :param int num_edges: Number of border layers to set to zero (either 1 or 2).
    :return: A flattened array of grid points.
  '''
  func_3d = func_3d.reshape(nx, ny, nz)
  if num_edges == 1:
    func_3d[0,:,:]  = 0
    func_3d[-1,:,:] = 0
    func_3d[:,0,:]  = 0
    func_3d[:,-1,:] = 0
    func_3d[:,:,0]  = 0
    func_3d[:,:,-1] = 0
  elif num_edges == 2:
    func_3d[0,:,:]  = 0
    func_3d[-1,:,:] = 0
    func_3d[:,0,:]  = 0
    func_3d[:,-1,:] = 0
    func_3d[:,:,0]  = 0
    func_3d[:,:,-1] = 0
    func_3d[1,:,:]  = 0
    func_3d[-2,:,:] = 0
    func_3d[:,1,:]  = 0
    func_3d[:,-2,:] = 0
    func_3d[:,:,1]  = 0
    func_3d[:,:,-2] = 0
  return func_3d.flatten()

def integ_3d(func_3d, dx, dy, dz):
  '''Integrates a 3D function over all defined space.

    :param np.array func_3d: 3D grid of points
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :return: Integrated 3D function as scalar value.

  '''
  return np.sum(func_3d * dx * dy * dz)

def norm_psi_and_den(e_vecs, occ_states, dx, dy, dz):
  '''Normalizes raw eigenvectors from the solver and finds electron density.

    :param np.array e_vecs: Array of eigenvectors from the eigenvalue solver.
    :param int occ_states: Number of occupied KS states (i.e. # electrons / 2).
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.

    :return:
        - **norm_psi** (*np.array*): Array of normalized eigenvectors.
        - **el_den** (*np.array*): Electron density as array.
  '''
  norm_psi = np.zeros_like(e_vecs)
  el_den   = np.zeros_like(e_vecs[:,0])
  for i in range(e_vecs.shape[1]):
      norm_psi[:,i] = e_vecs[:,i]/np.sqrt(integ_3d(e_vecs[:,i]**2, dx, dy, dz))
  for i in range(occ_states):
      el_den += 2* norm_psi[:,i]**2
  return norm_psi, el_den

def noninter_kin_e(norm_eigenvecs, occ_states, kin_mat, dx, dy, dz, nx, ny, nz):
  '''Finds noninteracting KS kinetic energy.

    :param np.array norm_eigenvecs: array of normalized eigenvectors as columns/rows
    :param int occ_states: number of occupied KS states (i.e. # electrons / 2)
    :param scipy.sparse.sparray kin_mat: kinetic energy operator as sparse matrix
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :param int nx: Number of grid points in each x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return: Scalar energy value in Ha.
  '''
  kin_energy_values = []
  for eig in norm_eigenvecs.T[:occ_states]:
    inner_prod   = eig*kin_mat.dot(eig)
    inner_prod   = edge_cleaner(inner_prod, nx, ny, nz, num_edges=1)
    orbital_k_en = integ_3d(inner_prod, dx, dy, dz)
    kin_energy_values.append(orbital_k_en)
  return sum(kin_energy_values)

def grid_density(l_x, l_y, l_z, plotting_lib='plotly'):
  '''Generates 3D scatter plot representing grid for DFT calculations.

    :param int l_x: Box length in x-direction.
    :param int l_y: Idem for y.
    :param int l_z: Idem for z.

    Example:   

    >>> # View numerical grid for a 16x8x3 Box
    >>> grid_density(16, 8, 3, plotting_lib='plotly')

    .. image:: figures/grid_density.png
       :align: center
       :scale: 60 %

  '''
  nx, ny, nz = (5 * l_x), (5 * l_y), (5 * l_z)
  gpoints    = (nx * ny * nz)
  xp, yp, zp = np.linspace(0, l_x, nx), np.linspace(0, l_y, ny), np.linspace(0, l_z, nz)
  X, Y, Z    = np.meshgrid(xp, yp, zp, indexing='ij')
  label      = '### <center> Grid Points (5 pts/bohr): '
  label     += '$x_p \\times y_p \\times z_p ='
  label     += f'{nx}' + '\\times' + f'{ny}' + '\\times' + f'{nz}' + ' = '+ f'{gpoints}$ <center/>'
  display(Markdown(label))
  if plotting_lib == 'ipyvol':
    fig1 = ipv.figure(title='PIB',width=450, height=450)
    fig1.camera.type = 'OrthographicCamera'
    ipv.quickscatter(X.flatten(), Y.flatten(), Z.flatten(), size=0.5, marker='box',description='grid points')
    ipv.squarelim()
    ipv.style.box_off()
    ipv.show()
  elif plotting_lib == 'plotly':
    scatter = go.Figure(data=[go.Scatter3d(x=X.flatten(), y=Y.flatten(), z=Z.flatten(),
                                       mode='markers',
                                       marker=dict(size=2,color='red',symbol='square',opacity=0.5))])
    scatter.update_layout(scene = dict(xaxis_title='x', yaxis_title='y', zaxis_title='z',
                  aspectmode='data'), width=900, height=500)
    scatter.show()

def exch_equation(den_number):
  '''Displays LaTeX equation of LDA exchange potential.

    :param float den_number: Value of electron density, n(r), at a specific position..

    Example:

    >>> # equation for electronic exchange at a point with a density of 0.50 e/Bohr**3
    >>> exch_equation(0.50)

    .. math::
      
      v_{\\text{X, n(r)=0.50}}^{\\text{LDA}}= - \\frac{3}{4} \left(\\frac{3}{\pi}\\right)^{1/3} (0.50)^{1/3} = -0.586\ \\text{Ha/e}

  '''
  exch_pot_number = np.round(-(3/4)*(3/np.pi)**(1/3)*den_number**(1/3), decimals=3)
  text = Markdown('## $v_{\\text{X},\ n(r)= \ ' + f'{den_number:.2f}' + '}^{\\text{LDA}} = -\\frac{3}{4}' +
                       '\left(\\frac{3}{\pi}\\right)^{1/3}' +
                       f'({den_number:.2f})' + '^{1/3}' +
                       f'={exch_pot_number:.3f}\ ' + '\\text{Ha/e}$')
  display(text)
  clear_output(wait=True)

def exch_pot_eq(den_number):
  '''Displays line plot of LDA exchange with scatterpoint at density value.

    :param float den_number: Value of electron density, n(r).
    
    Example:

    >>> # plot for lda electronic exchange with a point at a density of 0.50 e/Bohr**3
    >>> exch_pot_eq(0.50)

    .. image:: figures/LDA_exch.png

  '''
  exch_pot_number = np.round(-(3/4)*(3/np.pi)**(1/3)*den_number**(1/3), decimals=3)
  display(Markdown('<br>'))
  x = np.linspace(0, 10, 300)
  y = np.round(-(3/4)*(3/np.pi)**(1/3)*x**(1/3), decimals=3)
  fig = plt.figure()
  ax  = fig.add_subplot(1, 1, 1)
  plt.cla()
  plt.clf()
  clear_output(wait=True)
  plt.plot(x, y, label=r'Slater exchange')
  plt.hlines(exch_pot_number, 0, den_number, colors='black')
  plt.vlines(den_number, -2, exch_pot_number, colors='black')
  plt.scatter(den_number, exch_pot_number, marker="s", color='red', label='grid point')
  plt.title('LDA Exchange Potential', size=15)
  plt.xlabel('Density, $n(r)$', size=15)
  plt.ylabel('Potential, $v_{\mathrm{X}}^{\mathrm{LDA}}(n(r))$', size=15)
  plt.xlim(0,2)
  plt.ylim(-1.25,0)
  plt.rcParams["legend.markerscale"] = 1.5
  plt.legend(loc='upper right', fontsize=15)
  plt.show()

def LDA_c_display(ws_radius):
  '''Displays line plot of Chachiyo LDA correlation potential with a scatterpoint and LaTeX equation at the Wigner-Seitz radius (r_s) value.

    :param float ws_radius: Wigner-Seitz radius (inversely related to density).

    Example:

    >>> # display equation and diagram for LDA correlation with point @ r_s = 2.00
    >>> LDA_c_display(2.00)

    .. math::
     
      v_{\\text{C}, r_s=2.00}^{\\text{LDA}}=a\ln\left(1+\\frac{b}{(2.00)}+\\frac{b}{(2.00)^2}\\right) = -0.059\ \\text{Ha/e}

    .. image:: figures/LDA_cor.png
    
  '''
  # a, b, c fitting constants used in paper: 10.1063/1.4958669
  a, b, c = (np.log(2)-1)/(2*np.pi**2), 20.4562557, (4*np.pi/3)**(1/3)
  lda_expression = lambda rad: a*np.log(1 + b*c*1/(rad) + b*(c**2)*1/(rad))
  cor_pot_number = lda_expression(ws_radius)

  # LaTeX Equation
  text = Markdown('## $v_{\\text{C},\ r_s=' + f'{ws_radius:.2f}' + '}^{\\text{LDA}}' + '= a\cdot \\text{ln}\left(1+\\frac{b}{' +
                f'({ws_radius:.2f})' + '} + \\frac{b}{' + f'({ws_radius:.2f})^2' +
                '}\\right)=' + f'{cor_pot_number:.3f}\ ' + '\\text{Ha/e}$')
  display(text)

  # Matplotlib Diagram
  x = np.linspace(0.001, 10, 100)
  y = lda_expression(x)
  plt.plot(x, y, label='Chachiyo Correlation')
  plt.scatter(ws_radius, cor_pot_number, color='red', marker="s", label='grid point')
  plt.hlines(cor_pot_number, 0, ws_radius, colors='black')
  plt.vlines(ws_radius, -2, cor_pot_number, colors='black')
  plt.title('LDA Correlation Potential', size=15)
  plt.xlabel('Wigner-Seitz Radius, $r_s$', size=15)
  plt.ylabel('Potential, $v_{\mathrm{C}}^{\mathrm{LDA}}(r_s)$', size=15)
  plt.xlim(0,10)
  plt.ylim(-0.15,0)
  plt.rcParams["legend.markerscale"] = 1.5
  plt.legend(loc='upper right', fontsize=15)
  plt.show()

def GGAx_pot_eq(den_num, grad_den_num):
  '''Generates line plot of PBE exchange enhancement factor with a LaTeX equation and scatterpoint at the specified density + gradient value.

    :param float den_num: value of the electron density
    :param float grad_den_num: gradient of electron density at same point

    Example:

    >>> # Display GGA exchange equation and diagram for n(r) = 0.15 and grad n(r) = 1.00
    >>> GGAx_pot_eq(0.15, 1.00)

    .. math::

      s=\\frac{|1.00|}{2 \cdot 3^{1/3} \pi^{2/3} (0.150)^{4/3}}=2.028 \ \ \ F_{x,s=2.028} = 1 + \kappa - \\frac{\kappa}{\left(1 + \\frac{\mu (2.028)^2}{\kappa} \\right)} = 1.425

    .. image:: figures/GGA_exch_plot.png

  '''
  kappa, mu = 0.804, 0.21951 # constants used in the paper, 10.1103/PhysRevLett.77.3865
  # reduced density gradient, s:
  calculate_s   = lambda den, grad_den: abs(grad_den)/(2*3**(1/3)*np.pi**(2/3)*den**(4/3))
  # enhancement factor, F_s:
  calculate_F_s = lambda s: 1 + kappa - kappa/(1 + mu*s**2 / kappa)
  s_num   = calculate_s(den_num, grad_den_num)
  F_s_num = calculate_F_s(s_num)

  # LaTeX equation
  equation  = '## $s = \\frac{' + f'|{grad_den_num:.3f}|' + '}{2\cdot3' + '^{1/3}\pi^{2/3}' + f'({den_num:.3f})' + '^{4/3}}=' + f'{s_num:.3f}\ '
  equation += '\ \ \ \ F_{x,\ s=' + f'{s_num:.3f}' + '} = 1 + \kappa - \\frac{\kappa}{(1+\\frac{\mu' + f'({s_num:.3f})^2' + '}{\kappa})}=' + f'{F_s_num:.3f}\ ' + '$'
  equation = Markdown(equation)
  display(equation)
  display(Markdown('<br>'))

  # Matplotlib Diagram
  x_s = np.linspace(0, 10, 300)
  y   = calculate_F_s(x_s)
  plt.plot(x_s, y, label=r'PBE Exch. Factor')
  plt.hlines(F_s_num, 0, s_num, colors='black')
  plt.vlines(s_num, -2, F_s_num, colors='black')
  plt.scatter(s_num, F_s_num, color='red', marker="s", label='grid point')
  plt.title('GGA Enchancement Factor', size=15)
  plt.xlabel('RDG, $s$', size=15)
  plt.ylabel('$F_x(s)$', size=15)
  plt.xlim(0,10)
  plt.ylim(0.95, 2.1)
  plt.legend(loc='upper left', fontsize=15)
  plt.show()

def hartree_plotter(freq, solve_poisson):
  '''Display diagram solving for potential from model density in 1D.

    :param int freq: coefficient for frequency in trig function
    :param bool solve_poisson: display (True) potential from trig functionc

    Example:

    >>> # solve poisson equation for trig function with a frequency of 3
    >>> hartree_plotter(3, True)

    .. math::
    
      n(x) = \sin^2(2\pi(3)x)

      \\nabla^2 V(x) = -4\pi n(x) \\xrightarrow[\\text{linalg.cg}]{\\text{scipy}} V(x)
    
    .. image:: figures/1D_hartree.png

  '''
  # Data for density and Ax=b solver
  nx     = 300
  x      = np.linspace(0, np.pi, nx)
  y      = np.sin(freq*x)**2
  diag1x = np.ones(nx)/(x[1])
  D1x    = sparse.spdiags(np.array([-diag1x, diag1x]), np.array([0,1]), nx, nx)
  diagx  = np.ones(nx)/(x[1]**2)
  D2x    = sparse.spdiags(np.array([diagx, -2*diagx, diagx]), np.array([-1,0,1]), nx, nx)
  T      = -1/2 * D2x
  test   = sparse.linalg.cg(-2*T, -4.*np.pi*y) # solve for V

  # LaTeX equation
  equation1 = '## $n(x) = \sin^2( 2 \pi ' + f'({freq})' + 'x)$'
  equation2 = '## $\\nabla^2V(x)=-4\pi n(x) \\xrightarrow[\\text{linalg.cg}]{\\text{scipy}} V(x)$'
  display(Markdown(equation1))
  display(Markdown(equation2))

  # Matplotlib Diagram
  plt.gcf()
  plt.clf()
  plt.title('1D Hartree Potential', size=15)
  plt.plot(x, y, label=f'n(x) = $\sin^2(2 \pi ({freq}) x)$')
  if solve_poisson:
    plt.plot(x, test[0], label='V(x)')
  plt.xlabel('x', size=15)
  plt.ylabel('Amplitude', size=15)
  plt.legend(loc='upper right')
  plt.grid()
  plt.show()
  print() # a little extra space

def hamiltonian_display(functional, har, ex, cor):
  '''Generates LaTeX equation of effective single-particle Hamiltonian.

    :param str functional: XC potential (LDA or PBE)
    :param bool har: Hartree potential on or off
    :param bool ex: exchange term on or off
    :param bool cor: correlation term on or off

    Example:

    >>> # display effective hamiltonian with exch and hartree and without correlation
    >>> hamiltonian_display('LDA', True, True, False)
    
    .. math::

      \hat{h}_i = \hat{T}_{\\text{kin},i} + v_{\\text{Ha}}(n(r)) + v_{\\text{X}}^{\\text{LDA}}(n(r)) + 0

  '''
  ham = '## $$ \hat{h}_i = \hat{T}_{\\text{kin}, i} +'
  if har == True:
    ham += 'v_{\\text{Ha}}(n(r))'
  else:
    ham += '0'
  if ex == True:
    if functional == 'LDA':
      ham += '+ v_{\\text{X}}^{\\text{LDA}}(n(r))'
    elif functional == 'PBE':
      ham += '+ v_{\\text{X}}^{\\text{PBE}}(n(r), \\nabla n(r))'
  else:
    ham += '+ 0'
  if cor == True:
    if functional == 'LDA':
      ham += '+ v_{\\text{C}}^{\\text{LDA}}(n(r))'
    elif functional == 'PBE':
      ham += '+ v_{\\text{C}}^{\\text{PBE}}(n(r), \\nabla n(r))'
  else:
    ham += '+ 0'
  ham += '$$'
  display(Markdown(ham))

def energy_plot(den_log, ener_log, converge_state, show_fig=True, save_fig=False, filename=None):
  '''Generates iteration number vs total energy plot.

    :param list den_log: list of numpy arrays with electron density
    :param list ener_log: list of total energy values in Ha
    :param bool converge_state: SCF loop ended in a converged (True) or unconverged (False) state
    :param bool show_fig: display figure to output (True) or not (False)
    :param bool save_fig: save a .png file of diagram
    :param str filename: filename to save to (w/o extension)
  '''
  fig = plt.figure(figsize=(6, 4))
  ax = fig.add_subplot(1, 1, 1)
  plt.title('Convergence Plot')
  plt.scatter(0, energy_log[0], color='#F97306', label='noninter Energy')
  plt.plot(np.arange(1,len(energy_log[1:]) + 1), energy_log[1:], 'o-', label='DFT Energy')
  plt.legend(loc='upper right')
  plt.text(0.785, 0.800, f'iterations: {len(density_log) - 1}',
          horizontalalignment='center', verticalalignment='center',
          transform=ax.transAxes)
  plt.text(0.785, 0.730,f'            energy (Ha): {energy_log[-1]:.5f}',
          size=10.5, horizontalalignment='center', verticalalignment='center',
          transform= ax.transAxes)
  if converge_state == True:
    plt.text(0.79, 0.660,f'converged',
            size=10.5, horizontalalignment='center', verticalalignment='center',
            transform= ax.transAxes, weight='bold')
  elif converge_state == False:
    plt.text(0.79, 0.660,f'unconverged',
            size=10.5, horizontalalignment='center', verticalalignment='center',
            transform= ax.transAxes, weight='bold')
  plt.ylabel('Energy (Ha)')
  plt.xlabel('iteration #')
  plt.tight_layout()
  if save_fig:
    plt.savefig(f'{filename}.png', dpi = 800)
  if show_fig:
    plt.show()
  else:
    plt.close()

def hard_walls(potential, nx, ny, nz):
  '''Sets potential of outermost grid points to 1,000. Used for testing.

    :param np.array potential: KS (or other) potential as a flattened array
    :param int nx: Grid points in x-direction.
    :param int ny: Grid points in y-direction.
    :param int nz: Grid points in z-direction.
    :return: Flattened array of grid points.
  '''
  potential         = potential.reshape(nx, ny, nz)
  potential[0,:,:]  = 1000
  potential[-1,:,:] = 1000
  potential[:,0,:]  = 1000
  potential[:,-1,:] = 1000
  potential[:,:,0]  = 1000
  potential[:,:,-1] = 1000
  return potential.flatten()

def hartree(den, kin_oper, dx, dy, dz, nx, ny, nz):
  '''Uses Poisson's equation to find Hartree potential from density.

    :param np.array den: Electron density.
    :param np.array kin_oper: Kinetic energy operator as sparse matrix.
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return:
        -  **v_ha_flat** (*np.array*): Hartree potential as flattened array.
        - **v_ha_ener** (*float*): Hartree energy (in Ha units); needed to compute total energy
  '''
  clean_den = np.ma.array(den, mask= abs(den) < 0.000001) # Mask the low-density points
  clean_den = np.ma.filled(clean_den, fill_value=0.0) # Fill with zeros
  den       = clean_den
  den       = edge_cleaner(den, nx, ny, nz, num_edges=1) # Set edges to zero
  v_ha_flat = sparse.linalg.cg(-2*kin_oper,-4.*np.pi*den)[0]
  v_ha_flat = edge_cleaner(v_ha_flat, nx, ny, nz, num_edges=1)
  v_ha_ener = (1/2)*integ_3d(v_ha_flat*den, dx, dy, dz)
  return v_ha_flat, v_ha_ener

def lda_exchange(den, dx, dy, dz):
  '''Finds LDA exchange potential and energy from density.

    :param np.array den: Electron density.
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :return:
        - **exch_pot_lda** (*np.array*): LDA exchange potential.
        - **exch_ener_lda** (*float*): LDA exchnage energy (in Ha).
  '''
  exch_pot_lda = -(3/4)*(3/np.pi)**(1/3)*(den)**(1/3)
  clean_den = np.ma.array(den, mask= abs(den) < 0.000001)
  clean_den = np.ma.filled(clean_den, fill_value=0.0)
  den       = clean_den
  exch_ener_lda = -(3/4)*(3/np.pi)**(1/3)*integ_3d(den**(4/3), dx, dy, dz)
  return exch_pot_lda, exch_ener_lda

def lda_correlation(den, dx, dy, dz):
  '''Finds LDA correlation potential and energy from density.

    :param np.array den: Electron density.
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :return:
        - **corr_pot** (*np.array*): LDA correlation potential.
        - **corr_en** (*float*): LDA correlation energy (in Ha).
  '''
  # a, b, c fitting constants used in paper: 10.1063/1.4958669
  a, b, c   = (np.log(2)-1)/(2*np.pi**2), 20.4562557, (4*np.pi/3)**(1/3)
  corr_pot  = a*np.log(1 + b*c*den**(1/3) + b*(c**2)*den**(2/3))
  clean_den = np.ma.array(den, mask= abs(den) < 0.000001)
  clean_den = np.ma.filled(clean_den, fill_value=0.0)
  den       = clean_den
  corr_en   = integ_3d(den*corr_pot, dx, dy, dz)
  return corr_pot, corr_en

def RDG(den, der_1st, nx, ny, nz):
  '''Finds dimensionless reduced density gradient (needed for PBE exchange).

    :param np.array den: Electron density.
    :param np.array der_1st: 1st derivative operator as sparse matrix.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return:
      - **RDG** (*np.array*): Reduced density gradient as matrix.
  '''
  clean_den = np.ma.array(den, mask = abs(den) < 0.000001)
  den       = clean_den
  RDG       = (2*3**(1/3)*np.pi**(2/3))**(-1) * abs(der_1st.dot(den)) * den**(-4/3)
  RDG       = edge_cleaner(RDG, nx, ny, nz, num_edges=1) #set the edges to zero
  RDG       = np.ma.filled(RDG, fill_value=0.0) # return zeros where density is very small
  return RDG

def pbe_exchange(den, D1st, dx, dy, dz, nx, ny, nz):
  '''Finds PBE exchange potential and energy from density + gradient.

    :param np.array den: Electron density.
    :param np.array D1st: 1st derivative operator as sparse matrix.
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return:
        - **exch_pot_pbe** (*np.array*): PBE exchange potential.
        - **exch_ener_pbe** (*float*): PBE exchange energy (in Ha).
  '''
  # kappa and mu constants used in the paper, 10.1103/PhysRevLett.77.3865
  kappa, mu = 0.804, 0.2195149727645171
  s         = RDG(den, D1st, nx, ny, nz)
  F_xs      = 1 + kappa - kappa * (1 + mu * s**2 / kappa)**(-1) # exch enhancement factor
  exch_pot_pbe  = F_xs * -(3/4)*(3/np.pi)**(1/3)*((den)**(1/3))
  clean_den     = np.ma.array(den, mask= abs(den) < 0.000001)
  clean_den     = np.ma.filled(clean_den, fill_value=0.0)
  den           = clean_den
  exch_ener_pbe = integ_3d(den*exch_pot_pbe, dx, dy, dz)
  return exch_pot_pbe, exch_ener_pbe

def cor_den_grad(den, der_1st, nx, ny, nz):
  '''Finds dimensionless correlation density gradient (used in PBE correlation).

    :param np.array den: Electron density.
    :param np.array der_1st: 1st derivative operator as sparse matrix.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return:
        - **t** (*np.array*): Correlation density gradient.
  '''
  d_g       = abs(der_1st.dot(den))
  clean_den = np.ma.array(den, mask = abs(den) < 0.000001)
  den       = clean_den
  t         = (d_g*np.pi**(1/6))/(4*3**(1/6)*den**(7/6))
  t         = edge_cleaner(t, nx, ny, nz, num_edges=1)
  t         = np.ma.filled(t, fill_value=0.0)
  return t

def pbe_correlation(den, der_1st, dx, dy, dz, nx, ny, nz):
  '''Finds PBE correlation potential and energy from density + gradient.

    :param np.array den: Electron density.
    :param np.array der_1st: 1st derivative operator as sparse matrix.
    :param float dx: Differential volume element in x-direction.
    :param float dy: Idem for y.
    :param float dz: Idem for z.
    :param int nx: Number of grid points in x-direction.
    :param int ny: Idem for y.
    :param int nz: Idem for z.
    :return:
        - **cor_pot_pbe** (*np.array*): PBE correlation potential.
        - **cor_ener_pbe** (*float*): PBE correlation energy (in Ha).
  '''
  lda_c_pot = lda_correlation(den, dx, dy, dz)[0]
  # beta and gamma constants used in the paper, 10.1103/PhysRevLett.77.3865
  beta, gamma = 0.06672455060314922, 0.031090690869654894
  lda_c_pot = np.ma.array(lda_c_pot, mask = abs(lda_c_pot) < 0.000001)
  A = (beta/gamma)*((np.exp(-lda_c_pot/gamma)-1)**(-1))
  t = cor_den_grad(den, der_1st, nx, ny, nz)
  H = gamma*np.log(1+(beta/gamma)*t**2*((1+A*(t**2))/(1+A*(t**2)+(A**2)*(t**4))))
  cor_pot_pbe = lda_c_pot + H
  cor_ener_pbe = integ_3d(den*cor_pot_pbe, dx, dy, dz)
  return cor_pot_pbe, cor_ener_pbe

