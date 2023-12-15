# Functions used in Colab Notebook 1 (https://colab.research.google.com/github/tjz21/DFT_PIB_Code/blob/main/notebooks/NB1_3D_PIB.ipynb)

import numpy as np

def psi_reg(x, y, z, q_nx=1, q_ny=1, q_nz=1, lx=1, ly=1, lz=1):
    '''Numeric representation of the normalized 3D PIB wavefunction.

    :param float x: Cartesian spatial x variable.
    :param float y: Idem for y.
    :param float z: Idem for z.
    :param int q_nx: Quantum number specifying the state along x.
    :param int q_ny: Idem for y.
    :param int q_nz: Idem for z.
    :param int lx: Box dimension along x in Bohr.
    :param int ly: Idem for y.
    :param int lz: Idem for z.

    :return: Wavefunction evaluated at a point or array if input x, y, z are arrays.

    Example:

    >>> # value of the wavefunction at the center of a 3x3x3 box
    >>> # in the state 111
    >>> psi_reg(1.5, 1.5, 1.5, 1, 1, 1, 3, 3, 3)
    0.5443310539518174
    '''
    wvfn = np.sqrt(8/(lx*ly*lz)) * \
        np.sin((q_nx*np.pi*x)/lx) * \
        np.sin((q_ny*np.pi*y)/ly) * \
        np.sin((q_nz*np.pi*z)/lz)
    return wvfn

def psi_ener(qnx, qny, qnz, lx, ly, lz):
    '''Calculates energy of 3D PIB state.

    :param int qnx: Quantum number specifying the state along x.
    :param int qny: Idem for y.
    :param int qnz: Idem for z.
    :param int lx: Box dimension along x in Bohr.
    :param int ly: Idem for y.
    :param int lz: Idem for z.

    :return: Float energy value in Ha.

    Example:

    >>> # Energy of 111 state of an 8x8x3 box
    >>> psi_ener(1, 1, 1, 8, 8, 3)
    0.702524
    '''
    e_level = (4*np.pi**2/8)*((qnx/lx)**2 + (qny/ly)**2 + (qnz/lz)**2)
    return np.round(e_level, decimals=6)

def markdown_wvfn(q_nx, q_ny, q_nz, lx, ly, lz):
  '''Displays LaTeX equation of 3D wavefunction.

    :param int q_nx: Quantum number specifying the state along x.
    :param int q_ny: Idem for y.
    :param int q_nz: Idem for z.
    :param int lx: Box dimension along x in Bohr.
    :param int ly: Idem for y.
    :param int lz: Idem for z.

    Example: 

    >>> # Wavefunction of 321 state in 8x8x3 box
    >>> markdown_wvfn(3, 2, 1, 8, 8, 3)

    .. math::

      \\text{State:}\ \ \psi_{3,2,1}(x,y,z)=\sqrt{\\frac{8}{(8)(8)(3)}}\\sin\\bigg(\\frac{3\pi x}{8}\\bigg)\sin\\bigg(\\frac{2\pi y}{8}\\bigg)\sin\\bigg(\\frac{1\pi z}{3}\\bigg)

  '''
  subscript = f'{q_nx},{q_ny},{q_nz}'
  sqrt_denom = f'({lx})({ly})({lz})'
  display(Markdown('## $$ \\text{State: \ \ } \psi_{' + subscript + '}(x,y,z)=' +
                           '\sqrt{\\frac{8}{' + sqrt_denom + '}}' +
                           '\sin{\\Big(\\frac{' + f'{q_nx}' '\pi x}{' + f'{lx}' + '}\\Big)}' +
                           '\sin{\\Big(\\frac{' + f'{q_ny}' '\pi y}{' + f'{ly}' + '}\\Big)}' +
                           '\sin{\\Big(\\frac{' + f'{q_nz}' '\pi z}{' + f'{lz}' + '}\\Big)} \\newline $$'))

def markdown_ener(l_x, l_y, l_z):
  '''Displays LaTeX equation of energy.

    :param int l_x: Box dimension along x in Bohr.
    :param int l_y: Idem for y.
    :param int l_z: Idem for z.

    Example:

    >>> # Energy expression for a 16x8x3 Box
    >>> markdown_ener(16, 8, 3)

    .. math::

      E_{n_x,n_y,n_z}= \\frac{\pi^2}{2}\\bigg[\\bigg(\\frac{n_x}{16}\\bigg)^2 + \\bigg(\\frac{n_y}{8}\\bigg)^2 +\\bigg(\\frac{n_z}{3}\\bigg)^2\\bigg]

  '''
  display(Markdown('## $$ E_{n_x,n_y,n_z} = \\frac{\pi^2}{2} \\Bigl[ ' +
                         ' \\Bigl(\\frac{n_x}{' + f'{l_x}' + '}\\Bigl)^2 +' +
                         ' \\Bigl(\\frac{n_y}{' + f'{l_y}' + '}\\Bigl)^2 +' +
                         ' \\Bigl(\\frac{n_z}{' + f'{l_z}' + '}\\Bigl)^2\\Bigl] \\newline$$'))

def isoplotter(nx_val, ny_val, nz_val, lx, ly, lz, psi_square=False, plot_save=True, plotting_lib='plotly'):
  '''Displays and saves isosurface of 3D PIB wavefunction.

    :param int nx_val: Quantum number specifying the state along x.
    :param int ny_val: Idem for y.
    :param int nz_val: Idem for z.
    :param int lx: Box dimension along x in Bohr.
    :param int ly: Idem for y.
    :param int lz: Idem for z.
    :param bool psi_square: calculate prob. density (true) or wavefunction (false)
    :param bool plot_save: save a .png file of display
    :param str plotting_lib: 3d plotting library to use (plotly or ipyvol)
    
    Example:

    >>> # Render the 321 state of a 8x8x3 PIB system
    >>> isoplotter(3, 2, 1, 8, 8, 3, psi_square=False, plot_save=False, plotting_lib='plotly')
    
    .. image:: figures/321_wavefunction.png

    >>> # Render the 512 probability density of a 16x8x3 PIB system
    >>> isoplotter(5, 1, 2, 16, 8, 3, psi_square=True, plot_save=False, plotting_lib='plotly')
    
    .. image:: figures/512_prob_den.png

  '''
  # construct 3d grid of points
  nx_p, ny_p, nz_p = 7 * lx, 7 * ly, 7 * lz
  xp, yp, zp       = np.linspace(0, lx, nx_p), np.linspace(0, ly, ny_p), np.linspace(0, lz, nz_p)
  X, Y, Z          = np.meshgrid(xp, yp, zp, indexing='ij')
  psi              = psi_reg(X, Y, Z, nx_val, ny_val, nz_val, lx, ly, lz)
  norm_psi         = psi_reg(X, Y, Z, nx_val, ny_val, nz_val, lx, ly, lz)**2

  # ipyvolume potting commands
  if plotting_lib == 'ipyvol':
    ipv.clear()
    fig = ipv.figure(title='PIB',width=500, height=500)
    fig.camera.type = 'OrthographicCamera'
    if psi_square:
      norm_sur = ipv.pylab.plot_isosurface(norm_psi, color='red', level=norm_psi.mean(), controls=True,
                                          description='prob. density')
    else:
        pos_values = np.ma.array(psi, mask = psi < 0.0)
        if nx_val == ny_val == nz_val == 1:
          pos_sur = ipv.pylab.plot_isosurface(psi,color='red', level=np.sqrt(norm_psi.mean()), controls=True,
                                              description='positive')
        else:
          pos_sur = ipv.pylab.plot_isosurface(psi, color='red', level=np.sqrt(norm_psi.mean()), controls=True,
                                              description='positive')
          neg_sur = ipv.pylab.plot_isosurface(psi, color='blue', level=-np.sqrt(norm_psi.mean()), controls=True,
                                              description='negative')
    ipv.style.box_off()
    ipv.squarelim()
    ipv.view(0,-75)
    ipv.xyzlabel('lx', 'ly', 'lz')
    ipv.show()

  # plotly potting commands
  elif plotting_lib == 'plotly':
    global figure
    if psi_square:
      den_isoval = norm_psi.mean()
      figure = go.Figure(data=go.Isosurface(
      x=X.flatten(),
      y=Y.flatten(),
      z=Z.flatten(),
      value=norm_psi.flatten(),
      colorscale='BlueRed',
      isomin=-den_isoval,
      isomax=den_isoval,
      surface_count=2,
      showscale=False,
      caps=dict(x_show=False, y_show=False, z_show=False)
      ))
      figure.update_layout(scene = dict(
                      xaxis_title='x',
                      yaxis_title='y',
                      zaxis_title='z',
                      aspectmode='data'),
                    width=800,
                    height=500,
                    title_text=f'Prob. Den., isovalue = {den_isoval:.4f}')
      figure.show()
    else:
      wfn_isoval = np.sqrt(norm_psi.mean())
      figure = go.Figure(data=go.Isosurface(
      x=X.flatten(),
      y=Y.flatten(),
      z=Z.flatten(),
      value=psi.flatten(),
      colorscale='BlueRed',
      isomin=-wfn_isoval,
      isomax=wfn_isoval,
      surface_count=2,
      showscale=False,
      caps=dict(x_show=False, y_show=False, z_show=False)
      ))
      figure.update_layout(scene = dict(
                      xaxis_title='x',
                      yaxis_title='y',
                      zaxis_title='z',
                      aspectmode='data'),
                      width=800,
                      height=500,
                      title_text=f'Wavefunction, isovalue = {wfn_isoval:.4f}')
      figure.show()

def plot_saver(b):
  if lib_dropdown.value == 'ipyvol':
    ipv.savefig(f'{filename_text_ipv.value}.png', width=1200, height=1200)
  elif lib_dropdown.value == 'plotly':
    #plotly.offline.plot(figure, filename='test.html')
    figure.write_image(f'{filename_text_ipv.value}.png', width=800, height=800)


