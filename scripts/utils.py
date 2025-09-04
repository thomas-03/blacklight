import numpy as np

# Parameters
c = 2.99792458e10
gg_msun = 1.32712440018e26
pc = 9.69394202136e18 / np.pi
muas = np.pi / 180.0 / 60.0 / 60.0 / 1.0e6
data_format = np.float64

# Plotting parameters
cmap_i = 'inferno'
cmap_q = 'PuOr'
cmap_u = 'PiYG'
cmap_v = 'RdBu_r'
cmap_level = 'plasma'
cmap_crossings = 'plasma'
cmap_name = 'viridis'
nan_color = 'gray'
interpolation = 'none'
labelpad = 5
dpi = 300

# Import plotting module
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def read_data(**kwargs):

  # Verify input
  if sum((kwargs['stokes_q'], kwargs['stokes_u'], kwargs['stokes_v'], kwargs['name'] is not None,
      kwargs['refinement_level'])) > 1:
    raise RuntimeError('Can have at most one of Stokes Q/U/V, named quantity, or refinement level' \
        + ' selected.')

  # Prepare metadata
  width_rg = kwargs['width']
  mass_msun = kwargs['mass']
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read metadata
      if width_rg is None:
        width_rg = f['width'][0]
      elif not np.isclose(f['width'][0], width_rg):
        raise RuntimeError('Input width {0} does not match file value {1}.'.format(width_rg,
            f['width'][0]))
      if mass_msun is None:
        mass_msun = f['mass_msun'][0]
      elif not np.isclose(f['mass_msun'][0], mass_msun):
        raise RuntimeError('Input mass {0} does not match file value {1}.'.format(mass_msun,
            f['mass_msun'][0]))

      # Read root image
      if kwargs['refinement_level']:
        try:
          image = np.zeros_like(f['I_nu'][:])
        except KeyError:
          raise RuntimeError('No intensity data in file, so image cannot be refined.')
      elif kwargs['name'] is not None:
        try:
          image = f[kwargs['name']][:]
        except:
          raise RuntimeError('{0} not found in file.'.format(kwargs['name']))
      elif kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
        try:
          if kwargs['stokes_q']:
            image = f['Q_nu'][:]
          elif kwargs['stokes_u']:
            image = f['U_nu'][:]
          else:
            image = f['V_nu'][:]
        except KeyError:
          raise RuntimeError('No polarization data in file.')
      else:
        try:
          image = f['I_nu'][:]
        except KeyError:
          raise RuntimeError('No intensity data in file.')

      # Read adaptive image
      if max_level is None:
        max_level = f['adaptive_num_levels'][0]
      else:
        max_level = min(max_level, f['adaptive_num_levels'][0])
      
      num_blocks = {}
      block_locs = {}
      if max_level > 0:
        num_blocks = {}
        block_locs = {}
        image_adaptive = {}
        num_blocks[0] = f['adaptive_num_blocks'][0]
        for level in range(1, max_level + 1):
          num_blocks[level] = f['adaptive_num_blocks'][level]
          block_locs[level] = f['adaptive_block_locs_{0}'.format(level)][:]
          if kwargs['refinement_level']:
            image_adaptive[level] = \
                level * np.ones_like(f['adaptive_I_nu_{0}'.format(level)])
          elif kwargs['name'] is not None:
            image_adaptive[level] = f['adaptive_{0}_{1}'.format(kwargs['name'], level)][:]
          elif kwargs['stokes_q']:
            image_adaptive[level] = f['adaptive_Q_nu_{0}'.format(level)][:]
          elif kwargs['stokes_u']:
            image_adaptive[level] = f['adaptive_U_nu_{0}'.format(level)][:]
          elif kwargs['stokes_v']:
            image_adaptive[level] = f['adaptive_V_nu_{0}'.format(level)][:]
          else:
            image_adaptive[level] = f['adaptive_I_nu_{0}'.format(level)][:]
        return max_level, image_adaptive

      # Select frequency
      if kwargs['frequency_num'] is not None:
        if kwargs['name'] is not None and kwargs['name'] == 'time':
          raise RuntimeError('Cannot specify frequency for plot of ray times.')
        elif kwargs['name'] is not None and kwargs['name'] == 'length':
          raise RuntimeError('Cannot specify frequency for plot of ray lengths.')
        elif kwargs['name'] is not None and kwargs['name'] == 'crossings':
          raise RuntimeError('Cannot specify frequency for plot of plane crossings.')
        if len(image.shape) == 2 and kwargs['frequency_num'] != 1:
          raise RuntimeError('Only single frequency found in file.')
        if len(image.shape) == 3:
          image = image[kwargs['frequency_num']-1,...]
          for level in range(1, max_level + 1):
            image_adaptive[level] = image_adaptive[level][kwargs['frequency_num']-1,...]
      elif len(image.shape) != 2:
        raise RuntimeError('Must specify frequency_num.')

  # Read image data from .npy file
  elif kwargs['filename_data'][-4:] == '.npy':
    if kwargs['refinement_level']:
      raise RuntimeError('Adaptive refinement not supported by npy format.')
    if kwargs['frequency_num'] is not None and kwargs['frequency_num'] != 1:
      raise RuntimeError('Multiple frequencies not supported by npy format.')
    if kwargs['name'] is not None:
      raise RuntimeError('Arbitrary named quantities not supported by npy format.')
    image = np.load(kwargs['filename_data'])
    polarization = image.shape[0] == 4
    if not polarization and (kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']):
      raise RuntimeError('No polarization data in file.')
    if kwargs['stokes_q']:
      image = image[1,...]
    elif kwargs['stokes_u']:
      image = image[2,...]
    elif kwargs['stokes_v']:
      image = image[3,...]
    else:
      image = image[0,...]
    max_level = 0

  # Read image data from raw file
  else:
    if kwargs['refinement_level']:
      raise RuntimeError('Adaptive refinement not supported by raw format.')
    if kwargs['frequency_num'] is not None and kwargs['frequency_num'] != 1:
      raise RuntimeError('Multiple frequencies not supported by raw format.')
    if kwargs['name'] is not None:
      raise RuntimeError('Arbitrary named quantities not supported by raw format.')
    if kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
      raise RuntimeError('Polarization not supported by raw format.')
    image = np.fromfile(kwargs['filename_data'], dtype=data_format)
    pix_num = len(image)
    pix_res = int(pix_num ** 0.5)
    if pix_res ** 2 == pix_num:
      pass
    elif (pix_res + 1) ** 2 == pix_num:
      pix_res += 1
    else:
      raise RuntimeError('Image data not square.')
    image = np.reshape(image, (pix_res, pix_res))
    max_level = 0
  return max_level, num_blocks, block_locs, image

def get_grid(**kwargs):
  # Prepare metadata
  width_rg = kwargs['width']
  mass_msun = kwargs['mass']
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  max_level, num_blocks, block_locs, image = read_data(**kwargs)


  # Calculate root grid in pixels
  if kwargs['axes'] == 'pixel' or (kwargs['axes'] is None and width_rg is None):
    extent = np.array((-0.5, image.shape[0] - 0.5, -0.5, image.shape[1] - 0.5))
    x_label = r'$x$-pixel'
    y_label = r'$y$-pixel'

  # Calculate root grid in gravitational radii
  elif kwargs['axes'] == 'rg' or (kwargs['axes'] is None and mass_msun is None):
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    half_width = 0.5 * width_rg
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($GM/c^2$)'
      y_label = r'$y$ ($GM/c^2$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ GM/c^2$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ GM/c^2$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate root grid in cm
  elif kwargs['axes'] == 'cm' or (kwargs['axes'] is None and distance_pc is None):
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
    rg = gg_msun * mass_msun / c ** 2
    half_width = 0.5 * width_rg * rg
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($\mathrm{cm}$)'
      y_label = r'$y$ ($\mathrm{cm}$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{cm}$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{cm}$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate root grid in muas
  else:
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
    if distance_pc is None:
      raise RuntimeError('Must supply distance.')
    rg = gg_msun * mass_msun / c ** 2
    half_width = np.arctan(0.5 * width_rg * rg / (distance_pc * pc)) / muas
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($\mathrm{\mu as}$)'
      y_label = r'$y$ ($\mathrm{\mu as}$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{\mu as}$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{\mu as}$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate adaptive grid
  if max_level > 0:
    num_blocks_root_linear = image.shape[-1] / image[1].shape[-1]
    block_width = (extent[1] - extent[0]) / num_blocks_root_linear
    extent_adaptive = {}
    for level in range(1, max_level + 1):
      block_width_level = block_width / 2 ** level
      extent_adaptive[level] = np.empty((num_blocks[level], 4))
      for block in range(num_blocks[level]):
        x_loc = block_locs[level][block,1]
        y_loc = block_locs[level][block,0]
        extent_adaptive[level][block,0] = extent[0] + x_loc * block_width_level
        extent_adaptive[level][block,1] = extent[0] + (x_loc + 1) * block_width_level
        extent_adaptive[level][block,2] = extent[2] + y_loc * block_width_level
        extent_adaptive[level][block,3] = extent[2] + (y_loc + 1) * block_width_level
    return extent_adaptive
  return extent,x_label,y_label

def plot_data(**kwargs):
    max_level, num_blocks, block_locs, image = read_data(**kwargs)
    extent, x_label, y_label = get_grid(**kwargs)

    # Calculate refinement parameters
    if kwargs['refinement_level']:
        label = 'level'
        min_level = 0.0
        vmin = 0.0 if kwargs['vmin'] is None else round(kwargs['vmin'])
        vmax = float(max_level) if kwargs['vmax'] is None else round(kwargs['vmax'])
        num_ticks = int(vmax) - int(vmin) + 1
        tick_locs = np.linspace(vmin, vmax, num_ticks)
        vmin_eff = vmin if int(vmin) - int(min_level) <= 0 else vmin - 1.0
        vmax_eff = vmax if int(vmax) - int(max_level) >= 0 else vmax + 1.0
        cmap_chosen = cmap_level if kwargs['cmap'] is None else kwargs['cmap']
        cmap_continuous = plt.get_cmap(cmap_chosen)
        cmap_vals = [cmap_continuous((x - vmin_eff) / (vmax_eff - vmin_eff)) for x in tick_locs]
        cmap = LinearSegmentedColormap.from_list('levels', cmap_vals, num_ticks)
        cmap.set_under(cmap_continuous(0.0))
        cmap.set_over(cmap_continuous(1.0))
        vmin -= 0.5
        vmax += 0.5
        bounds = np.linspace(vmin, vmax, num_ticks + 1)

    # Calculate plane crossing parameters
    elif kwargs['name'] == 'crossings':
        label = 'crossings'
        min_crossings = np.min(image)
        max_crossings = np.max(image)
        for level in range(1, max_level + 1):
            for block in range(num_blocks[level]):
                max_crossings = max(max_crossings, np.max(image[level][block,...]))
                min_crossings = min(min_crossings, np.min(image[level][block,...]))
        vmin = 0.0 if kwargs['vmin'] is None else round(kwargs['vmin'])
        vmax = float(max_crossings) if kwargs['vmax'] is None else round(kwargs['vmax'])
        num_ticks = int(vmax) - int(vmin) + 1
        tick_locs = np.linspace(vmin, vmax, num_ticks)
        vmin_eff = vmin if int(vmin) - int(min_crossings) <= 0 else vmin - 1.0
        vmax_eff = vmax if int(vmax) - int(max_crossings) >= 0 else vmax + 1.0
        cmap_chosen = cmap_crossings if kwargs['cmap'] is None else kwargs['cmap']
        cmap_continuous = plt.get_cmap(cmap_chosen)
        cmap_vals = [cmap_continuous((x - vmin_eff) / (vmax_eff - vmin_eff)) for x in tick_locs]
        cmap = LinearSegmentedColormap.from_list('crossings', cmap_vals, num_ticks)
        cmap.set_under(cmap_continuous(0.0))
        cmap.set_over(cmap_continuous(1.0))
        vmin -= 0.5
        vmax += 0.5
        bounds = np.linspace(vmin, vmax, num_ticks + 1)

  # Calculate named quantity parameters
    elif kwargs['name'] is not None:
        label_mid = {'time': r'\mathrm{min}(t)', 'length': r'l', 'lambda': r'\lambda',
            'emission': r'\int j_\nu \nu^{-2}\ \mathrm{d}\lambda', 'tau': r'\tau', 'rho': r'\rho',
            'n_e': r'n_\mathrm{e}', 'p_gas': r'p_\mathrm{gas}', 'Theta_e': r'\Theta_\mathrm{e}',
            'B': r'B', 'sigma': r'\sigma', 'beta_inverse': r'\beta^{-1}'}
        label_unit = {'time': r' ($\mathrm{s}$)', 'length': r' ($\mathrm{cm}$)',
            'lambda': r' ($\mathrm{s\ cm}$)',
            'emission': r' ($\mathrm{erg\ s^2\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)', 'tau': '',
            'rho': r' ($\mathrm{g\ cm^{-3}}$)', 'n_e': r' ($\mathrm{cm^{-3}}$)',
            'p_gas': r' ($\mathrm{dyne\ cm^{-2}}$)', 'Theta_e': '', 'B': r' ($\mathrm{gauss}$)',
            'sigma': '', 'beta_inverse': ''}
        try:
            if kwargs['name'] in ('time', 'length', 'lambda', 'tau', 'emission'):
                key_mid = kwargs['name']
                label = r'$' + label_mid[kwargs['name']] + r'$'
            elif kwargs['name'][:11] == 'lambda_ave_':
                key_mid = kwargs['name'][11:]
                label = r'$\langle ' + label_mid[key_mid] + r' \rangle_\lambda$'
            elif kwargs['name'][:13] == 'emission_ave_':
                key_mid = kwargs['name'][13:]
                label = r'$\langle ' + label_mid[key_mid] + r' \rangle_{j_\nu / \nu^2}$'
            elif kwargs['name'][:8] == 'tau_int_':
                key_mid = kwargs['name'][8:]
                label = r'$\langle ' + label_mid[key_mid] + r' \rangle_\tau$'
            else:
                raise KeyError
            label += label_unit[key_mid]
        except KeyError:
            underscore = '_' if kwargs['notex'] else r'\_'
            label = ''.join([c if c != '_' else underscore for c in kwargs['name']])
        vmin = np.nanmin(image)
        vmax = np.nanmax(image)
        for level in range(1, max_level + 1):
            vmin_adaptive = np.nanmin(image[level])
            vmax_adaptive = np.nanmax(image[level])
            vmin = np.nanmin((vmin, vmin_adaptive))
            vmax = np.nanmax((vmax, vmax_adaptive))
        vmin = vmin if kwargs['vmin'] is None else kwargs['vmin']
        vmax = vmax if kwargs['vmax'] is None else kwargs['vmax']
        tick_locs = None
        cmap_chosen = cmap_name if kwargs['cmap'] is None else kwargs['cmap']
        cmap = plt.get_cmap(cmap_chosen)
        cmap.set_bad(nan_color)
        bounds = None

        # Calculate Stokes parameters
    else:
        vmax = np.nanmax(np.abs(image))
        for level in range(1, max_level + 1):
            vmax_adaptive = np.nanmax(np.abs(image[level]))
            vmax = np.nanmax((vmax, vmax_adaptive))
        if kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
            vmin = -vmax
        else:
            vmin = 0.0
        vmin = vmin if kwargs['vmin'] is None else kwargs['vmin']
        vmax = vmax if kwargs['vmax'] is None else kwargs['vmax']
        scale_exponent = int('{0:24.16e}'.format(vmax).split('e')[1])
        scale = 10.0 ** scale_exponent
        vmin /= scale
        vmax /= scale
        image /= scale
        for level in range(1, max_level + 1):
            image[level] /= scale
        if kwargs['stokes_q']:
            stokes_str = 'Q'
            cmap_chosen = cmap_q if kwargs['cmap'] is None else kwargs['cmap']
            cmap = plt.get_cmap(cmap_chosen)
        elif kwargs['stokes_u']:
            stokes_str = 'U'
            cmap_chosen = cmap_u if kwargs['cmap'] is None else kwargs['cmap']
            cmap = plt.get_cmap(cmap_chosen)
        elif kwargs['stokes_v']:
            stokes_str = 'V'
            cmap_chosen = cmap_v if kwargs['cmap'] is None else kwargs['cmap']
            cmap = plt.get_cmap(cmap_chosen)
        else:
            stokes_str = 'I'
            cmap_chosen = cmap_i if kwargs['cmap'] is None else kwargs['cmap']
            cmap = plt.get_cmap(cmap_chosen)
        if scale_exponent == 0:
            label = r'$' + stokes_str + r'_\nu$ ($\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
        elif scale_exponent == 1:
            label = r'$' + stokes_str + r'_\nu$ ($10\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
        else:
            label = r'$' + stokes_str + r'_\nu$ ($10^{' + repr(scale_exponent) \
            + r'}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
        tick_locs = None
        cmap.set_bad(nan_color)
        bounds = None
    
    fig, ax = plt.subplots(1, 1, figsize=kwargs['figsize'], dpi=dpi,facecolor='white')
    # Plot root image
    im = ax.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal', origin='lower', extent=extent,
      interpolation=interpolation)

    # Plot adaptive image
    for level in range(1, max_level + 1):
        for block in range(num_blocks[level]):
            im = ax.imshow(image[level][block,...], cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal',
            origin='lower', extent=extent[level][block,:], interpolation=interpolation)

    if kwargs['fontsize'] is None:
       fontsize = 8
    else:
       fontsize = kwargs['fontsize']

    # Make colorbar
    cb = fig.colorbar(im, ax=ax, ticks=tick_locs, boundaries=bounds,fraction=0.05)
    cb.set_label(label, labelpad=labelpad, fontsize=fontsize)
    cb.ax.tick_params(size=fontsize/2,labelsize=fontsize)
    if kwargs['refinement_level'] or kwargs['name'] == 'crossings':
        cb.ax.tick_params(length=0,labelsize=fontsize)

    # Adjust axes
    ax.set_xlim(extent[0], extent[1])
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylim(extent[2], extent[3])
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)

    # Save figure
    if 'filename_plot' in kwargs and kwargs['filename_plot'] is not None:
        plt.savefig(kwargs['filename_plot'], dpi=dpi)
    else:
       return fig,ax
