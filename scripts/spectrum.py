#! /usr/bin/env python

"""
Script for calculating total flux from outputs produced by Blacklight.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as cons
import astropy.units as u
h=cons.h.cgs.value


def get_flux(**kwargs):

  # Parameters
  c = 2.99792458e10
  gg_msun = 1.32712440018e26
  pc = 9.69394202136e18 / np.pi
  #jy = 1.0e-23
  eV = 1.60218e-12
  data_format = np.float64

  # Prepare metadata
  distance = kwargs['distance']
  mass_msun = kwargs['mass']
  width_rg = kwargs['width']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read metadata
      if mass_msun is None:
        mass_msun = f['mass_msun'][0]
      elif not np.isclose(f['mass_msun'][0], mass_msun):
        raise RuntimeError('Input mass {0} does not match file value {1}.'.format(mass_msun,
            f['mass_msun'][0]))
      if width_rg is None:
        width_rg = f['width'][0]
      elif not np.isclose(f['width'][0], width_rg):
        raise RuntimeError('Input width {0} does not match file value {1}.'.format(width_rg,
            f['width'][0]))
      try:
        f['Q_nu']
        polarization = True
      except KeyError:
        polarization = False
      multiple_frequencies = True if len(f['frequency']) > 1 else False
      frequencies = f['frequency'][:]
      if multiple_frequencies:
        if kwargs['frequency_num'] is None:
          raise RuntimeError('Must specify frequency_num.')
      else:
        if kwargs['frequency_num'] is not None and kwargs['frequency_num'] != 1:
          raise RuntimeError('Only single frequency found in file.')

      # Read root image
      try:
        i_nu = f['I_nu'][:]
      except KeyError:
        raise RuntimeError('No intensity data in file.')
      if polarization:
        if multiple_frequencies:
          i_nu = i_nu[kwargs['frequency_num']-1,...]
          q_nu = f['Q_nu'][kwargs['frequency_num']-1,...]
          u_nu = f['U_nu'][kwargs['frequency_num']-1,...]
          v_nu = f['V_nu'][kwargs['frequency_num']-1,...]
        else:
          q_nu = f['Q_nu'][:]
          u_nu = f['U_nu'][:]
          v_nu = f['V_nu'][:]
        image = np.vstack((i_nu[None,:,:], q_nu[None,:,:], u_nu[None,:,:], v_nu[None,:,:]))
      else:
        #if multiple_frequencies:
        #i_nu = i_nu[kwargs['frequency_num']-1,...]
        image = np.copy(i_nu[:,:,:])

      # Read adaptive image
      if max_level is None:
        max_level = f['adaptive_num_levels'][0]
      else:
        max_level = min(max_level, f['adaptive_num_levels'][0])
      if max_level > 0:
        num_blocks = {}
        block_locs = {}
        image_adaptive = {}
        num_blocks[0] = f['adaptive_num_blocks'][0]
        for level in range(1, max_level + 1):
          num_blocks[level] = f['adaptive_num_blocks'][level]
          block_locs[level] = f['adaptive_block_locs_{0}'.format(level)][:]
          if polarization:
            key_i = 'adaptive_I_nu_{0}'.format(level)
            key_q = 'adaptive_Q_nu_{0}'.format(level)
            key_u = 'adaptive_U_nu_{0}'.format(level)
            key_v = 'adaptive_V_nu_{0}'.format(level)
            if multiple_frequencies:
              i_nu = f[key_i][kwargs['frequency_num']-1,...]
              q_nu = f[key_q][kwargs['frequency_num']-1,...]
              u_nu = f[key_u][kwargs['frequency_num']-1,...]
              v_nu = f[key_v][kwargs['frequency_num']-1,...]
            else:
              i_nu = f[key_i][:]
              q_nu = f[key_q][:]
              u_nu = f[key_u][:]
              v_nu = f[key_v][:]
            image_adaptive[level] = \
                np.vstack((i_nu[None,:,:,:], q_nu[None,:,:,:], u_nu[None,:,:,:], v_nu[None,:,:,:]))
          else:
            key_i = 'adaptive_I_nu_{0}'.format(level)
            
            i_nu = f[key_i][:]
            image_adaptive[level] = np.copy(i_nu[:,:,:,:])

  

  # Calculate image size
  if distance is None:
    raise RuntimeError('Must supply distance.')
  if mass_msun is None:
    raise RuntimeError('Must supply mass.')
  if width_rg is None:
    raise RuntimeError('Must supply width.')
  rg = gg_msun * mass_msun / c ** 2
  width = 2.0 * np.arctan(0.5 * width_rg / (distance))

  # Prepare flag for NaN values
  nan_found = False
  flux = np.array([])
  # Calculate flux without adaptive refinement
  if max_level == 0:
    nan_found = np.any(np.isnan(image))
    for freq in range(len(frequencies)):
      #this is in ergs cm^-2 s^-1 Hz^-1
      #print(image[freq,:,:])
      tempImage = np.copy(image[freq,:,:])
      #I'm artifically setting a floor to the intensity 
      #tempImage[tempImage<1e-22]=1e-22
      print("frequency: "+frequencies[freq].__format__('.3e'))
      print(np.max(tempImage[~np.isnan(tempImage)]))
      print((width**2/eV).__format__('.2e'))
      print(np.average(tempImage[~np.isnan(tempImage)]).__format__('.2e'))
      flux = np.append(flux, (np.nanmean(tempImage) * width ** 2)/eV)
      #print(flux.shape)

  # Calculate flux with adaptive refinement
  else:

    # Disassemble root image into blocks
    block_res = image_adaptive[1].shape[-1]
    num_blocks_root_linear = image.shape[-1] / block_res
    image_adaptive[0] = np.reshape(image,
        (-1, num_blocks_root_linear, block_res, num_blocks_root_linear, block_res))
    image_adaptive[0] = np.swapaxes(image_adaptive[0], 2, 3)
    image_adaptive[0] = np.reshape(image_adaptive[0], (-1, num_blocks[0], block_res, block_res))

    # Describe disassembled image locations
    block_locs[0] = np.empty((num_blocks[0], 2), dtype=int)
    for block in range(num_blocks[0]):
      block_locs[0][block,1] = block % num_blocks_root_linear
      block_locs[0][block,0] = block / num_blocks_root_linear

    # Prepare flags indicating blocks to be counted
    flags = {}
    for level in range(max_level + 1):
      flags[level] = [True] * num_blocks[level]

    # Determine which blocks are masked by those at higher levels
    for level in range(1, max_level + 1):
      for block_fine in range(num_blocks[level]):
        x_loc_fine = block_locs[level][block_fine,1]
        y_loc_fine = block_locs[level][block_fine,0]
        x_loc_coarse = x_loc_fine / 2
        y_loc_coarse = y_loc_fine / 2
        for block_coarse in range(num_blocks[level-1]):
          x_loc_coarse_test = block_locs[level-1][block_coarse,1]
          y_loc_coarse_test = block_locs[level-1][block_coarse,0]
          if x_loc_coarse == x_loc_coarse_test and y_loc_coarse == y_loc_coarse_test:
            flags[level-1][block_coarse] = False

    # Prepare fluxes
    if polarization:
      flux = np.zeros(4)
    else:
      flux = np.zeros(1)

    # Calculate flux from non-masked blocks
    block_width_root = width / num_blocks_root_linear
    for level in range(max_level + 1):
      block_width = block_width_root / 2 ** level
      for block in range(num_blocks[level]):
        if flags[level][block]:
          if np.any(np.isnan(image_adaptive[level][:,block,:,:])):
            nan_found = True
          flux += np.nanmean(image_adaptive[level][:,block,:,:], axis=(1,2)) * block_width ** 2
    flux /= eV
    

  # Report results
  print('')
  if nan_found:
    print('Warning: ignoring NaN')
    print('')
  if multiple_frequencies:
    for freq in range(len(frequencies)):
      print('frequency: {0}'.format(frequencies[freq]))
      print('F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[freq])))
  else:
    if polarization:
      print('I: F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[0])))
      print('Q: F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[1])))
      print('U: F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[2])))
      print('V: F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[3])))
    else:
      print('F_nu = {0} erg cm^-2 s^-1 Hz^-1'.format(repr(flux[0])))
  print('')
  return flux, frequencies


def main(**kwargs):
  flux, frequencies = get_flux(**kwargs)
  plt.plot(frequencies,frequencies*flux)
  #plt.plot(frequencies, (10**10)*frequencies**3, label=r"$(10^{22.5})*\nu^3$")
  #plt.plot(frequencies, frequencies*1e25*np.exp(-frequencies*h/(2e-7)), label=r"$10^{45}*e^{-\nu/(2*10^{-7})}$")
  #plt.plot(frequencies, frequencies*1e25*np.exp(-frequencies*h/(2e-11)), label=r"$10^{45}*e^{-\nu/(2*10^{-11})}$")
  plt.xscale('log')
  plt.ylim(1e5, 1e21)
  plt.yscale('log')
  plt.xlabel('Frequency (Hz)')
  plt.ylabel('$\\nu F_\\nu (erg cm^{-2} s^{-1})$ ')
  plt.title('Flux vs Frequency for file '+kwargs['filename_data'].split('/')[-1])
  plt.legend()
  plt.grid()
  plt.show()


# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename_data', help='name of file containing raw image data')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in gravitational radii')
  parser.add_argument('-m', '--mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('-w', '--width', type=float,
      help='full width of image in gravitational radii')
  parser.add_argument('-f', '--frequency_num', type=int,
      help='index (1-indexed) of frequency to use if multiple present')
  parser.add_argument('-l', '--max_level', type=int,
      help='maximum adaptive level to use in calculation')
  args = parser.parse_args()
  main(**vars(args))
