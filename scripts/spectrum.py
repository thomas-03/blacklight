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


h_ev = 4.135667662e-15
eV = 1.60218e-12
h_erg = 6.626e-27
kB=cons.k_B.cgs.value
c = 2.99792458e10
gg_msun = 1.32712440018e26

def get_flux(**kwargs):

  # Parameters
  pc = 9.69394202136e18 / np.pi
  h = cons.h.cgs.value
  data_format = np.float64

  # Prepare metadata
  distance = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:
      mass_msun = f['mass_msun']
      width_rg = f['width']

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
      freq_num = len(frequencies)

      # Read root image
      try:
        i_nu = f['I_nu'][:]
      except KeyError:
        raise RuntimeError('No intensity data in file.')
      if polarization:
        q_nu = f['Q_nu'][:]
        u_nu = f['U_nu'][:]
        v_nu = f['V_nu'][:]
        image = np.vstack((i_nu[None,:,:], q_nu[None,:,:], u_nu[None,:,:], v_nu[None,:,:]))
      else:
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
  width =  2*np.arctan(0.5*width_rg /(distance))
  #width = width_rg/distance

  # Prepare flag for NaN values
  nan_found = False
  flux = np.array([])
  # Calculate flux without adaptive refinement
  if max_level == 0:
    nan_found = np.any(np.isnan(image))
    for freq in range(len(frequencies)):
      #this is in erg cm^-2 s^-1 Hz^-1
      tempImage = np.copy(image[freq,:,:])
      #print((width))
      flux = np.append(flux, (np.nanmean(tempImage[np.isfinite(tempImage)])*width**2))
      #print(flux.shape)
    
  #flux/=eV
  '''
  # Report results
  print('')
  if nan_found:
    print('Warning: ignoring NaN')
    print('')
  if multiple_frequencies:
    for freq in range(len(frequencies)):
      print('frequency: {0}'.format(frequencies[freq]))
      print('F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[freq])))
  else:
    if polarization:
      print('I: F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[0])))
      print('Q: F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[1])))
      print('U: F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[2])))
      print('V: F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[3])))
    else:
      print('F_nu = {0} eV cm^-2 s^-1 Hz^-1'.format(repr(flux[0])))
  print('')
  '''
  return flux, frequencies

def get_luminosity(**kwargs):
  flux, freqs = get_flux(**kwargs)
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:
      mass_msun = f['mass_msun']
  rg = gg_msun * mass_msun / c ** 2
  dA = 4*np.pi*(kwargs['distance']*rg)**2
  return flux*dA, freqs


def main(**kwargs):

  c = 2.99792458e10
  gg_msun = 1.32712440018e26

  if not kwargs['multiInc']:
    plt.figure(figsize=(8,6))
    if kwargs['luminosity']:
      files = kwargs['filename_data']
      if kwargs['files'] is not None:
        files.extend(kwargs['files'])
      for file in files:
        kwargs['filename_data'] = file
        lum,frequencies = get_luminosity(**kwargs)
        print(f"file:{file} luminosity:{lum}")
        if kwargs['labels'] is not None:
          plt.plot(frequencies*h_ev,frequencies*lum,label=kwargs['labels'][files.index(file)])
        else:
          plt.plot(frequencies*h_ev,frequencies*lum/(2*np.pi),label='Inclination {0} deg'.format(kwargs['inclination'][0]))
      #make it so that I can add in the line whether or not we want to compare against something else!!!
      #if kwargs['compare']:
      if kwargs['compare_file'] is not None:
        shaneResults = np.loadtxt(kwargs['compare_file'])
        plt.errorbar(shaneResults[:,0]*1e3,shaneResults[:,1],yerr=shaneResults[:,2],label='MC')
      
      if kwargs['compare_file2'] is not None:
        shaneResults2 = np.loadtxt(kwargs['compare_file2'])
        plt.errorbar(shaneResults2[:,0]*1e3,shaneResults2[:,1],yerr=shaneResults2[:,2],label='MC w/o screen')
      
      if kwargs['filename_data'][-4:] == '.npz':
        with np.load(kwargs['filename_data']) as f:
          mass_msun = f['mass_msun']
      rg = gg_msun * mass_msun / c ** 2
      #1e11 cm . 2rg = 5908253110111
      #B_nu = 2*h_erg*frequencies**3/c**2/(np.exp(h_erg*frequencies/(kB*1e5))-1)
      #get rid of the 2 because imu=sum so this basically returns flux
      #B_nu = h_erg/c**2*frequencies**3/(np.exp(h_erg*frequencies/(kB*1e5))-1)
      #plt.plot(frequencies*h_ev,frequencies*B_nu*4*np.pi*(1e11)**2,label='Blackbody at 10^5 K')
      #plt.errorbar(shaneResults[:,0]*1e3,shaneResults[:,1],yerr=shaneResults[:,2],label='MC Results')
      plt.xscale('log')
      plt.yscale('log')
      plt.xlabel('Frequency (eV)')
      plt.ylabel('$\\nu L_\\nu (erg s^{-1})$ ')
      #plt.title('Spectrum for file '+kwargs['filename_data'].split('/')[-1])
      plt.title("MC vs Blacklight")
      #plt.savefig('../plots/cbdisk/ff_only/i45spectrum_comparison.png',dpi=300)
    else:
      flux, frequencies = get_flux(**kwargs)
      plt.plot(frequencies,frequencies*flux)
      plt.xscale('log')
      plt.ylim(1e5, 1e25)
      plt.yscale('log')
      plt.xlabel('Frequency (Hz)')
      plt.ylabel('$\\nu F_\\nu (eV cm^{-2} s^{-1})$ ')
      plt.title('Flux vs Frequency for file '+kwargs['filename_data'].split('/')[-1])
  plt.legend()
  plt.grid()
  plt.savefig('/PellaShared/kcu8rf/blacklight/plots/spherical_thomson/temp_ff_thom.png')



# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename_data', nargs='+', help='name of file containing raw image data')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in gravitational radii')
  parser.add_argument('-m', '--mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('-w', '--width', type=float,
      help='full width of image in gravitational radii')
  parser.add_argument('-l', '--max_level', type=int,
      help='maximum adaptive level to use in calculation')
  parser.add_argument('--luminosity',type=bool,default=False,help='if true, plot luminosity instead of flux')
  parser.add_argument('--multiInc',type=bool,default=False,help='if true, plot multiple inclinations on one plot')
  parser.add_argument('--files',nargs='+',help='list of files to process',type=str)
  parser.add_argument('--inclination',nargs='+',type=float,default=0.0,help='inclination of image (degrees)')
  parser.add_argument('--labels',nargs='+',type=str,default=None,help='labels for the plot')
  parser.add_argument('--compare',type=bool,default=False,help='if true, compare against MC results')
  parser.add_argument('--compare_file',type=str,default=None,help='file containing comparison data')
  parser.add_argument('--compare_file2',type=str,default=None,help='file containing comparison data')
  args = parser.parse_args()
  main(**vars(args))
