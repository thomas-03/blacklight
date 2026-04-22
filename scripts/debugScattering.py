import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import athena_read
import h5py 

if __name__ == '__main__':
    #df = pd.read_csv('./debugOutput/scatteringDebug.csv',header=None,names=['x1','x2','x3','nu_mc','nu','scat','rho','kt'])
    
    #df['r'] = np.sqrt(df['x1']**2 + df['x2']**2 + df['x3']**2)
    #n_e = df['rho']/mp
    #plt.scatter(df['nu_mc'][::15],df['scat'][::15]/(n_e[::15]*sigma_t),label='blacklight',marker='x')
    #print(df[df['r']==df['r'].abs().min()])

    #jnu should be scattering/(ne sigma_t)
    #ne ~ rho/(mu*mp) 
    sigma_t = 6.65248e-25
    mp = 1.672e-24
    h_erg = 6.626e-27
    c=3e10
    kB = 1.38e-16
    #Bnu = h_erg/c**2*df['nu']**3/(np.exp(h_erg*df['nu']/(kB*1e5))-1)
    #plt.scatter(df['nu'][::7],Bnu[::7],c='red',label='blackbody')

    mcFreqs = [2.713029e+14,3.415501e+14 ,4.299860e+14,5.413204e+14 ,6.814820e+14 ,8.579350e+14 ,1.080076e+15 ,1.359735e+15 ,1.711805e+15,2.155035e+15 ,2.713029e+15,3.415501e+15 ,4.299860e+15 ,5.413204e+15 ,6.814820e+15 ,8.579350e+15 ,1.080076e+16,1.359735e+16 ,1.711805e+16 ,2.155035e+16]
    file = h5py.File('/PellaShared/swd8g/spherical_thomson/thomson/isoth.out2.00000.athdf','r')
    scat = file['mcscat'][()]
    file2 = h5py.File('/PellaShared/swd8g/spherical_thomson/thomson/isoth.out4.00000.athdf','r')
    print(file['x1v'].shape)
    #scat is in the shape nfreq, nMesh, nphi, ntheta, nr
    scat = np.average(scat[:,:,:,:,0],axis=(1,2,3))
    Bnu = 2*h_erg/c**2*np.array(mcFreqs)**3/(np.exp(h_erg*np.array(mcFreqs)/(kB*1e5))-1)
    #scat = scat[:,0,0,0,0]
    plt.plot(mcFreqs,scat/Bnu,c='red',label='MC/blackbody')
    plt.plot(mcFreqs,np.array(mcFreqs)*h_erg,c='blue',linestyle='--',label='$\\nu*h_{erg}$')
    #plt.scatter(mcFreqs,Bnu,c='red',label='Blackbody at 10^5 K')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\nu$')
    plt.ylabel('$ J_\\nu$')
    plt.legend()
    plt.savefig('./plots/spherical_thomson/MCvsBlackbodyRatio.png')