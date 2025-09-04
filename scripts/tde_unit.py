import numpy as np
from scipy.interpolate import griddata
import astropy.constants as cons
#Constants and units
#mp=1.6726e-24 #they use mh=1.6733e-24
mp=cons.m_p.cgs.value
#kB=1.3807e-16 #they use 1.380658e-16
kB=cons.k_B.cgs.value
#Lsun = 3.85e33
Lsun=cons.L_sun.cgs.value
#Msun=1.9891e33 #same
Msun=cons.M_sun.cgs.value
#Rsun=6.955e10
Rsun=cons.R_sun.cgs.value
#year = 3.155815e7
#pc = 3.0857e18
#AU=1.496e13
#c=2.99792458e10 #same
c=cons.c.cgs.value
#G=6.67259e-8 #same
G=cons.G.cgs.value
#a_rad = 7.5657e-15  #radiation constant they use: 7.5646e-15
a_rad=4*(cons.sigma_sb/cons.c).cgs.value
kes=0.34   #scattering opacity for solar composition
#Mbh=1e7*Msun   #blackhole mass
Mbh=3.0e6*Msun #black hole mass for injection at the boundary
Medd = 40*np.pi*G*Mbh/kes/c #Eddington accretion rate (g/s) (assuming eta=0.1)
v0=0.0005*c   #velocity unit
Crat=c/v0
rho0=1e-8  #density unit
r0=2*G*Mbh/c**2  #Schwarzchild radius
rs=r0
k0=1/(rho0*r0)  #scattering opacity unit
t0=r0/v0 #time unit
mu=0.6
T0=v0**2*mu*mp/kB #temperature unit
Prat=a_rad*pow(T0, 4)/(rho0*pow(v0, 2))
Meddcode=Medd/(r0**3*rho0/t0) #in code units
gamma=1.6666666666666667
Fr0=c*a_rad*T0**4  #radiation flux units
P0=rho0*v0**2 #gas pressure units
Er0=a_rad*T0**4 #radiation energy density units
Ledd=4*np.pi*G*Mbh*c/kes #Eddington luminosity
B0=np.sqrt(4*np.pi*rho0)*v0 #magnetic field units

#print('tstream (code units): ',1e5/T0)
#print('tfloor (cgs): ',2e-4*T0)


def scinot(x,sigfig=2):
    arg="{:."+str(sigfig)+"e}"
    return arg.format(x)

#print('prat:',scinot(Prat,10),' crat:',Crat,' Tunit:',scinot(T0,10),' lunit:',scinot(r0,10))

#print(scinot(P0,sigfig=4))
 
def cgs(var,varname):
    '''Converts variables in code units to cgs units'''
    '''fluxes refer to quantity/time'''
    varunit={'x1f':r0,'x1v':r0,'rho':rho0,'press':P0,'vel1':v0,'vel2':v0,'vel3':v0,
             'Er':Er0,'Fr1':Fr0,'Fr2':Fr0,'Fr3':Fr0,'Pr11':Er0, 'Pr22':Er0,'Pr33':Er0,
             'Pr12':Er0,'Pr13':Er0,'Pr23':Er0,'Pr21':Er0,'Pr31':Er0,'Pr32':Er0,
             'Er0':Er0,'Fr01':Fr0,'Fr02':Fr0,'Fr03':Fr0,'Sigma_s_0':1/r0,'Sigma_a_0':1/r0,
             'Sigma_p_0':1/r0,'time':t0,'Bcc1':B0,'Bcc2':B0,'Bcc3':B0,'divB':B0/r0,'wabsdivB':1,
             'wabsdivB2':1,'rhoghost':rho0,'pressghost':P0,'bcc1ghost':B0,'bcc2ghost':B0,'bcc3ghost':B0,
             'area':r0**2,'temp':T0,'lum':Fr0*r0**2,'KE':rho0*r0**3*v0**2,'massflux':rho0*v0*r0**2} 
    return var*varunit[varname]

def spherical_to_cartesian(Vr,Vtheta,Vphi,r,theta,phi,reduced=False): 
    '''Converts vectors from spherical polar coordinates to cartesian coordinates'''

    if reduced: #theta has been averaged previously
        phi_grid, r_grid = np.meshgrid(phi, r, indexing='ij')
        thetamid=theta[int(len(theta)/2)] #might not be okay? taking value closest to np.pi/2 (the previous averaging theta)
        Vx=Vr*np.sin(thetamid)*np.cos(phi_grid)-Vphi*np.sin(phi_grid)
        Vy=Vr*np.sin(thetamid)*np.sin(phi_grid)+Vphi*np.cos(phi_grid)
        Vz=Vr*np.cos(thetamid)
    else:
        phi_grid, theta_grid, r_grid = np.meshgrid(phi, theta, r, indexing='ij')
        Vx=Vr*np.sin(theta_grid)*np.cos(phi_grid)+Vtheta*np.cos(theta_grid)*np.cos(phi_grid)-Vphi*np.sin(phi_grid)
        Vy=Vr*np.sin(theta_grid)*np.sin(phi_grid)+Vtheta*np.cos(theta_grid)*np.sin(phi_grid)+Vphi*np.cos(phi_grid)
        Vz=Vr*np.cos(theta_grid)-Vtheta*np.sin(theta_grid)
    return Vx, Vy, Vz

def spherical_to_cylindrical(Vr,Vtheta,Vphi,r,theta,phi,reduced=False):
    '''Converts vectors from spherical polar coordinates to cylindrical coordinates'''

    phi_grid, theta_grid, r_grid = np.meshgrid(phi, theta, r, indexing='ij')

    Vrho=Vr*np.sin(theta_grid)+Vtheta*np.cos(theta_grid)
    Vz=Vr*np.cos(theta_grid)-Vtheta*np.sin(theta_grid)
    return Vrho, Vphi, Vz
    
def rhophiz(r,theta,phi):
    rho=r*np.sin(theta)
    phi=phi
    z=r*np.cos(theta)
    return rho,phi,z

def subsample_2d(q,u,x,y,step,average=True):
    """
    Subsampling of 2d array by factor of step, step must be integer
    """
   
    nx = x.shape[0]
    ny = x.shape[1]

    if (not average):
        x = x[step // 2:nx:step,step // 2:ny:step]
        y = y[step // 2:nx:step,step // 2:ny:step]
        q = q[step // 2:nx:step,step // 2:ny:step]
        u = u[step // 2:nx:step,step // 2:ny:step]
        return q,u,x,y
    else:
        # too lazy to work out a more pythony way of doing this
        xp = np.zeros((nx // step,ny // step))
        yp = np.zeros((nx // step,ny // step))
        qp = np.zeros((nx // step,ny // step))
        up = np.zeros((nx // step,ny // step))
        for i in range(nx // step):
            for j in range(ny // step):
               xp[i,j] = np.average(x[i*step:(i+1)*step,j*step:(j+1)*step])
               yp[i,j] = np.average(y[i*step:(i+1)*step,j*step:(j+1)*step])
               qp[i,j] = np.average(q[i*step:(i+1)*step,j*step:(j+1)*step])
               up[i,j] = np.average(u[i*step:(i+1)*step,j*step:(j+1)*step])

        return qp,up,xp,yp

def interp_even_grid(xc, yc, arr, coordout=0):
        '''
        read in 3D coordinate xc, yc that convereted from spherical polar,
        xc, yc are non-uniform Cartesian grids
        arr is the array of quantites we want to interpolate
        '''

        #for high resolution data, can skip some data points to save interpolation time
        step = 1
        xc_irr, yc_irr = xc[::step, ::step], yc[::step, ::step]
        coord_irr = np.append(xc_irr.reshape(-1,1),yc_irr.reshape(-1,1),axis=1)
        arr_irr = arr[::step, ::step].reshape(-1,1)
    
        #prepare uniform Cartesian mesh: x_rr, y_rr
        #face centered mesh xf_rr, yf_rr are not used for interpolation, only for plot
        xc_min, xc_max = np.min(xc_irr), np.max(xc_irr)
        yc_min, yc_max = np.min(yc_irr), np.max(yc_irr)
        nx1, nx2 = np.shape(xc_irr)[0], np.shape(xc_irr)[1]
        x_rr, y_rr = np.meshgrid(np.linspace(xc_min, xc_max, nx1, endpoint=True),\
                                 np.linspace(yc_min, yc_max, nx2, endpoint=True))
        xf_rr, yf_rr = np.meshgrid(np.linspace(xc_min, xc_max, nx1+1, endpoint=True),\
                                   np.linspace(yc_min, yc_max, nx2+1, endpoint=True))

        #griddata is the intepolation package
        #see documentation to pick different interpolation method
        grid_arr = griddata(coord_irr, arr_irr, (x_rr, y_rr), method='cubic')

        #return coordinate for plotting
        if coordout==0:
            return grid_arr[:,:,0]
        elif coordout==1:
            return x_rr, y_rr, xf_rr, yf_rr, grid_arr[:,:,0]
    
def load_vars(var,massweighted=False,photosphere=False,vector=False,var1=None,var2=None,var3=None,quiver=False,streamlines=False,qvar=None,ghost=False):
    qs=[]
    if massweighted:
        qs.append('rho')
    if photosphere:
        if kappa=='scat':
            qs.append('Sigma_s_0')
        elif kappa=='stara':
            qs.append('Sigma_s_0')
            qs.append('Sigma_a_0')
        elif kappa=='starp':
            qs.append('Sigma_s_0')
            qs.append('Sigma_p_0')
        else:
            qs.append('Sigma_s_0')
            qs.append('Sigma_a_0')
            qs.append('Sigma_p_0')
    if vector:
        qs.append(var1)
        qs.append(var2)
        qs.append(var3)
    elif var=='pbeta':
        if ghost:
            qs.append('pressghost')
            qs.append('bcc1ghost')
            qs.append('bcc2ghost')
            qs.append('bcc3ghost')
        else:
            qs.append('press')
            qs.append('Bcc1')
            qs.append('Bcc2')
            qs.append('Bcc3')
    elif var=='magflux':
        qs.append('rho')
        qs.append('Bcc1')
    elif var=='Tgas' or var=='specific_entropy':
        if ghost:
            qs.append('rhoghost')
            qs.append('pressghost')
        else:
            qs.append('rho')
            qs.append('press')
    elif var=='Trad':
        qs.append('Er')
    elif var=='lum':
        qs.append('Fr1')
    elif var=='Tratio':
        qs.append('rho')
        qs.append('press')
        qs.append('Er')
    elif var=='zeta':
        qs.append('Er')
        qs.append('Bcc1')
        qs.append('Bcc2')
        qs.append('Bcc3')
    elif var=='maglum':
        if ghost:
            qs.append('bcc1ghost')
        else:
            qs.append('Bcc1')
            qs.append('vel1')
    else:
        qs.append(var)
    if quiver or streamlines:
        qs.append(qvar[0])
        qs.append(qvar[1])
        qs.append(qvar[2])
    return list(set(qs))
