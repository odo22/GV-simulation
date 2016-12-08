from __future__ import division, print_function, absolute_import

from tmm_core import (inc_tmm, unpolarized_RT, ellips,
                       position_resolved, find_in_structure_with_inf)
import numpy as np
from numpy import pi, linspace, inf, array
import matplotlib.pyplot as plt

# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = pi/180
DEBUG = False

def main():
    lamb=linspace(400,800,num=21)
    angles=linspace(0,37,19)
    N = 20
    
    #Load experimental spectra so they can be plotted together with simulation
    expdata = np.genfromtxt("adultaverage.txt",delimiter="\t")
    expx = expdata[:,0]
    expy = expdata[:,1]

    RNAmeans = np.zeros((N, lamb.size))

    for n in range(N):
        #Function to read off the layers' positions and gray values from ImageJ profile.
        #Subtracts relative position values in the list to get absolute thickeness (in nm).
        t_list,gray_list = loadGrayvalues("Plot Values "+str(n+1)+".xls")
        n_list = generateNlist(gray_list,lamb)

        #Adds top and bottom inf layers to d list
        t_lists = np.append(t_list, inf)
        d_listss = np.insert(t_lists, 0, 21) #167.5
        d_lists = np.insert(d_listss, 0, inf)
        correction = 0.0
        d_list = d_lists*(1-correction)

        # C_list: same lenght as t_list. First and last value must be 'i'.
        a = ['c']*len(t_list)
        b = np.append(a, 'i')
        c_listt = np.insert(b, 0, 'c')
        c_list = np.insert(c_listt, 0, 'i')

        RNAmeans[n,:] = reflectionsim(lamb, angles, n_list, d_list, c_list)

        plt.figure(0)
        plt.plot(lamb,RNAmeans[n,:],'--')

    plt.plot(lamb,RNAmeans.mean(0),'blue',lw=2,label="Average")
    plt.plot(expx,expy,'purple',lw=2,label="Experimental")
    plt.legend()
    plt.xlabel('Wavelength (/nm)')
    plt.ylabel('Fraction reflected')
    plt.title('Simulated reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                'Experimental data (purple)')
    plt.axis([400, 800, 0, 0.35])

# Save simulation data in a text file made out of columns delimited by a space
    
    np.savetxt("/users/Olimpia/Desktop/datafile1.txt", RNAmeans, delimiter=" ")

    profile(n_list)
    plt.show()
    

def profile(n_list):
    if not DEBUG:
        return
    position = np.genfromtxt("Plot Values 20.xls",delimiter="\t")
    position = position[:,0]
    pos = np.delete(position, 0)
    poss = np.delete(pos,1)
    n_list_imaginary = n_list[1:-1].imag

    k = 0.002 #500nm
    plt.figure(1)
    #plt.plot(poss, n_list[1:-1].imag, 'red')
    plt.plot(poss, n_list[1:-1], 'red')
    
    plt.xlabel('Distance (/nm)')
    plt.ylabel('Im{refractive index}')
    #plt.title('Simulated reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                #'Experimental data (purple)')
    #plt.axis([0, 2000, -0.10, 0.10])
    plt.axis([0, 2000, 1.0, 2.0])

def loadGrayvalues(FIN):
    t_list = []
    gray_list = []
    f = open(FIN,'r')
    f.next()    #Skip the first line saying "X Y"
    f.next()    #Skip the first data point in 0.0
    for line in f:
        line = line.split()
        if len(t_list) == 0:
                thickness = float(line[0])
        else:
                for i in range (len(t_list),1):
                    thickness = t_list[i] - t_list[i-1] 
        t_list = t_list + [thickness]
        gray_list = gray_list + [float(line[1])]
    t_list = array(t_list)
    gray_list = array(gray_list)
    return t_list,gray_list


def generateNlist(gray_list,lamb):
    #Conversion from gray values to refractive index
    maximum = 67 #np.amax(gray_list) #71
    minimum = 40 #np.amin(gray_list) #40

    for k in 1/lamb:
        #N1R600 = 1.55 +0j #Real part of n1 at 600nm 1.45
        #N1I600 = 0 + 0.00j #Imaginary part of n1 at 600nm 0.03j
        #BR1 = 88000 # B Coefficient in Cauchy's equation (real)#8800
        #AR1 = N1R600 - BR1/(600**2) #Calculates Cauchy's A from given values of B and n at 600nm
        #BI1 = 0 + 5188j #B coefficient in Cauchy's equation (imaginary)
        #AI1 = N1I600 - BI1/(600**2) #Calculates Cauchy's A from given values of B and n at 600nm

        BR1 = 9464.8
        AR1 = 1.5145

        BR2 = 23700 # B Coefficient in Cauchy's equation (real)
        AR2 = 1.648 #Calculates Cauchy's A from given values of B and n at 600nm
        BI2 = 270 #B coefficient in exponential equation(imaginary part) #270
        AI2 = 0.56 #A coefficient in exponential equation (im part) #0.56

        #n1 = AR1 + BR1*(k**2) + 0*1j+ AI1 + BI1*(k**2)
        n1 = AR1 + BR1*(k**2)
        nmelanin = AR2 + BR2*(k**2) + (AI2*np.exp(-1/(BI2*k)))*1j 
        qpigm = 1
        n2 = n1*(1-qpigm) + nmelanin*qpigm
        
        conversion = (n2-n1)/(maximum-minimum)
        n_list_converted = n2 - conversion*(gray_list-minimum)
        
        wavel=1/k
        n1i = n1.imag
        n2i = n2.imag
        if DEBUG:
            plt.figure(3)
            plt.plot(wavel, n1, marker='o', ms = 10, alpha=1, color='b', label='Chitin layer')
            plt.plot(wavel, n2, marker='o', ms = 10, alpha=1, color='k', label='Melanin layer')
            plt.xlabel('Wavelength (/nm)')
            plt.ylabel('Refractive index')
            plt.title('Dispersion relations, imaginary part. Blue = Chitin layer. Black = Melanin layer')
        
    
#Adds bottom infinite layer at the end of the n list (average of ns)
    n_lists = np.append(n_list_converted, [n2])

#Adds air layer at the top of the n list

    extralayer = 1.6
    n_listt = np.insert(n_lists, 0, extralayer)
    air = 1.00029
    n_list = np.insert(n_listt, 0, air)

    return n_list
    
def reflectionsim(lamb, angles, n_list, d_list, c_list):

    RNA = np.zeros((lamb.size,angles.size))
    
    for i in range(lamb.size):
        print(lamb[i])
        for j in range(angles.size):
            RNA[i,j] = inc_tmm('s',n_list, d_list, c_list, angles[j]*degree, lamb[i])['R']
            RNA[i,j] += inc_tmm('p',n_list, d_list, c_list, angles[j]*degree, lamb[i])['R']
	    #both polarisations, averaged 
    RNA = RNA/2
    RNAmean = np.mean(RNA, axis=1)
    return RNAmean



if __name__ == "__main__":
    main()
