from scipy import signal
import numpy as np

def psd(x, dt, nseg):
    """
    Calculates the power spectral density of a given signal using the welch
    method. 

    Parameters
    ----------
    x 
        The time history of the signal.         
    dt
        The time step . 
    nseg
        The the number of segments to average the time series. 

    Returns
    -------
    freq, spectra
        Returns the frequencye and sepctra of the signal
    
    """
    x_no_mean = x - np.mean(x)
    freq, spectra = signal.welch(x_no_mean, fs=1.0/dt, nperseg=len(x_no_mean)/nseg)
       
    return freq[1:], spectra[1:]


def read_forces(file_name):
    """   
    Reads force data agregated over a surface from openfaom file and returns 
    origin(center of rotation), time, and the forces and moments vector for 
    each time step. 
    
    """
    origin = np.zeros(3)
    forces = []
    moments = []
    time = []
    
    #If pressure, viscous and porous are there
    #jump = 9
    
    #If only pressure and viscous
    jump = 6
    
    with open(file_name, "r") as f:
        for line in f:
            if line.startswith('#'): 
                #Read the origin where the force are itegrated
                if line.startswith('# CofR'): 
                    line = line.replace('(','')
                    line = line.replace(')','')
                    line = line.split()
                    origin[0] = line[3] # x-coordinate 
                    origin[1] = line[4] # y-coordinate          
                    origin[2] = line[5] # z-coordinate
                else:
                    continue
            # Read only the pressure part of force and moments. 
            # Viscous and porous are ignored
            else: 
                line = line.replace('(','')
                line = line.replace(')','')
                line = line.split()
                time.append(float(line[0]))
                forces.append([float(line[1]), float(line[2]), float(line[3])])
                moments.append([float(line[jump + 1]), float(line[jump + 2]), float(line[jump + 3])])
    
    time = np.asarray(time, dtype=np.float32)
    forces = np.asarray(forces, dtype=np.float32)
    moments = np.asarray(moments, dtype=np.float32)  
    
    return origin, time, forces, moments

def read_bin_forces(fileName):
    """   
    Reads binData measured at the center of each bin.

    Reads force data agregated over a surface from openfaom file and returns 
    bin heights, time, and the forces and moments vector on the bins for each 
    time step. 
    
    """
    forces = []
    moments = []
    time = []
    nbins = 0
    
    with open(fileName, "r") as f:
        for line in f:
            if line.startswith('#'): 
                #Read the origin where the force are itegrated
                if line.startswith('# bins'): 
                    line = line.replace('# bins', '')                    
                    line = line.replace(':', '')
                    line = line.split()
                    nbins = int(line[0])
                    coords = np.zeros((nbins, 3))
                elif line.startswith('# x co-ords'):
                    line = line.replace('# x co-ords', '')                    
                    line = line.replace(':', '')
                    line = line.split()
                    for i in range(nbins):
                        coords[i,0] = line[i]
                elif line.startswith('# y co-ords'):
                    line = line.replace('# y co-ords', '')                    
                    line = line.replace(':', '')
                    line = line.split()
                    for i in range(nbins):
                        coords[i,1] = line[i]
                elif line.startswith('# z co-ords'):
                    line = line.replace('# z co-ords', '')                    
                    line = line.replace(':', '')
                    line = line.split()
                    for i in range(nbins):
                        coords[i,2] = line[i]
                else:
                    continue
            # Read only the pressure part of force and moments. 
            # Viscous and porous are ignored
            else: 
                line = line.replace('(','')
                line = line.replace(')','')
                line = line.split()
                time.append(float(line[0]))
                story_force = np.zeros((nbins, 3))
                story_moments = np.zeros((nbins, 3))
                for i in range(nbins):
                    start = 18*i + 1

                    #Take only pressure part 
                    story_force[i,0] = line[start] 
                    story_force[i,1] = line[start + 1]
                    story_force[i,2] = line[start + 2]

                    story_moments[i,0] = line[start + 9] 
                    story_moments[i,1] = line[start + 9 + 1]
                    story_moments[i,2] = line[start + 9 + 2]

                forces.append(story_force)
                moments.append(story_moments)
    
    time = np.asarray(time, dtype=np.float32)
    forces = np.asarray(forces, dtype=np.float32)
    moments = np.asarray(moments, dtype=np.float32)  
    
    return coords, time, forces, moments
