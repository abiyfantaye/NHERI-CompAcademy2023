from scipy import signal
import numpy as np
import matplotlib.gridspec as gridspec    

##############################################################################
#  Sample plots for base load PSD functions 
##############################################################################
n_cmpt = 3
plt, fig = cwe.setup_plot(plt, 28, 24, 22)
gs = gridspec.GridSpec(ncols=1, nrows=n_cmpt, figure=fig)

B = 30.0
H = 450.0

markersize = 8.8
linewidth = 3.0
markeredgewidth = 1.25
titles  = ['(a)', '(b)', '(c)']
error_range = 0.2;
error_x = np.linspace(-10, 10, 3)
nseg = 8
ylabel = ['$fS_{M_x}(f)/\sigma_{M_x}^2$','$fS_{M_y}(f)/\sigma_{M_y}^2$','$fS_{M_z}(f)/\sigma_{M_z}^2$']
y_lim = [[1e-4, 10.0], [1e-4, 10.0], [1e-4, 10.0]]

for i in range(n_cmpt):
    ax = fig.add_subplot(gs[i])
    ax.tick_params(which='major', direction='in', size=10, width=1.5)
    ax.tick_params(which='minor', direction='in', size=6.5,  width=1.0)
    ax.grid(True, which='major', linestyle='-', alpha=0.5)
    ax.grid(True, which='minor', linestyle='--', alpha=0.25)
    
    ax.set_xlim([1e-3, 1.0])
    ax.set_ylim(y_lim[i])
    ax.set_axisbelow(True)

    ax.set_ylabel(ylabel[i])
        
    ax.set_xlabel('$fB/\overline{U}_H$')

    cfd_f, cfd_s = psd(moments[:,i], dt, nseg
    cfd_var = np.var(moments[:,i])

    ax.loglog(cfd_f*B/cfd_Uh, cfd_f*cfd_s/cfd_var, 'r--', linewidth=2)    

    # ax.plot(cfd_f, cfd_f*cfd_s/cfd_var, 'r--', linewidth=2)   

    # ax.legend(['Case-1','Case-2'], loc=0, framealpha=1.0,ncol=2,edgecolor='k')      

fig.set_size_inches(n_directions*22.5/2.54, n_cmpt*15/2.54)
plt.tight_layout()
plt.show()
plt.savefig('plots/base_moments_PSD_V.pdf')
plt.savefig('plots/base_moments_PSD_V.svg')
plt.savefig('plots/base_moments_PSD_V.png')


def setup_plot(plt, font_size=20, legend_font_size=20, axis_font_size=20):
    fig = plt.figure(facecolor='white')
    font = {'family' : 'Times New Roman','weight' : 'normal', 'size'   : font_size}
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10    
    plt.rcParams['xtick.direction'] = 'inout'
    plt.rcParams['ytick.direction'] = 'inout'
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['ytick.major.size'] = 8

    plt.rc('font', **font)    
    plt.rc('axes', labelsize=font_size)    # fontsize of the x and y labels
    plt.rc('axes', titlesize=font_size)  # fontsize of the axes title
    plt.rc('xtick', labelsize=axis_font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=axis_font_size)    # fontsize of the tick labels
    plt.rc('axes', linewidth=1.25)    
    plt.rc('legend', fontsize=legend_font_size)
    plt.rc('text', usetex=True)
    

    return plt, fig


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
