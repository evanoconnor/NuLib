import matplotlib
matplotlib.use('Agg')

import numpy as np
import pylab
import PyNuLib
import random
from scipy.interpolate import interp1d
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from collections import namedtuple
import os.path as path
from pylab import *
import sys
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import code
import cPickle as pickle
import math

from pynulib_shortcuts import *

def setfont(font='helvetica',unicode=True): 
    r""" 
    Set Matplotlibs rcParams to use LaTeX for font rendering. 
    Revert all changes by calling rcdefault() from matplotlib. 
    
    Parameters: 
    ----------- 
    font: string 
        "Helvetica" 
        "Times" 
        "Computer Modern" 
    
    usetex: Boolean 
        Use unicode. Default: False.     
    """ 
    
    # Use TeX for all figure text! 
    pylab.rc('text', usetex=True) 

    font = font.lower().replace(" ","") 
    if font == 'times': 
        # Times 
        font = {'family':'serif', 'serif':['Times']} 
        preamble  = r""" 
                       \usepackage{color} 
                       \usepackage{mathptmx} 
                    """ 
    elif font == 'helvetica': 
        # Helvetica 
        # set serif, too. Otherwise setting to times and then 
        # Helvetica causes an error. 
        font = {'family':'sans-serif','sans-serif':['Helvetica'], 
                'serif':['cm10']} 
        preamble  = r""" 
                       \usepackage{color} 
                       \usepackage[tx]{sfmath} 
                       \usepackage{helvet} 
                    """ 
    else: 
        # Computer modern serif 
        font = {'family':'serif', 'serif':['cm10']} 
        preamble  = r""" 
                       \usepackage{color} 
                    """ 
    
    if unicode: 
        # Unicode for Tex 
        #preamble =  r"""\usepackage[utf8]{inputenc}""" + preamble 
        # inputenc should be set automatically 
        pylab.rcParams['text.latex.unicode']=True 
    
    #print font, preamble 
    pylab.rc('font',**font) 
    pylab.rcParams['text.latex.preamble'] = preamble 

def input_parser(file,searchExp):
    input_file = open(file)
    for line in input_file:
        if searchExp in line:
            if searchExp == 'Directory':
                return line.split()[1].rstrip()
            if searchExp == 'Runs':
                files = []
                while True:
                    try:
                        line = next(input_file)
                        if line == '\n':
                            """ Break if empty line is 
                            found after reading last state """
                            break
                    except:
                        break
                    files.append(line.rstrip())
                return files    

def readfile(filename):
    return_data =[]
    file_data = []
    flg_read = False
    counter = 0.0        
    for line in open(filename):            
        counter += 1.0
        entries = line.split()
        try:
            entries = map(float,entries)
            file_data.append(entries)
        except:
            if flg_read == False:
                if entries[0] == 'Interpolant':
                    flg_read = True
                    continue
            else:
                entries = map(float,entries)
                file_data.append(entries[0:6])
    return_data.append(file_data)
    if len(return_data)>1:
        return return_data
    else:
        return file_data

def setticks(ax,xmajor=5,xminor=1,ymajor=2,yminor=0.5):

    # Set up ticks
    xmajorLocator   = MultipleLocator(xmajor)
    xmajorFormatter = FormatStrFormatter('%d')
    xminorLocator   = MultipleLocator(xminor)            
    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)
    ax.xaxis.set_minor_locator(xminorLocator)

    ymajorLocator   = MultipleLocator(ymajor)
    ymajorFormatter = FormatStrFormatter('%d')
    yminorLocator   = MultipleLocator(yminor)            
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)
    ax.yaxis.set_minor_locator(yminorLocator)

    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    for tick in ax.xaxis.get_ticklines():
        tick.set_markersize(2)
    for tick in ax.yaxis.get_ticklines():
        tick.set_markersize(2)
    for tick in ax.xaxis.get_ticklines(minor=True):
        tick.set_markersize(1)
    for tick in ax.yaxis.get_ticklines(minor=True):
        tick.set_markersize(1)


##########################################

def nan(shape):
    out = np.zeros(shape = shape)
    out.fill(np.NAN)
    return out

nuclei = namedtuple('nuclei',['ntime','Z','A','N','DYE','time'])

class nuclear_chart(object):
    """
    Reads in a table of binding energies as provided by Morten.
    Provides entries from (A,Z), (A,N), or (Z,N) requests.
    """
    def __init__(self):
        self.table = {} 
        self.ntimes = set()
        self.nuclei = set()
        # az = {}
        # nline = 0
        # times = [0.00,0.08,0.15,0.20]
        # for line in open(datafile):            
        #     entries = line.split()
        #     if len(entries)==16280:
        #         az = entries
        #     else:
        #         #if round(float(entries[0]),5) in times:
        #         if 1==1:
        #             print nline
        #             for i in range (0,len(entries)-1):
        #                 k = i*2
        #                 mylist = [int(nline),int(az[k+1]),int(az[k]),(int(az[k])-int(az[k+1])),float(entries[i+1]),float(entries[0])]
        #                 self.Addnuclei(mylist)      
        #         nline = nline + 1
        
    def AddNuclei(self,values):
        """Adds an isotope to the table.
        nucleis are indexed by their (Z,A) tuple."""
        ntime,Z,A,N,DYE,time = values
        ntime = int(ntime)
        Z = int(Z)
        A = int(A)
        N = int(N)
        DYE = float(DYE)
        time = float(time)
        self.table[(ntime,Z,A)] = nuclei(ntime,Z,A,N,DYE,time)   #Table is filled with lists of isotope information
        self.ntimes.add(ntime)
        self.nuclei.add((A,Z))
    def __getitem__(self,val):
        return self.table[val]
    def __iter__(self):
        return iter(self.table.values())
    def __contains__(self,val):
        return val in self.table
    def Get(self,ntime=None,A=None,Z=None,N=None):
        if ntime is not None and A is not None and Z is not None:
            return self[ntime,Z,A]
        elif ntime is not None and Z is not None and N is not None:
            return self[ntime,Z,Z+N]
        elif ntime is not None and A is not None and N is not None:
            return self[ntime,A-N,A]
        elif A is not None and Z is not None:
            return_list = []
            for i in self.ntimes:
                return_list.append(self[i,Z,A])
            return return_list
        else:
            raise Exception('Not enough information')


def colorbar_index(ncolors, cmap):
    cmap = cmap_discretize(cmap, ncolors)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(ncolors))

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
  
def do_plot_mask():
    plot_size = (70,120)
    test = nan(shape=plot_size)
    stable = path.dirname(path.abspath(__file__))
    for line in open(path.join(stable,'stableisotopes.dat')):
        entries = line.split()
        if int(entries[0])<plot_size[0] and int(entries[1])<plot_size[1]:
            test[int(entries[0]),int(entries[1])]=1000
    vals_masked = np.ma.array(test,mask=np.isnan(test))
    cmap = plt.get_cmap('gray')
    cmap.set_bad(color='w',alpha=1.0)
    chart = pylab.pcolor(vals_masked,edgecolor='k',linewidth=0.4,cmap='spectral_r',alpha=1.0)



def do_plot(vals,title=None,xlabel=None,ylabel=None):
    vals_masked = np.ma.array(vals,mask=np.isnan(vals))
    lower_bound = vals_masked < 1.0e-10
    vals_masked[lower_bound] = 0.0
    #vals_masked = vals_masked*0.0
    n_colours = 16
    new_cmap = truncate_colormap(cm.get_cmap('jet'), 0.0, 1.0, n=n_colours)
    #new_cmap.set_bad(color='#0000ff')
    new_cmap.set_bad(color='#cdcdcd')
    #chart = pylab.pcolor(vals_masked.T,cmap=new_cmap,edgecolor='w',linewidth=1.75,norm=LogNorm(vmin=1.0e-10,vmax=1.0))
    chart = pylab.pcolor(vals_masked.T,cmap=new_cmap,edgecolor='w',linewidth=0.75,norm=LogNorm(vmin=1.0e-10,vmax=1.0e0))
    #cbar = plt.colorbar(chart)

    if title is not None:
        pylab.title(title)
    if xlabel is not None:
        pylab.xlabel(xlabel)
    if ylabel is not None:
        pylab.ylabel(ylabel)

def make_chart(table,save=None):
    plot_size = (120,70)
    dyedt1 = nan(shape=plot_size)
    
    for iso in table:
#        if iso.ntime == 0:
 #           time1 = iso.time
        if iso.N<plot_size[0] and iso.Z<plot_size[1]:
            dyedt1[iso.N,iso.Z] = iso.DYE    

    pylab.figure(figsize=(14,6.6))
    pylab.subplot(1,1,1)
    do_plot(dyedt1,
            title="",
            xlabel="Neutron Number",
            ylabel="Proton Number")

    # if save is None:
    #     pylab.show()
    # else:
    #     pylab.savefig(save)
    #print "Done, saving..."
    #pickle.dump(table,open("dyedt_table.p","wb"))
    


top_contributers = []

#top_contributers = set()
def dyedt_top_integrated(table,ntop=5):
    
    dyesum = np.zeros([len(table.nuclei)])
    
    # loop over table, sums dyedt for all time steps
    for i,nuc in enumerate(table.nuclei):  
        dyesum[i] = sum([x.DYE for x in table.Get(A=nuc[0],Z=nuc[1])])

    max_dye_index = list(dyesum.argsort()[-ntop:][::-1])
    
    # loop over table, if index of nucleus is in the top 10/100/1000/ntop (in max_dyedt_index)
    # then include it in the top contributers list of tuples
    for i,nuc in enumerate(table.nuclei):  
        if i in max_dye_index:
            dyedt = [x.DYE for x in table.Get(A=nuc[0],Z=nuc[1])][0]
            top_contributers.append((nuc[0],nuc[1],dyedt))


def readM1xg(filename):
    data = []
    for i,line in enumerate(open(filename)):
        if 'Time' in line:
            if i!=0:
                data.append(data_tuple)
            line = line.split()
            data_tuple = {'time':float(line[2]),'prob_dist':[]}
        else:
            line = line.split()            
            if len(line)>1:
                data_tuple['prob_dist'].append([float(line[0]),float(line[1])])
    data.append(data_tuple) # get's the last tuple
    return data

def tbounce(dir_path):
    for line in open(dir_path+'tbounce.dat'):
        return float(line.split()[0])
    


if __name__=='__main__':
    setfont()
    rhofile = 'trajectory.dat'
    inputfile = './inputfiles.dat'
    data_dir = input_parser(inputfile,'Directory')
    out_dir = '/mnt/simulations/ceclub/tav/gr1d/'
    runs = input_parser(inputfile,'Runs')
    InitializeNuLib()
    species=GetNuclei()
    # index of trajectory point to plot
    # do_chart = [2000] 
    do_chart = []
    for nrun,run in enumerate(runs):
        filetoload = data_dir+run+'/'+rhofile
        M1_file = out_dir+run+'/M1_nue_enspectra_cen.xg'
        tb = tbounce(data_dir+run+'/')
        M1_data = readM1xg(M1_file)
        xtimes = []
        xtemps = []
        xrhos = []
        xyes = []
        ye = []
        profile = readfile(filetoload)
        for i in range(0,len(PyNuLib.pynulib.global_dyedt)):
            PyNuLib.pynulib.global_dyedt[i]=0.0
        PyNuLib.pynulib.global_dyedt_freep=0.0

        table = nuclear_chart()
        table_background = nuclear_chart()
        # loop over time steps of simulation
        for i,grid_pt in enumerate(profile):            

            dt = profile[i][0]-profile[i-1][0]

            #if i > 10: 
            #    continue
            #tfactor = 0.99
            tfactor = 0.01
            
            if grid_pt[0]>(tb-0.002)*tfactor:
                continue
            else:
                if i>0:
                    print "\nPercent of time to bounce: ",100.0*grid_pt[0]/((tb-0.002)*tfactor)," %   dt=",grid_pt[1],grid_pt[2],grid_pt[3]
            if len(do_chart)!=0:
                if i not in do_chart:
                    continue
            if i == 0:
                Ye_0 = grid_pt[3]
            rho = grid_pt[1]
            T = grid_pt[2]
            Ye = grid_pt[3]                    
            alpha = grid_pt[6]
            xtimes.append(grid_pt[0]-tb)
            xyes.append(Ye)
            xrhos.append(rho)
            
            # must be called before get_dyedt
            f_nu = np.asarray(M1_data[i]['prob_dist'])[:,1]
            PyNuLib.Pynulib.calculate_neutrino_blocking(rho,T,Ye,alpha,f_nu)
            # calculate dyedt (inside fortran routines)
            PyNuLib.Pynulib.get_dyedt(rho,T,Ye) 

            if i == 0:                
                # making an empyty background table for plotting the outline of the nuclear chart
                for k,nuc in enumerate(species):                
                    table_background.AddNuclei([i,nuc[1],nuc[0],nuc[0]-nuc[1],0.0,grid_pt[0]])
            else:
                # loop over all nuclei to calculate dyedt for each
                for k,nuc in enumerate(species):
                    dyedt = PyNuLib.pynulib.global_dyedt[k]
                    table.AddNuclei([i,nuc[1],nuc[0],nuc[0]-nuc[1],dyedt*dt,grid_pt[0]]) 
                
        # calculate dye
        dyesum = np.zeros([len(species)])
        #ntop = 2500
        ntop = 500
        for k,nuc in enumerate(species):
            #dyesum[nuc[0]-nuc[1],nuc[1]) = sum([timestep.DYE for timestep in table.Get(nuc[0],nuc[1])])
            dyesum[k] =  sum([timestep.DYE for timestep in table.Get(A=nuc[0],Z=nuc[1])])
        max_dye_index = list(dyesum.argsort()[-ntop:][::-1])      
       
        plot_size = (2*120,2*70)
        centroids = nan(shape=plot_size)                
        # loop over global list of the species with the largest value of dyedt in each time step
        for k,nuc in enumerate(species): # nuc = (A, Z, dyedt)
            if k in max_dye_index:
                if nuc[0]-nuc[1] < plot_size[0] and nuc[1] < plot_size[1]:
                    if centroids[nuc[0]-nuc[1],nuc[1]] != centroids[nuc[0]-nuc[1],nuc[1]]:  # if nan then 0 it
                        centroids[nuc[0]-nuc[1],nuc[1]]=0.0
                    if i != 0: # 
                        centroids[nuc[0]-nuc[1],nuc[1]] += dyesum[k] # += tdiff*dyedt - integrating
        vals_masked = np.ma.array(centroids,mask=np.isnan(centroids))
        make_chart(table_background)
        minimum = min([x for x in centroids.flatten() if x==x])
        maximum = max([x for x in centroids.flatten() if x==x])
        chart = pylab.pcolor(vals_masked.T,cmap='jet',edgecolor='w',linewidth=0.75,norm=LogNorm(vmin=minimum,vmax=maximum))
        cbar = pylab.colorbar(chart)
        pylab.xlim(0,120)
        pylab.ylim(0,70)
        pylab.savefig("/user/taverner/public_html/dyedt_top"+str(ntop)+"_"+str(tfactor)+"tb_overall_log.pdf")
        print "/dyedt_top"+str(ntop)+"_"+str(tfactor)+"tb_overall_log.pdf"
        code.interact(local=locals())



            

