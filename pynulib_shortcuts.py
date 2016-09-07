#This is the new version!
import PyNuLib
import numpy as np

def InitializeNuLib():
    parameters_filename = "./parameters"
    PyNuLib.Pynulib.standard_nulib_init()
    PyNuLib.Pynulib.weakrate_init()
    
def GetNuclei():
    species = []

    #print(PyNuLib.Pynulib.nuclei_a)
    
    for i, a in enumerate(PyNuLib.pynulib.nuclei_a):
        species.append((a,PyNuLib.pynulib.nuclei_z[i]))
    return species

def InitWeakRateLib():
    PyNuLib.Pynulib.weakrate_init()

if __name__=="__main__":
    InitializeNuLib()
#    InitWeakRateLib()
#    GetNuclei()
