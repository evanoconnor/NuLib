import h5py
import numpy as np

def interp3d(x,y,z,xt,yt,zt,data):
    shape = np.shape(data)
    nx = shape[0] # Ye
    ny = shape[1] # T
    nz = shape[2] # rho

    #!------  determine spacing parameters of (equidistant!!!) table

    dx    = (xt[nx-1] - xt[0]) / (nx-1)
    dy    = (yt[ny-1] - yt[0]) / (ny-1)
    dz    = (zt[nz-1] - zt[0]) / (nz-1)

    dxi   = 1.0 / dx
    dyi   = 1.0 / dy
    dzi   = 1.0 / dz

    dxyi  = dxi * dyi
    dxzi  = dxi * dzi
    dyzi  = dyi * dzi

    dxyzi = dxi * dyi * dzi

    #------- determine location in (equidistant!!!) table 

    ix = 1 + int( (x - xt[0] - 1e-10) * dxi )
    iy = 1 + int( (y - yt[0] - 1e-10) * dyi )
    iz = 1 + int( (z - zt[0] - 1e-10) * dzi )

    ix = max( 1, min( ix, nx ) )
    iy = max( 1, min( iy, ny ) )
    iz = max( 1, min( iz, nz ) )

    #------- set-up auxiliary arrays for Lagrange interpolation

    delx = xt[ix] - x
    dely = yt[iy] - y
    delz = zt[iz] - z

    corners = np.zeros(8)
    corners[0] = data[ix  , iy  , iz  ]
    corners[1] = data[ix-1, iy  , iz  ]
    corners[2] = data[ix  , iy-1, iz  ]
    corners[3] = data[ix  , iy  , iz-1]
    corners[4] = data[ix-1, iy-1, iz  ]
    corners[5] = data[ix-1, iy  , iz-1]
    corners[6] = data[ix  , iy-1, iz-1]
    corners[7] = data[ix-1, iy-1, iz-1]
    
    # coefficients
    
    a1 = corners[0]
    a2 = dxi   * ( corners[1] - corners[0] )       
    a3 = dyi   * ( corners[2] - corners[0] )       
    a4 = dzi   * ( corners[3] - corners[0] )       
    a5 = dxyi  * ( corners[4] - corners[1] - corners[2] + corners[0] )
    a6 = dxzi  * ( corners[5] - corners[1] - corners[3] + corners[0] )
    a7 = dyzi  * ( corners[6] - corners[2] - corners[3] + corners[0] )
    a8 = dxyzi * ( corners[7] - corners[0] + corners[1] + corners[2] +
                   corners[3] - corners[4] - corners[5] - corners[6] )

    return a1 +  a2 * delx \
        +  a3 * dely \
        +  a4 * delz \
        +  a5 * delx * dely \
        +  a6 * delx * delz \
        +  a7 * dely * delz \
        +  a8 * delx * dely * delz     


def interp2d(x, y, xt, yt, data):
    shape = np.shape(data)
    nx = shape[0] # eta
    ny = shape[1] # T

    #------  determine spacing parameters of (equidistant!!!) table

    dx    = (xt[nx-1] - xt[0]) / (nx-1)
    dy    = (yt[ny-1] - yt[0]) / (ny-1)
    
    dxi   = 1.0 / dx
    dyi   = 1.0 / dy
    
    dxyi  = dxi * dyi

    #------- determine location in (equidistant!!!) table 

    ix = 1 + int( (x - xt[0] - 1e-10) * dxi )
    iy = 1 + int( (y - yt[0] - 1e-10) * dyi )
    
    ix = max( 1, min( ix, nx ) )
    iy = max( 1, min( iy, ny ) )

    #------- set-up auxiliary arrays for Lagrange interpolation
    
    delx = xt[ix] - x
    dely = yt[iy] - y

    corners = np.zeros(4)
    corners[0] = data[ix  , iy  ]
    corners[1] = data[ix-1, iy  ]
    corners[2] = data[ix  , iy-1]
    corners[3] = data[ix-1, iy-1]

    #------ set up coefficients of the interpolation polynomial and

    a1 = corners[0]
    a2 = dxi   * ( corners[1] - corners[0] )       
    a3 = dyi   * ( corners[2] - corners[0] )       
    a4 = dxyi  * ( corners[3] - corners[1] - corners[2] + corners[0] )
    
    return a1 +  a2*delx + a3*dely + a4*delx*dely

    

def interpolate_eas(ig, s, rho, T, Ye, table, datasetname):
#!---------------------------------------------------------------------
#!
#!     purpose: interpolation of a function of three variables in an
#!              equidistant(!!!) table.
#!
#!     method:  8-point Lagrange linear interpolation formula          
#!
#!     ig       group number
#!     s        species number
#!     rho      density (g/ccm)
#!     T        temperature (MeV)
#!     Ye       electron fraction
#!
#!     table    h5py.File()
#!     datasetname  {"absorption_opacity", "emissivities", "scattering_opacity", "scattering_delta"}
#!---------------------------------------------------------------------
    data = table[datasetname]

    zt = np.log10(table["rho_points"])
    yt = np.log10(table["temp_points"])
    xt = np.array(table["ye_points"])

    if(rho<table["rho_points"][0] or rho>table["rho_points"][len(zt)-1] \
       or T<table["temp_points"][0] or T>table["temp_points"][len(yt)-1] \
       or Ye<table["ye_points"][0] or Ye>table["ye_points"][len(xt)-1]):
        print("Warning: outside table range")
        return 0

    x = Ye
    y = np.log10(T)
    z = np.log10(rho)

    return interp3d(Ye, np.log10(T), np.log10(rho), xt, yt, zt, data[ig,s,:,:,:])


def interpolate_kernel(ig_in, s, ig_out, eta, T, table, datasetname):
#!---------------------------------------------------------------------
#!
#!     purpose: interpolation of a function of three variables in an
#!              equidistant(!!!) table.
#!
#!     method:  8-point Lagrange linear interpolation formula          
#!
#!     ig_in    group number
#!     s        species number
#!     ig_out   group number
#!     eta      mue / T (dimensionless)
#!     T        temperature (MeV)
#!
#!     table    h5py.File()
#!     datasetname  {"inelastic_phi0", "inelastic_phi1"}
#!---------------------------------------------------------------------
    data = table[datasetname][ig_in,s,ig_out,:,:]

    xt = np.log10(table["eta_Ipoints"])
    yt = np.log10(table["temp_Ipoints"])

    if(eta<table["eta_Ipoints"][0] or eta>table["eta_Ipoints"][len(xt)-1] \
       or T<table["temp_Ipoints"][0] or T>table["temp_Ipoints"][len(yt)-1]):
        print("Warning: outside table range")
        return 0

    x = np.log10(eta)
    y = np.log10(T)

    return interp2d(x, y, xt, yt, data)

def interpolate_eos(rho, T, Ye, table, datasetname):
    data = table[datasetname]

    zt = table["logrho"]
    yt = table["logtemp"]
    xt = table["ye"]

    z = np.log10(rho)
    y = np.log10(T)
    x = Ye

    return interp3d(x, y, z, xt, yt, zt, data)
