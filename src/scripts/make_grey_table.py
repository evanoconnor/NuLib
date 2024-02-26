import numpy as np
from scipy import integrate
from scipy import interpolate
import h5py

def save_table_as_hf5f(table: dict, table_name: str = "table.h5") -> None:
    """
    Saves a table of data to an HDF5 file.

    This function iterates over the keys in a dictionary, saving each
    numpy array to an HDF5 file. The arrays are transposed before saving to match the
    expected layout in the HDF5 file format. 

    Parameters:
    - table (dict): The table to save, as a python dict. Each key will be saved as a hdf5 dataset.
    - table_name (str, optional): The name of the HDF5 file to save the data to. Defaults to "table.h5".

    Returns:
    - None: This function does not return a value but saves the data to an HDF5 file.
    """
    # Open or create the HDF5 file for writing
    hf = h5py.File(table_name, 'w')
    # Iterate over each key-value pair in the table dictionary
    for key in table.keys():
        # Create a dataset for each key in the HDF5 file and save the transposed numpy array
        hf.create_dataset(key, data=table[key].T)
    # Close the HDF5 file to ensure data is written properly and resources are freed
    hf.close()



def read_nulib_table(fil: str) -> dict:
    """
    Reads data from a NuLib table file and organizes it into a dictionary.

    This function opens a NuLib table file, reads
    temperature, density, electron fraction, scattering
    opacity, absorption opacity, neutrino energies, emissivities, and bin widths.
 
    Parameters:
    - fil (str): The file path to the NuLib table file, typically an HDF5 file.

    Returns:
    - dict: A dictionary containing arrays of neutrino energies, densities (rho),
      temperatures (temp), electron fractions (ye), scattering opacities,
      absorption opacities, emissivities, and bin widths. T

    The returned dictionary keys are: "energy", "rho", "temp", "ye", "scatter",
    "absorption", "emissivities", and "bin_width", each mapping to the corresponding
    numpy array extracted from the NuLib table file.
    """

    table = h5py.File(fil, 'r') 
    nulib_temp = np.asarray(table["temp_points"])
    nulib_rho = np.asarray(table["rho_points"])
    nulib_ye = np.asarray(table["ye_points"])
    scatter = np.asarray(table["scattering_opacity"])
    abso = np.asarray(table["absorption_opacity"])
    energy_bins = np.array(table["neutrino_energies"]) 
    emiss = np.asarray(table['emissivities'])
    bw = np.asarray(table['bin_widths'])

    nulib_table = {"energy" : energy_bins, "rho" : nulib_rho, "temp" : nulib_temp, "ye" : nulib_ye,
                   "scatter" : scatter, "absorption" : abso, "emissivities" : emiss, 'bin_width' : bw }

    return nulib_table



def read_eos_table(fil: str) -> dict:
    """
    Reads data from an EOS table file and organizes it into a dictionary.

    This function opens an EOS table file and extracts the 
    chemical potential density, temperature, and electron fraction.
    The chemical potentials are adjusted to represent the total muon chemical potential.

    Parameters:
    - fil (str): The file path to the EOS table file, typically an HDF5 file.

    Returns:
    - dict: A dictionary containing arrays for chemical potential ("mu"),
      log-scaled densities ("log_rho"), log-scaled temperatures ("log_temp"), and electron
      fractions ("ye").

    #IMPORTANT: This function is written for specific EOS tables. You might have to change it.
    """


    eosfile = h5py.File(fil, 'r')
    eos_mu = np.asarray(eosfile["mu_e"]) - np.asarray(eosfile["mu_n"]) + np.asarray(eosfile["mu_p"])
    eos_logden = np.asarray(eosfile["logrho"])
    eos_logtemp = np.asarray(eosfile["logtemp"])
    eos_ye = np.asarray(eosfile["ye"])

    eos_table = {"mu" : eos_mu, "log_rho" : eos_logden, "log_temp" : eos_logtemp, "ye" : eos_ye}
    return eos_table


def integrate_table(nulib_table: dict, eos_table: dict) -> np.ndarray:
    """
    Integrates neutrino interaction data over energy to compute interaction rates.
    
    This function calculates the neutrino absorption, emission, and scattering rates over a range of conditions
    specified by the NuLib and EOS tables. It is assumed that neutrinos follow a Fermi-Dirac distribution, and
    the energy dependence of the interaction rates is integrated out using the temperature data provided in the
    nulib_table.
    
    Parameters:
    - nulib_table (dict): Contains neutrino interaction data from the NuLib table, including neutrino energies,
      densities ("rho"), temperatures ("temp"), electron fractions ("ye"), scattering and absorption opacities.
    - eos_table (dict): Contains state of matter data from the EOS table, including chemical potentials ("mu"),
      log-scaled densities ("log_rho"), log-scaled temperatures ("log_temp"), and electron fractions ("ye").
    
    Returns:
    - np.ndarray: A multidimensional array containing integrated values for absorption, emission, and scattering rates,
      structured as [neutrino species, interaction type, density, temperature, electron fraction].
    
    Notes:
    - This function calculate the rates under the asumption that the neutrinos have the same temperature as the fluid.
    """
    # Factors to adjust the energy terms for different neutrino species
    mfac = np.array([1.0, -1.0, 0.0])
    # Coefficients for the distribution calculations
    fac = np.array([1.0, 1.0, 4.0])
    # Initialize the output array with zeros
    table = np.zeros([3, 4, len(nulib_table["rho"]), len(nulib_table["temp"]), len(nulib_table["ye"])])
    # Points for EOS interpolation
    points = (eos_table["ye"], eos_table["log_temp"], eos_table["log_rho"])
    # Neutrino energy array
    energy = np.array(nulib_table["energy"])

    # Generate meshgrids for interpolation
    xv, yv, zv = np.meshgrid(nulib_table["ye"], np.log10(nulib_table["temp"]), np.log10(nulib_table["rho"]), indexing="ij")
    # Interpolate the chemical potentials from the EOS table
    m1 = interpolate.interpn(points, eos_table["mu"], (xv, yv, zv))

    for n in range(0, 3):  # Loop over neutrino species
        # Calculate the eta factor for the Fermi-Dirac distribution
        eta = (energy[:, None, None, None] - mfac[n] * m1) / nulib_table["temp"][None, None, :, None]
        # Compute the distribution for neutrino interactions
        bb = fac[n] / (1.0 + np.exp(eta)) * (energy[:, None, None, None]) ** 3
        # Integrate over energy to get baseline interaction rates
        bb_int = np.trapz(bb, x=energy, axis=0)
        
        # Absorption calculations
        absorption = bb * nulib_table["absorption"][:, n, :, :, :]
        absorption_int = np.trapz(absorption, x=energy, axis=0)
        table[n, 0, :, :, :] = absorption_int.T
        # Normalized absorption rates
        table[n, 1, :, :, :] = np.divide(absorption_int, bb_int, out=np.zeros_like(absorption_int), where=bb_int != 0).T

        # Emission rate calculations
        emiss = np.zeros_like(absorption)
        for i in range(len(energy)):
            emiss[i, :, :, :] = absorption[i, :, :, :] / energy[i]
        emiss_int = np.trapz(emiss, x=energy, axis=0)
        table[n, 3, :, :, :] = emiss_int.T

        # Scattering rate calculations
        scatter = bb * nulib_table["scatter"][:, n, :, :, :]
        scatter_int = np.trapz(scatter, x=energy, axis=0)
        table[n, 2, :, :, :] = np.divide(scatter_int, bb_int, out=np.zeros_like(scatter_int), where=bb_int != 0).T

    return table



def integrate_table_n_temp(nulib_table: dict, eos_table: dict, ntemp) -> np.ndarray:
    """
    Integrates neutrino interaction data over energy to compute interaction rates for a given neutrino temperature.
    
    This function calculates the neutrino absorption, emission, and scattering rates across a range of conditions,
    assuming the neutrinos follow a Fermi-Dirac distribution. The energy dependence of the rates is integrated out
    using the specified neutrino temperature, ntemp.
    
    Parameters:
    - nulib_table (dict): Contains neutrino interaction data from the NuLib table, including energies, densities ("rho"),
      temperatures ("temp"), electron fractions ("ye"), and opacities.
    - eos_table (dict): Contains state of matter data from the EOS table, including chemical potentials ("mu"),
      log-scaled densities ("log_rho"), log-scaled temperatures ("log_temp"), and electron fractions ("ye").
    - ntemp (float): The neutrino temperature used to compute the Fermi-Dirac distribution.
    
    Returns:
    - np.ndarray: A multidimensional array containing integrated values for absorption, emission, and scattering rates,
      structured as [neutrino species, interaction type, density, temperature, electron fraction].

    Notes:
    - Unlike integrate_table(), this function allows for a neutrino temperature which differs from the fluid temperature.
    """
    # Multiplicative factors for energy adjustments based on neutrino species
    mfac = np.array([1.0, -1.0, 0.0])
    # Factors for the Fermi-Dirac distribution
    fac = np.array([1.0, 1.0, 4.0])
    # Initialize the output table with zeros
    table = np.zeros([3, 4, len(nulib_table["rho"]), len(nulib_table["temp"]), len(nulib_table["ye"])])
    # Reference points for EOS interpolation
    points = (eos_table["ye"], eos_table["log_temp"], eos_table["log_rho"])
    # Neutrino energy array
    energy = np.array(nulib_table["energy"])
    
    # Create meshgrids for interpolation
    xv, yv, zv = np.meshgrid(nulib_table["ye"], np.log10(nulib_table["temp"]), np.log10(nulib_table["rho"]), indexing="ij")
    # Interpolate chemical potentials from the EOS table
    m1 = interpolate.interpn(points, eos_table["mu"], (xv, yv, zv))

    for n in range(0, 3):
        # Compute the Fermi-Dirac distribution eta factor
        eta = (energy[:, None, None, None] - mfac[n] * m1) / ntemp
        # Calculate the distribution for neutrino interactions
        bb = fac[n] / (1.0 + np.exp(eta)) * (energy[:, None, None, None]) ** 3
        # Integrate over energy to get the baseline interaction rates
        bb_int = np.trapz(bb, x=energy, axis=0)
        
        # Calculate absorption rates
        absorption = bb * nulib_table["absorption"][:, n, :, :, :]
        absorption_int = np.trapz(absorption, x=energy, axis=0)
        table[n, 0, :, :, :] = absorption_int.T
        table[n, 1, :, :, :] = np.divide(absorption_int, bb_int, out=np.zeros_like(absorption_int), where=bb_int != 0).T

        # Calculate emission rates
        emiss = np.zeros_like(absorption)
        for i in range(len(energy)):
            emiss[i, :, :, :] = absorption[i, :, :, :] / energy[i]
        emiss_int = np.trapz(emiss, x=energy, axis=0)
        table[n, 3, :, :, :] = emiss_int.T

        # Calculate scattering rates
        scatter = bb * nulib_table["scatter"][:, n, :, :, :]
        scatter_int = np.trapz(scatter, x=energy, axis=0)
        table[n, 2, :, :, :] = np.divide(scatter_int, bb_int, out=np.zeros_like(scatter_int), where=bb_int != 0).T

    return table


def generate_and_save_table(nulib_file: str,eos_file: str,output_file: str,
                            num_temperatures: int = 12,temperature_range: tuple = (0.1, 16)) -> None:
    """
    Generates a 4D table for grey neutrino interactions over a range of temperatures and saves it to an HDF5 file.
    The table integrates neutrino interaction data (absorption, emission, and scattering rates) over energy, 
    assuming neutrinos follow a Fermi-Dirac distribution. It consists of two parts:
    
    1. Equilibrium interactions: Assumes the fluid temperature and the neutrino temperature are the same. 
    
    2. Non-equilibrium interactions: Allows for a range of neutrino temperatures, distinct from the fluid temperature. 
       This section is generated by sampling over specified neutrino temperatures, determined by num_temperatures 
       and temperature_range, to account for scenarios where neutrino and fluid temperatures differ.
    

    Parameters:
    - nulib_file (str): Path to the NuLib table file.
    - eos_file (str): Path to the EOS table file.
    - output_file (str): Path and filename for the output HDF5 file.
    - num_temperatures (int): Number of temperature points to generate.
    - temperature_range (tuple): The minimum and maximum temperatures (in MeV) for the logspace generation.

    Returns:
    - None: The function does not return a value but writes data to an HDF5 file.
    """
    # Read NuLib and EOS tables
    nulib_table = read_nulib_table(nulib_file)
    eos_table = read_eos_table(eos_file)

    # Generate log-spaced neutrino temperature array
    log_temperatures = np.logspace(np.log10(temperature_range[0]), np.log10(temperature_range[1]), num=num_temperatures)

    # Initialize the 4D table array
    table4d = np.zeros([3, 4, 82, 100, 56, len(log_temperatures)])

    # Integrate over temperatures
    tv4 = integrate_table(nulib_table, eos_table)  # Assuming integrate_table_v4 is similar to integrate_table
    for i, tt in enumerate(log_temperatures):
        table4d[:, :, :, :, :, i] = integrate_table_n_temp(nulib_table, eos_table, tt)
        #The emissivities are always the ones "of the fluid". 
        #It is done in this double counting way because it is easier down the line.
        #Change this if you worry about the size of your table.
        table4d[:, 0, :, :, :, i] = tv4[:, 0, :, :, :] 
        table4d[:, 3, :, :, :, i] = tv4[:, 3, :, :, :]

    # Construct the table dictionary
    table_dict = {
        "rho": nulib_table['rho'], "temp": nulib_table["temp"], "ye": nulib_table["ye"],
        "tnue": log_temperatures,
        "emissivities_eq": tv4[:, 0, :, :, :], "absorption_opacity_eq": tv4[:, 1, :, :, :],
        "scattering_opacity_eq": tv4[:, 2, :, :, :], "number_emissivities_eq": tv4[:, 3, :, :, :],
        "emissivities": table4d[:, 0, :, :, :, :], "absorption_opacity": table4d[:, 1, :, :, :, :],
        "scattering_opacity": table4d[:, 2, :, :, :, :], "number_emissivities": table4d[:, 3, :, :, :, :]
    }

    # Save the table to an HDF5 file
    save_table_as_hf5f(table_dict, table_name=output_file)


#Example use
#generate_and_save_table(nulib_file = 'NuLib_SRO_0.95.h5', eos_file = 'SRO.h5',output_path = 'NuLib_SRO_0.95_grey.h5')
