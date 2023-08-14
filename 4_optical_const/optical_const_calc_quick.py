from __future__ import annotations

import glob
import time
import os

import numpy as np

"""
Code reads in file with wave number in (1/cm) and absorbance [] for lab data
and calculates the following:
    the thickness (microns)
    Refractive index (no)
    Refractive index substrate (n2)
    MAPE ('mean absolute percentage error) cut for the given spectra at each wavelength
"""

print(
    """
     *--------------------------------------------------------------*
     |           NKABS - N and K from Absorbance Spectrum           |
     |                   Version: January, 2019                     |
     |                                                              |
     |           Astrochemisty and Astrophysics Group               |
     | LASA - Laboratorio de AStroquimica e Astrobiologia da Univap |
     |         Instituto de Pesquisa & Desenvolvimento              |
     |            Universidade do Vale do Paraiba                   |
     |  Authors: Dr. Will Robson and Dr. Sergio Pilling             |
     |                                                              |
     *--------------------------------------------------------------*

"""
)


def read_file(file_name: str) -> tuple[np.ndarray, np.ndarray]:
    """read file and return first two columns (wavelength & absorbance"""
    return np.loadtxt(file_name, delimiter=" ", dtype=float, unpack=True, usecols=(0, 1))


def nk_calc(xAb: np.ndarray, yAb: np.ndarray, d: float, n0: float, n2: float, error: float, dens: float, leng: np.ndarray):
    """calculate the optical constants (n-/ k-) from the given absorbance"""

    if len(xAb) != len(yAb): #check 1
        raise IndexError(f"xAb and yAb have different lengths! {len(xAb)=} {len(yAb)=}")

    leni = len(xAb)
    #############################################################################
    # 1. CONVERT Absorbance to experimental temperature
    #############################################################################
    xTe = xAb.copy()
    yTe = (np.power(10,yAb))**(-1)	#Rocha14, eq 4 #1.0 / (10**yAb)
    np.savetxt(out + Constants.Files.T_EXP, np.column_stack((xTe, yTe)))

    #############################################################################
    # 2. Calculate the value of optical constant k (imaginary part)
    #############################################################################
    init = np.ones(shape=xAb.shape, like=xAb)
    np.savetxt(out + Constants.Files.COEF, np.column_stack((xAb, init)))

    number_of_iterations = 0
    result = 1.0
    doubl = False # to check that right ice layer thickness is used (relevant for composite ices)
    while result > error:
        # convert Coef.dat to vector
        xC, yt01 = np.loadtxt(
            out + Constants.Files.COEF,
            delimiter=" ",
            dtype=float,
            unpack=True,
        )

        alpha = (1.0 / d) * (2.3 * yAb + 2.0 * (np.log(yt01)))
        imag = alpha / (12.5 * xC)

        #############################################################################
        # 3. Calculate optical constant n
        #############################################################################
        xk = xC.copy()
        yk = imag.copy()

        # open file for write n values
        # WARNING here I assume that what he wants to do is:
        # but there is a problem. If we actually do this, then h would be smaller than all other
        # variables.
        h = abs(xk[-1]-xk[0]) / leni #]*leni #np.linspace(xk[0],xk[-1],leni)#(xk[:-1] - xk[1:]) / np.arange(1, leni)
        # WARNING until we know how h should be defined, I use this placeholder
        print(h)
        xk_odd = xk[1::2]
        yk_odd = yk[1::2]

        xk_even = xk[::2]
        yk_even = yk[::2]
        print()
        print(f"{len(xk)=}")
        print(f"{len(yk)=}")
        g = number_of_iterations
        suma = []
        for g in range(1,leni,2): #iterate over all odd indexes
            #h = [(xk[leni - 1] - xk[0]) / leni]
            
            fodd = 0.5*(((yk_even)/(xk_even-xk[g]))+((yk_even)/(xk_even+xk[g]))) #oppposite logic, if value is even you use the odd array
            #if g2 % 2 != 0:
            feven = 0.5*(((yk_odd)/(xk_odd-xk[g-1]))+((yk_odd)/(xk_odd+xk[g-1])))

            suma.append(np.sum(feven))
            suma.append(np.sum(fodd))
        #print(h, suma)
        real = n0 + np.multiply((2.0 * np.pi**(-1)) * 2*h, suma)
        if 'inf' in real:
            print('Underflow error, optical constant n is nan. Please decrease ice thickness!')
            break
        yn = real.copy()

        # for loop just writing nk-file
        path = out + 'quick5_kn_constants.lnk'
        path2 = out2 + 'quick5_kn_constants.lnk' 
        f = open(path, 'w')
        g = open(path2, 'w')

        #adjust the line length for extended set
        step1 = int(abs(xk[-1] - 1e4/leng[0]) / h)  # xk[-1] because of reverse order
        step2 = int(abs(1e4/leng[1] - xk[0]) / h)   # abs(1e4/xk[1] - 1e4/xk[0])) #!!!! to few steps

        f.write("# First line is N_lam and rho [g/cm^3], \n# calculated optical constants n (REAL part) and k (IMAGINARY part) for their respective wavelengths \n# wavelength [um] \t k \t n \n")
        f.write('  ' + str(leni-1+step1+step2) + '  '+str(dens)+'\n')
        
        g.write("# First line is N_lam and rho [g/cm^3], \n# calculated optical constants n (REAL part) and k (IMAGINARY part) for their respective wavelengths \n# wavelength [um] \t k \t n \n")
        g.write('  ' + str(leni-1+step1+step2) + '  '+str(dens)+'\n')

        if leng[0] < 1e4 / xk[-1]: #generate left side of data points
            #print(leng[0],1e4/xk[-1],step1)
            for j in np.linspace(leng[0], 1e4 / xk[-1], step1):
                #print(int(leng[0]), int(1e4 / xk[-1]), step1)
                f.write("      {0:f}     {1:f}     {2:f} \n".format(j, real[-1], imag[-1]))
                g.write("      {0:f}     {1:f}     {2:f} \n".format(j, real[-1], imag[-1]))
        for i3 in range(1,leni):
            i3 = leni-i3
            #write output values in file
            if leng[0] <= 1e4/xk[i3] <= leng[1]: #check if selected output range is smaller than calculated range
                f.write("      {0:f}     {1:f}     {2:f} \n".format(1e4 / xk[i3], real[i3], imag[i3]))
                g.write("      {0:f}     {1:f}     {2:f} \n".format(1e4 / xk[i3], real[i3], imag[i3]))

        if leng[1] > 1e4/xk[0]:  #generate right side of data points
            for l in np.linspace(1e4/xk[0],leng[1],step2):
                #print(int(1e4/xk[0]),int(leng[1]),step2)
                f.write("      {0:f}     {1:f}     {2:f} \n".format(l, real[0], imag[0]))
                g.write("      {0:f}     {1:f}     {2:f} \n".format(l, real[0], imag[0]))


        #############################################################################
        # 4. CALCULATE THEORITICAL TRANSMITTANCE
        #############################################################################

        n = yn + 1j * yk
        t01 = 2 * n0 * (n0 + n) ** (-1)
        t12 = 2 * n * (n + n2) ** (-1)
        t02 = 2 * n0 * (n0 + n2) ** (-1)
        r01 = (n0 - n) * (n0 + n) ** (-1)
        r12 = (n - n2) * (n + n2) ** (-1)
        Ex1 = np.exp(-12.5 * d * xk * yk)
        Ex2 = np.exp(2j * (6.28 * xk * d * n))
        Trans = Ex1 * ((abs((t01 * t12 / t02) / (1.0 + r01 * r12 * Ex2))) ** 2)
        factor = abs((t01 * t12 / t02) * (1.0 + r01 * r12 * Ex2) ** (-1))

        np.savetxt(out + Constants.Files.T_TEO, np.column_stack((xk, Ex1, n, Trans)))
        np.savetxt(out + Constants.Files.COEF, np.column_stack((xk, factor)))

        # Convert Texp and Tteo to vector
        yTt = Trans.copy()

        #############################################################################
        # 5. CALCULATE CHI-SQUARE
        #############################################################################

        mape = np.abs((1.0 * leni**(-1)) * (yTe - yTt) / yTt) #Rocha14: eq 15
        chi = ((yTe - yTt) ** 2) / yTt

        result = np.sum(mape)
        if result == 'nan' and doub == False: #to correct for wrong input because of single counting the composite of ice
            d *= 2
            result = 1.0
            doub = True
            print('Increasing thinckness to: '+str(d))
        result1 = np.sum(chi)
        number_of_iterations += 1

    end_time = time.time()

    dt = end_time - start_time

    # INFO the number_of_iterations - 1 follows the convention of the previouse author
    print(
        f"""
======================================
             LIST OF RESULTS
The MAPE was:                 | {result * 100}%
The chi-square was:           | {result1}
The elapsed time was          | {dt:.2f} sec
with the number of iterations | {number_of_iterations - 1}
======================================
"""
    )

#%%%%%%%%%%%%%%%%%%%% personalize code from here on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Constants:
    """stores constants"""

    class Files:
        """different filenames"""

        COEF = "COEF.txt"
        KN_CONSTANTS = "kn_constants.lnk"
        T_EXP = "T_EXP.txt"
        T_TEO = "T_TEO.txt"
        K = "K.txt"

        # ! within your path
        LAB_DATA_FOLDER = "../lab_data"
        RESULTS_FOLDER = "./results"
        NEXT_STEP = "../5_optool/lnk_data"
        print(f"{LAB_DATA_FOLDER=}")
        print(f"{RESULTS_FOLDER=}")

print(
    "Please insert the thickness (microns), "
    "Refractive index (no), Refractive index substrate (n2) and "
    "MAPE ('mean absolute percentage error) cut!",
    end="\n\n",
)

print("In case of non-float / non-integer input: default values are used.", end="\n\n")

def input_data(inp_file):
    try:
        d = float(input("\n Thickness d (microns) (For example: 0.4):")) * 1e-4  # convert d to cm
        n0 = float(input("\n n0 (e.g. 1.22):"))
        n2 = float(input("\n n2 (Use e.g. 1.73 for CsI, 1.54 for KBr and 2.54 for ZnSe):"))
        error = float(input("\n MAPE (For example: 0.1):")) / 100  # to obtain fraction from percentage
    except ValueError:
        d = 2.5 * 1e-4 # converte parameter to cm
        n0 = 1.16
        n2 = 1.73
        error = 0.1 * 1e-3  # /100
        if 'CO_' in inp_file:
           n0 = 1.22
        elif 'H2O' in inp_file:
               d = 0.5 * 1e-4 #0.3 * 1e-4
    print('For the output file the ice density (default: 0.47 um) and the length of the wave number range in 1/cm (default: 1, 30) is required.')
    dens = 0.47
    leng = [1.,30.]
    print(
	f"""
got error: ValueError
the following default values will be used:
     {d=}
     {n0=}
     {n2=}
     MAPE {error=}
"""
    )
    return d,n0,n2,error,dens, leng
    
    
# Define own input file in .txt-format (wave number, absorbance)
# check all given .txt-files in directory and choose the first one
files = glob.glob(f"{Constants.Files.LAB_DATA_FOLDER}/*.txt")
print(f"Found the following files in the lab data folder: {files}")

# iterate through all files
for in_file in files:
    print(f"{in_file=}")
    d,n0,n2, error, dens, leng = input_data(in_file) #define input data

    # path where output data should be stored
    out_path = f"{Constants.Files.RESULTS_FOLDER}/d{d}_n0{n0}_n2{n2}_err{error}_"
    next_path= f"{Constants.Files.NEXT_STEP}/d{d}_n0{n0}_n2{n2}_err{error}_"
    start_time = time.time()  # to check timing
    xAb, yAb = read_file(in_file)
    # WARNING: xAb and yAb contain 1 value too much otherwise the calculation of h will go wrong
	# one we know what h is supposed to be this should be looked at
	# to avoid that len(h) != len(odd) + len(even), remove last value of when len(odd) != len(even)
    if len(xAb) % 2 != 0:
        xAb, yAb = xAb[:-1], yAb[:-1]

    # calculate n-/k- coefficients via the following routine
    out = out_path + os.path.basename(in_file)[:-4] + "_"
    out2 = next_path + os.path.basename(in_file)[:-4] + "_"
    nk_calc(xAb, yAb, d, n0, n2, error,dens, leng)
