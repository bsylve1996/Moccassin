import numpy as np
import pydl as pd
import math
import matplotlib.pyplot as plt
import savReaderWriter as srw
import statistics
import scipy.stats as sci
import SUPER_AMIY as am
import plot_dust_temp_amiy as pdta
import SUPER_KELSEY as sk
import os
import idlsave
from scipy import interpolate

def rebin(a, new_shape):
    M, N = a.shape
    m, n = new_shape
    if m < M:
        return a.reshape((m, M / m, n, N / n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m / M, axis = 0), n / N, axis = 1)

def planck(T, wave):
    """Planck function in wavelength.

    :INPUTS:
      T : scalar or array
        temperature in Kelvin

      wave : scalar or array
        wavelength in microns

    Value returned is in (nearly) cgs units: erg/s/cm^2/micron/sr
    """

    c = 29979245800. # cm/s
    nu = c/(wave/1e4) # Hz

    # Take advantage of the fact that (lam F_lam) = (nu F_nu):

    return (bnu(T, wave) * (nu / wave))

def bnu(T, wave):
    """Planck function in frequency.

    :INPUTS:
      T : scalar or array
        temperature in Kelvin

      wave : scalar or array
        wavelength in microns [but intensity will be per Hz]

    Value returned is in cgs units: erg/s/cm^2/Hz/sr
    """
    from numpy import exp

    c = 29979245800. # cm/s
    nu = c/(wave/1e4)
    h = 6.626068e-27 # cgs units
    k =  1.3806503e-16 # cgs units
    expo = h*nu/(k*T)
    nuoverc = 1./ (wave/1e4)
    return ((2*h*nuoverc**2 * nu)) / (exp(expo)-1)

def idl_tabulate(x, f, p=5) :
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        import scipy as sc
        weights = sc.integrate.newton_cotes(rn)[0]
        return (x[-1] - x[0]) / (x.shape[0] - 1) * np.dot(weights, f)
    ret = 0
    for idx in range(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret

def testbb(lum,temp,rstar,errorcheck="errorcheck"): #rstar is distance to object in pc
    if errorcheck.equals("errorcheck"):
        print,'testbbing'
    pc=3.086e18                  #1 parsec in cm
    const = 1.e36/pc^2


    npointsArray = [x for x in range(50000)] #npoints = 50000
    wave = []
    for x in npointsArray:
        wave.append((x+1)/50000*350000.0)# 50000 = npoints
    #wavelength in angstroms

    bbflux = planck(temp, wave)
    bbflux = float(bbflux)
    integral = idl_tabulate(wave, bbflux)
    bbflux = bbflux * (lum/integral)
    wave = wave/1.0e4 #convert to microns
    if errorcheck.equals("errorcheck"):
        print(integral)

    #stop
    yfinal = (bbflux/3e-13)*((wave)^2)
    yfinal = yfinal * const * math.pi / (4. * math.pi * rstar^2.)
    plt.plot(wave,yfinal/math.pi*1.0e3)

    temporary_var = []
    maxYFinal = max(yfinal)
    for i in yfinal:
        if (yfinal[i] == maxYFinal):
            temporary_var.append(yfinal[i])
    print("max flux at a wavelength of " + yfinal[temporary_var])
    print('corresponding to a TEST BB temp of  ' + 2897/yfinal[temporary_var])

def chi_squared(fObserved, fExpected, errors, p):
    if len(fObserved) != len(fExpected):
        print('interpolation failed')
    chi_sqr = 0.0
    for i in fObserved:
        chi_sqr += (fObserved[i] - fExpected[i]) ^ 2 / (errors[i]) ^ 2
    chi_sqr = chi_sqr / (float(len(fObserved)-p))
    return chi_sqr


def star_dust_mass(mass, chisqr, starname, runnum, username):
    name, dustmass, chisquare, id = np.loadtxt('/Users/' + username + '/mocassin-rw_changes/stardustmass.txt', unpack=True)

    i = 0
    starnum = -1

    while (i < len(name) and starnum == -1):
        if (name[i] == starname):
            starnum=i
        i += 1
    if (starnum >= 0):
        if (chisqr < chisquare[starnum] and chisqr >= 1):
            dustmass[starnum] = mass
            chisquare[starnum] = chisqr
            id[starnum] = runnum
    DAT = np.column_stack(name, dustmass, chisquare, id)
    np.savetxt('/Users/' + username + '/mocassin-rw_changes/stardustmass.txt', DAT, delimiter = " ", fmt = "%s")


def extra_moc_plot(flux_input, lam_input, flux_moc, lam_moc, error_flux_input, chi_sqr):
    srw.SavWriter('extra_moc_plot.sav',flux_input, lam_input, flux_moc, lam_moc, error_flux_input)

    #check input
    if min(lam_input) < 1: print('input fail2')
    if max(lam_input) > 45: print('input fail3')

#restrict the maximum and minimum values to be that of the IR
    minlam = min(lam_input)
    maxlam = max(lam_input)

    index = []
    for i in lam_moc:
        if (lam_moc[i] >= minlam and lam_moc[i] <= maxlam):
            index.append(lam_moc[i])

    flux_moc = flux_moc[index]
    lam_moc = lam_moc[index]

#interpolate and subtract the mocassin data from the input



    # Specifies the kind of interpolation as a string (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’ where
    # ‘zero’, ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline interpolation of zeroth, first, second or third order) or
    # as an integer specifying the order of the spline interpolator to use. Default is ‘linear’.
    interpfunc = interpolate.interp1d(lam_moc, flux_moc, kind='linear')
    flux_moc_int = interpfunc(lam_input) #idl equivalent flux_moc_int = interpol(flux_moc, lam_moc, lam_input)

#moc to many points.
    flux_subtracted = []
    for i in flux_input:
        flux_subtracted.append(flux_input[i] - flux_moc_int)

#calc and print useful things
    median = statistics.median(flux_subtracted)
    result =[np.mean(flux_subtracted),np.var(flux_subtracted),sci.skew(flux_subtracted),sci.kurtosis(flux_subtracted)]
    std_dev = np.std(flux_subtracted)
    chi_sqr = chi_squared(flux_input, flux_moc_int, error_flux_input, 10)


    plt.subplot(221)
    plt.gca().set_color_cycle(['red', 'black', 'green', 'blue','green', 'blue', 'black'])
#plot and stuff x, y and plot a line a y = 0 for refrence
    plt.plot(lam_input, flux_subtracted,'D', [.1,100], [0,0], linewidth = 1)
    plt.title = 'Real - Moc'
    plt.errorbar(lam_input, flux_subtracted, error_flux_input, 3)

    plt.subplot(222)
    plt.plot(lam_input, flux_input, '^', lam_moc, flux_moc, 'D',[.1, 100], [0, 0], linestyle = 1)

    plt.subplot(223)
    plt.plot(lam_input, flux_input, 's', lam_moc, flux_moc, 'D', [.1, 100], [0, 0], linestyle = 1)

def sm(tab):
    n_interm = len(np.reshape(tab[* 0, 0]))*5
    return(pd.smooth(rebin(np.reshape(tab[* 0, 0]), n_interm), 5))

def rerun_mocplot_amiy(username):
    #rerun_mocplot_amiy, 'bill'
    am.cd_amiy('/Users/' + username + '/mocassin-rw_changes/output')
    s = idlsave.read('mocplot_variables.sav')
    #mocplot_Clumpy(infile,paramfile)
    print('done plotting again!!!')

def mocplot_Clumpy(infile, paramfile, residual = "residual", errorcheck = "errorcheck", includePAHS = "includePAHS"):

    if (errorcheck == "errorcheck"): print('starting mocplot amiy 2')

    username = 'mocassin'
    starnum = -1

    starname, rin, rout, rho, cov_fac, diffuse, tstellar, lum, sil, carb, distance,donut, nova, theta1, phi1, theta2, phi2 = np.loadtxt(paramfile, unpack=True)

    am.cd_amiy('/Users/' + username + '/mocassin-rw_changes/output')

    theta1 = theta1 * (180 /math.pi)
    theta2 = theta2 * (180 /math.pi)

    constant = 0 #i dont know what constant to use

#COPY FILES OVER
    if diffuse:
        directoryname = '/Users/mocassin/mocassin-rw_changes/output/SN/' + starname + '_clumpy_' + 'f_' + str(cov_fac).strip() \
                        + '_Rin' + str(rin).strip() + '_Rout' + str(rout).strip() + '_Ls' + str(lum).strip() + '_Rho' + \
                        str(rho).strip() + '_sil' + str(math.floor(sil)).strip() + '_carb' + str(math.floor(carb)).strip()
    else:
        directoryname = 'RSG/' + id + '_' + starname + 'Co_' + str(math.floor(constant)).strip() + '_Rin' + str(rin).strip() \
                        + '_Rout' + str(rout).strip() + '_Ls' + str(lum).strip() + '_Rho' + str(math.floor(rho)).strip() \
                        + '_sil' + str(math.floor(sil)).strip() + '_carb' + str(math.floor(carb)).strip()

    outfoldername = directoryname

    os.system('mkdir ' + directoryname)
    os.system('cp dustGrid.out ' + directoryname + '/dustGrid.out.txt')

    if diffuse: os.system('cp equivalentTau.out '+directoryname+'/equivalentTau.out.txt')
    else: os.system('cp tauNu.out '+directoryname+'/tauNu_'+id+'.out.txt')

    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/input.in ' + directoryname + '/input.in.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/grainspecies.dat ' + directoryname + '/grainspecies.dat.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/ndust/nDUST ' + directoryname + '/nDUST.in.txt')


# ** ** ** ** ** MODIFY BELOW ** ** ** ** ** ** ** *
    mocassinfile = 'SED.out' #filename
    n = nova #inclination angles
    D = distance#3.36e6 #distance in pc

# ** ** ** ** ** END MODIFICATION ** ** ** ** ** *
    if (errorcheck == "errorcheck"): print(' SPOT 1')

#;GET THE MASS
    testarray1 = []
    os.system('grep "Total" dustGrid.out > test')
    lun1 = open('test', 'r')
    testarray1 = lun1.read()
    test = testarray1(0).index(':')#if does not exist code will throw a value error
    mass = float((float(testarray1(0)[test + 4, 13]) * 1.e30 * 1.e15) / 1.98892e33)
    lun1.close()

    os.system('rm -f test')

    sil = sk.sss(sil)
    carb = sk.sss(carb)

    if (errorcheck == "errorcheck" ): print('SPOT 1.5')

#TWO DIFFERENT WAYSTO CALCULATE TAU

#ONE FOR SN AND ONE FOR RSG
    if diffuse:
        print('Trying to set value of Tau')
        os.system('grep "      0.5623" /Users/mocassin/mocassin-rw_changes/output/equivalentTau.out > /Users/mocassin/mocassin-rw_changes/output/test.dat')

        a, wave, tau, d = np.loadtxt('/Users/mocassin/mocassin-rw_changes/output/equivalentTau.out', unpack=True)
        os.system('rm -f /Users/mocassin/mocassin-rw_changes/output/test.dat')
#;stop
    else:
        os.system('grep "0.5623" ./tauNu.out > test')
        lun1 = open('test', 'r')
        waveRSG = [] # because there are 3 coordinates
        tauRSG = [] # and each one makes a tau in that direction...
        waveRSG, tauRSG = np.loadtxt('test', unpack=True)
        wave = waveRSG[0]
        tau = tauRSG[0]
        lun1.close()
        os.system('rm -f test')

    if (errorcheck == "errorcheck"): print('SPOT 2')

    if (n == 0): col0, col1, col2 = np.loadtxt(mocassinfile, unpack=True)
    elif (n == 1): col0, col1, col2, col3 = np.loadtxt(mocassinfile, unpack=True)
    elif (n == 2): col0, col1, col2, col3, col4 = np.loadtxt(mocassinfile, unpack=True)

    data = [[0.0 for x in range(len(col1))]for x in range(3 + n)]#n is number of viewing angles
    for i in range(len(col0)):
        data[0,i] = col0[i]
    for i in range(len(col1)):
        data[1,i] = col1[i]
    for i in range(len(col2)):
        data[2,i] = col2[i]

    if (n > 0):
        for i in range(len(col3)):
            data[3, i] = col3[i]
    if (n == 2):
        for i in range(len(col4)):
            data[4, i] = col4[i]

    visibledatacheck = 0
    irsdatacheck = 0
    jhkcheck = 0
    temp = 'str'




#stop

    str = "mu"
    str1 = "lambda"
    str2 = "nu"

    nu1, fnu1, unu1, dataFlag = np.loadtxt(infile, unpack=True)
    upFlag = []
    datFl = 0
    for i in dataFlag:
        if(dataFlag[i] == 0):
            upFlag.append(dataFlag[i])
            datFl += 1
    if (upFlag[0] != -1):
        upNu = nu1[upFlag]
        upFnu = fnu1[upFlag]

#Check to make sure these are not 1D arrays
#Will not plot 1D!
        if len(upNu) < 2:
            upNu = [upNu, upNu]
            upFnu = [upFnu, upFnu]
#stop
    ymin = min(fnu1) - 0.35 * min(fnu1)
    ymax = max(fnu1) + 0.25 * max(fnu1)#EJM August 31, 2014
#stop
    xmax = 30
#PLOTS IN mJy CURRENTLY!!!
    D = distance #3.36e6
    plt.plot([0, 0], [0, 0])
    plt.axis([0.3, xmax, 0.002, 20.])
    plt.xlabel(str1 + ' [' + str + 'm]')
    plt.ylabel('F' + str2 + ' [mJy]')

    if (n > 0):
        y0dots = []
        for i in range(len(data[3])):
            y0dots.append(data[3][0]/(D^2)*1e3)
        plt.plot(data[1], y0dots, 'b-', linewidth = 3)
        if (n > 1):
            y1dots = []
            for i in range(len(data[4])):
                y0dots.append(data[4][0] / (D ^ 2) * 1e3)
            plt.plot(data[1], y1dots, 'g-', linewidth = 3)

    yndots = []
    for i in range(len(data[2])):
        yndots.append(data[2][n]/(4.9e7 ^ 2) * 1000)
    plt.plot(data[1], yndots, 'b')

    print('min = ' + min(data[2]) +'        max = '+ max(data[2]))

    nu1 = nu1[datFl]
    fnu1 = fnu1[datFl]
    unu1 = unu1[datFl]
    plt.plot(nu1, fnu1)
    plt.errorbar(nu1, fnu1, unu1)

    if (upFlag[0] != -1): plt.plot(upNu, upFnu, 'X')

    testbb(lum, tstellar, D)

    name2, vis, ir = np.loadtxt('/Users/mocassin/mocassin-rw_changes/accessories/staceyplots/lmcsmcfullspec2e.txt', unpack=True)

#read in a file with all of the SED data in a list called spectralist.txt
    if (starnum >= 0):
       if (os.path.isfile('/Users/mocassin/mocassin-rw_changes/accessories/staceyplots/IRspectra/' + ir[starnum])): visibledatacheck=1
    if (os.path.isfile('/Users/mocassin/mocassin-rw_changes/accessories/plots/text/' + starname + 'visibledata.txt')): irsdatacheck=1
    if (os.path.isfile('/Users/mocassin/mocassin-rw_changes/accessories/plots/text/' + starname + 'data.txt')): jhkcheck=1

#read in data to be plotted from RSG data.
    if (irsdatacheck): nu1, fnu1, irs_e_wave, irs_e_flux = np.loadtxt('/Users/mocassin/mocassin-rw_changes/accessories/staceyplots/IRspectra/'+ir[starnum], unpack=True)
    if (visibledatacheck): w, f = np.loadtxt('/Users/mocassin/mocassin-rw_changes/accessories/plots/text/'+starname+'visibledata.txt', unpack=True)
    if (jhkcheck): wdata, fdata = np.loadtxt('/Users/mocassin/mocassin-rw_changes/accessories/plots/text/'+starname+'data.txt', unpack=True)


    if (jhkcheck == 0):
        if os.path.isfile('/Users/' + username + '/rsgs/plots/text/' + starname + 'irsdata.txt'): visibledatacheck=1
        if os.path.isfile('/Users/' + username + '/rsgs/plots/text/' + starname + 'visibledata.txt'): irsdatacheck=1
        if os.path.isfile('/Users/' + username + '/rsgs/plots/text/' + starname + 'data.txt'): jhkcheck=1
        if irsdatacheck: nu1, fnu1, irs_e_wave, irs_e_flux = np.loadtxt('/Users/'+username+'/rsgs/plots/text/'+starname+'irsdata.txt', unpack=True)
        if visibledatacheck: w, f = np.loadtxt('/Users/'+username+'/rsgs/plots/text/'+starname+'visibledata.txt', unpack=True)
        if jhkcheck: wdata, fdata = np.loadtxt('/Users/'+username+'/rsgs/plots/text/'+starname+'data.txt', unpack=True)

    if visibledatacheck: plt.plot(w, f, 'r')#plot visible
    if irsdatacheck: plt.plot(nu1, fnu1, 'r')#irs
    if jhkcheck: plt.plot(wdata, fdata, 'rD')#jhk, iras mips data

    stringy = 'Covering Factor =  ' + cov_fac +'\nrin =\t' + rin + ' cm\nrout =\t'
    stringy += rout + ' cm\nL =\t' + lum + ' e36 erg/s\nT =\t' + tstellar.floor() + ' K'
    stringy += 'oss ' + sk.ssi(sil) + ' amC ' + sk.ssi(carb) + '\nrho = ' + rho.floor() + '\nDMass=' + mass

    for i in wave:
        if (wave[i] > 0.51 and wave[i] <= 0.55):
            idx = i
            break
    #idx = where(wave > 0.51 and wave <= 0.55)
    stringy += "\ntau = " + tau[idx] + ' at ' + wave[idx] + ' um'

#    if (nova > 0):
#        if (nova == 1):
#            if (theta1 == 0):
#                stringy += '\n' + thetaA + ' = ' + thetaangle1 + degree + ' , ' + phiA + ' = '+phiangle1 + degree + '. Face on'
#            else:
#                if (theta1 == 89):
#                    stringy += '\n' + thetaA + ' = ' + thetaangle1 + degree + ' , ' + phiA + ' = '+phiangle1 + degree + '. Edge on'
#        else:
#            if(nova == 2):
#                stringy += '\n' + thetaA + ' = ' + thetaangle1 + degree + ' , ' + phiA + ' = '+phiangle1 + degree
#                stringy += '\n' + thetaB + ' = ' + thetaangle2 + degree + ' , ' + phiB + ' = '+phiangle2 + degree

    if (donut == 0):
        stringy += 'Shape = Shell'
    else:
        stringy += 'Shape = Torus'

    plt.figtext(0.5, 0, stringy, wrap=True, horizontalalignment='center',fontsize=12)
    plt.savefig('/Users/' + username + '/mocassin-rw_changes/output/sn_smooth.eps')
#stop
    pdta.plot_dust_temp_amiy(username, diffuse)
    am.cd_amiy('/Users/' + username + '/mocassin-rw_changes/output')
    os.system('cp sn_smooth.eps ' + directoryname + '/' + starname + '_clumpy.eps')