import numpy as np
import pydl as pd
import math
import matplotlib.pyplot as plt
import savReaderWriter as srw
import statistics
import scipy.stats as sci
import plot_dust_temp_kelsey as pdtk
import SUPER_KELSEY as sk
import os
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline


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
    nu = []
    for i in range(len(wave)):
        nu.append(c/(wave[i]/1e4)) # Hz

    # Take advantage of the fact that (lam F_lam) = (nu F_nu):
    returnThis = []
    bnuArr = bnu(T, wave)
    for i in range(len(bnuArr)):
        returnThis.append(float(bnuArr[i] * (nu[i]/wave[i])))
    return returnThis

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
    nu = []
    for i in range(len(wave)):
        nu.append(c/(wave[i]/1e4))
    h = 6.626068e-27 # cgs units
    k =  1.3806503e-16 # cgs units
    expo = []
    for i in range(len(nu)):
        expo.append(h*nu[i]/(k*T))
    nuoverc = []
    for i in range(len(wave)):
        nuoverc.append(1./ (wave[i]/1e4))
    returningThis = []
    for i in range(len(nuoverc)):
        returningThis.append((2*h*(float(nuoverc[i])**2) * float(nu[i])) / (float(exp(float(expo[i])))-1))
    return returningThis

def testbb(lum,temp,rstar,errorcheck="errorcheck"): #rstar is distance to object in pc
    if errorcheck == "errorcheck":
        print('testbbing')
    pc=3.086e18                  #1 parsec in cm
    const = 1.e36/ pc ** 2


    npointsArray = [x for x in range(50000)] #npoints = 50000
    wave = []
    for x in npointsArray:
        wave.append((x+1)/50000*350000.0)# 50000 = npoints
    #wavelength in angstroms

    bbflux = planck(temp, wave)
    wave_x = np.array(wave)
    bbflux_y = np.array(bbflux)
    integral = InterpolatedUnivariateSpline(wave_x, bbflux_y, k = 5).integral(min(wave),max(wave))
    for i in range(len(bbflux)):
        bbflux[i] = float(bbflux[i] * (lum/integral))
    for i in range(len(wave)):
        wave[i] = wave[i]/1.0e4 #convert to microns
    if errorcheck == "errorcheck":
        print("integral is" + str(integral))

    #stop
    yfinal = []
    for i in range(len(bbflux)):
        yfinal.append(((bbflux[i]/3e-13)*((wave[i])**2))* const * math.pi / (4. * math.pi * rstar**2.))

    graphing_yfinal = []
    for i in range(len(yfinal)):
        graphing_yfinal.append(yfinal[i]/math.pi*1.0e3)
    plt.plot(wave,graphing_yfinal)

    maxYFinal = max(yfinal)
    for i in range(len(yfinal)):
        if (yfinal[i] == maxYFinal):
            temporary_var = i
    print("max flux at a wavelength of " + str(yfinal[temporary_var]))
    print('corresponding to a TEST BB temp of  ' + str(2897/yfinal[temporary_var]))

def chi_squared(fObserved, fExpected, errors, p):
    if len(fObserved) != len(fExpected):
        print('interpolation failed')
    chi_sqr = 0.0
    for i in fObserved:
        chi_sqr += (fObserved[i] - fExpected[i]) ^ 2 / (errors[i]) ^ 2
    chi_sqr = chi_sqr / (float(len(fObserved)-p))
    return chi_sqr


def star_dust_mass(mass, chisqr, starname, runnum, username):
    name = np.loadtxt('/Users/' + username + '/mocassin-rw_changes/stardustmass.txt', unpack=True, dtype='str', usecols=(0))
    dustmass, chisquare, id = np.loadtxt('/Users/' + username + '/mocassin-rw_changes/stardustmass.txt', unpack=True, usecols=(1,2,3))

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
    #srw.SavWriter('extra_moc_plot.sav',flux_input, lam_input, flux_moc, lam_moc, error_flux_input)
    f = open('extra_moc_plot.txt', 'w')
    f.write('flux_input = ' + str(flux_input))
    f.write('lam_input = ' + str(lam_input))
    f.write('flux_moc = ' + str(flux_moc))
    f.write('lam_moc = ' + str(lam_moc))
    f.write('error_flux_input = ' + str(error_flux_input))
    f.write('chi_sqr = ' + str(chi_sqr))
    f.close()

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

#def sm(tab):
#    n_interm = len(np.reshape(tab[* 0, 0]))*5
#    return(pd.smooth(rebin(np.reshape(tab[* 0, 0]), n_interm), 5))

def rerun_mocplot_kelsey(username):
    sk.cd_kelsey('/Users/'+username+'/mocassin-rw_changes/output')
    #restore,'mocplot_variables.sav'
    #mocplot_KELSEY(rin,rout,rho,lum,tstellar,starname,diffuse,username,outfoldername,starnum,donut,sil,carb,name_gs,filename_gs,percent_gs,totaltime)
    print('done plotting again!!!')

def mocplot_KELSEY(rin, rout, rho, lum, tstellar, starname, diffuse, distance, symmetric,username, outfoldername, starnum, donut, sil, carb, name_gs, filename_gs, percent_gs, totaltime,infile, nduststr, distributionFile, residual = "residual", errorcheck = "errorcheck", includePAHS = "includePAHS"):

    if (errorcheck == "errorcheck"): print('starting mocplot kelsey')

    sk.cd_kelsey('/Users/' + username + '/mocassin-rw_changes/output')
    #srw.SavWriter('mocplot_variables.sav',rin, rout, rho, lum, tstellar, starname, diffuse, distance, symmetric,username, outfoldername, starnum, donut, sil, carb, name_gs, filename_gs, percent_gs, totaltime,infile, nduststr, distributionFile)
    f = open('mocplot_variables.txt', 'w')
    f.write("rin = " + str(rin))
    f.write("rout = " + str(rout))
    f.write("rho = " + str(rho))
    f.write("lum = " + str(lum))
    f.write("tstellar = " + str(tstellar))
    f.write("starname = " + str(starname))
    f.write("diffuse = " + str (diffuse))
    f.write("distance = " + str(distance))
    f.write("symmetric = " + str(symmetric))
    f.write("username = " + str(username))
    f.write("outfoldername = " + str(outfoldername))
    f.write("starnum = " + str(starnum))
    f.write("donut = " + str(donut))
    f.write("sil = "  + str(sil))
    f.write("carb = " + str(carb))
    f.write("name_gs = " + str(name_gs))
    f.write("filename_gs = " + str(filename_gs))
    f.write("perscent_gs = " + str(percent_gs))
    f.write("total time = " + str(totaltime))
    f.write("infile = " + str(infile))
    f.write("nduststr = " + str(nduststr))
    f.write("distributionFile = " + str(distributionFile))
    numbers = []
    with open('/Users/' + username + '/mocassin-rw_changes/KELSEY_number.txt') as r:
        for line in r:
            number_str = line.split()
            numbers = [int(x) for x in number_str]
    if not numbers:
        KELSEY_number = 1
    else:
        KELSEY_number = int(numbers[0]) + 1
    id = int(KELSEY_number)
    r.close()
    writer = open('/Users/' + username + '/mocassin-rw_changes/KELSEY_number.txt', 'w')
    writer.write(str(KELSEY_number))
    writer.close()
    sil = np.floor(sil)
    rho2 = np.floor(rho)
    carb2 = np.floor(carb)
#COPY FILES OVER
    if diffuse:
        directoryname = 'SN/' + str(id) + '_' + str(starname) + '_Rin' + str(rin).strip() + '_Rout' + str(rout).strip() + '_Ls' + \
                        str(lum).strip() + '_Rho' + str(rho).strip() + '_sil' + str(sil).strip() + '_carb' + \
                        str(carb2).strip()
    else:
        directoryname = 'RSG/' + str(id) + '_' + str(starname) + '_Rin' + str(rin).strip() + '_Rout' + str(rout).strip() + '_Ls' + \
                        str(lum).strip() + '_Rho' + str(rho2).strip() + '_sil' + str(sil).strip() + \
                        '_carb' + str(carb2).strip()

    outfoldername = directoryname

    os.system('mkdir ' + directoryname)
    os.system('cp dustGrid.out ' + directoryname + '/dustGrid_' + str(id) + '.out.txt')
    os.system('mv runinfo.txt ' + directoryname + '/runinfo_' + str(id) + '.txt')
    os.system('cp SED.out ' + directoryname + '/SED_' + str(id) + '.out.txt')
    os.system('cp equivalentTau.out ' + directoryname + '/equivalentTau_' + str(id) + '.out.txt')
    os.system('cp tauNu.out ' + directoryname + '/tauNu_' + str(id) + '.out.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/input.in ' + directoryname + '/input_' + str(id) + '.in.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/grainspecies.dat ' + directoryname + '/grainspecies_' + str(id) + '.dat.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/' + str(nduststr) + ' ' + directoryname + '/nDUST_' + str(id) + '.in.txt')

# ** ** ** ** ** MODIFY BELOW ** ** ** ** ** ** ** *
    mocassinfile = 'SED.out' #file name
    n = 0 #Inclination set by grid file from Roger's website
#so n does not need to set by defunct nova keyword
#hard coded to be 0 - -EJM Dec 7, 2015
    D = distance #distance in pc
# ** ** ** ** ** END MODIFICATION ** ** ** ** ** *
    if (errorcheck == "errorcheck"): print(' SPOT 1')

#GET THE MASS
    testarray1 = []
    os.system('grep "Total" dustGrid.out > test')#keep this make it total instead of test and remove stringy everywhere and stupid extras
    #os.system('rm -f test')

    sil = sk.sss(sil)
    carb = sk.sss(carb)

    if (errorcheck == "errorcheck"): print('SPOT 1.5')

    if (errorcheck == "errorcheck"): print('SPOT 2')

    #col0, col1, col2 = np.loadtxt(mocassinfile, unpack=True)
    col0 = []
    col1 = []
    col2 = []

    with open(mocassinfile, "r") as f:
        for i in range(4):
            f.readline()#ignore beginning
        data = f.readlines()
        for line in data:
            try:
                arr = line.split()
                if isinstance(float(arr[0]), float):
                    col0.append(float(arr[0]))
                    col1.append(float(arr[1]))
                    col2.append(float(arr[2]))
            except:
                pass
  
    
    data = [col0, col1, col2]
    visibledatacheck = 0
    irsdatacheck = 0
    jhkcheck = 0
    temp = 'str'

    str0 = "mu"
    str1 = "lambda"
    str2 = "nu"
    
    nu1 = []
    fnu1 = []
    unu1 = []
    dataFlag = []

    with open(infile, "r") as f:
        allLines = f.readlines()
        for line in allLines:
            arr = line.split()
            if len(arr) == 0:
                break
            nu1.append(float(arr[0]))
            fnu1.append(float(arr[1]))
            unu1.append(float(arr[2]))
            dataFlag.append(int(arr[3]))
    upFlag = []
    datFl = 0
    for i in range(len(dataFlag)):
        if(dataFlag[i] == 0):
            upFlag.append(dataFlag[i])
            datFl += 1
    if (len(upFlag) != 0 and upFlag[0] != -1):
        upNu = nu1[upFlag]
        upFnu = fnu1[upFlag]
#Check to make sure these are not 1D arrays
#Will not plot 1D!
        if len(upNu) < 2:
            upNu = [upNu, upNu]
            upFnu = [upFnu, upFnu]


    ymin = min(fnu1) - 0.35 * min(fnu1)
    ymax = max(fnu1) + 0.25 * max(fnu1) #EJM August 31, 2014
    xmax = 30
#PLOTS IN mJy CURRENTLY!!!

    D = distance

    plt.plot([0, 0], [0, 0])
    plt.axis([0.3,xmax, ymin, ymax])
    plt.xlabel(str1 + ' [' + str0 + 'm]')
    plt.ylabel('F' + str2 + ' [mJy]')

    yndots = []
    for i in range(len(data[2])):
        yndots.append(float(data[2][n]) / (4.9e7 ** 2) * 1000)
    plt.plot(np.array(data[1]), np.array(yndots), 'b')

    print('min = ' + str(min(data[2])) + '        max = ' + str(max(data[2])))

    nu1 = nu1[datFl]
    fnu1 = fnu1[datFl]
    unu1 = unu1[datFl]
    plt.plot(nu1, fnu1)
    plt.errorbar(nu1, fnu1, unu1)

    if (len(upFlag) != 0 and upFlag[0] != -1):
        plt.plot(upNu, upFnu,marker = 'X')

    testbb(lum, tstellar, D)

    irs = np.loadtxt('/Users/mocassin/dust_modelling/accessories/staceyplots/lmcsmcfullspec2e.txt', unpack=True, usecols=(2), skiprows=1)

#read in a file with all of the SED data in a list called spectralist.txt
    if (starnum >= 0):
        if (os.path.isfile('/Users/mocassin/dust_modelling/accessories/staceyplots/IRspectra/' + str(irs[starnum]))): visibledatacheck=1
    if (os.path.isfile('/Users/mocassin/dust_modelling/accessories/plots/text/' + str(starname) + 'visibledata.txt')): irsdatacheck=1
    if (os.path.isfile('/Users/mocassin/dust_modelling/accessories/plots/text/' + str(starname) + 'data.txt')): jhkcheck=1

#read in data to be plotted from RSG data.
    if (irsdatacheck): nu1, fnu1, irs_e_wave, irs_e_flux = np.loadtxt('/Users/mocassin/dust_modelling/accessories/staceyplots/IRspectra/'+ str(irs[starnum]), unpack=True)
    if (visibledatacheck):  w, f = np.loadtxt('/Users/mocassin/dust_modelling/accessories/plots/text/' + str(starname) + 'visibledata.txt', unpack=True)
    if (jhkcheck):  wdata, fdata = np.loadtxt('/Users/mocassin/dust_modelling/accessories/plots/text/' + str(starname) + 'data.txt', unpack=True)

#if mocplot cant find the file from bill, then it will search in the user's directory... IF IT IS THERE!!!
    if (jhkcheck == 0):#if the JHK data is not in bill's directory, search user's
        if (os.path.isfile('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'irsdata.txt')): visibledatacheck=1
        if (os.path.isfile('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'visibledata.txt')): irsdatacheck=1
        if (os.path.isfile('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'data.txt')): jhkcheck=1
        if (irsdatacheck):  nu1, fnu1, irs_e_wave, irs_e_flux = np.loadtxt('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'irsdata.txt', unpack=True)
        if (visibledatacheck):  w, f = np.loadtxt('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'visibledata.txt', unpack=True)
        if (jhkcheck): wdata, fdata = np.loadtxt('/Users/' + str(username) + '/rsgs/plots/text/' + str(starname) + 'data.txt', unpack=True)

    if (visibledatacheck): plt.plot(w, f) #plot visible
    if (irsdatacheck): plt.plot(nu1, fnu1) #irs
    if (jhkcheck): plt.plot(wdata, fdata) #jhk, iras mips data

    index = []
    numspecies = 0
    for i in range(len(percent_gs)):
        for j in range(len(percent_gs[i])):
            if percent_gs[i][j] > 0:
                index.append(percent_gs[i])
                numspecies += 1

#i fixed this part ish
"""
    if (residual == "residual"):
        if irsdatacheck:
            #FIND TAU
            wave, tau = np.loadtxt('/Users/mocassin/mocassin-rw_changes/output/tauNu.out', unpack=True)
            idx = 0
            for i in range(len(wave)):
                if (wave[i] > 0.52 and wave[i] <= 0.54):
                    idx = i
                    break
            wav = wave[idx]
            tau_print = tau[idx]
            extra_moc_plot(fnu1, nu1, data[2] / (D**2), data[1], irs_e_flux, ch1)
            star_dust_mass(mass, ch1, starname, id, username)
            if n > 0:
                extra_moc_plot(fnu1, nu1, data[3] / (D**2), data[1], irs_e_flux, ch2)
                star_dust_mass(mass, ch2, starname, id, username)
                if n > 1:
                    extra_moc_plot(fnu1, nu1, data[4] / (D**2), data[1], irs_e_flux, ch3)
                    star_dust_mass(mass, ch3, starname, id, username)
"""

    pdtk.plot_dust_temp_kelsey(username, distributionFile, symmetric)
    sk.cd_kelsey('/Users/' + username + '/mocassin-rw_changes/output')
    os.system('cp sn_smooth.eps ' + str(directoryname) + '/' + str(starname) + '_' + str(id) + '.eps')