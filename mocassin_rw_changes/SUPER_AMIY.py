#version 7

import numpy as np
import os
import sys
from os import walk as w
import time
import savReaderWriter as srw
import math

def make_input_amiy(j, torus, diffuse, n, luminosity, tstellar, nphotons,lstar, numiterations, convpercent, rin, rout, nova,theta1, theta2, phi1, phi2, username, nduststr, includePAHS = "includePAHS"):
#make the input. in file

    filenumber = open('/Users/' + username + '/dust_modelling/mocassin.2.02.70/input/input.in')

#super input
#torus must subtract one from the grid size
    if (isinstance(torus,list)):
        if torus[j]:
            nvalue = n[j] - 1
        else:
            nvalue= n[j]
        if torus[j] == 0:
            filenumber.write('symmetricXYZ')
    else:
        if torus:
            nvalue = n[j] - 1
        else:
            nvalue = n[j]
        if torus == 0:
            filenumber.write('symmetricXYZ')

#FOR ALL
    if (isinstance(diffuse, list)):
        if diffuse[j]:
#inputcfor SN
            filenumber.write('autoPackets 0.20 2. 900000000')
            filenumber.write('writeGrid 10.')
            filenumber.write('convLimit 0.01')
            filenumber.write('nPhotons 0')
            filenumber.write('LStar 0')
            filenumber.write('TStellar 0')
            filenumber.write("dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'")
            filenumber.write('diffuseSource ' + sss(luminosity[j]) + ' '+sss(tstellar) + " 'blackbody' " + sss(nphotons[j]) + " 1")
    #1 is for smooth
            filenumber.write('getEquivalentTau')

        else:
#input for the RSGs
            filenumber.write('autoPackets 0.10 3. 1000000000') #.20 5. 100000000
    #when convergence increase is less than first number percent increase, increase the
    #number of photons by 2nd number, to a maximum of the third.
            filenumber.write('writeGrid 80.')
            filenumber.write('convLimit 0.05')
            filenumber.write('nPhotons ' + sss(nphotons[j]))
            filenumber.write('LStar ' + sss(lstar))
            filenumber.write('TStellar ' + sss(tstellar))
    #IF THE PAHS are turned on...
            if includePAHS == "inclusePAHS":
                filenumber.write("dustFile 'input/grainspecies.dat' 'input/PAH_sizes.dat'")#input / PAH_sizes.dat;;input / mrn.dat
                filenumber.write('quantumHeatGrain .005 90')#size, convergence;.04 is 400 angstroms...
                print('Qheating turned ON!!!')
            else:
                filenumber.write("dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'")
    else:
        if diffuse:
            # inputcfor SN
            filenumber.write('autoPackets 0.20 2. 900000000')
            filenumber.write('writeGrid 10.')
            filenumber.write('convLimit 0.01')
            filenumber.write('nPhotons 0')
            filenumber.write('LStar 0')
            filenumber.write('TStellar 0')
            filenumber.write("dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'")
            filenumber.write(
                'diffuseSource ' + sss(luminosity) + ' ' + sss(tstellar) + " 'blackbody' " + sss(nphotons) + " 1")
            # 1 is for smooth
            filenumber.write('getEquivalentTau')

        else:
            # input for the RSGs
            filenumber.write('autoPackets 0.10 3. 1000000000')  # .20 5. 100000000
            # when convergence increase is less than first number percent increase, increase the
            # number of photons by 2nd number, to a maximum of the third.
            filenumber.write('writeGrid 80.')
            filenumber.write('convLimit 0.05')
            filenumber.write('nPhotons ' + sss(nphotons))
            filenumber.write('LStar ' + sss(lstar))
            filenumber.write('TStellar ' + sss(tstellar))
            # IF THE PAHS are turned on...
            if includePAHS == "inclusePAHS":
                filenumber.write(
                    "dustFile 'input/grainspecies.dat' 'input/PAH_sizes.dat'")  # input / PAH_sizes.dat;;input / mrn.dat
                filenumber.write('quantumHeatGrain .005 90')  # size, convergence;.04 is 400 angstroms...
                print('Qheating turned ON!!!')
            else:
                filenumber.write("dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'")

#all
    filenumber.write('output')


    filenumber.write('contShape  blackbody')
    filenumber.write('nebComposition noGas')
    filenumber.write('maxIterateMC ' + ssi(numiterations[j]) + ' ' + ssf(convpercent[j]))
    filenumber.write("Ndust file '" + sss(nduststr) + "'")
    filenumber.write('Rin ' + sss(rin[j]))
    filenumber.write('Rout ' + sss(rout[j]))
    filenumber.write('nx ' + sss(nvalue))#1 was already subtracted if it is a torus
    filenumber.write('ny ' + sss(nvalue))
    filenumber.write('nz ' + sss(nvalue))

#edit viewing angles
    if (isinstance(nova, list)):
        if nova[j] == 1:
            filenumber.write('inclination '+ssi(nova[j])+' '+ssf(theta1[j])+' '+ssf(phi1[j]))
        elif nova[j] == 2:
            filenumber.write('inclination '+ssi(nova[j])+' '+ssf(theta1[j])+' '+ssf(phi1[j])+' ' + ssf(theta2[j]) + ' ' + ssf(phi2[j]))
    else:
        if nova == 1:
            filenumber.write('inclination ' + ssi(nova[j]) + ' ' + ssf(theta1[j]) + ' ' + ssf(phi1[j]))
        elif nova == 2:
            filenumber.write('inclination ' + ssi(nova[j]) + ' ' + ssf(theta1[j]) + ' ' + ssf(phi1[j]) + ' ' + ssf(theta2[j]) + ' ' + ssf(phi2[j]))


def make_grainspecies_amiy(name, file, percent, silicatepercent, AMCpercent):
    filenumber = open('input/grainspecies.dat', "w")

    numspecies = 0
    index = []
    for i in percent:
        if (percent[i] > 0):
            numspecies += 1
            index.append(i)

    if (numspecies == 0):
        if (AMCpercent == 0):
            filenumber.write(ssi(1))
            filenumber.write("'dustData/sil-oss1.nk' " + sss(silicatepercent / 100.0))
        else:
            filenumber.write(ssi(2))
            filenumber.write("'dustData/sil-oss1.nk' " + sss(silicatepercent / 100.0))
            filenumber.write("'dustData/amC-hann.nk' " + sss(AMCpercent / 100.0))
            if (silicatepercent / 100.0 + AMCpercent / 100.0 != 1):
                print('error in grainspecies!!!')
    else:
        filenumber.write(ssi(numspecies))
        # THIS NUMBER MUST MATCH THE NUMBER OF dustData FILES!!!!
        for i in range(0, numspecies - 1):
            filenumber.write("'dustData/" + sss(file[index[i]]) + "' " + sss(percent[index[i]] / 100.0))

    filenumber.write('')
    # blank space at the end!
    filenumber.close()


def makeGrainSizeDistribution(*args):
#makeGrainSizeDistribution, 20, 3.5, 5e-3, .25, 'mrn_2.dat'
    if (len(args) != 6):
        nsizes = input('how many bins?')
        slope = input('what is the exponential giving the slope of the distribution? (p)')
        amin = input('what is amin in microns?')
        amax = input('what is amax in microns?')
        outfile = input('what should the output file be called?')
        username = input('what is your username??')
    else:
        arguments = []
        for i in args:
            arguments.add(i)
        nsizes = arguments[0]
        slope = arguments[1]
        amin = arguments[2]
        amax = arguments[3]
        outfile = arguments[4]
        username = arguments[5]

    astep = (np.log10(amax) - np.log10(amin)) / float(nsizes - 1)

    lun = open('/Users/' + username + '/mocassin-rw_changes/input/' + outfile)

    lun.write(nsizes + ' sizes ')
    for i in range(1, nsizes):
        lun.write(i + "   " + 10.0 ^ (np.alog10(amin) + (i - 1) * astep) + "   " (10.0 ^ (np.alog10(amin) + (i - 1) * astep)) ^ (-slope))

    lun.write('')
    lun.write('')
    lun.write('slope (p): ' + slope)
    lun.close()

    print('Done making new grain size distribution!!!!!')

def mocassin_fail_amiy(j, username, diffuse, directoryname, outfoldername, starname):
    print("RUN FAILED! Writing output.")
    print("Failed on line number" + (j + 1) + "of AMIY_input.txt")

    with srw.SavReader('/Users/' + username + '/mocassin-rw_changes/AMIY_number.sav') as reader:
        AMIY_number = reader.next()
    id = ssi(AMIY_number)
    AMIY_number += 1
    srw.SavWriter('/Users/'+username+'/mocassin-rw_changes/AMIY_number.sav',AMIY_number)

    if(diffuse[j]): type='SN'
    else: type='RSG'

    directoryname = "/Users/" + username + "/mocassin-rw_changes/output/" + type + "/" + id + '_' + starname + '_FAIlED'
    os.system("mkdir " + directoryname)
    outfoldername = type + "/" + id + '_' + starname + '_FAIlED'

    os.chdir('/Users/' + username + '/mocassin-rw_changes/output')

    os.system('cp dustGrid.out ' + directoryname + '/dustGrid_' + id + '.out.txt')
    os.system('cp runinfo.txt ' + directoryname + '/runinfo_' + id + '.txt')
    os.system('cp SED.out ' + directoryname + '/SED_' + id + '.out.txt')
    if (diffuse[j]):
        os.system('cp equivalentTau.out ' + directoryname + '/equivalentTau_' + id + '.out.txt')
    else:
        os.system('cp tauNu.out ' + directoryname + '/tauNu_' + id + '.out.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/input.in ' + directoryname + '/input_' + id + '.in.txt')
    os.system('cp /Users/' + username + '/mocassin-rw_changes/input/ndust/nDUST ' + directoryname + '/nDUST_' + id + '.in.txt')



def check_temp_AMIY(temperature, temp, starnum, diffuse, errorcheck = "errorcheck"):
#first check for user input temp
#next search using starname...
#finally go to a default temp...
#default for RSG! AND SN are different!
    tstellar = 0
    if (temperature > 0):
        tstellar = temperature
        if (errorcheck == "errorcheck"):
            print("user supplied a temperature, so we are using temperature " + str(tstellar))
    elif (tstellar <= 0):
        if (starnum > -1):  #if the star is in the list...
            tstellar = temp[starnum]
            if (errorcheck == "errorcheck"):
                print('found temperature from file, that temperature is ' + str(tstellar))
        else: #if not in list use a default
            if diffuse:
                default=6000.0
            else:
                default=3500.0
        tstellar = default
        if (errorcheck == "errorcheck"):
            print("deafulting to temp of " + str(tstellar))

    if (starnum > -1):
        if (temperature > 0 and temperature != temp[starnum]):
            if (errorcheck == "errorcheck"):
                print("Temp in kelsey_input is different than RSG_info.txt (in Bill's folder)")
                print("TStellar= " + str(tstellar))
                print("temp from RSG_info.txt is " + temp[starnum])
    return(tstellar)


def check_lstar_AMIY(luminosity, luminositycalc, starnum, diffuse, errorcheck = "errorcheck"):
    lstar = 0
    if (luminosity > 0):
        lstar = luminosity
        if (errorcheck == "errorcheck"):
            print("user supplied a luminosity, so we are using luminosity " + str(lstar))
    else:
        if (starnum > -1):
            lstar = luminositycalc[starnum]
            if (errorcheck == "errorcheck"):
                print('found luminosity from file, that luminosity is ' + str(lstar))
        else:      #SN DEFAULT........RSG DEFAULT
            default=650.0
            if (errorcheck == "errorcheck"):
                print('Luminosity is defaulting...')
            lstar = default

    return(lstar)

def makenuryd(diffuse, old, errorcheck = "errorcheck"):
    c = 2.9979250e10
    ryd = 3.28984e15
    if (os.path.isfile('dustData/wavelength_resolution.txt')):
        lambdaAstro = np.loadtxt('dustData/wavelength_resolution.txt', unpack=True)
        nu = []
        for i in lambdaAstro:
            nu.append(c/lambdaAstro[i] * 1.e-4 * ryd)

        ct = 0
        index = []
        for i in lambdaAstro:
            if (lambdaAstro[i] == .5623):
                ct += 1
                index.append(i)

        if (ct == 0):
            print('FAIL!  You need a lambda of .5623 to get tau out silly pants!')

        if diffuse:
            np.savetxt('dustData/nuRyd.dat', nu)
        else:
            np.savetxt('dustData/nuDustRyd.dat', nu)

    if (old == 1): #restores the old NU file...
        if (os.path.isfile('dustData/nuDustRyd_default.txt')):
            nu = np.loadtxt('dustData/nuDustRyd_default.txt', unpack=True)
            #nu = c / (lambda * 1.e-4 * ryd)  ## ALREADY IN rydbergs.
            print('restoring OLD')
        if diffuse:
            np.savetxt('dustData/nuRyd.dat', nu)
        else:
            np.savetxt('dustData/nuDustRyd.dat', nu)


def cd_AMIY(new_dir, errorcheck = "errorcheck"):
    old_dir = '.'
    os.chdir(old_dir)
    if (errorcheck == "errorcheck"):
        print('old dir =', old_dir)
        print('new dir =', new_dir)
    if (old_dir != new_dir):
        os.chdir(new_dir)

from math import floor

def ssf(number):
    return(str(number).replace(" ", ""))
def ssi(number):
    return(str(int(floor(number))).replace(" ", ""))
def sss(number):
    return(str(number).replace(" ", ""))

def SUPER_AMIY(infile, dontcheckinput = "dontcheckinput", errorcheck = "errorcheck"):

    includePAHS = 0 #1 includes pahs... 0 is for normal run...
    errorcheck = 0 #1 prints status statements... 0 does not!
    dontcheckinput = 1 #1 prevents error checking... 0 checks for errors
    old = 1 #1 = old nuDustRyd.dat, 0 uses wavelength from dustdata / wavelength_resolution.txt

    ans = 'str'
    username = 'mocassin'

    cd_AMIY('/Users/' + username + '/mocassin-rw_changes')

    # read in the input from the user
    data = np.loadtxt('AMIY_input.txt', unpack=True, skiprows=3, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22))
    starname = np.loadtxt('KELSEY_input.txt', unpack=True, skiprows=3, usecols=(0), dtype='str')
    n = data[0]
    rin = data[1]
    rout = data[2]
    p = data[3]
    constant = data[4]
    nphotons = data[5]
    numiterations = data[6]
    convpercent = data[7]
    diffuse = data[8]
    temperature = data[9]
    luminosity = data[10]
    silicatepercent = data[11]
    AMCpercent = data[12]
    distance = data[13]
    torus = data[14]
    dens1 = data[15]
    dens2 = data[16]
    nova = data[17]
    theta1 = data[18]
    phi1 = data[19]
    theta2 = data[20]
    phi2 = data[21]

    if isinstance(starname, list):
        numrows = len(starname)
    else:
        numrows = 1

    #NAME, FILE, PERCENT
    if (os.path.isfile('AMIY_input_grainspecies.txt')):
        data = np.loadtxt('AMIY_input_grainspecies.txt', unpack=True, dtype='str', usecols=(0,1,3,4,6,7,9,10))
        data2 = np.loadtxt('AMIY_input_grainspecies.txt', unpack=True, usecols=(2,5,8,11))
        gs1 = data[0]
        file1 = data[1]
        gs2 = data[2]
        file2 = data[3]
        gs3 = data[4]
        file3 = data[5]
        gs4 = data[6]
        file4 = data[7]
        percent1 = data2[0]
        percent2 = data2[1]
        percent3 = data2[2]
        percent4 = data2[3]
    else:
        percent1 = np.repeat(0, numrows)
        percent2 = np.repeat(0, numrows)
        percent3 = np.repeat(0, numrows)
        percent4 = np.repeat(0, numrows)
        gs1 = np.repeat('nothing', numrows)
        gs2 = np.repeat('nothing', numrows)
        gs3 = np.repeat('nothing', numrows)
        gs4 = np.repeat('nothing', numrows)
        file1 = np.repeat('nothing', numrows)
        file2 = np.repeat('nothing', numrows)
        file3 = np.repeat('nothing', numrows)
        file4 = np.repeat('nothing', numrows)

        # correct so that if the user inputs a file
        # that does not have enough rows, kelsey wont crash
    if isinstance(percent1, list):
        if len(percent1) < numrows:
            percent1 = np.repeat(0, numrows)
            percent2 = np.repeat(0, numrows)
            percent3 = np.repeat(0, numrows)
            percent4 = np.repeat(0, numrows)
            gs1 = np.repeat('nothing', numrows)
            gs2 = np.repeat('nothing', numrows)
            gs3 = np.repeat('nothing', numrows)
            gs4 = np.repeat('nothing', numrows)
            file1 = np.repeat('nothing', numrows)
            file2 = np.repeat('nothing', numrows)
            file3 = np.repeat('nothing', numrows)
            file4 = np.repeat('nothing', numrows)
            if includePAHS == 1:
                import pymsgbox as p
                p.alert(text='ERROR', title='REMEMBER TO CHANGE AMIY_INPUT_GRAINSPECIES.TXT', button='OK')
                endOfProgram(1)
                sys.exit()
    else:
        if len(percent1) < numrows:
            percent1 = np.repeat(0, numrows)
            percent2 = np.repeat(0, numrows)
            percent3 = np.repeat(0, numrows)
            percent4 = np.repeat(0, numrows)
            gs1 = np.repeat('nothing', numrows)
            gs2 = np.repeat('nothing', numrows)
            gs3 = np.repeat('nothing', numrows)
            gs4 = np.repeat('nothing', numrows)
            file1 = np.repeat('nothing', numrows)
            file2 = np.repeat('nothing', numrows)
            file3 = np.repeat('nothing', numrows)
            file4 = np.repeat('nothing', numrows)
            if includePAHS == 1:
                import pymsgbox as p
                p.alert(text='ERROR', title='REMEMBER TO CHANGE AMIY_INPUT_GRAINSPECIES.TXT', button='OK')
                endOfProgram(1)
                sys.exit()


    #temperature and luminosity are optional, set to 0 if unknown
    #for RSG, luminosity will be derived from the temp and the radius
    #both temp and radius will be stolen from lmcsmcfullspec2e.txt
    #as long as the star name is correct, it will find the temperature, radius, and luminosity.

    #if AMIY can't find the star's name as an EXACT match to the name in lmcsmcfullspec2e.txt,
    #then it will try to use the user input temperature and luminosity from AMIY_input.txt
    #if it STILL can't find these variables, it will default to some value....

    #temp VS temperature VS Tstellar
    #temperature is what the USER has supplied us with in AMIY_input.txt {THIS IS AN ARRAY!!}
    #temp is what we get from LMCSMCfullspec2e.txt
    #tstellar is what we finally give MOCASSIN as an input(for both RSG and SN)

    #luminosity Vs.lstar
    #luminosity is from the user {THIS IS AN ARRAY}
    #lstar is what we give to MOCASSIN as an input
    #if luminosity is not given, then AMIY will try to calculate it from radius and temperature.

    theta1 = theta1 * (math.pi / 180.)
    theta2 = theta2 * (math.pi / 180.)

    rsg_file = '/Users/mocassin/mocassin-rw_changes/accessories/RSG_info.txt'
    name = np.loadtxt(rsg_file, unpack=True, skiprows=1, dtype='str', usecols=(0))
    temp = np.loadtxt(rsg_file, unpack=True, skiprows=1, usecols=(1))
    luminositycalc = np.loadtxt(rsg_file, unpack=True, skiprows=1, usecols=(5))

    j = 0
    name_gs = []
    percent_gs = []
    filename_gs = []
    while (j < numrows):
        if (isinstance(starname, list)):
            print("Beginning run number" + str(j + 1) + "  for  " + starname[j])
        else:
            print("Beginning run number" + str(j + 1) + "  for  " + str(starname))

        successtest = 0  # check for if mocassian was successful
        lstar = 0  # luminosity used by mocassin
        tstellar = 0  # tstellar "     "   "
        starnum = -1  # number of the star in rsg_info

        # grain species from the file
        name_gs.append([])
        name_gs[j].append(gs1[j])
        name_gs[j].append(gs2[j])
        name_gs[j].append(gs3[j])
        name_gs[j].append(gs4[j])
        percent_gs.append([])
        percent_gs[j].append(percent1[j])
        percent_gs[j].append(percent2[j])
        percent_gs[j].append(percent3[j])
        percent_gs[j].append(percent4[j])
        filename_gs.append([])
        filename_gs[j].append(file1[j])
        filename_gs[j].append(file2[j])
        filename_gs[j].append(file3[j])
        filename_gs[j].append(file4[j])

        # searchfor the starnum
        i = 0
        while (i < len(name) and starnum == -1):
            if (isinstance(starname, list)):
                if (name[i] == starname[j]):
                    starnum = i
            else:
                if (name[i] == starname):
                    starnum = i
            i += 1
        if (isinstance(temperature, list)):
            if (isinstance(diffuse, list)):
                tstellar = check_temp_AMIY(temperature[j], temp, starnum, diffuse[j])
                if (isinstance(luminosity, list)):
                    lstar = check_lstar_AMIY(luminosity[j], luminositycalc, starnum, diffuse[j])
                else:
                    lstar = check_lstar_AMIY(luminosity, luminositycalc, starnum, diffuse[j])
            else:
                tstellar = check_temp_AMIY(temperature[j], temp, starnum, diffuse)
                if (isinstance(luminosity, list)):
                    lstar = check_lstar_AMIY(luminosity[j], luminositycalc, starnum, diffuse)
                else:
                    lstar = check_lstar_AMIY(luminosity, luminositycalc, starnum, diffuse)
        else:
            if (isinstance(diffuse, list)):
                tstellar = check_temp_AMIY(temperature, temp, starnum, diffuse[j])
                if (isinstance(luminosity, list)):
                    lstar = check_lstar_AMIY(luminosity[j], luminositycalc, starnum, diffuse[j])
                else:
                    lstar = check_lstar_AMIY(luminosity, luminositycalc, starnum, diffuse[j])
            else:
                tstellar = check_temp_AMIY(temperature, temp, starnum, diffuse)
                if (isinstance(luminosity, list)):
                    lstar = check_lstar_AMIY(luminosity[j], luminositycalc, starnum, diffuse)
                else:
                    lstar = check_lstar_AMIY(luminosity, luminositycalc, starnum, diffuse)
        if (errorcheck == "errorcheck"):
            print('ts=', tstellar)
            print('ls=', lstar)

    # change dust percent composition
    # make dust shell, returns the name of the ndust file...
        nduststr = np.loadtxt('nDust.list', dtype='str', unpack=True)

        cd_AMIY('/Users/' + username + '/dust_modelling/mocassin.2.02.70')

        if (isinstance(n,list)):
            n[j] = floor(n[j])
        else:
            n = floor(n)

    #MAKE THE GRAINSPECIES FILES!!!
        if (errorcheck == "errorcheck"):
            print('calling grainspecies function')
        if (isinstance(silicatepercent, list)):
            if(isinstance(AMCpercent, list)):
                make_grainspecies_amiy(name_gs, filename_gs, percent_gs, silicatepercent[j], AMCpercent[j])
            else:
                make_grainspecies_amiy(name_gs, filename_gs, percent_gs, silicatepercent[j], AMCpercent)
        else:
            if (isinstance(AMCpercent, list)):
                make_grainspecies_amiy(name_gs, filename_gs, percent_gs, silicatepercent, AMCpercent[j])
            else:
                make_grainspecies_amiy(name_gs, filename_gs, percent_gs, silicatepercent, AMCpercent)

    #MAKE THE INPUT.IN FILES!!!
        if (errorcheck == "errorcheck"):
            print('calling input.in making function')
        if(isinstance(nduststr,list)):
            make_input_amiy(j, torus, diffuse, n, luminosity, tstellar, nphotons,lstar, numiterations, convpercent, rin, rout, nova,theta1, theta2, phi1, phi2, username, nduststr[j])
        else:
            make_input_amiy(j, torus, diffuse, n, luminosity, tstellar, nphotons, lstar, numiterations, convpercent, rin, rout, nova, theta1, theta2, phi1, phi2, username, nduststr)

        # MAKE THE WAVELENGTH RESOLUTION FILE...
        if (isinstance(diffuse, list)):
            makenuryd(diffuse[j], old)
        else:
            makenuryd(diffuse, old)

        if (includePAHS == 1):
            makeGrainSizeDistribution(20, 3.5, 3.548e-04, .25, 'PAH_sizes.dat', username)

        start_time = time.time()
        runcounter = 0
        while (successtest != 1 and runcounter <= 1):
            if (runcounter == 1):
                print("running MOCASSIN again!!!")
            print("beginning to run MOCASSIN!!!")

            os.system('mpiexec -n 4 ./mocassin >output/runinfo.txt') #CMD TO RUN MOCASSIN


    #test if the run was successful
            os.system("grep -c '! mocassin: MCIterationDriver done' /Users/" + username + "/dust_modelling/mocassin.2.02.70/output/runinfo.txt > test1.txt")
            successtest = np.loadtxt('test1.txt', unpack=True)
            if (successtest == 1):
                print('Mocassin was a great success!!!')
            else:
                print('Mocassin = epic fail!')
            os.system("rm test1.txt")
            if (isinstance(starname,list)):
                if (successtest == 0):
                    mocassin_fail_amiy(j, username, diffuse, os.getcwd(), outfoldername, starname[j])
            else:
                if (successtest == 0):
                    mocassin_fail_amiy(j, username, diffuse, os.getcwd(), outfoldername, starname)
            runcounter += 1

        if (successtest == 1):
            totaltime = time.time() - start_time
        print('Total run time is ' + str(ssi(totaltime)) + ' seconds (' + str(ssi((totaltime) / 60.0)) + ')')
        outfoldername = 'str'

        cd_AMIY('/Users/' + username + '/dust_modelling/mocassin.2.02.70/')

    #make an AMIY output with what AMIY did.
        filenumber = open('/Users/' + username + '/dust_modelling/mocassin.2.02.70/output/' + outfoldername + '/AMIY_output.txt','w')
        filenumber.write("all AMIY_input.txt variables...")

        if (isinstance(starname,list)):
            filenumber.write(starname[j], n[j], rin[j], rout[j], p[j], constant[j], nphotons[j], numiterations[j],convpercent[j], diffuse[j], temperature[j], luminosity[j], silicatepercent[j], AMCpercent[j], torus[j],dens1[j], dens2[j], nova[j], theta1[j], phi1[j], theta2[j], phi2[j])

            filenumber.write("inputs to mocplot")
            filenumber.write(rin[j], rout[j], constant[j], p[j], lstar, tstellar, ' ' + starname[j], diffuse[j],' ' + username, ' ' + outfoldername, starnum, nova[j], theta1[j], phi1[j], theta2[j], phi2[j], torus[j],silicatepercent[j], AMCpercent[j], ' ' + name_gs[0], ' ' + name_gs[1], ' ' + name_gs[2], ' ' + name_gs[3], ' ' + filename_gs[0],' ' + filename_gs[1], ' ' + filename_gs[2], ' ' + filename_gs[3], percent_gs[0], percent_gs[1], percent_gs[2], percent_gs[3])

            filenumber.write("starname        ", starname[j])
            filenumber.write("n               ", n[j])
            filenumber.write("rin             ", rin[j])
            filenumber.write("rout            ", rout[j])
            filenumber.write("p               ", p[j])
            filenumber.write("constant        ", constant[j])
            filenumber.write("nphotons        ", nphotons[j])
            filenumber.write("numiterations   ", numiterations[j])
            filenumber.write("convpercent     ", convpercent[j])
            filenumber.write("diffuse         ", diffuse[j])
            filenumber.write("temperature     ", temperature[j])
            filenumber.write("luminosity      ", luminosity[j])
            filenumber.write("silicatepercent ", silicatepercent[j])
            filenumber.write("AMCpercent      ", AMCpercent[j])
            filenumber.write("inclination     ", nova[j], theta1[j], phi1[j], theta2[j], phi2[j])
            filenumber.write("                ")
            filenumber.write("starnum         ", starnum)
            filenumber.write("tstellar        ", tstellar)
            filenumber.write("Lstar           ", lstar)
            filenumber.write("  ")
            filenumber.write('namegs', name_gs)
            filenumber.write('filenamegs', filename_gs)
            filenumber.write('percent', percent_gs)
        else:
            filenumber.write(starname, n, rin, rout, p, constant, nphotons, numiterations, convpercent, diffuse, temperature, luminosity, silicatepercent, AMCpercent, torus, dens1, dens2, nova, theta1, phi1, theta2, phi2)

            filenumber.write("inputs to mocplot")
            filenumber.write(rin, rout, constant, p, lstar, tstellar, ' ' + starname, diffuse, ' ' + username, ' ' + outfoldername, starnum, nova, theta1, phi1, theta2, phi2, torus, silicatepercent, AMCpercent, ' ' + name_gs[0], ' ' + name_gs[1],' ' + name_gs[2], ' ' + name_gs[3], ' ' + filename_gs[0], ' ' + filename_gs[1],' ' + filename_gs[2], ' ' + filename_gs[3], percent_gs[0], percent_gs[1], percent_gs[2],percent_gs[3])

            filenumber.write("starname        ", starname)
            filenumber.write("n               ", n)
            filenumber.write("rin             ", rin)
            filenumber.write("rout            ", rout)
            filenumber.write("p               ", p)
            filenumber.write("constant        ", constant)
            filenumber.write("nphotons        ", nphotons)
            filenumber.write("numiterations   ", numiterations)
            filenumber.write("convpercent     ", convpercent)
            filenumber.write("diffuse         ", diffuse)
            filenumber.write("temperature     ", temperature)
            filenumber.write("luminosity      ", luminosity)
            filenumber.write("silicatepercent ", silicatepercent)
            filenumber.write("AMCpercent      ", AMCpercent)
            filenumber.write("inclination     ", nova, theta1, phi1, theta2, phi2)
            filenumber.write("                ")
            filenumber.write("starnum         ", starnum)
            filenumber.write("tstellar        ", tstellar)
            filenumber.write("Lstar           ", lstar)
            filenumber.write("  ")
            filenumber.write('namegs', name_gs)
            filenumber.write('filenamegs', filename_gs)
            filenumber.write('percent', percent_gs)
        j += 1

def endOfProgram(fatalerror=0):
    print("PROGRAM KELSEY IS OVER!")
    print("PROGRAM KELSEY IS OVER!")
    print("PROGRAM KELSEY IS OVER!")
    # spawn, 'rm /Users/mocassin/dust_modelling/mocassin.2.02.70/output/runinfo.txt'
    if (fatalerror == 1):
        print("(because of a fatal error!!!)")

def main(argv):
    SUPER_AMIY(argv[0],argv[1])

if __name__ == "__main__":
    main(sys.argv[1:])