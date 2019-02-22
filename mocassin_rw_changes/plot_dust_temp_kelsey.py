#Komparing Effectively Like SEDs Efficiently....Yeah
#modified from plot_dust_temp_amiy.pro
#this program plots the dust tempratures obtained
# for a simulation with N grain sizes and 2 grain
# species.

import SUPER_KELSEY as sk
import numpy as np
from Scanner import Scanner
import math
import matplotlib.pyplot as plt

def ssf(number):
    return(str(number).replace(" ", ""))

def plot_dust_temp_kelsey(username, distributionFile, symmetric, errorcheck = "errorcheck", includePAHS = "includePAHS"):

    sk.cd_kelsey('/Users/' + username + '/mocassin-rw_changes/')
    #if (errorcheck == "errorcheck"):
    print('plotting dust')

    lun = Scanner(file = 'output/grid0.out')
    #lun = open('output/grid0.out')
    lun.next()
    nx = lun.next_int()
    ny = lun.next_int()
    nz = lun.next_int()
    lun.next_int()
    lun.next_int()
    rout = lun.next_float()
    class MyStructConver():
        def __init__(self, sizeofn):
            self.n = [0 for x in range(sizeofn)]
    convergence = [[[MyStructConver(3)for col in range(nz)]for row in range(ny)] for x in range(nx)]
    x = []
    for xl in range(nx):
        x.append(lun.next_float)
    y = []
    for yl in range(ny):
        y.append(lun.next_float)
    z = []
    for zl in range(nz):
        z.append(lun.next_float)
    for x1 in range(nx):
        for y1 in range(ny):
            for z1 in range(nz):
                for t1 in range(3):
                    convergence[x1][y1][z1].n[t1] = lun.next_float()

    sizefile = str(distributionFile)
    #if (includePAHS == "includePAHS"):
    #    sizefile = 'input/PAH_sizes.dat'
    lun = Scanner(file=sizefile)
    nsizes = lun.next_int()
    print(sizefile)
    #class MyStructDustSizes():
     #   def __init__(self, index, radius, weight):
      #      self.index = index
       #     self.radius = radius
        #    self.weight = weight

    #dustsizes = []
    #for i in range(nsizes):
     #   dustsizes.append(MyStructDustSizes())

    speciesfile = 'input/grainspecies.dat' #'primary_grainspecies.dat'
    lun = Scanner(file=speciesfile)
    nspecies =  lun.next_int()
    dustnames = []
    dustabundances = []
    while(lun.has_next()):
        dustnames.append(lun.next())
        dustabundances.append(lun.next())

    #if errorcheck.equals("errorcheck"):
    print('dnames, dabundances', dustnames, dustabundances)

    mocassinfile = 'output/dustGrid.out'

    class MyStructDustProp():
        def __init__(self, sizeofn, sizeoft, sizeoftrow):
            self.n = [0 for x in range(sizeofn)]
            self.t = [[0 for row in range(sizeoft)] for x in range(sizeoftrow)]

    dustprop = [[[MyStructDustProp(1,nspecies+1,nsizes+1) for col in range(nx)]for row in range(ny)] for x in range(nz)]
    lun = Scanner(file=mocassinfile)
    for x1 in range(nx):
        for y1 in range(ny):
            for z1 in range(nz):
                dustprop[x1][y1][z1].n[0] = lun.next_float()
                for t1 in range(nspecies+1):
                    for t2 in range(nsizes+1):
                        dustprop[x1][y1][z1].t[t1][t2] = lun.next_float()
    r = [[[0 for col in range(nz)]for row in range(ny)] for x in range(nx)]
    celltot = nx * ny * nz
    r1 = [0 for col in range(celltot)]
    r2 = [0 for col in range(celltot)]
    r3 = [0 for col in range(celltot)]

    class MyStructT():
        def __init__(self, sizeoft):
            self.t = [0 for x in range(sizeoft)]
    temp = [[[MyStructT(nsizes)for col in range(nz)]for row in range(ny)] for x in range(nx)]
    t1 = [MyStructT(nsizes) for x in range(nx)]
    t2 = [MyStructT(nsizes) for x in range(nx)]
    t3 = [MyStructT(nsizes) for x in range(nx)]

    ii = 1
    jj = 1
    kk = 1
    for i in range(nx-1):
        for j in range(nx-1):
            for k in range(nx-1):
                r[i, j, k] = math.sqrt(x[i] ^ 2 + x[j] ^ 2 + x[k] ^ 2)
                if (convergence[i, j, k].n[1] == 1):
                    if (math.sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] > 0.707107):
                        r1[ii] = r[i, j, k]
                    elif (math.sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] < 0.707107):
                        r2[jj] = r[i, j, k]
                    r3[kk] = r[i, j, k]

                for l in range(nsizes-1):
                    for m in range(nspecies-1):
                        temp[i, j, k].t[l] = temp[i, j, k].t[l] + dustprop[i, j, k].t[m + 1, l + 1] * dustabundances[m]
                        if (math.sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] > 0.707107 and convergence[i, j, k].n[1] == 1):
                            t1[ii] = temp[i, j, k]
        if (convergence[i, j, k].n[1] == 1):
            if (math.sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] > 0.707107 ):
                ii=ii+1
            if (math.sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] < 0.707107):
                jj=jj+1
            kk=kk+1

#This should start at the origin and move outward radially through
#the dust shell along the x - axis(to account for torus models!)
    plt.figure(1)
    if (symmetric == 0):
        plotRout = r[15:nx - 1, 15, 15]
        plotTemp_1 = temp[15:nx - 1, 15, 15].t[nsizes - 1] #covers the largest grain
        plotTemp_2 = temp[15:nx - 1, 15, 15].t[0] #covers the smallest grain

        plt.subplot(221)
        plt.plot(plotRout,plotTemp_1, marker = '+')
        plt.plot(plotRout, plotTemp_2, marker = '+')
        plt.axis([0.,max(plotRout) + (0.25 * max(plotRout)),0, max(plotTemp_1) + (0.25 * max(plotTemp_1))])
        plt.xlabel('r [cm]')
        plt.ylabel('Temp (k)')
        plt.title('Temperature Vs. Radius')

        #stop
        plt.subplot(222)
        plt.plot(plotRout, plotTemp_1)
        plt.title('Temperature Vs. Radius2')
        index = []
        for i in len(r):
            if (r[i,15,15] == math.median(r[i,15,15])):
                index.append(i)
        ymax = -1e308
        for i in len(temp):
            for j in len(temp[i,15,15].t):
                if(temp[i,15,15].t[j] > ymax):
                    ymax = temp[i,15,15].t[j]
        ymin = ymax
        for i in len(temp):
            for j in len(temp[i, 15, 15].t):
                if (temp[i, 15, 15].t[j] < ymin):
                    ymin = temp[i, 15, 15].t[j]
        elements = 0
        for i in range(nx - 1):
            if (r[i, 15, 15].exists()):
                elements += 1
        plt.subplot(223)
        plt.plot(dustsizes.radius, temp[index, 15, 15].t, marker = 'D', label = 'a')
        plt.plot(dustsizes.radius, temp[2, 15, 15].t, marker = '^', label = 'b')
        plt.plot(dustsizes.radius, temp[elements - 2, 15, 15].t, marker = 's', label = 'c')
        plt.title('Temp vs size')
        plt.xlabel('size in microns')
        plt.ylabel('Temp (k)')
        darr = np.array(dustsizes.radius)
        plt.axis([darr.min(),darr.max(),ymin - 15,ymax + 15])
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)

        #xyouts, xstart, ystart + N * yint, 'Square:   ' + ssf(r[elements - 2, 15, 15]), / normal & N + +
        #xyouts, xstart, ystart + N * yint, 'Diamond:  ' + ssf(r[where(r[*, 15, 15] eq median(r[*, 15, 15]))]   ), / normal & N + +
        #xyouts, xstart, ystart + N * yint, 'Triangle: ' + ssf(r[2, 15, 15]), / normal & N + +
    else:
        plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])
        plt.subplot(221)
        plt.plot(r[0, 0, 0:nx-1], temp[0, 0, 0:nx-1].t[nsizes - 1], '+', r[0, 0, 0:nx-1], temp[0, 0, 0,nx-1].t[0], '+', linewidth = 1)
        maxTemp = -1e308
        for i in range (nx - 1):
            if (temp[0,0,i].t19 > maxTemp):
                maxTemp = temp[0,0,i].t[19]
        plt.axis([0., rout, 0., maxTemp])
        plt.xlabel('r [cm]')
        plt.ylabel('Temp (k)')
        plt.title('Temperature Vs. Radius')

        plt.subplot(222)
        plt.plot(r[0, 0,0:nx-1], temp[0, 0,0:nx-1].t[nsizes - 1])
        plt.title = 'Temperature Vs. Radius2'

        index = []
        for i in len(r):
            if (r[i, 0, 0] == math.median(r[i, 0, 0])):
                index.append(i)
        ymax = -1e308
        for i in len(temp):
            for j in len(temp[i, 0, 0].t):
                if (temp[i, 0, 0].t[j] > ymax):
                    ymax = temp[i, 0, 0].t[j]
        ymin = ymax
        for i in len(temp):
            for j in len(temp[i, 0, 0].t):
                if (temp[i, 0, 0].t[j] < ymin):
                    ymin = temp[i, 0, 0].t[j]
        elements = 0
        for i in range(nx - 1):
            if (r[i, 0, 0].exists()):
                elements += 1
        plt.subplot(223)
        plt.plot(dustsizes.radius, temp[index, 0, 0].t,marker = 'D')
        plt.plot(dustsizes.radius, temp[2, 0, 0].t,marker = '^')
        plt.plot(dustsizes.radius, temp[elements - 2, 0, 0].t,marker = 's')
        plt.title('Temp vs size')
        plt.xlabel('size in microns')
        plt.ylabel('Temp (k)')
        darr = np.array(dustsizes.radius)
        plt.axis([darr.min(), darr.max(), ymin - 15, ymax + 15])
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)

        #xyouts, xstart, ystart + N * yint, 'Square:   ' + ssf(r[elements - 2, 0, 0]), / normal & N + +
        #xyouts, xstart, ystart + N * yint, 'Diamond:  ' + ssf(r[where(r[*, 0, 0] eq median(r[*, 0, 0]))]   ), / normal & N + +
        #xyouts, xstart, ystart + N * yint, 'Triangle: ' + ssf(r[2, 0, 0]), / normal & N + +

    if (includePAHS == "includePAHS"):
        plt.figure(2)
        elements = 0
        for i in range(nx - 1):
            if (r[i, 0, 0].exists()):
                elements += 1
        plt.plot(dustsizes.radius, temp[index, 0, 0].t, marker = 'D')
        plt.plot(dustsizes.radius, temp[2, 0, 0].t, marker = '^')
        plt.plot(dustsizes.radius, temp[elements - 2, 0, 0].t, marker = 's')
        plt.title('Temp vs size')
        plt.xlabel('size in microns')
        plt.ylabel('Temp (k)')
        plt.axis([min(dustsizes.radius), .01, ymin - 15, ymax + 15])
    plt.show()
    if (errorcheck == "errorcheck"):
        print('Done plotting dust Temperatures!')
