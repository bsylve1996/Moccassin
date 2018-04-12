"""
this program plots the dust tempratures obtained
for a simulation with 20 grain sizes and 2 grain
species.
this was a multichemistry model and the results
are also plotted by sector.
needs the following input files: dustGrid.out,
grainsize file, grainspecies file, grid0.out.
"""

def ssf(number):
    return(str(number).replace(" ", ""))

def plot_dust_temp_amiy(username, diffuse, errorcheck = "errorcheck", includePAHS = "includePAHS"):

    cd_amiy, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/'
    if keyword_set(errorcheck) then print, 'plotting dust'


    nx = 0
    ny = 0
    nz = 0
    gridfile = 'output/grid0.out'
    openr, lun, gridfile, / get_lun
    ncells = 0
    motherP = 0
    rout = 0.
    readf, lun, motherP
    readf, lun, nx, ny, nz, ncells, motherP, rout
    convergence = replicate({n: intarr(3)}, nx, ny, nz)
    x = fltarr(nx)
    y = fltarr(ny)
    z = fltarr(nz)
    readf, lun, x
    readf, lun, y
    readf, lun, z
    readf, lun, convergence
    close, lun
    free_lun, lun

    if diffuse then sizefile='input/mrnsmall.dat' $
    else IF KEYWORD_SET(includePAHS) then sizefile='input/PAH_sizes.dat' $
    else sizefile = 'input/mrnsmall.dat'

    nsizes = 0
    openr, lun, sizefile, / get_lun
    readf, lun, nsizes
    dustsizes = replicate({index: intarr(1), radius: fltarr(1), weight: fltarr(1)}, nsizes)
    readf, lun, dustsizes
    close, lun
    free_lun, lun

    speciesfile = 'input/grainspecies.dat' #'primary_grainspecies.dat'
    nspecies = 0
    openr, lun, speciesfile, / get_lun
    readf, lun, nspecies
    close, lun
    free_lun, lun

    readcol, 'input/grainspecies.dat', dustnames, dustabundances, format = 'A,F', skipline = 1, / silent


    if keyword_set(errorcheck) then print, 'dnames, dabundances', dustnames, dustabundances

    mocassinfile = 'output/dustGrid.out'
    dustprop = replicate({n: fltarr(1), t: fltarr(nspecies + 1, nsizes + 1)}, nx, ny, nz)
    openr, lun, mocassinfile, / get_lun
    readf, lun, dustprop
    close, lun
    free_lun, lun

    r = fltarr(nx, ny, nz)
    celltot = nx * ny * nz
    r1 = fltarr(celltot)
    r2 = fltarr(celltot)
    r3 = fltarr(celltot)

    temp = replicate({t: fltarr(nsizes)}, nx, ny, nz)
    t1 = replicate({t: fltarr(nsizes)}, celltot)
    t2 = replicate({t: fltarr(nsizes)}, celltot)
    t3 = replicate({t: fltarr(nsizes)}, celltot)

    ii = 1
    jj = 1
    kk = 1
    for i = 0, nx-1 do
        for j = 0, nx-1 do
            for k =0, nx-1 do begin
                r[i, j, k] = sqrt(x[i] ^ 2 + x[j] ^ 2 + x[k] ^ 2)
                if (sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] gt 0.707107 and convergence[i, j, k].n[1] eq 1) then r1[ii] = r[i, j, k]
                if (sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] lt 0.707107 and convergence[i, j, k].n[1] eq 1) then r2[jj] = r[i, j, k]
                if (convergence[i, j, k].n[1] eq 1) then r3[kk] = r[i, j, k]

                for l = 0, nsizes-1 do begin
                    for m = 0, nspecies-1 do begin
                        temp[i, j, k].t[l] = temp[i, j, k].t[l] + dustprop[i, j, k].t[m + 1, l + 1] * dustabundances[m]
                        if (sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] gt 0.707107 and convergence[i, j, k].n[1] eq 1) then t1[ii] = temp[i, j, k]

                if (sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] gt 0.707107 and convergence[i, j, k].n[1] eq 1) then ii=ii+1
                if (sqrt(r[i, j, k] ^ 2 - x[i] ^ 2) / r[i, j, k] lt 0.707107 and convergence[i, j, k].n[1] eq 1) then jj=jj+1
                if (convergence[i, j, k].n[1] eq 1) then kk=kk+1


#This should start at the origin and move outward radially through
#the dust shell along the x - axis(to account for torus models!)
    diffuse = 1
    if diffuse eq 1 then begin
        plotRout = r[15:nx - 1, 15, 15]
        plotTemp_1 = temp[15:nx - 1, 15, 15].t[19]
        plotTemp_2 = temp[15:nx - 1, 15, 15].t[0]
        plot, plotRout, plotTemp_1, psym = -1, xrange = [0., max(PlotRout) + (0.25 * max(PlotRout))], xtitle = 'r [cm]',$;, yrange =
            [0., max(temp[0, 0, *].t[19])], yrange = [0, max(plotTemp_1) + (0.25 * max(plotTemp_1))], / ysty, ytitle = 'Temp (k)', / xsty, title = 'Temperature Vs. Radius'

        oplot, plotRout, plotTemp_2, psym = -1, col =!col.green #big granis
#stop
        plot, plotRout, plotTemp_1, title = 'Temperature Vs. Radius2'

        index = where(r[*, 15, 15] eq median(r[*, 15, 15]))
        ymax = max(temp[*, 15, 15].t[*])
        ymin = min(temp[1: *, 15, 15].t[*])
        plot, dustsizes.radius, temp[index, 15, 15].t[*], psym = -4, / ys, / xs, title = 'Temp vs size', $
            xtitle = 'size in microns', ytitle = 'Temp (k)', yrange = [ymin - 15, ymax + 15]#diamond
        oplot, dustsizes.radius, temp[2, 15, 15].t[*], psym = -5#triange

        elements = n_elements(r[*, 15, 15])
        oplot, dustsizes.radius, temp[elements - 2, 15, 15].t[*], psym = -6#square

        xyouts, xstart, ystart + N * yint, 'Square:   ' + ssf(r[elements - 2, 15, 15]), / normal & N + +
        xyouts, xstart, ystart + N * yint, 'Diamond:  ' + ssf(r[where(r[*, 15, 15] eq median(r[*, 15, 15]))]   ), / normal & N + +
        xyouts, xstart, ystart + N * yint, 'Triangle: ' + ssf(r[2, 15, 15]), / normal & N + +
    else begin
        plot, r[0, 0, *], temp[0, 0, *].t[19], psym = -1, xrange = [0., rout], xtitle = 'r [cm]',$;, yrange =
                [0., max(temp[0, 0, *].t[19])] ytitle = 'Temp (k)', / xsty, title = 'Temperature Vs. Radius'

        oplot, r[0, 0, *], temp[0, 0, *].t[0], psym = -1, col =!col.green #big granis
        plot, r[0, 0, *], temp[0, 0, *].t[19], title = 'Temperature Vs. Radius2'

        index = where(r[*, 0, 0] eq median(r[*, 0, 0]))
        ymax = max(temp[*, 0, 0].t[*])
        ymin = min(temp[1: *, 0, 0].t[*])
        plot, dustsizes.radius, temp[index, 0, 0].t[*], psym = -4, / ys, / xs, title = 'Temp vs size', $
            xtitle = 'size in microns', ytitle = 'Temp (k)', yrange = [ymin - 15, ymax + 15]#diamond
        oplot, dustsizes.radius, temp[2, 0, 0].t[*], psym = -5#triange

        elements = n_elements(r[*, 0, 0])
        oplot, dustsizes.radius, temp[elements - 2, 0, 0].t[*], psym = -6#square


        xyouts, xstart, ystart + N * yint, 'Square:   ' + ssf(r[elements - 2, 0, 0]), / normal & N + +
        xyouts, xstart, ystart + N * yint, 'Diamond:  ' + ssf(r[where(r[*, 0, 0] eq median(r[*, 0, 0]))]   ), / normal & N + +
        xyouts, xstart, ystart + N * yint, 'Triangle: ' + ssf(r[2, 0, 0]), / normal & N + +

    IF KEYWORD_SET(includePAHS)
        plot, dustsizes.radius, temp[index, 0, 0].t[*], psym = -4, / ys, / xs, title = 'Temp vs size', $
            xtitle = 'size in microns', ytitle = 'Temp (k)', yrange = [ymin - 15, ymax + 15], xrange = [min(dustsizes.radius), .01]#diamond
        oplot, dustsizes.radius, temp[2, 0, 0].t[*], psym = -5#triange

        elements = n_elements(r[*, 0, 0])
        oplot, dustsizes.radius, temp[elements - 2, 0, 0].t[*], psym = -6#square

    if keyword_set(errorcheck) then print, 'Done plotting dust Temperatures!'