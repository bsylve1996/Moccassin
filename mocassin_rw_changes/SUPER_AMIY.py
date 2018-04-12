#version 7

def make_input_amiy(j, torus, diffuse, n, luminosity, tstellar, nphotons,lstar, numiterations, convpercent, rin, rout, nova,theta1, theta2, phi1, phi2, username, nduststr, includePAHS = "includePAHS"):
#make the input. in file

    get_LUN, filenumber
    openw, filenumber, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/input/input.in'

#super input
#torus must subtract one from the grid size
    if torus[j] then nvalue=n[j]-1 else nvalue=n[j]
    if torus[j] eq 0 then printf, filenumber, 'symmetricXYZ'

#FOR ALL
    if diffuse[j] then begin
#inputcfor SN
        printf, filenumber, 'autoPackets 0.20 2. 900000000'
        printf, filenumber, 'writeGrid 10.'
        printf, filenumber, 'convLimit 0.01'
        printf, filenumber, 'nPhotons 0'
        printf, filenumber, 'LStar 0'
        printf, filenumber, 'TStellar 0'
        printf, filenumber, "dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'"
        printf, filenumber, 'diffuseSource ' + sss(luminosity[j]) + ' '+sss(tstellar) + " 'blackbody' " + sss(nphotons[j]) + " 1"
    #1 is for smooth
        printf, filenumber, 'getEquivalentTau'

    else begin
#input for the RSGs
        printf, filenumber, 'autoPackets 0.10 3. 1000000000' #.20 5. 100000000
    #when convergence increase is less than first number percent increase, increase the
    #number of photons by 2nd number, to a maximum of the third.
        printf, filenumber, 'writeGrid 80.'
        printf, filenumber, 'convLimit 0.05'
        printf, filenumber, 'nPhotons ' + sss(nphotons[j])
        printf, filenumber, 'LStar ' + sss(LStar)
        printf, filenumber, 'TStellar ' + sss(tstellar)
    #IF THE PAHS are turned on...
        if keyword_set(includePAHS) then begin
            printf, filenumber, "dustFile 'input/grainspecies.dat' 'input/PAH_sizes.dat'"#input / PAH_sizes.dat;;input / mrn.dat
            printf, filenumber, 'quantumHeatGrain .005 90'#size, convergence;.04 is 400 angstroms...
            print, 'Qheating turned ON!!!'
        else
            printf, filenumber, "dustFile 'input/grainspecies.dat' 'input/mrnsmall.dat'"


#all
    printf, filenumber, 'output'

    printf, filenumber, 'contShape  blackbody'
    printf, filenumber, 'nebComposition noGas'
    printf, filenumber, 'maxIterateMC ' + ssi(numiterations[j]) + ' ' + ssf(convpercent[j])
    printf, filenumber, "Ndust file '" + sss(nduststr) + "'"
    printf, filenumber, 'Rin ' + sss(rin[j])
    printf, filenumber, 'Rout ' + sss(rout[j])
    printf, filenumber, 'nx ' + sss(nvalue)#1 was already subtracted if it is a torus
    printf, filenumber, 'ny ' + sss(nvalue)
    printf, filenumber, 'nz ' + sss(nvalue)

#edit viewing angles
    if nova[j] eq 1 then printf, filenumber, 'inclination '+ssi(nova[j])+' '+ssf(theta1[j])+' '+ssf(phi1[j])
    if nova[j] eq 2 then printf, filenumber, 'inclination '+ssi(nova[j])+' '+ssf(theta1[j])+' '+ssf(phi1[j])+' ' + ssf(theta2[j]) + ' ' + ssf(phi2[j])

    close, filenumber
    free_lun, filenumber

def make_grainspecies_amiy(name, file, percent, silicatepercent, AMCpercent):

    get_LUN, filenumber
    openw, filenumber, 'input/grainspecies.dat'














index = where(percent
gt
0, numspecies)
if numspecies eq 0 then BEGIN
if AMCpercent eq 0 then begin
printf, filenumber, ssi(1)
printf, filenumber, "'dustData/sil-oss1.nk' " + sss(silicatepercent / 100.0)
endif else begin
printf, filenumber, ssi(2)
printf, filenumber, "'dustData/sil-oss1.nk' " + sss(silicatepercent / 100.0)
printf, filenumber, "'dustData/amC-hann.nk' " + sss(AMCpercent / 100.0)
if silicatepercent / 100.0 + AMCpercent / 100.0 ne 1 then message, 'error in grainspecies!!!'
;printf, filenumber, "'dustData/sil-dlee.nk' " + sss(silicatepercent[j] / 100.0)
endelse
endif else begin
printf, filenumber, ssi(numspecies);
THIS
NUMBER
MUST
MATCH
THE
NUMBER
OF
dustData
FILES!!!!
for i=0, numspecies-1, 1 do printf, filenumber, "'dustData/"+sss(file[index[i]])+"' "+sss(percent[index[i]] / 100.0)
endelse

printf, filenumber, '';
blank
space
at
the
end!
close, filenumber
free_lun, filenumber

end

pro
makeGrainSizeDistribution, nsizes, slope, amin, amax, outfile, username
;makeGrainSizeDistribution, 20, 3.5, 5e-3, .25, 'mrn_2.dat'
FORWARD_FUNCTION
ssi
FORWARD_FUNCTION
sss
FORWARD_FUNCTION
ssf

if n_params() ne 6 then begin
print, 'how many bins?'
read, nsizes
print, 'what is the exponential giving the slope of the distribution? (p)'
read, slope
print, 'what is amin in microns?'
read, amin
print, 'what is amax in microns?'
read, amax
print, 'what should the output file be called?'
read, outfile
print, 'what is your username??'
read, outfile
endif

astep = (alog10(amax) - alog10(amin)) / float(nsizes - 1)

get_lun, lun
openw, lun, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/input/' + outfile

printf, lun, nsizes, ' sizes '
for i=1, nsizes, 1 do begin
printf, lun, i, 10.0 ^ (alog10(amin) + (i - 1) * astep), (10.0 ^ (alog10(amin) + (i - 1) * astep)) ^ (-slope)
endfor

printf, lun, ''
printf, lun, ''
printf, lun, 'slope (p): ', slope
close, lun
free_lun, lun
print, 'Done making new grain size distribution!!!!!'
end

pro
mocassin_fail_amiy, j, username, diffuse, directoryname, outfoldername, starname
FORWARD_FUNCTION
ssi
FORWARD_FUNCTION
sss
FORWARD_FUNCTION
ssf

print, "RUN FAILED! Writing output."
print, "Failed on line number", j + 1, "of AMIY_input.txt"
restore, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/AMIY_number.sav'
id = ssi(AMIY_number)
AMIY_number + +
save, AMIY_number, filename = '/Users/' + username + '/dust_modelling/mocassin.2.02.70/AMIY_number.sav'

directoryname = 'str'
outfoldername = 'str'
if diffuse[j] then type='SN' else type='RSG'
directoryname = "/Users/" + username + "/dust_modelling/mocassin.2.02.70/output/" + type + "/"$
+id + '_' + starname + '_FAIlED'
spawn, "mkdir " + directoryname
outfoldername = type + "/" + id + '_' + starname + '_FAIlED'

cd_amiy, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/output'
spawn, 'cp dustGrid.out ' + directoryname + '/dustGrid_' + id + '.out.txt'
spawn, 'cp runinfo.txt ' + directoryname + '/runinfo_' + id + '.txt'
spawn, 'cp SED.out ' + directoryname + '/SED_' + id + '.out.txt'
if diffuse[j] then $
spawn, 'cp equivalentTau.out ' + directoryname + '/equivalentTau_' + id + '.out.txt' $
else spawn, 'cp tauNu.out ' + directoryname + '/tauNu_' + id + '.out.txt'
spawn, 'cp /Users/' + username + '/dust_modelling/mocassin.2.02.70/input/input.in ' + directoryname + '/input_' + id + '.in.txt'
spawn, 'cp /Users/' + username + '/dust_modelling/mocassin.2.02.70/input/ndust/nDUST ' + directoryname + '/nDUST_' + id + '.in.txt'

end

function
check_temp_amiy, temperature, temp, starnum, diffuse, errorcheck = errorcheck
;first
check
for user input temp
;next
search
using
starname...
;finally go
to
a
default
temp...
;default
for RSG! AND SN are different!

tstellar = 0
if temperature gt 0 then begin
tstellar = temperature
if keyword_set(errorcheck) then print, "user supplied a temperature, so we are using temperature ", tstellar
endif

if tstellar LE 0 then begin
if starnum gt -1 then begin; if the star is in the list...
tstellar = temp[starnum]
if keyword_set(errorcheck) then print, 'found temperature from file, that temperature is ', tstellar
endif else begin; if not in list
use
a
default
if diffuse then default=6000.0 else default=3500.0
tstellar = default
if keyword_set(errorcheck) then print, "deafulting to temp of ", tstellar
endelse
endif

if starnum gt -1 then begin
if temperature gt 0 and temperature ne temp[starnum] then begin
if keyword_set(errorcheck) then begin
print, "Temp in amiy_input is different than RSG_info.txt (in Bill's folder)"
print, "TStellar= ", Tstellar
print, "temp from RSG_info.txt is ", temp[starnum]
endif
endif
endif

return, tstellar

end

function
check_lstar_amiy, luminosity, luminositycalc, starnum, diffuse, errorcheck = errorcheck

lstar = 0
if luminosity gt 0 then begin
lstar = luminosity
if keyword_set(errorcheck) then print, "user supplied a luminosity, so we are using luminosity ", lstar
endif else begin
if starnum gt -1 then begin
lstar = luminositycalc[starnum]
if keyword_set(errorcheck) then print, 'found luminosity from file, that luminosity is ', lstar
endif else begin;
SN
DEFAULT........RSG
DEFAULT
if diffuse then default=650.0 else default=650.0
if keyword_set(errorcheck) then print, 'Luminosity is defaulting...'
lstar = default
endelse
endelse

return, lstar
end

pro
makenuryd, diffuse, old, errorcheck = errorcheck
c = 2.9979250e10
ryd = 3.28984e15

if file_test('dustData/wavelength_resolution.txt') then begin
readcol, 'dustData/wavelength_resolution.txt', lambda , / silent
nu = c / ( lambda * 1.e-4 * ryd)
index=where( lambda eq .5623, ct)
if ct eq 0 then print, 'FAIL!  You need a lambda of .5623 to get tau out silly pants!'

if diffuse then forprint, nu, textout='dustData/nuRyd.dat', / nocomment, format='F', / silent $
else forprint, nu, textout='dustData/nuDustRyd.dat', / nocomment, format='F', / silent
endif

if old eq 1 then begin; restores the old NU file...
if file_test('dustData/nuDustRyd_default.txt') then begin
readcol, 'dustData/nuDustRyd_default.txt', NU, / silent
;nu = c / (lambda * 1.e-4 * ryd);; ALREADY IN rydbergs.
                  print, 'restoring OLD'
                  if diffuse then forprint, nu, textout='dustData/nuRyd.dat', / nocomment, format='F', / silent $
                  else forprint, nu, textout='dustData/nuDustRyd.dat', / nocomment, format='F', / silent
                  endif
                  endif



                  end

                  pro cd_amiy, new_dir, errorcheck=errorcheck
old_dir = ''
cd, current = old_dir
if keyword_set(errorcheck) then print, 'old dir =', old_dir
if keyword_set(errorcheck) then print, 'new dir =', new_dir
if old_dir ne new_dir then cd, new_dir
end

function
ssf, number
return, strcompress(string(number), / remove_all)
end
function
ssi, number
return, strcompress(string(fix(number)), / remove_all)
end
function
sss, number
return, strcompress(string(number), / remove_all)
end
pro
SUPER_AMIY, infile, dontcheckinput = dontcheckinput, errorcheck = errorcheck

includePAHS = 0;
1
includes
pahs...
0 is
for normal run...
    errorcheck = 0;
    1
    prints
    status
    statements...
    0
    does
    not!
    dontcheckinput = 1;
    1
    prevents
    error
    checking...
    0
    checks
    for errors
        old = 1;
        1 = old
        nuDustRyd.dat, 0
        uses
        wavelength
        from dustdata / wavelength_resolution.txt

    FORWARD_FUNCTION
    makeRSG_AMIY
    FORWARD_FUNCTION
    maketorus_AMIY
    FORWARD_FUNCTION
    mocplot_AMIY2
    FORWARD_FUNCTION
    amiy_check
    FORWARD_FUNCTION
    testbb
    FORWARD_FUNCTION
    get_color
    ;FORWARD_FUNCTION
    ssf
    ;FORWARD_FUNCTION
    ssi
    ;FORWARD_FUNCTION
    sss
    close, / all

    ans = 'str'
    username = 'mocassin'

    cd_amiy, '/Users/' + username + '/dust_modelling/mocassin.2.02.70'

    ;read in the
    input
    from the user

    fmt = 'A,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'
    readcol, 'AMIY_input.txt', format = fmt, starname, n, rin, rout, p, constant, nphotons, numiterations,$
    convpercent, diffuse, temperature, luminosity, silicatepercent, AMCpercent, distance, torus, dens1, dens2, nova, theta1,$
    phi1, theta2, phi2, / silent, comment = '#'

    numrows = n_elements(starname)

    ;;NAME, FILE, PERCENT
    if file_test('AMIY_input_grainspecies.txt') then begin
    readcol, 'AMIY_input_grainspecies.txt', format = 'X,A,A,F,A,A,F,A,A,F,A,A,F',$
    gs1, file1, percent1, gs2, file2, percent2, gs3, file3, percent3, gs4, file4, percent4, / silent
    endif else begin
    percent1 = replicate(0, numrows)
    percent2 = replicate(0, numrows)
    percent3 = replicate(0, numrows)
    percent4 = replicate(0, numrows)
    gs1 = replicate('nothing', numrows)
    gs2 = replicate('nothing', numrows)
    gs3 = replicate('nothing', numrows)
    gs4 = replicate('nothing', numrows)
    file1 = replicate('nothing', numrows)
    file2 = replicate('nothing', numrows)
    file3 = replicate('nothing', numrows)
    file4 = replicate('nothing', numrows)
    endelse

    ;correct
    so
    that if the
    user
    inputs
    a
    file
    ;that
    does
    not have
    enough
    rows, amiy
    wont
    crash
    if n_elements(percent1) lt numrows then begin
    percent1 = replicate(0, numrows)
    percent2 = replicate(0, numrows)
    percent3 = replicate(0, numrows)
    percent4 = replicate(0, numrows)
    gs1 = replicate('nothing', numrows)
    gs2 = replicate('nothing', numrows)
    gs3 = replicate('nothing', numrows)
    gs4 = replicate('nothing', numrows)
    file1 = replicate('nothing', numrows)
    file2 = replicate('nothing', numrows)
    file3 = replicate('nothing', numrows)
    file4 = replicate('nothing', numrows)
    if includePAHS eq 1 then message, 'ERROR, REMEMBER TO CHANGE AMIY_INPUT_GRAINSPECIES.TXT'
    endif

    ;print, 'p1', percent1
    ;percent1 = 0

    ;temperature and luminosity
    are
    optional, set
    to
    0 if unknown
    ;for RSG, luminosity will be derived from the temp and the radius
    ;both
    temp and radius
    will
    be
    stolen
    from lmcsmcfullspec2e.txt
    ;as long as the
    star
    name is correct, it
    will
    find
    the
    temperature, radius, and luminosity.

    ;if AMIY can't find the star's name as an EXACT match to the name in lmcsmcfullspec2e.txt,
    ;then
    it
    will
    try to use the user input temperature and luminosity from AMIY_input.txt
    ;if it STILL can't find these variables, it will default to some value....

    ;temp
    VS
    temperature
    VS
    Tstellar
    ;temperature is what
    the
    USER
    has
    supplied
    us
    with in AMIY_input.txt {THIS IS AN ARRAY!!}
    ;temp is what
    we
    get
    from LMCSMCfullspec2e.txt
    ;tstellar is what
    we finally give
    MOCASSIN as an
    input(
    for both RSG and SN)

    ;luminosity
    Vs.lstar
    ;luminosity is
    from the user

    {THIS
    IS
    AN
    ARRAY}
    ;lstar is what
    we
    give
    to
    MOCASSIN as an
    input
    ;if luminosity is not given, then AMIY will try to calculate it from radius and temperature.

    ;convert
    to
    radians
    theta1 = theta1 * (!pi / 180.)
    theta2 = theta2 * (!pi / 180.)

    fmt = 'A,F,F,F,F'
    readcol, '/Users/mocassin/dust_modelling/accessories/RSG_info.txt', format = fmt, name, temp, apparentlum, dist, mbol, luminositycalc, / silent
    mbol = double(mbol)

    ;begin
    input
    checking...
    fatalerror = 0
    If
    ~KEYWORD_SET(dontcheckinput)
    then
    begin
    AMIY_check, username, fatalerror
    if fatalerror then goto, theend
    ENDIF

    j = 0
    while j lt numrows do begin

    PRINT, "Beginning run number", j + 1, "  for  ", starname[j]

    successtest = 0;
    check
    for if mocassian was successful
    lstar = 0;
    luminosity
    used
    by
    mocassin
    tstellar = 0;
    tstellar
    "     "   "
    starnum = -1;
    number
    of
    the
    star in rsg_info
    name_gs = [gs1[j], gs2[j], gs3[j], gs4[j]];
    grain
    species
    from the file

    percent_gs = [percent1[j], percent2[j], percent3[j], percent4[j]]
    filename_gs = [file1[j], file2[j], file3[j], file4[j]];

    ;set
    dens1 and dens2
    according
    to
    constant and rho in cases
    where
    torus is used and homogenous
    dust
    density
    ;if torus[j] then begin
    ; dens1 = constant / 10000.
    ; if p[j] eq 0 then begin
    ;dens2 = dens1
    ;endif
    ;endif

    ;;search
    for the starnum
        i = 0
    while i lt n_elements(name) and starnum eq -1 do begin
    if name[i] eq starname[j] then starnum=i
    i + +
    endwhile

    tstellar = check_temp_amiy(temperature[j], temp, starnum, diffuse[j], errorcheck=errorcheck)
    lstar = check_lstar_amiy(luminosity[j], luminositycalc, starnum, diffuse[j], errorcheck=errorcheck)
    if keyword_set(errorcheck) then print, 'ts=', tstellar
    if keyword_set(errorcheck) then print, 'ls=', lstar

    ;change
    dust
    percent
    composition
    ;make
    dust
    shell, returns
    the
    name
    of
    the
    ndust
    file...
    nduststr = 'string'
    if torus[j] then readcol, 'nDust.list', nduststr, format='(a)' $
    else nduststr=makeRSG_AMIY(n[j], rin[j], rout[j], p[j], constant[j], username, errorcheck=errorcheck)

    cd_amiy, '/Users/' + username + '/dust_modelling/mocassin.2.02.70'

    n[j] = fix(n[j])

    ;;MAKE
    THE
    GRAINSPECIES
    FILES!!!
    if keyword_set(errorcheck) then print, 'calling grainspecies function'
    make_grainspecies_amiy, name_gs, filename_gs, percent_gs, silicatepercent[j], AMCpercent[j]

    ;;MAKE
    THE
    INPUT.IN
    FILES!!!
    if keyword_set(errorcheck) then print, 'calling input.in making function'
    make_input_amiy, j, torus, diffuse, n, luminosity, tstellar, nphotons,$
    lstar, numiterations, convpercent, rin, rout, nova,$
    theta1, theta2, phi1, phi2, username, nduststr[j], includePAHS = includePAHS

    ;;MAKE
    THE
    WAVELENGTH
    RESOLUTION
    FILE...
    makenuryd, diffuse[j], old

    ;makeGrainSizeDistribution, nsizes, slope, amin, amax, outfile
    ;makeGrainSizeDistribution, 20, 3.5, 5e-3, .25, 'mrn.dat', username
    if includePAHS eq 1 then $
    makeGrainSizeDistribution, 20, 3.5, 3.548e-04, .25, 'PAH_sizes.dat', username;;COMMENT
    IN
    FOR
    PAH
    'S
    ;makeGrainSizeDistribution, 20, 3.5, 3.548e-04, 1.000e+01, 'PAH_sizes.dat', username

    start_time = SYSTIME(1, / seconds)
    runcounter = 0
    while successtest ne 1 and runcounter le 1 do begin
    if runcounter ge 1 then print, "runniing MOCASSIN again!!!"
    print, "beginning to run MOCASSIN!!!"
    spawn, 'mpiexec -n 4 ./mocassin >output/runinfo.txt';
    CMD
    TO
    RUN
    MOCASSIN
    ;spawn, './mocassin >output/runinfo.txt'
    ;tail - f / Users / bill / MOCASSIN / newest_compilation / mocassin
    .2
    .02
    .56 / output / runinfo.txt

    ;test if the
    run
    was
    successful
    spawn, "grep -c '! mocassin: MCIterationDriver done' /Users/" + username + "/dust_modelling/mocassin.2.02.70/output/runinfo.txt > test1.txt"
    readcol, 'test1.txt', format = 'I', successtest, / silent
    if successtest eq 1 then print, 'Mocassin was a great success!!!' else print, 'Mocassin = epic fail!'
    spawn, "rm test1.txt"
    ;spawn, 'mv /Users/mocassin/dust_modelling/mocassin.2.02.70/output/runinfo.txt /Users/mocassin/dust_modelling/mocassin.2.02.70/output/runinfo_' + strtrim(
        string(j + 1), 2) + '.txt'

    if successtest eq 0 then mocassin_fail_amiy, j, username, diffuse, directoryname, outfoldername, starname[j]
    runcounter + +
    endwhile

    if successtest eq 1 then BEGIN
    totaltime = fix(SYSTIME(1, / seconds) - start_time)
    print, 'Total run time is ' + ssi(totaltime) + ' seconds (' + ssi((totaltime) / 60.0) + ')'
    outfoldername = 'str'
    ;    print, 'mocplot_AMIY2,', rin[j], rout[j], constant[j], p[j], lstar, tstellar, starname[j], diffuse[j],$
    ;                  username, outfoldername, starnum, nova[j], theta1[j], phi1[j], theta2[j], phi2[j], torus[j], $
    ;                  silicatepercent[j], AMCpercent[j], name_gs, filename_gs, percent_gs, 1
    if keyword_set(errorcheck) then print, 'ENTERING AMIY MOCPLOT'
    mocplot_AMIY2, rin[j], rout[j], constant[j], p[j], lstar, tstellar, starname[j], diffuse[j], distance[j],$
    username, outfoldername, starnum, nova[j], theta1[j], phi1[j], theta2[j], phi2[j], torus[j], $
    silicatepercent[j], AMCpercent[j], name_gs, filename_gs, percent_gs, totaltime, infile, nduststr[
        j], / residual, errorcheck = errorcheck, includePAHS = includePAHS
    if keyword_set(errorcheck) then print, 'leaving mocplot'

    endif
    cd_amiy, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/'

    ;make
    an
    AMIY
    output
    with what AMIY did.

    get_LUN, filenumber
    openw, filenumber, '/Users/' + username + '/dust_modelling/mocassin.2.02.70/output/' + outfoldername + '/AMIY_output.txt'

    printf, filenumber, "all AMIY_input.txt variables..."

    printf, filenumber, format = '(A,I,2E,2F,E,I,F,2I,F,3I,2F,5I)', starname[j], n[j], rin[j], rout[j], p[j], constant[
        j], nphotons[j], numiterations[j],$
    convpercent[j], diffuse[j], temperature[j], luminosity[j], silicatepercent[j], AMCpercent[j], torus[j],$
    dens1[j], dens2[j], nova[j], theta1[j], phi1[j], theta2[j], phi2[j]

    printf, filenumber, "inputs to mocplot"
    printf, filenumber, format = '(E,E,F,F,I,I,A,I,A,A,I,I,F,F,F,F,I,I,I,4A,4A,4I)', rin[j], rout[j], constant[j], p[
        j], lstar, tstellar, ' ' + starname[j], diffuse[j],$
    ' ' + username, ' ' + outfoldername, starnum, nova[j], theta1[j], phi1[j], theta2[j], phi2[j], torus[j],$
    silicatepercent[j], AMCpercent[j], ' ' + name_gs[0], ' ' + name_gs[1], ' ' + name_gs[2], ' ' + name_gs[3], ' ' +
    filename_gs[0],$
    ' ' + filename_gs[1], ' ' + filename_gs[2], ' ' + filename_gs[3], percent_gs[0], percent_gs[1], percent_gs[2], \
    percent_gs[3]

    printf, filenumber, "starname        ", starname[j]
    printf, filenumber, "n               ", n[j]
    printf, filenumber, "rin             ", rin[j]
    printf, filenumber, "rout            ", rout[j]
    printf, filenumber, "p               ", p[j]
    printf, filenumber, "constant        ", constant[j]
    printf, filenumber, "nphotons        ", nphotons[j]
    printf, filenumber, "numiterations   ", numiterations[j]
    printf, filenumber, "convpercent     ", convpercent[j]
    printf, filenumber, "diffuse         ", diffuse[j]
    printf, filenumber, "temperature     ", temperature[j]
    printf, filenumber, "luminosity      ", luminosity[j]
    printf, filenumber, "silicatepercent ", silicatepercent[j]
    printf, filenumber, "AMCpercent      ", AMCpercent[j]
    printf, filenumber, "inclination     ", nova[j], theta1[j], phi1[j], theta2[j], phi2[j]
    printf, filenumber, "                "
    printf, filenumber, "starnum         ", starnum
    printf, filenumber, "tstellar        ", tstellar
    printf, filenumber, "Lstar           ", lstar
    printf, filenumber, "  "
    printf, filenumber, 'namegs', name_gs
    printf, filenumber, 'filenamegs', filename_gs
    printf, filenumber, 'percent', percent_gs

    close, filenumber
    free_lun, filenumber

    j + +
    endwhile

    theend:
    close, / all
    print, "PROGRAM AMIY IS OVER!"
    print, "PROGRAM AMIY IS OVER!"
    print, "PROGRAM AMIY IS OVER!"
    ;spawn, 'rm /Users/mocassin/dust_modelling/mocassin.2.02.70/output/runinfo.txt'
    if fatalerror then print, "(because of a fatal error!!!)"
    end
