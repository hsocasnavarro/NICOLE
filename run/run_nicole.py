#!/usr/bin/env python

# Script to prepare inputs for use with nicole, run nicole and prepare outputs
#  for human use



# Validate input files. Check for unrecognized entries. First build a
# database of allowed entries by scanning the source file for
# lines of the form variable=get_value(config,'xxxxx',default,file_to_check)
# filestocheck and keysinfile might contain values on input of additional
# allowed keys
# Note: ignores everything after #, = or [
def validate(source,filestocheck,keysinfile,quiet=0):
    import re
    import sys

    badlines=list()
    try:
        length=len(filestocheck)
        findfiles=0
    except:
        findfiles=1
        filestocheck=list()
    try:
        length=len(keysinfile)
    except:
        keysinfile=list()

    #   Compile list of files being read and valid keys
    f=open(source,'r')
    source_lines=f.readlines()
    f.close()

    for l in source_lines:
        l=re.sub('#.*','',l)
        if re.search('(?i)=get_value\(',l):
            l=re.sub("(?i).*=get_value\(",'',l)
#            l=re.sub('\)[^\,*]','',l)
            l=re.sub("\n",'',l)
            l=re.sub("'",'',l)
            l=re.sub('"','',l)
            l=re.sub('\)$','',l)
            filename=l.split(',')[3]
            key=l.split(',')[1].lower()
            
            try:
                ifile=filestocheck.index(filename)
                if len(keysinfile[ifile]) == 0: keysinfile[ifile]=list()
                keysinfile[ifile].append(key)
            except:
                if (findfiles == 1):
                    ifile=len(filestocheck)
                    filestocheck.append(filename)
                    keysinfile.append('')
                    keysinfile[ifile]=list()
                    keysinfile[ifile].append(key)


    #   Loop through files and check lines against keys

    for ifile in range(len(filestocheck)):
        if quiet != 1: print 'Checking syntax in file:'+filestocheck[ifile]
        errors=0
        try:
            f=open(filestocheck[ifile],'r')
            source_lines=f.readlines()
            f.close()
        except:
            source_lines=''
        indline=0
        for line in source_lines:
            indline=indline+1
            l=line.strip()
            l=l.lower()
            l=re.sub('#.*','',l)
            l=re.sub('=.*','',l)
            l=re.sub('^ *','',l)
            l=re.sub(' *$','',l)
            l=re.sub('\[.*','',l)
            l=re.sub('\n','',l)
            if l != '':
                try:
                    ind=keysinfile[ifile].index(l)
                except:
                    badlines.append('\nError in file '+filestocheck[ifile]+'. Offending line ['+str(indline)+']:'+line)
                    errors=1
                    import difflib
                    matches=difflib.get_close_matches(l,keysinfile[ifile])
                    if (len(matches) > 0):
                        if (len(matches) == 1):
                            badlines.append('Did you mean the following instead?')
                        else:
                            badlines.append('Did you mean one of the following instead?')
                        for suggestion in matches:
                            badlines.append(' '+suggestion+'=')

#        if errors == 0:
#            print '  ... no errors found'

    return badlines


#        l=re.sub('#$','',l)
#        l=re.sub('[$','',l)
 #       l=re.sub('=$','',l)

# Get value from configuration file or stop with an error if not found
# (unless default is specified, in which case the default is returned)
def get_value(config,key,default,filename,section='',subsection=''):
    import sys
    import re

    try:
        if section == '':
            value=config[key.lower()]
        elif subsection == '':
            value=config[section.lower()][key.lower()]
        else:
            value=config[section.lower()][subsection.lower()][key.lower()]
    except KeyError:
        if default == '':
            readstr=''
            if section != '': readstr=readstr+', section: \"'+section+'\"'
            if subsection != '': readstr=readstr+', subsection: \"'+subsection+'\"'
            print 'Error in ',filename+readstr
            print 'Required field \"'+key+'\" was not found'
            sys.exit(1)
        else:
            value=default
    if value != None:
        value=re.sub("^'",'',value)
        value=re.sub("'$",'',value)
        value=re.sub('^"','',value)
        value=re.sub('"$','',value)
    return value
# Convert string to lower case up to the first occurrence of a separator
def lower_to_sep(string, separator='='):
    list=string.split(separator)
    nlist=len(list)
    if nlist == 1: 
        return string.lower()
    retstring=list[0].lower()
    for i in range(nlist-1):
        retstring=retstring+'='+list[i+1]
    return retstring

# Main program
from configobj import ConfigObj
import sys
import os
import re
import struct
import subprocess
import getopt
from model_prof_tools import *

elements=['h','he','li','be','b','c','n','o','f','ne', 
          'na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr', 
          'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr', 
          'rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in', 
          'sn','sb','te','i','xe','cs','ba','la','ce','pr','nd','pm', 
          'sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w', 
          're','os','ir','pt','au','hg','tl','pb','bi','po','at','rn', 
          'fr','ra','ac','th','pa','u']

# Get command-line arguments
nicolecommand_argv=''
inputmodel_argv=''
outprof_argv=''
try:
    opts, args = getopt.getopt(sys.argv[1:],"", ["nicolecommand=","modelin=","profout="])
except:
    print 'Command-line option not recognized'
    print 'Usage:'
    print "run_nicole.py [--nicolecommand='command' --modelin='inputmodel' --profout='outputprofiles']"
    sys.exit(2)
for o, a in opts:
    if o == '--nicolecommand':
        nicolecommand_argv=a
    if o == '--modelin':
        inputmodel_argv=a
    if o == '--profout':
        outprof_argv=a

if ((inputmodel_argv <> '') and (outprof_argv <> '')): print 'Command-line options, model ', inputmodel_argv,', profiles output ',outprof_argv

# Check sizes of default types
[int4f,intf,flf]=check_types()
 
# Validate LINES file
badlines=validate('run_nicole.py',['LINES'],[''])
#
# Parse first NICOLE.input
#
if os.path.isfile('NICOLE.input_1'):
    f=open('NICOLE.input_1','r')
else:
    f=open('NICOLE.input','r')
NICOLE_input_lines=f.readlines()
f.close()
NICOLE_input_lower=['']
for l in NICOLE_input_lines:
    NICOLE_input_lower.append(lower_to_sep(l)) # Convert keys to lowercase
config=ConfigObj(NICOLE_input_lower)
ncycles=get_value(config,'Cycles','1','NICOLE.input')
cycle0=get_value(config,'Start cycle','1','NICOLE.input')

for icycle in range(int(ncycles)):
    files=list()
    keys=['']
    keys[0]=list()
    keys[0][:]=elements

    if ncycles > 1: print 'Preparing cycle',icycle+1
    suffix='_'+str(icycle+1)
    if os.path.isfile('NICOLE.input'+suffix):
        f=open('NICOLE.input'+suffix,'r')
        files2=list(files) # Make a copy
        files2.append('NICOLE.input') # Use NICOLE.input to get keys
        void=validate('run_nicole.py',files2,keys,quiet=1) # get keys
        files.append('NICOLE.input'+suffix)
        badlines=badlines+validate('run_nicole.py',files,keys)
    else:
        f=open('NICOLE.input','r')
        files.append('NICOLE.input')
        badlines=badlines+validate('run_nicole.py',files,keys)
    if len(badlines) >= 1:
        print 'The following errors have been found in the input files'
        for l in badlines: print l
        sys.exit(1)
    else:
        print '  ... no errors found\n'
    NICOLE_input_lines=f.readlines()
    f.close()
    NICOLE_input_lower=['']
    for l in NICOLE_input_lines:
        NICOLE_input_lower.append(lower_to_sep(l)) # Convert keys to lowercase
    config=ConfigObj(NICOLE_input_lower)
    if nicolecommand_argv == '':
        nicolecommand=get_value(config,'Command','../main/nicole','NICOLE.input')
    else:
        nicolecommand=nicolecommand_argv
    mode=get_value(config,'Mode','','NICOLE.input')
#    ncycles=get_value(config,'Cycles','1','NICOLE.input')
    if inputmodel_argv == '':
        inputmodel=get_value(config,'Input model','','NICOLE.input')
    else:
        inputmodel=inputmodel_argv
    inputmodel2=get_value(config,'Input model 2',None,'NICOLE.input')
    if outprof_argv == '':
        outprof=get_value(config,'Output profiles','','NICOLE.input')
    else:
        outprof=outprof_argv
    helio=get_value(config,'Heliocentric angle','None','NICOLE.input')
    xang=get_value(config,'Heliocentric X-angle','None','NICOLE.input')
    yang=get_value(config,'Heliocentric Y-angle','None','NICOLE.input')
    dx=get_value(config,'X scale','1E10','NICOLE.input')
    dy=get_value(config,'Y scale','1E10','NICOLE.input')
    obsprof=get_value(config,'Observed profiles',None,'NICOLE.input')
    outputmodel=get_value(config,'Output model',None,'NICOLE.input')
    outputmodel2=get_value(config,'Output model 2',None,'NICOLE.input')
    formal=get_value(config,'Formal solution method','0','NICOLE.input')
    boundarycond=get_value(config,'Formal solution boundary condition','Normal','NICOLE.input')
    stray=get_value(config,'Stray light file',None,'NICOLE.input')
    printout=get_value(config,'Printout detail','1','NICOLE.input')
    maxiters=get_value(config,'Maximum inversion iterations','25','NICOLE.input')
    noise=get_value(config,'Noise level','1e-3','NICOLE.input')
    maxinv=get_value(config,'Maximum number of inversions','5','NICOLE.input')
    acceptchisq=get_value(config,'Acceptable Chi-square','-1','NICOLE.input')
    if acceptchisq == -1:
        acceptchisq=get_value(config,'Acceptable Chisquare','-1','NICOLE.input')
    if acceptchisq == -1:
        acceptchisq=get_value(config,'Acceptable Chi square','-1','NICOLE.input')
    if acceptchisq == -1:
        acceptchisq=get_value(config,'Acceptable Chi-sq','-1','NICOLE.input')
    if acceptchisq == -1:
        acceptchisq=get_value(config,'Acceptable Chi2','5','NICOLE.input')
    speed=get_value(config,'Speed optimization','1','NICOLE.input')
    always_compute_der=get_value(config,'Always compute derivatives','Y','NICOLE.input')
    cent_der=get_value(config,'Centered derivatives','0','NICOLE.input')
    gravity=get_value(config,'Gravity','2.7414e+4','NICOLE.input')
    regul=get_value(config,'Regularization','1.e-3','NICOLE.input')
    update_opac=get_value(config,'Update opacities every','10','NICOLE.input')
    negligible_opac=get_value(config,'Negligible opacity','0','NICOLE.input')
    contref=get_value(config,'Continuum reference','1','NICOLE.input')
    contval=get_value(config,'Continuum value','1','NICOLE.input')
    sethydro=get_value(config,'Impose hydrostatic equilibrium','Y','NICOLE.input')
    setnH=get_value(config,'Compute Hydrogen populations','Y','NICOLE.input')
    depcoef=get_value(config,'Depcoef behavior','1','NICOLE.input')
    write_depcoef=get_value(config,'Write depcoef','N','NICOLE.input')
    inputdens=get_value(config,'Input density','Pel','NICOLE.input')
    keep_el_p=get_value(config,'Keep El_p',' -1 ','NICOLE.input')
    keep_gas_p=get_value(config,'Keep Gas_p',' -1 ','NICOLE.input')
    keep_rho=get_value(config,'Keep Rho',' -1 ','NICOLE.input')
    keep_nH=get_value(config,'Keep nH',' -1 ','NICOLE.input')
    keep_nHminus=get_value(config,'Keep nHminus',' -1 ','NICOLE.input')
    keep_nHplus=get_value(config,'Keep nHplus',' -1 ','NICOLE.input')
    keep_nH2=get_value(config,'Keep nH2',' -1 ','NICOLE.input')
    keep_nh2plus=get_value(config,'Keep nh2plus',' -1 ','NICOLE.input')
    restart=get_value(config,'Restart','-1','NICOLE.input')
    inputdens=inputdens.lower()
    hscale=get_value(config,'Height scale','t','NICOLE.input')
    hscale=hscale[0]
    hscale=hscale.lower()
    opacities=get_value(config,'Opacity package','wittmann','NICOLE.input')
    opacities=opacities.lower()
    opacitiesUV=get_value(config,'Opacity package UV','-','NICOLE.input')
    opacitiesUV=opacitiesUV.lower()
    if opacitiesUV == '-':
        if opacities == 'sopa' or opacities == 'sopas' or opacities == 'shchukina': 
            opacitiesUV = 'sopa'
        else:
            opacitiesUV = 'top'
    eqstate=get_value(config,'Eq state','0','NICOLE.input')
    eqstate=eqstate.lower()
    if (eqstate == '0' or eqstate == 'nicole'):
        eqstate = '0'
    elif (eqstate == '1' or eqstate == 'ann' or eqstate == 'simpleann'):
        eqstate = '1'
    elif (eqstate == '2' or eqstate == 'wittmann' or eqstate == 'wittman'):
        eqstate = '2'
    else:
        print 'Unknown value for "Eq state"'
        print 'Must be NICOLE, ANN or Wittmann'
        sys.exit(1)
    eqstateH=get_value(config,'Eq state for H','0','NICOLE.input')
    eqstateH=eqstateH.lower()
    if (eqstateH == '0' or eqstateH == 'nicole'):
        eqstateH = '0'
    elif (eqstateH == '1' or eqstateH == 'asensio' or eqstateH == 'andres' or 
          eqstateH == 'asensio2' or eqstateH == 'andres2'):
        eqstateH = '1'
    elif (eqstateH == '2' or eqstateH == 'asensio273' or eqstateH == 'andres273'):
        eqstateH = '2'
    elif (eqstateH == '3' or eqstateH == 'wittmann' or eqstateH == 'wittman'):
        eqstateH = '3'
    else:
        print 'Unknown value for "Eq state for H"'
        print 'Must be NICOLE, Asensio, Asensio273 or Wittmann'
        sys.exit(1)
    peconsistency=get_value(config,'Pe consistency','1e-1','NICOLE.input')
    debug=get_value(config,'Debug mode','0','NICOLE.input')
    interp=get_value(config,'Optimize grid','-','NICOLE.input')
    outputpop=get_value(config,'Output populations','0','NICOLE.input')
    outputcontop=get_value(config,'Output continuum opacity','0','NICOLE.input')
    outputNLTEsf=get_value(config,'Output NLTE source function','0','NICOLE.input')
    pixx1=get_value(config,'Start X position','1','NICOLE.input')
    pixy1=get_value(config,'Start Y position','1','NICOLE.input')
    pixx2=get_value(config,'End X position','10000000','NICOLE.input')
    pixy2=get_value(config,'End Y position','10000000','NICOLE.input')
    irec0=get_value(config,'Start irec','1','NICOLE.input')

    if stray == None: stray=''
    if obsprof == None: obsprof=''
    if inputmodel2 == None: inputmodel2=''
    if outputmodel == None: outputmodel=''
    if outputmodel2 == None: outputmodel2=''
    if len(inputmodel) > 256 :
        print inputmodel,' is too long. Must be under 256 characters'
    if len(inputmodel2) > 256 :
        print inputmodel2,' is too long. Must be under 256 characters'
    if len(outprof) > 256 :
        print outprof,' is too long. Must be under 256 characters'
    if len(obsprof) > 256 :
        print obsprof,' is too long. Must be under 256 characters'
    if len(outputmodel) > 256 :
        print outputmodel,' is too long. Must be under 256 characters'
    if len(outputmodel2) > 256 :
        print outputmodel2,' is too long. Must be under 256 characters'
    if len(stray) > 256 :
        print stray,' is too long. Must be under 256 characters'

    mode=mode.lower()
    mode=mode[0:1]
    if mode != 's' and mode != 'i' and mode != 'c':
        print 'Error in NICOLE.input'
        print 'Unknown mode. Must be either Synthesis, Inversion or Convert'
        sys.exit(1)
    if inputdens != 'pel' and inputdens != 'pgas' and inputdens != 'nel' and inputdens != 'dens':
        print 'Error in NICOLE.input. Input density must be either Pel, Pgas, Nel or Dens'
        sys.exit(1)
    if formal < '0' or formal > '9':
        print 'Error in NICOLE.input. Invalid formal solution method'
        sys.exit(1)
    boundarycond=boundarycond.lower()
    if boundarycond == 'normal' or boundarycond == 'difussion':
        boundarycond = '0'
    elif boundarycond == 'zero':
        boundarycond = '1'
    else:
        print 'Error in NICOLE.input. Invalid formal solution boundary condition'
        print 'Must be Normal, Difusion or Zero'
        sys.exit(1)    
    if hscale != 't' and hscale != 'z':
        print 'Error in NICOLE.input. Input density must be either tau or z'
        sys.exit(1)
    if opacities == 'natasha' or opacities == 'sopas' or opacities == 'shchukina': opacities = 'sopa'
    if opacities == 'andres' or opacities == 'asensio': opacities = 'andres'
    if opacities != 'andres' and opacities != 'sopa' and opacities != 'wittmann':
        print 'Error in NICOLE.input. Unknown opacity package:',opacities
        sys.exit(1)
    if opacitiesUV != 'top' and opacitiesUV != 'dm' and opacitiesUV != 'sopa':
        print 'Error in NICOLE.input. Unknown opacity package UV:',opacitiesUV
        sys.exit(1)
    if opacities == 'wittmann': opacities = '1'
    if opacities == 'andres': opacities = '2'
    if opacities == 'sopa': opacities = '3'
    if opacitiesUV == 'top': opacitiesUV = '1'
    if opacitiesUV == 'dm': opacitiesUV = '2'
    if opacitiesUV == 'sopa': opacitiesUV = '3'
    if opacities == '3' and opacitiesUV != '3':
        print 'If you select SOPA as opacity package, you cannot select a different'
        print 'opacity package for the UV. Leave it blank or set it to SOPA as well'
        sys.exit(1)
    always_compute_der=always_compute_der.lower()
    always_compute_der=always_compute_der[0:1]
    if always_compute_der != 'y' and always_compute_der != 'n':
        print 'Error in NICOLE.input'
        print 'Always compute derivatives must be either Y or N'
        sys.exit(1)
    if always_compute_der == 'y': 
        always_compute_der='1' 
    else: 
        always_compute_der='0'
    interp=interp.lower()
    interp=interp[0:1]
    if interp == '-': # Default
        if mode == 's': interp='1'
        if mode == 'i' or mode == 'c': interp='0'
    if interp == 'y': interp='1'
    if interp == 'n': interp='0'
    sethydro=sethydro.lower()
    sethydro=sethydro[0:1]
    if sethydro != 'y' and sethydro != 'n':
        print 'Error in NICOLE.input'
        print 'Impose hydrostatic equilibrium must be either Y or N'
        sys.exit(1)
    if sethydro == 'y': 
        sethydro = 'T'
    else: 
        sethydro = 'F'
    if restart != '-1' and restart != '0' and restart !='1':
        print 'Error in NICOLE.input'
        print 'Restart must be either -1, 0 or 1'
        sys.exit(1)
    setnH=setnH.lower()
    setnH=setnH[0:1]
    if setnH != 'y' and setnH != 'n':
        print 'Error in NICOLE.input'
        print 'Compute Hydrogen populations must be either Y or N'
        sys.exit(1)
    if setnH == 'y': 
        setnH = 'T'
    else: 
        setnH = 'F'
    if depcoef != '1' and depcoef != '2':
        print 'Incorrect choice for depcoef behavior:',depcoef
        print 'Must be 1 or 2'
        sys.exit(1)
    write_depcoef=write_depcoef[0:1]
    write_depcoef=write_depcoef.lower()
    if write_depcoef != 'y' and write_depcoef != 'n':
        print 'Incorrect "Write depcoef":',depcoef
        print 'Must be Y or N'
        sys.exit(1)
    if write_depcoef == 'y':
        write_depcoef = 'T'
    else:
        write_depcoef = 'F'
    if mode == 'i' and (obsprof == None or outputmodel == None):
        print 'Error in NICOLE.input'
        print 'When Mode is Inversion, the fields "Observed profiles" and'
        print '"Ouptut model" are mandatory'
        sys.exit(1)
    # Parse spectral information
    nlines=0
    line=''
    while line != None:
        nlines=nlines+1
        line=get_value(config,'Line',None,'NICOLE.input','Line '+str(nlines),)
    nlines=str(nlines-1)
    nregions=0
    region=''
    while region != None:
        nregions=nregions+1
        region=get_value(config,'First wavelength',None,'NICOLE.input','Region '+str(nregions),)
    nregions=str(nregions-1)
    nlam=0
    for iregion in range(int(nregions)):
        tmp3=get_value(config,'Number of wavelengths','','NICOLE.input'
                      ,'Region '+str(iregion+1))
        nlam=nlam+int(tmp3)
    if (nlam == 0):
        print 'Error in NICOLE.input. At least one wavelength is required'
        sys.exit(1)
    # 
    # Parse abundance data now because it may be needed
    #
    abund_set=get_value(config,'Abundance set','grevesse_sauval_1998','NICOLE.input','Abundances')
    abund_set=abund_set.strip()
    elements=['H','HE','LI','BE','B','C','N','O','F','NE', 
              'NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR', 
              'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR', 
              'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN', 
              'SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM', 
              'SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W', 
              'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN', 
              'FR','RA','AC','TH','PA','U']
    for i in range(len(elements)):
        elements[i]=elements[i].lower()
    if abund_set.lower() == 'grevesse_1984':
        abundances=[ 12.00,11.00,1.00,1.15,2.60,8.55,7.99,8.77,4.56,8.00,6.18, 
                  7.48, 6.40,7.55,5.45,7.21,5.50,6.58,5.12,6.36,3.10,5.02, 
                  4.00, 5.67,5.45,7.50,4.92,6.25,4.21,4.60,2.88,3.63,2.39, 
                  3.35, 2.63,3.21,2.60,2.90,2.24,2.56,2.10,1.92,0.00,1.84, 
                  1.12, 1.69,0.94,1.86,1.66,2.00,1.00,2.25,1.51,2.19,1.12, 
                  2.13, 1.22,1.55,0.71,1.34,0.00,0.80,0.51,1.12,0.20,1.10, 
                  0.26, 0.93,0.00,1.08,0.76,0.88,-.09,1.11,0.26,1.45,1.35, 
                  1.80, 1.13,1.27,0.90,1.90,0.71,-8.0,-8.0,-8.0,-8.0,-8.0, 
                  -8.00, 0.02,-8.0 ,-.47]
    elif abund_set.lower() == 'thevenin_1989': # Actually, taken from Holweger 1979
        abundances=[ 12.00, 11.00, 1.00, 1.15, 2.60, 8.69, 7.99, 8.91, 4.56, 8.00, 6.280,
            7.53, 6.43, 7.50, 5.45, 7.21, 5.50, 6.58, 5.05, 6.36, 2.99, 4.880,
            3.91, 5.61, 5.47, 7.46, 4.85, 6.18, 4.24, 4.60, 2.88, 3.57, 2.390,
            3.35, 2.63, 3.21, 2.60, 2.93, 2.18, 2.46, 1.46, 2.10, 0.00, 1.780,
            1.10, 1.69, 0.94, 1.86, 1.66, 2.00, 1.00, 2.25, 1.51, 2.19, 1.120,
            2.18, 1.07, 1.58, 0.76, 1.40, 0.00, 0.88, 0.48, 1.13, 0.20, 1.070,
            0.26, 0.93, 0.00, 1.08, 0.76, 0.88, -0.09, 0.98, 0.26, 1.45, 1.360,
            1.80, 1.13, 1.27, 0.90, 1.90, 0.71, -8.00, -8.00, -8.00, -8.00, -8.00,
           -8.00, 0.02, -8.00, -0.470]
    elif abund_set.lower() == 'asplund_et_al_2009':
        abundances=[ 12.00, 10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
                     6.24,7.60,6.45,7.51,5.41,7.12,5.50,6.4,5.03,6.34,3.15,4.95,
                     3.93,5.64,5.43,7.50,4.99,6.22,4.19,4.56,3.04,3.65,-10,-10,-10,
                     3.25,2.52,2.87,2.21,2.58,1.46,1.88,-10,
                     1.75,0.91,1.57,0.94,-10,0.80,2.04,-10,-10,-10,2.24,-10,2.18,
                     1.10,1.58,0.72,1.42,-10,0.96,0.52,1.07,0.30,1.10,0.48,0.92,
                     0.10,0.84,0.10,0.85,-10,0.85,-10,1.40,1.38,-10,0.92,-10,0.90,
                     1.75,-10,-10,-10,-10,-10,-10,-10,.02,-10,-10]
    elif abund_set.lower() == 'grevesse_sauval_1998':
        abundances=[ 12.00,10.93,1.10,1.40,2.55,8.52,7.92,8.83,4.56,8.08,6.33,7.58,
                     6.47,7.55,5.45,7.33,5.5,6.40,5.12,6.36,3.17,5.02,4.00,5.67,
                     5.39,7.50,4.92,6.25,4.21,4.60,2.88,3.41,-10,-10,-10,-10,
                     2.60,2.97,2.24,2.60,1.42,
                     1.92,-10,1.84,1.12,1.69,0.94,1.77,1.66,2.0,1.0,-10,-10,-10,-10,
                     2.13,1.17,1.58,0.71,1.50,-10,1.01,0.51,1.12,-0.1,1.14,0.26,
                     0.93,0.00,1.08,0.06,0.88,-10,1.11,-10,1.45,1.35,1.8,1.01,-10,
                     0.9,1.95,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10]
    elif abund_set.lower() != 'model':
        print 'Abundance set not known in NICOLE.input, section Abundances'
        print 'Must be one of Grevesse_Sauval_1998, Asplund_et_al_2009,'
        print 'Thevenin_1989, or Grevesse_1984'
        sys.exit(1)
    if abund_set.lower() == 'model': 
        abundances=[0]*len(elements)
        abundances[0]=-1 # Flag to use abundance set in the model
    if len(abundances) != len(elements):
        print 'Error in abundances section of run_nicole.py'
        sys.exit(1)
    abundances=dict(zip(elements,abundances))
    abund_file=get_value(config,'Abundance file',None,'NICOLE.input','Abundances')
    if abund_file != None:
        print 'Reading abundance file:',abund_file
        try:
            fa=open(abund_file,'r')
        except IOError:
            print 'Error reading abundance file:',abund_file
            print 'in NICOLE.input, section: Abundances'
            sys.exit(1)
        ab_lines=fa.readlines()
        fa.close()
        ab_lines_lower=['']
        for l in ab_lines:
            ab_lines_lower.append(lower_to_sep(l))
        config_ab=ConfigObj(ab_lines_lower)
        ablist=get_value(config_ab,'abundances',None,abund_file)
        if ablist != None:
            if len(ablist) != 92:
                print 'Error in ',abundfile
                print 'Need 92 values for abundances. Found:',len(ablist)
                sys.exit(1)
            for iel in range(len(elements)): abundances[elements[iel]]=ablist[iel]
        for element in abundances:
            abundances[element]=get_value(config_ab,element,str(abundances[element]),abund_file)
    for element in elements:
        override=get_value(config,element,str(abundances[element]),'NICOLE.input','Abundances','Override')
        if re.search('(?i)ppm',override): # Convert ppm to log scale
            ppm=re.sub('(?i)ppm','',override) # Remove text
            import math
            override=math.log10(float(ppm)*1e-6) + 12
        abundances[element]=override
    numab=list()
    for el in elements: numab.append(float(abundances[el]))
    # Elements to ignore in background opac
    elneglectopac=get_value(config,'Elements to ignore in background opacities',' ','NICOLE.input')
    elneglectopac=elneglectopac.lower()
    elneglectopac=re.sub(',',' ',elneglectopac)
    elneglectopac=elneglectopac.split()
    elneglectopacidx=list()
    for el in elements:
        el=el.lower()
        if elneglectopac.count(el) >= 1:
            elneglectopacidx.append(int([elements.index(el)][0])+1)
    #
    # NLTE
    #
    elim1=get_value(config,'Elim','-','NICOLE.input','NLTE')
    if elim1 == '-':
        elim1=get_value(config,'Populations relative change','-','NICOLE.input','NLTE')
    if elim1 == '-': 
        elim1='1e-3'
        if mode == 'i': elim1='1e-4'
    isum=get_value(config,'isum','1','NICOLE.input','NLTE')
    istart=get_value(config,'istart','1','NICOLE.input','NLTE')
    cper=get_value(config,'cper','1.0','NICOLE.input','NLTE')
    usecolswitch=get_value(config,'Use collisional switching','-','NICOLE.input','NLTE')
    usecolswitch=usecolswitch.lower()
    if usecolswitch == '-': usecolswitch='n'
    if usecolswitch != 'y' and usecolswitch != 'n':
        print 'Error, unknown value for Use Collisional Switching'
        sys.exit(1)
    if usecolswitch == 'y':
        usecolswitch = '1'
    else:
        usecolswitch = '0'
    nmu=get_value(config,'nmu','-','NICOLE.input','NLTE')
    if nmu == '-': nmu='3'
    qnorm=get_value(config,'qnorm','10','NICOLE.input','NLTE')
    nlteformalsolution=get_value(config,'formal solution','-','NICOLE.input','NLTE')
    if nlteformalsolution == '-':
        nlteformalsolution=get_value(config,'nlte formal solution','-','NICOLE.input','NLTE')
    if nlteformalsolution=='-': nlteformalsolution='1'
    if nlteformalsolution != '1' and nlteformalsolution != '2':
        print 'Error in NLTE Formal Solution'
        sys.exit(1)
    velfree=get_value(config,'vel free','-','NICOLE.input','NLTE')
    if velfree == '-':
        velfree=get_value(config,'velocity free','-','NICOLE.input','NLTE')
    if velfree == '-': velfree='y'
    if velfree != 'y' and velfree != 'n':
        print 'Error, unknown value for Velocity Free'
        sys.exit(1)
    if velfree == 'y':
        velfree = 'T'
    else:
        velfree = 'F'
    ngacc=get_value(config,'ngacc','-','NICOLE.input','NLTE')
    if ngacc == '-':
        ngacc=get_value(config,'ng acc','-','NICOLE.input','NLTE')
    if ngacc == '-':
        ngacc=get_value(config,'ng acceleration','-','NICOLE.input','NLTE')
    if ngacc == '-': ngacc='y'
    if ngacc != 'y' and ngacc != 'n':
        print 'Error, unknown value for NG Acceleration'
        sys.exit(1)
    if ngacc == 'y':
        ngacc = 'T'
    else:
        ngacc = 'F'
    nlteoptthin=get_value(config,'optically thin','1e-3','NICOLE.input','NLTE')
    nlteoptthick=get_value(config,'optically thick','1e3','NICOLE.input','NLTE')
    nltelinear=get_value(config,'linear formal solution','0','NICOLE.input','NLTE')
    nltemaxiters=get_value(config,'max iters','-','NICOLE.input','NLTE')
    if nltemaxiters == '-': nltemaxiters='500'
    lambdaiters=get_value(config,'lambda iterations','-','NICOLE.input','NLTE')
    if lambdaiters == '-': 
        lambdaiters=get_value(config,'number of lambda iterations','-','NICOLE.input','NLTE')    
    if lambdaiters == '-': lambdaiters='3'
    ltepop=get_value(config,'ltepop','nicole','NICOLE.input','NLTE')
    ltepop=ltepop.lower()
    ltepop=ltepop[0:1]
    if ltepop != 'n' and ltepop != 'm':
        print 'Unknown value for ltepop in LTE section of NICOLE.input'
        print 'Needs to be either NICOLE or MULTI'
        sys.exit(1)
    if ltepop == 'n': 
        ltepop='1' 
    else: 
        ltepop='2'
    #
    # Nodes
    #
    nodesT=get_value(config,'Temperature','-1','NICOLE.input','Nodes')
    if nodesT == '-1': nodesT=get_value(config,'T','-1','NICOLE.input','Nodes')
    if nodesT == '-1': nodesT=get_value(config,'Temp','-1','NICOLE.input','Nodes')
    nodesTx=''
    nod=nodesT
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesT=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesTx=locnod
    nodesv=get_value(config,'Velocity','-1','NICOLE.input','Nodes')
    if nodesv == '-1': nodesv=get_value(config,'v','-1','NICOLE.input','Nodes')
    if nodesv == '-1': nodesv=get_value(config,'vel','-1','NICOLE.input','Nodes')
    if nodesv == '-1': nodesv=get_value(config,'vlos','-1','NICOLE.input','Nodes')
    if nodesv == '-1': nodesv=get_value(config,'v_los','-1','NICOLE.input','Nodes')
    nodesvx=''
    nod=nodesv
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesv=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesvx=locnod
    nodesvmic=get_value(config,'Microturbulence','-1','NICOLE.input','Nodes')
    if nodesvmic == '-1': nodesvmic=get_value(config,'Micro','-1','NICOLE.input','Nodes')
    if nodesvmic == '-1': nodesvmic=get_value(config,'vmic','-1','NICOLE.input','Nodes')
    if nodesvmic == '-1': nodesvmic=get_value(config,'v_mic','-1','NICOLE.input','Nodes')
    nodesvmicx=''
    nod=nodesvmic
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesvmic=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesvmicx=locnod
    nodesvmac=get_value(config,'Macroturbulence','-1','NICOLE.input','Nodes')
    if nodesvmac == '-1': nodesvmac=get_value(config,'Macro','-1','NICOLE.input','Nodes')
    if nodesvmac == '-1': nodesvmac=get_value(config,'vmac','-1','NICOLE.input','Nodes')
    if nodesvmac == '-1': nodesvmac=get_value(config,'v_mac','-1','NICOLE.input','Nodes')
    nodesblong=get_value(config,'B_long','-1','NICOLE.input','Nodes')
    if nodesblong == '-1': nodesblong=get_value(config,'Blong','-1','NICOLE.input','Nodes')
    if nodesblong == '-1': nodesblong=get_value(config,'Bl','-1','NICOLE.input','Nodes')
    if nodesblong == '-1': nodesblong=get_value(config,'B_z','-1','NICOLE.input','Nodes')
    if nodesblong == '-1': nodesblong=get_value(config,'Bz','-1','NICOLE.input','Nodes')
    nodesblongx=''
    nod=nodesblong
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesblong=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesblongx=locnod
    nodesbx=get_value(config,'B_x','-1','NICOLE.input','Nodes')
    if nodesbx == '-1': nodesbx=get_value(config,'Bx','-1','NICOLE.input','Nodes')
    nodesbxx=''
    nod=nodesbx
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesbx=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesbxx=locnod
    nodesby=get_value(config,'B_y','-1','NICOLE.input','Nodes')
    if nodesby == '-1': nodesby=get_value(config,'By','-1','NICOLE.input','Nodes')
    nodesbyx=''
    nod=nodesby
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesby=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesbyx=locnod
    nodesstray=get_value(config,'Stray','-1','NICOLE.input','Nodes')
    if nodesstray == '-1': nodesstray=get_value(config,'Stray light','-1','NICOLE.input','Nodes')
    nodesexp=get_value(config,'ffactor','-1','NICOLE.input','Nodes')
    if nodesexp == '-1': nodesexp=get_value(config,'filling factor','-1','NICOLE.input','Nodes')
    nodesab=get_value(config,'Abundances','-1','NICOLE.input','Nodes')
    if nodesab == '-1': nodesab=get_value(config,'Number of abundances','-1','NICOLE.input','Nodes')
# Nodes for the second component
    nodesT2=get_value(config,'Temperature','-1','NICOLE.input','Nodes 2')
    if nodesT2 == '-1': nodesT2=get_value(config,'T','-1','NICOLE.input','Nodes 2')
    if nodesT2 == '-1': nodesT2=get_value(config,'Temp','-1','NICOLE.input','Nodes 2')
    nodesTx2=''
    nod=nodesT2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesT2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesTx2=locnod
    nodesv2=get_value(config,'Velocity','-1','NICOLE.input','Nodes 2')
    if nodesv2 == '-1': nodesv2=get_value(config,'v','-1','NICOLE.input','Nodes 2')
    if nodesv2 == '-1': nodesv2=get_value(config,'vel','-1','NICOLE.input','Nodes 2')
    if nodesv2 == '-1': nodesv2=get_value(config,'vlos','-1','NICOLE.input','Nodes 2')
    if nodesv2 == '-1': nodesv2=get_value(config,'v_los','-1','NICOLE.input','Nodes 2')
    nodesvx2=''
    nod=nodesv2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesv2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesvx2=locnod
    nodesvmic2=get_value(config,'Microturbulence','-1','NICOLE.input','Nodes 2')
    if nodesvmic2 == '-1': nodesvmic2=get_value(config,'Micro','-1','NICOLE.input','Nodes 2')
    if nodesvmic2 == '-1': nodesvmic2=get_value(config,'vmic','-1','NICOLE.input','Nodes 2')
    if nodesvmic2 == '-1': nodesvmic2=get_value(config,'v_mic','-1','NICOLE.input','Nodes 2')
    nodesvmicx2=''
    nod=nodesvmic2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesvmic2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesvmicx2=locnod
    nodesblong2=get_value(config,'B_long','-1','NICOLE.input','Nodes 2')
    if nodesblong2 == '-1': nodesblong2=get_value(config,'Blong','-1','NICOLE.input','Nodes 2')
    if nodesblong2 == '-1': nodesblong2=get_value(config,'Bl','-1','NICOLE.input','Nodes 2')
    if nodesblong2 == '-1': nodesblong2=get_value(config,'B_z','-1','NICOLE.input','Nodes 2')
    if nodesblong2 == '-1': nodesblong2=get_value(config,'Bz','-1','NICOLE.input','Nodes 2')
    nodesblongx2=''
    nod=nodesblong2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesblong2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesblongx2=locnod
    nodesbx2=get_value(config,'B_x','-1','NICOLE.input','Nodes 2')
    if nodesbx2 == '-1': nodesbx2=get_value(config,'Bx','-1','NICOLE.input','Nodes 2')
    nodesbxx2=''
    nod=nodesbx2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesbx2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesbxx2=locnod
    nodesby2=get_value(config,'B_y','-1','NICOLE.input','Nodes 2')
    if nodesby2 == '-1': nodesby2=get_value(config,'By','-1','NICOLE.input','Nodes 2')
    nodesbyx2=''
    nod=nodesby2
    if re.search(',',nod): # Node locations
        locnod=[float(x) for x in nod.split(',')]
        if (len(locnod) > 99): sys.exit(1)
        locnod.sort()
        nodesby2=str(len(locnod))
        for i in range(99-len(locnod)): locnod.append(0.)
        nodesbyx2=locnod
    nodesab2=get_value(config,'Abundances','-1','NICOLE.input','Nodes 2')
    if nodesab2 == '-1': nodesab2=get_value(config,'Number of abundances','-1','NICOLE.input','Nodes 2')
    f=open('nodelocations.dat'+suffix,'wb')
    scratch=list()
    if (len(nodesTx) > 0):
        scratch.append(nodesT)
        scratch.extend(nodesTx)
    else:
        scratch.append(-1)
    if (len(nodesvx) > 0):
        scratch.append(nodesv)
        scratch.extend(nodesvx)
    else:
        scratch.append(-1)
    if (len(nodesvmicx) > 0):
        scratch.append(nodesvmic)
        scratch.extend(nodesvmicx)
    else:
        scratch.append(-1)
    if (len(nodesblongx) > 0):
        scratch.append(nodesblong)
        scratch.extend(nodesblongx)
    else:
        scratch.append(-1)
    if (len(nodesbxx) > 0):
        scratch.append(nodesbx)
        scratch.extend(nodesbxx)
    else:
        scratch.append(-1)
    if (len(nodesbyx) > 0):
        scratch.append(nodesby)
        scratch.extend(nodesbyx)
    else:
        scratch.append(-1)
    if (len(nodesTx2) > 0):
        scratch.append(nodesT2)
        scratch.extend(nodesTx2)
    else:
        scratch.append(-1)
    if (len(nodesvx2) > 0):
        scratch.append(nodesv2)
        scratch.extend(nodesvx2)
    else:
        scratch.append(-1)
    if (len(nodesvmicx2) > 0):
        scratch.append(nodesvmic2)
        scratch.extend(nodesvmicx2)
    else:
        scratch.append(-1)
    if (len(nodesblongx2) > 0):
        scratch.append(nodesblong2)
        scratch.extend(nodesblongx2)
    else:
        scratch.append(-1)
    if (len(nodesbxx2) > 0):
        scratch.append(nodesbx2)
        scratch.extend(nodesbxx2)
    else:
        scratch.append(-1)
    if (len(nodesbyx2) > 0):
        scratch.append(nodesby2)
        scratch.extend(nodesbyx2)
    else:
        scratch.append(-1)
    for data in scratch: f.write(struct.pack('<'+flf,float(data)))
    f.close()
    #
    # Read input data files and create native NICOLE data files
    #
    if mode == 'i': # In inversion mode
        [filemodeprof,nxout,nyout,nlout]=check_prof(obsprof)
        if filemodeprof == 'inexistent': 
            if icycle == 0: 
                sys.exit(1)
            else:
                useorigfileprof=1
        else:
            nxprof=nxout ; nyprof=nyout ; nlam2=nlout
        nx=nxprof
        ny=nyprof
        if nlam2 != nlam:
            print 'Error. Wavelength grid in NICOLE.input has ',nlam,' wavelengths.'
            print 'But there are ',nlam2,' in the observed profile file:',obsprof
            sys.exit(1)
        useorigfileprof=0
        if filemodeprof == 'nicole2.3' or icycle >= 1:
            useorigfileprof=1
        else:
            f=open('__inputprof.bin'+suffix,'wb')
            f.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.3bp     ',nxprof,nyprof,nlam)) # First record
            for i in range(nlam*4-16/8-1-1): f.write(struct.pack('<'+flf,0.)) # Fill record
            percent=-1
            seq=0
            for ix in range(nxprof):
                for iy in range(nyprof):
                    data=read_prof(obsprof, filemodeprof, nx, ny, nlam, ix, iy, sequential=seq)
                    seq=1
                    for d in data:
                        d=d/float(contval) # Normalization value
                        f.write(struct.pack('<'+flf,d))
                    if ((ix*ny+iy+1)*100./nx/ny > percent):
                        percent=int((ix*ny+iy+1)*1./nx/ny*100)
                        sys.stdout.write('\r'+'Preparing file with observed profiles...'+str(percent)+'%')
                        sys.stdout.flush()
            f.close() # Done with input profiles
        print ''
        [filemodemod,nxout,nyout,nzout]=check_model(inputmodel)
        useorigfilemod=0
        if filemodemod == 'inexistent': 
            if icycle == 0: 
                sys.exit(1)
            else:
                useorigfilemod=1
        else:
            nxmod=nxout ; nymod=nyout ; nz=nzout
        if inputmodel2 !='':
            [filemodemod2,nxmod2,nymod2,nz2]=check_model(inputmodel2)
            if (filemodemod2 != 'nicole2.6' or nxmod2 != nxmod or nymod2 != nymod
                or nz2 != nz): 
                print 'Error. Model 2 is not in native binary format or is not'
                print 'compatible with model 1'
                print 'Offending file:'+inputmodel2
                sys.exit(1)
        padding=0
        if (nxprof != nxmod or nyprof != nymod) and useorigfilemod == 0:
            print 'Warning! Observed profile file:'+obsprof
            print 'has dimensions:',nxprof,'x',nyprof,'=',nxprof*nyprof,'pixels'
            print 'The input model has:',nxmod,'x',nymod,'=',nxmod*nymod,'pixels'
            print 'Adopting geometry in the profile file'
            if nxmod*nymod > nxprof*nyprof:
                print 'The last models in the input model file will be ignored'
            if nxprof*nyprof > nxmod*nymod:
                print 'The input model will be padded by repeating the last model'
            padding=1
        if padding == 0 and (filemodemod == 'nicole2.6' or icycle >= 1):
            useorigfilemod=1
        else:
            useorigfilemod=0
            f=open('__inputmodel.bin'+suffix,'wb')
            f.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.6bm     ',nxprof,nyprof,nz)) # First record
            for i in range(22*nz+3+8+92-16/8-1-1): f.write(struct.pack('<'+flf,0.)) # Fill record
            percent=-1
            seq=0
            for ix in range(nxprof):
                for iy in range (nyprof):
                    ix2=ix
                    iy2=iy
                    if ix2 > nxmod-1: 
                        ix2=nxmod-1
                        seq=0
                    if iy2 > nymod-1: 
                        iy2=nymod-1
                        seq=0
                    try:
                        data=read_model(inputmodel, filemodemod, nx, ny, nz, ix2, iy2, sequential=seq)
                    except:
                        data=[x*0. for x in range(22*nz+3+8+92)]
                        print 'error in file ',inputmodel
                    if filemodemod != 'nicole2.6': # Old format, need to add abundances
                        for ab in numab:
                            data.append(float(ab))
                        if len(abundances) != 92:
                            print 'Number of elements in the abundance set != 92'
                            sys.exit(2)
                    seq=1

                    for d in data: f.write(struct.pack('<'+flf,d))
                    if ((ix*ny+iy+1)*100./nx/ny > percent):
                        percent=int((ix*ny+iy+1)*1./nx/ny*100)
                        sys.stdout.write('\r'+'Preparing file with input model...'+str(percent)+'%')
                        sys.stdout.flush()
            f.close() # Done with input model
        print ''
        if stray != '':
            [filemodestray,nxout,nyout,nlout]=check_prof(stray)
            useorigfilestray=0
            if filemodestray == 'inexistent': 
                if icycle == 0: 
                    sys.exit(1)
                else:
                    useorigfilestray=1
            else:
                nxstray=nxout ; nystray=nyout ; nlam2=nlout
            if nlam2 != nlam:
                print 'Error. Wavelength grid in NICOLE.input has ',nlam,' wavelengths.'
                print 'But there are ',nlam2,' in the stray light profile file:',stray
                sys.exit(1)
            if useorigfilestray == 0 and (nxprof != nxstray or nyprof != nystray):
                print 'Error! Stray light profile file:'+stray
                print 'has ',nxstray,'x',nystray,'=',nxstray*nystray,' profiles.'
                print 'The input profile has ',nxprof*nyprof
                sys.exit(1)
            if filemodestray == 'nicole2.3' or icycle >= 1:
                useorigfilestray=1
            else:
                f=open('__strayprof.bin'+suffix,'wb')
                f.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.3bp     ',nxprof,nyprof,nlam)) # First record
                for i in range(nlam*4-16/8-1-1): f.write(struct.pack('<'+flf,0.)) # Fill record
                percent=-1
                seq=0
                for ix in range(nxprof):
                    for iy in range(nyprof):
                        ix2=ix
                        iy2=iy
                        if ix2 > nxstray-1: 
                            ix2=nxstray-1
                            seq=0
                        if iy2 > nystray-1: 
                            iy2=nystray-1
                            seq=0
                        data=read_prof(stray, filemodestray, nxstray, nystray, nlam, ix2, iy2, sequential=seq)
                        seq=1
                        for d in data:
                            d=d/float(contval) # Normalization value
                            f.write(struct.pack('<'+flf,d))
                        if ((ix*ny+iy+1)*100./nx/ny > percent):
                            percent=int((ix*ny+iy+1)*1./nx/ny*100)
                            sys.stdout.write('\r'+'Preparing file with stray light...'+str(percent)+'%')
                            sys.stdout.flush()
                f.close() # Done with stray light profiles
            print ''
    elif mode == 's' or mode == 'c': # In Synthesis mode
        [filemodemod,nxout,nyout,nzout]=check_model(inputmodel)
        if filemodemod == 'inexistent': 
            if icycle == 0: 
                sys.exit(1)
            else:
                useorigfilemod=1
        else:
            nxmod=nxout ; nymod=nyout ; nz=nzout
        nx=nxmod
        ny=nymod
        useorigfileprof=0
        if filemodemod == 'nicole2.6' or icycle >= 1:
            useorigfilemod=1
        else:
            useorigfilemod=0
            f=open('__inputmodel.bin'+suffix,'wb')
            f.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.6bm     ',nxmod,nymod,nz)) # First record
            for i in range(22*nz+92+3+8-16/8-1-1): f.write(struct.pack('<'+flf,0.)) # Fill record
            percent=-1
            seq=0
            for ix in range(nxmod):
                for iy in range(nymod):
                    data=read_model(inputmodel, filemodemod, nx, ny, nz, ix, iy, sequential=seq)
                    if filemodemod != 'nicole2.6': # Old format, need to add abundances
                        for ab in numab:
                            data.append(ab)
                        if len(abundances) != 92:
                            print 'Number of elements in the abundance set != 92'
                            sys.exit(2)
                    for d in data: f.write(struct.pack('<'+flf,d))
                    if (int((ix*ny+(iy+1))*100./nx/ny) > percent):
                        percent=int((ix*ny+(iy+1))*100./nx/ny)
                        sys.stdout.write('\r'+'Preparing input model...'+str(percent)+'%')
                        sys.stdout.flush()
            f.close() # Done with input model
        print ''
        if stray != '':
            [filemodestray,nxout,nyout,nlout]=check_prof(stray)
            useorigfilestray=0
            if filemodestray == 'inexistent': 
                if icycle == 0: 
                    sys.exit(1)
                else:
                    useorigfilestray=1
            else:
                nxstray=nxout ; nystray=nyout ; nlam2=nlout
            if nlam2 != nlam:
                print 'Error. Wavelength grid in NICOLE.input has ',nlam,' wavelengths.'
                print 'But there are ',nlam2,' in the stray light profile file:',stray
                sys.exit(1)
            if useorigfilestray == 0 and (nxmod != nxstray or nymod != nystray):
                print 'Warning! Stray light profile file:'+stray
                print 'has ',nxstray,'x',nystray,'=',nxstray*nystray,' profiles. The input model file has ',nxmod,'x',nymod,'=',nxmod*nymod
                print 'Adopting the model geometry'
                if nxstray*nystray > nxmod*nymod:
                    print 'The last profiles in the stray light file will be ignored'
                if nxmod*nymod > nxstray*nystray:
                    print 'The stray light file will be padded by repeating the last profile'
            if filemodestray == 'nicole2.3' or icycle >= 1:
                useorigfilestray=1
            else:
                f=open('__strayprof.bin'+suffix,'wb')
                f.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.3bp     ',nxmod,nymod,nlam)) # First record
                for i in range(nlam*4-16/8-1-1): f.write(struct.pack('<'+flf,0.)) # Fill record
                percent=-1
                for ix in range(nxmod):
                    for iy in range(nymod):
                        ix2=ix
                        iy2=iy
                        if ix2 > nxstray-1: 
                            ix2=nxstray-1
                            seq=0
                        if iy2 > nystray-1: 
                            iy2=nystray-1
                            seq=0
                        data=read_prof(stray, filemodestray, nx, ny, nlam, ix2, iy2)
                        for d in data:
                            d=d/float(contval) # Normalization value
                            f.write(struct.pack('<'+flf,d))
                        if (int((ix*ny+(iy+1))*100./nx/ny) > percent):
                            percent=int((ix*ny+(iy+1))*100./nx/ny)
                            sys.stdout.write('\r'+'Preparing file with stray light profile...'+str(percent)+'%')
                            sys.stdout.flush()
                f.close() # Done with stray light profiles
            print ''
#
# Write __input.dat file for NICOLE
    f=open('__input.dat'+suffix,'w')
    f.write(str(ncycles)+'        ! ncycles \n') # debug
    f.write(str(cycle0)+'        ! cycle0 \n') # debug
    f.write(str(nx)+'        ! nx \n')
    f.write(str(ny)+'        ! ny \n')
    f.write(str(nz)+'        ! nz \n')
    f.write(pixx1+'   '+pixy1+'      ! pixx1, pixy1 \n')
    f.write(pixx2+'   '+pixy2+'      ! pixx2, pixy2 \n')
    f.write(irec0+'          ! irec0 \n')
    f.write(mode+'           ! mode. Below: modelin,modelin2, profin, profout, modelout,modelout2 \n')
    if useorigfilemod == 1:
        f.write(inputmodel+'\n')
    else:
        f.write('__inputmodel.bin'+suffix+'  \n')
    f.write(inputmodel2+'\n')
    if useorigfileprof == 1:
        f.write(obsprof+'\n')
    else:
        f.write('__inputprof.bin'+suffix+'  \n')
    f.write(outprof+'\n')
    f.write(outputmodel+'\n')
    f.write(outputmodel2+'\n')
    f.write(formal+' '+boundarycond+'      ! formal boundarycond \n')
    f.write(helio+' '+'   ! helio \n')
    f.write(printout+'    ! printout \n')
    f.write(maxiters+'    ! maxiters. Below: stray profile \n')
    if stray != '':
        if useorigfilestray == 1:
            f.write(stray+'\n')
        else:
            f.write('__strayprof.bin'+suffix+'  \n')
    else:
        f.write('\n')
    f.write(noise+'    ! noise \n')
    f.write(maxinv+'    ! maxinv \n')
    f.write(acceptchisq+'     ! acceptchisq\n')
    f.write(speed+'     ! speed \n')
    f.write(always_compute_der+ ' '+cent_der+'     ! always_compute_der, centered_derivatives \n')
    f.write(gravity+'     ! gravity \n')
    f.write(regul+'     ! regularization \n')
    f.write(update_opac+'    ! update_opac \n')
    f.write(negligible_opac+'    ! negligible_opac \n')
    f.write(contref+'       !  contref \n')
    f.write(inputdens+'       ! inputdens \n')
    f.write(keep_el_p+' '+keep_gas_p+' '+keep_rho+' '+keep_nH+' '+ 
            keep_nHminus+' '+keep_nHplus+' '+keep_nH2+' '+keep_nh2plus+
            '       ! Keep: el_p, gas_p, rho, nH, nH-, nH+, nH2, nH2+ \n')
    f.write(hscale+'          ! hscale \n')
    f.write(opacities+'  '+opacitiesUV+'  ! opac, opacUV \n')
    l=len(elneglectopacidx)
    f.write(str(l)+' ! elements to ignore in opac\n')
    if l > 0:
        for el in elneglectopacidx: f.write(str(el)+'\n')
    f.write(eqstate+'  '+eqstateH+' '+peconsistency+'          ! Eq state, Eq state for H, Pe_consistency \n')
    f.write(sethydro+'     '+setnH+'   '+restart+'   ! sethydro, setnH, restart \n')
    f.write(depcoef+' '+write_depcoef+'   ! depcoef mode, write \n')
    f.write(nlines+'   ! nlines \n')
    f.write(nregions+'    ! nregions \n')
    for iregion in range(int(nregions)):
        tmp1=get_value(config,'First wavelength','','NICOLE.input',
                       'Region '+str(iregion+1))
        tmp2=get_value(config,'Wavelength step','','NICOLE.input'
                      ,'Region '+str(iregion+1))
        if re.search('mA',tmp2): # Convert mA to A
            tmp2=re.sub('[^0-9.-]','',tmp2) # Remove text
            step=float(tmp2)*1e-3
            tmp2=str(step)
        tmp2=re.sub('[^0-9.-]','',tmp2) # Remove text
        tmp3=get_value(config,'Number of wavelengths','','NICOLE.input'
                      ,'Region '+str(iregion+1))
        tmp4=get_value(config,'Macroturbulence enhancement','1.0','NICOLE.input'
                       ,'Region '+str(iregion+1))
        tmp5=get_value(config,'Layer','1','NICOLE.input'
                      ,'Region '+str(iregion+1))
        tmp6=get_value(config,'Opacity enhancement','1.0','NICOLE.input'
                      ,'Region '+str(iregion+1))
        tmp7=get_value(config,'Observations additive constant','0.0','NICOLE.input'
                      ,'Region '+str(iregion+1))
        tmp8=get_value(config,'Observations multiplicative constant','1.0','NICOLE.input'
                      ,'Region '+str(iregion+1))
        tmp9=get_value(config,'Gaussian profile sigma','0.0','NICOLE.input'
                      ,'Region '+str(iregion+1))
        f.write(tmp1+' '+tmp2+' '+tmp3+' '+tmp4+' '+tmp5+' '+tmp6+' '+tmp7+' '+tmp8+' '+tmp9+'\n')
    # Parse spectral line database file
    fl=open('LINES')
    LINES_lines=fl.readlines()
    fl.close()
    LINES_lines_lower=['']
    for l in LINES_lines:
        LINES_lines_lower.append(lower_to_sep(l)) # Convert keys to lowercase
    config_lines=ConfigObj(LINES_lines_lower)
    for iline in range(int(nlines)):
        line_label=get_value(config,'Line','','NICOLE.input','Line '+str(iline+1))
        # Line data in LINES file database
        elem=get_value(config_lines,'Element','','LINES',line_label)
        ion_stage=get_value(config_lines,'Ionization stage','','LINES',line_label)
        if int(ion_stage) < 1 or int(ion_stage) > 3:
            print 'Ionization stage '+ion_stage+' not supported'
            print 'File LINES, section:',line_label
            sys.exit()
        wlength=get_value(config_lines,'Wavelength','','LINES',line_label)
        energy_low=get_value(config_lines,'Excitation potential','','LINES',line_label)
        if re.search('(?i)cm',energy_low): # Convert cm-1 to eV
            energy_low=re.sub('[^0-9.]','',energy_low) # Remove text
            ev=float(energy_low)*0.000123980262342235
            energy_low=str(ev)
        energy_low=re.sub('[^0-9.]','',energy_low) # Remove text
        loggf=get_value(config_lines,'Log(gf)','','LINES',line_label)
        string=get_value(config_lines,'Term (lower)','','LINES',line_label)
        mult_low=string[0]
        des_low=string[1]
        j_low=string[2:]
        if re.search('[^0-9]',mult_low):
            print 'Wrong multiplicity in file LINES, section:',line_label
            print 'Found:'+mult_low+' in string:'+string
            print 'Multiplicity must be 2s+1 (s is the spin moment). The string must be of the'
            print 'form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^0-9.]',j_low):
            print 'Wrong total angular momentum j in file LINES, section:',line_label
            print 'Found:'+j_low+' in string:'+string
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^SPDFGHIJKLMNOQpfhkmortuvw]',des_low):
            print 'Wrong angular momentum L in file LINES, section:',line_label
            print 'Found:'+des_low+' in string:'+string
            print 'Must be one of S, P, D, F, G, H, I... (p, f, h, k for semi-integers)'
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        string=get_value(config_lines,'Term (upper)','','LINES',line_label)
        mult_up=string[0]
        des_up=string[1]
        j_up=string[2:]
        if re.search('[^0-9]',mult_up):
            print 'Wrong multiplicity in file LINES, section:',line_label
            print 'Found:'+mult_up+' in string:'+string
            print 'Multiplicity must be 2s+1 (s is the spin moment). The string must be of the'
            print 'form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^0-9.]',j_up):
            print 'Wrong total angular momentum j in file LINES, section:',line_label
            print 'Found:'+j_up+' in string:'+string
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^SPDFGHIJKLMNOQpfhkmortuvw]',des_up):
            print 'Wrong angular momentum L in file LINES, section:',line_label
            print 'Found:'+des_up+' in string:'+string
            print 'Must be one of S, P, D, F, G, H, I... (p, f, h, k for semi-integers)'
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        collisions=get_value(config_lines,'Collisions','1','LINES',line_label)
        collisions=collisions.lower()
        if collisions[0] == 'u': collisions='1'
        if collisions[0] == 'b': collisions='2'
        if collisions[0] == 'g': collisions='3'
        bark_sigma=get_value(config_lines,'Damping sigma','-1','LINES',line_label)
        bark_alpha=get_value(config_lines,'Damping alpha','-1','LINES',line_label)
        g_rad=get_value(config_lines,'Gamma Radiative','-1','LINES',line_label)
        g_strk_12=get_value(config_lines,'Gamma Stark','-1','LINES',line_label)
        g_vdW_16=get_value(config_lines,'Gamma van der Waals','-1','LINES',line_label)
        if collisions != '1' and collisions != '2' and collisions != '3':
            print 'Collisions must be either Unsold, Barklem et al or Gamma'
            print 'Or alternatively a number between 1 and 3'
            print 'File LINES, section:'+line_label
            sys.exit(1)
        if collisions == 1:
            bark_sigma=-1
            bark_alpha=-1
        if collisions == 3:
            g_rad=-1
            g_strk_12=-1
            g_vdW_16=-1
        vdw_enh=get_value(config_lines,'Damping enhancement','1.0','LINES',line_label)
        width=get_value(config_lines,'Width','2.0','LINES',line_label)
        linemode=get_value(config_lines,'Mode','lte','LINES',line_label)
        linemode=linemode.lower()
        transition='0' ; nlratio='1.0' ; nuratio='1.0'
        if linemode != 'lte' and linemode != 'nlte':
            print 'Mode must be either LTE or NLTE'
            print 'in file LINES, section:',line_label
            sys.exit(1)
        if linemode == 'nlte':
            transition=get_value(config_lines,'Transition index in model atom','','LINES',line_label)
            nlratio=get_value(config_lines,'Lower level population ratio','1.0','LINES',line_label)
            nuratio=get_value(config_lines,'Upper level population ratio','1.0','LINES',line_label)
        hf=get_value(config_lines,'Hyperfine structure','n','LINES',line_label)
        hf=hf.lower()
        if hf == 'n' or hf == 'no' or hf == 'false' or hf == '.false.' or hf == 'f' or hf == '.f.': hf='F'
        if hf == 'y' or hf == 'yes' or hf == 'true' or hf == '.true.' or hf == 't' or hf == '.t.': hf='T'
        alow='';blow='';aup='';bup='';spini=''
        if hf == 'T':
            alow=get_value(config_lines,'Hyperfine Alow','','LINES',line_label)
            blow=get_value(config_lines,'Hyperfine Blow','','LINES',line_label)
            aup=get_value(config_lines,'Hyperfine Aup','','LINES',line_label)
            bup=get_value(config_lines,'Hyperfine Bup','','LINES',line_label)
            spini=get_value(config_lines,'Nuclear spin','','LINES',line_label)
        else:
            alow='0' ; blow='0' ; aup='0' ; bup='0' ; spini='0'
        extra_vmic=get_value(config_lines,'Extra v_mic','0','LINES',line_label)
        if extra_vmic == '0':
            extra_vmic=get_value(config_lines,'Extra micro','0','LINES',line_label)
        if extra_vmic == '0':
            extra_vmic=get_value(config_lines,'Extra microturbulence','0','LINES',line_label)
        # Override with values from NICOLE.input?
        elem=get_value(config,'Element',elem,'NICOLE.input','Line '+str(iline+1))
        ion_stage=get_value(config,'Ionization stage',ion_stage,'NICOLE.input','Line '+str(iline+1))
        if int(ion_stage) < 1 or int(ion_stage) > 3:
            print 'Ionization stage '+ion_stage+' not supported'
            print 'File NICOLE.input, section:','Line '+str(iline+1)
            sys.exit()
        wlength=get_value(config,'Wavelength',wlength,'NICOLE.input','Line '+str(iline+1))
        energy_low=get_value(config,'Excitation potential',energy_low,'NICOLE.input','Line '+str(iline+1))
        if re.search('(?i)cm',energy_low): # Convert cm-1 to eV
            energy_low=re.sub('[^0-9.]','',energy_low) # Remove text
            ev=float(energy_low)*0.000123980262342235
            energy_low=str(ev)
        energy_low=re.sub('[^0-9.]','',energy_low) # Remove text
        loggf=get_value(config,'Log(gf)',loggf,'NICOLE.input','Line '+str(iline+1))
        string=get_value(config,'Term (lower)',mult_low+des_low+j_low,'NICOLE.input','Line '+str(iline+1))
        mult_low=string[0]
        des_low=string[1]
        j_low=string[2:]
        if re.search('[^0-9]',mult_low):
            print 'Wrong multiplicity in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+mult_low+' in string:'+string
            print 'Multiplicity must be 2s+1 (s is the spin moment). The string must be of the'
            print 'form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^0-9.]',j_low):
            print 'Wrong total angular momentum j in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+j_low+' in string:'+string
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^SPDFGHIJKLMNOQpfhkmortuvw]',des_low):
            print 'Wrong angular momentum L in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+des_low+' in string:'+string
            print 'Must be one of S, P, D, F, G, H, I... (p, f, h, k for semi-integers)'
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        string=get_value(config,'Term (upper)',mult_up+des_up+j_up,'NICOLE.input','Line '+str(iline+1))
        mult_up=string[0]
        des_up=string[1]
        j_up=string[2:]
        if re.search('[^0-9]',mult_up):
            print 'Wrong multiplicity in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+mult_up+' in string:'+string
            print 'Multiplicity must be 2s+1 (s is the spin moment). The string must be of the'
            print 'form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^0-9.]',j_up):
            print 'Wrong total angular momentum j in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+j_up+' in string:'+string
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        if re.search('[^SPDFGHIJKLMNOQpfhkmortuvw]',des_up):
            print 'Wrong angular momentum L in file NICOLE.input, section:','Line '+str(iline+1)
            print 'Found:'+des_up+' in string:'+string
            print 'Must be one of S, P, D, F, G, H, I... (p, f, h, k for semi-integers)'
            print 'The string must be of the form 2s+1Lj. For example: 3S0.5'
            sys.exit(1)
        collisions=get_value(config,'Collisions',collisions,'NICOLE.input','Line '+str(iline+1))
        collisions=collisions.lower()
        if collisions[0] == 'u': collisions='1'
        if collisions[0] == 'b': collisions='2'
        if collisions[0] == 'g': collisions='3'
        bark_sigma=get_value(config,'Damping sigma',bark_sigma,'NICOLE.input','Line '+str(iline+1))
        bark_alpha=get_value(config,'Damping alpha',bark_alpha,'NICOLE.input','Line '+str(iline+1))
        g_rad=get_value(config,'Gamma Radiative',g_rad,'NICOLE.input','Line '+str(iline+1))
        g_strk_12=get_value(config,'Gamma Stark',g_strk_12,'NICOLE.input','Line '+str(iline+1))
        g_vdW_16=get_value(config,'Gamma van der Waals',g_vdW_16,'NICOLE.input','Line '+str(iline+1))
        if collisions != '1' and collisions != '2' and collisions != '3':
            print 'Collisions must be either Unsold, Barklem et al or Gamma'
            print 'Or alternatively a number between 1 and 3'
            print 'File LINES, section:'+line_label
            sys.exit(1)
        if collisions == 1:
            bark_sigma=-1
            bark_alpha=-1
        if collisions == 3:
            g_rad=-1
            g_strk_12=-1
            g_vdW_16=-1
        vdw_enh=get_value(config,'Damping enhancement',vdw_enh,'NICOLE.input','Line '+str(iline+1))
        width=get_value(config,'Width',width,'NICOLE.input','Line '+str(iline+1))
        linemode=get_value(config,'Mode',linemode,'NICOLE.input','Line '+str(iline+1))
        linemode=linemode.lower()
        if linemode != 'lte' and linemode != 'nlte':
            print 'Mode must be either LTE or NLTE'
            print 'in file NICOLE.input, section:','Line '+str(iline+1)
            sys.exit(1)
        if linemode == 'lte':
            transition='0' ; nlratio='1.0' ; nuratio='1.0'
        if linemode == 'nlte':
            transition=get_value(config,'Transition index in model atom',transition,'NICOLE.input','Line '+str(iline+1))
            nlratio=get_value(config,'Lower level population ratio',nlratio,'NICOLE.input','Line '+str(iline+1))
            nuratio=get_value(config,'Upper level population ratio',nuratio,'NICOLE.input','Line '+str(iline+1))
        hf=get_value(config,'Hyperfine structure',hf,'NICOLE.input','Line '+str(iline+1))
        hf=hf.lower()
        if hf == 'n' or hf == 'no' or hf == 'false' or hf == '.false.' or hf == 'f' or hf == '.f.': hf='F'
        if hf == 'y' or hf == 'yes' or hf == 'true' or hf == '.true.' or hf == 't' or hf == '.t.': hf='T'
        if hf == 'T':
            alow=get_value(config,'Hyperfine Alow',alow,'NICOLE.input','Line '+str(iline+1))
            blow=get_value(config,'Hyperfine Blow',blow,'NICOLE.input','Line '+str(iline+1))
            aup=get_value(config,'Hyperfine Aup',aup,'NICOLE.input','Line '+str(iline+1))
            bup=get_value(config,'Hyperfine Bup',bup,'NICOLE.input','Line '+str(iline+1))
            spini=get_value(config,'Nuclear spin',spini,'NICOLE.input','Line '+str(iline+1))
        else:
            alow='0' ; blow='0' ; aup='0' ; bup='0' ; spini='0'


        f.write(width+' '+elem.upper()+' '+ion_stage+' '+wlength+' '+vdw_enh+' '+energy_low+
                ' '+loggf+' '+des_low+' '+mult_low+' '+j_low+' '+des_up+' '+ 
                mult_up+' '+j_up+' '+collisions+' '+bark_sigma+' '+bark_alpha+
                ' '+g_rad+' '+g_strk_12+' '+g_vdW_16+' '+
                ' '+transition+' '+nlratio+' '+nuratio+' '+hf+' '+alow+' '+
                blow+' '+aup+' '+bup+' '+spini+' '+extra_vmic+'\n')
    # Abundances
    for element in elements:
        f.write(str(abundances[element])+' ')
    f.write('\n')
    # Nodes
    f.write(nodesT+' ! Nodes in T \n')
    f.write(nodesv+' ! Nodes in v_los \n')
    f.write(nodesvmic+' ! Nodes in v_mic \n')
    f.write(nodesblong+' ! Nodes in B_long \n')
    f.write(nodesbx+' ! Nodes in B_x \n')
    f.write(nodesby+' ! Nodes in B_y \n')
    f.write(nodesstray+' ! Nodes in Stray light \n')
    f.write(nodesvmac+' ! Nodes in v_mac \n')
    f.write(nodesexp+' ! Nodes in ffactor \n')
    nab=0
    if nodesab != '-1' and nodesab != '0': 
        nodesab=nodesab.lower()
        nodesab=nodesab.split(" ")
        nab=len(nodesab)
        for iab in range(nab): # In case user lists element convert to Z
            if (nodesab[iab] in abundances):
                nodesab[iab]=str(elements.index(nodesab[iab])+1)
    f.write(str(nab)+' ! Number of abundances to invert \n')
    if nab >= 1: f.write('  '.join(nodesab)+'  ! Abundances to invert \n')
    f.write(nodesT2+' ! Nodes in T 2\n')
    f.write(nodesv2+' ! Nodes in v_los 2\n')
    f.write(nodesvmic2+' ! Nodes in v_mic 2\n')
    f.write(nodesblong2+' ! Nodes in B_long 2\n')
    f.write(nodesbx2+' ! Nodes in B_x 2\n')
    f.write(nodesby2+' ! Nodes in B_y 2\n')
    nab=0
    if nodesab2 != '-1': 
        nodesab2=nodesab2.lower()
        nodesab2=nodesab2.split(" ")
        nab=len(nodesab2)
        for iab in range(nab): # In case user lists element convert to Z
            if (nodesab2[iab] in abundances):
                nodesab2[iab]=str(elements.index(nodesab2[iab])+1)
    f.write(str(nab)+' ! Number of abundances to invert 2\n')
    if nab >= 1: f.write('  '.join(nodesab2)+'  ! Abundances to invert 2\n')

    outputpop=str(outputpop == '1')[0]
    outputcontop=str(outputcontop == '1')[0]
    outputNLTEsf=str(outputNLTEsf == '1')[0]
    f.write(debug+' '+interp+' '+outputpop+' '+outputcontop+' '+outputNLTEsf+
            ' ! debug, optimize grid, outputpop, outputcontop, outputNLTEsf')
    f.write('\n')
    # NLTE
    f.write(elim1+' Elim 1 \n')
    f.write(isum+' '+istart+' ISUM, ISTART \n')
    f.write(nmu+' '+usecolswitch+' '+nlteformalsolution+' NMU, usecolswitch, formalsolution \n')
    f.write(qnorm+' '+cper+' QNORM, CPER \n')
    f.write(velfree+' '+ngacc+' '+lambdaiters+' velfree, ngacc,lambdaiters \n')
    f.write(nlteoptthin+' '+nlteoptthick+' '+nltelinear+' nlte opt thin, opt thick, linear \n')
    f.write(nltemaxiters+' nltemaxiters \n')
    f.write(ltepop+' NLTE ltepop \n')
    f.close()

#
# Run NICOLE
#

#print ' READY\nYou can now run NICOLE'

print ''
print 'Starting code execution'
status=subprocess.call(nicolecommand.split())
try:
    import calccorrect
    calccorrect.message(correct= status == 0)
except:
    pass
