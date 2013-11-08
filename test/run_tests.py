#!/usr/bin/env python

# Script to run the test cases and check results

# To use in non-interactive mode (e.g. to test parallel mode on a
# queued job system), first launch manually the tests by going into
# each directory (syn1, syn2, ..., inv1, inv2, ..., etc) and launch
# each test with:
#   ./run_nicole.py --nicolecommand='xxxxx' Then run this
# script with the check-only option to analyze the results:
#   ./run_tests.py --check-only
#
# You can use --sequential to pause after launching NICOLE and before checking 
#  the results
# You can specify the command to pass to run_nicole.py with the 
#  --nicolecommand='xxxx' switch

# Get command line arguments
def get_args():
    import getopt
    import sys

# Get command-line arguments
    nicolecommand='../../main/nicole'
    interactive=1
    sequential=0
    clean=0
    try:
        opts, args = getopt.getopt(sys.argv[1:],"", ["nicolecommand=","check-only","sequential","clean"])
    except:
        print 'Command-line option not recognized'
        print 'Usage:'
        print "run_nicole.py [--check-only] [--nicolecommand='command'] [--sequential] [--clean]"
        sys.exit(2)
    for o, a in opts:
        if o == '--nicolecommand':
            nicolecommand=a
        if o == '--check-only':
            interactive=0
        if o == '--sequential':
            sequential=1
        if o == '--clean':
            clean=1
    return [interactive,nicolecommand,sequential,clean]

# Run one test
def test_nicole(dir,message,interactive=1,nicolecommand='../../main/nicole'):
    import time
    import os

    print '\n\n'
    print '*******************   '+dir+'   *********************'

    result=1
    print message
    try:
        os.remove('Chisq.dat_1')
    except:
        donothing=1
    try:
        import glob
        for filename in glob.glob('profiling*.txt') :
            os.remove( filename ) 
    except:
        donothing=1
    try:
        import glob
        for filename in glob.glob('core*.txt') :
            os.remove( filename ) 
    except:
        donothing=1
    if interactive == 1:
        print 'Starting run. Output will be kept in '+dir+'/log.txt'
        print '     (Starting at '+str(time.localtime()[3])+':'+ \
            str(time.localtime()[4])+':'+str(time.localtime()[5])+')'
        logfile=open('log.txt','w')
        pipe=subprocess.Popen(['./run_nicole.py','--nicolecommand='+nicolecommand],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err=pipe.communicate()
        out=out.splitlines()
        if ' DONE' in out:
            result=0
        if result == 0:
            print '    The run has completed normally'
        else:
            print '     ** ERROR!! NICOLE crashed during this test'

        print '     (Finished at '+str(time.localtime()[3])+':'+ \
            str(time.localtime()[4])+':'+str(time.localtime()[5])+')'
    else:
        result=0
    return result

# Main program
import sys
import os
import re
import struct
import subprocess
from model_prof_tools import *

# Command-line arguments
[interactive,nicolecommand,sequential,clean]=get_args()

[int4f,intf,flf]=check_types()
 
print '\n\n                  Testing NICOLE'

# Check if numpy is available
numpy=0
try:
    import numpy
    numpy=1
except:
    print '\nPython module numpy has not been found. You will not be able to read'
    print 'IDL save files with NICOLE on this machine'

# Run tests
cwd=os.getcwd()
success=1

# conv1
dir='conv1'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('modelout.model')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will convert electron to gas pressure in the HSRA model',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_model('modelout.model')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,1,95]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_model('modelout.model',mode,1,1,95,0,0)
        m_pg=data[95*3:95*3+95]
        data2=read_model('hsra_full.bin','nicole2.3',1,1,95,0,0)
        hsra_pg=data2[95*3:95*3+95]

#        print 'pres=',m_pg[80],hsra_pg[80] # debug

        dif1=abs(m_pg[80]-hsra_pg[80])/hsra_pg[80]
        dif2=abs(m_pg[60]-hsra_pg[60])/hsra_pg[60]
        if dif1 > 0.35 or dif2 > 0.35:
            print '    ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'

# conv2
dir='conv2'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('modelout.model')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will convert gas to electron pressure in the HSRA model',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_model('modelout.model')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,1,95]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_model('modelout.model',mode,1,1,95,0,0)
        m_elp=data[95*5:95*5+95]
        data2=read_model('hsra_full.bin','nicole2.3',1,1,95,0,0)
        hsra_elp=data2[95*5:95*5+95]

#        print 'pres 1=',m_elp[80],hsra_elp[80] # debug
#        print 'pres 2=',m_elp[60],hsra_elp[60] # debug

        dif1=abs(m_elp[80]-hsra_elp[80])/hsra_elp[80]
        dif2=abs(m_elp[60]-hsra_elp[60])/hsra_elp[60]
        if dif1 > 0.5 or dif2 > 0.5:
            print '    ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'

# syn1
dir='syn1'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('hsra_mag.pro')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will test a simple synthesis in LTE with a magnetic atmosphere',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_prof('hsra_mag.pro')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,1,200]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_prof('hsra_mag.pro',mode,1,1,200,0,0)
        i=data[::4]
        q=data[1::4]
        u=data[2::4]
        v=data[3::4]

#        print 'max(i,q,u,v)=',max(i),max(q),max(u),max(v) # debug
#        print 'min(i,q,u,v)=',min(i),min(q),min(u),min(v) # debug

        if abs(max(i)-1.000)>0.03 or abs(min(i)-0.230)>0.03 or \
                abs(max(q)-0.)>2e-3 or abs(min(q)+.001)>2e-3 or \
                abs(max(u)-0.0036)>2e-3 or abs(min(u))>2e-3 or \
                abs(max(v)-0.029)>2e-3 or abs(min(v)+0.0295)>2e-3:
            print '     ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'


# syn2
dir='syn2'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('manyhsra_mag.pro')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will test many syntheses in LTE with a magnetic atmosphere',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_prof('manyhsra_mag.pro')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,10,200]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
#        data=read_prof('manyhsra_mag.pro',mode,1,1000,200,0,0)
        data=read_prof('manyhsra_mag.pro',mode,1,10,200,0,0)
        i=data[::4]
        q=data[1::4]
        u=data[2::4]
        v=data[3::4]
 
        limi=[0.99999, 0.219]
        limq=[0.0, -0.002]
        limu=[0.0068, -0.000274]
        limv=[0.0491, -0.0440]

#        print 'max(i,q,u,v)=',max(i),max(q),max(u),max(v) # debug
#        print 'min(i,q,u,v)=',min(i),min(q),min(u),min(v) # debug

        if abs(max(i)-limi[0])>0.01 or abs(min(i)-limi[1])>0.01 or \
                abs(max(q)-limq[0])>2e-3 or abs(min(q)-limq[1])>2e-3 or \
                abs(max(u)-limu[0])>2e-3 or abs(min(u)-limu[1])>2e-3 or \
                abs(max(v)-limv[0])>2e-3 or abs(min(v)-limv[1])>2e-3:
            print '     ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'

# syn3
dir='syn3'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('valc_mag.pro')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will test a simple synthesis in NLTE with a magnetic atmosphere',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_prof('valc_mag.pro')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,1,400]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_prof('valc_mag.pro',mode,1,1,400,0,0)
        data0=data # For use with syn4
        i=data[::4]
        q=data[1::4]
        u=data[2::4]
        v=data[3::4]
        i=i[200:]
        q=q[200:]
        u=u[200:]
        v=v[200:]


        limi=[ 0.9724 , 0.2539 ]
        limq=[ 1.8653e-05 , -1.03727e-06 ]
        limu=[ 2.48e-06 , -1.481088e-05 ]
        limv=[ 0.00270 , -0.0035475 ]

#        print 'max(i,q,u,v)=',max(i),max(q),max(u),max(v) # debug
#        print 'min(i,q,u,v)=',min(i),min(q),min(u),min(v) # debug

        if abs(max(i)-limi[0])>0.03 or abs(min(i)-limi[1])>0.03 or \
                abs(max(q)-limq[0])>2e-3 or abs(min(q)-limq[1])>2e-3 or \
                abs(max(u)-limu[0])>2e-3 or abs(min(u)-limu[1])>2e-3 or \
                abs(max(v)-limv[0])>2e-3 or abs(min(v)-limv[1])>2e-3:
            print '     ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'

# inv1
dir='inv1'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('inversion.pro')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will test a simple LTE inversion',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nlam]=check_model('inversion.mod')
    if [mode[0:6],nx,ny,nlam] != ['nicole',1,1,56]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_model('inversion.mod',mode, 1,1,56, 0,0)
        f=open('Chisq.dat_2')
        chisq=struct.unpack('<1'+flf,f.read(8))
        if chisq < 10 and \
                abs(data[56*2+26]-4718) > 250: # Temperature at ltau_500=-1.5
            print '     ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'
        
# inv2
dir='inv2'
os.chdir(cwd)
os.chdir(dir)

if interactive != 0 or clean == 1:
    try:
        os.remove('log.txt')
        os.remove('inversion.pro')
    except:
        do_nothing=1
if clean == 1:
    print '\n'+dir+': Removing output files...\n'
    result = 1
else:
    result=test_nicole(dir,'This will test two LTE inversions',interactive,nicolecommand)
if sequential == 1:
    print 'Press ENTER to continue'
    input=raw_input()
if result != 0:
    success=0
else:
    print '    Checking the results produced...'
    [mode,nx,ny,nz]=check_model('inversion.mod_2')
    if [mode[0:6],nx,ny,nz] != ['nicole',1,1,95]:
        print '     ** ERROR!! NICOLE produced an incorrect file'
        success=0
    else:
        data=read_model('inversion.mod_2',mode, 1,1,95, 0,0)
        f=open('Chisq.dat_2')
        chisq=struct.unpack('<1'+flf,f.read(8))
        if chisq < 10 and \
                abs(data[95*2+65]-4718) > 250: # Temperature at ltau_500=-1.5
            print '     ** ERROR!! NICOLE produced inaccurate results'
            success=0
        else:
            print '    Results appear to be correct'


# Finished. Print summary
if clean == 1:
    sys.exit(0)

if success == 1:
    print '\n\n         =========================================='
    print 'Congratulations! Your build of NICOLE has passed all tests successfully :)\n\n'
    if numpy == 0:
        print ' (Note, however, that you will not be able to use IDL save files on this machine)'
    sys.exit(0)
else:
    print '\n\n         =========================================='
    print 'Unfortunately your build of NICOLE has failed one or more tests :(\n\n'
    sys.exit(1)
