#!/usr/bin/env python

# Script to analyze current system and produce a suitable makefile.
#
# Instructions for USERS:
#
# Use the -h or --help options for usage instructions
#
#
# Instructions for DEVELOPERS:
#
# To adapt this script to a different code, follow these instructions:
#
# 1)You may want to have different command-line switches and arguments.
# Change these in the code that follows the getopt call. Remember to update
# the usage() function. Remember also to initialize the variables properly
# at the beginning of the main program
#
# 2)Edit the code-specific block. This is where you do all the preprocessing
# that is specific to your code, check certain compiler functionalities and,
# if necessary, edit the source code. Start by defining treetop and mains.
# Treetop points to the top of the source code tree. For example, if the tree 
# is under the current directory, treetop should be "./" (or even better
# '.'+dirsep, see dirsep below).
# Mains is a list of the files that contain main programs. These will be the
# targets built when a "make all" command is issued.
# Ignoredeps is a list of dependencies that is ok not to have in the source
# tree, so this script will not crash in the event that they are not found
# In doing all of this, use the variable dirsep to refer to the directory
# separator. This is usually the slash (/) character, except in Windows-like
# systems where it's the backslash (\). This script detects the operating
# system and sets the correct directory separation character in dirsep so
# use it to make your code more portable.
#
# 3)If necessary, change the scratch_dir variable to the name of a directory
# that can be used for temporary files. Make sure that this name does not 
# interfere with your code. 
#
# 4)Finally, check the makefile section to make any necessary tweaks required
# by your code
#

# Replicate UNIX which functionality in a cross-platform implementation
def which(program):
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

# Print usage information
def usage():
    print ""
    print "This program will automatically generate a makefile for your"
    print "NICOLE distribution. If invoked without arguments it will search"
    print "the default path in your system for one of the known F90 compilers."
    print ""
    print "Alternatively you may specify which compiler you want to use"
    print "but then you may also need to specify certain common switches"
    print "to control the compiler behavior."
    print ""
    print "Optional arguments:"
    print " -h, --help: Prints this help and exits"
    print " -v, --version: Prints version information and exits"
    print " -y: Answer yes to all questions"
    print " --mpi: Create makefile for MPI parallel version. If this option"
    print "   is not present, the makefile is for the serial version (nicole)."
    print " --recl: Manually specify record length for a 64-bit real in your "
    print "   compiler/architecture. If not present, detect automatically"
    print " --quiet: Produces less output (default)"
    print " --verbose: Produces more output"
    print " --showtree: Shows the full dependency tree of all source files"
    print " --ignore: Continue creating makefile even if some dependencies are not"
    print "   found in the source tree."
    print " --compiler=: Specify F90 compiler to use (may be preceded by the path"
    print "   if not in the search path). Use --compiler==list to print a list of"
    print "   known compilers. If your compiler is not one of these you will also"
    print "   need to specify the following options. Use double quotes to enclose a"
    print "   complex string. Examples: --compiler=ifort, --compiler=\"/usr/bin/ifort\" "
    print " --cswitch=: Command-line argument passed to the compiler to invoke it in "
    print "   compile-only (no linking) mode. Usually this is --cswitch=-c"
    print " --modpath=: Command-line argument passed to the compiler to specify the"
    print "   module search path. Usually this is --modpath=-M, --modpath=-I or --modpath=\"-module \""
    print "   (note the trailing space in this last example)."
    print " --modsuf=: Suffix of modules produced by the compiler. This is usually"
    print "   --modsuf=\".mod\", or in some older compilers --modsuf=\".M\""
    print " --otherflags=: Other flags that will be passed to the compiler, e.g. for"
    print "   optimization options, debugging, etc. Examples: --otherflags=-O3"
    print "   --otherflags=-g"
    print ""
    print ""

# Print program version
def printversion():
    print ""
    print "Version 1.0"
    print ""

# List known compilers
def listknowncompilers():
    print "Known Fortran90 compilers"
    print "gfortran"
    print "ifort"
    print "pgf90"
    print "f95"
    print "xlf90"

# Recursively find all files under a given directory
def walktree (top = ".", depthfirst = True):
    import os
    import stat

    names = os.listdir(top)
    if not depthfirst:
        yield top, names
    for name in names:
        try:
            st = os.lstat(os.path.join(top, name))
        except os.error:
            continue
        if stat.S_ISDIR(st.st_mode):
            for (newtop, children) in walktree (os.path.join(top, name), depthfirst):
                yield newtop, children
    if depthfirst:
        yield top, names

# Recursive function used to print out dependencies (used if showtree option 
# is selected)
def printsubtree(source_files,requires,ifile,irecurs):
    if irecurs > 20:
        print 'Stopping. Too many recursion levels in printsubtree'
        exit(1)
#    print 'd',requires[ifile]
#    print 'd',len(requires[ifile])
    for required in requires[ifile]:
        for i in range(irecurs):
            sys.stdout.write('....')
        print required
        irecurs2=irecurs+1
        ifile2=source_files.index(required)
        printsubtree(source_files,requires,ifile2,irecurs2)

# Main program

import sys
import os
import shutil
import getopt
import datetime
import re
import subprocess
import shlex

compiler=None
cswitch=None
modpath=None
modsuf=None
otherflags=''
ignore=0
quiet=1
recl=-1
yes=0
mpi=0
sopa=0
showtree=0

scratch_dir='000crmak' # Directory for temporary scratch files
shutil.rmtree(scratch_dir,ignore_errors=True)
cwd=os.getcwd()

fnull = open(os.devnull, 'w') # Null file

# Reuse previous flags?
sysarg=''
for arg in sys.argv[1:]:
    sysarg=sysarg+' "'+arg+'"'
if "--keepflags" in sys.argv: # Read previous flags from makefile
    f=open('makefile','r')
    lines=f.readlines()
    f.close()
    fline=''
    for line in lines:
        if re.search('# Flags:',line) != None:
            fline=line
    if fline == '':
        sysarg=''
    else:
        fline=re.sub('# Flags:','',fline)
        sysarg=re.sub('\n','',fline)
    print 'Using previous flags:',sysarg

# Look for compiler and flags in the command line
try:
    sysarglist=shlex.split(sysarg)
    opts, args = getopt.getopt(sysarglist, "hvy", ["help", "version", 
                                                    "ignore","quiet","verbose",
                                                     "mpi","showtree","sopa",
                                                     "keepflags","recl=",
                                                 "cswitch=", "compiler=", 
                                              "modpath=", "modsuf=",
                                                    "otherflags="]) 

except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)



for o,a in opts:
    if o == '-h' or o == '--help':
        usage()
        sys.exit(0)
    if o == '-v' or o == '--version':
        printversion()
        sys.exit(0)
    if o == '-y':
        yes=1
    if o == '--quiet':
        quiet=1
    if o == '--verbose':
        quiet=0
    if o == '--ignore':
        ignore=1
    if o == '--showtree':
        showtree=1
    if o == '--mpi':
        mpi=1
    if o == '--sopa':
        sopa=1
    if o == '--compiler':
        compiler=a
        if a == 'list':
            listknowncompilers()
            sys.exit(0)
    if o == '--recl':
        recl=a
    if o == '--cswitch':
        cswitch=a
    if o == '--modpath':
        modpath=a
    if o == '--modsuf':
        modsuf=a
    if o == '--otherflags':
        otherflags=a

# OS-dependent variables
if os.name == 'posix':
    rm='rm -f '
    cp='cp '
    dirsep='/'
elif os.name == 'mac':
    rm='rm -f '
    cp='cp '
    dirsep='/'
elif os.name == 'nt' or os.name == 'ce':
    rm='del /f '
    cp='copy '
    dirsep='\\'
else:
    rm='rm -f '
    cp='cp '
    dirsep='/'
    print 'Operating system not recognized:',os.name
    if (not ignore):
        print 'Use the --ignore switch to continue creating the makefile'
        print 'But the remove file commands in the makefile will likely be wrong'
        print 'You may need to edit it manually'
        sys.exit(1)

scratch_dir=scratch_dir+dirsep


if (compiler == None):
    # Search for F90 compiler
    print 'User did not specify a F90 compiler. Searching for available options...'
    stri='f95'
    if (which(stri) != None):
        compiler='f95'
        print 'Found '+stri
    else:
        print stri+' not found'
    stri='gfortran'
    if (which(stri) != None):
        compiler='gfortran'
        print 'Found '+stri
    else:
        print stri+' not found'
    stri='pgf90'
    if (which(stri) != None):
        compiler='pgf90'
        print 'Found '+stri
    else:
        print stri+' not found'
    stri='xlf90'
    if (which(stri) != None):
        compiler='xlf90'
        print 'Found '+stri
    else:
        print stri+' not found'
    stri='ifort'
    if (which(stri) != None):
        compiler='ifort'
        print 'Found '+stri
    else:
        print stri+' not found'
    stri='mpif90'
    if (which(stri) != None and mpi == 1):
        compiler='mpif90'
        print 'Found '+stri
    else:
        print stri+' not found'

    if (compiler == None):
        print 'Compiler auto-detect failed'
        sys.exit(1)

# Choose options for compiler
if (compiler == 'pgf90'):
    if cswitch == None:
        cswitch='-c'
    if modpath == None:
        modpath='-module '
    if modsuf == None:
        modsuf='.mod'
    if otherflags=='':
        otherflags=otherflags+' -O3 '
if (compiler == 'f95'):
    if cswitch == None:
        cswitch='-c'
    if modpath == None:
        modpath='-I'
    if modsuf == None:
        modsuf='.mod'
    if otherflags=='':
        otherflags=otherflags+' -O3 '
if (compiler == 'gfortran'):
    if cswitch == None:
        cswitch='-c'
    if modpath == None:
        modpath='-I'
    if modsuf == None:
        modsuf='.mod'
    if otherflags=='':
        otherflags=otherflags+' -O3 '
if (compiler == 'xlf90'):
    if cswitch == None:
        cswitch='-c'
    if modpath == None:
        modpath='-I'
    if modsuf == None:
        modsuf='.mod'
    if otherflags=='':
        otherflags=otherflags+' -O3 '
if (compiler == 'ifort'):
    if cswitch == None:
        cswitch='-c'
    if modpath == None:
        modpath='-I '
#        modpath='-module '
    if modsuf == None:
        modsuf='.mod'
    if otherflags=='':
        otherflags=otherflags+' -O3 '

if (compiler == None):
    print 'Error. You need to have a F90 compiler installed in your system'
    usage()
    sys.exit(1)

if cswitch == None: # Try to autodetect cswitch option
    os.mkdir(scratch_dir)
    f=open(scratch_dir+'test.f90','w')
    f.write('end\n')
    f.close()    
    os.chdir(scratch_dir)
    try:
        retcode=subprocess.call([compiler,'-c','test.f90'],stdout=fnull,stderr=fnull)
    except:
        retcode=10
    os.chdir(cwd)
    if retcode == 0 and os.path.exists(scratch_dir+'test.o'):
        cswitch='-c'
    else:
        print 'Failed to auto-detect cswitch flag for compiler:'+compiler
        print "Usually it is -c but this didn't seem to work"
        print 'Please, specify it manually by invoking this program with the --cswitch option'
    shutil.rmtree(scratch_dir,ignore_errors=True)
    print 'Autodetected cswitch:',cswitch
            
if modsuf == None: # Try to autodetect modsuf option
    shutil.rmtree(scratch_dir,ignore_errors=True)
    os.mkdir(scratch_dir)
    f=open(scratch_dir+'test.f90','w')
    f.write('module test\nend module test\n')
    f.close()    
    os.chdir(scratch_dir)
    try:
        retcode=subprocess.call([compiler,cswitch,'test.f90'],stdout=fnull,stderr=fnull)
    except:
        retcode=10
    os.chdir(cwd)
    testfor=['.mod','.M','.module'] # Possible options for modsuf
    for option in testfor:
        if os.path.exists(scratch_dir+'test'+option):
            modsuf=option
            os.remove(scratch_dir+'test'+option)
            print 'Autodetected modsuf:',modsuf
    if modsuf == None or retcode != 0:
        print 'Failed to auto-detect modsuf flag for compiler:'+compiler
        print "Usually it is .mod or .M but these didn't seem to work"
        print 'Please, specify it manually by invoking this program with the --modsuf option'
    shutil.rmtree(scratch_dir,ignore_errors=True)

if modpath == None: # Try to autodetect modpath option
    shutil.rmtree(scratch_dir,ignore_errors=True)
    os.mkdir(scratch_dir)
    os.mkdir(scratch_dir+dirsep+'modu')
    f=open(scratch_dir+dirsep+'modu/modu.f90','w')
    f.write('module modu\nend module modu\n')
    f.close()    
    f=open(scratch_dir+'test.f90','w')
    f.write('use modu\nend\n')
    f.close()    
    os.chdir(scratch_dir+dirsep+'modu')
    try:
        retcode=subprocess.call([compiler,cswitch,'modu.f90'],stdout=fnull,stderr=fnull)
    except:
        retcode=10
    options=['-I','-module','-L']
    for option in options:
        os.chdir(cwd)
        os.chdir(scratch_dir)
        try:
            retcode2=subprocess.call([compiler,option,'modu','test.f90'],stdout=fnull,stderr=fnull)
        except:
            retcode2=10
        if os.path.exists('a.out') and retcode == 0 and retcode2 == 0:
            modpath=option+' '
            print 'Autodetected modpath:',modpath
            break
    os.chdir(cwd)
    if modpath == None:
        print 'Failed to auto-detect modpath flag for compiler:'+compiler
        print "Usually it is -module or -I but these didn't seem to work"
        print 'Please, specify it manually by invoking this program with the --modpath option'
    shutil.rmtree(scratch_dir,ignore_errors=True)


if (cswitch == None or modpath == None or modsuf == None):
    print 'Error autodetecting compiler flags'
    print 'For cswitch found ',cswitch
    print 'For modpath found ',modpath
    print 'For modsuf found ',modsuf
    usage()
    sys.exit(1)

print '**** Using '+compiler+' ****'
print ''

# Endianness
print 'Testing system byte order (endianness)'
endianness=sys.byteorder
print 'Found '+endianness+' endian system'
if endianness != 'little' and endianness != 'big':
    print 'sys.byteorder returned an unknown value for machine endianness'
    sys.exit(1)

# Check that compiler/flags work
print ''
print 'Testing compiler and flags with Hello World program'
shutil.rmtree(scratch_dir,ignore_errors=True)
os.mkdir(scratch_dir)
f=open(scratch_dir+'test.f90','w')
f.write('print *,"Hello World"\nend\n')
f.close()
os.chdir(scratch_dir)
cmd=list()
cmd.append(compiler)
cmd.extend(otherflags.split())
cmd.append('test.f90')
print 'Command:'," ".join(cmd)
try:
    retcode=subprocess.call(cmd,stdout=fnull,stderr=fnull)
    pipe=subprocess.Popen('./a.out',stdout=subprocess.PIPE)
except:
    print '*** ERROR! Your compiler did not produce a suitable executable!'
    os.chdir(cwd)
    shutil.rmtree(scratch_dir,ignore_errors=True)
    sys.exit(1)
print 'Ok. Compiler and flags seem to work'
os.chdir(cwd)
shutil.rmtree(scratch_dir,ignore_errors=True)


# Find recl for current architecture
print ''
print 'Trying to determine recl parameter'
if recl == -1:
    shutil.rmtree(scratch_dir,ignore_errors=True)
    os.mkdir(scratch_dir)
    f=open(scratch_dir+'test.f90','w')
    f.write('integer::i\nreal(kind=8)::number\nnumber=1.\ndo i=1,257\nopen (35,file=\'test.dat\',access=\'direct\',recl=i)\nwrite (35,rec=1,err=100) number\nclose (35)\nprint *,i\nstop\n100  close (35)\nend do\nprint *,-1\nend\n')
    f.close()
    os.chdir(scratch_dir)
    cmd=list()
    cmd.append(compiler)
    cmd.extend(otherflags.split())
    cmd.append('test.f90')
    retcode=subprocess.call(cmd,stdout=fnull,stderr=fnull)
    pipe=subprocess.Popen('./a.out',stdout=subprocess.PIPE)
    recl=pipe.stdout.read()
    recl=recl[0:len(recl)-1] # Drop the newline at the end
    os.chdir(cwd)
    shutil.rmtree(scratch_dir,ignore_errors=True)
    if recl==-1:
        print '**** Error, impossible to determine recl automatically'
        print 'Proceeding assuming recl=4'
        print 'Note that this could result in your code not being able to properly'
        print 'read and write binary files and may lead to crashes!!'
        print 'You are strongly advised to use the --recl switch to specify it manually'
        recl=4
print 'Using recl = ',recl

# Start code-specific block
# Files that contain main programs
treetop='..'+dirsep
mains=[treetop+'main'+dirsep+'nicole.f90']

ignoredeps=['mpif.h'] # To allow compilation of mpi version
#ignoredeps.append('NRTYPE;') # Avoid crash in svdcmp.f90 because we can't handle ;

# Check Numerical Recipes routines
print 'Checking for availability of required routines from the Numerical Recipes'
files=['svdcmp','svbksb','pythag','nr','nrutil','nrtype','convlv','twofft','four1','realft','fourrow']
num_recipes=1
for file in files:
    if not os.path.isfile(treetop+'numerical_recipes'+dirsep+file+'.f90'):
        print 'Missing Numerical Recipes routine: '+treetop+file
        num_recipes=0
if num_recipes != 1:
    print 'NICOLE cannot be compiled without the required routines'
    print 'from the Numerical Recipes. Aborting...'
    sys.exit(1)
else:
    print'  ...Ok'

# Modify Numerical Recipes routines to suit NICOLE
print 'Modifying Numerical Recipes routines'
for file in files:
    source=open(treetop+'numerical_recipes'+dirsep+file+'.f90').read()
    source=re.sub(';','\n',source)
    f=open(treetop+'numerical_recipes'+dirsep+file+'.f90','w')
    f.write(source)
    f.close()
source=open(treetop+'numerical_recipes'+dirsep+'svdcmp.f90').read()
source=re.sub('USE nr,','USE Debug_module \n USE   nr,',source)
source=re.sub('call nrerror(.*)','then \n Debug_errorflags(flag_svdcmp)=1 \n   return \n   end if',source)
f=open(treetop+'numerical_recipes'+dirsep+'svdcmp.f90','w')
f.write(source)
f.close()
source=open(treetop+'numerical_recipes'+dirsep+'ludcmp.f90').read()
source=re.sub('SUBROUTINE ludcmp','SUBROUTINE ludcmp_dp',source)
source=re.sub('USE nrtype;','USE nrtype; USE Debug_module ; USE nrutil2;',source)
source=re.sub('REAL\(SP\)','REAL(DP)',source)
source=re.sub('\+imaxloc','+imaxloc_dp',source)
source=re.sub('call swap','call swap_dp',source)
source=re.sub('call nrerror(.*)','then \n Debug_errorflags(flag_ludcmp)=1 \n   return \n   end if',source)
f=open(treetop+'numerical_recipes'+dirsep+'ludcmp_dp.f90','w')
f.write(source)
f.close()
source=open(treetop+'numerical_recipes'+dirsep+'lubksb.f90').read()
source=re.sub('SUBROUTINE lubksb','SUBROUTINE lubksb_dp',source)
source=re.sub('REAL\(SP\)','REAL(DP)',source)
f=open(treetop+'numerical_recipes'+dirsep+'lubksb_dp.f90','w')
f.write(source)
f.close()

# Check if the FLUSH statement is available
print 'Testing availability of FLUSH statement...'
flush=0
shutil.rmtree(scratch_dir,ignore_errors=True)
os.mkdir(scratch_dir)
f2=open(scratch_dir+'test.f90','w')
f2.write('call flush(1)\nend')
f2.close()
os.chdir(scratch_dir)
cmd=list()
cmd.append(compiler)
cmd.extend(otherflags.split())
cmd.append('test.f90')
pipe=subprocess.Popen(cmd, stdout=subprocess.PIPE, \
                          stderr=subprocess.PIPE)
pipe.wait()
os.chdir(cwd)
if os.path.exists(scratch_dir+'a.out'):
    flush=1
    os.remove(scratch_dir+'a.out')
    print '   Flush is available'
else:
    print '   Flush is not available'
shutil.rmtree(scratch_dir,ignore_errors=True)

f=open(treetop+'misc'+dirsep+'myflush.f90','w')
source=open(treetop+'misc'+dirsep+'myflush.presource').read()
if flush == 0:
    source=re.sub("(?i)call flush\ *\(.*\)","",source)
f.write(source)
f.close()

if sopa == 0: # Do not compile SOPA (safer, default option)
    f=open(treetop+'forward'+dirsep+'sopa'+dirsep+'sopa.f90','w')
    source=open(treetop+'forward'+dirsep+'sopa'+dirsep+'nosopa.presource').read()
    f.write(source)
    f.close()
else: # Compile with SOPA
    f=open(treetop+'forward'+dirsep+'sopa'+dirsep+'sopa.f90','w')
    source=open(treetop+'forward'+dirsep+'sopa'+dirsep+'sopa.presource').read()
    f.write(source)
    f.close()

if mpi == 0: # Produce serial version
    print '*** Producing serial version ***'
    f=open(treetop+'main'+dirsep+'nicole.f90','w')
    source=open(treetop+'main'+dirsep+'nicole.presource').read()
    source=re.sub("!NOMPI","",source)
    source=re.sub("!MPI.*\n","",source)
    f.write(source)
    f.close()
    f=open(treetop+'time_code'+dirsep+'profiling.f90','w')
    source = open(treetop+'time_code'+dirsep+'profiling.presource').read()
    source=re.sub("!NOMPI","",source)
    source=re.sub("!MPI.*\n","",source)
    f.write(source)
    f.close()
else: # Produce MPI parallel version
    print '*** Producing MPI parallel version ***'
    f=open(treetop+'main'+dirsep+'nicole.f90','w')
    source = open(treetop+'main'+dirsep+'nicole.presource').read()
    source=re.sub("!MPI","",source)
    source=re.sub("!NOMPI.*\n","",source)
    f.write(source)
    f.close()
    f=open(treetop+'time_code'+dirsep+'profiling.f90','w')
    source = open(treetop+'time_code'+dirsep+'profiling.presource').read()
    source=re.sub("!MPI","",source)
    source=re.sub("!NOMPI.*\n","",source)
    f.write(source)
    f.close()
# Set correct RealBytes and Endian parameters in param_struct.f90
source=open(treetop+'main'+dirsep+'param_struct.f90').read()
source=re.sub("  Integer, Parameter :: RealBytes=.* !","  Integer, Parameter :: RealBytes="+recl+" !",source)
if endianness == 'little': boolean='.True.'
if endianness == 'big': boolean='.False.'
source=re.sub("Logical, Parameter :: LittleEndian=.* !","Logical, Parameter :: LittleEndian="+boolean+" !",source)
f=open(treetop+'main'+dirsep+'param_struct.f90','w')
f.write(source)
f.close()
# End of code-specific block

print ''
print '**** Analyzing source tree for dependences and modules ****'
# Create list of all source files
source_files=['']
for (basepath, children) in walktree(treetop,False):
    for child in children:
        if child.endswith('.f90'):
            source_files.append(os.path.join(basepath, child))
source_files=source_files[1:]
directories=['']
for file in source_files: # Directories in the source tree
    fpath, fname = os.path.split(file)
    if not fpath in directories:
        directories.append(fpath)
directories=directories[1:]

# Analyze dependencies and modules

# First find which files contain which modules
contains=['']
for file in source_files:
    contains.append([''])
contains=contains[1:len(contains)]
modulename=['']
modulefile=['']
for ifile in range(len(source_files)):
    file=source_files[ifile]
    f = open(file,'r')
    while True:
        line=f.readline()
        if len(line) == 0: break
        line=line.split('!')[0] # Remove comments
        line=line.split('\n')[0] # Remove end-of-line CR
        line=line.replace(',',' ') # Replace commas with spaces
        line=line.upper() # Upper case
        if re.match('[ ,\t]*END[ ,\t]',line):
            line='' # Remove END statements
        if re.match('[ ,\t]*MODULE[ ,\t]PROCEDURE[ ,\t]',line):
            line='' # Remove MODULE PROCEDURE statements
        tmp=''
        if re.match('[ ,\t]*MODULE[ ,\t]',line):
            tmp=line.split('MODULE ')
        tmp=line.split('MODULE ')
        if len(tmp) > 1 and re.search('[A-z,0-9]',tmp[0]) == None:
            module_defined=tmp[1].split()[0]
            contains[ifile].append(module_defined)
    contains[ifile]=contains[ifile][1:len(contains[ifile])]
    if (len (contains[ifile]) > 0):
        if (not quiet): print 'File: '+file+' contains the following module(s):'
        for module in contains[ifile]:
            if (not quiet): print '  '+module
            if module in modulename:
                print '      ** Error: Duplicated module *************'
                print '    Module name:',module
                print '    Found in files:'+modulefile[modulename.index(module)]+' and '+file
                print '      *****************************************'
                sys.exit(1)
            modulename.append(module)
            modulefile.append(file)
    f.close()

# Now figure out dependencies
requires=['']
missing=['']
for file in source_files:
    requires.append([''])
requires=requires[1:len(requires)]
for ifile in range(len(source_files)):
    file=source_files[ifile]
    f = open(file,'r')
    while True:
        line=f.readline()
        if len(line) == 0: break
        line=line.split('!')[0] # Remove comments
        line=line.split('\n')[0] # Remove end-of-line CR
        line=line.replace(',',' ') # Replace commas with spaces
        origline=line
        line=line.upper() # Upper case
        if re.match('[ ,\t]*END[ ,\t]',line):
            line='' # Remove END statements
        tmp=''
        if re.match('[ ,\t]*USE[ ,\t]',line):
            tmp=line.split('USE ')
        if len(tmp) > 1 and re.search('[A-z,0-9]',tmp[0]) == None:
            for tmp2 in tmp[1:]:
                if re.search(';',tmp2):
                    tmp2=re.sub(';','',tmp2)
                required_module=tmp2.split()[0]
                if not required_module in modulename:
                    if not required_module in missing and \
                            not required_module in ignoredeps:
                        print '     ** Error: Required module not found ******'
                        print '   Module name:',required_module
                        print '   Required by file:',file
                        missing.append(required_module)
                    if (not ignore) and \
                            not required_module in ignoredeps:
                        print '   Use the --ignore switch to continue creating the makefile'
                        sys.exit(1)
                else:
                    required_file=modulefile[modulename.index(required_module)]
                    if not required_file in requires[ifile] and \
                            required_file != file:
                        requires[ifile].append(required_file)
        tmp=''
        if re.match('(?i)[ ,\t]*INCLUDE[ ,\t]',origline):
            tmp=re.split('(?i)INCLUDE ',origline)
        if len(tmp) > 1 and re.search('[A-z,0-9]',tmp[0]) == None:
            tmp[1]=tmp[1].replace('"',' ')
            tmp[1]=tmp[1].replace("'",' ')
            required_module=tmp[1].split()[0]
            if not required_module in modulename:
                if not required_module in missing and \
                        not required_module in ignoredeps:
                    print '     ** Error: Required included file not found ******'
                    print '   Required file:',required_module
                    print '   Required by file:',file
                    missing.append(required_module)
                if (not ignore) and \
                        not required_module in ignoredeps: 
                    print '   Use the --ignore switch to continue creating the makefile'
                    sys.exit(1)
            else:
                required_file=modulefile[modulename.index(required_module)]
                if not required_file in requires[ifile] and \
                        required_file != file:
                    requires[ifile].append(required_file)
    requires[ifile]=requires[ifile][1:len(requires[ifile])]
    if (len (requires[ifile]) > 0):
        if (not quiet): print 'File: '+file+' depends on the following file(s):'
        for dependency in requires[ifile]:
            if (not quiet): print '  '+dependency
    f.close()
# If requested, show dependency tree
if showtree:
    for ifile in range(len(source_files)):
        print 'File: ',source_files[ifile]
        printsubtree(source_files,requires,ifile,1)

# Finished analyzing dependencies in the source tree
print '**** Finished analyzing dependences in source tree'
print ''
missing=missing[1:len(missing)]
if len(missing) > 0:
    print '***** WARNING: The following modules or files were not found in'
    print 'the source tree:'
    for m in missing: print '   '+m
    print 'Continuing anyway...'
    print ''


# Create makefile
ans=''
if (not yes):
    while (ans != 'y' and ans != 'n'):
        ans=raw_input("About to overwrite your makefile. Is this ok? [y/n] ")
if ans == 'n':
    sys.exit(0)

f = open('makefile', 'w')
f.write('# ************** AUTOGENERATED MAKEFILE ******************* \n')
f.write('#    Created on: '+datetime.datetime.now().ctime()+' \n')
f.write('# OS: '+os.name+'\n')
f.write('# Options used in create_makefile.py (some may have been autodetected): \n')
f.write('# --compiler='+compiler+'\n')
f.write('# --cswitch='+cswitch+'\n')
f.write('# --modpath='+modpath+'\n')
f.write('# --modsuf='+modsuf+'\n')
f.write('# --otherflags='+otherflags+'\n')
f.write('# --mpi='+str(mpi)+'\n')
f.write('# --recl='+recl+'\n')
f.write('# Flags:'+sysarg+'\n')
f.write('# ********************************************************* \n')
f.write('F90comp='+compiler+'\n')
f.write('CSwitch='+cswitch+'\n')
f.write('ModPath='+modpath+'\n')
f.write('ModSuf='+modsuf+'\n')
f.write('OPTIONS='+otherflags)
for dir in directories:
    f.write('  $(ModPath)'+dir)
f.write('\n')
f.write('RM='+rm+'\n')
f.write('CP='+cp+'\n')
f.write('\n')
f.write('.SUFFIXES: .f90 .o'+'\n')
f.write('\n')
f.write('.f90.o:'+'\n')
f.write('\t $(F90comp) $(CSwitch) $(OPTIONS) $< -o $@'+'\n')
f.write('\n')
tmp='all: '
for main in mains:
    fpath, fname=os.path.split(main)
    tmp=tmp+fname[0:len(fname)-4]+' '
f.write(tmp+' \n')
f.write('\n')
f.write('clean: \n')
allobj=directories[0]+dirsep+'*.o '
for dir in directories[1:]:
    allobj=allobj+' '+dir+dirsep+'*.o '
allmod=directories[0]+dirsep+'*$(ModSuf) '
for dir in directories[1:]:
    allmod=allmod+' '+dir+''+dirsep+'*$(ModSuf) '
f.write('\t $(RM) nicole '+allobj+allmod+' \n')
f.write('\n')
for main in mains:
    fpath, fname=os.path.split(main)
    tmp=fname[0:len(fname)-4]+': object \n'
    f.write(tmp)
    f.write('\t $(F90comp) $(OPTIONS) '+main+' '+allobj+' -o '+main[0:len(main)-4]+' \n\t $(RM) '+main[0:len(main)-4]+'.o \n')
    f.write('\n')
f.write('\n')
f.write('\n')
f.write('object: \\\n')
for file in source_files:
    if file not in mains:
        fileobj=file[0:len(file)-4]+'.o'
        f.write(fileobj+'    \\\n')
f.write('\n\n#Dependencies of individual files:\n')
for ifile in range(len(source_files)):
    file=source_files[ifile]
    if file not in mains:
        fileobj=file[0:len(file)-4]+'.o'
        f.write(fileobj+': ')
        for dep in requires[ifile]:
            depobj=dep[0:len(dep)-4]+'.o'
            f.write('  '+depobj+'  ')
        f.write('\n')
f.close()

print '\n'
print 'Done. The new makefile has been created\n\n'

fnull.close()
