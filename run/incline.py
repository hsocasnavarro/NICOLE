#!/usr/bin/env python
#
# Incline the model to put the line of sight along the z-direction
#
import struct
import sys
import numpy as n
import getopt
import model_prof_tools
import math
import sys

def usage():
    print 'Usage:'
    print "incline.py [-h, --help] [-q, --quiet] [-v, --verbose] [--xang=xxx] [--yang=xxx] [--dx=xxx] [--dy=xxx] [--model_in=filename] [--model_out=filename]"
    print "By default xang=0, yang=0"
    print "dx and dy in the same units as z"
    print "z scale is mandatory!!"
    print "modelin and modelout both default to __inputmodel.bin"

"""
NICOLE output tools
AUTHOR: Jaime de la Cruz Rodriguez (ITA-UiO 2011)

DESCRIPTION: The package is divided in two classes (model and profile) that read the output model
             and profiles from NICOLE.

USAGE:
       MODEL:
                import nicoletools as nic
                model = nic.model('modelout.mod', nx = 256, ny = 128) # nx and ny are optional keywords.
                model.read()
                # check the variables: model.z, model.tau, model.t, model.vz, model.bz, ...
       PROFILE:
                import nicoletools as nic
                prof = nic.profile('profiles.pro', nx = 256, ny = 128)
                prof.read()
                # check profiles: prof.i, prof.q, prof.u, prof.v

DEPENDENCIES: numpy 
--------------
MODIFICATIONS:
               2014-04-02, JdlCR: added support for model 2.6 and corrected buggy behaviour in tilt_model. 
"""
#
# MODEL CLASS
# 
class model:
    def __init__(self, name, nx = 0, ny = 0, verbose = 0):
        self.filename = name
        self.xy = n.int8(0)
        #
        # Read header
        self.s1,self.s2,self.npix,self.nz  = n.fromfile(name, dtype = 'int64', count = 4)
        if self.s1 == 3328834590979877230 and self.s2 == 2314885530823516723: self.xy+=1
        if self.s1 == 3328834590979877230 and self.s2 == 2314885530823516726: self.xy+=2
        self.s1,self.s2,self.s1,self.s2,self.nx, self.ny  = n.fromfile(name, dtype = 'int32', count = 6)
        self.npix=self.nx*self.ny
        [nx,ny,npix,nz]=[self.nx,self.ny,self.npix,self.nz]
        #
        #
        if verbose:
            print 'class::model::init : nz =   ', self.nz
            print 'class::model::init : npix = ', self.npix
        if (self.xy >= 1 ):
            if verbose:
                print 'class::model::init : nx = ', self.nx
                print 'class::model::init : ny = ', self.ny
                print 'class::model::init : mtype = ', self.xy
            self.reshape = [self.nx, self.ny, self.nz]
            reshape2 = [self.nx, self.ny]
        else:
            self.reshape = [self.npix, self.nz]
            reshape2 = [self.npix]
        #
        # create arrays for model 
        self.z = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.tau = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.t = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.el_p = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.gas_p = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.rho = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.vlos = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.vx = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.vy = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.vz = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.bz = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.bx = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.by = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.blx = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.bly = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.blz = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.v_mic = n.zeros(self.nz * self.npix, dtype = 'float32').reshape(self.reshape)
        self.v_mac = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.stray = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.ffactor = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.nH = n.zeros(self.nz *self.npix, dtype = 'float32').reshape(self.reshape)
        self.nHminus =  n.zeros(self.nz *self.npix, dtype = 'float32').reshape(self.reshape)
        self.nHplus =  n.zeros(self.nz *self.npix, dtype = 'float32').reshape(self.reshape)
        self.nH2 =  n.zeros(self.nz *self.npix, dtype = 'float32').reshape(self.reshape)
        self.nh2plus =  n.zeros(self.nz *self.npix, dtype = 'float32').reshape(self.reshape)
        self.keep_el_p = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_gas_p = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_rho = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_nH = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_nHminus = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_nHplus = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_nH2 = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        self.keep_nH2plus = n.zeros(self.npix, dtype = 'float32').reshape(reshape2)
        if(self.xy >= 1):
            self.abundance = n.zeros(self.npix*92, dtype = 'float32').reshape([self.nx,self.ny,92])
        else:
            self.abundance = n.zeros(self.npix*92, dtype = 'float32').reshape([self.npix,92])

#        self.abundance = n.zeros
#            

    def check_types(self):

# Check size of default types. Need 32-bit integers, 64-bit integers and floats
# 32-bit integers
        if struct.calcsize('<i') == 4:
            int4f='i'
        elif struct.calcsize('<l') == 4:
            int4f='l'
        elif struct.calcsize('<q') == 4:
            int4f='q'
        else:
            print 'This architecture has a non-standard integer size'
            print 'which is not supported in this version. Please, contact the'
            print 'author to get support for this platform.'
            sys.exit(1)
# 64-bit integers
        if struct.calcsize('<i') == 8:
            intf='i'
        elif struct.calcsize('<l') == 8:
            intf='l'
        elif struct.calcsize('<q') == 8:
            intf='q'
        else:
            print 'This architecture has a non-standard integer size'
            print 'which is not supported in this version. Please, contact the'
            print 'author to get support for this platform.'
            sys.exit(1)
# 64-bit floats
        if struct.calcsize('<f') == 8:
            flf='f'
        elif struct.calcsize('<d') == 8:
            flf='d'
        else:
            print 'This architecture has a non-standard float size'
            print 'which is not supported in this version. Please, contact the'
            print 'author to get support for this platform.'
            sys.exit(1)
        return [int4f,intf,flf]

    def read(self,verbose=1):
        [int4f,intf,flf]=self.check_types()
        if(self.xy == 1): sizerec = 17 * self.nz + 3
        if(self.xy == 2): sizerec = 22 * self.nz + 11 + 92
            
        ndat  = sizerec
        form = str(ndat)+'f8'
        #
        if verbose:
            print 'class::model::read : reading ', self.filename
        #
        fd = open(self.filename, 'rb')
        fd.seek(0)
        #
        # discard header
 
        
        tmp = struct.unpack('<'+str(sizerec)+flf,fd.read(sizerec*8))
        # loop and read
        if(self.xy == 1):
            for jj in n.arange(self.nx):
                for ii in n.arange(self.ny):
                    #
                    # Read pixel
                    tmp = struct.unpack('<'+str(sizerec)+flf,fd.read(sizerec*8))
                    
                    #
                    # Copy to vars
                    self.z[jj,ii,:] =       tmp[0 : self.nz]
                    self.tau[jj,ii,:] =     tmp[1 * self.nz : 2 * self.nz]
                    self.t[jj,ii,:] =       tmp[2 * self.nz : 3 * self.nz]
                    self.gas_p[jj,ii,:] =   tmp[3 * self.nz : 4 * self.nz]
                    self.rho[jj,ii,:] =     tmp[4 * self.nz : 5 * self.nz]
                    self.el_p[jj,ii,:] =    tmp[5 * self.nz : 6 * self.nz]
                    self.vlos[jj,ii,:] =    tmp[6 * self.nz : 7 * self.nz]
                    self.v_mic[jj,ii,:] =   tmp[7 * self.nz : 8 * self.nz]
                    self.bz[jj,ii,:] =      tmp[8 * self.nz : 9 * self.nz]
                    self.bx[jj,ii,:] =      tmp[9 * self.nz : 10 * self.nz]
                    self.by[jj,ii,:] =      tmp[10 * self.nz: 11 * self.nz]
                    self.blx[jj,ii,:] =     tmp[11 * self.nz: 12 * self.nz]
                    self.bly[jj,ii,:] =     tmp[12 * self.nz: 13 * self.nz]
                    self.blz[jj,ii,:] =     tmp[13 * self.nz: 14 * self.nz]
                    self.vx[jj,ii,:] =      tmp[14 * self.nz: 15 * self.nz]
                    self.vy[jj,ii,:] =      tmp[15 * self.nz: 16 * self.nz]
                    self.vz[jj,ii,:] =      tmp[16 * self.nz: 17 * self.nz]
                    self.v_mac[jj,ii] =     tmp[17 * self.nz    ]
                    self.stray[jj,ii] =     tmp[17 * self.nz + 1]
                    self.ffactor[jj,ii] =   tmp[17 * self.nz + 2]
        elif(self.xy == 2):
         for jj in n.arange(self.nx):
                for ii in n.arange(self.ny):
                    #
                    # Read pixel
                    tmp = struct.unpack('<'+str(sizerec)+flf,fd.read(sizerec*8))
                    
                    #
                    # Copy to vars
                    self.z[jj,ii,:] =       tmp[0 : self.nz]
                    self.tau[jj,ii,:] =     tmp[1 * self.nz : 2 * self.nz]
                    self.t[jj,ii,:] =       tmp[2 * self.nz : 3 * self.nz]
                    self.gas_p[jj,ii,:] =   tmp[3 * self.nz : 4 * self.nz]
                    self.rho[jj,ii,:] =     tmp[4 * self.nz : 5 * self.nz]
                    self.el_p[jj,ii,:] =    tmp[5 * self.nz : 6 * self.nz]
                    self.vlos[jj,ii,:] =    tmp[6 * self.nz : 7 * self.nz]
                    self.v_mic[jj,ii,:] =   tmp[7 * self.nz : 8 * self.nz]
                    self.bz[jj,ii,:] =      tmp[8 * self.nz : 9 * self.nz]
                    self.bx[jj,ii,:] =      tmp[9 * self.nz : 10 * self.nz]
                    self.by[jj,ii,:] =      tmp[10 * self.nz: 11 * self.nz]
                    self.blx[jj,ii,:] =     tmp[11 * self.nz: 12 * self.nz]
                    self.bly[jj,ii,:] =     tmp[12 * self.nz: 13 * self.nz]
                    self.blz[jj,ii,:] =     tmp[13 * self.nz: 14 * self.nz]
                    self.vx[jj,ii,:] =      tmp[14 * self.nz: 15 * self.nz]
                    self.vy[jj,ii,:] =      tmp[15 * self.nz: 16 * self.nz]
                    self.vz[jj,ii,:] =      tmp[16 * self.nz: 17 * self.nz]
                    self.nH[jj,ii,:] =      tmp[17 * self.nz: 18 * self.nz]
                    self.nHminus[jj,ii,:] = tmp[18 * self.nz: 19 * self.nz]
                    self.nHplus[jj,ii,:] =  tmp[19 * self.nz: 20 * self.nz]
                    self.nH2[jj,ii,:] =     tmp[20 * self.nz: 21 * self.nz]
                    self.nh2plus[jj,ii,:] = tmp[21 * self.nz: 22 * self.nz]
                    self.v_mac[jj,ii] =     tmp[22*self.nz]
                    self.stray[jj,ii] =     tmp[22*self.nz+1]
                    self.ffactor[jj,ii] =   tmp[22*self.nz+2]
                    self.keep_el_p[jj,ii] = tmp[22*self.nz+3]
                    self.keep_gas_p[jj,ii]= tmp[22*self.nz+4]
                    self.keep_rho[jj,ii] =  tmp[22*self.nz+5]
                    self.keep_nH[jj,ii] =   tmp[22*self.nz+6]
                    self.keep_nHminus[jj,ii] = tmp[22*self.nz+7]
                    self.keep_nHplus[jj,ii] =  tmp[22*self.nz+8]
                    self.keep_nH2[jj,ii] =     tmp[22*self.nz+9]
                    self.keep_nH2plus[jj,ii] = tmp[22*self.nz+10]
                    self.abundance[jj,ii,:] = tmp[22*self.nz+11:22*self.nz+11+92]
                    
        else:
            for ii in n.arange(self.npix):
                #
                # Read pixel
                tmp = ((n.core.records.fromfile(fd, formats = form, shape = 1))[0])[0]
                #
                self.z[ii,:] =       tmp[0 : self.nz]
                self.tau[ii,:] =     tmp[    self.nz : 2 * self.nz]
                self.t[ii,:] =       tmp[2 * self.nz : 3 * self.nz]
                self.gas_p[ii,:] =   tmp[3 * self.nz : 4 * self.nz]
                self.rho[ii,:] =     tmp[4 * self.nz : 5 * self.nz]
                self.el_p[ii,:] =    tmp[5 * self.nz : 6 * self.nz]
                self.vz[ii,:] =      tmp[6 * self.nz : 7 * self.nz]
                self.v_mic[ii,:] =   tmp[7 * self.nz : 8 * self.nz]
                self.blz[ii,:] =      tmp[8 * self.nz : 9 * self.nz]
                self.blx[ii,:] =      tmp[9 * self.nz : 10 * self.nz]
                self.bly[ii,:] =      tmp[10 * self.nz : 11 * self.nz]
                self.v_mac[ii] =     tmp[13 * self.nz    ]
                self.stray[ii] =     tmp[13 * self.nz + 1]
                self.expansion[ii] = tmp[13 * self.nz + 2]
#
                self.blong[ii,:]=self.blz[ii,:]
                self.bx[ii,:]=self.blx[ii,:]
                self.by[ii,:]=self.bly[ii,:]
                self.vlos[ii,:]=self.vz[ii,:]
#                self.vy[ii,:]=self.vx[ii,:]
#                self.vz[ii,:]=self.vy[ii,:]
        #
        # Close file
        fd.close()

    def write(self,filename,verbose=1):

        [int4f,intf,flf]=self.check_types()
        [nx,ny,nz]=[self.nx,self.ny,self.nz]
        sizerec = 22 * nz + 11 + 92

        fd = open(filename, 'wb')
        ndat  = (17 * self.nz + 3)
        fd.write(struct.pack('<16s'+int4f+int4f+intf,'nicole2.6bm     ',nx,ny,nz))
        for i in range(sizerec-16/8-1-1): 
                fd.write(struct.pack('<'+flf,0.)) # Fill record
        #
        Vector=range(sizerec)
        percent=-1
        for ix in range(nx):
            for iy in range(ny):
                Vector[0:nz]=self.z[ix,iy,0:nz]
                Vector[nz:2*nz]=self.tau[ix,iy,0:nz]
                Vector[2*nz:3*nz]=self.t[ix,iy,0:nz]
                Vector[3*nz:4*nz]=self.gas_p[ix,iy,0:nz]
                Vector[4*nz:5*nz]=self.rho[ix,iy,0:nz]
                Vector[5*nz:6*nz]=self.el_p[ix,iy,0:nz]
                Vector[6*nz:7*nz]=self.vlos[ix,iy,0:nz]
                Vector[7*nz:8*nz]=self.v_mic[ix,iy,0:nz]
                Vector[8*nz:9*nz]=self.bz[ix,iy,0:nz]
                Vector[9*nz:10*nz]=self.bx[ix,iy,0:nz]
                Vector[10*nz:11*nz]=self.by[ix,iy,0:nz]
                Vector[11*nz:12*nz]=self.blx[ix,iy,0:nz]
                Vector[12*nz:13*nz]=self.bly[ix,iy,0:nz]
                Vector[13*nz:14*nz]=self.blz[ix,iy,0:nz]
                Vector[14*nz:15*nz]=self.vx[ix,iy,0:nz]
                Vector[15*nz:16*nz]=self.vy[ix,iy,0:nz]
                Vector[16*nz:17*nz]=self.vz[ix,iy,0:nz]
                Vector[17 * self.nz: 18 * self.nz] =      self.nH[ix,iy,:] 
                Vector[18 * self.nz: 19 * self.nz] =      self.nHminus[ix,iy,:] 
                Vector[19 * self.nz: 20 * self.nz] =      self.nHplus[ix,iy,:] 
                Vector[20 * self.nz: 21 * self.nz] =      self.nH2[ix,iy,:] 
                Vector[21 * self.nz: 22 * self.nz] =      self.nh2plus[ix,iy,:] 
                Vector[22*self.nz] =		       self.v_mac[ix,iy] 
                Vector[22*self.nz+1] =		       self.stray[ix,iy] 
                Vector[22*self.nz+2] =		       self.ffactor[ix,iy] 
                Vector[22*self.nz+3] =		       self.keep_el_p[ix,iy] 
                Vector[22*self.nz+4] = 		       self.keep_gas_p[ix,iy]
                Vector[22*self.nz+5] =		       self.keep_rho[ix,iy] 
                Vector[22*self.nz+6] =		       self.keep_nH[ix,iy] 
                Vector[22*self.nz+7] =		       self.keep_nHminus[ix,iy]
                Vector[22*self.nz+8] =		       self.keep_nHplus[ix,iy] 
                Vector[22*self.nz+9] =		       self.keep_nH2[ix,iy]
                Vector[22*self.nz+10] =		       self.keep_nH2plus[ix,iy]
                Vector[22*self.nz+11:22*self.nz+11+92] =  self.abundance[ix,iy,:] 
                for d in Vector: fd.write(struct.pack('<'+flf,d))
        if verbose:
          if (int((ix*ny+(iy+1))*100./nx/ny) > percent):
            percent=int((ix*ny+(iy+1))*100./nx/ny)
            sys.stdout.write('\r'+'Writing model...'+str(percent)+'%')
            sys.stdout.flush()
        fd.close() # Done with input model  
        print ''


#
# PROFILE CLASS
#
class profile:
    def __init__(self, name, nx = 0 , ny = 0):
        self.filename = name
        #
        # Read header
        self.s1, self.s2, self.npix, self.nwav, self.nx, self.ny  = n.fromfile(name, dtype = 'int64', count = 6)
        self.xy = n.int8(0)
        #
        # 
        if(((nx*ny) != 0) and ((self.nx * self.ny) ==0) and ((nx * ny) == self.npix)):
            self.nx = nx
            self.ny = ny
        #
        print 'class::profile::init : nwav = ', self.nwav
        print 'class::profile::init : npix = ', self.npix
        if((self.nx != 0 ) and (self.ny != 0)):
            print 'class::profile::init : nx = ', self.nx
            print 'class::profile::init : ny = ', self.ny
            self.xy+= 1
            self.reshape = [self.ny, self.nx, self.nwav]
        else:
            self.reshape = [self.npix, self.nwav]
        #
        # create arrays for model 
        self.i = n.zeros(self.nwav * self.npix, dtype = 'float32').reshape(self.reshape)
        self.q = n.zeros(self.nwav * self.npix, dtype = 'float32').reshape(self.reshape)
        self.u = n.zeros(self.nwav * self.npix, dtype = 'float32').reshape(self.reshape)
        self.v = n.zeros(self.nwav * self.npix, dtype = 'float32').reshape(self.reshape)
    #
    def read(self,verbose=1):
        ndat = 4 * self.nwav
        form = str(ndat)+'f8'
        nwav4 = self.nwav * 4
        #
        # open file
        fd = open(self.filename, 'r')
        fd.seek(0)
        if verbose:
            print 'class::profile::read : reading ', self.filename

        #
        # discard header
        tmp = n.core.records.fromfile(fd, formats = form, shape = 1)
        
        #
        # Read data
        if(self.xy !=0): 
            for jj in n.arange(self.ny):
                for ii in n.arange(self.nx):
                    #
                    # Read pixel
                    tmp = ((n.core.records.fromfile(fd, formats = form, shape = 1))[0])[0]
                    #
                    # Copy vars
                    self.i[jj, ii, :] = tmp[0:nwav4:4]
                    self.q[jj, ii, :] = tmp[1:nwav4:4]
                    self.u[jj, ii, :] = tmp[2:nwav4:4]
                    self.v[jj, ii, :] = tmp[3:nwav4:4]
        else:
            for ii in n.arange(self.npix):
                #
                # Read pixel
                tmp = ((n.core.records.fromfile(fd, formats = form, shape = 1))[0])[0]
                
                #
                # Copy vars
                self.i[ii, :] = tmp[0:nwav4:4]
                self.i[ii, :] = tmp[1:nwav4:4]
                self.i[ii, :] = tmp[2:nwav4:4]
                self.i[ii, :] = tmp[3:nwav4:4]
        #
        fd.close()

#
def tilt_model(m, xang = 0.0, yang = 0.0, dx = 0, dy = 0, verbose = 1):
    dtor = 0.017453292521655169 #3.14159265389793 / 180.
    #
    # X axis tilt
    # Modify depth scale along new rays
    z = m.z[0,0,:]
    x = n.arange(0,m.nx)*dx
    y = n.arange(0,m.ny)*dy
    mmax = n.max(x)
    mmay = n.max(y)
    
    mu = n.cos(n.sqrt(xang**2 + yang**2) * dtor)
    xmu = n.cos(dtor * xang)
    ymu = n.cos(dtor * yang)
    if verbose: print 'funct::tilt_model : mu = ', mu
    print 'funct::tilt_model : tilting model assuming periodic boundary conditions!'
    #
    zx2 = (z - z[0]) / abs(xmu) + z[0]
    zy2 = (z - z[0]) / abs(ymu) + z[0]
    #
    xx0, yy0 = n.meshgrid(x / mmax, y / mmay)
    points = (yy0.flatten(), xx0.flatten())
    percent = int(-1)
    ntot = 100. / (m.nz - 1.0 + 1e-4)
    #
    # Temporary storage for vector fields
    tmpz = n.copy(m.vlos)
    tmpx = n.copy(m.vx)
    tmpy = n.copy(m.vy)
    tmpz1 = n.copy(m.bz)
    tmpx1 = n.copy(m.bx)
    tmpy1 = n.copy(m.by)
    #
    for k in range(m.nz):
        #
        # X-axis grid (compute displacement of current layer)
        #
        ddx = -n.sign(xang) * n.abs(zx2[k] - zx2[0]) * \
            n.sqrt(1.0 - xmu * xmu) / mmax * (m.nx - 1.0)
        ddy = -n.sign(yang) * n.abs(zy2[k] - zy2[0]) * \
            n.sqrt(1.0 - ymu * ymu) / mmay * (m.ny - 1.0)
        #                
        # Interpolate variables
        #
        m.t[:,:,k] = bilint(m.t[:,:,k], ddx, ddy)
        #
        #if(n.min(m.rho[:,:,k]) > 0.0):
        #    m.rho[:,:,k] = n.exp(bilint(n.log(m.rho[:,:,k]), dx, dy))
        #else:
        m.rho[:,:,k] = bilint(m.rho[:,:,k], ddx, ddy)
        #
        #if(n.min(m.el_p[:,:,k]) > 0.0):
        #    m.el_p[:,:,k] = n.exp(bilint(n.log(m.el_p[:,:,k]), dx, dy))
        #else:
        m.el_p[:,:,k] = bilint(m.el_p[:,:,k], ddx, ddy)
        #
        #if(n.min(m.gas_p[:,:,k]) > 0.0):
        #    m.gas_p[:,:,k] = n.exp(bilint(n.log(m.gas_p[:,:,k]), dx, dy))
        #else:
        m.gas_p[:,:,k] = bilint(m.gas_p[:,:,k], ddx, ddy)
        #

        tmpz[:,:,k] = bilint(tmpz[:,:,k], ddx, ddy)
        tmpx[:,:,k] = bilint(tmpx[:,:,k], ddx, ddy)
        tmpy[:,:,k] = bilint(tmpy[:,:,k], ddx, ddy)
        
        tmpz1[:,:,k] = bilint(tmpz1[:,:,k], ddx, ddy)
        tmpx1[:,:,k] = bilint(tmpx1[:,:,k], ddx, ddy)
        tmpy1[:,:,k] = bilint(tmpy1[:,:,k], ddx, ddy)

        
        m.v_mic[:,:,k] = bilint(m.v_mic[:,:,k], ddx, ddy)
        #
        # counter
        #
        if(int(k*ntot) > percent) and verbose:
            percent = int(k * ntot)
            sys.stdout.write('\r'+'funct::tilt_model : slanting 3D model: '+str(percent)+'%')
            sys.stdout.flush()
    print(' ')
    sys.stdout.flush()
    
    # Proyect vectors
    print 'funct::tilt_model : projecting velocities and the magnetic field to the new LOS'
    xsign = n.sign(xmu)
    ysign = n.sign(ymu)
    #
    ymu2 = n.sqrt(1.0 - ymu * ymu)
    xmu2 = n.sqrt(1.0 - xmu * xmu)
    #
    m.bz[:,:,:] = tmpz1 * xmu * ymu - ysign * tmpy1 * xmu * ymu2 - xsign * tmpx1 * xmu2
    m.by[:,:,:] = tmpy1 * ymu + tmpz1 * ymu2 * ysign
    m.bx[:,:,:] = tmpx1 * xmu + (tmpz1 * ymu - ysign * tmpy1 * ymu2) * xmu2 * xsign
    #
    m.vlos[:,:,:] = tmpz * xmu * ymu - ysign * tmpy * xmu * ymu2 - xsign * tmpx * xmu2
    m.vy[:,:,:] =   tmpy * ymu + tmpz * ymu2 * ysign
    m.vx[:,:,:] =   tmpx * xmu + (tmpz * ymu - ysign * tmpy * ymu2) * xmu2 * xsign
    
    #
    # new height scale
    #
    for k in n.arange(m.nz):
        m.z[:,:,k] = (z[k] - z[0]) / abs(mu) + z[0]
    #
    return(m)
#
"""
Function bilint

PURPOSE: Fast bilinear interpolation of 2D data assuming periodic bourdary conditions.
         The data is assumed to be in a regular grid.

INPUT:
        var: input 2D array to be interpolated
         dx: shift in the x-axis (float)
         dy: shift in the y-axis (float)
         
AUTHOR: Jaime de la Cruz Rodriguez (ITA-UiO 2011)

"""
def bilint(var, dx, dy):
    dtor = 0.017453292521655169 #3.14159265389793 / 180.
    #
    dim = n.shape(var)
    img_int = n.zeros(dim)
    #
    var2= n.zeros([dim[0] + 1, dim[1] + 1])
    var2[0:-1, 0:-1] = var
    var2[0:-1,dim[1]] = var[:, dim[1]-1]
    var2[dim[0], 0:-1] = var[dim[0]-1, :]
    var2[dim[0], dim[1]] = var[0,0]
    #
    # Separate shift into integer and float parts.
    #
    if(dx < 0.0):
        idx = int(dx - 1.0)
        ddx = 1.0 + (dx - int(dx))
    else:
        idx = int(dx)
        ddx = dx - int(dx)
    if(dy < 0.0):
        idy = int(dy - 1.0)
        ddy = 1.0 + (dy - int(dy))
    else:
        idy = int(dy)
        ddy = dy - int(dy)
    #
    # compute interpolation coeffs at dx - int(dx) and dy-int(dy),
    # the non-integer part.
    #
    img_int[:,:] = var2[0:-1, 0:-1] * (1.0 - ddx) * (1 - ddy) + \
                   var2[1:: , 0:-1] * (1.0 - ddx) * ddy + \
                   var2[0:-1, 1:: ] * (1.0 - ddy) * ddx + \
                   var2[1:: , 1:: ] * ddx * ddy
                   
    #
    # Shift the array by the integer part of the shifts.
    #
    if(idx != 0):
        img_int = n.roll(img_int, -idx, axis = 1)
    if(idy != 0):
        img_int = n.roll(img_int, -idy, axis = 0)
    #
    #
    return(img_int)

# Get command-line arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],"hqv", ["help","quiet",
                      "verbose","xang=","yang=","dx=","dy=",
                      "model_in=","model_out="])
except:
    print 'Command-line option not recognized'
    usage()
    sys.exit(2)

verbose=1
xang=0.
yang=0.
dx=1
dy=1
modelin='__inputmodel.bin'
modelout='__inputmodel.bin'
for o,a in opts:
    if o == '-h' or o == '--help':
        usage()
        sys.exit(0)
    if o == '-q' or o == '--quiet':
        verbose=0
    if o == '-v' or o == '--verbose':
        verbose=2
    if o == '--xang':
        xang=float(a)
    if o == '--yang':
        yang=float(a)
    if o == '--dx':
        dx=float(a)
    if o == '--dy':
        dy=float(a)
    if o == '--model_in':
        modelin=a
    if o == '--model_out':
        modelout=a

if verbose == 2:
    print 'Inclining. Angles: x=',xang,' y=',yang
    print '           Horizontal scales: dx=',dx,' dy=',dy

m_in=model(modelin,verbose=verbose)
m_in.read()
tmdl=tilt_model(m_in, xang, yang, dx, dy,verbose=verbose)
tmdl.write(modelout)

sys.exit(0)
