#!/usr/bin/env python

# Tools to read/write NICOLE models and profiles

def check_types():
    import struct

# Check size of default types. Need 32-bit integers, 64-bit integers and floats
# 32-bit integers
    if struct.calcsize('<i') == 4:
        int4f='i'
    elif struct.calcsize('<l') == 4:
        int4f='l'
    elif struct.calcsize('<q') == 4:
        int4f='q'
    else:
        print ('This architecture has a non-standard integer size')
        print ('which is not supported in this version. Please, contact the')
        print ('author to get support for this platform.')
        sys.exit(1)
# 64-bit integers
    if struct.calcsize('<i') == 8:
        intf='i'
    elif struct.calcsize('<l') == 8:
        intf='l'
    elif struct.calcsize('<q') == 8:
        intf='q'
    else:
        print ('This architecture has a non-standard integer size')
        print ('which is not supported in this version. Please, contact the')
        print ('author to get support for this platform.')
        sys.exit(1)
# 64-bit floats
    if struct.calcsize('<f') == 8:
        flf='f'
    elif struct.calcsize('<d') == 8:
        flf='d'
    else:
        print ('This architecture has a non-standard float size')
        print ('which is not supported in this version. Please, contact the')
        print ('author to get support for this platform.')
        sys.exit(1)

    return [int4f,intf,flf]


import re
import os
import struct
import sys

def check_model(filename):

    def check_nicole_binary_format(header, filename):
        # Logic for checking NICOLE binary format
        string = header[0:11].decode('utf-8')  # Decode from bytes to string and remove padding zeros
        if string == 'nicole1.6bm':
            filetype = 'nicole1.6'
            nx = 1
            [_, ny, nz] = struct.unpack('<16s' + intf + intf, header)
            expected_size = (13 * nz + 3) * (nx * ny + 1) * 8
        elif string == 'nicole2.3bm':
            filetype = 'nicole2.3'
            [_, nx, ny, nz] = struct.unpack('<16s' + int4f + int4f + intf, header)
            expected_size = (17 * nz + 3) * (nx * ny + 1) * 8
        elif string == 'nicole2.6bm':
            filetype = 'nicole2.6'
            [_, nx, ny, nz] = struct.unpack('<16s' + int4f + int4f + intf, header)
            expected_size = (22 * nz + 3 + 8 + 92) * (nx * ny + 1) * 8
        elif string == 'nicole18.04' or string == 'nicole1804m':
            filetype = 'nicole18.04'
            [_, nx, ny, nz] = struct.unpack('<16s' + int4f + int4f + intf, header)
            expected_size = (22 * nz + 13 + 92) * (nx * ny + 1) * 8

        else:
            return None

        filesize = os.path.getsize(filename)
        if filesize != expected_size:
            print('Incorrect size of model file:', filename)
            print('The file is probably corrupted. Proceeding anyway...')

        return [filetype, nx, ny, nz]

    def check_ascii_format(filename):
        try:
            with open(filename, 'rb') as f:
                content = f.read()

            # Attempt to decode content
            try:
                content = content.decode('utf-8')
            except UnicodeDecodeError:
                content = content.decode('utf-8', errors='replace')

            if not content:
                print('Empty file!')
                print('filename:', filename)
                sys.exit(1)

            lines = content.split('\n')
            if re.search('(?i)format version', lines[0]):
                filetype = 'ascii' + lines[0].split()[-1]
                if filetype not in ['ascii1.0', 'ascii2.3']:
                    print('Model file:', filename)
                    print('seems to be an ASCII file of unsupported format')
                    print(lines[0])
                    sys.exit(1)

                nz = len(lines) - 2 - lines.count('')  # Subtract 1 for the first line

                return [filetype, 1, 1, nz]  # ASCII files are one-dimensional
        except Exception as e:
            print("Error processing ASCII format:", e)
            return None

    def check_idl_format(filename):
        try:
            with open(filename, 'rb') as f:
                stri = struct.unpack('10s', f.read(10))
                bytes = list(stri[0])

            if bytes == [83, 82, 0, 4, 0, 0, 0, 10, 0, 0]:
                try:
                    import idlsave
                    idl = idlsave.read(filename, verbose=0)
                    shape = idl.t.shape
                    ls = list(shape)
                    while len(ls) < 3:
                        ls.append(1)
                    idl.t = idl.t.reshape(ls)
                    nz = idl.nz if 'nz' in idl else idl.t.shape[0]
                    ny = idl.ny if 'ny' in idl else idl.t.shape[1]
                    nx = idl.nx if 'nx' in idl else idl.t.shape[2]
                    return ['idl', nx, ny, nz]
                except ImportError:
                    print('Error. Model file is an IDL save file:', filename)
                    print('Python modules IDLSave and Numpy are needed to read this kind of files')
                    print("The modules don't appear to be installed in your system")
                    sys.exit(1)
        except Exception as e:
            print("Error processing IDL format:", e)
            return None

    [int4f, intf, flf] = check_types()

    try:
        with open(filename, 'rb') as f:
            header = f.read(16 + 4 + 4 + 8)
            format_check = check_nicole_binary_format(header, filename)
            if format_check:
                return format_check

            format_check = check_ascii_format(filename)
            if format_check:
                return format_check

            format_check = check_idl_format(filename)
            if format_check:
                return format_check

    except FileNotFoundError:
        print('Could not find model file:', filename)
        return ['inexistent', -1, -1, -1]

    print('Unknown model file type:', filename)
    sys.exit(1)




import re
import os
import struct
import sys

def check_prof(filename):
    
    def check_ascii_format(filename):
        with open(filename, 'r') as f:
            content = f.read()

        if not content:
            print('Empty file!')
            print('filename:', filename)
            sys.exit(1)

        nprintable = sum(1 for char in content if char in string.printable or char in string.whitespace)
        if float(nprintable) / float(len(content)) > 0.95:
            content = re.sub('!.*\n', '\n', content)
            content = re.sub('(?i).*format version.*\n', '\n', content)
            lines = content.split('\n')
            nlam = len(lines) - 1 - lines.count('')
            return ['ascii', 1, 1, nlam]
        return None

    def check_nicole_binary_format(header, filename):
        readstr = header[0:11].decode('utf-8')  # Decode and remove padding zeros
        if readstr == 'nicole1.6bp':
            filetype = 'nicole1.6'
            nx = 1
            [_, ny, nlam] = struct.unpack('<16s' + intf + intf, header)
            expected_size = (4 * nlam) * (nx * ny + 1) * 8
        elif readstr == 'nicole2.3bp':
            filetype = 'nicole2.3'
            [_, nx, ny, nlam] = struct.unpack('<16s' + int4f + int4f + intf, header)
            expected_size = (4 * nlam) * (nx * ny + 1) * 8
        else:
            return None

        filesize = os.path.getsize(filename)
        if filesize != expected_size:
            print('Incorrect size of profile file:', filename)
            print('The file is probably corrupted. Proceeding anyway...')
        return [filetype, nx, ny, nlam]

    def check_idl_format(filename):
        with open(filename, 'rb') as f:
            stri = struct.unpack('10s', f.read(10))
            bytes = list(stri[0])

        if bytes == [83, 82, 0, 4, 0, 0, 0, 10, 0, 0]:
            try:
                import idlsave
                idl = idlsave.read(filename, verbose=0)
                shape = idl.stki.shape
                ls = list(shape)
                while len(ls) < 3:
                    ls.append(1)
                idl.stki = idl.stki.reshape(ls)
                nlam = getattr(idl, 'nlam', idl.stki.shape[0])
                ny = getattr(idl, 'ny', idl.stki.shape[1])
                nx = getattr(idl, 'nx', idl.stki.shape[2])
                return ['idl', nx, ny, nlam]
            except ImportError:
                print('Error. Profile file is an IDL save file:', filename)
                print('Python modules IDLSave and Numpy are needed to read this kind of files')
                print("The modules don't appear to be installed in your system")
                sys.exit(1)
        return None

    [int4f, intf, flf] = check_types()

    try:
        with open(filename, 'rb') as f:
            header = f.read(16 + 8 + 8)

        format_check = check_nicole_binary_format(header, filename)
        if format_check:
            return format_check

        format_check = check_ascii_format(filename)
        if format_check:
            return format_check

        format_check = check_idl_format(filename)
        if format_check:
            return format_check

    except FileNotFoundError:
        print('Could not find profile file:', filename)
        return ['inexistent', -1, -1, -1]

    print('Unknown profile file type:', filename)
    sys.exit(1)








import re
import os
import struct
import sys

def check_ascii_format(filename):
    try:
        with open(filename, 'rb') as f:
            content = f.read()

        # Attempt to decode content; use 'replace' or 'ignore' for error handling
        try:
            content = content.decode('utf-8')
        except UnicodeDecodeError:
            content = content.decode('utf-8', errors='replace')

        if not content:
            print('Empty file!')
            print('filename:', filename)
            sys.exit(1)

        nprintable = sum(1 for char in content if char in string.printable or char in string.whitespace)
        if float(nprintable) / float(len(content)) > 0.95:
            content = re.sub('!.*\n', '\n', content)
            content = re.sub('(?i).*format version.*\n', '\n', content)
            lines = content.split('\n')
            nlam = len(lines) - 1 - lines.count('')
            return ['ascii', 1, 1, nlam]
    except Exception as e:
        print("Error processing ASCII format:", e)
        return None

    def check_nicole_binary_format(header, filename):
        readstr = header[0:11].decode('utf-8')  # Decode and remove padding zeros
        if readstr == 'nicole1.6bp':
            filetype = 'nicole1.6'
            nx = 1
            [_, ny, nlam] = struct.unpack('<16s' + intf + intf, header)
            expected_size = (4 * nlam) * (nx * ny + 1) * 8
        elif readstr == 'nicole2.3bp':
            filetype = 'nicole2.3'
            [_, nx, ny, nlam] = struct.unpack('<16s' + int4f + int4f + intf, header)
            expected_size = (4 * nlam) * (nx * ny + 1) * 8
        else:
            return None

        filesize = os.path.getsize(filename)
        if filesize != expected_size:
            print('Incorrect size of profile file:', filename)
            print('The file is probably corrupted. Proceeding anyway...')
        return [filetype, nx, ny, nlam]

    def check_idl_format(filename):
        with open(filename, 'rb') as f:
            stri = struct.unpack('10s', f.read(10))
            bytes = list(stri[0])

        if bytes == [83, 82, 0, 4, 0, 0, 0, 10, 0, 0]:
            try:
                import idlsave
                idl = idlsave.read(filename, verbose=0)
                shape = idl.stki.shape
                ls = list(shape)
                while len(ls) < 3:
                    ls.append(1)
                idl.stki = idl.stki.reshape(ls)
                nlam = getattr(idl, 'nlam', idl.stki.shape[0])
                ny = getattr(idl, 'ny', idl.stki.shape[1])
                nx = getattr(idl, 'nx', idl.stki.shape[2])
                return ['idl', nx, ny, nlam]
            except ImportError:
                print('Error. Profile file is an IDL save file:', filename)
                print('Python modules IDLSave and Numpy are needed to read this kind of files')
                print("The modules don't appear to be installed in your system")
                sys.exit(1)
        return None

    [int4f, intf, flf] = check_types()

    try:
        with open(filename, 'rb') as f:
            header = f.read(16 + 8 + 8)

        format_check = check_nicole_binary_format(header, filename)
        if format_check:
            return format_check

        format_check = check_ascii_format(filename)
        if format_check:
            return format_check

        format_check = check_idl_format(filename)
        if format_check:
            return format_check

    except FileNotFoundError:
        print('Could not find profile file:', filename)
        return ['inexistent', -1, -1, -1]

    print('Unknown profile file type:', filename)
    sys.exit(1)

def read_model(filename, filetype, nx, ny, nz, ix, iy, sequential=0):
# Read a particular model from a model file 
# Note that ix,iy=0 represents the first model
# Use the sequential keyword to specify that this file is going to be
# read sequentially rather than accessing specific records. In that
# case, the values of ix and iy are irrelevant. The first record
# needs to be read with sequential=0 and ix,iy = 0
# Model is returned WITHOUT abundances
    import struct
    import re
    import sys
    global idl, irec, f # Save values between calls

    [int4f,intf,flf]=check_types()
    if ix < 0 or iy < 0: 
        print ('Error in call to read_model')
        sys.exit(1)
    if filetype[0:5] == 'ascii':
        if ix != 0 or iy != 0:
            print ('Error in read_model')
            print ('Attempt to read ix,iy different from 0 from ascii file:',filename)
            print ('Requested ix,iy:',ix,iy)
            sys.exit(1)
        f=open(filename,'r')
        readstr=f.read()
        readstr=re.sub('!.*\n','\n',readstr)
        readstr=re.sub('(?i).*format version.*\n','\n',readstr)
        lines=readstr.rsplit('\n')
        for l in range(lines.count('')):
            lines.remove('')        
        f.close()
        l=lines[0].split()
        if len(lines)-1 != nz:
            print ('Error reading model file:',filename)
            print ('Incorrect number of depth points',nz,len(lines))
            sys.exit(1)
        vmac=0.
        stray_frac=0.
        expansion=0.
        if len(l) >= 1: vmac=float(l[0])
        if len(l) >= 2: stray_frac=float(l[1])
        if len(l) >= 3: expansion=float(l[2])       
        data=list()
        for icol in range(22): # Loop in model variables
            for l in lines[1:]: # First line is vmac, stray, expansion
                row=l.split()
                if filetype == 'ascii1.0':
                    if icol == 0:
                        data.append( 0. ) # z
                    elif icol == 1:
                        data.append( row[0] ) # tau
                    elif icol == 2:
                        data.append( row[1] ) # temp
                    elif icol == 3 or icol == 4:
                        data.append( row[2] ) # gas_p or rho
                    elif icol == 5:
                        data.append( row[2] ) # el_p
                    elif icol == 6:
                        data.append( row[5] ) # v_los
                    elif icol == 7:
                        data.append( row[3] ) # v_mic
                    elif icol == 8:
                        data.append( row[4] ) # b_long
                    elif icol == 9:
                        data.append( row[6] ) # b_x
                    elif icol == 10:
                        data.append( row[7] ) # b_y
                    elif icol == 11:
                        data.append( row[6] ) # bl_x
                    elif icol == 12:
                        data.append( row[7] ) # bl_y
                    elif icol == 13:
                        data.append( row[5] ) # by_z
                    elif icol == 14:
                        data.append( 0. ) # v_x
                    elif icol == 15:
                        data.append( 0. ) # v_y
                    elif icol == 16:
                        data.append( row[5] ) # v_z
                    elif icol >= 17:
                        data.append( 0. ) # nH, nH-, nH+, nH2, nH2+
                if filetype == 'ascii2.3':
                    if icol == 0:
                        data.append( row[0] ) # z
                    elif icol == 1:
                        data.append( row[1] ) # tau
                    elif icol == 2:
                        data.append( row[2] ) # temp
                    elif icol == 3:
                        data.append( row[3] ) # gas_p
                    elif icol == 4:
                        data.append( row[4] ) # rho
                    elif icol == 5:
                        data.append( row[5] ) # el_p
                    elif icol == 6:
                        data.append( row[6] ) # v_los
                    elif icol == 7:
                        data.append( row[7] ) # v_mic
                    elif icol == 8:
                        data.append( row[8] ) # b_long
                    elif icol == 9:
                        data.append( row[9] ) # b_x
                    elif icol == 10:
                        data.append( row[10] ) # b_y
                    elif icol == 11:
                        data.append( row[11] ) # bl_x
                    elif icol == 12:
                        data.append( row[12] ) # bl_y
                    elif icol == 13:
                        data.append( row[13] ) # bl_z
                    elif icol == 14:
                        data.append( row[14] ) # vl_x
                    elif icol == 15:
                        data.append( row[15] ) # vl_y
                    elif icol == 16:
                        data.append( row[16] ) # vl_z
                    elif icol >= 17:
                        data.append( 0. ) # nH, nH-, nH+, nH2, nH2+
        data.append(vmac)
        data.append(stray_frac)
        data.append(expansion)
        for i in range(8): data.append(0) # Keep el_p,gas_p,rho,nH,nH-,nH+,nH2,nH2+
        for ind in range(2): data.append(0.) # chrom_x,chrom_y
        for ind in range(92): data.append(0.) # abundances
        for i in range(len(data)): data[i]=float(data[i])
        return data
    elif filetype == 'nicole1.6':
        irec=iy+ix*ny
        f=open(filename,'rb')
        sizerec=13*nz+3 # Floats (multiply by 8 to convert to bytes)
        f.seek(sizerec*8*(irec+1)) # Skip header and previous records
        data=struct.unpack('<'+str(sizerec)+flf,f.read(sizerec*8))
        data2=range(22*nz+3+8)
        data2[0:sizerec]=data[0:sizerec]
        for ind in range(nz):
            data2[6*nz+ind]=data[6*nz+ind] # v_los
        for ind in range(nz):
            data2[8*nz+ind]=data[8*nz+ind] # b_long
        for ind in range(nz):
            data2[9*nz+ind]=data[9*nz+ind] # b_x
        for ind in range(nz):
            data2[10*nz+ind]=data[10*nz+ind] # b_y
        for ind in range(nz):
            data2[11*nz+ind]=data[9*nz+ind] # bl_x=b_x
        for ind in range(nz):
            data2[12*nz+ind]=data[10*nz+ind] # bl_y=b_y
        for ind in range(nz):
            data2[13*nz+ind]=data[8*nz+ind] # bl_z=b_long
        for ind in range(nz):
            data2[14*nz+ind]=0. # vl_x
        for ind in range(nz):
            data2[15*nz+ind]=0. # vl_y
        for ind in range(nz):
            data2[16*nz+ind]=data[6*nz+ind] # bl_y=b_y
        for ivar in range(5):
            for ind in range(nz):
                data2[(17+ivar)*nz+ind]=0. 
        data2[22*nz:22*nz+3]=data[13*nz:13*nz+3] # v_mac, stray_frac, expansion
        for ind in range(8): data2[22*nz+3+ind]=0. # keep
        for ind in range(2): data2.append(0.) # chrom_x,chrom_y
        for ind in range(92): data2.append(0.) # abundances
        f.close()
        return data2
    elif filetype == 'nicole2.3':
        if (sequential == 0):
            irec=iy+ix*ny
        else:
            irec=irec+1
        sizerec=17*nz+3 # Floats (multiply by 8 to convert to bytes)
        if (sequential == 0):
            f=open(filename,'rb')
            f.seek(sizerec*8*(irec+1)) # Skip header and previous records
        data=struct.unpack('<'+str(sizerec)+flf,f.read(sizerec*8))
        data2=list()
        for ind in range(17*nz): data2.append(data[ind])
        for ivar in range(5):
            for ind in range(nz):
                data2.append(0.)
        for ind in range(3): data2.append(data[13*nz+ind]) # v_mac, stray_frac, expansion
        for ind in range(8): data2.append(0.) # keep
        for ind in range(2): data2.append(0.) # chrom_x,chrom_y
        for ind in range(92): data2.append(0.) # abundances
        for ind in range(len(data2)): data2[ind]=float(data2[ind])
        return data2
    elif filetype == 'nicole2.6':
        if (sequential == 0):
            irec=iy+ix*ny
        else:
            irec=irec+1
        sizerec=22*nz+11+92 # Floats (multiply by 8 to convert to bytes)
        if (sequential == 0):
            f=open(filename,'rb')
            f.seek(sizerec*8*(irec+1)) # Skip header and previous records
        data=struct.unpack('<'+str(sizerec)+flf,f.read(sizerec*8))
        data=list(data)
        data.insert(22*nz+11,1.)
        data.insert(22*nz+11,-5.)
        return data
    elif filetype == 'nicole18.04':
        if (sequential == 0):
            irec=iy+ix*ny
        else:
            irec=irec+1
        sizerec=22*nz+13+92 # Floats (multiply by 8 to convert to bytes)
        if (sequential == 0):
            f=open(filename,'rb')
            f.seek(sizerec*8*(irec+1)) # Skip header and previous records
        data=struct.unpack('<'+str(sizerec)+flf,f.read(sizerec*8))
        data=list(data)
        return data
    elif filetype == 'idl':
        import idlsave
        try: 
            a=idl.z[0,0,0]
        except:
            idl=idlsave.read(filename,verbose=0)
            try:
                idl.z=idl.z.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.tau=idl.tau.reshape(nz,ny,nx)
            except: donothing=0
            idl.t=idl.t.reshape(nz,ny,nx)
            try:
                idl.gas_p=idl.gas_p.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.rho=idl.rho.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.el_p=idl.el_p.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_los=idl.v_los.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_mic=idl.v_mic.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_long=idl.b_los_z.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_tx=idl.b_los_x.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_ty=idl.b_los_y.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_x=idl.b_x.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_y=idl.b_y.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.b_z=idl.b_z.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_x=idl.v_x.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_y=idl.v_y.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_z=idl.v_z.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.nH=idl.nH.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.nHminus=idl.nHminus.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.nHplus=idl.nHplus.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.nH2=idl.nH2.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.nH2plus=idl.nH2plus.reshape(nz,ny,nx)
            except: donothing=0
            try:
                idl.v_mac=idl.v_mac.reshape(ny,nx)
            except: donothing=0
            try:
                idl.stray_frac=idl.stray_frac.reshape(ny,nx)
            except: donothing=0
            try:
                idl.expansion=idl.expansion.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_el_p=idl.keep_el_p.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_gas_p=idl.keep_gas_p.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_rho_p=idl.keep_rho_p.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_nH=idl.keep_nH.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_nHminus=idl.keep_nHminus.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_nHplus=idl.keep_nHplus.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_nH2=idl.keep_nH2.reshape(ny,nx)
            except: donothing=0
            try:
                idl.keep_nH2plus=idl.keep_nH2plus.reshape(ny,nx)
            except: donothing=0
            try:
                idl.abundance=idl.abundance.reshape(92,ny,nx)
            except: donothing=0
        data=list()
        try:
            for d in idl.z[:,iy,ix]: data.append(d)
        except: 
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.tau[:,iy,ix]: data.append(d)
        except: 
            for d in range(nz): data.append( 0. )
        for d in idl.t[:,iy,ix]: data.append(d)
        try:
            for d in idl.gas_p[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.rho[:,iy,ix]: data.append(d)
        except: 
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.el_p[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.v_los[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.v_mic[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_long[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_tx[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_ty[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_x[:,iy,ix]: data.append(d) 
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_y[:,iy,ix]: data.append(d) 
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.b_z[:,iy,ix]: data.append(d) 
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.v_x[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.v_y[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.v_z[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.nH[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.nHminus[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.nHplus[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.nH2[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            for d in idl.nH2plus[:,iy,ix]: data.append(d)
        except:
            for d in range(nz): data.append( 0. )
        try:
            data.append(idl.v_mac[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.stray_frac[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.expansion[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_el_p[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_gas_p[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_rho[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_nH[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_nHminus[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_nHplus[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_nH2[iy,ix])
        except:
            data.append( 0. )
        try:
            data.append(idl.keep_nH2plus[iy,ix])
        except:
            data.append( 0. )
#        try:
#            data.append(idl.abundance[:,iy,ix])
#        except:
#            data.append( 0. )
        for i in range(len(data)): data[i]=float(data[i])
        return data
    else:
        print ('Unknown file type')
        sys.exit(1)


def read_prof(filename, filetype, nx, ny, nlam, ix, iy, sequential=0):
# Read a particular profile from a profile file 
# Note that ix,iy=0 represents the first profile
# Use the sequential keyword to specify that this file is going to be
# read sequentially rather than accessing specific records. In that
# case, the values of ix and iy are irrelevant. The first record
# needs to be read with sequential=0 and ix,iy = 0
    import struct
    import re
    import sys
    global idl, irec, f, lastfilename # Save values between calls

    [int4f,intf,flf]=check_types()
# Read a particular profile from a profile file 
# Note that ix, iy=0 represents the first profile
    if ix < 0 or iy < 0: 
        print ('Error in call to read_prof')
        sys.exit(1)
    if filetype == 'ascii':
        if ix != 0 and iy != 0:
            print ('Error in read_prof')
            print ('Attempt to read ix,iy different from 0 from ascii file:',filename)
            print ('Requested ix,iy:',ix,iy)
            sys.exit(1)
        f=open(filename,'r')
        readstr=f.read()
        readstr=re.sub('!.*\n','\n',readstr)
        lines=readstr.rsplit('\n')
        for l in range(lines.count('')):
            lines.remove('')        
        f.close()
        if len(lines) != nlam:
            print ('Error reading profile file:',filename)
            print ('Incorrect number of wavelengths')
            sys.exit(1)
        data=list()
        for l in lines:
            row=l.split()    
            for istk in range(4):
                data.append(row[istk+1])
        for i in range(len(data)): data[i]=float(data[i])
        return data
    elif filetype[0:6] == 'nicole':
        if (sequential == 0):
            irec=iy+ix*ny
        else:
            irec=irec+1
        sizerec=4*nlam # Floats (multiply by 8 to convert to bytes)
        if (sequential == 0):
            f=open(filename,'rb')
            f.seek(sizerec*8*(irec+1)) # Skip header and previous records
        data=struct.unpack('<'+str(sizerec)+flf,f.read(sizerec*8))
        data=list(data)
        return data
    elif filetype == 'idl':
        import idlsave
        try:
            a=lastfilename
        except:
            lastfilename=''
        if filename != lastfilename:
            idl=idlsave.read(filename,verbose=0)
            idl.stki=idl.stki.reshape(nlam,ny,nx)
            idl.stkq=idl.stkq.reshape(nlam,ny,nx)
            idl.stku=idl.stku.reshape(nlam,ny,nx)
            idl.stkv=idl.stkv.reshape(nlam,ny,nx)
        lastfilename=filename
        try: 
            a=idl.stki[0,0,0]
        except:
            idl=idlsave.read(filename,verbose=0)
            idl.stki=idl.stki.reshape(nlam,ny,nx)
            idl.stkq=idl.stkq.reshape(nlam,ny,nx)
            idl.stku=idl.stku.reshape(nlam,ny,nx)
            idl.stkv=idl.stkv.reshape(nlam,ny,nx)
        data=list()
        for ilam in range(nlam): 
            data.append(idl.stki[ilam,iy,ix])
            data.append(idl.stkq[ilam,iy,ix])
            data.append(idl.stku[ilam,iy,ix])
            data.append(idl.stkv[ilam,iy,ix])
        for i in range(len(data)): data[i]=float(data[i])
        return data
    else:
        print ('Unknown file type')
        sys.exit(1)
