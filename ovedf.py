#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 12:00:00 2016

@author: ovallcorba

---------------
Changelog
Current:
 - Calculation of 2theta considering tilt/rot fit2d convention
 - Reading EDF and INP files. IntegOpts.

TODO:
 - 
"""

import numpy as np
import os, time, datetime
import argparse  
import struct
import logging as log
import matplotlib.pyplot as plt

class EDFdata:
    def __init__(self):
        self.filename = ''
        self.dimx = None    # cols
        self.dimy = None    # rows
        self.intenXY = None   # np array inten [y][x] (for matplotlib)
        self.pixSXmm = 0.079
        self.pixSYmm = 0.079
        self.distMDmm = 180
        self.wavelA = 0.4246
        self.xcen = 1024
        self.ycen = 1024
        # these data will be retrieved from the edf header
        
    def readFile(self, edfFilename):
        """Reads a EDF file and populates the data"""
        
        fedf = open(edfFilename,'r')
        self.filename = fedf.name
        size = -1
        lcount=0
        maxheaderLines=60
        try:
            for line in fedf:
                if lcount>=maxheaderLines:break
                if line.startswith("}"):break  # end of header
                
                log.info(line)
                try:
                    iigual = line.index("=")
                except:
                    lcount=lcount+1
                    continue
                    
                if line.startswith("Size"):
                    size = int(line[iigual+1:len(line)-2]) # -2 for the ;
                    log.debug("size="+str(size))
                    log.debug("typesize="+str(type(size)))                
                if line.startswith("Dim_1"):
                    self.dimx = int(line[iigual+1:len(line)-2])
                    log.debug("Dim_1="+str(self.dimx))
                if line.startswith("Dim_2"):
                    self.dimy = int(line[iigual+1:len(line)-2])
                    log.debug("Dim_2="+str(self.dimy))
                if line.startswith("beam_center_x"):
                    self.xcen = int(line[iigual+1:len(line)-2])
                    log.debug("beam_center_x="+str(self.xcen))
                if line.startswith("beam_center_y"):
                    self.ycen = int(line[iigual+1:len(line)-2])
                    log.debug("beam_center_y="+str(self.ycen))
                if line.startswith("pixelsize_x"):
                    self.pixSXmm = int(line[iigual+1:len(line)-2])/1000.
                    log.debug("pixelsize_x="+str(self.pixSXmm))
                if line.startswith("pixelsize_y"):
                    self.pixSYmm = int(line[iigual+1:len(line)-2])/1000.
                    log.debug("pixelsize_y="+str(self.pixSYmm))
                if line.startswith("ref_distance"):
                    self.distMDmm = int(line[iigual+1:len(line)-2])
                    log.debug("ref_distance="+str(self.distMDmm))
                if line.startswith("ref_wave"):
                    self.wavelA = int(line[iigual+1:len(line)-2])
                    log.debug("ref_wave="+str(self.wavelA))
                
                lcount=lcount+1
        finally:
            log.debug("finally executed")
            fedf.close()
        
        # now the binary part
        filesize = os.path.getsize(edfFilename)
        log.info("filesize="+str(filesize))
        headersize = filesize - size
        log.info("headersize="+str(headersize))
        
        fedf = open(edfFilename, "rb")
        self.intenXY = np.empty([self.dimx,self.dimy])
        row = []
        header = fedf.read(headersize)
        log.info(header)
        try:
            for y in range(self.dimy):
                for x in range(self.dimx):
                    # edf two bytes unsigned short LowByteFirst (uint8 little endian)
                    self.intenXY[y][x] =  struct.unpack("<H",fedf.read(2))[0]
                log.debug(self.intenXY[x])
        finally:
            fedf.close()

class IntegOptions:
    def __init__(self):
        self.xcen = -1.0 # if set this will override image header info
        self.ycen = -1.0
        self.distMDmm = -1.0
        self.wavelA = -1.0
        self.tiltRotation = 0.0
        self.angleTilt = 0.0
        self.subadu = 0.0
        self.startAzim = 0
        self.endAzim = 360
        self.inRadi = 0
        self.outRadi = 1000
        self.xbin = 1
        self.ybin = 1
        self.azimBins = 1
        self.radialBins = 1000
        # TODO:MASK implementation
        
        # this inside here to speed up
        self.costilt = None
        self.sintilt = None
        self.cosrot = None
        self.sinrot = None

        self.dist = None
        self.x0 = None
        self.phi = None
        self.tilt = None

    def readINPFile(self, inpFilename):
        """Reads an INP file with calibration info"""
        finp = open(inpFilename,'r')
        try:
            for line in finp:
                log.info(line)
                try:
                    iigual = line.index("=")
                except:
                    continue
                if line.startswith("X-BEAM"):
                    self.xcen = float(line[iigual+1:len(line)-3]) 
                    log.debug("X-BEAM="+str(self.xcen))
                if line.startswith("Y-BEAM"):
                    self.ycen = float(line[iigual+1:len(line)-3])
                    log.debug("Y-BEAM="+str(self.ycen))
                if line.startswith("DISTANCE"):
                    self.distMDmm = float(line[iigual+1:len(line)-3])
                    log.debug("DISTANCE="+str(self.distMDmm))
                if line.startswith("WAVELE"):
                    self.wavelA = float(line[iigual+1:len(line)-3])
                    log.debug("WAVELENGTH="+str(self.wavelA))
                if line.startswith("TILT"):
                    self.tiltRotation = float(line[iigual+1:len(line)-3])
                    log.debug("TILT ROTATION="+str(self.tiltRotation))
                if line.startswith("ANGLE"):
                    self.angleTilt = float(line[iigual+1:len(line)-3])
                    log.debug("ANGLE OF TILT="+str(self.angleTilt))
                if line.startswith("SUBADU"):
                    self.subadu = float(line[iigual+1:len(line)-3])
                    log.debug("SUBADU="+str(self.subadu))
                if line.startswith("START"):
                    self.startAzim = float(line[iigual+1:len(line)-3])
                    log.debug("START AZIM="+str(self.startAzim))                 
                if line.startswith("END"):
                    self.endAzim = float(line[iigual+1:len(line)-3])
                    log.debug("END AZIM="+str(self.endAzim))
                if line.startswith("INNER"):
                    self.inRadi = float(line[iigual+1:len(line)-3])
                    log.debug("INNER RADIUS="+str(self.inRadi))                 
                if line.startswith("OUTER"):
                    self.outRadi = float(line[iigual+1:len(line)-3])
                    log.debug("OUTER RADIUS="+str(self.outRadi))

                # TODO read xbin ybin azimbins radialbins

        finally:
            log.debug("finally executed")
            finp.close()
        
        # calculation here to speed up the process
        self.costilt = math.cos(math.radians(integOpts.angleTilt))
        self.sintilt = math.sin(math.radians(integOpts.angleTilt))
        self.cosrot = math.cos(math.radians(integOpts.tiltRotation))
        self.sinrot = math.sin(math.radians(integOpts.tiltRotation))

        self.rotM_rotZ = np.array(([self.cosrot,-self.sinrot,0.0],[self.sinrot,self.cosrot,0.0],[0.0,0.0,1.0]),dtype=np.float32)
        self.rotM_tiltX = np.array(([1.0,1.0,0.0],[0.0,self.costilt,-self.sintilt],[0.0,self.sintilt,self.costilt]),dtype=np.float32)
        
class D1Pattern:
    def __init__(self):
        self.t2 = []     # array 2Theta
        self.Iobs = []   # array Yobs
        self.esd = []    # array esd
        self.npix = []   # array contributing pixels
        
######################################################## GENERAL FUNCTIONS

def calc2tDeg(integOpts, edfdata, rowY,colX):

    costilt = integOpts.costilt
    sintilt = integOpts.sintilt
    cosrot = integOpts.cosrot
    sinrot = integOpts.sinrot
    
    # vector centre-pixel
    vPCx=(float)(colX)-integOpts.xcen
    vPCy=integOpts.ycen-(float)(rowY)
    
    sum1 = ((vPCx*cosrot + vPCy*sinrot)*(vPCx*cosrot + vPCy*sinrot))+((vPCy*cosrot - vPCx*sinrot)*(vPCy*cosrot - vPCx*sinrot))
    num = costilt*costilt*sum1
    den = (integOpts.distMDmm/edfdata.pixSXmm) + sintilt*(vPCx*cosrot+vPCy*sinrot)
    t2p = math.atan(math.sqrt(num/(den*den)))
    
    return math.degrees(t2p)



def calc2tDegF2D(integOpts, edfdata, rowY,colX):
    
    costilt = integOpts.costilt
    sintilt = integOpts.sintilt
    cosrot = integOpts.cosrot
    sinrot = integOpts.sinrot
    
    # vector centre-pixel
    vPCx=colX-integOpts.xcen;
    vPCy=integOpts.ycen-rowY;
    
    num = costilt*costilt*((vPCx*cosrot+vPCy*sinrot)*(vPCx*cosrot+vPCy*sinrot))+((vPCy*cosrot-vPCx*sinrot)*(vPCy*cosrot-vPCx*sinrot))
    den = (integOpts.distMDmm/edfdata.pixSXmm) + sintilt*(vPCx*cosrot+vPCy*sinrot)
    t2p = math.atan(math.sqrt(num/(den*den)))
    
    return math.degrees(t2p)


def calc2tGS(integOpts, edfdata, rowY,colX):
    
    costilt = integOpts.costilt
    sintilt = integOpts.sintilt
    cosrot = integOpts.cosrot
    sinrot = integOpts.sinrot

    # vector centre-pixel
    vPCx=colX-integOpts.xcen
    vPCy=integOpts.ycen-rowY
    
    vPCx = np.array(colX-integOpts.xcen,dtype=np.float32)
    vPCy = np.array(integOpts.ycen-rowY,dtype=np.float32)
    X = np.array(([vPCx,vPCy,0.0]),dtype=np.float32).T
    X = np.dot(X,integOpts.rotM_rotZ)
    Z = np.dot(X,integOpts.rotM_tiltX).T[2]
    tth = math.atan(np.sqrt(vPCx**2+vPCy**2-Z**2)/(integOpts.distT-Z))
    
    return math.degrees(tth)


def makeMat(Angle,Axis):
    '''Make rotation matrix from Angle and Axis

    :param float Angle: in degrees
    :param int Axis: 0 for rotation about x, 1 for about y, etc.
    '''
    cs = math.cos(Angle)
    ss = math.sin(Angle)
    M = np.array(([1.,0.,0.],[0.,cs,-ss],[0.,ss,cs]),dtype=np.float32)
    return np.roll(np.roll(M,Axis,axis=0),Axis,axis=1)

######################################################## MAIN PROGRAM
if __name__=="__main__":

    __version__ = '161017'

    ts = time.time()
    st = '[' + datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S') + '] '
    welcomeMsg = st+'=== OV EDF READING ver.'+__version__+' ==='

    parser = argparse.ArgumentParser(description="EDF file reading and integration")
    parser.add_argument('filename', nargs=1, help='edf file')

    parser.add_argument("-i", "--input", default=None, nargs=1, help='calibration input file')
    parser.add_argument("-o", "--output", default=None, nargs=1, help='1D output file')
    parser.add_argument("-pl", "--plot", action="store_true", help='plot image')
    parser.add_argument("-d", "--debug", action="count", help='debug mode')    

    args = parser.parse_args()
    
    print args.debug
    if args.debug==0:
        log.basicConfig(level=log.CRITICAL)
    if args.debug==1:
        log.basicConfig(level=log.INFO)
    if args.debug==2:
        log.basicConfig(level=log.DEBUG)

    log.info('parser.parse_args()= %s',args)
    
    print ''
    print welcomeMsg
    print ''
    
    log.info("Input file: %s"%args.filename[0])
    print "  Input file: ",args.filename[0]
    edf = EDFdata()
    edf.readFile(args.filename[0])

    # IN CASE WE WANT TO PLOT THE IMAGE
    if (args.plot):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.imshow(edf.intenXY)
        plt.colorbar()
        numrows, numcols = edf.intenXY.shape
        def format_coord(x, y):
            col = int(x+0.5)
            row = int(y+0.5)
            if col>=0 and col<numcols and row>=0 and row<numrows:
                z = edf.intenXY[row,col]
                return 'col_x=%.1f, row_y=%.1f, Intensity=%d'%(x, y, z)
            else:
                return 'col_x=%.1f, row_y=%.1f'%(x, y)
        ax.format_coord = format_coord
        plt.show()
    
    # NOW READ THE OPTIONS
    integOpts = IntegOptions()
    if args.input is not None:
        log.info("reading inp file:"+args.input[0])
        integOpts.readINPFile(args.input[0])
        
        if integOpts.xcen <= 0:
            integOpts.xcen = edf.xcen
        if integOpts.ycen <= 0:
            integOpts.ycen = edf.ycen
        if integOpts.distMDmm <= 0:
            integOpts.distMDmm = edf.distMDmm
        if integOpts.wavelA <= 0:
            integOpts.wavelA = edf.wavelA
        
        print ""
        print "Integration options"
        print "Center X,Y = %.3f %.3f"%(integOpts.xcen,integOpts.ycen)
        print "Distance Wavelengh = %.3f %.4f"%(integOpts.distMDmm,integOpts.wavelA)
        print "tiltRot AngTilt = %.3f %.3f"%(integOpts.tiltRotation,integOpts.angleTilt)
        print ""

    
    # NOW WE SHOULD INTEGRATE
    
    # FAST TEST to check
    step = 0.023
    t2in = 1
    t2fin = 23
    
    size = (int)((t2fin-t2in)/step +2)
    patt = D1Pattern()
    for i in range(size):
        #log.info(i)
        patt.t2.append(t2in + i*step)
        patt.esd.append(0)
        patt.Iobs.append(0)
        patt.npix.append(0)
    
    integOpts.distT = (integOpts.distMDmm/edf.pixSXmm)/integOpts.costilt

    for y in range(edf.dimy):
        for x in range(edf.dimx):
            #TODO: check for excluded zones
            t2p = calc2tGS(integOpts, edf, y,x)
            if (t2p>t2in) and (t2p<t2fin):
                #position to the vector
                p = (int)(t2p/step + 0.5)-(int)(t2in/step + 0.5)
                
                #TODO: CORRECTIONS!!
                inten = edf.intenXY[y][x]
                patt.Iobs[p] = patt.Iobs[p] + inten - (int)(integOpts.subadu)
                patt.npix[p] = patt.npix[p] + 1
        #print "line %d"%y
    
    #write a xy file
    fout = open("test.xy", 'w')
    fout.write("#TEST PATTERN \n")
    for i in range(size):
        if patt.npix[i] <= 0:continue
        fout.write(" %.5f %5f \n"%(patt.t2[i], (patt.Iobs[i]/patt.npix[i])))
    fout.close()
    log.info("xy file written")
