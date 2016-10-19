#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 12:00:00 2016

@author: ovallcorba

---------------
Changelog
Current (last change 161019 19.15h):
 - Azimuth range, bins
 - Added PD lorentz correction (at the writting of the dat file)
 - Speed up calculation by numpy operations
 - Calculation of 2theta considering tilt/rot fit2d convention
 - Reading EDF and INP files. IntegOpts.

TODO:
 - xbin, ybin, azimbins, masks
 - geometrical corrections needed?
 - esd
 - support for multiple files (batch). It is easy but for testing I 
   prefer to keep 1 argument (image) only
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
        
        self.t2XY = None   #np array of t2
        self.azim = None
        self.lorfac = None  #np array of lorentz factors

    def readFile(self, edfFilename):
        """Reads a EDF file and populates the data"""
        t0 = time.time()        
        fedf = open(edfFilename,'r')
        self.filename = fedf.name
        size = -1
        lcount=0
        maxheaderLines=60
        try:
            for line in fedf:
                if lcount>=maxheaderLines:break
                if line.startswith("}"):break  # end of header
                
                log.debug(line)
                try:
                    iigual = line.index("=")
                    ipcoma = line.index(";")
                except:
                    lcount=lcount+1
                    continue
                    
                if line.startswith("Size"):
                    size = int(line[iigual+1:ipcoma].strip())
                    log.debug("size="+str(size))
                    log.debug("typesize="+str(type(size)))                
                if line.startswith("Dim_1"):
                    self.dimx = int(line[iigual+1:ipcoma].strip())
                    log.debug("Dim_1="+str(self.dimx))
                if line.startswith("Dim_2"):
                    self.dimy = int(line[iigual+1:ipcoma].strip())
                    log.debug("Dim_2="+str(self.dimy))
                if line.startswith("beam_center_x"):
                    self.xcen = float(line[iigual+1:ipcoma].strip())
                    log.debug("beam_center_x="+str(self.xcen))
                if line.startswith("beam_center_y"):
                    self.ycen = float(line[iigual+1:ipcoma].strip())
                    log.debug("beam_center_y="+str(self.ycen))
                if line.startswith("pixelsize_x"):
                    self.pixSXmm = float(line[iigual+1:ipcoma].strip())/1000.
                    log.debug("pixelsize_x="+str(self.pixSXmm))
                if line.startswith("pixelsize_y"):
                    self.pixSYmm = float(line[iigual+1:ipcoma].strip())/1000.
                    log.debug("pixelsize_y="+str(self.pixSYmm))
                if line.startswith("ref_distance"):
                    self.distMDmm = float(line[iigual+1:ipcoma].strip())
                    log.debug("ref_distance="+str(self.distMDmm))
                if line.startswith("ref_wave"):
                    self.wavelA = float(line[iigual+1:ipcoma].strip())
                    log.debug("ref_wave="+str(self.wavelA))
                lcount=lcount+1
        finally:
            log.debug("finally executed")
            fedf.close()
        
        # now the binary part
        filesize = os.path.getsize(edfFilename)
        log.debug("filesize="+str(filesize))
        headersize = filesize - size
        log.debug("headersize="+str(headersize))

        #init
        self.intenXY = np.empty([self.dimx,self.dimy])
        self.t2XY = np.empty([self.dimx,self.dimy])
        self.azim = np.empty([self.dimx,self.dimy])
        self.lorfac = np.empty([self.dimx,self.dimy])
                
        fedf = open(edfFilename, "rb")
        row = []
        header = fedf.read(headersize)
        log.debug(header)
        try:
            for y in range(self.dimy):
                for x in range(self.dimx):
                    self.intenXY[y,x] =  struct.unpack("<H",fedf.read(2))[0] #edf I as unsigned short little endian
                log.debug(self.intenXY[x])
        finally:
            fedf.close()
        log.info(" EDF reading: %.4f sec"%(time.time() - t0))

class IntegOptions:
    def __init__(self):
        #if set these parameters will override image header info
        self.xcen = -1.0   
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
                log.debug(line)
                try:
                    iigual = line.index("=")
                    ipcoma = line.index(";")
                except:
                    continue
                if line.startswith("X-BEAM"):
                    self.xcen = float(line[iigual+1:ipcoma].strip()) 
                    log.debug("X-BEAM="+str(self.xcen))
                if line.startswith("Y-BEAM"):
                    self.ycen = float(line[iigual+1:ipcoma].strip())
                    log.debug("Y-BEAM="+str(self.ycen))
                if line.startswith("DISTANCE"):
                    self.distMDmm = float(line[iigual+1:ipcoma].strip())
                    log.debug("DISTANCE="+str(self.distMDmm))
                if line.startswith("WAVELE"):
                    self.wavelA = float(line[iigual+1:ipcoma].strip())
                    log.debug("WAVELENGTH="+str(self.wavelA))
                if line.startswith("TILT"):
                    self.tiltRotation = float(line[iigual+1:ipcoma].strip())
                    log.debug("TILT ROTATION="+str(self.tiltRotation).strip())
                if line.startswith("ANGLE"):
                    self.angleTilt = float(line[iigual+1:ipcoma].strip())
                    log.debug("ANGLE OF TILT="+str(self.angleTilt).strip())
                if line.startswith("SUBADU"):
                    self.subadu = float(line[iigual+1:ipcoma].strip())
                    log.debug("SUBADU="+str(self.subadu))
                if line.startswith("XBIN"):
                    self.xbin = float(line[iigual+1:ipcoma].strip())
                    log.debug("XBIN="+str(self.xbin))
                if line.startswith("YBIN"):
                    self.ybin = float(line[iigual+1:ipcoma].strip())
                    log.debug("YBIN="+str(self.ybin))
                if line.startswith("START"):
                    self.startAzim = float(line[iigual+1:ipcoma].strip())
                    log.debug("START AZIM="+str(self.startAzim))                 
                if line.startswith("END"):
                    self.endAzim = float(line[iigual+1:ipcoma].strip())
                    log.debug("END AZIM="+str(self.endAzim))
                if line.startswith("INNER"):
                    self.inRadi = float(line[iigual+1:ipcoma].strip())
                    log.debug("INNER RADIUS="+str(self.inRadi))                 
                if line.startswith("OUTER"):
                    self.outRadi = float(line[iigual+1:ipcoma].strip())
                    log.debug("OUTER RADIUS="+str(self.outRadi))
                if line.startswith("AZIMUTH BINS"):
                    self.azimBins = float(line[iigual+1:ipcoma].strip())
                    log.debug("AZIMUTH BINS="+str(self.azimBins))
                if line.startswith("RADIAL BINS"):
                    self.radialBins = float(line[iigual+1:ipcoma].strip())
                    log.debug("RADIAL BINS="+str(self.radialBins))
                #TODO MASK

        finally:
            log.debug("finally executed")
            finp.close()
        
        #calculation here to speed up the process
        self.costilt = math.cos(math.radians(self.angleTilt))
        self.sintilt = math.sin(math.radians(self.angleTilt))
        self.cosrot = math.cos(math.radians(self.tiltRotation))
        self.sinrot = math.sin(math.radians(self.tiltRotation))
        
class D1Pattern:
    def __init__(self):
        self.t2 = []     # array 2Theta
        self.Iobs = []   # array Yobs
        self.esd = []    # array esd
        self.npix = []   # array contributing pixels
        
######################################################## GENERAL FUNCTIONS

def calc2tDeg(integOpts, edfdata, rowY,colX):
    """single pixel calculation"""
    costilt = integOpts.costilt
    sintilt = integOpts.sintilt
    cosrot = integOpts.cosrot
    sinrot = integOpts.sinrot
    distPix = (integOpts.distMDmm/edfdata.pixSXmm)/integOpts.costilt
    
    # vector centre-pixel
    vPCx=(float)(colX)-integOpts.xcen
    vPCy=integOpts.ycen-(float)(rowY)
    
    zcomponent = (-vPCx*sinrot + vPCy*cosrot)*(-sintilt)
    tiltedVecMod = math.sqrt(vPCx**2+vPCy**2-zcomponent**2)
    t2p = math.atan(tiltedVecMod/(distPix-zcomponent))
    
    return math.degrees(t2p)

def getAzim(integOpts,edfdata,rowY,colX):
    if rowY == int(edfdata.ycen+0.5) and colX == int(edfdata.xcen+0.5): return 0
    
    #vector centre-pixel
    vCPx=(float)(colX)-integOpts.xcen
    vCPy=integOpts.ycen-(float)(rowY)
    
    #angle --- NOW is defined FROM 12h clockwise +
    azim = np.atan2(vCPx,vCPy)
    
    if (azim<0): azim = azim + 2*math.pi
    return math.degrees(azim)



def calc2tDegFull(integOpts, edfdata): #and azimuth
    """full image calculation (to speed up using numpy operations)"""
    #first fill up useful arrays (e.g. vPCx**2+vPCy**2-zcomponent**2)
    t0 = time.time()
    #inits
    vCPx = np.empty([edfdata.dimx,edfdata.dimy])
    vCPy = np.empty([edfdata.dimx,edfdata.dimy])
    vCPz = np.empty([edfdata.dimx,edfdata.dimy])
    distPix = (integOpts.distMDmm/edfdata.pixSXmm)/integOpts.costilt    
    #horizontal vector (for the azimuth)
    horX = 1.0
    horY = 0.0
    for y in range(edfdata.dimy):
        cpy = integOpts.ycen-(float)(y)
        for x in range(edfdata.dimx):
            cpx = (float)(x)-integOpts.xcen
            vCPx[y,x] = cpx
            vCPy[y,x] = cpy
            vCPz[y,x] = (-cpx*integOpts.sinrot + cpy*integOpts.cosrot)*(-integOpts.sintilt)
    log.info(" vector calculation: %.4f sec"%(time.time() - t0))
    
    #now the square root
    t0 = time.time()
    edfdata.t2XY = vCPx**2 + vCPy**2 - vCPz**2
    edfdata.t2XY = np.sqrt(edfdata.t2XY)
    log.info(" sqrt calculation: %.4f sec"%(time.time() - t0))
    
    #now the arctan
    t0 = time.time()
    edfdata.t2XY = edfdata.t2XY/ (distPix-vCPz)
    edfdata.t2XY = np.arctan(edfdata.t2XY)
    log.info(" atan calculation: %.4f sec"%(time.time() - t0))
    
    #TODO: GEOM CORRECTIONS (Lorentz,...)
    t0 = time.time()
    #edfdata.lorfac = 1/((np.sin(edfdata.t2XY))*((np.sin(edfdata.t2XY/2))))
    #edfdata.lorfac = 1/((np.cos(edfdata.t2XY/2))*((np.sin(edfdata.t2XY/2))**2))
    #edfdata.lorfac = ((np.sin(edfdata.t2XY))**2)/((np.cos(edfdata.t2XY/2)))
    #edfdata.lorfac = 1/((np.cos(edfdata.t2XY))**2)
    log.info(" lorfac calculation: %.4f sec"%(time.time() - t0))
    
    #t2 to degrees
    t0 = time.time()
    edfdata.t2XY = np.degrees(edfdata.t2XY)
    log.info(" to deg calculation: %.4f sec"%(time.time() - t0))
    
    #Azim calc
    t0 = time.time()
    edfdata.azim = np.degrees(np.arctan2(vCPx,vCPy))
    log.debug(" (y,x)=1200,1200 azim = %.2f deg"%(edfdata.azim.item(1200,1200)))
    edfdata.azim=np.where(edfdata.azim < 0, edfdata.azim+360,edfdata.azim)
    log.info(" azim calculation: %.4f sec"%(time.time() - t0))
    log.debug(" (y,x)=1200,1200 azim = %.2f deg"%(edfdata.azim.item(1200,1200)))
    log.debug(" (y,x)=1200,800 azim = %.2f deg"%(edfdata.azim.item(1200,800)))
    log.debug(" (y,x)=800,1200 azim = %.2f deg"%(edfdata.azim.item(800,1200)))
    log.debug(" (y,x)=800,800 azim = %.2f deg"%(edfdata.azim.item(800,800)))


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

    
    # NOW WE INTEGRATE
    
    #Integration options
    t2in = calc2tDeg(integOpts, edf, edf.ycen,edf.xcen+integOpts.inRadi)
    t2fin = calc2tDeg(integOpts, edf, edf.ycen,edf.xcen+integOpts.outRadi)
    #TODO azimbins
    size = (int)(integOpts.radialBins)
    step = (t2fin-t2in)/size
    
    #put start and end azim in range 0-360 starting at 12o'clock direction clockwise !!
    #Inserted as fit2d 3o'clock=zero and positive clockwise (because it is y inverted!)
    startAzim = integOpts.startAzim
    if (startAzim < 0): startAzim = startAzim + 360
    endAzim = integOpts.endAzim
    if (endAzim < 0): endAzim = endAzim + 360
    #now +90 to put at 12 (remember + clockwise)
    startAzim = startAzim +90
    endAzim = endAzim +90

    log.info("t2in=%f t2f=%f size=%d step=%f"%(t2in,t2fin,size,step))
    
    print ""
    print "Integration parameters and options:"
    print " Center X,Y = %.3f %.3f"%(integOpts.xcen,integOpts.ycen)
    print " Distance Wavelengh = %.3f %.4f"%(integOpts.distMDmm,integOpts.wavelA)
    print " tiltRot AngTilt = %.3f %.3f"%(integOpts.tiltRotation,integOpts.angleTilt)
    print " t2ini t2end step = %.4f %.4f %.4f"%(t2in,t2fin,step)
    print " startAzim endAzim = %.4f %.4f (ref. 12h CW+ %.4f %.4f)"%(integOpts.startAzim,integOpts.endAzim,startAzim,endAzim)
    print ""    
        
    patt = D1Pattern()
    for i in range(size+1):
        #log.info(i)
        patt.t2.append(t2in + i*step)
        patt.esd.append(0)
        patt.Iobs.append(0)
        patt.npix.append(0)
    
    calc2tDegFull(integOpts,edf)

    #IN CASE WE WANT TO PLOT THE IMAGE (I moved it after calculations)
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
                inten = edf.intenXY.item(row,col)
                #t2 = calc2tDeg(integOpts, edf, row,col)
                #azim = getAzim(integOpts,edf,row,col)
                t2 = edf.t2XY.item(row,col)
                azim = edf.azim.item(row,col)
                return 'col_x=%.1f, row_y=%.1f, Intensity=%d, t2=%.3f, azim=%.2f'%(x, y, inten,t2,azim)
            else:
                return 'col_x=%.1f, row_y=%.1f'%(x, y)
        ax.format_coord = format_coord
        plt.show()

    t0 = time.time()

    #TODO: APPLY GEOM CORRECTIONS HERE
    #edf.intenXY = edf.intenXY*edf.lorfac
    #log.info(" apply lor corr: %.4f sec"%(time.time() - t0))
    #t0 = time.time()
    
    for y in range(edf.dimy):
        for x in range(edf.dimx):
            #TODO: check for excluded zones
 
            azim = edf.azim.item(y,x)
            if (endAzim<startAzim):
                if (azim < startAzim) and (azim > endAzim): continue
            else:
                if (azim < startAzim) or (azim > endAzim): continue
            
            t2p = edf.t2XY.item(y,x)
            if (t2p<t2in):continue
            if (t2p>t2fin):continue

            #position in the histogram
            p = (int)(t2p/step + 0.5)-(int)(t2in/step + 0.5)
            
            inten = edf.intenXY.item(y,x)
            patt.Iobs[p] = patt.Iobs[p] + inten + integOpts.subadu
            patt.npix[p] = patt.npix[p] + 1
    
    log.info(" histogram filling: %.4f sec"%(time.time() - t0))
    
    #write a dat file
    #TODO: HEADER
    fout = open("test.dat", 'w')
    fout.write("#TEST PATTERN \n")
    for i in range(size+1):
        if patt.npix[i] <= 0:continue
        lorfact = (1./((math.cos(math.radians(patt.t2[i])))**2))
        icorr = patt.Iobs[i] * lorfact
        if (icorr<0): icorr = 0
        fout.write(" %.5f %5f %5f\n"%(patt.t2[i], icorr/patt.npix[i], math.sqrt(icorr)/patt.npix[i]))
    fout.close()
    log.info("dat file written")
