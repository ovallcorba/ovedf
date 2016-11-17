#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 12:00:00 2016

@author: ovallcorba

---------------
Changelog
Current (last change 161117 17.30h):
 - Mask in BIN format (image with all I=0 except mask pixels with I=-1)
 - D2Dplot or Fit2d convention (argument option)
 - Added counters on header
 - Speed optimization, only histogram filling is pending
 - direct calculation of ESD using npix
 - Azimuthal
 - Support for multiple files (batch)
 - Header as in ffops generated files
 - ESD (as in ffops)
 - Azimuth range, bins
 - Added PD lorentz correction (at the writting of the dat file)
 - Speed up calculation by numpy operations
 - Calculation of 2theta considering tilt/rot fit2d convention
 - Reading EDF and INP files. IntegOpts.

TODO:
 - Speed up histogram filling
 - maybe apply mask before to the whole image? would it be faster?

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
        self.comment = ''
        self.counter_mne = []
        self.counter_pos = []
        # these data will be retrieved from the edf header
        
        self.t2XY = None   #np array of t2
        self.azim = None   #np array of azim

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
                if line.strip().startswith("}"):break  #end of header

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
                if line.startswith("Comment"):
                    self.comment = line[iigual+1:ipcoma].strip()
                    log.debug("Comment="+self.comment)
                if line.startswith("counter_mne"):
                    self.counter_mne = line[iigual+1:ipcoma].strip().split()
                    log.debug("counter_mne="+str(self.counter_mne))
                if line.startswith("counter_pos"):
                    self.counter_pos = line[iigual+1:ipcoma].strip().split()
                    log.debug("counter_pos="+str(self.counter_pos))
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

        fedf = open(edfFilename, "rb")
        header = fedf.read(headersize)
        log.debug(header)
        try:
            self.intenXY = np.fromfile(fedf,np.uint16)
            self.intenXY = self.intenXY.reshape(self.dimx,self.dimy)
        finally:
            fedf.close()
        log.info(" EDF reading: %.4f sec"%(time.time() - t0))

    def readBINFile(self, binFilename):
        """Reads a BIN file and populates the data"""
        t0 = time.time()

        #binary 60 bytes Header, then intensity as signed int*2 as (1,1), (2,1), (3,1),....
        #header contains:
        # Int*4 NXMX(cols)
        # Int*4 NYMX(rows)
        # Real*4 SCALE
        # Real*4 CENTX
        # Real*4 CENTY
        # Real*4 PIXLX
        # Real*4 PIXLY
        # Real*4 DISTOD
        # Real*4 WAVEL
        # Real*4 OME/PHI ini (degrees)
        # Real*4 OME/PHI final (degrees)
        # REAL*4 ACQTIME
        headersize = 60

                
        fbin = open(binFilename, "rb")
        
        self.dimx = struct.unpack("<i",fbin.read(4))[0]
        self.dimy = struct.unpack("<i",fbin.read(4))[0]
        fbin.read(4) #scale
        self.xcen = struct.unpack("<f",fbin.read(4))[0]
        self.ycen = struct.unpack("<f",fbin.read(4))[0]
        self.pixSXmm = struct.unpack("<f",fbin.read(4))[0]
        self.pixSYmm = struct.unpack("<f",fbin.read(4))[0]
        self.distMDmm = struct.unpack("<f",fbin.read(4))[0]
        self.wavelA = struct.unpack("<f",fbin.read(4))[0]
        #up to here 36 bytes readed, 24 remaining for the header
        fbin.read(24)
        
        log.info("dimx=%d dimy=%d xcen=%.2f ycen=%.2f pixSXmm=%.2f pixSYmm=%.2f distMDmm=%.2f wavelA=%.4f"%(self.dimx,self.dimy,self.xcen,self.ycen,self.pixSXmm,self.pixSYmm,self.distMDmm,self.wavelA))
        
        #init
        self.intenXY = np.empty([self.dimx,self.dimy])

        try:
            self.intenXY = np.fromfile(fbin,np.int16)
            self.intenXY = self.intenXY.reshape(self.dimx,self.dimy)
        finally:
            fbin.close()
        log.info(" BIN reading: %.4f sec"%(time.time() - t0))

class IntegOptions:
    def __init__(self):
        #if set these parameters will override image header info
        self.xcen = -1.0   
        self.ycen = -1.0
        self.distMDmm = -1.0
        self.wavelA = -1.0
        self.tiltRotation = 0.0
        self.angleTilt = 0.0
        self.tiltRotationIN = 0.0
        self.subadu = 0.0
        self.startAzim = 0
        self.endAzim = 360
        self.inRadi = 0
        self.outRadi = 1000
        self.xbin = 1
        self.ybin = 1
        self.azimBins = 1
        self.radialBins = 1000
        self.fit2d = False
        self.maskf = None
        
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
                if line.strip().startswith("#"): continue
                try:
                    iigual = line.index("=")
                except:
                    continue
                try:
                    ipcoma = line.index(";")
                except:
                    ipcoma = len(line)-1

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
                    self.tiltRotationIN = float(line[iigual+1:ipcoma].strip())
                    log.debug("TILT ROTATION="+str(self.tiltRotationIN).strip())
                if line.startswith("ANGLE"):
                    self.angleTilt = float(line[iigual+1:ipcoma].strip())
                    log.debug("ANGLE OF TILT="+str(self.angleTilt).strip())
                if line.startswith("SUBADU"):
                    self.subadu = float(line[iigual+1:ipcoma].strip())
                    log.debug("SUBADU="+str(self.subadu))
                if line.startswith("XBIN"):
                    self.xbin = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("XBIN="+str(self.xbin))
                if line.startswith("YBIN"):
                    self.ybin = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("YBIN="+str(self.ybin))
                if line.startswith("START"):
                    self.startAzim = float(line[iigual+1:ipcoma].strip())
                    log.debug("START AZIM="+str(self.startAzim))                 
                if line.startswith("END"):
                    self.endAzim = float(line[iigual+1:ipcoma].strip())
                    log.debug("END AZIM="+str(self.endAzim))
                if line.startswith("INNER"):
                    self.inRadi = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("INNER RADIUS="+str(self.inRadi))                 
                if line.startswith("OUTER"):
                    self.outRadi = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("OUTER RADIUS="+str(self.outRadi))
                if line.startswith("AZIMUTH BINS"):
                    self.azimBins = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("AZIMUTH BINS="+str(self.azimBins))
                if line.startswith("RADIAL BINS"):
                    self.radialBins = int(float(line[iigual+1:ipcoma].strip()))
                    log.debug("RADIAL BINS="+str(self.radialBins))
                if line.startswith("MASK"):
                    self.maskf = line[iigual+1:ipcoma].strip()
                    log.debug("MASK="+self.maskf)

        finally:
            log.debug("finally executed")
            finp.close()

        if self.fit2d:
            #rot convention from fit2d (0 at 3 CW+ when corrected flip H) to 0 at 12h CW+
            self.tiltRotation = self.tiltRotationIN + 90
            if (self.tiltRotation>=180):
                self.tiltRotation = self.tiltRotationIN%(-180)
            self.angleTilt = -self.angleTilt
        else:
            #direct value
            self.tiltRotation=self.tiltRotationIN

        #calculation here to speed up the process
        self.costilt = math.cos(math.radians(self.angleTilt))
        self.sintilt = math.sin(math.radians(self.angleTilt))
        self.cosrot = math.cos(math.radians(self.tiltRotation))
        self.sinrot = math.sin(math.radians(self.tiltRotation))
        
        log.debug('tilt=%.6f rot=%.6f sT=%.6f cT=%.6f sR=%.6f cR=%.6f'%(self.tiltRotation,self.angleTilt,self.sintilt,self.costilt,self.sinrot,self.cosrot))

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
    
    zcomponent = (vPCx*sinrot + vPCy*cosrot)*(sintilt)
    tiltedVecMod = math.sqrt(vPCx**2+vPCy**2-zcomponent**2)
    t2p = math.atan(tiltedVecMod/(distPix-zcomponent))
    
    return math.degrees(t2p)

def getAzim(integOpts,edfdata,rowY,colX):
    if rowY == int(edfdata.ycen+0.5) and colX == int(edfdata.xcen+0.5): return 0
    
    #vector centre-pixel
    vCPx=(float)(colX)-integOpts.xcen
    vCPy=integOpts.ycen-(float)(rowY)
    
    #angle --- defined FROM 12h clockwise +
    azim = np.arctan2(vCPx,vCPy)
    
    if (azim<0): azim = azim + 2*math.pi
    return math.degrees(azim)

def calc2tDegFull(integOpts, edfdata): #and azimuth
    """full image calculation (to speed up using numpy operations)"""
    t0 = time.time()
    #inits
    vCPy,vCPx = np.mgrid[0:edfdata.dimy,0:edfdata.dimx]
    vCPz = np.empty([edfdata.dimx,edfdata.dimy])
    distPix = (integOpts.distMDmm/edfdata.pixSXmm)/integOpts.costilt
    
    #vectors calculation
    vCPx = vCPx - integOpts.xcen
    vCPy = integOpts.ycen - vCPy
    vCPz = (vCPx*integOpts.sinrot + vCPy*integOpts.cosrot)*(integOpts.sintilt)
    log.info(" vector calculation: %.4f sec"%(time.time() - t0))
    
    #Azim calc
    t0 = time.time()
    horX = 1.0  #horizontal vector
    horY = 0.0
    edfdata.azim = np.degrees(np.arctan2(vCPx,vCPy))
    edfdata.azim=np.where(edfdata.azim < 0, edfdata.azim+360,edfdata.azim)
    log.info(" azim calculation: %.4f sec"%(time.time() - t0))

    #T2 CALC (modulus, arctan, geomCorr, todegrees)
    t0 = time.time()
    edfdata.t2XY = vCPx**2 + vCPy**2 - vCPz**2
    edfdata.t2XY = np.sqrt(edfdata.t2XY)
    log.info(" sqrt calculation: %.4f sec"%(time.time() - t0))

    t0 = time.time()
    edfdata.t2XY = edfdata.t2XY/ (distPix-vCPz)
    edfdata.t2XY = np.arctan(edfdata.t2XY)
    log.info(" atan calculation: %.4f sec"%(time.time() - t0))

    t0 = time.time()
    edfdata.t2XY = np.degrees(edfdata.t2XY)
    log.info(" to deg calculation: %.4f sec"%(time.time() - t0))

######################################################## MAIN PROGRAM
if __name__=="__main__":

    __version__ = '161117'
    
    ts = time.time()
    st = '[' + datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S') + '] '
    welcomeMsg = '=== MSPD EDF IMAGE INTEGRATION ver.'+__version__+' by OV === '+st

    #argument parsing
    parser = argparse.ArgumentParser(description="EDF image integration")
    parser.add_argument('filenames', nargs='*', help='edf file')
    parser.add_argument("-i", "--par", default=None, help='calibration input file')
    parser.add_argument("-f2d", "--fit2d", action="store_true", help='fit2D convention for tilt/rot (default is d2Dplot convention)')
    parser.add_argument("-o", "--output", default="", help='1D output file only when processing one file (default: filename.dat')
    parser.add_argument("-s", "--suffix", default="", help='suffix for the output names (before the dat extension)')
    parser.add_argument("-ne", "--noesd", default=False, help='do not calculate esd\'s')
    parser.add_argument("-pl", "--plot", action="store_true", help='plot image')
    parser.add_argument("-d", "--debug", action="count", help='debug mode')
    args = parser.parse_args()
    
    #print args.debug
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
    
    if len(sys.argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit()

    #Entered filenames to process and output filenames to generate:
    files = []
    for i in range(len(args.filenames)):
        for filename in glob.glob(args.filenames[i]):
            files.append(filename)
    log.info("Input file(s): %s"%files)
    print "  Input file(s): %s"%(files)

    outFiles = []
    if (len(args.filenames)==1):
        if args.output == "":
            ipunt = args.filenames[0].rfind(".")
            outFiles.append(args.filenames[0][0:ipunt]+args.suffix+".dat")
        else:
            outFiles.append(args.output)
    else:
        #same names as input but with dat extension
        for fname in files:
            ipunt = fname.rfind(".")
            outFiles.append(fname[0:ipunt]+args.suffix+".dat")
    log.info("output file(s): %s"%outFiles)
    print "  Output file(s): %s"%(outFiles)
   
    #now we read the calibration parameter file
    integOpts = IntegOptions()
    if args.fit2d:
        #enable fit2d option
        integOpts.fit2d = True

    if args.par is not None: #otherwise default parameters (image headers) will be used
        log.info("reading inp file:"+args.par)
        integOpts.readINPFile(args.par)
        
        if integOpts.fit2d:
            #put start and end azim in range 0-360 starting at 12o'clock direction clockwise !!
            #Inserted as fit2d 3o'clock=zero and positive clockwise (because it is y inverted!)
            startAzim = integOpts.startAzim
            if (startAzim < 0): startAzim = startAzim + 360
            endAzim = integOpts.endAzim
            if (endAzim < 0): endAzim = endAzim + 360
            #now +90 to put at 12 (remember + clockwise)
            startAzim = startAzim +90
            endAzim = endAzim +90
        else:
            #direct values (d2dplot convention)
            startAzim = integOpts.startAzim
            endAzim = integOpts.endAzim

    isMask = False
    if integOpts.maskf is not None:
        #read the mask file
        maskbin = EDFdata()
        maskbin.readBINFile(integOpts.maskf)
        isMask = True

    for i in range(len(files)):
        fname = files[i]
        edf = EDFdata()
        edf.readFile(fname)
        
        #check for existence calib parameters (otherwise use header and default pars: tilt/rot=0, full azim,...)
        if integOpts.xcen <= 0:
            integOpts.xcen = edf.xcen
        if integOpts.ycen <= 0:
            integOpts.ycen = edf.ycen
        if integOpts.distMDmm <= 0:
            integOpts.distMDmm = edf.distMDmm
        if integOpts.wavelA <= 0:
            integOpts.wavelA = edf.wavelA
    
        #NOW WE INTEGRATE
        
        #Integration options
        t2in = calc2tDeg(integOpts, edf, integOpts.ycen,integOpts.xcen+integOpts.inRadi)
        t2fin = calc2tDeg(integOpts, edf, integOpts.ycen,integOpts.xcen+integOpts.outRadi)
        size = (int)(integOpts.radialBins)
        step = (t2fin-t2in)/size
        
        log.info("Processing %s"%(fname))
        log.info("t2in=%f t2f=%f size=%d step=%f"%(t2in,t2fin,size,step))
        
        print ""
        print "File: %s"%(fname)
        print " CenX(px)=%.3f CenY(px)=%.3f Dist(mm)=%.3f Wave(A)=%.4f PixSX(mm)=%.4f"%(integOpts.xcen,integOpts.ycen,integOpts.distMDmm,integOpts.wavelA,edf.pixSXmm)
        if integOpts.fit2d:
            print " TiltRot(º)=%.3f AngTilt(º)=%.3f startAzim(º)=%.4f endAzim(º)=%.4f [fit2d convention: 0 at 3h CW+]"%(integOpts.tiltRotationIN,integOpts.angleTilt,integOpts.startAzim,integOpts.endAzim)
        else:
            print " TiltRot(º)=%.3f AngTilt(º)=%.3f startAzim(º)=%.4f endAzim(º)=%.4f [d2dplot convention: 0 at 12h CW+]"%(integOpts.tiltRotationIN,integOpts.angleTilt,integOpts.startAzim,integOpts.endAzim)
        print " Inner/Outer Radius(Px)=%d %d (º) x/y bin= %d %d azimBins=%d radialBins=%d"%(integOpts.inRadi,integOpts.outRadi,integOpts.xbin,integOpts.xbin,integOpts.azimBins,integOpts.radialBins)
        print " T2ini(º)=%.4f step(º)=%.4f T2end(º)=%.4f"%(t2in+step/2,step,t2fin-step/2)
        if(isMask):print " Mask file= %s"%(integOpts.maskf)

        log.info(" tiltRot internally used = %f"%integOpts.tiltRotation)
        log.info(" startAzim endAzim internally used (ref. 12h CW+) = %f %f"%(startAzim,endAzim))
    
        #HEAVY CALCULATION
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
                    t2 = edf.t2XY.item(row,col)
                    azim = edf.azim.item(row,col)
                    return 'col_x=%.1f, row_y=%.1f, Intensity=%d, t2=%.3f, azim=%.2f'%(x, y, inten,t2,azim)
                else:
                    return 'col_x=%.1f, row_y=%.1f'%(x, y)
            ax.format_coord = format_coord
            plt.show()
        
        t0 = time.time()
    
        #TODO: IF NECESSARY APPLY GEOM CORRECTIONS HERE?
        #edf.intenXY = edf.intenXY*edf.gcorr
        #log.info(" apply gcorr: %.4f sec"%(time.time() - t0))
        #t0 = time.time()
    
        #to implement azimBins more than one pattern is needed (in order to do it in one image loop)
        #we have to divide the azumutal range in azimBins and relate it to the vector position
        azimRange = endAzim - startAzim
        if (endAzim<startAzim):
            azimRange = 360 + azimRange
        azimInc = (int) (azimRange / integOpts.azimBins)
        
        patts = []
        azimIni = []
        azimFin = []
        for k in range(integOpts.azimBins):
            ini = startAzim+azimInc*k
            if ini > 360: ini = ini - 360
            fin = ini + azimInc
            if fin > 360: fin = fin - 360
            azimIni.append(ini)
            azimFin.append(fin)

            patts.append(D1Pattern())
            for j in range(size+1):
                patts[k].t2.append(t2in + j*step)
                patts[k].esd.append(0)
                patts[k].Iobs.append(0)
                patts[k].npix.append(0)

        for y in range(edf.dimy):
            for x in range(edf.dimx):
                #first if there is a mask readed we check it
                if isMask:
                    if (maskbin.intenXY.item(y,x)<0):continue
    
                azim = edf.azim.item(y,x)

                if (endAzim<startAzim):
                    if (azim < startAzim) and (azim > endAzim): continue
                else:
                    if (azim < startAzim) or (azim > endAzim): continue
                
                #which position of the azim vector in case azimBins >1:
                azimPos = 0
                if (integOpts.azimBins>1):
                    for k in range(integOpts.azimBins):
                        if (azimFin[k]<azimIni[k]):
                            if (azim < azimIni[k]) and (azim > azimFin[k]): continue
                        else:
                            if (azim < azimIni[k]) or (azim > azimFin[k]): continue
                        azimPos = k
                        break
                    if azimPos < 0: 
                        log.info("error in azimPos")
                        continue
                
                t2p = edf.t2XY.item(y,x)
                if (t2p<t2in):continue
                if (t2p>t2fin):continue
    
                #position in the histogram
                p = (int)(t2p/step-t2in/step)
                inten = edf.intenXY.item(y,x)
                patts[azimPos].Iobs[p] = patts[azimPos].Iobs[p] + inten + integOpts.subadu
                patts[azimPos].npix[p] = patts[azimPos].npix[p] + 1
        
        log.info(" histogram filling: %.4f sec"%(time.time() - t0))
    
        #write the output dat(s) file
        for k in range(integOpts.azimBins):
            if (integOpts.azimBins>1):
                if k==0: iniName = outFiles[i]
                ipunt = iniName.rfind(".")
                outFiles[i] = iniName[0:ipunt]+"_bin%02d.dat"%(k)
            fout = open(outFiles[i], 'w')
            fout.write("#%s  COM = %s\n"%(fname,edf.comment))
            fout.write("#I vs. 2Theta [deg] Azim: %.2f %.2f Radial %d %d Tilts: %.2f %.2f Dist: %.2f Wave: %.4f Cen: %.2f %.2f Pix: %.2f Bin: %d\n"
                    %(integOpts.startAzim,integOpts.endAzim,integOpts.inRadi,integOpts.outRadi,integOpts.tiltRotationIN,integOpts.angleTilt,integOpts.distMDmm,integOpts.wavelA,integOpts.xcen,integOpts.ycen,edf.pixSXmm*1000,integOpts.azimBins))
            if (integOpts.azimBins>1):
                fout.write("# Azim Bin %d: from %.2f to %.2f\n"%(k,azimIni[k],azimFin[k]))
            line = ''
            for j in range(len(edf.counter_mne)):
              line = line + ' ' + edf.counter_mne[j] + ' ' + edf.counter_pos[j]
            fout.write("# %s\n"%(line))
            fout.write("# %d \n"%(size)) #number of points
            for j in range(size):
                if patts[k].npix[j] <= 0:continue
                
                #correction of t2, e.g. 0 it is really 0+(step/2) to be at the center of the bin
                t2 = patts[k].t2[j] + (step/2)
                
                #lorentz correction
                lorfact = (1./((math.cos(math.radians(t2)))**3))
                icorr = patts[k].Iobs[j] * lorfact
                if (icorr<0): icorr = 0
                inten = icorr/patts[k].npix[j]
                
                #esd (direct calculation using npix)
                if (args.noesd == False):
                    esd = math.sqrt(inten/patts[k].npix[j])
                    #write the line
                    fout.write(" %.7E  %.7E %.7E\n"%(round(t2,4), inten, esd))  #round to 4 decimals (ffops uses 3)
                else:
                    fout.write(" %.7E  %.7E\n"%(round(t2,4), inten)) #no esd
            fout.close()
            print "==> Output DAT file (xye format): %s"%(outFiles[i])
            log.info("%s dat file written"%(outFiles[i]))
        
    print ""
    print "=== END EXECUTION (total time %.4f sec) ==="%(time.time()-ts)
    print ""
    log.info(" Total execution time: %.4f sec"%(time.time() - ts))
