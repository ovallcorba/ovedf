#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 12:00:00 2016

@author: ovallcorba

---------------
Changelog
Current:
 - 1st version

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
        self.pixSXmm = None
        self.pixSYmm = None
        self.distMDmm = None
        self.wavelA = None

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
                    size = int(line[iigual+1:len(line)-2]) # for the ;
                    log.debug("size="+str(size))
                    log.debug("typesize="+str(type(size)))                
                if line.startswith("Dim_1"):
                    self.dimx = int(line[iigual+1:len(line)-2])
                    log.debug("Dim_1="+str(self.dimx))
                if line.startswith("Dim_2"):
                    self.dimy = int(line[iigual+1:len(line)-2])
                    log.debug("Dim_2="+str(self.dimy))
                    
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
        
######################################################## MAIN PROGRAM

if __name__=="__main__":

    __version__ = '161017'

    ts = time.time()
    st = '[' + datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S') + '] '
    welcomeMsg = st+'=== OV EDF READING ver.'+__version__+' ==='

    parser = argparse.ArgumentParser(description="EDF file reading and integration")
    parser.add_argument('filename', nargs=1, help='edf file')

    parser.add_argument("-d", "--debug", action="count", help='debug mode')
    parser.add_argument("-pl", "--plot", action="store_true", help='plot image')
    parser.add_argument("-out", "--outspr", action="store_true", help='out spr file (testing)')

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
    
    # debug
    log.info("1024,847 intensity (should be 129)="+str(edf.intenXY[1024][847]))
    log.info("847,1024 intensity (should be 129)="+str(edf.intenXY[847][1024]))
    log.info("1200,1015 intensity (should not be 0)="+str(edf.intenXY[1200][1015]))

    if (args.plot):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # ax.imshow(edf.intenXY)
        
        numrows, numcols = edf.intenXY.shape
        # data = edf.intenXY
        def format_coord(x, y):
            col = int(x+0.5)
            row = int(y+0.5)
            if col>=0 and col<numcols and row>=0 and row<numrows:
                z = edf.intenXY[row,col]
                return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
            else:
                return 'x=%1.4f, y=%1.4f'%(x, y)
        
        ax.format_coord = format_coord
        plt.imshow(edf.intenXY)
        plt.colorbar()
        plt.show()
    
    if (args.outspr):
        fout = open("test.spr", 'w')
        fout.write("  2048  2048 \n")
        for y in range(edf.dimy):
            for x in range(edf.dimx):
                fout.write(" %.5e"%(edf.intenXY[x][y]))
            fout.write("\n")
            log.info("line "+str(y))
        fout.close()
        log.info("spr file written")