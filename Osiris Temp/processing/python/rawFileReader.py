#from rawEventBuilder import eventBuilder
import rawEventBuilder
import struct
import datetime
import os
import importlib

def readTimeStamp(tsData):
    ts = struct.iter_unpack("I", tsData)
    tStamp = []
    for tsWord in ts:
        tStamp.append(tsWord[0]) 
    return datetime.datetime.fromtimestamp(tStamp[0]+tStamp[1]/1000000000)

class fileReader():
    def __init__(self,filename,activeTDCs=[tdc for tdc in range(5)]):
        self.fname = filename
        importlib.reload(rawEventBuilder)
        self.evtBuilder = rawEventBuilder.eventBuilder(activeTDCs=activeTDCs)
        self.wordSize = 4
        self.data = open(self.fname,'rb')
        self.data.seek(0, os.SEEK_END)
        self.fsizeBytes = self.data.tell()
        self.data.seek(0)
        self.bytesRead = 0
        stDat = self.data.read(2*self.wordSize)
        self.st = readTimeStamp(stDat)
        self.leadWords = []

    def doneReading(self):
        return self.bytesRead==(self.fsizeBytes-2*self.wordSize)
        
    def hasEvents(self):
        return len(self.evtBuilder.events)>0

    def getEvents(self):
        evts = []
        for event in range(len(self.evtBuilder.events)):
            evts.append(self.evtBuilder.events.pop(0))
        return evts

    def readBlock(self,p=False):
       tsDat = self.data.read(2*self.wordSize)
       thisTime = readTimeStamp(tsDat)
       headDat = self.data.read(self.wordSize)
       headWord = 0
       try:
           headWord = struct.unpack("I",headDat)
       except:
           return False
       thisTDC = headWord[0]>>24
       nWords = headWord[0]&0xffffff
       if thisTDC>5:
           print("Bad TDC - number is",thisTDC,", header word was:",hex(leadWords[2]))
           return False
       if self.bytesRead+nWords>(self.fsizeBytes-2*self.wordSize):
           print("Going to over-read the file. Corrupted number of bytes? Header is", hex(leadWords[2]))
           return False
       tdcReadData = self.data.read(nWords*self.wordSize)
       self.bytesRead = self.bytesRead+self.wordSize*(nWords+3)
       #Skip the fifth TDC for now
       if thisTDC<5:
           self.evtBuilder.addTDCRead(thisTDC, thisTime, tdcReadData, p)
       return True
