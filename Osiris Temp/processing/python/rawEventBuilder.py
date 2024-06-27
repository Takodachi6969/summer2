import struct
import math

class tdcEvent():
    def __init__(self, header, time, words = [], eob=None, qual=0):
        self.header = header
        self.words = words
        self.time = time
        self.EOB = eob
        self.qual = qual

class proAnubEvent():
    def __init__(self, tdcEvents):
        self.tdcEvents = tdcEvents

class eventBuilder():
    def __init__(self, activeTDCs = [0,1,2,3,4]):
        self.events = []
        self.tdcEventBuffer = [[] for tdc in activeTDCs]
        self.activeTDCs = activeTDCs
        self.missItr = [[0,0] for tdc in activeTDCs]
        self.eventCounts = [-1 for tdc in activeTDCs]
        self.headerMask = 0x400000
        self.EOBMask = 0x200000 
        self.HEADER = 0
        self.EOB = 1
        self.DATA = 2
        self.CORRUPTHEADER = 3
        # self.processedEventCount = 0    
        # self.alignmentCheckInterval = 100
        # self.order = [[0,1], [1,2], [2,4]]
        # self.chunk = []
        # self.adjustment = [0 for i in range(len(self.order))]
        # self.lastWasBad = [False for i in range(len(self.order))]
        # self.alignment_event = []

        
    def addTDCRead(self, thisTDC, thisTime, tdcReadData, p=True):
        if thisTDC not in self.activeTDCs:
            return
        readData = struct.iter_unpack("I", tdcReadData)
        thisRead = []
        splitEvtBuf = []
        lastWordType = None
        missedEvts = 0
        maybeMissed = []
        startCount = self.eventCounts[thisTDC]
        lastHead = 0
        for word in readData:
            thisWord = word[0]
            if p:
                print("TDC",thisTDC,hex(thisWord))
            wordType = self.getWordType(thisWord)
            if wordType==self.HEADER:
                #If we get a header, check the previous word to throw out extra headers
                if lastWordType==self.HEADER or lastWordType==self.CORRUPTHEADER:
                    thisRead[-1].qual = thisRead[-1].qual|0x10
                    lastWordType=self.EOB
                    continue
                else:
                    #Otherwise, create a new event.
                    #If we haven't had any events yet, start the counter with the first header in case they don't begin with zero
                    if self.eventCounts[thisTDC]<0:
                        self.eventCounts[thisTDC]=(thisWord&0x3ff)+1
                    else:
                        self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                    #Store the difference in TDC header number and the expected event number to check for missed events
                    evtDiff = int(thisWord&0x3ff)-int(self.eventCounts[thisTDC]%(0x3ff+1))+1
                    if p:
                        print("Header from TDC",thisTDC,"evtDiff:",evtDiff,"Saw:",thisWord&0x3ff,"Expected:", self.eventCounts[thisTDC]%(0x3ff+1)-1)
                    if evtDiff<0:
                        evtDiff = evtDiff+0x3ff
                    maybeMissed.append(evtDiff) 
                    thisRead.append(tdcEvent(thisWord, thisTime,[]))
            elif wordType==self.CORRUPTHEADER:
                #Can treat corrupted headers the same way, but as we know the event number is bad we don't check for missed events
                if lastWordType==self.HEADER or lastWordType==self.CORRUPTHEADER:
                    continue
                else:
                    self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                    thisRead.append(tdcEvent(thisWord, thisTime, [], qual=0x1))
            #For non-headers, need to check if we've started a word yet to see if we've got an event that's split over two reads.
            #Need to put split events into a temp buffer, as we don't yet know whether there were missed events in between.
            elif len(thisRead)==0:
                splitEvtBuf.append(thisWord)
            elif wordType==self.EOB:
                #End the current TDC event if we get an EOB, flag the event if we get two in a row.
                if lastWordType==self.EOB:
                    thisRead[-1].qual = thisRead[-1].qual+0x2
                else:
                    thisRead[-1].EOB = thisWord
            else:
                #Have a data word. Should add them into the event unless we've already got an EOB in the last event, suggesting that we missed a header. 
                if thisRead[-1].EOB is not None:
                    #If we did miss a header, make a fake one with a quality flag set to show we made one up
                    thisRead.append(tdcEvent((self.eventCounts[thisTDC]+1)%(0x3ff+1), thisTime, [], qual=0x5))
                    self.eventCounts[thisTDC]=self.eventCounts[thisTDC]+1
                thisRead[-1].words.append(thisWord)
            lastWordType = wordType
        
        #Now that we've read in the new events, want to check if we missed any in between TDC reads by taking the mode of maybeMissed
        if len(maybeMissed)>0:
            missedEvts = max(set(maybeMissed), key=maybeMissed.count)
            
        if missedEvts>0:
            #Add to the miss counter if more than one read was made 
            #and the offset occurred at least twice
            if len(maybeMissed)>1 and maybeMissed.count(missedEvts)>1:
                if p:
                    print("Inserting "+str(missedEvts)+" into TDC",thisTDC)
                for missedEvt in range(missedEvts):
                    if missedEvts<100:
                        thisRead.insert(0, tdcEvent((startCount+missedEvt)%(0x3ff+1), thisTime, [], qual=0x8))
                    self.eventCounts[thisTDC] = self.eventCounts[thisTDC]+1
                    self.missItr[thisTDC] = [0,0]
            else:
                #If there was only one TDC read, store it in an iterator and wait for another event to confirm the offset.
                if self.missItr[thisTDC][1] == missedEvts:
                    if p:
                        print("Inserting",missedEvts,"into TDC",thisTDC)
                    for evt in range(missedEvts):
                        if missedEvts<100:
                            self.tdcEventBuffer[thisTDC].insert(self.missItr[thisTDC][0],tdcEvent(0xffff, thisTime, [], qual=0x8))
                        self.missItr[thisTDC] = [0,0]
                        self.eventCounts[thisTDC] = self.eventCounts[thisTDC]+1
                else:
                    self.missItr[thisTDC] = [len(self.tdcEventBuffer[thisTDC]), missedEvts]
        
        #If there were no missed events, tack any crossover words onto the proper event from the previous TDC read
        if missedEvts==0:
            self.missItr[thisTDC]=[0,0]
            if len(splitEvtBuf)>0:
                for word in splitEvtBuf:
                    wordType = self.getWordType(word)
                    if len(self.tdcEventBuffer[thisTDC])>0:
                        if wordType==self.EOB:
                            self.tdcEventBuffer[thisTDC][-1].EOB=word
                        else:
                            self.tdcEventBuffer[thisTDC][-1].words.append(word)
                       
        #Now store the tdc events from the new read into the full buffer, and check whether we can pull out any complete proANUBIS events
        for event in thisRead:
            self.tdcEventBuffer[thisTDC].append(event)
        self.buildFullEvents()
        thisRead = []
    
    def buildFullEvents(self):
        while self.checkBufferForEvents():
            fullEvent = []
            for tdc in range(5):
                if tdc in self.activeTDCs:
                    fullEvent.append(self.tdcEventBuffer[tdc].pop(0))
                else:
                    fullEvent.append(tdcEvent(0,0,[],qual=0xff))
            self.events.append(proAnubEvent(fullEvent))
    #         self.chunk.append(proAnubEvent(fullEvent))
    #         self.processedEventCount = self.processedEventCount + 1
    #         if self.processedEventCount % self.alignmentCheckInterval == 0:
    #             self.realigner
    #         self.chunk = []
    #     return  
    
    # def realigner(self):
    #     aligned_check = []
    #     for idx, item in enumerate(self.order):
    #         i, j = item
    #         x, y, l, m = self.find_tdc_alignment_metric(i, j)
    #         aligned, update = self.doRealign(self.chunk, x, y, l, m, i, j, adjustment=self.adjustment[idx], processedEvents=self.processedEventCount)
    #         if update:
    #             self.alignment_event[idx].append(self.processedEventCount)
    #         aligned_check.append(aligned)
    #         print(f"Alignment check for pair {idx}: {aligned}, lastWasBad: {self.lastWasBad[idx]}, adjustment: {self.adjustment[idx]}")
            
    #         if self.lastWasBad[idx] and not aligned_check[idx]:
    #             self.adjustment[idx] += 1
    #         elif not aligned_check[idx]:
    #             self.lastWasBad[idx] = True
    #         else:
    #             self.lastWasBad[idx] = False
    #             self.adjustment[idx] = 0
    #     print("Alignment check results:", aligned_check)
                
    # def testAlign(self, rpc1Hits, rpc2Hits):
    #     minTimes = [300,300]
    #     minChans = [-1,-1]
    #     if len(rpc1Hits)<1 or len(rpc2Hits)<1:
    #         return -1
    #     for hit in rpc1Hits:
    #         if hit.time<minTimes[0]:
    #             minTimes[0]=hit.time
    #             minChans[0]=hit.channel
    #     for hit in rpc2Hits:
    #         if hit.time<minTimes[1]:
    #             minTimes[1]=hit.time
    #             minChans[1]=hit.channel
    #     return abs(minChans[1]-minChans[0])

    # def calcAvgAlign(self, eventList,offSet=0, i = 1, j = 2, k = 0, l = 2, tdc1 =0, tdc0 = 1, processedEvents = 0):
    #     mets = []
    #     for idx, event in enumerate(eventList):
    #         etaHits = [[] for rpc in range(6)]
    #         phiHits = [[] for rpc in range(6)]
    #         if (idx+abs(offSet))<len(eventList):
    #             if offSet<=0:
    #                 oneIdx = idx+abs(offSet)
    #                 twoIdx = idx
    #             else:
    #                 oneIdx = idx
    #                 twoIdx = idx+offSet
    #             for word in eventList[oneIdx].tdcEvents[tdc1].words:
    #                 rpc, thisHit = self.tdcChanToRPCHit(word,tdc1, processedEvents + idx)
    #                 if thisHit.eta:
    #                     etaHits[rpc].append(thisHit)
    #                 # elif thisHit.eta == False and tdc1 == 2 and thisHit.channel < 31:
    #                 #     continue
    #                 else:
    #                     phiHits[rpc].append(thisHit)
    #             for word in eventList[twoIdx].tdcEvents[tdc0].words:
    #                 rpc, thisHit = self.tdcChanToRPCHit(word,tdc0, processedEvents + idx)
    #                 if thisHit.eta:
    #                     etaHits[rpc].append(thisHit)
    #                 # elif thisHit.eta == False and tdc1 == 2 and thisHit.channel < 31:
    #                 #     continue
    #                 else:
    #                     phiHits[rpc].append(thisHit)
                        
    #             if i != -1:  
    #                 etOff = self.testAlign(etaHits[i],etaHits[j])
    #                 phOff = self.testAlign(phiHits[k],phiHits[l])
    #                 if etOff>=0 and phOff>=0:
    #                     mets.append(math.sqrt(etOff*etOff+phOff*phOff))
    #             else:
    #                 phOff = self.testAlign(phiHits[k],phiHits[l])
    #                 if phOff>=0:
    #                     mets.append(math.sqrt(phOff*phOff))
    #         if len(mets)>0:
    #             return sum(mets)/len(mets)
    #         else:
    #             return 100
    
    # def doRealign(self, events, i=1, j=2, k=0, l=2, tdc1 = 0, tdc0 = 1, adjustment = 0, max_tdc = 4, processedEvents = 0):
    #     aligned = True
    #     update = []
    #     alignMet = self.calcAvgAlign(events, offSet=0, i=i, j=j, k=k, l=l, tdc1 = tdc1, tdc0 = tdc0)
    #     if alignMet>18 and alignMet<100:
    #         aligned = False
    #         generated_list = [p for o in range(1, (4 + adjustment)) for p in (o, -o)]

    #         for testOffset in generated_list:
    #             testAlignMet = self.calcAvgAlign(events,testOffset, i=i, j=j, k=k, l=l, tdc1 = tdc1, tdc0 = tdc0, processedEvents=processedEvents)
    #             if testAlignMet<17:
    #                 if testOffset>0:
    #                     for tdc in range(tdc1 + 1):
    #                         update.append(tdc)
    #                         for fakeEvent in range(testOffset):
    #                             self.insertFakeEvent(tdc=tdc)
    #                 else:
    #                     for tdc in range(max_tdc + 1):
    #                         if tdc >= tdc0:
    #                             update.append(tdc)
    #                             for fakeEvent in range(abs(testOffset)):
    #                                 if tdc >= tdc0:
    #                                     self.insertFakeEvent(tdc=tdc)            
    #                 aligned = True
    #                 print("Found a new alignment, offsetting by",testOffset, "idx is", processedEvents, "updated TDC", update)
    #                 break
    #     return aligned, update
    
    # def find_tdc_alignment_metric(self, tdc0, tdc1):
    #     if tdc0 > tdc1:
    #         tdc0, tdc1 = tdc1, tdc0
    #     i, j, k, l = None, None, None, None
    #     if tdc0 == 0:
    #         if tdc1 == 1:
    #             i, j, k, l = 1, 2, 0, 1
    #         if tdc1 == 2:
    #             i, j, k, l = 3, 0, 3, 0
    #         if tdc1 == 3:
    #             i, j, k, l = 4, 0, 4, 0
    #         if tdc1 == 4:
    #             i, j, k, l = -1, -1, 5, 0
    #     if tdc0 == 1:
    #         if tdc1 == 2:
    #             i, j, k, l = 3, 2, 3, 1
    #         if tdc1 == 3:
    #             i, j, k, l = 5, 2, 4, 1
    #         if tdc1 == 4:
    #             i, j, k, l = -1, -1, 5, 1
    #     if tdc0 == 2:
    #         if tdc1 == 3:
    #             i, j, k, l = 4, 3, 4, 3
    #         if tdc1 == 4:
    #             i, j, k, l = -1, -1, 5, 3
    #     if tdc0 == 3:
    #         if tdc1 == 4:
    #             i, j, k, l = -1, -1, 5, 4
            
    #     return i, j, k, l
    
    
    
    def checkBufferForEvents(self):
        haveAnEvent = True
        for tdc in self.activeTDCs:
            if len(self.tdcEventBuffer[tdc])==0:
                haveAnEvent = False
            #If the last event is missing an EOB, wait for an extra read to see if it got split
            elif len(self.tdcEventBuffer[tdc])==1 and self.tdcEventBuffer[tdc][-1].EOB==None:
                haveAnEvent = False
            #Also don't write events if we might need to insert missing ones
            elif self.missItr[tdc][1] > 0:
                haveAnEvent = False
        return haveAnEvent
           
    def insertFakeEvent(self, tdc):
        self.tdcEventBuffer[tdc].insert(0,tdcEvent(0xffff, 0, [], qual=0x8))
 
    def getWordType(self, word):
        if word&self.headerMask:
            if word&self.EOBMask:
                return self.CORRUPTHEADER
            else:
                return self.HEADER
        elif word&self.EOBMask:
            return self.EOB
        else:
            return self.DATA
