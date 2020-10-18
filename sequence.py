import numpy as np
import matplotlib.pyplot as plt
#import inst as inst


#Input variables
MAXREPS=255 #Max repititions of a sequence element, used for smart wait

#Define Phase Class
class Phase:
    def __init__(self, freq, ampli, angle):
        self.freq = freq
        self.amp = ampli
        self.angle = angle

     #equality overwrite
    def __eq__(self,other):
        if type(self)!=type(other):
            return False
        if int(self.freq)==int(other.freq):
            if int(self.amp*10000)==int(other.amp*10000):
                if int(self.angle*100)==int(other.angle*100):
                    return True
        return False

    def __ne__(self,other):
        return not self==other

#Define Phase Tracking Class used for phase stability within chunk
#Pass prevChunk to make phase stable between chunks as well
class PhaseTrack:

    #phase to track and start time in points
    def __init__(self, phase, start,prevChunk=None, samplingRate=50*10**9):
        self.phase=phase
        self.start=start #in points
        self.prevAngle=0
        self.samplingRate = samplingRate
        if not prevChunk is None:
            if not prevChunk.getPhaseTrack(phase) is None:
                prevChunkreps=1 #takes into account repititions of a chunk when setting phase
                if not prevChunk.meta.reps == 'ONCE':
                    prevChunkreps=prevChunk.meta.reps
                prevTrack=prevChunk.getPhaseTrack(phase)
                self.prevAngle=prevTrack.getAngle(prevChunk.getDurationGroomed()*prevChunkreps)

    #get angle a certain number of points into the chunk
    #does not include the angle of the phase object itself
    def getAngle(self, t):
        freq=self.phase.freq
        return ((360)*(t-self.start)*freq/self.samplingRate)%360+self.prevAngle

    #Check exact phase including angle and amplitude
    def hasPhase(self, phase):
        return phase==self.phase

    #check only frequency
    def hasFreq(self,phase):
        return self.phase.freq==phase.freq

#Define Data structure that carries all metadata for a chunk
#Name, reps, id, next,jump  next, flags, trigger
class MetaChunk:
    def __init__(self,idNum,name=None,nextNum=None,flags=None,reps=None,jumpNum=None,trigger=None):
        if name is None:
            name=str(idNum)
        if nextNum is None:
            nextNum=0
        if flags is None:
            flags=['LOW','LOW','LOW','LOW']
        if reps is None:
            reps='ONCE'
        if jumpNum is None:
            jumpNum=0
        if trigger is None:
            trigger='OFF'

        self.idNum=idNum
        self.name=name
        self.nextNum=nextNum+1
        self.flags=flags
        self.reps=reps
        self.jumpNum=jumpNum+1
        self.trigger=trigger


#Define Chunk Class
#phase takes a Phase object
#dur: duration in ns
#off: wait time at the start of the pulse in ns
#reps: times chunk is repeated
#idNum: target ID number in AWG
#nextNum: next sequence element to run after pulse is completed
#jumpNum: next sequence element to run after pulse is completed if jump is activated
class Chunk:

    def __init__(self, phase, dur, off,meta, prevChunk=None, samplingRate=50*10**9):
        self.meta=meta
        self.groomedOutputs=[]
        self.trackedPhases=[]
        self.outputs=[]
        self.prevChunk=prevChunk
        self.samplingRate=samplingRate
        #loads all the tracked phases from prev Chunk into the new chunk
        if not prevChunk is None:
            if not prevChunk.trackedPhases is None:
                for trPh in prevChunk.trackedPhases:
                    newTrPh=PhaseTrack(trPh.phase,0,prevChunk)
                    self.trackedPhases.append(newTrPh)
        waveOC=[]
        #if only one set of inputs not stored in an array make it an array
        if np.size(off)==1:
            if type(off) is float or type(off) is int:
                dur=[dur]
                phase=[phase]
                off=[off]

        # once inputs are in an array

        for i,ph in enumerate(phase): #go through list of inputs

            #code to calculate the extra phase required for phase coherence for a given frequency within a chunk
            offLength=round(off[i]*self.samplingRate/10**9)
            prevAngle=0
            hasPhase=0

            #Does a phase track exist for this phase?
            if np.size(self.trackedPhases)>0:

                for tPh in self.trackedPhases:
                    if tPh.hasFreq(ph): #if yes grab the old angle
                        hasPhase=1
                        prevAngle=prevAngle+tPh.getAngle(np.size(self.outputs)+offLength)
            if hasPhase==0: #if no generate a new copy of the phase and append it
                newPT=PhaseTrack(ph,np.size(self.outputs)+offLength,self.prevChunk)
                self.trackedPhases.append(newPT)
                prevAngle=newPT.getAngle(np.size(self.outputs)+offLength)

            #generate sin with the correct phase and add to outputs
            waveT=np.linspace(0,round(dur[i]*self.samplingRate/10**9)-1,round(dur[i]*self.samplingRate/10**9))
            waveO=ph.amp*np.sin(2*np.pi*(waveT*ph.freq/self.samplingRate+(ph.angle+prevAngle)/360))
            if off[i]>0:
                offO=np.zeros(offLength)
                waveOC = np.concatenate((waveOC,offO,waveO))
                self.outputs=waveOC
            else:
                waveOC = np.concatenate((waveOC,waveO))
                self.outputs=waveOC
        self.groomChunk(0) #groom after every change



    def appendPulse(self, phase, dur, off):

        #if only one set of inputs and not stored in an array make it an array


        waveOC=[]

        if np.size(off)==1:
            if type(off) is float or type(off) is int:
                dur=[dur]
                phase=[phase]
                off=[off]
        for i,ph in enumerate(phase):

            #code to calculate the extra phase required for phase coherence for a given frequency within a chunk
            offLength=round(off[i]*self.samplingRate/10**9)
            prevAngle=0
            hasPhase=0

            if np.size(self.trackedPhases)>0:

                for tPh in self.trackedPhases:

                    if tPh.hasFreq(ph):

                        hasPhase=1
                        prevAngle=prevAngle+tPh.getAngle(np.size(self.outputs)+offLength)
            if hasPhase==0:
                newPT=PhaseTrack(ph,np.size(self.outputs)+offLength,self.prevChunk)
                self.trackedPhases.append(newPT)
                prevAngle=newPT.getAngle(np.size(self.outputs)+offLength)
            waveT=np.linspace(0,round(dur[i]*self.samplingRate/10**9)-1,round(dur[i]*self.samplingRate/10**9))
            waveO=ph.amp*np.sin(2*np.pi*(waveT*ph.freq/self.samplingRate+(ph.angle+prevAngle)/360))
            if off[i]>0:
                offO=np.zeros(round(off[i]*self.samplingRate/10**9))
                waveOC = np.concatenate((offO,waveO))
                self.outputs=np.concatenate((self.outputs,waveOC))
            else:
                waveOC = waveO
                self.outputs=np.concatenate((self.outputs,waveOC))
        self.groomChunk(0)

    #Returns duration in points
    def getDuration(self):
        return np.size(self.outputs)

    #Returns duration in points
    def getDurationGroomed(self):
        return np.size(self.groomedOutputs)

    #Grooms Chunk in accordance with Groom Waveforms VI
    #Rules:
    #1) Add 0s corresponding to Channel Delay (ns) to beginning of Chunk
    #2) Add 0s to end of Chunk corresponding to difference between max delay and Ch Delay
    #3) Add 0s to end of Chunk to ensure that Chunk is mimumum of 4800 points and is a multiple of 240 points
    def groomChunk(self, chDelay=0, maxDelay=None):
        minSize=4800
        mult=240
        if maxDelay is None:
            maxDelay = chDelay
        chDelayOutput = np.zeros(int(chDelay*self.samplingRate/10**9))
        extraDelayOutput = np.zeros(int((maxDelay-chDelay)*self.samplingRate/10**9))
        testOutput = np.concatenate((chDelayOutput,self.outputs,extraDelayOutput))
        if np.size(testOutput)%mult > 0:
            multSup = np.zeros(240-np.size(testOutput)%mult)
            testOutput = np.concatenate((testOutput,multSup))

        if np.size(testOutput) < minSize:
            sizeSup=np.zeros(minSize-np.size(testOutput))
            testOutput = np.concatenate((testOutput,sizeSup))

        self.groomedOutputs=testOutput

    #looks for and returns phaseTrack assosciated with the frequency of the given phase if it exists
    def getPhaseTrack(self, phase):
        for tPh in self.trackedPhases:
            if tPh.hasFreq(phase):
                return tPh
        return None


class Sequence:
    def __init__(self, maxReps=255):
        self.maxReps = maxReps
        self.chunkList=[]

    def overwriteChunk(self, phase, dur, off,meta,prevChunk):
        idNum=meta.idNum
        if idNum < self.getNextID():
            self.chunkList[idNum] = Chunk(phase, dur, off,meta,prevChunk)
            return self.getC(idNum)
        else:
            print("Warning: attempted to overwrite sequence at too high an id number")
            print("Chunk appended at ID number " + str(self.getNextID()))
            meta.idNum=self.getNextID()
            self.appendChunk(phase, dur, off,meta,prevChunk)
            return self.getC(idNum)


    def appendChunk(self, phase, dur, off,meta=None,prevChunk=None):
        if meta is None:
            meta=MetaChunk(self.getNextID())
        self.chunkList.append(Chunk(phase, dur, off,meta,prevChunk))
        return self.chunkList[meta.idNum]

    # generates conditional xy8-N for N given by reps and adds it to the queue as a single chunk
    # FinalPi is -1,0,1. If 1 uses phase1, if 0 deletes final pi, if -1 inverts the direction of the final pi/2
    def appendXY8(self, phase1, phase2, pi1, pi2, tau, num, meta=None, firstPi=None, finalPi=None, prevChunk=None):
        if finalPi is None:
            finalPi=0
        if firstPi is None:
            firstPi=0
        if meta is None:
            meta=MetaChunk(self.getNextID())

        if firstPi==1 :
            xyChunk=self.appendChunk(phase1,pi1/2,0,meta, prevChunk)
        else:
            xyChunk=self.appendChunk(phase1,0,0,meta,prevChunk)

        N=np.linspace(1,num,num)
        for n in N:
            xyChunk.appendPulse(phase1, pi1, tau/2-pi1/2)
            xyChunk.appendPulse(phase2, pi2, tau-pi1/2-pi2/2)
            xyChunk.appendPulse(phase1, pi1, tau-pi1/2-pi2/2)
            xyChunk.appendPulse(phase2, pi2, tau-pi1/2-pi2/2)
            xyChunk.appendPulse(phase2, pi2, tau-pi2/2-pi2/2)
            xyChunk.appendPulse(phase1, pi1, tau-pi1/2-pi2/2)
            xyChunk.appendPulse(phase2, pi2, tau-pi1/2-pi2/2)
            xyChunk.appendPulse(phase1, pi1, tau-pi1/2-pi2/2)
            #last wait segment
            xyChunk.appendPulse(phase1, 0, tau/2-pi1/2)
        if finalPi == 1:
            xyChunk.appendPulse(phase1,pi1/2,0)
        if finalPi == -1:
            xyChunk.appendPulse(Phase(phase1.freq,phase1.amp,phase1.angle+180),pi1/2,0)
        return xyChunk

    # generates unconditional xy8-N for N given by reps and adds it to the queue as a single chunk
    # FinalPi is -1,0,1. If 1 uses phase1, if 0 deletes final pi, if -1 inverts the direction of the final pi/2
    def appendUXY8(self, phase1X, phase1Y, phase2X, phase2Y, pi1, pi2, tau, num, meta=None, finalPi=None, prevChunk=None):
        if finalPi is None:
            finalPi=1
        if meta is None:
            meta=MetaChunk(self.getNextID())
        xyChunk=self.appendChunk(phase1X,pi1/2,0,meta, prevChunk)
        xyChunk.appendPulse(phase2X,pi2/2,0)

        combPi=pi1+pi2
        N=np.linspace(1,num,num)
        for n in N:
            xyChunk.appendPulse(phase1X, pi1, tau/2-combPi/2)
            xyChunk.appendPulse(phase2X, pi2, 0)

            xyChunk.appendPulse(phase1Y, pi1, tau-combPi)
            xyChunk.appendPulse(phase2Y, pi2, 0)

            xyChunk.appendPulse(phase1X, pi1, tau-combPi)
            xyChunk.appendPulse(phase2X, pi2, 0)

            xyChunk.appendPulse(phase1Y, pi1, tau-combPi)
            xyChunk.appendPulse(phase2Y, pi2, 0)

            xyChunk.appendPulse(phase1Y, pi1, tau-combPi)
            xyChunk.appendPulse(phase2Y, pi2, 0)

            xyChunk.appendPulse(phase1X, pi1, tau-combPi)
            xyChunk.appendPulse(phase2X, pi2, 0)

            xyChunk.appendPulse(phase1Y, pi1, tau-combPi)
            xyChunk.appendPulse(phase2Y, pi2, 0)

            xyChunk.appendPulse(phase1X, pi1, tau-combPi)
            xyChunk.appendPulse(phase2X, pi2, 0)
            #last wait segment
            xyChunk.appendPulse(phase1X, 0, tau/2-combPi/2)
        if finalPi == 1:
            xyChunk.appendPulse(phase1X,pi1/2,0)
            xyChunk.appendPulse(phase2X,pi2/2,0)
        if finalPi == -1:
            xyChunk.appendPulse(Phase(phase1X.freq,phase1X.amp,phase1X.angle+180),pi1/2,0)
            xyChunk.appendPulse(Phase(phase2X.freq,phase2X.amp,phase2X.angle+180),pi2/2,0)
        return xyChunk


    # Unconditional Pi Echo with optional RF pulses between every second pi pulse
    # phase 1 and 2 are electron phases for unconditional pi
    # dur RF is the duration of the RF pulse. Will be shortened to fit within the tau window if too long
    def appendUnEchoNPi(self, phase1, phase2, phaseRF, pi1, pi2, durRF, tau, num, meta=None, prevChunk=None):
        if meta is None:
            meta=MetaChunk(self.getNextID())

        if durRF > tau - pi1/2 - pi2/2:
            durRF = tau - pi1/2 - pi2/2

        tauRed=(tau-durRF-pi1/2-pi2/2)/2
        #generate first echo component
        echoChunk=self.appendChunk(phase1,pi1,0,meta, prevChunk)
        echoChunk.appendPulse(phase2,pi2,0)
        echoChunk.appendPulse(phaseRF,durRF,tauRed)
        echoChunk.appendPulse(phase1,pi1,tauRed)
        echoChunk.appendPulse(phase2,pi2,0)
        echoChunk.appendPulse(phase1,0,tau-pi1/2-pi2/2)

        if num > 1:
            N=np.linspace(1,num-1,num-1)
            for n in N:
                echoChunk.appendPulse(phase1,pi1,0)
                echoChunk.appendPulse(phase2,pi2,0)
                echoChunk.appendPulse(phaseRF,durRF,tauRed)
                echoChunk.appendPulse(phase1,pi1,tauRed)
                echoChunk.appendPulse(phase2,pi2,0)
                echoChunk.appendPulse(phase1,0,tau-pi1/2-pi2/2)

        return echoChunk



    # 3 CPis, Electron, Nuclear, Electron with off in between them
    def appendSwap(self, phase1, phase2, pi1, pi2, off, meta=None, prevChunk=None):
        if meta is None:
            meta=MetaChunk(self.getNextID())
        swapChunk=self.appendChunk(phase1,pi1,off,meta,prevChunk)
        swapChunk.appendPulse(phase2, pi2, off)
        swapChunk.appendPulse(phase1, pi1, off)
        return swapChunk

    def appendWait(self,dur,meta=None, prevChunk=None):
        if meta is None:
            meta=MetaChunk(self.getNextID())
        waitChunk=self.appendChunk(Phase(1,1,1),0,dur, meta, prevChunk)

        return waitChunk

    #Find the optimal duration/repition combination for a given wait duration, ns error tolerance, and self.maxReps
    #Then append this "smart" wait to the sequence
    def appendSmartWait(self,dur,tol,meta=None,prevChunk=None):
        if meta is None:
            meta=MetaChunk(self.getNextID())

        repList=np.linspace(1,self.maxReps,self.maxReps)
        error=abs(np.round(dur/repList)*repList-dur)
        repList=np.where(error<tol,repList,0)
        reps=int(np.max(repList))
        meta.reps=reps

        smChunk=self.appendChunk(Phase(1,1,1),0, int(dur/self.maxReps), meta, prevChunk)

        return smChunk

    #Print a table of the current chunks
    def tabulate(self):
        elements=np.linspace(0,np.size(self.chunkList)-1,np.size(self.chunkList))
        print('Name     Dur    DurGroom   Reps     Next     Jump')
        print('__________________________________________________')
        for i in elements:
            i=int(i)
            print(self.getC(i).meta.name + '       ' + str(self.getC(i).getDuration()) + '     ' + str(self.getC(i).getDurationGroomed()) + '     '
                  + str(self.getC(i).meta.reps)+ '       ' + str(self.getC(i).meta.nextNum)+ '        ' + str(self.chunkList[i].meta.jumpNum))

    #Get a chunk from the list
    def getC(self,n):
        for ch in self.chunkList:
            if ch.meta.idNum==n:
                return ch

    #Plot a chunk from the list
    def plotChunk(self,n):
        plt.plot(np.linspace(0,np.size(self.getC(n).outputs)-1,np.size(self.getC(n).outputs)),self.getC(n).outputs)
        return self.getC(n)

    #Plot a groomed chunk from the list
    def plotGroomedChunk(self,n):

        gp=plt.plot(np.linspace(0,np.size(self.getC(n).groomedOutputs)-1,np.size(self.getC(n).groomedOutputs)),self.getC(n).groomedOutputs)
        return gp

    #Applies grooming procedure to all chunks- call as part of the loading sequence
    def groomAll(self, chDelay, maxDelay=None):
        if maxDelay is None:
            maxDelay = chDelay
        for ch in self.chunkList:
            ch.groomChunk(chDelay, maxDelay)

    def getNextID(self):
        return np.size(self.chunkList)

    #load sequence into instrument- should be a inst object as defined in inst.py
    #If not run in notebook form this may need to be modified to indicate where the AWG interface commands are
    def loadSequence(self,instrument,seqName=None):
        if seqName is None:
            seqName='ExpSequence'
        self.groomAll(0)
        print(np.size(self.chunkList))
        instrument.configureNewSequence(seqName,np.size(self.chunkList))
        for ch in self.chunkList:
            name=str(ch.meta.idNum)+'_'+str(ch.meta.name)
            instrument.loadsinglewaveform(name,ch.groomedOutputs)
            meta=ch.meta
            instrument.addSeqElements(seqName,name,'OFF',
                              meta.flags[0],meta.flags[1],meta.flags[2],meta.flags[3]
                           ,meta.reps,meta.trigger,meta.nextNum,meta.jumpNum)

