import pyvisa
import numpy as np
import time


class Driver():

    def __init__(self, gpib_address):
        """Initializes the AWG
        
        :gpib_address: (str) Ressource address of AWG, e.g.
        TCPIP0::140.247.189.152::inst0::INSTR

        Can be found via
                rm = pyvisa.ResourceManager()
                rm.list_resources()
        """
        rm = pyvisa.ResourceManager()
        self.instrument = rm.open_resource('TCPIP0::140.247.189.152::inst0::INSTR')

    def initialize(self):
        
        self.instrument.write('*IDN?')
        return self.instrument.read()
    
    def prepareforsequenceloading(self):
        
        self.instrument.write('WLIS:WAV:DEL ALL')
        frequency = 100000000
        self.instrument.write('SOUR:FREQ {frequency}'.format(frequency=frequency))
        print('SOUR:FREQ {frequency}'.format(frequency=frequency))
        self.instrument.write('*OPC?')
        print(self.instrument.read())
        ratein=50*10**9
        self.instrument.write('CLOCK:SRATE {rate}'.format(rate=ratein))
        self.instrument.write('CLOCK:SRATE?')
        out=self.instrument.read()
        print(out)
        self.instrument.write('SOUR:FREQ?')
        out = self.instrument.read()
        return out

    def importwaveformfromrawdata(self,waveform_name,waveform):
    
        size = 1024
        
        self.instrument.write('WLIS:WAV:NEW "{name}", {size}'.format(name=waveform_name,size=size))
        self.instrument.write('SYST:ERR?')
        print(self.instrument.read())
    
        #self.instrument.write('WLIS:WAV:DATA "{name}",0,{sizef},#44096{waveform}'.format(name=waveform_name, 
        #                                                                            sizef = size, waveform=waveform))
    
        #self.instrument.write_binary_values('WLIS:WAV:DATA "Test1",0,{sizef},#44096'.format(sizef= size), waveform)

        self.instrument.write_binary_values('WLIS:WAV:DATA "Test2",0,', waveform)
    
        #Test raw labveiw input:
        #self.instrument.write(visaWrite)
    
        #print('WLIS:WAV:DAT "{name}",0,{size},#72001600{waveform}'.format(name=waveform_name, size=size, waveform = waveform))
        #WLIST:WAVEFORM:DATA "TestWfm",0,1024,#44096
    
        self.instrument.write('*OPC?')
        print(self.instrument.read())
        self.instrument.write('SYST:ERR?')
        print(self.instrument.read())
        #:A "%s",0,%d,#%d%d%s;
        #:WLIS:WAV:MARK:DATA "%s",0,%d,#%d%d%s;
    
        return
    
    def importwaveformfromrawdata(self,waveform_name,waveform):
        
        size = len(waveform)
    
        #Define a new waveform
        self.instrument.write('WLIS:WAV:NEW "{name}", {size}'.format(name=waveform_name,size=size))
    
        #self.instrument.write('SYST:ERR?')
        #print(self.instrument.read())

        #Load data to waveform
        self.instrument.write_binary_values('WLIS:WAV:DATA "{name}",0,'.format(name = waveform_name), waveform.tolist())  

        self.instrument.write('SYST:ERR?')

        return self.instrument.read()


    def loadwaveforms(self,names,waveforms):
        
        #Delete previous waveforms
        self.instrument.write('WLIS:WAV:DEL ALL')
        self.instrument.write('SLIS:SEQ:DEL ALL')
    
        #Load list of waveforms
        for i in np.arange(np.size(names)):
            importwaveformfromrawdata(self.instrument,names[i],waveforms[i])
        
        self.instrument.write('SYST:ERR?')
        return self.instrument.read()


    def loadsinglewaveform(self,names,waveforms):
        
        self.instrument.write("WLIS:LIST?")
        out = self.instrument.read()
        out=','+out[1:-2]+','
        if out.find(','+names+',')>-1:
            self.instrument.write('WLIS: DEL '+names)
        self.instrument.write('WLIS:WAV:DEL "'+str(names)+'"')
        #Load list of waveforms

        self.importwaveformfromrawdata(names,waveforms)
        
        self.instrument.write('SYST:ERR?')
        return self.instrument.read()

    def configuresequenceinchunks(self,sequence_name,element_names,wait_triggers,
                              flagsA,flagsB,flagsC,flagsD,repeats,jump_triggers,gotos,jump_steps):
        
        step_num = np.size(element_names)
    
        #Delete old sequences
        self.instrument.write('SLIS:SEQ:DEL ALL')
    
        #Add a new sequence
        self.instrument.write('SLIS:SEQ:NEW  "{name}", {size}'.format(name=sequence_name,size=step_num))
    
        #Whether jump happens immediately or at end of step, default = end
        self.instrument.write('SLIS:SEQ:EVEN:JTIM "{name}", END'.format(name=sequence_name))
    
        #Whether flags toggle on repeat, default = off
        self.instrument.write('SLIS:SEQ:RFL "{name}", OFF'.format(name=sequence_name))
    
    
        for step in np.arange(step_num):
        
            #Add a step for each wavform
            waveform = element_names[step]
            self.instrument.write('SLIS:SEQ:STEP{step}:TASS1:WAV "{sequence}", "{waveform}"'.format(step=step+1,sequence=sequence_name,waveform = waveform))
    
            #Whether to wait for A or B trigger, default = off
            wtrigger = wait_triggers[step] #'ATR','BTR','OFF'
            self.instrument.write('SLIS:SEQ:STEP{step}:WINP "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=wtrigger))
    
            #Set flags
            flagA = flagsA[step] #HIGH,LOW,PULSE,TOGG,NCH
            flagB = flagsB[step]
            flagC = flagsC[step]
            flagD = flagsD[step]
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:AFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagA))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:BFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagB))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:CFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagC))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:DFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagD))
    
            #Set repeat count
            repeat = repeats[step] #'INF','ONCE'
            self.instrument.write('SLIS:SEQ:STEP{step}:RCO "{sequence}",{repeat}'.format(step=step+1,sequence=sequence_name,repeat=repeat))
    
            #Whether to jump on A or B trigger, default = off
            jtrigger = jump_triggers[step] #'ATR','BTR','OFF'
            self.instrument.write('SLIS:SEQ:STEP{step}:EJIN "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=jtrigger))
    
            #Sets Goto
            goto = gotos[step] #'NEXT','FIRS','LAST','END',#
            self.instrument.write('SLIS:SEQ:STEP{step}:GOTO "{sequence}",{goto}'.format(step=step+1,sequence=sequence_name,goto=goto))
    
            #sets step to jump to
            jumpto = jump_steps[step] #'NEXT','FIRS','LAST','END',#
            self.instrument.write('SLIS:SEQ:STEP{step}:EJUM "{sequence}",{jumpto}'.format(step=step+1,sequence=sequence_name,jumpto=jumpto))

        return

    def load_seq_csv(self,file):
    
        seq_table = np.loadtxt(file,dtype='str',skiprows=1,delimiter=',',encoding='utf8')

        step_names = seq_table[:,1]

        wait_triggers  = ['OFF']*np.size(step_names)

        jump_triggers = seq_table[:,7]

        jump_steps = seq_table[:,8]

        gotos = seq_table[:,9]

        repeats = seq_table[:,6]

        event_jumps = np.zeros_like(step_names)

        flagsA = seq_table[:,2]
        flagsB = seq_table[:,3]
        flagsC = seq_table[:,4]
        flagsD = seq_table[:,5]

        return [step_names,wait_triggers,flagsA,flagsB,flagsC,flagsD,repeats,jump_triggers,gotos,jump_steps]
    
    def generateDJtable(self, seqname, pattern,step):
        
        #Check number of steps
        self.instrument.write('SLIST:SEQUENCE:LENGTH? "{seq}"'.format(seq=seqname))
        a = self.instrument.read()
        print(a)
        steps = int(a) 
    
    
        #Deleted two commas after JSTROBE and SEDG (One Each) Formerly
        #self.instrument.write('AWGCONTROL:PJUMP:JSTROBE, ON')
        #self.instrument.write('AWGC:PJUM:SEDG, RIS')
        #Turn strobe jump on, regardless of address pattern change
        self.instrument.write('AWGCONTROL:PJUMP:JSTROBE ON')
    
        #Trigger on rising strobe
        self.instrument.write('AWGC:PJUM:SEDG RIS')

        #Enable pattern jump for current sequence
        self.instrument.write('SLIS:SEQ:EVEN:PJUM:ENAB "{seq}", ON'.format(seq=seqname))

        #If only one pattern jump, convert to array
        if np.size(step) == 1:
            step=[step]
            pattern=[pattern]
    
        #Iterate through requested pattern jumps and write only if they point to an existing step
        for i, stepE in enumerate(step):
            patternE=pattern[i]
            if stepE < steps+1:
                self.instrument.write('SLIS:SEQ:EVEN:PJUM:DEF "{seq}", {pat}, {stp}'.format(seq=seqname,pat=patternE,stp=int(stepE)))

        return

    def assignsequence(self,seqname):
        
        self.instrument.write('CASS:SEQ "{seqname}", 1'.format(seqname=seqname))
        time.sleep(3)
        self.instrument.write('*OPC?')
        print(self.instrument.read())
    
        return

    def configureNewSequence(self,sequence_name,step_num):
        
        #Delete old sequences
        self.instrument.write('SLIS:SEQ:DEL ALL')
    
        #Add a new sequence
        self.instrument.write('SLIS:SEQ:NEW  "{name}", {size}'.format(name=sequence_name,size=step_num))
        
        #Whether jump happens immediately or at end of step, default = end
        self.instrument.write('SLIS:SEQ:EVEN:JTIM "{name}", END'.format(name=sequence_name))
        
        #Whether flags toggle on repeat, default = off
        self.instrument.write('SLIS:SEQ:RFL "{name}", OFF'.format(name=sequence_name))
        
        return

    def addSeqElements(self,sequence_name,element_names,wait_triggers,
                                flagsA,flagsB,flagsC,flagsD,repeats,jump_triggers,gotos,jump_steps,stepTar=None):
        
        step_num = np.size(jump_steps)
        if step_num==1:
            if stepTar==None:
                undsc=element_names.find('_')
                stepTar=int(element_names[0:int(undsc)])
            step=stepTar
            #Add a step for each wavform
            waveform = element_names
            self.instrument.write('SLIS:SEQ:STEP{step}:TASS1:WAV "{sequence}", "{waveform}"'.format(step=step+1,sequence=sequence_name,waveform = waveform))
        
            #Whether to wait for A or B trigger, default = off
            wtrigger = wait_triggers #'ATR','BTR','OFF'
            self.instrument.write('SLIS:SEQ:STEP{step}:WINP "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=wtrigger))
        
            #Set flags
            flagA = flagsA #HIGH,LOW,PULSE,TOGG,NCH
            flagB = flagsB
            flagC = flagsC
            flagD = flagsD
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:AFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagA))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:BFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagB))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:CFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagC))
            self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:DFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagD))
        
            #Set repeat count
            repeat = repeats #'INF','ONCE'
            self.instrument.write('SLIS:SEQ:STEP{step}:RCO "{sequence}",{repeat}'.format(step=step+1,sequence=sequence_name,repeat=repeat))
        
            #Whether to jump on A or B trigger, default = off
            jtrigger = jump_triggers #'ATR','BTR','OFF'
            self.instrument.write('SLIS:SEQ:STEP{step}:EJIN "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=jtrigger))
        
            #Sets Goto
            goto = gotos #'NEXT','FIRS','LAST','END',#
            self.instrument.write('SLIS:SEQ:STEP{step}:GOTO "{sequence}",{goto}'.format(step=step+1,sequence=sequence_name,goto=goto))
        
            #sets step to jump to
            jumpto = jump_steps #'NEXT','FIRS','LAST','END',#
            self.instrument.write('SLIS:SEQ:STEP{step}:EJUM "{sequence}",{jumpto}'.format(step=step+1,sequence=sequence_name,jumpto=jumpto))
        else:
            for step in np.arange(step_num):
            
                #Add a step for each wavform
                waveform = element_names[step]
                self.instrument.write('SLIS:SEQ:STEP{step}:TASS1:WAV "{sequence}", "{waveform}"'.format(step=step+1,sequence=sequence_name,waveform = waveform))
        
                #Whether to wait for A or B trigger, default = off
                wtrigger = wait_triggers[step] #'ATR','BTR','OFF'
                self.instrument.write('SLIS:SEQ:STEP{step}:WINP "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=wtrigger))
        
                #Set flags
                flagA = flagsA[step] #HIGH,LOW,PULSE,TOGG,NCH
                flagB = flagsB[step]
                flagC = flagsC[step]
                flagD = flagsD[step]
                self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:AFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagA))
                self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:BFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagB))
                self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:CFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagC))
                self.instrument.write('SLIS:SEQ:STEP{step}:TFL1:DFL "{sequence}",{flag}'.format(step=step+1,sequence=sequence_name,flag=flagD))
        
                #Set repeat count
                repeat = repeats[step] #'INF','ONCE'
                self.instrument.write('SLIS:SEQ:STEP{step}:RCO "{sequence}",{repeat}'.format(step=step+1,sequence=sequence_name,repeat=repeat))
        
                #Whether to jump on A or B trigger, default = off
                jtrigger = jump_triggers[step] #'ATR','BTR','OFF'
                self.instrument.write('SLIS:SEQ:STEP{step}:EJIN "{sequence}",{trigger}'.format(step=step+1,sequence=sequence_name,trigger=jtrigger))
        
                #Sets Goto
                goto = gotos[step] #'NEXT','FIRS','LAST','END',#
                self.instrument.write('SLIS:SEQ:STEP{step}:GOTO "{sequence}",{goto}'.format(step=step+1,sequence=sequence_name,goto=goto))
        
                #sets step to jump to
                jumpto = jump_steps[step] #'NEXT','FIRS','LAST','END',#
                self.instrument.write('SLIS:SEQ:STEP{step}:EJUM "{sequence}",{jumpto}'.format(step=step+1,sequence=sequence_name,jumpto=jumpto))

        return    