'''
Test
'''


import numpy as np
import re
import os

class picomax:
    '''
    A wrapper class for using picomax in python
    '''
    # susceptibility matrix size
    nchi:int = 1

    # planewave cutoff
    encut:float = 200
    # symmetrized (spherical) planewave grid
    gsym:int = 1

    # number of electronic bands
    nband:int = 8
    # number of valence electronic bands
    nvalence:int = 4

    # number of q-points
    nqpt:int = 0
    iqpt = np.array([0],dtype=np.int32)
    # qpt = np.array([[0,0,0]],dtype=np.float64)
    qnum = []
    qpath:str = ""
    # gammazero = [0,0,0]

    # number of frequency points
    omega = np.array([],dtype=np.float64)
    nomega:int = 1
    domega:float = 0.1
    omega2:float = 0
    epsilon:float = 0.1
    kk_transform:int = 0
    delta:int = 0

    # BZ integration (k-point)
    kptfile:str
    kptopt:int = 0
    kptorder:int = 0


    # refpoint = [0,0,0]

    debug:int = 0
    wdir:str = '/tmp'
    outputfile:str = 'tmpepsilon'

    picomax_exe:str = '/home/mjh92/proj/picomax/src/picomax'

    crystal:str = 'zincblende'
    a:float
    f:float = -1
    u:float = -1
    vg = np.array([],dtype=np.float32)

    Xij = []

    def __init__(self):
        self._refpoint = 'G'
        self._gammazero = [0,0,0]
        self._qpt = np.array([[0,0,0]],dtype=np.float64)

    @property
    def refpoint(self):
        return self._refpoint
    @refpoint.setter
    def refpoint(self,v):
        if len(v)==1 and v.isalpha() and v.isupper():
            self._refpoint = v
        elif len(v)==3:
            self._refpoint = v
        else:
            raise ValueError("Check the argument: refpoint")
        
    @property
    def gammazero(self):
        return self._gammazero
    @gammazero.setter
    def gammazero(self,v):
        if len(v)!=3:
            raise ValueError("Check the argument: gammazero")
        self._gammazero = v

    @property
    def qpt(self):
        return self._qpt
    @qpt.setter
    def qpt(self,v):
        if isinstance(v,list):
            v = np.array(v,dtype=np.float64)
        
        if isinstance(v,np.ndarray) and len(v.shape)==1:
            v = np.array([v],dtype=np.float64)
        
        if v.shape[1]!=3:
            raise ValueError("Check the argument: qpt")
        
        self._qpt = v

    def set(self, **kwargs):
        for k,v in kwargs.items():
            setattr(self, k, v)

    def run(self,switch:str):
        '''
        Run PicoMax
        '''
        if switch == 'eband':
            cmd = self.gencmd_eband()
        elif switch == 'Xij':
            cmd = self.gencmd_Xij()
        elif switch == 'Xxyz':
            cmd = self.gencmd_Xij()
        elif switch == 'XLL':
            cmd = self.gencmd_Xij()
        else:
            print("Check out the run option")

        # run picomax through system terminal
        print(cmd)
        os.system(cmd)

        # load results
        filename = self.wdir + '/' + self.outputfile
        if switch == 'eband':
            self.load(filename+'.qpt','qpt')
            self.load(filename+'.dat','eband')
        elif switch == 'Xij':
            self.load(filename+'.qpt','qpt')
            self.load(filename+'.gpt','gpt')
            self.load(filename+'.fpt','fpt')
            self.load(filename+'.dat','Xij')
        elif switch == 'Xxyz':
            self.load(filename+'.qpt','qpt')
            self.load(filename+'.qpt','gpt')
            self.load(filename+'.qpt','fpt')
            self.load(filename+'.dat','Xij')
        elif switch == 'XLL':
            self.load(filename+'.qpt','qpt')
            self.load(filename+'.qpt','gpt')
            self.load(filename+'.qpt','fpt')
            self.load(filename+'.dat','Xij')
        print("Finished importing the calculation results!")
    
    def load(self,filename:str,switch:str):
        '''
        Load PicoMax
        '''
        if switch == 'eband':
            self.load_eband(filename)
        elif switch == 'Xij':
            self.load_Xij(filename)
        elif switch == 'Xxyz':
            self.load_Xij(filename)
        elif switch == 'XLL':
            self.load_Xij(filename)
        elif switch == 'gpt':
            self.load_gpt(filename)
        elif switch == 'kpt':
            self.load_gpt(filename)
        elif switch == 'fpt':
            self.load_fpt(filename)
        elif switch == 'qpt':
            self.load_qpt(filename)
        else:
            print("Check out the load option")

    def load_eband(self,filename:str):
        file = open(filename,'rb')
        data = np.fromfile(file, dtype=np.float64)
        nqpt = self.nqpt
        nband = int(len(data)/nqpt)
        self.eband = data.reshape((nqpt,nband))

    def load_Xij(self,filename:str):
        file = open(filename,'rb')
        data = np.fromfile(file, dtype=np.float64)
        nomega = self.nomega
        nqpt = self.nqpt
        ne = int(self.nchi * (self.nchi+1) / 2)
        N = 3*3*ne*nqpt*nomega
        self.Xij = np.transpose(data[0:N].reshape((3,3,nqpt,ne,nomega)),(4,2,3,0,1))\
                +1j*np.transpose(data[N:2*N].reshape((3,3,nqpt,ne,nomega)),(4,2,3,0,1))

    def load_gpt(self,filename:str):
        file = open(filename,'rb')
        data = np.fromfile(file, dtype=np.float64)
        self.nchi = int(len(data)/3)
        self.gpt = data.reshape((self.nchi,3))

    def load_qpt(self,filename:str):
        file = open(filename,'rb')
        data = np.fromfile(file, dtype=np.float64)
        self.nqpt = int(len(data)/3)
        self.qpt = data.reshape((self.nqpt,3))
        self.iqpt = np.arange(0,self.nqpt,1)

    def load_fpt(self,filename:str):
        file = open(filename,'rb')
        data = np.fromfile(file, dtype=np.float64)
        self.nomega = int(len(data))
        self.omega = data

    def gencmd_Xij(self)->str:
        cmd = self.gencmd_base() \
            + '-switch Xij ' \
            + self.gencmd_sys() \
            + self.gencmd_qpt() \
            + self.gencmd_kpt() \
            + self.gencmd_fpt()
        return cmd

    def gencmd_eband(self)->str:
        cmd = self.gencmd_base() \
            + '-switch eband ' \
            + self.gencmd_sys() \
            + self.gencmd_qpt()
        return cmd

    def gencmd_base(self)->str:
        return "%s -wdir %s -outputfile %s " % (self.picomax_exe,self.wdir,self.outputfile)

    def gencmd_kpt(self)->str:
        return "-kptopt %d -kptorder %d " % (self.kptopt,self.kptorder)

    def gencmd_sys(self)->str:
        cmd = "-crystal %s -a %f -vg %s "\
                    % (self.crystal,self.a,arr2csl(self.vg))
        if self.f != -1:
            cmd += "-f %f " % (self.f)

        if self.u != -1:
            cmd += "-u %f " % (self.u)

        return cmd + "-encut %d -gsym %d -nchi %d -nband %d -nvalence %d "\
                        % (self.encut,self.gsym,self.nchi,self.nband,self.nvalence)

    def gencmd_fpt(self)->str:
        return "-nomega %d -domega %f -omega2 %f -kk_transform %d -epsilon %f -delta %d "\
                    % (self.nomega,self.domega,self.omega2,self.kk_transform,self.epsilon,self.delta)

    def gencmd_qpt(self)->str:
        if ((len(self.qnum) == 0) or (len(self.qpath) == 0)):
            cmd = "-qvec %s " % (arr2csl(self.qpt.flatten()))
            self.nqpt = self.qpt.shape[1]
        else:
            cmd = "-qnum %s -qpath %s " % (arr2csl(self.qnum),self.qpath)
            self.nqpt = np.sum(self.qnum)+1
        
        return cmd + "-refpoint %s -gammazero %s "\
                        % (arr2csl(self.refpoint),arr2csl(self.gammazero))

def arr2csl(v:np.array)->str:
    '''
    Converts a numpy array to a comma-separated-list string
    '''
    return re.sub(r'\s+', ',', ' '.join(map(str,v)))


