import argparse
import os
import math

#JOULE_PER_CAL = 4.184

class Atom:
    """An atom class"""
    def __init__(self, charge, lj):
        self.charge = charge
        self.epsilon = lj[0]
        self.sigma = lj[1]

#    def __eq__(self, other):
#        return (self.charge == other.charge and 
#                self.epsilon == other.epsilon and
#                self.sigma == other.sigma)

class Molecule:
    """A molecule class"""
    def __init__(self, atoms):
        self.sites = atoms

#    def __eq__(self, other):
#        return self.sites == other.sites

class System:
    """A system class"""
    def __init__(self, volume, solute, solvent, slvNum):
        self.volume = volume
        self.solute = solute
        self.solvent = solvent
        self.slvNum = dict(zip(solvent, slvNum))
        self.slvNumDensity = dict(zip(solvent, [self.slvNum[slv] / volume for slv in solvent]))

    def getElrc(self, rs, rc):
        elrc=0.0
        for slv in self.solvent:
            C0 = 16.0 * math.pi * self.slvNumDensity[slv]
            if rs != rc:
                C3_0 = (rcutoff**2 - rswitch**2)**3
                C3_2 = 3.0 * rswitch**2 + 3.0 * rcutoff**2
                C3_3 = 6.0 * rcutoff**2 * rswitch**2
                C3_4 = rcutoff**6 - 3.0 * rcutoff**4 * rswitch**2

            for slvSite in slv.sites:
                for sltSite in self.solute.sites:
                    eps = math.sqrt(slvSite.epsilon * sltSite.epsilon)
                    sig_2 = slvSite.sigma * sltSite.sigma    
                    sig_6 = sig_2**3
                    sig_12 = sig_6 * sig_6 

                    term1 = sig_12 / (9.0 * rswitch**9)
                    term2 = -sig_6 / (3.0 * rswitch**3)

                    if (rswitch == rcutoff): 
                        term3 = 0.
                        term4 = 0.
                    else:
                        term3_1 = -2.0 / 3.0 * (1/rcutoff**3 - 1/rswitch**3)
                        term3_2 = C3_2 / 5.0 * (1/rcutoff**5 - 1/rswitch**5)
                        term3_3 = -C3_3 / 7.0 * (1/rcutoff**7 - 1/rswitch**7)
                        term3_4 = -C3_4 / 9.0 * (1/rcutoff**9 - 1/rswitch**9)
                        term3 = -sig_12 / C3_0 * (term3_1 + term3_2 + term3_3 + term3_4)
                         
                        term4_1 = -2.0 / 3.0 * (rcutoff**3 - rswitch**3)
                        term4_2 = C3_2 * (rcutoff - rswitch)
                        term4_3 = C3_3 * (1./rcutoff - 1./rswitch)
                        term4_4 = C3_4 / 3 * (1./rcutoff**3 - 1./rswitch**3)
                        term4 = -sig_6 / C3_0 * (term4_1 + term4_2 + term4_3 + term4_4)
                    
                    elrc += C0 * eps * (term1 + term2 + term3 + term4)
        return elrc
                

def readGroLog(groLogFile):
    isAverage = False
    isVolumeNext = False
    pars = dict.fromkeys(['volume', 'rswitch', 'rcutoff'])
    for line in groLogFile:
        if "rvdw_switch" in line.split():
            pars['rswitch'] = float(line.split()[2])
            continue
        if "rvdw" in line.split():
            pars['rcutoff'] = float(line.split()[2])
            continue
        if "<====  A V E R A G E S  ====>" in line:
            isAverage = True
            continue
        if isAverage and "Volume" in line:
            isVolumeNext = True
            idxVolume = line.split().index("Volume")
            continue
        if isVolumeNext == True:
             pars['volume'] = float(line.split()[idxVolume])
             isVolumeNext = False
             isAverage = False
    return pars


parser = argparse.ArgumentParser(description='Calculate LJ long-range correction')
#parser.add_argument('-rs', '--rswitch', type=float, dest='rswitch', required=True,
#    help='r_switch in nm')
#parser.add_argument('-rc', '--rcutoff', type=float, dest='rcutoff', required=True,
#    help='r_cutoff in nm')
parser.add_argument('-l', '--log', type=argparse.FileType('r'), required=True,
    help='Gromacs log file, for obtaining average volume.')
parser.add_argument('-d', '--dir', default=os.getcwd(),
    help='directory where MDinfo and SltInfo are put (default is current working dir)')

args = parser.parse_args()
pars = readGroLog(args.log)

MDinfo = open(args.dir + "/MDinfo", 'r')
SltInfo = open(args.dir + "/SltInfo", 'r')
rswitch = pars["rswitch"]
rcutoff = pars["rcutoff"]



lineCounter = 1
for line in MDinfo:
    if lineCounter == 1:
        numTotalType = int(line.split()[1])
#        numSlvType = numTotalType - 1 #there is always only one solute type
    if lineCounter == 2:
        molNum = [int(line.split()[i]) for i in range(len(line.split()))]
    if lineCounter == 3:
        siteNum = [int(line.split()[i]) for i in range(len(line.split()))]
    lineCounter += 1

lineCounter = 1
sltAtoms = []
for line in SltInfo:
    record = [float(line.split()[i]) for i in range(len(line.split())) if i > 1]
    sltAtoms.append(Atom(record[0], record[1:3])) 

solute = Molecule(sltAtoms)

for i in range(1, numTotalType):
    MolPrm = [open("MolPrm"+str(i), 'r')]

solvent = []
for file in MolPrm:
    slvAtoms = []
    for line in file:
        record = [float(line.split()[i]) for i in range(len(line.split())) if i > 1]
        slvAtoms.append(Atom(record[0], record[1:3]))
    solvent.append(Molecule(slvAtoms))

system = System(pars['volume'], solute, solvent, molNum[1:])
print('DIR = ' + args.dir)
print('LOG = ' + args.log.name)
print('rswitch = ' + str(rswitch))
print('rcutoff = ' + str(rcutoff))
print('average volume = ' + str(pars['volume']))
print('-------------------------------')
print('elrc = ' + str(system.getElrc(rswitch, rcutoff)))

