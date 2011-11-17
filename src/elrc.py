#!/home1/kmtu/local/Python-3.2/bin/python3.2

import argparse
import os
import math
import sys

JOULE_PER_CAL = 4.184

class Atom:
    """An atom class"""
    def __init__(self, charge, lj):
        self.charge = charge
        self.epsilon = lj[0]
        self.sigma = lj[1]

    def __str__(self):
        return ('charge' + str(self.charge) + '\n' 
                'epsilon' + str(self.epsilon) + '\n'
                'sigma' + str(self.sigma))

    def display(self):
        print('charge:' + str(self.charge) + ' ' 
              'epsilon:' + str(self.epsilon) + ' '
              'sigma:' + str(self.sigma))

#    def __eq__(self, other):
#        return (self.charge == other.charge and 
#                self.epsilon == other.epsilon and
#                self.sigma == other.sigma)

class Molecule:
    """A molecule class"""
    def __init__(self, atoms):
        self.sites = atoms

    def __str__(self):
        return ('sites ' + str(len(self.sites)) + '\n' +
                 str(self.sites))

    def display(self):
        print('*** start print molecule ***')
        print('sites', str(len(self.sites)))
        for atom in self.sites:
            atom.display()
        print('*** end print molecule **')
        
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

    def __str__(self):
        return ('volume ' + str(self.volume) + '\n'
                'solute\n' +
                str(self.solute) + '\n'
                'solvent ' + str(len(self.solvent)) + '\n' +
                str(self.solvent) + '\n' +
                'solvent number' + '\n' +
                str(self.slvNum) + '\n' +
                'solvent numDensity\n' + 
                str(self.slvNumDensity))
     
    def display(self):
        print('***** start print system *****')
        print('volume ', self.volume)
        print('solute')
        self.solute.display()
        print('solvent ', len(self.solvent))
        for slv in self.solvent:
            print('number', self.slvNum[slv])
            print('numDensity', self.slvNumDensity[slv])
            slv.display()
        print('***** end print system *****')

    #obsolete (incomplete)
    def getElrc_Karino(self, rs, rc):
        elrc = {}
        for slv in self.solvent:
            for slvSite in slv.sites:
                for sltSite in self.solute.sites:
                    ep = math.sqrt(slvSite.epsilon * sltSite.epsilon)


    def getElrc(self, rs, rc):
        elrc = {}
        if rs != rc:
            C3_0 = (rc**2 - rs**2)**3
            C3_2 = 3.0 * rs**2 + 3.0 * rc**2
            C3_3 = 6.0 * rc**2 * rs**2
            C3_4 = rc**6 - 3.0 * rc**4 * rs**2

        for slv in self.solvent:
            C0 = 16.0 * math.pi * self.slvNumDensity[slv]
            elrc[slv] = 0.0
            for slvSite in slv.sites:
                for sltSite in self.solute.sites:
                    eps = math.sqrt(slvSite.epsilon * sltSite.epsilon)
                    sig_2 = slvSite.sigma * sltSite.sigma    
                    sig_6 = sig_2**3
                    sig_12 = sig_6 * sig_6 

                    term1 = sig_12 / (9.0 * rs**9)
                    term2 = -sig_6 / (3.0 * rs**3)

                    if (rs == rc): 
                        term3 = 0.
                        term4 = 0.
                    else:
                        term3_1 = -2.0 / 3.0 * (1./rc**3 - 1./rs**3)
                        term3_2 = C3_2 / 5.0 * (1./rc**5 - 1./rs**5)
                        term3_3 = -C3_3 / 7.0 * (1./rc**7 - 1./rs**7)
                        term3_4 = -C3_4 / 9.0 * (1./rc**9 - 1./rs**9)
                        term3 = -sig_12 / C3_0 * (term3_1 + term3_2 + term3_3 + term3_4)
                         
                        term4_1 = -2.0 / 3.0 * (rc**3 - rs**3)
                        term4_2 = C3_2 * (rc - rs)
                        term4_3 = C3_3 * (1./rc - 1./rs)
                        term4_4 = C3_4 / 3.0 * (1./rc**3 - 1./rs**3)
                        term4 = -sig_6 / C3_0 * (term4_1 + term4_2 + term4_3 + term4_4)
                    
                    elrc[slv] += eps * (term1 + term2 + term3 + term4)
            elrc[slv] *= C0 / JOULE_PER_CAL
        elrc['total'] = sum(list(elrc.values())) 
        return elrc
                

def readGroLog(groLogFile, volume):
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
    if pars['volume'] == None:
        if volume == None:
            sys.exit("Unable to find volume in file: " + str(groLogFile.name) + '\n' +
                     "If it is an NVT simulation, please input the volume manually, with -v option")
        else:
            pars['volume'] = volume
    else:
        if volume != None:
            print("Warning: two sources of volume were provided.")    
            print("         I choose the value of -v option only: ", volume)
            pars['volume'] = volume
    return pars


parser = argparse.ArgumentParser(description='Calculate LJ long-range correction')
parser.add_argument('-l', '--log', type=argparse.FileType('r'), required=True,
    help='Gromacs log file, for obtaining average volume.')
parser.add_argument('-d', '--dir', default=os.getcwd(),
    help='directory where MDinfo and SltInfo are put (default is current working dir)')
parser.add_argument('-v', '--volume', type=float, default=None, 
    help='Volume info for NVT simulations')

args = parser.parse_args()
pars = readGroLog(args.log, args.volume)

MDinfo = open(args.dir + "/MDinfo", 'r')
SltInfo = open(args.dir + "/SltInfo", 'r')



lineCounter = 1
for line in MDinfo:
    if lineCounter == 1:
        numTotalType = int(line.split()[1])
    if lineCounter == 2:
        molNum = [int(i) for i in line.split()]
    if lineCounter == 3:
        siteNum = [int(i) for i in line.split()]
    lineCounter += 1

lineCounter = 1
sltAtoms = []
for line in SltInfo:
    record = [float(line.split()[i]) for i in range(len(line.split())) if i > 1]
    sltAtoms.append(Atom(record[0], record[1:3])) 

solute = Molecule(sltAtoms)

MolPrm = [open(args.dir+"/MolPrm"+str(i), 'r') for i in range(1,numTotalType)]
solvent = []
for file in MolPrm:
    slvAtoms = []
    for line in file:
        record = [float(line.split()[i]) for i in range(len(line.split())) if i > 1]
        slvAtoms.append(Atom(record[0], record[1:3]))
    solvent.append(Molecule(slvAtoms))

system = System(pars['volume'], solute, solvent, molNum[1:])
elrc = system.getElrc(pars['rswitch'], pars['rcutoff'])

print()
print('dir = ' + args.dir)
print('log = ' + args.log.name)
print('rswitch = ' + str(pars['rswitch']))
print('rcutoff = ' + str(pars['rcutoff']))
print('average volume = ' + str(pars['volume']))
print('number density (nm^-3):')
for i in range(len(system.solvent)):
    print('  solvent %i = ' % int(i+1), system.slvNumDensity[system.solvent[i]])
print('-------------------------------')
print('elrc (kcal/mol): ')
for i in range(len(system.solvent)):
    print('  solvent %i = ' % int(i+1), elrc[system.solvent[i]])
print('  total = ', elrc['total'])
print()

#debug
#print("\n\n**** debug ****")
#system.display()

