import numpy as np
import constants
import sys
#
class polarizable_embedding:
    #
    """Polarizable embedding class object"""
    #
    def __init__(self, force_field = 'fq', atomtypes = [], chi = [], eta = [], alpha = [], \
                 Rq = [], Rmu = [], pqeq = False):
        #
        """Initialization of the polarizable embedding class object"""
        """force_field = fq,fqfmu
           atomtypes = list of str with the name of the atomtypes 
           chi       = list of electronegativities
           eta       = list of chemical hardnesses
           alpha     = list of atomic polarizabilities
           Rq        = list of charge radii (only if pqeq)
           Rmu       = list of dipoles radii (only if pqeq)
           pqeq      = if pqeq computation
        """
        #
        # Sanity checks and initialization
        #
        self.check_parameters(force_field, atomtypes, chi, eta, alpha, Rq, Rmu, pqeq)
        #
    #
    def print_info(self, file_=''):
        #
        # Print info of the forcefield, in case you pass a file, write in there
        #
        original_stdout = sys.stdout 
        if (file_ != ''):
            sys.stdout = file_
        else:
            print('---polarizable force field---')
        #
        print(' force_field : ' + self.force_field)
        infostring = ''
        for ii,i in enumerate(self.atomtypes):
            infostring += "'" + i + "'"
            if ii < len(self.atomtypes) -1:
                infostring += ' , '
        print(' atomtypes   : ' + infostring)
        #
        infostring = ''
        for ii,i in enumerate(self.chi):
            infostring += str(i)
            if ii < len(self.chi) -1:
                infostring += ' , '
        print(' chi         : ' + infostring)
        #
        infostring = ''
        for ii,i in enumerate(self.eta):
            infostring += str(i)
            if ii < len(self.eta) -1:
                infostring += ' , '
        print(' eta         : ' + infostring)
        #
        if (self.force_field == 'fqfmu' or self.force_field == 'fqfmu_pqeq'):
            infostring = ''
            for ii,i in enumerate(self.alpha):
                infostring += str(i)
                if ii < len(self.alpha) -1:
                    infostring += ' , '
            print(' alpha       : ' + infostring)
        #
        if (self.force_field == 'fq_pqeq' or self.force_field == 'fqfmu_pqeq'):
            infostring = ''
            for ii,i in enumerate(self.Rq):
                infostring += str(i)
                if ii < len(self.Rq) -1:
                    infostring += ' , '
            print(' Rq          : ' + infostring)
        #
        if (self.force_field == 'fqfmu_pqeq'):
            infostring = ''
            for ii,i in enumerate(self.Rmu):
                infostring += str(i)
                if ii < len(self.Rmu) -1:
                    infostring += ' , '
            print(' Rmu         : ' + infostring)
        #
        sys.stdout = original_stdout
        #print('pqeq:       : ' + str(self.pqeq))
        #

    #
    #
    #
    def check_parameters(self,force_field = '',  atomtypes = [], chi = [], eta = [], alpha = [], \
                         Rq = [], Rmu = [], pqeq = True):
        #
        """Procedure to check the given inputs to the class"""
        #
        # Force Field
        #
        if (force_field != 'fq' and force_field != 'fqfmu'):
            if (force_field == 'fq_pqeq' and pqeq == True):
                self.force_field = 'fq_pqeq'
                self.pqeq = True
            elif(force_field == 'fq_pqeq' and pqeq == False):
                print('ERROR: forcefield fq_pqeq and pqeq keyword set to False')
                sys.exit()
            elif(force_field == 'fqfmu_pqeq' and pqeq == True):
                self.force_field = 'fqfmu_pqeq'
                self.pqeq = True
            elif(force_field == 'fqfmu_pqeq' and pqeq == False):
                print('ERROR: forcefield fqfmu_pqeq and pqeq keyword set to False')
                sys.exit()
            else:
                print('ERROR: forcefield ' + str(force_field) + ' not recognized')
                sys.exit()
        else:
            self.force_field = force_field
        #
        # Atomtypes
        #
        if (type(atomtypes) != list):
            atomtypes = [atomtypes]
        for i in atomtypes:
            if (type(i) != str):
                print('ERROR: polarizable embedding with non-character atomtypes provided')
                sys.exit()
        self.atomtypes = atomtypes
        #
        if (len(self.atomtypes) != len(set(atomtypes))):
            print('ERROR: polarizable embedding with repeated atomtypes')
            sys.exit()
        #
        # Chi
        #
        if (type(chi) != list):
            chi = [chi]
        for i in chi:
            if (type(i) != float):
                print('ERROR: polarizable embedding with non-float electronegativities provided')
                sys.exit()
        self.chi = chi
        #
        # Eta
        #
        if (type(eta) != list):
            eta = [eta]
        for i in eta:
            if (type(i) != float):
                print('ERROR: polarizable embedding with non-float chemical hardness provided')
                sys.exit()
        self.eta = eta
        #
        # pqeq
        #
        if type(pqeq) != bool:
            print('ERROR: polarizable embedding with non-boolean pqeq option provided')
            sys.exit()
        self.pqeq = pqeq
        #
        # Rq
        #
        if (self.pqeq):
            if (type(Rq) != list):
                Rq = [Rq]
            for i in Rq:
                if (type(i) != float):
                    print('ERROR: polarizable embedding with non-float Rq provided')
                    sys.exit()
            self.Rq = Rq
        #
        # Alpha
        #
        if (self.force_field == 'fqfmu' or self.force_field == 'fqfmu_pqeq'):
            if (type(alpha) != list):
                alpha = [alpha]
            for i in alpha:
                if (type(i) != float):
                    print('ERROR: polarizable embedding with non-float polarizabilities provided')
                    sys.exit()
            self.alpha = alpha
            #
            # Rmu
            #
            if (self.pqeq):
                if (type(Rmu) != list):
                    Rmu = [Rmu]
                for i in Rmu:
                    if (type(i) != float):
                        print('ERROR: polarizable embedding with non-float Rmu provided')
                        sys.exit()
                self.Rmu = Rmu
        #
        # Check that the length of the different lists are the same and initialize uninitialize variables
        #
        if (self.force_field == 'fq'):
            if(len(self.atomtypes) != len(self.chi)   or \
               len(self.atomtypes) != len(self.eta)):
               print('ERROR: fq uncorrectly initilized (different length of the lists)')
               sys.exit()
            else:
                self.alpha = ['None']
                self.Rq = ['None']
                self.Rmu = ['None']
        #
        elif (self.force_field == 'fq_pqeq'):
            if(len(self.atomtypes) != len(self.chi)   or \
               len(self.atomtypes) != len(self.eta)   or \
               len(self.atomtypes) != len(self.Rq)):
               print('ERROR: fq_pqeq uncorrectly initilized (different length of the lists)')
               sys.exit()
            else:
                self.alpha = ['None']
                self.Rmu = ['None']
        #
        elif (self.force_field == 'fqfmu'):
            if(len(self.atomtypes) != len(self.chi)   or \
               len(self.atomtypes) != len(self.eta)   or \
               len(self.atomtypes) != len(self.alpha)):
               print('ERROR: fqfmu uncorrectly initilized (different length of the lists)')
               sys.exit()
            else:
                self.Rq = ['None']
                self.Rmu = ['None']
        #
        elif (self.force_field == 'fqfmu_pqeq'):
            if(len(self.atomtypes) != len(self.chi)   or \
               len(self.atomtypes) != len(self.eta)   or \
               len(self.atomtypes) != len(self.alpha) or \
               len(self.atomtypes) != len(self.Rq)    or \
               len(self.atomtypes) != len(self.Rmu)):
               print('ERROR: fqfmu_pqeq uncorrectly initilized (different length of the lists)')
               sys.exit()
        else:
            print('ERROR: polarizable embedding unrecognized embedding.')
        #
        return
           



