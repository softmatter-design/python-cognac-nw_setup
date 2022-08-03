#
# File: CognacTrajectoryAnalysis.py
#
#   Copyright (c) 2001-2020
#   OCTA Licensing Committee 
#   All right reserved.                     
#
#              2017/4/2 T.Aoyagi, OCTA User's group
#             
#Description: Python script for analysis of trajectory data
#Need: 
#    numpy
#    CognacUtility.dll
#    CognacBasicAnalysis.py
#
#Usage:
#    put this file to the directory where is included in python path 
#
#Class:
#    CognacTrajectoryAnalysis (subclass of CognacBasicAnalysis)
#
#Constructor:
#     CognacTrajectoryAnalysis(udfname)
#        udfname ... UDF file name
#
#Method:
#
#    molMsd(molname,start_record,end_record) ... return mean square displacement of center of mass 
#                    of molecules
#        molname ... molecular name
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, msd, (msd_x,msd_y,msd_z))
#
#    atomMsd(molname,atomIndex,start_record,end_record) ... return mean square displacement of specified atoms
#        molname ... molecular name
#        atomIndex ... list of Index of atom
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, msd, (msd_x,msd_y,msd_z))
#
#    normalCoordinate(molname, p,start_record,end_record) ... return autocorrelation of Normal coordinate, Cp(t)
#        molname ... molecular name
#        p ... p-th mode  (default = 1)
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cp, (Cp_x, Cp_y, Cp_z))
#
#    vectorAutoCorrelation(molname, atomIndex,start_record,end_record) ... return auto correlation of specfied vector
#        molname ... molecular name
#        atomIndex ... list of atom Index. e.g.[0,10]
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cvec, (Cvec_x, Cvec_y, Cvec_z))
#
# This function is a contribution of Dr. M. Malvaldi, Universita' di Pisa. 2004/8/28
# The name of function is changed. 2010/4/16 T.Aoyagi
#    atomVelocityAutoCorrelation(atom_name,start_record,end_record) ... return auto correlation of velocity of atoms of specified atom type
#        atom_name ... atom name
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cvel, (Cvel_x, Cvel_y, Cvel_z))
#
#    molVelocityAutoCorrelation(mol_name,start_record,end_record) ... return auto correlation of velocity of center of mass of molecules of specified molecular name
#        mol_name ... molecular name
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cvel, (Cvel_x, Cvel_y, Cvel_z))
#
#    forceAutoCorrelation(atom_name,start_record,end_record) ... return auto correlation of force of atoms of specified atom type
#        atom_name ... atom name
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cforce, (Cforce_x, Cforce_y, Cforce_z))
#
#    stressAutoCorrelation(start_record,end_record) ... return auto correlation of stress of the system
#        start_record ... start record for sampling (default=first)
#        end_record ... end record for sampling (default=last)
#        return object ... list of (time, Cpress, (Cforce_xx, Cforce_yy, Cforce_zz), (Cforce_xy, Cforce_yz, Cforce_zx))
#
import numpy as np
import CognacUtility as cu
from CognacBasicAnalysis import CognacBasicAnalysis

class CognacTrajectoryAnalysis (CognacBasicAnalysis):
    def __init__(self,udfname):
        CognacBasicAnalysis.__init__(self,udfname)
        # setup timeRecord[] and timeList[]
        startTime = 0.0        # output UDF always starts at time=0
        self.timeRecord = [0]
        self.timeList = [startTime]
        for i in range(0,self.totalRecord()):
            self.jump(i)
            time = self.get("Time")
            if time == startTime:
                self.timeRecord[0] = i
            else:
                self.timeRecord.append(i)
                self.timeList.append(time - startTime)

    #----- utility functions -----

    def _recList(self, start_record, end_record):
        '''list of records to be used for analysis

        records in the range start <= rec <= end will be used

        Args:
            start_record (str): 'first' or str(record_number)
            end_record (str): 'last' or str(record_number)

        Returns:
            list of record numbers
        '''
        if start_record == "first" and end_record == "last":
            return self.timeRecord
        else:
            if start_record == "first":
                beginRec = 0
            else:
                beginRec = eval(start_record)    
            if end_record == "last":
                endRec = self.totalRecord()
            else:
                endRec = eval(end_record) + 1
            return range(beginRec,endRec)

    def _molList(self, molname, suffix=None):
        '''list of molecules having the specifiled molname

        Args:
            molname (str): name of the molecules to be selected
            suffix (str): suffix to be added to the 'key'

        Returns:
            list [(molIndex, key), ...]
            where key = 'molname:molIndex:suffix'
        '''
        key = '{}:{}'
        if suffix != None:
            key += ':' + suffix
        molList = []
        for m in range(self.totalmol):
            if self.get("Set_of_Molecules.molecule[].Mol_Name",[m]) == molname:
                molList.append((m, key.format(molname, m)))
        if len(molList) == 0:
            print('no molecule with Mol_Name: ' + molname)
        return molList

    # def _atomList(self, atomname):
    #     '''list of atoms having the specified atomname

    #     Args:
    #         atomname (str): name of the atoms to be selected

    #     Returns:
    #         list [ (molIndex, atomIndex, key), ...],
    #         where key = 'atomname:molIndex:atomIndex'
    #     '''
    #     udfpath = 'Set_of_Molecules.molecule[].atom[]'
    #     atomList = []
    #     for m in range(self.totalmol):
    #         for a in range(self.size(udfpath, [m])):
    #             if self.get(udfpath+'.Atom_Name', [m,a]) == atomname:
    #                 atomList.append((m, a, '{}:{}:{}'.format(atomname, m, a)))
    #     if len(atomList) == 0:
    #         print('no atom with Atom_Name: ' + atomname)
    #     return atomList

    # def _mol_atomList(self, molname, atomList=None, suffix=None):
    #     '''list of selected atoms

    #     Args:
    #         molname (str): select molecules with this name.
    #         atomList (None, list of atomIndex, or a single atomIndex):
    #                 None: select all the atoms in the molecule
    #                 list of atomIndex: selct atoms in the list
    #                 single atomIndex: select the atom

    #     Returns:
    #         list [(molIndex, atomIndex, key), ...],
    #         where key = 'molname:molIndex:atomIndex:suffix'
    #     '''
    #     key = '{}:{}:{}'
    #     if suffix != None:
    #         key += ':' + suffix
    #     m_aList = []
    #     for m in range(self.totalmol):
    #         if self.get("Set_of_Molecules.molecule[].Mol_Name",[m]) == molname:
    #             if atomList==None:
    #                 natom = self.size("Set_of_Molecules.molecule[].atom[]",[m])
    #                 atoms = range(natom)
    #             elif type(atomList) == type([]):
    #                 atoms = atomList
    #             else:
    #                 atoms = [atomList]
    #             for a in atoms:
    #                     m_aList.append((m, a, key.format(molname, m, a)))
    #     if len(m_aList) == 0:
    #         print('no molecule with Mol_Name: ' + molname)
    #     return m_aList

    def _resultsList(self, vec):
        '''normalize vec and create the results list

        Arg:
            vec: [ [x0,y0,z0], [x1,y1,z1], ...]

        Returns:
            list [ (t0, s0, (xn0,yn0,zn0)), (t1, s1, (xn1,yn1,zn1)), ... ]
            where
                ti = dt*i,
                si = (xi+yi+zi)/S    (S = x0+y0+z0),
                xni = xi/x0, yni = yi/y0, zni = zi/z0.
        '''
        S = np.sum(vec[0])
        x0, y0, z0 = vec[0]
        results = [ (self.timeList[i], np.sum(vec[i])/S,
                        (vec[i][0]/x0, vec[i][1]/y0, vec[i][2]/z0))
                        for i in range(len(vec)) ]
        return results

    #----- public methods -----

    # def molMsd(self, molname, start_record="first", end_record="last",
    #                                                 cancel_trans=False):
    #     # example usage of timer()
    #     self.timer_off()        # use timer_on() to enable timing

    #     molList = self._molList(molname, 'com')
    #     if len(molList) == 0:
    #         return 0

    #     self.timer('create molList')    # do nothing if timer is off

    #     cu.clearVectorMap()
    #     for rec in self._recList(start_record, end_record):
    #         self.jump(rec)
    #         if cancel_trans:
    #             Gx, Gy, Gz = self.system_centerofmass()
    #         for m, key in molList:
    #             x, y, z = self.centerofmass(m)
    #             if cancel_trans:
    #                 x -= Gx
    #                 y -= Gy
    #                 z -= Gz
    #             cu.pushVector(key, (x,y,z))

    #     self.timer('calculate centers of mass')

    #     msdAve = cu.msd()

    #     self.timer('calculate MSD')

    #     return [ (self.timeList[i], np.sum(msdAve[i]), tuple(msdAve[i]))
    #                 for i in range(1,len(msdAve)) ]

    # def atomMsd(self, molname, atomIndex=None,
    #             start_record="first", end_record="last", cancel_trans=False):
    #     m_aList = self._mol_atomList(molname, atomIndex)
    #     if len(m_aList) == 0:
    #         return 0

    #     cu.clearVectorMap()
    #     for rec in self._recList(start_record, end_record):
    #         self.jump(rec)
    #         if cancel_trans:
    #             Gx, Gy, Gz = self.system_centerofmass()
    #         for m, a, key in m_aList:
    #             x, y, z = self.position([m,a])
    #             if cancel_trans:
    #                 x -= Gx
    #                 y -= Gy
    #                 z -= Gz
    #             cu.pushVector(key, (x,y,z))

    #     msdAve = cu.msd()

    #     return [ (self.timeList[i], np.sum(msdAve[i]), tuple(msdAve[i]))
    #                 for i in range(1,len(msdAve)) ]

    def normalCoordinate(self, molname, p=1,
                                    start_record="first", end_record="last"):
        molList = self._molList(molname, 'Xp')
        
        cu.clearVectorMap()
        for rec in self._recList(start_record, end_record):
            self.jump(rec)
            for m, key in molList:
                cu.pushVector(key, self.Xp(m,p))

        CpAve = cu.vectorCorrelation()
        return self._resultsList(CpAve)

    def Xp(self, molIndex, p=1):
        '''normal mode of a molecule

        Args:
            molIndex: int: specify the molecule
            p: int: mode number (1, 2, ..., numAtoms-1)

        Returns:
            tuple (Xpx, Xpy, Xpz)
        '''
        pos = np.array(self.position([molIndex]))
        print(pos)

        N = len(pos)
        k = np.pi*p/N
        xp = np.zeros(3)
        for n in range(N):
            xp += np.cos(k*(n+0.5))*pos[n]
        return tuple(xp/N)

    def position(self, mol_atom):
        '''position of an atom, or list of positions of atoms in a molecule

        position([molIndex, atomIndex])
            Returns: tuple (x,y,z): position of the atom

        position([molIndex])
            Returns: ([x,y,z], [x,y,z], ...): positions of atoms in the molecule
        '''
        return tuple(self.get("Structure.Position.mol[].atom[]", mol_atom))

#     def vectorAutoCorrelation(self, molname, atomIndex,
#                                     start_record="first", end_record="last"):
#         molList = self._molList(molname, 'vec')
#         if len(molList) == 0:
#             return 0

#         cu.clearVectorMap()
#         for rec in self._recList(start_record, end_record):
#             self.jump(rec)
#             for m, key in molList:
#                 cu.pushVector(key,
#                             self.vector([m,atomIndex[0]],[m,atomIndex[1]]))

#         vecCorrAve = cu.vectorCorrelation()
#         return self._resultsList(vecCorrAve)

# # This function is a contribution of Dr. M. Malvaldi, Universita' di Pisa. 2004/8/28
# # Atom Name is used to specify the set of atoms. 2009/5/1 T.Aoyagi
# # The name of function is changed for the extention of molecular velocity auto correlation. 
# #  2010/3/19 T.Aoyagi
# #
#     def atomVelocityAutoCorrelation(self, atom_name,
#                                     start_record="first", end_record="last"):
#         atomList = self._atomList(atom_name)
#         if len(atomList) == 0:
#             return 0

#         cu.clearVectorMap()
#         for rec in self._recList(start_record, end_record):
#             self.jump(rec)
#             for m, a, key in atomList:
#                 cu.pushVector(key, self.velocity([m,a]))

#         velocityCorrAve = cu.vectorCorrelation()
#         return self._resultsList(velocityCorrAve)

#     def molVelocityAutoCorrelation(self, mol_name,
#                                     start_record="first", end_record="last"):
#         molList = self._molList(mol_name)
#         if len(molList) == 0:
#             return 0

#         cu.clearVectorMap()
#         for rec in self._recList(start_record, end_record):
#             self.jump(rec)
#             for m, key in molList:
#                 cu.pushVector(key, self.transvelofmol(m))

#         velocityCorrAve = cu.vectorCorrelation()
#         return self._resultsList(velocityCorrAve)

#     def forceAutoCorrelation(self, atom_name,
#                                     start_record="first", end_record="last"):
#         atomList = self._atomList(atom_name)
#         if len(atomList) == 0:
#             return 0

#         cu.clearVectorMap()
#         for rec in self._recList(start_record, end_record):
#             self.jump(rec)
#             for m, a, key in atomList:
#                 cu.pushVector(key, self.force([m,a]))

#         forceCorrAve = cu.vectorCorrelation()
#         return self._resultsList(forceCorrAve)

#     def stressAutoCorrelation(self, start_record="first", end_record="last"):
#         stress = []
#         for rec in self._recList(start_record, end_record):
#             self.jump(rec)
#             if self.get('Steps') != 0:
#                 stress.append(
#                         self.get("Statistics_Data.Stress.Total.Instantaneous"))

#         scorr = self.stressCorrelation(stress)
#         nRec = len(scorr)

#         for i in range(1,nRec):
#             for k in range(0,7):
#                 scorr[i][k] /= scorr[0][k]

#         for k in range(0,7):
#             scorr[0][k] = 1.0
            
#         return [ (self.timeList[i], scorr[i][6],
#                         (scorr[i][0],scorr[i][1],scorr[i][2]),
#                         (scorr[i][3],scorr[i][4],scorr[i][5])
#                  ) for i in range(nRec) ]

#     def stressCorrelation(self,stress):
#         nrec = len(stress)

#         press = [ (stress[i][0]+stress[i][1]+stress[i][2])/3.0
#                         for i in range(nrec) ]

#         corr = [ [0.0]*7 for i in range(nrec) ]
#         for i in range(0,nrec):
#             for j in range(i,nrec):
#                 for k in range(0,6):
#                     corr[j-i][k] += stress[i][k]*stress[j][k]
#                 corr[j-i][6] += press[i]*press[j]
#         for i in range(0,nrec):
#             for j in range(0,7):
#                 corr[i][j] /= float(nrec-i)

#         return corr
