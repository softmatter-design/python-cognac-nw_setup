import numpy as np
import CognacUtility as cu
from CognacBasicAnalysis import CognacBasicAnalysis

class XpCalc (CognacBasicAnalysis):
    def __init__(self,udfname):
        CognacBasicAnalysis.__init__(self,udfname)
        # # setup timeRecord[] and timeList[]
        # startTime = 0.0        # output UDF always starts at time=0
        # self.timeRecord = [0]
        # self.timeList = [startTime]
        # for i in range(0,self.totalRecord()):
        #     self.jump(i)
        #     time = self.get("Time")
        #     if time == startTime:
        #         self.timeRecord[0] = i
        #     else:
        #         self.timeRecord.append(i)
        #         self.timeList.append(time - startTime)

    #----- utility functions -----
    def time(self):
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



    def normalCoordinate(self, molname, p=1,
                                    start_record="first", end_record="last"):
        self.time()
        molList = self._molList(molname, 'Xp')
        print(molList)
        
        cu.clearVectorMap()
        print(self._recList(start_record, end_record))
        for rec in self._recList(start_record, end_record):
            self.jump(rec)
            pos_list = tuple(self.get("Structure.Position.mol[].atom[]"))
            for m, key in molList:
                pos = np.array(pos_list[m])
                cu.pushVector(key, self.Xp(pos,p))

        CpAve = cu.vectorCorrelation()
        return self._resultsList(CpAve)

    def Xp(self, pos, p=1):
        N = len(pos)
        k = np.pi*p/N
        xp = np.zeros(3)
        for n in range(N):
            xp += np.cos(k*(n+0.5))*pos[n]
        return tuple(xp/N)

    # def position(self, mol_atom):
    #     '''position of an atom, or list of positions of atoms in a molecule

    #     position([molIndex, atomIndex])
    #         Returns: tuple (x,y,z): position of the atom

    #     position([molIndex])
    #         Returns: ([x,y,z], [x,y,z], ...): positions of atoms in the molecule
    #     '''
    #     return tuple(self.get("Structure.Position.mol[].atom[]", mol_atom))
