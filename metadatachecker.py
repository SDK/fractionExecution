__author__ = 'sagonzal'
from sacm import *
import math as mymath
import itertools

class AsdmCheck:
    #uid = ''
    #asdmDict = dict()
    #check = dict()
    #toc = ''
    #main = ''
    #antennas = ''
    #source = ''
    #scan = ''
    #field = ''
    #syscal = ''


    def __init__(self):
        self.uid = ''
        self.asdmDict = {}
        self.check = {}
        self.toc = ''
        self.main = ''
        self.antennas = ''
        self.source = ''
        self.scan = ''
        self.field = ''
        self.syscal = ''


    def setUID(self,uid=None):
        self.uid = uid
        try:
            asdmList, asdm , self.toc = getASDM(self.uid)
            for i in asdmList:
                self.asdmDict[i[0].strip()] = i[3].strip()
            self.main = getMain(self.asdmDict['Main'])
            self.antennas = getAntennas(self.asdmDict['Antenna'])
            self.source = getSource(self.asdmDict['Source'])
            self.scan = getScan(self.asdmDict['Scan'])
            self.field = getField(self.asdmDict['Field'])
            self.syscal = getSysCal(self.asdmDict['SysCal'])


            return True
        except Exception as e:
            print 'There is a problem with the uid: ', uid
            print e


    def isNullState(self):
        """
        Checks the Main.xml table for "null" states
        Sets in the self.check dictionary the value 'NullState' with True or False
        True: There is no state with null values
        False: There is at least one state in the Main.xml table with null state.
        :return:
        """
        try:

            self.main['null'] = self.main.apply(lambda x: True if 'null' in x['stateId'] else False, axis = 1)
            if len(self.main['null'].unique()) > 1:
                self.check['NullState'] = False
            else:
                self.check['NullState'] = True
        except Exception as e:
            print e
            return False
        return True


    def isValidUID(self):
        """

        :return:
        """
        import re
        regex = re.compile("^uid\:\/\/A00.\/X[a-zA-Z0-9]+\/X[a-zA-Z0-9]+")
        try:
            for k,v in self.asdmDict.iteritems():
                if '/X0/X0/X0' in v:
                    self.check['ValidUID'] = False
                    return True
                if regex.match(v) is None:
                    self.check['ValidUID'] = False
                    return True
            self.check['ValidUID'] = True
            return True
        except Exception as e:
            print e
            return False

    def isSyscaltimestamp(self):
        try:
            dfa = self.syscal[['spectralWindowId','antennaId','timeInterval']]
            spw = dfa.groupby('spectralWindowId')
            for name,group in spw:
                df = group
                df['time'] = df.apply(lambda x: int(x['timeInterval'].strip().split(' ')[0]) / 1000000000.0, axis = 1 )
                df['interval'] = df.apply(lambda x: int(x['timeInterval'].strip().split(' ')[1]) / 1000000000.0, axis = 1 )
                df['timestamp'] = df.apply(lambda x: x.time - x.interval / 2, axis = 1)
                t0 = 86400.0 * mymath.floor(df.timestamp.min() / 86400.0)
                df['utimes'] = df.apply(lambda x: x['time'] - t0, axis =1)
                nT = df.utimes.nunique()
                df['utimestamp'] = df.apply(lambda x: mymath.floor(x['timestamp']) - t0 , axis =1)
                nTS = df.utimestamp.nunique()
                #print name,nT, nTS
                #print (group)
                if nT != nTS:
                    self.check['SysCalTimes'] = True
                    return True
            self.check['SysCalTimes'] = False
            return True
        except Exception as e:
            return False


    def iscsv2555(self):
        try:
            src = self.source[['sourceId', 'sourceName']]
            src['sourceName1'] = src.apply(lambda x: x['sourceName'].strip(), axis = 1)
            src = src.drop_duplicates()
            fld = self.field[['sourceId', 'fieldName']]
            fld['fieldName1'] = fld.apply(lambda x: x['fieldName'].strip(), axis = 1)
            fld = fld.drop_duplicates()
            a = pd.merge(src,fld,left_on = 'sourceId',right_on='sourceId',how='outer')
            a['csv2555'] = a.apply(lambda x: True if x['sourceName1'] == x['fieldName1'] else False, axis = 1)
            if a['csv2555'].nunique() == 1 and a['csv2555'].unique()[0] is True:
                self.check['CSV2555'] = True
            else:
                self.check['CSV2555'] = False
            return True
        except Exception as e:
            return False


    def isfixplanets(self):
        try:
            df = self.source[['sourceName','direction']].drop_duplicates()
            df['coordinate'] = df.apply(lambda x: True if float(arrayParser(x['direction'].strip() , 1 )[0]) == 0.0 else False , axis = 1)
            df['coordinate2'] = df.apply(lambda x: True if float(arrayParser(x['direction'].strip() , 1 )[1]) == 0.0 else False , axis = 1)
            if df['coordinate'].unique()[0] is False and df['coordinate'].nunique() == 1 and df['coordinate2'].unique()[0] is False and df['coordinate2'].nunique() == 1:
                self.check['FixPlanets'] = True
            else:
                self.check['FixPlanets'] = False
            return True
        except Exception as e:
            return False

    def ict4871(self):
        try:
            nant = self.antennas.antennaId.nunique()
            sys_nant = self.syscal.antennaId.nunique()
            sc = self.scan[['scanNumber','startTime','scanIntent']]
            problem_list = list()
            if nant == sys_nant:
                df = self.syscal[['timeInterval','antennaId','spectralWindowId']]
                df['start'] = df.apply(lambda x: int(x['timeInterval'].strip().split(' ')[0]) - int(x['timeInterval'].strip().split(' ')[1])/2, axis = 1)
                df2 = pd.merge (df,sc, left_on='start',right_on='startTime',copy=False,how='inner')
                antlist = [x.strip(' ') for x in self.antennas.antennaId.unique().tolist()]
                df3 = df2.groupby(['antennaId','spectralWindowId','scanNumber'])
                spw_list = self.syscal.spectralWindowId.unique().tolist()
                scan_list = df2.scanNumber.unique().tolist()
                fu = list(itertools.product(antlist,spw_list,scan_list))
                for i in fu:
                    try:
                        df3.groups[i]
                    except KeyError as k:
                        problem_list.append(i)

                if len(problem_list) > 0:
                    df = pd.DataFrame(problem_list, columns= ['antennaId','spectralWindowId','scanNumber'])
                    popular_scan = df.scanNumber.mode().values[0]
                    ant = df.antennaId.nunique()
                    self.check['MissingTsys'] = False
                    self.check['MissingTsys_explain'] = 'Scan: '+str(popular_scan)+' Antennas Affected: '+str(ant)
                else:
                    self.check['MissingTsys'] = True

            else:
                self.check['AntennasMissing'] = "Number of antennas are diferent between Antenna.xml ("+nant+") and SysCal.xml ("+sys_nant+")"
                self.check['MissingTsys'] = False
                self.check['MissingTsys_explain'] = 'Some antenna is completly missing from Syscal.xml table'


        except Exception as e:
            print e
            return False

    def doCheck(self):
        try:
            self.iscsv2555()
            self.isSyscaltimestamp()
            self.isValidUID()
            self.isfixplanets()
            self.isNullState()
            self.ict4871()
            return True
        except Exception as e:
            print e
            return False

