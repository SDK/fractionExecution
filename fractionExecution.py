__author__ = 'sagonzal'
from metadatachecker import *
import pandas as pd
import sys
import Sensitivity.sensitivity as sensitivity
from waitress import serve
from pyramid.config import Configurator
##############################

#configDir = '/home/sagonzal/druva/PycharmProjects/MetaData/Sensitivity'
configDir = '/home/sagonzal/PycharmProjects/fractionExecution/Sensitivity'

def fractionExecution(request):

    try:
        asdmUID = request.params.get('uid',None)
        asdm = AsdmCheck()
        asdm.setUID(asdmUID)
        sb = getSBSummary(asdm.asdmDict['SBSummary'])
    except Exception as e:
        print e
        sys.exit(1)
    sbUID = sb.values[0][0]
    sbtype = getSBType(sbUID)

    c3 = datetime.datetime(2015,10,1,0,0,0)
    c2 = datetime.datetime(2014,6,3,0,0,0)
    c1 = datetime.datetime(2013,10,12,0,0,0)

    antennas4cycle = {
        'TWELVE-M':{
            'c1': {
                'NAntOT' : 32.
            },
            'c2': {
                'NAntOT' : 34.
            },
            'c3':{
                'NAntOT' : 36.
            },
            'c4': {
                'NAntOT' : 40.
            }
        },
        'SEVEN-M':{
            'c1': {
                'NAntOT' : 9.
            },
            'c2': {
                'NAntOT' : 9.
            },
            'c3':{
                'NAntOT' : 10.
            },
            'c4': {
                'NAntOT' : 10.
            }

        },
        'TP-Array':{
            'c1': {
                'NAntOT' : 2.
            },
            'c2': {
                'NAntOT' : 2.
            },
            'c3':{
                'NAntOT' : 3.
            },
            'c4': {
                'NAntOT' : 3.
            }
        }
    }

    if asdm.toc >= c3:
        cycle = 'c3'
    elif asdm.toc >= c2:
        cycle = 'c2'
    elif asdm.toc >= c1:
        cycle = 'c1'
    else:
        cycle = 'c0'



    ############################

    freq = float(sb.values[0][3])*1e9
    band = sb.values[0][4]

    for i in sb.values[0][6].split('"'):
        if i.find('maxPWVC') >= 0: maxPWV = float(i.split('=')[1].split(' ')[1])

    maxPWV = returnMAXPWVC(maxPWV)


    print sbUID

    science = getSBScience(sbUID)
    field = getSBFields(sbUID)
    target = getSBTargets(sbUID)

    df1 = pd.merge(science,target, left_on='entityPartId',right_on='ObsParameter',how='inner')
    df2 = pd.merge(df1,field,left_on='FieldSource',right_on='entityPartId',how='inner')

    dec = float(df2['latitude'].values[0])
    ToSOT = float(df2['integrationTime'].values[0])
    ###############################

    #############################

    s = sensitivity.SensitivityCalculator(config_dir=configDir)
    result =  s.calcSensitivity(maxPWV,freq,dec=dec,latitude=-23.029, N=antennas4cycle[sbtype][cycle]['NAntOT'], BW=7.5e9, mode='image', N_pol=2,returnFull=True)
    TsysOT = result['Tsys']

    ################################
    ant = getAntennas(asdm.asdmDict['Antenna'])
    NantEB = float(ant.antennaId.count())


    #print "EB:", a.uid
    scan = getScan(asdm.asdmDict['Scan'])
    subscan = getSubScan(asdm.asdmDict['Subscan'])
    scan['target'] = scan.apply(lambda x: True if str(x['scanIntent']).find('OBSERVE_TARGET') > 0 else False ,axis = 1)
    targets = map(unicode.strip,list(scan[scan['target'] == True].sourceName.values))
    subscan['target'] = subscan.apply(lambda x: True if str(x['fieldName']).strip() in targets else False, axis = 1)
    subscan['onsource'] = subscan.apply(lambda x: True if str(x['subscanIntent']) == 'ON_SOURCE'  and x['target'] == True else False, axis = 1)
    subscan['delta'] = subscan.apply(lambda x:  (gtm2(x['endTime']) - gtm2(x['startTime'])) ,axis = 1)

    ToSEB = float(subscan['delta'][subscan['onsource'] == True].sum().total_seconds())

    ################################
    target = list(scan[subscan['target'] == True]['sourceName'].unique())
    scan['atm'] = scan.apply(lambda x: True if str(x['scanIntent']).find('CALIBRATE_ATMOSPHERE') > 0 and x['sourceName'] in target else False,axis =1 )

    ################################
    syscal = getSysCal (asdm.asdmDict['SysCal'])
    syscal['startTime'] = syscal.apply(lambda x: int(x['timeInterval'].split(' ')[1]) - int(x['timeInterval'].split(' ')[2])/2 ,axis=1 )


    ###################################
    spw = getSpectralWindow(asdm.asdmDict['SpectralWindow'])
    spw['repWindow'] = spw.apply(lambda x: findChannel(float(x['chanFreqStart']),float(x['chanFreqStep']), freq, int(x['numChan'])), axis = 1)

    ################################


    df1 = syscal[syscal['startTime'].isin(scan[scan['atm'] == True]['startTime'])]
    spwList = list(spw[spw['repWindow'] != 0]['spectralWindowId'].values)
    spwList = map(unicode.strip, spwList)
    df2 = df1[df1.spectralWindowId.isin(spwList)]
    #df2 = df1[df1['spectralWindowId'] == spw[spw['repWindow'] != 0]['spectralWindowId'].values[0].strip()]
    matches = list(df2.spectralWindowId.unique())
    spw['spectralWindowId'] = spw.apply(lambda x: unicode.strip(x['spectralWindowId']), axis = 1)
    channels = spw[spw.spectralWindowId.isin(matches)].repWindow.values

    data = df2.apply(lambda x: arrayParser(x['tsysSpectrum'],2), axis = 1)
    data2 = pd.DataFrame (data)
    data2.columns = ['hola']
    x = pd.concat([pd.DataFrame(v,index=np.repeat(k,len(v))) for k,v in data2.hola.to_dict().items()])
    x = x.convert_objects(convert_numeric=True)

    TsysEB = list()
    for i in channels:
        TsysEB.append(x[int(i)].median())

    result = {
        'fractionExec' : ((TsysOT/max(TsysEB))**2)*((NantEB*(NantEB-1.))/(antennas4cycle[sbtype][cycle]['NAntOT']*(antennas4cycle[sbtype][cycle]['NAntOT']-1.)))*(ToSEB/ToSOT),
        'TimeOnSource_EB':ToSEB,
                  'TimeOnSource_OT':ToSOT,
                  'Tsys_OT':TsysOT,
                  'Tsys_EB':max(TsysEB),
                  'NAntennas_EB':NantEB,
                  'NAntennas_OT': antennas4cycle[sbtype][cycle]['NAntOT']}

    return result


if __name__ == '__main__':
    config = Configurator()
    config.add_route('fractionExecution', '/')
    config.add_view(fractionExecution, route_name='fractionExecution',renderer='string')
    app = config.make_wsgi_app()
    serve(app, host='0.0.0.0', port=8080)