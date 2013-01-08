import os

for simFlags in simFlagsList[:1]:
    simFlags['errorQuantities']=['u']
    simFlags['errorTypes']= ['numericalSolution'] #compute error in soln and glob. mass bal
    simFlags['errorNorms']= ['L2'] #compute L2 norm in space
    simFlags['errorTimes']= ['All']#['Last']
    simFlags['echo']=True
    simFlags['dataFile']       = simFlags['simulationName']+'_results'
    simFlags['dataDir']        = os.getcwd()+'/results'
    simFlags['storeQuantities']= ['simulationData','errorData'] #include errorData for mass bal
    simFlags['storeTimes']     = ['All']#['Last']
    #simFlags['plotQuantities'].append('u_exact')
    #simFlags['plotQuantities'].append('velocity')
start
quit
