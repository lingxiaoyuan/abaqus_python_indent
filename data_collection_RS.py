#!/user/bin/python
#-*-coding:UTF-8-*-
from abaqus import *
from abaqusConstants import *
from odbAccess import *
import visualization
import displayGroupOdbToolset as dgo
import math
import os.path
import xlsxwriter
import xlwt
import xlrd

EY=[2200]
ks=[0,0.9]
#postfix='_RS.xlsx'
postfix='RS2.xlsx'
forcebook=xlsxwriter.Workbook('E:\Matlab Files\FORCE'+postfix)
depthbook=xlsxwriter.Workbook('E:\Matlab Files\DEPTH'+postfix)
radiusbook=xlsxwriter.Workbook('E:\Matlab Files\RADIUS'+postfix)
#pdepthbook=xlsxwriter.Workbook('E:\Matlab Files\PDEPTH'+postfix)
pillupbook=xlsxwriter.Workbook('E:\Matlab Files\PILLUP'+postfix)
forcesheet=[]
depthsheet=[]
radiussheet=[]
#pdepthsheet=[]
pillupsheet = []

def datacollection(EYi,ksi):
    strEY=filter(str.isdigit,str(EY[EYi]))
    strks=filter(str.isdigit,str(ks[ksi]))
    if ks[ksi]<0:
        filename='RM'+strks+'C'+'-EY'+strEY
    elif ks[ksi]>0:
        filename='RM'+strks+'T'+'-EY'+strEY
    else:
        filename='RM0'+'-EY'+strEY

    forcesheet[EYi].write(0,ksi,'ks'+str(ks[ksi]))
    radiussheet[EYi].write(0,ksi,'ks'+str(ks[ksi]))
    depthsheet[EYi].write(0,ksi,'ks'+str(ks[ksi]))
    #pdepthsheet[EYi].write(0,ksi,'ks'+str(ks[ksi]))
    pillupsheet[EYi].write(0,ksi,'ks'+str(ks[ksi]))
    #创建数据文件
    if not os.path.exists(filename+'.odb'):
        return
    odb=openOdb(path=filename+'.odb')
    Move=odb.steps['Loading']
    #写入接触总压力
    contactForce=Move.historyRegions['Node INDENTER-1.1'].historyOutputs['RF2'].data
    i=1
    for time,contactForceValue in contactForce :
        forcesheet[EYi].write(i,ksi,-contactForceValue)
        i=i+1
    forcesheet[EYi].write(1,ksi,0)
    '''
    #写入pill-up高度
    i=1
    for frame in Move.frames:
        U=frame.fieldOutputs['U'].values
        maxvalue = 0.0
        for valuepoint in U:
            Ycoord=valuepoint.data[1]
            if Ycoord>maxvalue:
                maxvalue = Ycoord
        pillupsheet[EYi].write(i,ksi,maxvalue)
        i=i+1
        '''
    #写入接触半径
    Ydisdata=[0,]
    i=1
    for frame in Move.frames:
        edgeCoord=0
        edgeValue=None
        cpress=frame.fieldOutputs['CPRESS'].values
        U=frame.fieldOutputs['U']
        for cpressValue in cpress:
            cLabel=cpressValue.nodeLabel
            cpressCoord=cpressValue.instance.nodes[cLabel-1].coordinates[0]
            if frame.frameId!=0 and cpressValue.data!=0 and cpressCoord>edgeCoord: #从第二个frame开始，对所有接触压力不为零的点，找到横坐标最大的点，即为接触点
                edgeValue=cpressValue
                edgeCoord=cpressCoord
        if edgeValue:
            label=edgeValue.nodeLabel
            Xcoord=edgeValue.instance.nodes[label-1].coordinates[0]   #node[507]的nodeLabel为508，故需要-1
            Xdis=U.getSubset(region=edgeValue.instance.nodes[label-1]).values[0].data[0] #只有一个点，第一个点就是Value[0]
            Ydis=U.getSubset(region=edgeValue.instance.nodes[label-1]).values[0].data[1]
            contactRadius=Xcoord+Xdis
            Ydisdata.append(Ydis)
            radiussheet[EYi].write(i,ksi,contactRadius)
        i=i+1
    radiussheet[EYi].write(1,ksi,0)

    #写入压入深度
    penetration=Move.historyRegions['Node INDENTER-1.1'].historyOutputs['U2'].data
    i=1
    for time,penetrationValue in penetration :
        penetrationValue=-(penetrationValue)
        depthsheet[EYi].write(i,ksi,penetrationValue)
        #pdepth=penetrationValue+Ydisdata[i-1]
        #pdepthsheet[EYi].write(i,ksi,pdepth)
        i=i+1

for EYi in range(len(EY)):
    strEY=filter(str.isdigit,str(EY[EYi]))
    forcesheet.append(forcebook.add_worksheet('EY'+strEY))
    depthsheet.append(depthbook.add_worksheet('EY'+strEY))
    radiussheet.append(radiusbook.add_worksheet('EY'+strEY))
    #pdepthsheet.append(pdepthbook.add_worksheet('EY'+strEY))
    pillupsheet.append(pillupbook.add_worksheet('EY'+strEY))
    for ksi in range(len(ks)):
        datacollection(EYi,ksi)
forcebook.close()
depthbook.close()
radiusbook.close()
#pdepthbook.close()
pillupbook.close()


