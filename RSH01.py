#!/user/bin/python
#-*-coding:UTF-8-*-
#!/user/bin/python
#-*-coding:UTF-8-*-
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import math
import displayGroupMdbToolset as dgm
executeOnCaeStartup()
Mdb()

# only have the mesh of thinPart refined

# <editor-fold desc="parameter">
kh=0.1 ; kd = 0.02 ; ks=0.5;
mdb.saveAs(pathName='C:/ABAQUSF/Residual/RSH01')   #InitialConditionMethod
mdb.Job(name='RSH01RS001T', model='Model-1', description='', type=ANALYSIS,
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4,
    numDomains=4, numGPUs=0)
#l=height&weight;E=elastic modulus;miu=poisson ratio;athin=expansion ratio of thinPart;RS=Residual Stress
kEY = 275
l = 40.0  ; R=2; miu=0.3 ;Y=276; E = Y*kEY*(1-miu**2); RS=ks*E
tr = 5
h = kh*R    #thickness of residuelayer
Disp = -kd*R
mes = h/(2**tr)      #minElementSize
datanum = 100  #data of numintervals(output of field and history)
biamag=10    #lagist element size=biamag*mesf

alpha = 0.001
deltaT = 0.5*(1-miu)/alpha/(1-miu**2)/kEY
print 'RS=',RS
print 'mes=',mes
print 'l=',l
print 'R=',R
print 'Disp=',Disp
# </editor-fold>

#create Part
# <editor-fold desc="Part">
#create rigid indenter
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=1000.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -500.0), point2=(0.0, 500.0))
s.FixedConstraint(entity=g[2])
s.ArcByCenterEnds(center=(0.0, R), point1=(R, R), point2=(0.0,
    0.0), direction=CLOCKWISE)
s.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
p = mdb.models['Model-1'].Part(name='indenter', dimensionality=AXISYMMETRIC,
    type=ANALYTIC_RIGID_SURFACE)
p = mdb.models['Model-1'].parts['indenter']
p.AnalyticRigidSurf2DPlanar(sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
p.ReferencePoint(point=v1[0])

#create elastic Body
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=2000.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.sketchOptions.setValues(viewStyle=AXISYM)
s1.setPrimaryObject(option=STANDALONE)
s1.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s1.FixedConstraint(entity=g[2])
s1.rectangle(point1=(0.0, 0.0), point2=(l,-h))
s1.CoincidentConstraint(entity1=v[0], entity2=g[2], addUndoState=False)
p = mdb.models['Model-1'].Part(name='thinPart', dimensionality=AXISYMMETRIC,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['thinPart']
p.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['thinPart']
del mdb.models['Model-1'].sketches['__profile__']
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=2000.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -1000.0), point2=(0.0, 1000.0))
s.FixedConstraint(entity=g[2])
s.rectangle(point1=(0.0, -h), point2=(l, -l))
s.CoincidentConstraint(entity1=v[0], entity2=g[2], addUndoState=False)
p = mdb.models['Model-1'].Part(name='thickPart', dimensionality=AXISYMMETRIC,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['thickPart']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
# </editor-fold>

#creat Propriety
# <editor-fold desc="Propriety">
mdb.models['Model-1'].Material(name='M-thick')
mdb.models['Model-1'].materials['M-thick'].Elastic(table=((E, miu),))
mdb.models['Model-1'].materials['M-thick'].Plastic(table=((Y, 0.0), ))
mdb.models['Model-1'].materials['M-thick'].Expansion(table=((alpha, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='S-thick',
    material='M-thick', thickness=None)
p = mdb.models['Model-1'].parts['thickPart']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p.SectionAssignment(region=region, sectionName='S-thick', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].Material(name='M-thin')
mdb.models['Model-1'].materials['M-thin'].Elastic(table=((E, miu),))
mdb.models['Model-1'].materials['M-thin'].Plastic(table=((Y, 0.0), ))
mdb.models['Model-1'].materials['M-thin'].Expansion(table=((alpha, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='S-thin',
    material='M-thin', thickness=None)
p = mdb.models['Model-1'].parts['thinPart']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p.SectionAssignment(region=region, sectionName='S-thin', offset=0.0,
    offsetType=MIDDLE_SURFACE, offsetField='',
    thicknessAssignment=FROM_SECTION)
# </editor-fold>

#create Assembly
# <editor-fold desc="Assembly">
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0),
    point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
p = mdb.models['Model-1'].parts['indenter']
a.Instance(name='indenter-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['thickPart']
a.Instance(name='thickPart-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['thinPart']
a.Instance(name='thinPart-1', part=p, dependent=ON)
# </editor-fold>

#Set&Surf
# <editor-fold desc="creat set&urf">
#set
a = mdb.models['Model-1'].rootAssembly
edges1 = a.instances['thinPart-1'].edges.findAt(((l/2,0.0,0.0),))
a.Set(edges=edges1, name='contactSet')
r1 = a.instances['indenter-1'].referencePoints
refPoints1=(r1[2], )
a.Set(referencePoints=refPoints1, name='indenterRef')
edges1 = a.instances['thinPart-1'].edges.findAt(((0.0,-h/2,0.0),),((0.0,-2*h,0.0),))
edges2 = a.instances['thickPart-1'].edges.findAt(((l,-h/2,0.0),),((l,-2*h,0.0),))
a.Set(edges=edges1, name='X-fix-left')
edges1 = a.instances['thickPart-1'].edges.findAt(((l/2,-l,0.0),))
a.Set(edges=edges1, name='Y-fix-bottom')
faces1 = a.instances['thinPart-1'].faces.getSequenceFromMask(mask=('[#1 ]', ), )
a.Set(faces=faces1, name='Set-thinPart')
faces1 = a.instances['thickPart-1'].faces.getSequenceFromMask(mask=('[#1 ]', ), )
a.Set(faces=faces1, name='Set-thickPart')
faces1 = a.instances['thinPart-1'].faces.getSequenceFromMask(mask=('[#1 ]', ), )
faces2 = a.instances['thickPart-1'].faces.getSequenceFromMask(mask=('[#1 ]', ), )
a.Set(faces=faces1+faces2, name='Set-whole')
#Surface
side1Edges1 = a.instances['thinPart-1'].edges.findAt(((0.01,0.0,0.0),))
a.Surface(side1Edges=side1Edges1, name='contactSurf')
s1 = a.instances['thinPart-1'].edges
side1Edges1 = s1.findAt(((0.01,-h,0.0),))
a.Surface(side1Edges=side1Edges1, name='tieOFthinPart')
s1 = a.instances['thickPart-1'].edges
side1Edges1 = s1.findAt(((0.01,-h,0.0),))
a.Surface(side1Edges=side1Edges1, name='tieOFthickPart')
# </editor-fold>

#create Step
# <editor-fold desc="Step">
mdb.models['Model-1'].StaticStep(name='Move', previous='Initial', nlgeom=ON,maxNumInc=100000)
mdb.models['Model-1'].StaticStep(name='Expansion', previous='Initial')
#create Output
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValuesInStep(
    stepName='Move', numIntervals=datanum)
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValuesInStep(
    stepName='Move', numIntervals=datanum)
regionDef=mdb.models['Model-1'].rootAssembly.sets['contactSet']
mdb.models['Model-1'].FieldOutputRequest(name='F-Contact',
    createStepName='Move', variables=('CSTRESS', ), numIntervals=datanum,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
regionDef=mdb.models['Model-1'].rootAssembly.sets['indenterRef']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Contact',
    createStepName='Move', variables=('U2', 'RF2', 'CF2'), numIntervals=datanum,
    region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
regionDef=mdb.models['Model-1'].rootAssembly.sets['contactSet']
mdb.models['Model-1'].HistoryOutputRequest(name='H-contactSet',
    createStepName='Move', variables=('CFN1', 'CFN2', 'CAREA'),
    numIntervals=datanum, region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
# </editor-fold>

#create Interaction
# <editor-fold desc="Interaction">
mdb.models['Model-1'].ContactProperty('IntProp-1')
a = mdb.models['Model-1'].rootAssembly
side1Edges1 = a.instances['indenter-1'].edges.getSequenceFromMask(mask=('[#1 ]', ), )
region1=regionToolset.Region(side1Edges=side1Edges1)
region2=a.surfaces['contactSurf']
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1',
    createStepName='Move', master=region1, slave=region2, sliding=FINITE,
    enforcement=NODE_TO_SURFACE, thickness=OFF,
    interactionProperty='IntProp-1', surfaceSmoothing=NONE, adjustMethod=NONE,
    smooth=0.2, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
region1=a.surfaces['tieOFthickPart']
region2=a.surfaces['tieOFthinPart']
mdb.models['Model-1'].Tie(name='Constraint-1', master=region1, slave=region2,
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

# </editor-fold>

#create Load
# <editor-fold desc="Load">
a = mdb.models['Model-1'].rootAssembly
region = a.sets['X-fix-left']
mdb.models['Model-1'].XsymmBC(name='BC-X-fix', createStepName='Initial',
    region=region, localCsys=None)
region = a.sets['Y-fix-bottom']
mdb.models['Model-1'].YsymmBC(name='BC-Y-fix', createStepName='Initial',
    region=region, localCsys=None)
region = a.sets['indenterRef']
mdb.models['Model-1'].DisplacementBC(name='BC-Move', createStepName='Initial',
    region=region, u1=SET, u2=SET, ur3=SET, amplitude=UNSET,
    distributionType=UNIFORM, fieldName='', localCsys=None)
mdb.models['Model-1'].TabularAmplitude(name='Amp-RSnull', timeSpan=STEP,
    smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (0.1, 0.01), (0.2, 0.04), (0.3,
    0.09), (0.4, 0.16), (0.5, 0.25), (0.6, 0.36), (0.7, 0.5), (0.8, 0.64), (
    0.9, 0.81), (1.0, 1.0)))
mdb.models['Model-1'].boundaryConditions['BC-Move'].setValuesInStep(
    stepName='Move', u2=Disp,amplitude='Amp-RSnull')
#Initial conditions
a = mdb.models['Model-1'].rootAssembly
region = a.sets['Set-thinPart']
mdb.models['Model-1'].Temperature(name='Predefined Field-T',
    createStepName='Initial', region=region, distributionType=UNIFORM,
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Move')
mdb.models['Model-1'].predefinedFields['Predefined Field-T'].setValuesInStep(
    stepName='Expansion', magnitudes=(-deltaT, ))
mdb.models['Model-1'].predefinedFields['Predefined Field-T'].resetToPropagated(
    stepName='Move')
# </editor-fold>

#create mesh
# <editor-fold desc="Mesh">
# <editor-fold desc="mesh thin part">
#hy:y coordinate of origin point of the part
def MeshRefineA(lengOld,lengNew,round,hy=0.0,meshPart='thinPart'):
    ratio1=2**round
    ratio2=2**(round-1)
    p = mdb.models['Model-1'].parts[meshPart]
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=1414.21, gridSpacing=mes*ratio2, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    l1=s.Line(point1=(lengNew, -hy), point2=(lengNew, -lengNew-hy))
    l2=s.Line(point1=(lengNew, -lengNew-hy), point2=(0.0, -lengNew-hy))
    l3=s.Line(point1=(lengNew+mes*ratio1, -hy), point2=(lengNew+mes*ratio1, -lengNew-mes*ratio1-hy))
    l4=s.Line(point1=(lengNew+mes*ratio1, -lengNew-mes*ratio1-hy), point2=(0.0, -lengNew-mes*ratio1-hy))
    l5=s.Line(point1=(lengNew, -lengNew-hy), point2=(lengNew+mes*ratio1, -lengNew-hy))
    l6=s.Line(point1=(lengNew, -lengNew-hy), point2=(lengNew, -lengNew-mes*ratio1-hy))

    pattern1=s.linearPattern(geomList=(l5, ), vertexList=(), number1=1, spacing1=141.421,
                             angle1=0.0, number2=int(ceil(lengNew/(ratio1*mes))), spacing2=ratio1*mes, angle2=90.0)
    pattern2=s.linearPattern(geomList=(l6, ), vertexList=(), number1=int(ceil(lengNew/(ratio1*mes))),
                             spacing1=ratio1*mes,angle1=180.0, number2=1, spacing2=ratio1*mes, angle2=90.0)
    l7=s.Line(point1=(lengOld+mes*ratio2,-lengOld-mes*ratio2-hy),point2=(lengNew,-lengOld-mes*ratio2-hy))
    l8=s.Line(point1=(lengOld+mes*ratio2,-lengOld-mes*ratio2-hy),point2=(lengOld+mes*ratio2,-lengNew-hy))
    pickedFaces = f.findAt(((lengNew,-hy,0.0),))
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    # creat seed
    i=2
    y=-hy
    num2re=[]  #the region of whose elements number is 2
    while y>=-lengNew-hy:
        num2re.append(e.findAt(((lengNew+mes,y,0.0),),((lengNew,y-mes,0.0),),((lengNew,y+mes,0.0),)))
        y=-hy-ratio1*mes*i
        i=i+2
    x=0
    i=2
    while x<=lengNew:
        num2re.append(e.findAt(((x,-lengNew-hy-mes,0.0),),((x+mes,-lengNew-hy,0.0),),((x-mes,-lengNew-hy,0.0),)))
        x=ratio1*mes*i
        i=i+2
    for i in num2re:
        p.seedEdgeByNumber(edges=i, number=2, constraint=FINER)
    p.seedEdgeBySize(edges=e.findAt(((lengNew-mes,-hy,0.0),),((0.0,-lengNew-mes,0.0),)),
                     size=mes*ratio2,constraint=FINER)

def FinalRefineA(lengf,roundf,hy=0.0,meshPart='thinPart'):
    ratiof=2**roundf
    p = mdb.models['Model-1'].parts[meshPart]
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=1414.21, gridSpacing=ratiof*mes, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    l7=s.Line(point1=(lengf+mes*ratiof,-lengf-mes*ratiof-hy),point2=(l,-lengf-mes*ratiof-hy))
    l8=s.Line(point1=(lengf+mes*ratiof,-lengf-mes*ratiof-hy),point2=(lengf+mes*ratiof,-l))
    pickedFaces = f.findAt(((l,-lengf-mes*ratiof-hy-mes,0.0),))
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    if meshPart=='thickPart':
        pickedEdges2=e.findAt(((l-mes,-hy,0.0),),((l-mes,-lengf-mes*ratiof-hy,0.0),),((l,-l+mes,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
        pickedEdges1=e.findAt(((0,-l+mes,0.0),),((lengf+mes*ratiof,-l+mes,0.0),),((l-mes,-l,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
    else:
        pickedEdges2=e.findAt(((l-mes,-hy,0.0),),((l-mes,-lengf-mes*ratiof-hy,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
        pickedEdges1=e.findAt(((l-mes,-h,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
        p.seedEdgeBySize(edges=e.findAt(((0,-h+mes,0.0),)),size=mes*ratiof,constraint=FINER)

def meshRefineB(lengNew,round):
    ratio=2**round
    p = mdb.models['Model-1'].parts['thinPart']
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=1414.21, gridSpacing=mes, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    l1=s.Line(point1=(lengNew, 0.0), point2=(lengNew, -h))
    l2=s.Line(point1=(lengNew+mes*ratio, 0.0), point2=(lengNew+mes*ratio, -h))
    l3=s.Line(point1=(lengNew,-h),point2=(lengNew+mes*ratio,-h))
    pattern1=s.linearPattern(geomList=(l3, ), vertexList=(), number1=1, spacing1=141.421,
            angle1=0.0, number2=int(ceil(h/(mes*ratio))), spacing2=mes*ratio, angle2=90.0)
    pickedFaces = f.findAt(((l,0.0,0.0),))
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    i=2
    y=0
    num2re=[]  #the region of whose elements number is 2
    while y>=-h:
        num2re.append(e.findAt(((lengNew+mes,y,0.0),),((lengNew,y-mes,0.0),),((lengNew,y+mes,0.0),)))
        y=-ratio*mes*i
        i=i+2
    for i in num2re:
        p.seedEdgeByNumber(edges=i, number=2, constraint=FINER)
    if round==1:
        p.seedEdgeBySize(edges=e.findAt(((0.0,-h+mes,0.0),)),size=mes,constraint=FINER)
    p.seedEdgeBySize(edges=e.findAt(((lengNew-mes,0.0,0.0),)),size=mes*(ratio/2),constraint=FINER)

def FinalRefineB(roundf,hy=0.0,meshPart='thinPart'):
    ratiof=2**roundf
    p = mdb.models['Model-1'].parts[meshPart]
    f, e, d = p.faces, p.edges, p.datums
    if meshPart=='thickPart':
        pickedEdges2=e.findAt(((l-mes,-hy,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
        pickedEdges1=e.findAt(((l-mes,-l,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
    else:
        pickedEdges2=e.findAt(((l-mes,-hy,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=pickedEdges2, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)
        pickedEdges1=e.findAt(((l-mes,-h,0.0),))
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, minSize=ratiof*mes,
            maxSize=biamag*ratiof*mes, constraint=FINER)

def TranMesh(lengNew,round,meshPart='thinPart'):
    ratio1=2**round
    p = mdb.models['Model-1'].parts[meshPart]
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0.0, 0.0, 0.0))
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=1414.21, gridSpacing=mes*ratio1, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    l8=s.Line(point1=(lengNew+mes*ratio1,-lengNew-mes*ratio1),point2=(lengNew+mes*ratio1,-h))
    pickedFaces = f.findAt(((l,0.0,0.0),))
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    p.seedEdgeBySize(edges=e.findAt(((0,-h+mes,0.0),),((lengNew+mes*ratio1,-h+mes,0.0),)),
                     size=mes*ratio1,constraint=FINER)

leng=[]
leng.append(-mes)
#range(2,11):[1,2,3...10]
start = 1
end = tr-1
for i in range(start,end):
    leng.append(leng[i-1]+2*(2**(i-2)*mes)+62*(2**(i-1)*mes))
    mesf=2**i*mes
    roundf=i
    lengf=leng[i]
    if 0.5*mesf+mes > leng[1]:
        meshRefineB(lengNew=leng[i],round=i)
    else:
        MeshRefineA(lengOld=leng[i-1],lengNew=leng[i],round=i)

    if i== end-1:
        FinalRefineB(roundf=i,hy=0)

# </editor-fold>

#mesh thickPart
# <editor-fold desc="mesh thickPart">

mes = 2**(end-2)*mes
leng=[]
leng.append(-mes)
#range(2,11):[1,2,3...10]
start = 1
end = 4
for i in range(start,end):
    if i == start:
        leng.append(leng[i-1]+2*(2**(i-2)*mes)+62*(2**(i-1)*mes))
    else:
        leng.append(leng[i-1]+2*(2**(i-2)*mes)+30*(2**(i-1)*mes))
    roundf=i
    lengf=leng[i]
    MeshRefineA(lengOld=leng[i-1],lengNew=leng[i],round=i,hy=h,meshPart='thickPart')
    if i == end-1:
        FinalRefineA(lengf=leng[i],roundf=i,hy=h,meshPart='thickPart')


'''leng.append(-2**(end1-1)*mes)
leng.append(62*(2**(end1-1)*mes))
start2 = end1
end2 = tr+2
for i in range(start2,end2):
    leng. append(leng[i]+2*(2**(i-1)*mes)+30*(2**(i)*mes))
    MeshRefineA(lengOld=leng[i-1],lengNew=leng[i],round=i,hy=h,meshPart='thickPart')
MeshRefineA(lengOld=leng[end2-1],lengNew=leng[end2],round=end2,hy=h,meshPart='thickPart')'''


'''for i in range(3,5):
    leng.append(leng[i-1]+2*(2**(i+2)*mes)+14*(2**(i+3)*mes))
    MeshRefineA(lengOld=leng[i-1],lengNew=leng[i],round=i,hy=h,meshPart='thickPart')
    if i==4:
        FinalRefineA(lengf=leng[i],roundf=i,hy=h,meshPart='thickPart')'''


# </editor-fold>



