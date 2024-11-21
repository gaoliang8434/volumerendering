#
#   Various utilities for quickly getting volumes and renders set up
#

import sys
import os
import calendar
import datetime
import time
import math
import random
from bishop import *
#from scheduleArchive import *

#archiveDataLocation = os.getenv("LAUGHING_ARCHIVE_DATA_PATH")
#print "Archive location: " + str(archiveDataLocation)


def formattedFrame( f ):
    fstr = str(f)
    if f < 1000:
        fstr = "0" + fstr
    if f < 100:
        fstr = "0" + fstr
    if f < 10:
        fstr = "0" + fstr
    return fstr

def formattedFrameName( parameters, extension ):
    fname = parameters["-name"] + parameters["-case"] + "-" + parameters["-version"] + "." + formattedFrame(parameters["-turntableframe"]) + "." + extension
    return fname





def formattedTime() :
    T = time.localtime(time.time())
    year  = str(T[0])
    month = str(T[1])
    if T[1] < 10 :
        month = "0" + month
    day   = str(T[2])
    if T[2] < 10 :
        day = "0" + day
    hour = str(T[3])
    if T[3] < 10:
        hour = "0" + hour
    minute = str(T[4])
    if T[4] < 10 :
        minute = "0" + minute
    sec = str(T[5])
    if T[5] < 10 :
        sec = "0" + sec
    formatted = year + "-" + month + "-" + day + " " + hour + ":" + minute + ":" + sec
    return formatted


def MergeParameters( p1, p2 ):
    return dict( list(p1.items()) + list(p2.items()) )


def LogIt( tag ):
    log = "bishop"
    user = str( os.getenv("USER")  )
    logFileName = "/group/dpa/logs/" + log + ".log"
    pid = str(os.getpid()) + " "
    currentTime = str( formattedTime() )
    message = pid + currentTime + " " + str(tag) + " "
    for item in sys.argv:
        message = message + str(item) + " "
    if os.path.exists(logFileName ):
        logFile = open( logFileName, 'a' )
        logFile.write(message)
        logFile.close()
    else:
        print("LOG MESSAGE: " + message)


def DeclareIt( tag ):
    user = str( os.getenv("USER")  )
    pid = str(os.getpid()) + " "
    currentTime = str( formattedTime() )
    message = pid + currentTime + " " + str(tag)
    print("BISHOP: " + message)







CmdLineEntries = []
CmdLineDefault = []
CmdLineActual= []



def CmdLineFindIndex( tag ):
    for i in range(len(sys.argv)):
        if sys.argv[i] == tag:
            return i
    return -1

def CmdLineFind( tag, defaultvalue ):
    CmdLineEntries.append( tag )
    CmdLineDefault.append( defaultvalue )
    i = CmdLineFindIndex(tag)
    if i > 0:
        if i < len(sys.argv)-1:
            return sys.argv[i+1]
    return defaultvalue

def CmdLineFindTuple( tag, nbElements ):
    i = CmdLineFindIndex(tag)
    CmdLineEntries.append( tag )
    CmdLineDefault.append( str(nbElements) + " items" )
    defstring = ""
    result = []
    if i > 0:
        if i < len(sys.argv)- nbElements+1:
            for j in range(0,nbElements):
                result.append(sys.argv[i+1+j])
                defstring = defstring + str(sys.argv[i+1+j]) + " "
    CmdLineActual.append( defstring )
    return result


def CmdLineFindIndices( tag ):
    indices = []
    for i in range(len(sys.argv)):
        if sys.argv[i] == tag:
            indices.append(i)
    return indices

def CmdLineFindArray( tag, defaultvalue ):
    CmdLineEntries.append( tag )
    CmdLineDefault.append( defaultvalue )
    clfList = CmdLineFindIndices(tag)
    clfArray = []
    for i in clfList:
        if i < len(sys.argv)-1:
            clfArray.append( sys.argv[i+1] )
    return clfArray

def CmdLineFindTuples( tag, nbElements ):
    CmdLineEntries.append( tag )
    CmdLineDefault.append( str(nbElements) + " items" )
    idx = CmdLineFindIndices(tag)
    result = []
    for i in idx:
        if i < len(sys.argv) - nbElements+1:
            elem = []
            for j in range(0,nbElements):
                elem.append(sys.argv[i+1+j])
            result.append( elem )
    return result

def CmdLineHelp( tag ):
    i = CmdLineFindIndex(tag)
    if i > 0:
        if len(CmdLineEntries) > 0:
            nbchmax = 0
            for parm in CmdLineEntries:
                nbch = len(parm)
                if nbchmax < nbch:
                    nbchmax = nbch
            nbextra = nbchmax - len("Option") 
            extra = "Option"
            for b in range(0,nbextra):
                extra = extra + " "
            print(extra + " -------> Default value")
            parameters = {}
            for j in range(0,len(CmdLineEntries) ):
                parameters = MergeParameters( parameters, { CmdLineEntries[j]:CmdLineDefault[j] } )
            keys = list(parameters.keys())
            keys.sort()
            for j in range(0,len(keys) ):
                extra = ""
                nbextra = nbchmax - len(keys[j])
                for b in range(0,nbextra):
                    extra = extra + " "
                print(keys[j] + extra + " -------> " + str( parameters[keys[j]] ))
        sys.exit()

def PrintParameters( choices ):
    nbchmax = 0
    for parm, value in choices.items():
        nbch = len(parm)
        if nbchmax < nbch:
            nbchmax = nbch
    nbextra = nbchmax - len("Parameters") 
    extra = "Parameters"
    for b in range(0,nbextra):
        extra = extra + " "
    print("============================================================================")
    print(extra + " -------> Value")
    parameters = list(choices.keys())
    parameters.sort()
    for parm in parameters:
        value = choices[parm]
        nbch = len(parm)
        nbextra = nbchmax - nbch
        extra = ""
        for b in range(0,nbextra):
            extra = extra + " "
        print(str(parm) + extra + " -------> " + str(value))
    print("============================================================================")


def ShowParameters( choices ):
    result = ""
    nbchmax = 0
    for parm, value in choices.items():
        nbch = len(parm)
        if nbchmax < nbch:
            nbchmax = nbch
    nbextra = nbchmax - len("Parameters") 
    extra = "Parameters"
    for b in range(0,nbextra):
        extra = extra + " "
    result = result + "============================================================================\n"
    print(extra + " -------> Value")
    parameters = list(choices.keys())
    parameters.sort()
    for parm in parameters:
        value = choices[parm]
        nbch = len(parm)
        nbextra = nbchmax - nbch
        extra = ""
        for b in range(0,nbextra):
            extra = extra + " "
        result = result + str(parm) + extra + " -------> " + str(value) + "\n"
    result = result + "============================================================================\n"
    return result


def CmdLineFindSet( choices ):
    result = {}
    for parm, value in choices.items():
        result[parm] = value
        if parm.count("---")>0:
            newvalue = CmdLineFindTuple( parm, 3 )
            if len(newvalue) == 3:
                result[parm] = newvalue
        elif parm.count("--")>0:
            newvalue = CmdLineFindTuple( parm, 2 )
            if len(newvalue) == 2:
                result[parm] = newvalue
        elif parm.count("-") >0:
            newvalue = CmdLineFind( parm, value  )
            result[parm] = newvalue
        else:
            newvalue = CmdLineFindIndex( parm )
            if newvalue > 0:
                result[parm] = True
    return result


def CmdLineFindValues( choices ):
    result = {}
    for parm, value in choices.items():
        result[parm] = value
        if isinstance( value, str ):
            newValue = CmdLineFind( parm, value )
            result[parm] = newValue
        else:
            try:
                if len(value) > 1:
                    nbArgs = len(value)
                    newvalue = CmdLineFindTuple( parm, nbArgs )
                    if len(newvalue) == nbArgs:
                        result[parm] = newvalue
            except:
                newValue = CmdLineFind( parm, value )
                result[parm] = newValue
    #put full command line into a parameter
    fullcommandline = sys.argv[0]
    for item in sys.argv[1:]:
        fullcommandline = fullcommandline + " " + str(item)
    result["Full Command Line"] = fullcommandline
    return result




vrVolumesContainer = []

def SaveIt( volume ):
    vrVolumesContainer.append( volume )

def KeepVolume( volume ):
    vrVolumesContainer.append( volume )
    return volume

def ClearVolumes():
    vrVolumesContainer = []


def CameraParameters():
    parms = { "-eye":[0,0,10], "-view":[0,0,0], "-up":[0,1,0], "-fov":60.0, "-aspect":(16.0/9.0), "-far":1.0e6, "-near":0.0   }
    return parms

def ImageParameters():
    parms = { "-size":[1920, 1080], "-name":"image.exr", "-location":"", "-version":"1", "-case":"baseline", "-bright":1.0, "-gamma":1.0, "-backgroundcolor":[0.0,0.0,0.0,0.0] }
    return parms

def RenderDataParameters():
    parms = { "-nbdensitysamples":1, "-maxpathlength":30.0, "-ds":0.1, "-scatter":1.0, "-color":[1,1,1], "-ambient":[0.1,0.1,0.1], "-renderllc":[-15,-15,-15], "-renderurc":[15,15,15],  "ScatteredDensityField":"", "ScatteredColorField":"", "EmissiveDensityField":"", "EmissiveColorField":"" }
    return parms

def RenderLoopParameters():
    parms = { "-nbpatches":100, "-nbthreads":1, "-subsamples":1, "-subsampleseed":0 }
    return parms

def DSMParameters():
    parms = { "-nodsm":0, "-dsmcell":[0.1, 0.1, 0.1], "-dsmlength":[10, 10, 10], "-dsmllc":[-5, -5, -5], "-dsmsamples":1,  "-lightstyle":"whitekeyfill",
          "-lightP":[ [0,100,0 ] ], "-lightCd": [ [ 1, 1, 1 ] ], "-dsmrange":-1.0 }
    return parms

def SimulationParameters():
    parms = { "-frame":1, "-nbframes":100, "-t":0.0, "-dt":1.0/24.0, "-substeps":1 }
    return parms

def TurntableParameters():
    parms = { "-nbturntableframes":120, "-rotationaxis":[0,1,0], "-rotationorigin":[0,0,0], "-turntableframe":1 }
    return parms

def VolumeLightParameters():
    parms = { "-volumelightthreshold":1.0, "-volumelightsamples":4, "-volumelightpicks":100000 }
    return parms


def NoiseParameters():
    parms = {
               "-frequency":1,
               "-translate":[0,0,0],
               "-octaves":1.0,
               "-amplitude":1,
               "-offset":0,
               "-fjump":2,
               "-roughness":0.5,
               "-radius":1.0,
               "-capradius":0.1,
               "-pscale":1.0,
               "-gamma":1.0,
               "-time":0.0,
               "-fftLowCutoff":0.01,
               "-fftHighCutoff":1.0,
               "-fftPower":3.5,
               "-fftNbGridPoints":128,
               "-fftLength":10.0,
               "-lognormalmean":1.0,
               "-gaussianstandarddeviation":1.0,
               "-seed":12345,
               "-tangent":[0,0,1],
               "-normal":[0,1,0],
               "-binormal":[1,0,0],
               "-axis":[0,1,0],
               "-angle":0.0,
               "-P":[0,0,0],
               "-v":[0,0,0],
               "-A":[0,0,0],
               "-age":0.0,
               "-lifeTime":1.0,
               "-shutter":0.5,
               "-frameRate":1.0/24.0,
               "-falloff":1.0,
               "-noiseCd":[1,1,1,0],
               "-nbWisps":0,
               "-wispOctaves":1.0,
               "-wispFreq":1.0,
               "-wispTranslate":[0,0,0],
               "-wispOffset":0.0,
               "-wispFjump":2.13,
               "-wispRoughness":0.5,
               "-wispCorrelation":0.0,
               "-wispRadialGroup":1.0/3.0,
               "-wispDisplacementScale":1.0
            }
    return parms


def GenerateStandardBasename( parameters ):
    fnameparts = os.path.splitext( parameters["-name"] )
    name = str(parameters["-location"]) + str(fnameparts[0])
    name = name + str(parameters["-case"]) 
    name = name + "-" + str(int(parameters["-version"])).zfill(4) 
    return name

def GenerateStandardFilename( parameters ):
    fnameparts = os.path.splitext( parameters["-name"] )
    name = GenerateStandardBasename( parameters )
    #name = str(parameters["-location"]) + str(fnameparts[0])
    #name = name + str(parameters["-case"]).capitalize() 
    #name = name + "-" + str(int(parameters["-version"])).zfill(4) 
    name = name + "." + str(int(parameters["-turntableframe"])).zfill(4)
    name = name + fnameparts[1]
    return name


def AllParameters( parms ):
    allparms = { "versionString":versionString() }
    allparms = MergeParameters( allparms, CameraParameters() )
    allparms = MergeParameters( allparms, ImageParameters() )
    allparms = MergeParameters( allparms, RenderDataParameters() )
    allparms = MergeParameters( allparms, RenderLoopParameters() )
    allparms = MergeParameters( allparms, DSMParameters() )
    allparms = MergeParameters( allparms, TurntableParameters() )
    allparms = MergeParameters( allparms, SimulationParameters() )
    allparms = MergeParameters( allparms, VolumeLightParameters() )
    allparms = MergeParameters( allparms, NoiseParameters() )
    allparms = MergeParameters( allparms, parms )
    return allparms


def Turntable( parameters ):
    if float(parameters["-turntableframe"]) == 1.0:
        return parameters
    axisTuple = parameters["-rotationaxis"]
    originTuple = parameters["-rotationorigin"]
    axis = Vector( float(axisTuple[0]), float(axisTuple[1]), float(axisTuple[2]) )
    axis.unitvector()
    origin = Vector( float(originTuple[0]), float(originTuple[1]), float(originTuple[2]) )
    eyeTuple = parameters["-eye"]
    eye = Vector( float(eyeTuple[0]), float(eyeTuple[1]), float(eyeTuple[2]) )
    viewTuple = parameters["-view"]
    view = Vector( float(viewTuple[0]), float(viewTuple[1]), float(viewTuple[2]) )
    upTuple = parameters["-up"]
    up = Vector( float(upTuple[0]), float(upTuple[1]), float(upTuple[2]) )
    theta = float( parameters["-turntableframe"] ) - 1.0
    theta = theta * 2.0 * 3.14159265/ float( parameters["-nbturntableframes"] )
    ctheta = math.cos(theta)
    stheta = math.sin(theta)
    rotatedView = view * ctheta + axis * (axis*view) * (1.0-ctheta) + (axis^view) * stheta
    rotatedUp = up * ctheta + axis * (axis*up) * (1.0-ctheta) + (axis^up) * stheta
    rotatedEye = origin + (eye-origin) * ctheta + axis * (axis*(eye-origin)) * (1.0-ctheta) + (axis^(eye-origin)) * stheta
    parameters["-eye"] = [ rotatedEye.X(), rotatedEye.Y(), rotatedEye.Z() ]
    parameters["-view"] = [ rotatedView.X(), rotatedView.Y(), rotatedView.Z() ]
    parameters["-up"] = [ rotatedUp.X(), rotatedUp.Y(), rotatedUp.Z() ]
    return parameters


def MakeCamera( parameters ):
    camEyeTuple = parameters[ "-eye"]
    camViewTuple = parameters[ "-view"]
    camUpTuple = parameters[ "-up" ]
    camFov = float(parameters["-fov"])
    camAspectRatio = float( parameters[ "-aspect"] )
    camNear = float( parameters[ "-near" ] )
    camFar  = float( parameters[ "-far" ] )
    camEye = Vector( 0, 0, 10 )
    if len(camEyeTuple) == 3:
        camEye = Vector( float(camEyeTuple[0]), float(camEyeTuple[1]), float(camEyeTuple[2]) )
    camView = Vector( 0, 0, 0 ) - camEye
    if len(camViewTuple) == 3:
        camView = Vector( float(camViewTuple[0]), float(camViewTuple[1]), float(camViewTuple[2]) ) - camEye
    camUp = Vector(0, 1, 0)
    if len(camUpTuple) == 3:
        camUp = Vector( float(camUpTuple[0]), float(camUpTuple[1]), float(camUpTuple[2]) )
    cam = Camera()
    cam.setEyeViewUp( camEye, camView, camUp )
    cam.setFov( camFov )
    cam.setAspectRatio( camAspectRatio )
    cam.setNearPlane( camNear )
    cam.setFarPlane( camFar )
    return cam


def MakeImage( parameters ):
    NXNY = parameters[ "-size" ]
    depth = int(parameters.get('-depth', 4))
    if len(NXNY) < 2:
        NXNY = [ 1920, 1080 ]
    img = Image()
    img.reset( int(NXNY[0]), int(NXNY[1]), depth )
    pixvalue = FloatArray()
    for v in parameters["-backgroundcolor"]:
        pixvalue.push_back( float(v) )
    for j in range(0,int(NXNY[1])):
        for i in range(0,int(NXNY[0])):
            setPixel( img, i, j, pixvalue )
    return img


#def MakeRenderData( density, denColor, ambColor, parameters ):
#    rd = RenderData()
#    rd.nbDensitySamples = int( parameters["-nbdensitysamples"] )
#    rd.maxPathlength = float( parameters["-maxpathlength"] )
#    rd.ds = float( parameters["-ds" ] )
#    rd.scatterCoefficient = float( parameters["-scatter"] )
#    SetDensityField( rd, density )
#    SetColorField( rd, denColor )
#    SetAmbientColorField( rd, ambColor )
#    return rd

def MakeSparseRenderData( grid, density, denColor, ambColor, parameters ):
    rd = RenderData()
    rd.nbDensitySamples = int( parameters["-nbdensitysamples"] )
    rd.maxPathlength = float( parameters["-maxpathlength"] )
    rd.ds = float( parameters["-ds" ] )
    rd.scatterCoefficient = float( parameters["-scatter"] )
    SetDensityField( rd, density )
    SetColorField( rd, denColor )
    SetAmbientColorField( rd, ambColor )
    if float(parameters["-dsmrange"]) > 0.0:
        SetDSMRange( rd, float( parameters["-dsmrange"] ) )
    SetSparseGrid( rd, grid )
    return rd


def MakeRenderData( density, ambientdensity, denColor, ambColor, parameters ):
    rd = RenderData()
    rd.nbDensitySamples = int( parameters["-nbdensitysamples"] )
    rd.maxPathlength = float( parameters["-maxpathlength"] )
    rd.ds = float( parameters["-ds" ] )
    rd.scatterCoefficient = Color(float( parameters["-scatter"]),float( parameters["-scatter"]),float( parameters["-scatter"]),float( parameters["-scatter"]) )
    SetDensityField( rd, density )
    SetAmbientDensityField( rd, ambientdensity )
    SetColorField( rd, denColor )
    SetAmbientColorField( rd, ambColor )
    if float(parameters["-dsmrange"]) > 0.0:
        SetDSMRange( rd, float( parameters["-dsmrange"] ) )
    parameters["ScatteredDensityField"] = str(density)
    parameters["ScatteredColorField"] = str(denColor)
    parameters["EmissiveDensityField"] = str(ambientdensity)
    parameters["EmissiveColorField"] = str(ambColor)
    return rd

def MakeSparseRenderData( grid, density, ambientdensity, denColor, ambColor, parameters ):
    rd = RenderData()
    rd.nbDensitySamples = int( parameters["-nbdensitysamples"] )
    rd.maxPathlength = float( parameters["-maxpathlength"] )
    rd.ds = float( parameters["-ds" ] )
    rd.scatterCoefficient = float( parameters["-scatter"] )
    SetDensityField( rd, density )
    SetAmbientDensityField( rd, ambientdensity )
    SetColorField( rd, denColor )
    SetAmbientColorField( rd, ambColor )
    if float(parameters["-dsmrange"]) > 0.0:
        SetDSMRange( rd, float( parameters["-dsmrange"] ) )
    SetSparseGrid( rd, grid )
    return rd



def RenderLoop( camera, image, renderdata, parameters ):
    nbPixels = image.Width() * image.Height()
    nbPatches = parameters["-nbpatches"]
    count = int( nbPixels/nbPatches )
    nbthreads = int( parameters["-nbthreads"])
    subsampleSeed = int( parameters["-subsampleseed"])
    subsamples = int( parameters["-subsamples"])
    setNbCores( nbthreads )
    metr = ProgressMeter( nbPatches, "Render" )
    for pixSet in range(0,nbPatches):
        random.seed( pixSet + subsampleSeed )
        renderdata.startPosition.clear()
        renderdata.startDirection.clear()
        for p in range(0,count):
            pixelIndex = p + pixSet * count
            pixely = int( pixelIndex/image.Width() )
            pixelx = int( pixelIndex - pixely*image.Width() )
            D = camera.view( float(pixelx)/float(image.Width()), float(pixely)/float(image.Height()) )
            P = camera.eye() + D*camera.nearPlane()
            renderdata.startPosition.append( P )
            renderdata.startDirection.append( D )
            if subsamples > 1:
                for sam in range(0,subsamples-1):
                    x = float(pixelx) + random.random()
                    y = float(pixely) + random.random()
                    D = camera.view( x/float(image.Width()), y/float(image.Height()) )
                    P = camera.eye() + D*camera.nearPlane()
                    renderdata.startPosition.append( P )
                    renderdata.startDirection.append( D )
        output = ColorArray()
        ssRayMarchAccumulation( renderdata, output )
        pixelIndex = -1
        for p in range(0,len(output),subsamples):
            pixelIndex += 1
            pixely = int( pixelIndex/image.Width() )
            pixelx = int( pixelIndex - pixely*image.Width() )
            red = 0.0
            green = 0.0
            blue = 0.0
            alpha = 0.0
            for pp in range(0,subsamples):
                red = red + output( p + pp ).X()
                green = green + output( p + pp ).Y()
                blue = blue + output( p + pp ).Z()
                alpha = alpha + output( p + pp ).W()
            red /= subsamples
            green /= subsamples
            blue /= subsamples
            alpha /= subsamples
            pixvalue = FloatArray()
            pixvalue.push_back( red )
            pixvalue.push_back( green )
            pixvalue.push_back( blue )
            pixvalue.push_back( alpha )
            setPixel( image, pixelx, pixely, pixvalue )
        metr.update()



def RenderVolume( renderdata, parameters, callback=None ):
    camera = MakeCamera( parameters )
    image = MakeImage( parameters )
    nbPixels = image.Width() * image.Height()
    nbPatches = parameters["-nbpatches"]
    count = int( nbPixels/nbPatches )
    nbthreads = int( parameters["-nbthreads"])
    # code for subsampling
    subsampleSeed = int( parameters["-subsampleseed"])
    subsamples = int( parameters["-subsamples"])
    setNbCores( nbthreads )
    metr = callback or ProgressMeter( nbPatches, "Render" )
    for pixSet in range(0,nbPatches):
        random.seed( pixSet + subsampleSeed )
        renderdata.startPosition.clear()
        renderdata.startDirection.clear()
        for p in range(0,count):
            pixelIndex = p + pixSet * count
            pixely = int( pixelIndex/image.Width() )
            pixelx = int( pixelIndex - pixely*image.Width() )
            D = camera.view( float(pixelx)/float(image.Width()), float(pixely)/float(image.Height()) )
            P = camera.eye() + D*camera.nearPlane()
            renderdata.startPosition.append( P )
            renderdata.startDirection.append( D )
            if subsamples > 1:
                for sam in range(0,subsamples-1):
                    x = float(pixelx) + random.random()
                    y = float(pixely) + random.random()
                    D = camera.view( x/float(image.Width()), y/float(image.Height()) )
                    P = camera.eye() + D*camera.nearPlane()
                    renderdata.startPosition.append( P )
                    renderdata.startDirection.append( D )
        output = ColorArray()
        ssRayMarchAccumulation( renderdata, output )
        pixelIndex = pixSet * count - 1
        for p in range(0,len(output),subsamples):
            pixelIndex += 1
            pixely = int( pixelIndex/image.Width() )
            pixelx = int( pixelIndex - pixely*image.Width() )
            red = 0.0
            green = 0.0
            blue = 0.0
            alpha = 0.0
            for pp in range(0,subsamples):
                red = red + output[ p + pp ].X()
                green = green + output[ p + pp ].Y()
                blue = blue + output[ p + pp ].Z()
                alpha = alpha + output[ p + pp ].W()
            red /= subsamples
            green /= subsamples
            blue /= subsamples
            alpha /= subsamples
            currentpixel = image.pixel( pixelx, pixely )
            red   =  red   + currentpixel[0]*(1.0-alpha)
            green =  green + currentpixel[1]*(1.0-alpha)
            blue  =  blue  + currentpixel[2]*(1.0-alpha)
            alpha =  alpha + currentpixel[3]*(1.0-alpha)
            pixvalue = FloatArray()
            pixvalue.push_back( red )
            pixvalue.push_back( green )
            pixvalue.push_back( blue )
            pixvalue.push_back( alpha )
            setPixel( image, pixelx, pixely, pixvalue )
#    nbthreads = int( parameters["-nbthreads"])
#    setNbCores( nbthreads )
#    metr = ProgressMeter( nbPatches, "Render" )
#    for pixSet in range(0,nbPatches):
#        renderdata.startPosition.clear()
#        renderdata.startDirection.clear()
#        for p in range(0,count):
#            pixelIndex = p + pixSet * count
#            pixely = int( pixelIndex/image.Width() )
#            pixelx = int( pixelIndex - pixely*image.Width() )
#            D = camera.view( float(pixelx)/float(image.Width()), float(pixely)/float(image.Height()) )
#            P = camera.eye() + D*camera.nearPlane()
#            renderdata.startPosition.append( P )
#            renderdata.startDirection.append( D )
#        output = ColorArray()
#        ssRayMarchAccumulation( renderdata, output )
#        for p in range(0,len(output)):
#            pixelIndex = p + pixSet * count
#            pixely = int( pixelIndex/image.Width() )
#            pixelx = int( pixelIndex - pixely*image.Width() )
#            pixvalue = FloatArray()
#            pixvalue.push_back( output[p].X() )
#            pixvalue.push_back( output[p].Y() )
#            pixvalue.push_back( output[p].Z() )
#            pixvalue.push_back( output[p].W() )
#            setPixel( image, pixelx, pixely, pixvalue )
        metr.update()
    return image



def ComputeLightFrustum( eye, Origin, Length ):
    verts = [
           Vector(  float(Origin[0]), float(Origin[1]), float(Origin[2]) ),
               Vector(  float(Origin[0]) + float(Length[0]), float(Origin[1]), float(Origin[2]) ),
               Vector(  float(Origin[0]), float(Origin[1]) + float(Length[1]), float(Origin[2]) ),
               Vector(  float(Origin[0]), float(Origin[1]), float(Origin[2]) + float(Length[2]) ),
               Vector(  float(Origin[0]) + float(Length[0]), float(Origin[1]) + float(Length[1]), float(Origin[2]) ),
               Vector(  float(Origin[0]), float(Origin[1]) + float(Length[1]), float(Origin[2]) + float(Length[2]) ),
               Vector(  float(Origin[0]) + float(Length[0]), float(Origin[1]), float(Origin[2]) + float(Length[2]) ),
               Vector(  float(Origin[0]) + float(Length[0]), float(Origin[1]) + float(Length[1]), float(Origin[2]) + float(Length[2]) )
            ]
    center =  Vector(  float(Origin[0]) + 0.5*float(Length[0]), float(Origin[1]) + 0.5*float(Length[1]), float(Origin[2]) + 0.5*float(Length[2]) )
    origin =  Vector(  float(Origin[0]), float(Origin[1]), float(Origin[2]) )
    view = (center-eye).unitvector()
    up = (origin - eye).unitvector()
    up = (up - view*dot_product(up,view)).unitvector()
    fov = 0.0
    near = 100000000.0
    far = 0.0
    for v in verts:
        d = (v-eye)
        dmag = d.magnitude()
        if dmag < near:
            near = dmag
        if dmag > far:
            far = dmag
        cos = dot_product(d,view)/dmag
        fv = math.acos(cos) * 180.0/3.14159265
        if fv > fov:
            fov = fv
    fov = 2.0*fov
    if fov >= 180.0:
        fov = 179.0
    aspect = 1.0
    camera = Camera()
    camera.setEyeViewUp( eye, view, up )
    camera.setFov( fov )
    camera.setAspectRatio( aspect )
    camera.setNearPlane( near )
    camera.setFarPlane( far )
    return camera




def MakeDSMs( renderdata, parameters ):
    dsmGridList = []
    dsmCell = parameters[ "-dsmcell" ]
    if len(dsmCell) == 0:
        return dsmGridList    
    dsmLength = parameters[ "-dsmlength"]
    if len(dsmLength) == 0:
        return dsmGridList    
    dsmOrigin = parameters[ "-dsmllc" ]
    if len(dsmOrigin) == 0:
        return dsmGridList    
    lightPs = []
    lightCds = []
    lightstyle = parameters["-lightstyle"]
    if lightstyle == "keyrimfill" :
        lightCds.append( Vector( 0,0,1 ) )
        lightCds.append( Vector( 1,0,0 ) )
        lightCds.append( Vector( 0,1,0 ) )
        lightPs.append( Vector( 0, 950, 200 ) )
        lightPs.append( Vector( 0, -950, 0 ) )
        lightPs.append( Vector( 0, 0, -950 ) )
        #lightPs.append( Vector( 5, 30, 5 ) )
        #lightPs.append( Vector( 0, -30, 0 ) )
        #lightPs.append( Vector( 0, 0, -30 ) )
    elif lightstyle == "keyfill":
        lightCds.append( Vector( 0,0,1 ) )
        lightCds.append( Vector( 0.3,0,0 ) )
        lightPs.append( Vector( 0, 950, 200 ) )
        lightPs.    append( Vector( 0, -950, 0 ) )
    elif lightstyle == "whitekeyrimfill":
        lightCds.append( Vector( 1,1,1 ) )
        lightCds.append( Vector( 0.3,0.3,0.3 ) )
        lightCds.append( Vector( 1,1,1 ) )
        lightPs.append( Vector( 0, 950, 200 ) )
        lightPs.append( Vector( 0, -950, 0 ) )
        lightPs.append( Vector( 0, 0, -950 ) )
    elif lightstyle == "whitekeyfill":
        lightCds.append( Vector( 1,1,1 ) )
        lightCds.append( Vector( 0.3,0.3,0.3 ) )
        lightPs.append( Vector( 0, 950, 200 ) )
        lightPs.append( Vector( 0, -950, 0 ) )
    elif lightstyle == "custom":
        cmdlinelightpositions = parameters["-lightP"]
        cmdlinelightcolors = parameters["-lightCd"]
        nblights = len(cmdlinelightpositions)
        if nblights > len(cmdlinelightcolors):
            nblights = len(cmdlinelightcolors)
	print "\nUsing custom lights"
        print "Light   Position            Color"
        for i in range(0,nblights):
            print(str(i) + "       " + str( cmdlinelightpositions[i] ) + "           " + str( cmdlinelightcolors[i] ))
            lightPs.append( Vector( float( cmdlinelightpositions[i][0]), float(cmdlinelightpositions[i][1]), float(cmdlinelightpositions[i][2]) ) )
            lightCds.append( Vector( float( cmdlinelightcolors[i][0]), float(cmdlinelightcolors[i][1]), float(cmdlinelightcolors[i][2]) ) )
    nblights = len(lightPs)
    if nblights > len(lightCds):
        nblights = len(lightCds)
    for i in range(0,nblights):
    	if parameters["-nodsm"] == 0:
        	X0 = Vector( float(dsmOrigin[0]),  float(dsmOrigin[1]),  float(dsmOrigin[2]) )
        	XU = Vector( float(dsmOrigin[0]) + float(dsmLength[0]),  float(dsmOrigin[1]) + float(dsmLength[1]),  float(dsmOrigin[2]) + float(dsmLength[2]) )
        	Cell = Vector(  float(dsmCell[0]),  float(dsmCell[1]),  float(dsmCell[2]) )
        	lightnx = int( float(dsmLength[0])/float(dsmCell[0]) )
        	lightny = int( float(dsmLength[1])/float(dsmCell[1]) )
        	lightnz = int( float(dsmLength[2])/float(dsmCell[2]) )
		dsm_aabb = AABB( X0, XU )
		if isInside( dsm_aabb, lightPs[i] ):
			lightcam = MakeCamera( parameters )
			print("\nLight " + str(i) + " is inside dsm box.\nDSM Camera:" + str(lightcam))
        		fgb = makeGridBox( X0, XU, Cell )
        		dsmGrid = makeGrid( fgb, 0.0 )
        		#fgb = makeFrustumBox( lightnx, lightny, lightnz, lightcam )
        		#dsmFGrid = makeFrustumGrid( fgb, 0.0 )
        		RayMarchVisibleDSMAccumulation( renderdata, lightPs[i], lightcam, dsmGrid )
        		#nbsamples = int( parameters["-dsmsamples"])
        		#RayMarchDSMAccumulation( renderdata.densityField, dsmFGrid, nbsamples )
        		AddDSM( renderdata, gridded(dsmGrid) )
        		dsmGridList.append( dsmGrid )
		else:
        		lightcam = ComputeLightFrustum( lightPs[i], dsmOrigin, dsmLength );
        		print("Light " + str(i) + " is outside dsm box.\nDSM Camera: " + str(lightcam))
        		fgb = makeFrustumBox( lightnx, lightny, lightnz, lightcam )
        		dsmFGrid = makeFrustumGrid( fgb, 0.0 )
        		nbsamples = int( parameters["-dsmsamples"])
        		RayMarchDSMAccumulation( renderdata.densityField, dsmFGrid, nbsamples )
        		AddDSM( renderdata, gridded(dsmFGrid) )
        		dsmGridList.append( dsmFGrid )
        renderdata.lightColor.append( Color( lightCds[i].X(), lightCds[i].Y(), lightCds[i].Z(), 1 ) )
        renderdata.lightPosition.append( lightPs[i] )
    return dsmGridList


def PackageDSMs( dsms, renderdata ):
    griddedVolumeList = []
    for dsm in dsms:
        griddedVolumeList.append( gridded(dsm) )
    return griddedVolumeList

def WriteImage( img, parameters ):
    fname = GenerateStandardFilename(parameters)
    keys = StringArray()
    values = StringArray()
    for key,value in parameters.items():
        keys.push_back(str(key))
        values.push_back(str(value))
    writeOIIOImage( fname, img, keys, values, 1.0, 1.0 )
    print("Image written to file " + str(fname))

def OutputImage( img, parameters ):
    WriteImage( img, parameters)
    #fname = GenerateStandardFilename(parameters)
    #scheduleArchiveAction( [fname], archiveDataLocation )

def ReadImage( parameters ):
    fname = GenerateStandardFilename(parameters)
    img = MakeImage(parameters)
    readOIIOImage( fname, img )
    print("Image read from file " + str(fname))
    return img



def RenderSetup( volum, parms ):
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    return data


def RenderFromData( data, parms ):
    print("Executing Standard Render")
    renderStartTime = formattedTime()
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )


def StandardRender( volum, parms ):
    print("Executing Standard Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )

def StandardIntervalRender( volum, parms, intervalset ):
    print("Executing Standard Interval Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    #AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    SetIntervalTree( data, intervalset )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )

def StandardHGIntervalRender( volum, parms, intervalset, hg ):
    print("Executing Standard Henyey-Greenstein Interval Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    SetHenyeyGreensteinPhaseFunction( data, hg )
    SetIntervalTree( data, intervalset )
    #llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    #urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    #AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    #packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )

def VolumeHGIntervalRender( volum, parms, intervalset, hg, light_density, light_color, lightbnds ):
    print("Executing Volume Light Henyey-Greenstein Interval Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    SetHenyeyGreensteinPhaseFunction( data, hg )
    SetIntervalTree( data, intervalset )
    vl = makeVolumeLight( lightbnds, light_density, light_color, float(parms["-volumelightthreshold"]) )
    SetSampleLimit( vl, int( parms["-volumelightpicks"] ) )
    AddVolumeLight( data, vl )
    SetVolumeLightSamples( data, int( parms["-volumelightsamples"] ) )
    #llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    #urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    #AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    #packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )



def StandardColorRender( volum, colorfield, parms ):
    print("Executing Standard Color Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  colorfield
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )



def StandardSparseRender( grid, volum, parms ):
    print("Executing Standard Sparse Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeSparseRenderData( grid, volum, volum, litcolor, ambient, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )


def StandardTwoPartRender( emitted, scattered, emittedCd, scatteredCd, parms ):
    print("Executing Standard Two Part Render")
    renderStartTime = formattedTime()
    data = MakeRenderData( scattered, emitted, scatteredCd, emittedCd, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )

def StandardTwoPartIntervalRender( emitted, scattered, emittedCd, scatteredCd, parms, intervalset ):
    print("Executing Standard Two Part Interval Render")
    renderStartTime = formattedTime()
    data = MakeRenderData( scattered, emitted, scatteredCd, emittedCd, parms )
    PrintParameters(parms)
    SetIntervalTree( data, intervalset )
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    #packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )

def StandardTwoPartHGIntervalRender( emitted, scattered, emittedCd, scatteredCd, parms, intervalset, hg ):
    print("Executing Standard Two Part Henyey-Greenstein Interval Render")
    renderStartTime = formattedTime()
    data = MakeRenderData( scattered, emitted, scatteredCd, emittedCd, parms )
    PrintParameters(parms)
    SetIntervalTree( data, intervalset )
    SetHenyeyGreensteinPhaseFunction( data, hg )
    #llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    #urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    #AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    #packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )



def StandardTwoPartSparseRender( grid, emitted, scattered, emittedCd, scatteredCd, parms ):
    print("Executing Standard Two Part Sparse Render")
    renderStartTime = formattedTime()
    data = MakeSparseRenderData( grid, scattered, emitted, scatteredCd, emittedCd, parms )
    PrintParameters(parms)
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    AddBoundingBox( data, llc, urc )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    image = RenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )




def TurntableRender( volum, parms ):
    print("Executing Turntable Render")
    renderStartTime = formattedTime()
    ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    data = MakeRenderData( volum, volum, litcolor, ambient, parms )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    eye0 = parms["-eye"]
    nbframes = parms["-nbturntableframes"]
    for i in range(0,nbframes):
        renderStartTime = formattedTime()
        print("\n\n\tFrame " + str(i+1) + "\n\n")
        parms["-turntableframe"] = i+1
        parms["-eye"] = eye0
        parms = Turntable( parms )
        image = RenderVolume( data, parms )
        renderEndTime = formattedTime()
        nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
        OutputImage( image, nparms )




def TurntableTwoPartRender( emitted, scattered, emittedCd, scatteredCd, parms ):
    print("Executing Turntable Render")
    data = MakeRenderData( scattered, emitted, scatteredCd, emittedCd, parms )
    dsms = MakeDSMs( data, parms )
    packagedDsms = PackageDSMs( dsms, data )
    eye0 = parms["-eye"]
    nbframes = parms["-nbturntableframes"]
    for i in range(0,nbframes):
        renderStartTime = formattedTime()
        print("\n\n\tFrame " + str(i+1) + "\n\n")
        parms["-turntableframe"] = i+1
        parms["-eye"] = eye0
        parms = Turntable( parms )
        image = RenderVolume( data, parms )
        renderEndTime = formattedTime()
        nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
        OutputImage( image, nparms )










totalexecutiontimemeter = ProgressMeter(1, "TOTAL EXECUTION TIME")

def beginJob():
    print(versionString())
    print(LabLogo())
    LogIt( "OPEN" )

def endJob():
    totalexecutiontimemeter.update()
    LogIt("CLOSE")



class BackPlaneRenderer:
    def __init__(self, bgcolor, aabb):
        self.bgcolor = bgcolor
        self.aabb = aabb
    
    def intersection( self, P, D ):
        distance = self.aabb.farIntersection( P, D )
        if distance<=0:
            return Color(0,0,0,0)
        hitP = P + D*distance
        normal = self.aabb.normal( hitP )
        pseudoD = ( hitP - (aabb.llc()+aabb.urc())*0.5 )
        pseudoD.normalize()
        colorMod = math.fabs( normal * pseudoD )
        return self.bgcolor*colorMod

    def shadowed_intersection( self, P, D, sf, ds, scatter, lightPs, lightCds ):
        distance = self.aabb.farIntersection( P, D )
        Ci = Color(0,0,0,0)
        if distance <= 0:
            return Ci
        startP = P + D*distance
        normal = self.aabb.normal( P + D*distance )
        for i in range(0, len(lightPs) ):
            direction = lightPs[i] - startP
            direction.normalize()
            t = 1.0
            endDistance = self.aabb.farIntersection( startP + direction*0.00001, direction )
            if endDistance > 0:
                endP = startP + direction*endDistance
                t = transmissivity( sf, startP, endP, ds, scatter )
            Ci = Ci + lightCds[i] * ( t * math.fabs( normal*direction ) )
        return Ci * self.bgcolor



