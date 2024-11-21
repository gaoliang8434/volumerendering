#
#  Various utilities for quickly getting volumes and renders set up
#  CUDA Version
#

import sys
import os
import calendar
import datetime
import time
import math
import random

from bishop import *
from bishopUtils import *

def cuMakeRenderData( density, ambientdensity, denColor, ambColor, parameters ):
    rd = cu_RenderData()
    rd.nbDensitySamples = int( parameters["-nbdensitysamples"] )
    rd.maxPathlength = float( parameters["-maxpathlength"] )
    rd.ds = float( parameters["-ds" ] )
    rd.scatterCoefficient = float( parameters["-scatter"] )
    rd.clampv = float( parameters["-clamp"] )
    cu_SetDensityField( rd, density )
    cu_SetAmbientDensityField( rd, ambientdensity )
    cu_SetColorField( rd, denColor )
    cu_SetAmbientColorField( rd, ambColor )
    parameters["ScatteredDensityField"] = str(density)
    parameters["ScatteredColorField"] = str(denColor)
    parameters["EmissiveDensityField"] = str(ambientdensity)
    parameters["EmissiveColorField"] = str(ambColor)
    return rd


def cuRenderVolume( renderdata, parameters, callback=None ):
    camera = MakeCamera( parameters )
    image = MakeImage( parameters )
    nbPixels = image.Width() * image.Height()
    nbPatches = parameters["-nbpatches"]
    count = int( nbPixels/nbPatches )
    #print("Pixels in patch: " + str(count))
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
        cu_ssRayMarchAccumulation( renderdata, output )
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
        metr.update()
    return image


def cuMakeDSMs( renderdata, parameters ):
    dsmGridList = []
    if parameters["-nodsm"] > 0:
        return dsmGridList
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
        lightPs.append( Vector( 5, 30, 5 ) )
        lightPs.append( Vector( 0, -30, 0 ) )
        lightPs.append( Vector( 0, 0, -30 ) )
    elif lightstyle == "keyfill":
        lightCds.append( Vector( 0,0,1 ) )
        lightCds.append( Vector( 0.3,0,0 ) )
        lightPs.append( Vector( 0, 950, 200 ) )
        lightPs.append( Vector( 0, -950, 0 ) )
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
        print("\nLight   Position            Color")
        for i in range(0,nblights):
            print(str(i) + "       " + str( cmdlinelightpositions[i] ) + "           " + str( cmdlinelightcolors[i] ))
            lightPs.append( Vector( float( cmdlinelightpositions[i][0]), float(cmdlinelightpositions[i][1]), float(cmdlinelightpositions[i][2]) ) )
            lightCds.append( Vector( float( cmdlinelightcolors[i][0]), float(cmdlinelightcolors[i][1]), float(cmdlinelightcolors[i][2]) ) )
    nblights = len(lightPs)
    if nblights > len(lightCds):
        nblights = len(lightCds)
    for i in range(0,nblights):
        X0 = Vector( float(dsmOrigin[0]),  float(dsmOrigin[1]),  float(dsmOrigin[2]) )
        XU = Vector( float(dsmOrigin[0]) + float(dsmLength[0]),  float(dsmOrigin[1]) + float(dsmLength[1]),  float(dsmOrigin[2]) + float(dsmLength[2]) )
        Cell = Vector( float(dsmCell[0]),  float(dsmCell[1]),  float(dsmCell[2]) )
        lightcam = ComputeLightFrustum( lightPs[i], dsmOrigin, dsmLength );
        print("DSM Camera: " + str(lightcam))
        lightnx = int( float(dsmLength[0])/float(dsmCell[0]) )
        lightny = int( float(dsmLength[1])/float(dsmCell[1]) )
        lightnz = int( float(dsmLength[2])/float(dsmCell[2]) )
        fgb = makeFrustumBox( lightnx, lightny, lightnz, lightcam )
        dsmFGrid = makeFrustumGrid( fgb, 0.0, int(parameters["-dsmpartition"]) )
        nbsamples = int( parameters["-dsmsamples"])
        volum = gridded(renderdata.densityField)
        volum = clamp( volum / constant(renderdata.clampv), 0.0, 1.0 )
        RayMarchDSMAccumulation( volum, dsmFGrid, nbsamples )
        initCUDA(dsmFGrid)
        cu_AddDSM( renderdata, dsmFGrid )
        dsmGridList.append( dsmFGrid )
        #gb = makeGridBox( X0, XU, Cell )
        #dsmGrid = makeGrid( gb, 0.0 )
        #RayMarchDSMAccumulation( renderdata.densityField, lightPs[i], renderdata.ds, dsmGrid )
        #AddDSM( renderdata, gridded(dsmGrid) )
        #dsmGridList.append( dsmGrid )
        renderdata.lightColor.append( Color( lightCds[i].X(), lightCds[i].Y(), lightCds[i].Z(), 1 ) )
        renderdata.lightPosition.append( lightPs[i] )
    return dsmGridList


def cuStandardRender( grid, parms ):
    #print("Executing CUDA Standard Render")
    renderStartTime = formattedTime()

    # Need to convert everything to grids for CUDA render
    llc = Vector( float(parms["-dsmllc"][0]), float(parms["-dsmllc"][1]), float(parms["-dsmllc"][2]) )
    NND = Vector( float(parms["-dsmlength"][0]), float(parms["-dsmlength"][1]), float(parms["-dsmlength"][2]) )
    res = Vector( float(parms["-dsmcell"][0]), float(parms["-dsmcell"][1]), float(parms["-dsmcell"][2]) )
    gb = makeGridBox(llc, llc+NND, res)

    # ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
    amb_grid = makeGrid(gb, Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), int(parms["-dsmpartition"])))
    initCUDA(amb_grid)

    # litcolor =  constant( Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), 0.0 ) )
    litC_grid = makeGrid(gb, Color( float(parms["-color"][0]), float(parms["-color"][1]), float(parms["-color"][2]), int(parms["-dsmpartition"])))
    initCUDA(litC_grid)

    data = cuMakeRenderData( grid, grid, litC_grid, amb_grid, parms )
    llc = Vector( float(parms["-renderllc"][0]), float(parms["-renderllc"][1]), float(parms["-renderllc"][2]) )
    urc = Vector( float(parms["-renderurc"][0]), float(parms["-renderurc"][1]), float(parms["-renderurc"][2]) )
    cu_AddBoundingBox( data, llc, urc )
    dsms = cuMakeDSMs( data, parms )
    #packagedDsms = PackageDSMs( dsms, data )
    # cu_InitFields( data )

    image = cuRenderVolume( data, parms )
    renderEndTime = formattedTime()
    nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
    OutputImage( image, nparms )
