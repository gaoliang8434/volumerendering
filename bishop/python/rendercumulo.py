#!/usr/bin/python


from bishoputils import *
from ImplicitShapes import *


def CumuloParameters():
	return {  "-nbGenerations":2, "-fjump":3.2, "-rough":0.25, "-noisescale":0.5, "-noisegamma":0.3  }




beginJob()

cumuloParameters = CumuloParameters()
parms = AllParameters( MergeParameters(cumuloParameters, ImplicitJackParameters()) )
CmdLineHelp("-h")
parms["-scatter"] = 30.0
parms["-dsmllc"] = [ -2, -2, -2 ]
parms["-dsmlength"] = [ 4,4,4 ]
parms["-dsmNxNyNz"] = [ 400, 400, 400 ]
parms["-near"] = 6.0
parms["-maxpathlength"] = 6.0
parms["-noisescale"] = 0.5
parms["-fov"] = 35.0
parms["-ambient"] = [ 0, 0, 0 ]
parms["-ds"] = 0.007
parms["-name"] = "cumulo.exr"
parms["-color"] = [10.0, 10.0, 10.0]
parms = CmdLineFindValues( parms )
PrintParameters( parms )



def CumuloNoiseLayers( parameters, noise ):
	noiseArray = FloatVolumeArray()
	freqScale = 1.0
	ampScale = 1.0
	fjump = float(parameters["-fjump"])
	rough = float(parameters["-rough"])
	nbGenerations = int(parameters["-nbGenerations"])
        for i in range(0,nbGenerations):
		nextNoise = scale( noise, freqScale )
		#nextNoise = KeepVolume(ScaleVolume( noise, freqScale ))
		nextnextNoise = nextNoise*constant(ampScale) 
		#nextnextNoise = KeepVolume( MultiplyVolume( nextNoise, ampScale ) )
		noiseArray.push_back( nextnextNoise )
		ampScale = ampScale * fjump
		freqScale = freqScale * rough
	return noiseArray



objSDF = Sphere( Vector(0,0,0), 1.0 )

noiseparm = Noise_t()

noisegamma = float( parms["-noisegamma"] )


basenoise = perlin()
basenoise.getParameters( noiseparm )

noiseparm.wavelength = 2.245356
noiseparm.octaves = 3.0
noiseparm.roughness = 0.7
noiseparm.amplitude = 2.0


basenoise.setParameters(noiseparm)
basenoisevolume = noise( basenoise )
noiseScale = float( parms["-noisescale"] )
cumulonoise0 = basenoisevolume*constant( noiseparm.amplitude )
cumulonoise1 = gamma( cumulonoise0, noisegamma )
cumulonoise = cumulonoise1*constant(noiseScale )

Xf = ImplicitPointVectorVolume( objSDF, 0.01, 6, 1.0 )
identity = IdentityVectorVolume()
dXf = SubtractVectorVolume( Xf, identity )
#Grid it
Xgrid = SparseVectorGrid()
#Xgrid = VectorVolumeGrid()
Xgrid.init( 100, 100, 100, 12, 12, 12, Vector(-6,-6,-6) )
Sample( Xgrid, dXf )
#VectorSample( Xgrid, dXf )
Xgg = GriddedVectorVolume( Xgrid )
X = AddVectorVolume( identity, Xgg )
mappedcumulonoise = WarpVolume( cumulonoise, X )
cloud = AddVolume( objSDF , mappedcumulonoise ) 

theshold = 1.0;
mcloud = MultiplyVolume( cloud, theshold )
density = ClampVolume( mcloud, 0.0, 1.0 )
parms["-name"] = "cumulo_one_generation.exr"
StandardRender( density, parms )

### 2nd Generation

basenoise.getParameters(noiseparm)
noiseparm.wavelength = noiseparm.wavelength*2.2
noiseparm.octaves = 2.0
noiseparm.roughness = 0.5
noiseScale = noiseScale * 0.5
nextnoisegamma = 1.0



nextbasenoise = PerlinFractalSum()
nextbasenoise.setParameters(noiseparm)
nextbasenoisevolume = NoiseVolume( nextbasenoise )
nextcumulonoise0 = GammaVolume( nextbasenoisevolume, nextnoisegamma )
nextcumulonoise = MultiplyVolume( nextcumulonoise0, noiseScale )

nextXf = ImplicitPointVectorVolume( cloud, 0.01, 6, 0.5 )
nextdXf = SubtractVectorVolume( nextXf, identity )
nextXgrid = SparseVectorGrid()
#nextXgrid = VectorVolumeGrid()
nextXgrid.init( 200, 200, 200, 12, 12, 12, Vector(-6,-6,-6) )
Sample( nextXgrid, nextdXf )
#VectorSample( nextXgrid, nextdXf )
nextXgg = GriddedVectorVolume( nextXgrid )
nextX = AddVectorVolume( identity, nextXgg )


nextmappedcumulonoise = WarpVolume( nextcumulonoise, nextX )
nextcloud = AddVolume( cloud, nextmappedcumulonoise ) 


mcloud = MultiplyVolume( nextcloud, theshold )
density = ClampVolume( mcloud, 0.0, 1.0 )
parms["-name"] = "cumulo_two_generations_notcleared.exr"
StandardRender( density, parms )



clearanceThickness = 0.5
clearanceGamma = 2.0
clearanceField0 = MultiplyVolume( mappedcumulonoise, 1.0/clearanceThickness )
clearanceField1 = ClampVolume( clearanceField0, 0.0001, 1.0 )
clearanceField2 = GammaVolume( clearanceField1, clearanceGamma )
nextclearedcumulonoise = MultiplyVolume( nextmappedcumulonoise, clearanceField2 )

nextnextcloud = AddVolume( cloud, nextclearedcumulonoise ) 


mcloud = MultiplyVolume( nextnextcloud, theshold )
density = ClampVolume( mcloud, 0.0, 1.0 )
parms["-name"] = "cumulo_two_generations_cleared.exr"
StandardRender( density, parms )


endJob()
