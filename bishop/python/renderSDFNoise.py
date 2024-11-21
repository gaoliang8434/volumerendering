#!/usr/bin/python


from bishoputils import *
from ImplicitShapes import *
from Cumulo import *


def NoiseParameters():
	return { "-freq":0.645356, "-octaves":3.0, "-roughness":0.5, "-noisescale":0.4  }



beginJob()




parms = AllParameters( MergeParameters( NoiseParameters(), MergeParameters( ImplicitJackParameters(),  TurntableParameters() )  )  )
CmdLineHelp("-h")
parms["-ds"] = 0.01
parms["-dsmNxNyNz"] = [ 150, 150, 150 ]
parms["-ambient"] = [ 0, 0, 0 ]
parms["-size"] = [ 960, 540 ]
parms["-version"] = 4
parms["-scatter"] = 5
parms["-case"] = "LocalNormal"
parms = CmdLineFindValues( parms )
PrintParameters( parms )




noiseparm = Noise_t()

noiseparm.wavelength = float( parms["-freq"] )
noiseparm.octaves = float( parms["-octaves"] )
noiseparm.roughness = float( parms["-roughness"] )
noiseScale = float( parms["-noisescale"] )


basenoise = perlin(noiseparm)
basenoisevolume = SFNoise( basenoise )
absbasenoisevolume = abs( basenoisevolume )
cumulonoise =  absbasenoisevolume * constant(noiseScale)


sdf = ImplicitJackSDF( parms )
sdfnormal =  grad(sdf) / pow( grad(sdf)*grad(sdf) + constant(0.0000001), 0.5 )
#cloud = sdf + warp(cumulonoise, sdfnormal )
cloud = sdf + warp(cumulonoise,  ImplicitSurfacePoint( sdf, 0.1, 3 ) )
density = clamp( cloud/constant(0.1), 0.0, 1.0 )


StandardRender( density, parms )

endJob()
