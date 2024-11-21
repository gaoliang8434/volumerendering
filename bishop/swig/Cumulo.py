from bishoputils import *

def CumuloParameters():
    return {  "-cumuloNbGenerations":2, 
              "-cumuloFjump":(1.0/3.2), 
              "-cumuloRoughness":0.75, 
              "-cumuloNoisescale":1.0, 
              "-cumuloNoiseGamma":0.3, 
              "-cumuloNxNyNz":[100, 100, 100], 
              "-cumulollc":[-5,-5,-5],
              "-cumulolength":[10,10,10],
              "-cumuloClearance":0.0  }




def CumuloNoiseLayers( parameters, basenoise ):
    noiseArray = []
    freqScale = 1.0
    ampScale = 1.0
    fjump = float(parameters["-cumuloFjump"])
    rough = float(parameters["-cumuloRoughness"])
    gamma = float(parameters["-cumuloNoiseGamma"])
    nbGenerations = int(parameters["-cumuloNbGenerations"])
    for i in range(0,nbGenerations):
        print "Cumulo Noise Generation: " + str(i) + " freq=" + str(freqScale) + " amp=" + str(ampScale)
        nextNoise = scale( basenoise, Vector(freqScale,freqScale,freqScale) )
        nextNoise = pow( nextNoise, gamma )
        nextNoise = nextNoise * constant(ampScale)
        noiseArray.append( nextNoise )
        ampScale = ampScale * rough
        freqScale = freqScale * fjump
    return noiseArray


def ApplyCumulo( baseSDF, noiseArray, parameters ):
    sdf = baseSDF
    llc = Vector(  float(parameters["-cumulollc"][0]), float(parameters["-cumulollc"][1]),float(parameters["-cumulollc"][2]) )
    urc = llc + Vector(  float(parameters["-cumulolength"][0]), float(parameters["-cumulolength"][1]),float(parameters["-cumulolength"][2]) )
    dx = Vector(float(parameters["-cumulolength"][0]) / float(parameters["-cumuloNxNyNz"][0]), 
                float(parameters["-cumulolength"][1]) / float(parameters["-cumuloNxNyNz"][1]), 
                float(parameters["-cumulolength"][2]) / float(parameters["-cumuloNxNyNz"][2])   )
    gb = makeGridBox( llc, urc, dx )
    for noiz in noiseArray:
        print "Cumulo Pyro Generation: " + str(noiz)
        cptgrid = makeGrid( gb, Vector(0,0,0) )
        stamp( cptgrid, sdf*grad(sdf), 1 )
        cpt = identity() - gridded(cptgrid)
        wnoise = warp( noiz, cpt )
        sdf = sdf + wnoise
        print "Cumulo Field: " + str(sdf)
    return sdf
