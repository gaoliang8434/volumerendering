#!/usr/bin/python

from bishoputils import *


def YAdvectField( ufield, timestep, gridbox ):
	gu = grad(ufield)
	egu = exp(gu*constant(timestep))
	Y = identity() - ufield*constant(timestep)
	dX = Vector( dx(gridbox), 0, 0 )
	dY = Vector( 0, dy(gridbox), 0 )
	dZ = Vector( 0, 0, dz(gridbox) )
	wpx = identity() - constant( dX )
	wpy = identity() - constant( dY )
	wpz = identity() - constant( dZ )
	upx = identity() - constant( dX*2.0 )
	upy = identity() - constant( dY*2.0 )
	upz = identity() - constant( dZ*2.0 )
	wmx = identity() + constant( dX )
	wmy = identity() + constant( dY )
	wmz = identity() + constant( dZ )
	umx = identity() + constant( dX*2.0 )
	umy = identity() + constant( dY*2.0 )
	umz = identity() + constant( dZ*2.0 )
	for i in range(0,25):
		print "Y advect iteration " + str(i)
		homogeneous = warp(Y,upx) + warp(Y,upy) + warp(Y,upz) + warp(Y,umx) + warp(Y,umy) + warp(Y,umz)
		source = constant(dX)*warp(egu,warp(Y,wmx)) + constant(dY)*warp(egu,warp(Y,wmy)) + constant(dZ)*warp(egu,warp(Y,wmz)) - constant(dX)*warp(egu,warp(Y,wpx)) - constant(dY)*warp(egu,warp(Y,wpy)) - constant(dZ)*warp(egu,warp(Y,wpz))
		Ynew = homogeneous*constant(1.0/6.0) + source*constant(1.0/3.0)
		Ygrid = makeGrid( gridbox, Vector(0,0,0) )
		stamp( Ygrid, Ynew-identity() )
		Y = gridded(Ygrid) + identity()
	return Y


def StandardColorRender( grid, volum, litcolor, parms ):
	print "Executing Standard Color Render"
	renderStartTime = formattedTime()
	ambient = constant( Color( float(parms["-ambient"][0]), float(parms["-ambient"][1]), float(parms["-ambient"][2]), 0.0 ) )
	data = MakeSparseRenderData( grid, volum, volum, litcolor, ambient, parms )
	dsms = MakeDSMs( data, parms )
	packagedDsms = PackageDSMs( dsms, data )
	image = RenderVolume( data, parms )
	renderEndTime = formattedTime()
	nparms = MergeParameters( parms, { "RenderStartTime":renderStartTime, "RenderEndTime":renderEndTime } )
	OutputImage( image, nparms )



parms = AllParameters( {} )
parms["-size"] = [960, 540]
parms["-aspect"] = float( parms["-size"][0] ) / float( parms["-size"][1] )
parms["-near"] = 5.0
parms["-maxpathlength"] = 15.0
parms["-ds"] = 0.01
parms["-fov"] = 60.0
parms["-view"] = [0,4,0]
parms["-eye"] = [0,4,20]
parms["-dsmllc"] = [-1,-1,-1]
parms["-dsmlength"] = [2,8,2]
parms["-dsmsize"] = [64,256,64]
parms = CmdLineFindValues( parms )
PrintParameters(parms)


nx = 128
nz = nx
ny = 4*nx

simllc = Vector(-1,-1,-1)
simsize = Vector(2,8,2)
simurc = simllc + simsize
simres = Vector( simsize.X()/nx, simsize.Y()/ny, simsize.Z()/nz )
simgb = makeGridBox( simllc, simurc, simres )



radius = 0.3
x0 = (identity() - constant( Vector(0,0,0) ))/constant(radius)
initialdensity = clamp( ( constant(1.0) - (x0*x0) )/constant(0.03), 0.0, 1.0 )


slavelocity = constant(Vector(0,0,0))
selmavelocity = constant(Vector(0,0,0))
gavelocity = constant(Vector(0,0,0))
yvelocity = constant(Vector(0,0,0))

X = identity()
XY = identity()

sladensity = initialdensity
selmadensity = initialdensity
gadensity = initialdensity
ydensity = initialdensity

gravity = constant(Vector(0,-10,0))

slacolor = Color(1,0,0,1)
selmacolor = Color(0,1,0,1)
gacolor =Color(0,0,1,1)
ycolor =Color(0,0,1,1)

nbframes = 100
dt = 0.1
for i in range(0,nbframes+1):
	######### SLA ##########
	sladensity = advect( sladensity, slavelocity, dt )
	densityGrid = makeGrid( simgb, 0.0 )
	stamp( densityGrid, sladensity )
	sladensity = gridded(densityGrid)
	slavelocity = advect( slavelocity, slavelocity, dt )
	slavelocity += -gravity*sladensity*constant(dt)
	slavelocity = FFTDivFree( simgb, slavelocity )
	slacf = constant( slacolor ) * sladensity
	########## SELMA #######
	selmadensity = warp( initialdensity, X )
	selmavelocity = advect( selmavelocity, selmavelocity, dt )
	selmavelocity += -gravity*selmadensity*constant(dt)
	selmavelocity = FFTDivFree( simgb, selmavelocity )
	selmacf = constant( selmacolor ) * selmadensity
	X = advect( X, selmavelocity, dt )
	Xgrid = makeGrid( simgb, Vector(0,0,0) )
	stamp( Xgrid, X-identity())
	X = identity() + gridded(Xgrid)
	########## Y #######
	Y = YAdvectField( yvelocity, dt, simgb )
	XY = warp( XY, Y )
	XYgrid = makeGrid( simgb, Vector(0,0,0) )
	stamp( XYgrid, XY-identity())
	XY = identity() + gridded(XYgrid)
	ydensity = warp( initialdensity, XY )
	yvelocity = warp( yvelocity, Y )
	yvelocity += -gravity*ydensity*constant(dt)
	yvelocity = FFTDivFree( simgb, yvelocity )
	ycf = constant( ycolor ) * ydensity
	########## GA ##########
	#gadensity = advect(gadensity, gavelocity, dt )
	#gavelocity = advect( gavelocity, gavelocity, dt )
	#gavelocity += -gravity*gadensity*constant(dt)
	#gavelocity = FFTDivFree( simgb, gavelocity )
	#gacf = constant( gacolor ) * gadensity
	####### Assemble ######
	parms["-turntableframe"] = i+1
	densityfield = sladensity + selmadensity + ydensity
	dengrid = makeGrid( simgb, 0.0 )
	stamp(dengrid, densityfield )
	colorfield = slacf + selmacf + ycf
	StandardColorRender( dengrid, densityfield, colorfield, parms )
