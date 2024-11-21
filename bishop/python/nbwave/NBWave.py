
#bishoputils is the volume modeling system I wrote to assist teaching students volume methods
import bishoputils as bsp
import math



#Some utilities for looking at slices and statistics of data
def slice_state( state ):
	gb = state["gridbox"]
	dx = bsp.dx(gb)
	nx = bsp.nx(gb)
	llc = bsp.llc(gb)
	urc = bsp.urc(gb)
	center = (urc+llc)/2.0
	xx = llc.X()
	for i in range(0,nx):
		P = bsp.Vector( xx, center.Y(), center.Z() )
		phi = bsp.evaluate( state["phi"], P )
		xe = bsp.evaluate( state["XE"], P )
		se = bsp.evaluate( state["SE"], P )
		gphi = bsp.evaluate( state["Gphi"], P )
		print str(xx) + "     " + str(phi) + "     " + str(xe) + "     " + str(se) + "     " + str(gphi)
		xx = xx + dx


def slice_data( data, gb, offset ):
	dx = bsp.dx(gb)
	nx = bsp.nx(gb)
	llc = bsp.llc(gb)
	urc = bsp.urc(gb)
	center = (urc+llc)/2.0
	xx = llc.X()
	for i in range(0,nx):
		P = bsp.Vector( xx, center.Y(), center.Z() ) + offset
		v = bsp.evaluate( data, P )
		print str(xx) + "     " + str(v)
		xx = xx + dx



def statistics_state( state ):
	gb = state["gridbox"]
	print "\nPhi statistics:\n"
	bsp.fieldStatistics( state["phi"], gb )
	print "\nSE statistics:\n"
	bsp.fieldStatistics( state["SE"], gb )
	print "\nXE statistics:\n"
	bsp.fieldStatistics( state["XE"], gb )
	print "\nGphi statistics:\n"
	bsp.fieldStatistics( state["Gphi"], gb )
	print "\n\n\n"
		


#Very simple erosion generator
def cheesy_levelset_filter( ls, gb, erode ):
	lsremove = erode * bsp.dx(gb)
	thinls = ls - bsp.constant(lsremove)
	dr = float( int( lsremove / bsp.dx(gb) ) + 0.5 ) * bsp.dx(gb)
	fatgb = bsp.makeGridBox( bsp.llc(gb), bsp.urc(gb), bsp.Vector(dr,dr,dr) )
	fatgrid = bsp.makeGrid( fatgb, 0.0 )
	bsp.stamp( fatgrid, thinls, 1 )
	lsgrid = bsp.makeGrid( gb, 0.0 )
	bsp.stamp( lsgrid, bsp.gridded(fatgrid) + bsp.constant(lsremove), 1 )
	return bsp.gridded(lsgrid)



#Input parameters for nbwave simulation
def nbwave_parameters():
	return { "-grad_width":3, 
                 "-dx":0.1, 
                 "-dampen":0.0, 
                 "-solver":"garrett1", 
                 "-output":1, 
                 "-esize":1, 
                 "-igthreshold":0.1,
                 "-amplify":1.0,
                 "-blurs":0,
                 "-st_sigma":0.0,
                 "-st_alpha":1.0,
                 "-erode":0
               }


#creates a dictionary with the state data for nbwave simulation
def build_state( phi, xe, se, sb, gb, gw, grav ):
	gsb = bsp.fdboundedgrad( sb, int(gw), gb )
	nbgrid = bsp.makeGrid( gb, bsp.Vector(0,-1,0) )
	bsp.stamp( nbgrid, gsb, 1 )
	nb = bsp.gridded(nbgrid)
	gsb = bsp.fdboundedgrad( sb+se, int(gw), gb )
	nbgrid2 = bsp.makeGrid( gb, bsp.Vector(0,-1,0) )
	bsp.stamp( nbgrid2, gsb, 1 )
	nbe = bsp.gridded(nbgrid2)
	gphi = bsp.fdinteriorgrad(phi, int(gw), bsp.dx(gb), bsp.dy(gb), bsp.dz(gb), sb+se )
	gphi = bsp.RMIncompressibleGradient(gphi, gb, sb+se, nbe, nb, 1.0,  int(gw), 1 )
	state = { "phi":phi, 
                  "XE":xe, 
                  "SE":se, 
                  "SB":sb, 
	          "NB":nb,
                  "gridbox":gb, 
                  "dampen":0.0,
                  "grad_width":gw, 
                  "Gphi":gphi,
                  "gravity":grav
                }
	return state



#After lots of procedural math, make sure everything has been stamped back into grid form
def gridify_state( state ):
	gb = state["gridbox"]
	phi = bsp.makeGrid( gb, 0.0 )
	bsp.stamp( phi, state["phi"], 1 )
	state["phi"] = bsp.gridded(phi)
	xe = bsp.makeGrid( gb, bsp.Vector(0,0,0) )
	bsp.stamp( xe, state["XE"], 1 )
	state["XE"] = bsp.gridded(xe)
	se = bsp.makeGrid( gb, 0.0 )
	bsp.stamp( se, state["SE"], 1 )
	state["SE"] = bsp.gridded(se)
	return state


#Generates a field that is 1 interior of grid, 0 outside, and ramps between 0 and 1 for a thin region
def boundary_shell( gb, thickness ):
	llc = bsp.llc(gb) 
	urc = bsp.urc(gb) - bsp.Vector( bsp.dx(gb), bsp.dy(gb), bsp.dz(gb) )
	bs = bsp.HardBox( llc, urc )
	bs = bsp.clamp( bs/bsp.constant(thickness), 0.0, 1.0 )
	return bs


#Shell region for surface tension
def smooth_shell( f, dx ):
	outside = bsp.constant(0.0)
	test = bsp.constant(1.0) - bsp.abs( f / bsp.constant(dx) )
	inside =  bsp.constant(0.5)*( bsp.constant(1.0) + bsp.cos( f * bsp.constant( 3.14159265*0.5/dx ) ) )
	ss = bsp.which( inside, outside, test )
	return ss

#Construct surface tension force field
def force_surface_tension( sb, se, sigma, alpha, dx, gb, gw ):
	if sigma == 0.0:
		return bsp.constant(0.0)
	nb =  bsp.unitvector(bsp.fdboundedgrad(sb, gw, gb))
	nbe =  bsp.unitvector(bsp.fdboundedgrad(sb+se, gw, gb) )
	f = bsp.constant(sigma) * bsp.fdboundeddiv(nbe - bsp.constant(alpha)*nb, gw, gb) * smooth_shell(sb+se, dx)
	fgrid = bsp.makeGrid( gb, 0.0 )
	bsp.stamp( fgrid, f, 1 )
	f = bsp.gridded( fgrid )
	return f


#
#  Coefficients for the Garrett solvers come from
#
#  "Fast Polynomial Approximations to Sine and Cosine"
#  by  Charles K Garrett, appendices A and A.2
#
def garrett_solver_1( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [- 9.33876972831853619134e-02, + 8.56983327795249506462e-01]
	ccoeff = [- 2.30984600730397541756e-01, + 7.59908877317533285829e-01]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )


def garrett_solver_2( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [+ 5.64311797634681035370e-03, - 1.55271410633428644799e-01, + 9.87862135574673806965e-01]
	ccoeff = [+ 2.61598203821710351549e-02, - 4.52287810766610989686e-01, + 9.78326390892394782255e-01]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )
	

def garrett_solver_3( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [- 1.47740880797318521837e-04,
                  + 7.99858143743132551201e-03,
                  - 1.65838452698030873892e-01,
                  + 9.99450193893262089505e-01
                 ]
	ccoeff = [- 9.92863295193013173583e-04,
                  + 3.95223221293306431394e-02,
                  - 4.96248679451054559990e-01,
                  + 9.98987171037332669123e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )
	


def garrett_solver_4( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [+ 2.17326217498596729611e-06,
                  - 1.93162796407356830500e-04,
                  + 8.31238887417884598346e-03,
                  - 1.66632595072086745320e-01,
                  + 9.99984594193494365437e-01
                 ]
	ccoeff = [+ 1.90652668840074246305e-05,
                  - 1.34410769349285321733e-03,
                  + 4.15223086250910767516e-02,
                  - 4.99837602272995734437e-01,
                  + 9.99971094606182687341e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )



def garrett_solver_5( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [- 2.05342856289746600727e-08,
                  + 2.70405218307799040084e-06,
                  - 1.98125763417806681909e-04,
                  + 8.33255814755188010464e-03,
                  - 1.66665772196961623983e-01, 
                  + 9.99999707044156546685e-01
                 ]
	ccoeff = [- 2.21941782786353727022e-07,
                  + 2.42532401381033027481e-05,
                  - 1.38627507062573673756e-03,
                  + 4.16610337354021107429e-02,
                  - 4.99995582499065048420e-01,
                  + 9.99999443739537210853e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )




def garrett_solver_6( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [+ 1.35333825545218599272e-10,
                  - 2.47016425480527869032e-08,
                  + 2.75322955330449911163e-06,
                  - 1.98403112669018996690e-04,
                  + 8.33331451433080749755e-03,
                  - 1.66666650437066346286e-01,
                  + 9.99999995973569972699e-01
                 ]
	ccoeff = [+ 1.73691489450821293670e-09,
                  - 2.71133771940801138503e-07,
                  + 2.47734245730930250260e-05,
                  - 1.38879704270452054154e-03,
                  + 4.16665243677686230461e-02,
                  - 4.99999917728614591900e-01,
                  + 9.99999992290827491711e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )


def garrett_solver_7( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [- 6.58075489175121657026e-13,
                  + 1.58850004791504823423e-10,
                  - 2.50368914392103083120e-08,
                  + 2.75565598752102704008e-06,
                  - 1.98412483604340805859e-04,
                  + 8.33333301181570639096e-03,
                  - 1.66666666451352167974e-01,
                  + 9.99999999958141380079e-01
                 ]
	ccoeff = [- 9.77507131527006498114e-12,
                  + 2.06207503915813519567e-09,
                  - 2.75369918573799545860e-07,
                  + 2.48006913718665260256e-05,
                  - 1.38888674687691339750e-03,
                  + 4.16666641590361985136e-02,
                  - 4.99999998886526927002e-01,
                  + 9.99999999919365479957e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )




def garrett_solver_8( state, dt, threshfactor, nbblurs, sigma, alpha ):
	scoeff = [+ 2.45928524290153002259e-15,
                  - 7.58106260511361554811e-13,
                  + 1.60521984385153059172e-10,
                  - 2.50516861359706378210e-08,
                  + 2.75573034843986111280e-06,
                  - 1.98412694971242118241e-04,
                  + 8.33333332926687803703e-03,
                  - 1.66666666664489411560e-01,
                  + 9.99999999999659411867e-01
                 ]
	ccoeff = [+ 4.14869721869947572436e-14,
                  - 1.13600777795958675706e-11,
                  + 2.08661897358261903687e-09,
                  - 2.75567298437160383039e-07,
                  + 2.48015679993921751541e-05,
                  - 1.38888885344276371809e-03,
                  + 4.16666666341518636873e-02,
                  - 4.99999999988560571910e-01,
                  + 9.99999999999340745485e-01
                 ]
	return garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha )


	
def garrett_solver_N( state, dt, threshfactor, scoeff, ccoeff, nbblurs, sigma, alpha ):
	phi0 = state["phi"]
	xe0 = state["XE"]
	se0 = state["SE"]
	gb = state["gridbox"]
	surface = state["SB"] + se0
	gsb = bsp.fdboundedgrad( surface, int(state["grad_width"]), gb )
	nbgrid = bsp.makeGrid( gb, bsp.Vector(0,-1,0) )
	bsp.stamp( nbgrid, gsb, 1 )
	nb = bsp.gridded(nbgrid)
	nb = state["NB"]
	surface = state["SB"]
	gw = int( state["grad_width"])
	gravity = state["gravity"]
	threshold = 9.8 * dt * threshfactor
	N = len(scoeff)
	for i in range(0,N):
		filtercoeff = threshold*math.fabs(scoeff[i])
		if i==0:
			state["phi"] = phi0*bsp.constant(ccoeff[i]) - bsp.constant(dt*scoeff[i])*gravity*xe0 
			state["XE"]  = xe0 *bsp.constant(ccoeff[i]) - bsp.constant(dt*scoeff[i])*state["Gphi"]
		else:
			state["phi"] = phi0*bsp.constant(ccoeff[i]) - bsp.constant(dt)*gravity*(xe0*bsp.constant(scoeff[i]) + bsp.constant(dt)*state["Gphi"]) 
			aa = bsp.fdboundedgrad(   bsp.constant(scoeff[i])*phi0 + bsp.constant(dt)*gravity*state["XE"],
                                           int(state["grad_width"]), gb )
			aa = bsp.RMIncompressibleGradient(aa, gb, surface, nb, state["NB"], filtercoeff,  int(state["grad_width"]), nbblurs )
			state["XE"]  = xe0 *bsp.constant(ccoeff[i]) - bsp.constant(dt)*aa
		gphi = bsp.fdboundedgrad( state["phi"], int(state["grad_width"]), gb  )
		state["Gphi"] = bsp.RMIncompressibleGradient(gphi, gb, surface, nb, state["NB"], filtercoeff, int(state["grad_width"]), nbblurs )
		state = gridify_state(state)
	if sigma != 0.0:
		print "\n\tSURFACE TENSION\n"
		dx = bsp.dx(gb) * 2.5
		fst = force_surface_tension( state["SB"], se0, sigma, alpha, dx, gb, gw )
        	state["phi"] += bsp.constant(dt)*fst	
		gphi = bsp.fdboundedgrad( state["phi"], int(state["grad_width"]), gb  )
		state["Gphi"] = bsp.RMIncompressibleGradient(gphi, gb, surface, nb, state["NB"], filtercoeff, int(state["grad_width"]), nbblurs )
		print "\n\tDONE SURFACE TENSION\n"
	state["SE"]  = se0 - bsp.constant(dt)*(state["NB"]*state["Gphi"])
	state = gridify_state(state)
	return state



#Selects the correct solver based on command line parameters, and handles subframe stepping
def solve( state, dt, parameters, nbsubstep ):
	ddt = dt/float(nbsubstep)
	igthreshold = float(parameters["-igthreshold"])
	nbblurs = int(parameters["-blurs"])
	sigma = float(parameters["-st_sigma"])
	alpha = float(parameters["-st_alpha"])
	for i in range(0,nbsubstep):
		if parameters["-solver"] == "garrett1":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold)
			state = garrett_solver_1( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett2":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold) 
			state = garrett_solver_2( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett3":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold)
			state = garrett_solver_3( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett4":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold)
			state = garrett_solver_4( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett5":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold) 
			state = garrett_solver_5( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett6":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold)
			state = garrett_solver_6( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett7":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold) 
			state = garrett_solver_7( state, ddt, igthreshold, nbblurs, sigma, alpha )
		elif parameters["-solver"] == "garrett8":
			print "Solver: " + parameters["-solver"] + " with threshold: " + str(igthreshold) 
			state = garrett_solver_8( state, ddt, igthreshold, nbblurs, sigma, alpha )
		else:
			print "Solver: " + parameters["-solver"] + " is not a valid choice."
	return state



