
import math
import bishop
import bishoputils as bu




# From 
# Back and forth error compensation and correction mehtods for semi-lagrangian schemes with application to level set interface computations
# TODD F. DUPONT AND YINGJIE LIU
def advect_bfecc( field, velocity, timestep ):
	step1 = bishop.advect( field, velocity, timestep )
	step2 = bishop.advect( step1, velocity, -timestep )
	step3 = field * bishop.constant(1.5) - step2 * constant(0.5)
	step4 = bishop.advect( step3, velocity, timestep )
	return step4


# From 
# An Unconditionally Stable MacCormack Method
# Andrew Selle	Ronald Fedkiw	ByungMoon Kim Yingjie Liu	Jarek Rossignac
def advect_modified_maccormack( field, velocity, timestep ):
	step1 = bishop.advect( field, velocity, timestep )
	step2 = bishop.advect( step1, velocity, -timestep )
	step3 = step1 + (field - step2) * bishop.constant(0.5)
	return step3




def nested_warp_field( base_warp, nb ):
    nested_warp = base_warp
    for i in range(0,nb):
        nested_warp = bishop.warp( nested_warp, nested_warp )
    return nested_warp




def log_advect( field, velocity, dt, nb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.advect(bishop.identity(), velocity, ddt )
    log_warp_field = nested_warp_field( base_warp, nb )
    return bishop.warp( field, log_warp_field )

def log_advect_bfecc( field, velocity, dt, nb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_bfecc(bishop.identity(), velocity, ddt )
    log_warp_field = nested_warp_field( base_warp, nb )
    return bishop.warp( field, log_warp_field )

def log_advect_modified_maccormack( field, velocity, dt, nb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_modified_maccormack(bishop.identity(), velocity, ddt )
    log_warp_field = nested_warp_field( base_warp, nb )
    return bishop.warp( field, log_warp_field )

def log_advect_gradient_stretch( field, velocity, dt, nb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.gradientStretchCM(velocity, ddt, 1 )
    log_warp_field = nested_warp_field( base_warp, nb )
    return bishop.warp( field, log_warp_field )






def gridded_nested_warp_field( base_warp, nb, gb ):
    #bu.DeclareIt( "Gridded Nested Warp Field computation STARTED" )
    griddata = bishop.makeGrid( gb, bishop.Vector(0,0,0) )
    bishop.stamp( griddata, base_warp-bishop.identity(), 1 )
    nested_warp = bishop.gridded(griddata) + bishop.identity()
    for i in range(0,nb):
        nested_warp = bishop.warp( nested_warp, nested_warp )
        griddata = bishop.makeGrid( gb, bishop.Vector(0,0,0) )
        bishop.stamp( griddata, nested_warp-bishop.identity(), 1 )
        nested_warp = bishop.gridded(griddata) + bishop.identity()
    #bu.DeclareIt("Gridded Nested Warp Field computation DONE" )
    return nested_warp



def gridded_log_advect( field, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.advect(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    return bishop.warp( field, log_warp_field )

def gridded_log_advect_bfecc( field, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_bfecc(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    return bishop.warp( field, log_warp_field )

def gridded_log_advect_modified_maccormack( field, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_modified_maccormack(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    return bishop.warp( field, log_warp_field )

def gridded_log_advect_gradient_stretch( field, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.gradientStretchCM(velocity, ddt, 1 )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    return bishop.warp( field, log_warp_field )





def gridded_log_multi_advect( fields, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.advect(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    advected_fields = []
    for f in fields:
        advected_fields.append( bishop.warp( f, log_warp_field ) )
    return advected_fields

def gridded_log_multi_advect_bfecc( fields, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_bfecc(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    advected_fields = []
    for f in fields:
        advected_fields.append( bishop.warp( f, log_warp_field ) )
    return advected_fields

def gridded_log_multi_advect_modified_maccormack( fields, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = advect_modified_maccormack(bishop.identity(), velocity, ddt )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    advected_fields = []
    for f in fields:
        advected_fields.append( bishop.warp( f, log_warp_field ) )
    return advected_fields

def gridded_log_multi_advect_gradient_stretch( fields, velocity, dt, nb, gb ):
    ddt = float(dt) / math.pow( 2, int(nb) )
    base_warp = bishop.gradientStretchCM(velocity, ddt, 1 )
    log_warp_field = gridded_nested_warp_field( base_warp, nb, gb )
    advected_fields = []
    for f in fields:
        advected_fields.append( bishop.warp( f, log_warp_field ) )
    return advected_fields



