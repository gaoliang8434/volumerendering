
docDocumentationList = []


def newDocumentationItem( item, ret, sig, desc, seealso ):
	docItem = {
	"item":item,
	"return":ret,
	"signature":sig,
	"description":desc,
	"seealso":seealso
	}
	docDocumentationList.append( docItem.copy() )

def bishopdocLongList( docitem ):
	listing = docitem["return"] + " " + docitem["item"] + "( "
	signature = docitem["signature"]
	for i in range(0,len(signature)):
		listing += signature[i] + " (" + str(i+1) + ")"
		if i < len(signature)-1:
			listing += ", "
	listing +=  " )"
	listing += "\n\n"
	listing += docitem["description"] + "\n\n"
	if docitem["seealso"] != "":
		listing += "See also: " + docitem["seealso"] + "\n\n"
	return listing




def bishophelp( item ):
	nbFoundItems = 0
	for docitem in docDocumentationList:
		if str(docitem["item"]).lower().find( str(item).lower() ) >= 0:
			nbFoundItems += 1
			listing = "--------------------------------------------------------------\n\n"
			listing += bishopdocLongList( docitem )
			print(listing)
	if nbFoundItems > 0:
		print("--------------------------------------------------------------\n\n")



def bishopsearch( item ):
	nbFoundItems = 0
	for docitem in docDocumentationList:
		if str(docitem["item"]).lower().find( str(item).lower() ) >= 0 or str(docitem["description"]).lower().find( str(item).lower() ) >= 0:
			nbFoundItems += 1
			listing = "--------------------------------------------------------------\n\n"
			listing += bishopdocLongList( docitem )
			print(listing)
	if nbFoundItems > 0:
		print("--------------------------------------------------------------\n\n")





IN = "integer"
FL = "float"
VE = "Vector"
VO = "No Retrun"
SF = "ScalarField"
VF = "VectorField"
MF = "MatrixField"
CF = "ColorField"
FF = "FormField"
CO = "Color"
MA = "Matrix"
FO = "Form"
ST = "String"
SG = "ScalarGrid"
VG = "VectorGrid"
CG = "ColorGrid"
MG = "MatrixGrid"
GB = "GridBox"


newDocumentationItem( "evaluate", FL, [SF, VE], 
          "Evaluates the scalar field (1) at the spatial position (2) and returns the value.",
          ""
          )

newDocumentationItem( "evaluate", VE, [VF, VE], 
          "Evaluates the vector field (1) at the spatial position (2) and returns the value.",
          ""
          )


newDocumentationItem( "evaluate", CO, [CF, VE], 
          "Evaluates the color field (1) at the spatial position (2) and returns the value.",
          ""
          )

newDocumentationItem( "evaluate", MA, [MF, VE], 
          "Evaluates the matrix field (1) at the spatial position (2) and returns the value.",
          ""
          )

newDocumentationItem( "evaluate", FO, [FF, VE], 
          "Evaluates the form field (1) at the spatial position (2) and returns the value.",
          ""
          )

newDocumentationItem( "grad", VF, [SF], 
          "Returns a vector field of the gradient of the scalar field (1)",
          "curl, div"
          )

newDocumentationItem( "grad", MF, [VF], 
          "Returns a matrix field of the gradient of the vector field (1)",
          ""
          )

newDocumentationItem( "grad", FF, [FF], 
          "Returns a form field of the exterior derivative of the form field (1)",
          ""
          )

newDocumentationItem( "constant", SF, [FL], 
          "Returns a scalar field that always evaluates to the value of the float (1)",
          ""
          )

newDocumentationItem( "constant", VF, [VE], 
          "Returns a vector field that always evaluates to the value of the vector (1)",
          ""
          )

newDocumentationItem( "constant", MF, [MA], 
          "Returns a matrix field that always evaluates to the value of the matrix (1)",
          ""
          )

newDocumentationItem( "constant", CF, [CO], 
          "Returns a color field that always evaluates to the value of the color (1)",
          ""
          )

newDocumentationItem( "constant", FF, [FO], 
          "Returns a form field that always evaluates to the value of the form (1)",
          ""
          )

newDocumentationItem( "scale", SF, [SF,VE], 
          "Generates a scalar field using a scaling transformation on the scalar field (1).\n" +
	  "The evaluation position is scaled in each direction by the scaling parameters in (2).",
          "translate, rotate"
          )

newDocumentationItem( "scale", VF, [VF,VE], 
          "Generates a vector field using a scaling transformation on the vector field (1).\n" +
	  "The evaluation position is scaled in each direction by the scaling parameters in (2).\n"+
	  "Inverse scaling is applied to the value of the vector field.",
          "translate, rotate"
          )

newDocumentationItem( "scale", CF, [CF,VE], 
          "Generates a color field using a scaling transformation on the color field (1).\n" +
	  "The evaluation position is scaled in each direction by the scaling parameters in (2).",
          "translate, rotate"
          )

newDocumentationItem( "scale", FF, [FF,VE], 
          "Generates a form field using a scaling transformation on the form field (1).\n" +
	  "The evaluation position is scaled in each direction by the scaling parameters in (2).",
          "translate, rotate"
          )


newDocumentationItem( "translate", SF, [SF,VE], 
          "Generates a scalar field using a translation by (2) on the scalar field (1).",
          "scale, rotate"
          )

newDocumentationItem( "translate", VF, [VF,VE], 
          "Generates a vector field using a translation by (2) on the vector field (1).",
          "scale, rotate"
          )

newDocumentationItem( "translate", CF, [CF,VE], 
          "Generates a color field using a translation by (2) on the color field (1).",
          "scale, rotate"
          )

newDocumentationItem( "translate", FF, [FF,VE], 
          "Generates a form field using a translation by (2) on the form field (1).",
          "scale, rotate"
          )




newDocumentationItem( "rotate", SF, [SF,VE], 
          "Generates a scalar field using a rotation by (2) on the scalar field (1).\n" +
	  "The direction of the Vector (2) is the axis of rotation, and its magnitude\n" +
	  "is the angle of rotation in radians.",
          "scale, translate"
          )

newDocumentationItem( "rotate", VF, [VF,VE], 
          "Generates a vector field using a rotation by (2) on the vector field (1).\n" +
	  "The direction of the Vector (2) is the axis of rotation, and its magnitude\n" +
	  "is the angle of rotation in radians.",
          "scale, translate"
          )

newDocumentationItem( "rotate", CF, [CF,VE], 
          "Generates a color field using a rotation by (2) on the color field (1).\n" +
	  "The direction of the Vector (2) is the axis of rotation, and its magnitude\n" +
	  "is the angle of rotation in radians.",
          "scale, translate"
          )

newDocumentationItem( "rotate", FF, [FF,VE], 
          "Generates a form field using a rotation by (2) on the form field (1).\n" +
	  "The direction of the Vector (2) is the axis of rotation, and its magnitude\n" +
	  "is the angle of rotation in radians.",
          "scale, translate"
          )


newDocumentationItem( "rotation", MF, [VE], 
          "Generates a matrix field of of a rotation by vector field in (1).\n" +
	  "The direction of the vector field (1) is the axis of rotation, and\n" +
	  "its magnitude is the angle of rotation in radians.",
          "scale, translate, rotate"
          )


newDocumentationItem( "exp", SF, [SF], 
          "Generates a scalar field that is the exponential of the values of the scalar field (1).",
          ""
          )

newDocumentationItem( "exp", MF, [MF], 
          "Generates a matrix field that is the exponential of the values of the matrix field (1).",
          ""
          )


newDocumentationItem( "det", SF, [MF], 
          "Generates a scalar field that is the determinant of the values of the matrix field (1).",
          ""
          )


newDocumentationItem( "report", SF, [SF,ST], 
          "Generates a scalar field that returns the same value as the scalar field (1),\n" +
	  "but also prints out the string (2), location of evaluation, and the evaluated value.",
          ""
          )

newDocumentationItem( "report", VF, [VF,ST], 
          "Generates a vector field that returns the same value as the vector field (1),\n" +
	  "but also prints out the string (2), location of evaluation, and the evaluated value.",
          ""
          )

newDocumentationItem( "report", CF, [CF,ST], 
          "Generates a color field that returns the same value as the color field (1),\n" +
	  "but also prints out the string (2), location of evaluation, and the evaluated value.",
          ""
          )

newDocumentationItem( "report", FF, [FF,ST], 
          "Generates a form field that returns the same value as the form field (1),\n" +
	  "but also prints out the string (2), location of evaluation, and the evaluated value.",
          ""
          )


newDocumentationItem( "negate", SF, [SF], 
          "Generates a scalar field that returns the negative of the scalar field (1).",
          ""
          )

newDocumentationItem( "negate", VF, [VF], 
          "Generates a vector field that returns the negative of the vector field (1).",
          ""
          )

newDocumentationItem( "negate", CF, [CF], 
          "Generates a color field that returns the negative of the color field (1).",
          ""
          )

newDocumentationItem( "negate", MF, [MF], 
          "Generates a matrix field that returns the negative of the matrix field (1).",
          ""
          )

newDocumentationItem( "negate", FF, [FF], 
          "Generates a form field that returns the negative of the form field (1).",
          ""
          )


newDocumentationItem( "abs", SF, [SF], 
          "Generates a scalar field that returns the absolute value of the scalar field (1).",
          ""
          )

newDocumentationItem( "abs", SF, [VF], 
          "Generates a scalar field that returns the magnitude of the vector field (1).",
          ""
          )

newDocumentationItem( "which", SF, [SF,SF,SF], 
          "Generates a scalar field that evaluates to the value of either scalar field (2)\n" +
	  "or scalar field (3) depending on the value of scalar field (1). A positive value\n" +
	  "for (1) chooses the value from (2), negative chooses the value from (3). The unchosen\n" +
	  "field is not evaluated.",
          ""
          )

newDocumentationItem( "which", VF, [SF,VF,VF], 
          "Generates a vector field that evaluates to the value of either vector field (2)\n" +
	  "or vector field (3) depending on the value of scalar field (1). A positive value\n" +
	  "for (1) chooses the value from (2), negative chooses the value from (3). The unchosen\n" +
	  "field is not evaluated.",
          ""
          )

newDocumentationItem( "which", CF, [SF,CF,CF], 
          "Generates a color field that evaluates to the value of either color field (2)\n" +
	  "or color field (3) depending on the value of scalar field (1). A positive value\n" +
	  "for (1) chooses the value from (2), negative chooses the value from (3). The unchosen\n" +
	  "field is not evaluated.",
          ""
          )

newDocumentationItem( "which", FF, [SF,FF,FF], 
          "Generates a form field that evaluates to the value of either form field (2)\n" +
	  "or form field (3) depending on the value of scalar field (1). A positive value\n" +
	  "for (1) chooses the value from (2), negative chooses the value from (3). The unchosen\n" +
	  "field is not evaluated.",
          ""
          )



newDocumentationItem( "multiply", SF, [SF,FL], 
          "Generates a scalar field that returns the product of the scalar field (1)\n" +
	  "and the float value (2).",
          ""
          )

newDocumentationItem( "multiply", SF, [SF,SF], 
          "Generates a scalar field that returns the product of the scalar field (1)\n" +
	  "and the scalar field (2).",
          ""
          )

newDocumentationItem( "multiply", VF, [VF,FL], 
          "Generates a vector field that returns the product of the vector field (1)\n" +
	  "and the float value (2).",
          ""
          )

newDocumentationItem( "multiply", VF, [VF,SF], 
          "Generates a vector field that returns the product of the vector field (1)\n" +
	  "and the scalar field (2).",
          ""
          )

newDocumentationItem( "multiply", CF, [CF,FL], 
          "Generates a color field that returns the product of the color field (1)\n" +
	  "and the float value (2).",
          ""
          )

newDocumentationItem( "multiply", CF, [CF,SF], 
          "Generates a color field that returns the product of the color field (1)\n" +
	  "and the scalar field (2).",
          ""
          )

newDocumentationItem( "multiply", CF, [CF,CF], 
          "Generates a color field that returns the product of the color field (1)\n" +
	  "and the color field (2). This is a component-by-component multiplication.",
          ""
          )


newDocumentationItem( "multiply", MF, [MF,FL], 
          "Generates a matrix field that returns the product of the matrix field (1)\n" +
	  "and the float value (2).",
          ""
          )

newDocumentationItem( "multiply", MF, [MF,SF], 
          "Generates a matrix field that returns the product of the matrix field (1)\n" +
	  "and the scalar field (2).",
          ""
          )

newDocumentationItem( "multiply", FF, [FF,FL], 
          "Generates a form field that returns the product of the form field (1)\n" +
	  "and the float value (2).",
          ""
          )

newDocumentationItem( "multiply", FF, [FF,SF], 
          "Generates a form field that returns the product of the form field (1)\n" +
	  "and the scalar field (2).",
          ""
          )


newDocumentationItem( "divide", SF, [SF,FL], 
          "Generates a scalar field that returns the scalar field (1) divided by\n" +
	  "the float value (2).",
          ""
          )

newDocumentationItem( "divide", SF, [SF,SF], 
          "Generates a scalar field that returns the scalar field (1) divided by\n" +
	  "the scalar field (2).",
          ""
          )


newDocumentationItem( "divide", VF, [VF,FL], 
          "Generates a vector field that returns the vector field (1) divided by\n" +
	  "the float value (2).",
          ""
          )

newDocumentationItem( "divide", FF, [VF,SF], 
          "Generates a vector field that returns the vector field (1) divided by\n" +
	  "the scalar field (2).",
          ""
          )


newDocumentationItem( "divide", CF, [CF,FL], 
          "Generates a color field that returns the color field (1) divided by\n" +
	  "the float value (2).",
          ""
          )

newDocumentationItem( "divide", CF, [CF,SF], 
          "Generates a color field that returns the color field (1) divided by\n" +
	  "the scalar field (2).",
          ""
          )

newDocumentationItem( "divide", MF, [MF,FL], 
          "Generates a matrix field that returns the matrix field (1) divided by\n" +
	  "the float value (2).",
          ""
          )

newDocumentationItem( "divide", MF, [MF,SF], 
          "Generates a matrix field that returns the matrix field (1) divided by\n" +
	  "the scalar field (2).",
          ""
          )

newDocumentationItem( "divide", FF, [FF,FL], 
          "Generates a form field that returns the form field (1) divided by\n" +
	  "the float value (2).",
          ""
          )

newDocumentationItem( "divide", FF, [FF,SF], 
          "Generates a form field that returns the form field (1) divided by\n" +
	  "the scalar field (2).",
          ""
          )


newDocumentationItem( "add", SF, [SF,SF], 
          "Generates a scalar field that returns the sum of scalar fields (1) and (2).",
          ""
          )

newDocumentationItem( "add", VF, [VF,VF], 
          "Generates a vector field that returns the sum of vector fields (1) and (2).",
          ""
          )

newDocumentationItem( "add", CF, [CF,CF], 
          "Generates a color field that returns the sum of color fields (1) and (2).",
          ""
          )

newDocumentationItem( "add", MF, [MF,MF], 
          "Generates a matrix field that returns the sum of matrix fields (1) and (2).",
          ""
          )

newDocumentationItem( "add", FF, [FF,FF], 
          "Generates a form field that returns the sum of form fields (1) and (2).",
          ""
          )


newDocumentationItem( "subtract", SF, [SF,SF], 
          "Generates a scalar field that returns the difference of scalar fields (1) and (2).",
          ""
          )

newDocumentationItem( "subtract", VF, [VF,VF], 
          "Generates a vector field that returns the difference of vector fields (1) and (2).",
          ""
          )

newDocumentationItem( "subtract", CF, [CF,CF], 
          "Generates a color field that returns the difference of color fields (1) and (2).",
          ""
          )

newDocumentationItem( "subtract", MF, [MF,MF], 
          "Generates a matrix field that returns the difference of matrix fields (1) and (2).",
          ""
          )

newDocumentationItem( "subtract", FF, [FF,FF], 
          "Generates a form field that returns the difference of form fields (1) and (2).",
          ""
          )

newDocumentationItem( "Sphere", SF, [VE,FL], 
          "Generates a scalar field for the implicit function of a sphere centered at\n" +
	  "position (1) with radius (2).",
          ""
          )

newDocumentationItem( "Ellipse", SF, [VE,VE,FL,FL], 
          "Generates a scalar field for the implicit function of an ellipse centered at\n" +
	  "position (1) with radii (3) and (4).  The axis (2) is the orientation of the\n" +
	  "radius (3).",
          ""
          )

newDocumentationItem( "CsgBox", SF, [VE,FL,FL], 
          "Generates a scalar field for the implicit function of a rounded-corner box\n" +
	  "centered at position (1) with radius (2). The sharpness of the corners depends\n" +
	  "on the value of (3).  If (3) has the value 2, the box is a sphere. Higher values\n" +
	  "better define the corners of the box.",
          "CsgRectangularBox, HardBox"
          )

newDocumentationItem( "CsgRectangularBox", SF, [VE,FL,VE,FL], 
          "Similar to CsgBox, except the sides may be unequal with aspect ratios (3). The box is\n" +
	  "centered at position (1) with radius (2). The sharpness of the corners depends\n" +
	  "on the value of (4).  If (4) has the value 2, the box is an ellipse. Higher values\n" +
	  "better define the corners of the box.",
          "CsgBox, HardBox"
          )

newDocumentationItem( "HardBox", SF, [VE,VE], 
          "Generates a scalar field of the implicit function of a sharp-edged box from the\n" +
	  "intersection of six planes. The bounds of the box are the \"lower left corner\" (1)\n" +
	  "and the \"upper right corner\" (2).\n",
          "CsgBox, CsgRectangularBox, Plane"
          )

newDocumentationItem( "Cone", SF, [VE,VE,FL,FL], 
          "Generates a scalar field of the implicit function of a cone with the tip of the cone\n" +
	  "at the position (1), cone axis (2), cone length (3), and angular width (4).",
          ""
          )

newDocumentationItem( "Plane", SF, [VE,VE], 
          "Generates a scalar field of the implicit function of an infinite plane which has\n" +
	  "the position (1) in the plane, and (2) is the normal to the plane.",
          ""
          )

newDocumentationItem( "Torus", SF, [VE,VE,FL,FL], 
          "Generates a scalar field of the implicit function of a torus centered at the\n" +
	  "position (1) with axis (2), major radius (3) and minor radius (4).",
          ""
          )

newDocumentationItem( "SteinerPatch", SF, [], 
          "Generates a scalar field of the implicit function of a Steiner Patch.",
          ""
          )

newDocumentationItem( "Icosahedron", SF, [], 
          "Generates a scalar field of the implicit function of an icosahedron.",
          ""
          )

newDocumentationItem( "Cylinder", SF, [VE,FL], 
          "Generates a scalar field of the implicit function of an infinitely long cylinder\n" +
	  "centered at the origin, oriented along axis (1) with radius (2).",
          "CappedCylinder"
          )

newDocumentationItem( "CappedCylinder", SF, [VE,VE,FL,FL], 
          "Generates a scalar field of the implicit function of an cylinder with capped ends,\n" +
	  "centered at (1), oriented along axis (2) with radius (4) and length (3).",
          "Cylinder"
          )

newDocumentationItem( "Shell", SF, [SF,FL], 
          "Generates a scalar field of the implicit function of shell around scalar field (1)\n" +
	  "with thickness (2).",
          ""
          )


newDocumentationItem( "mask", SF, [SF], 
          "Generates a scalar field mask of the scalar field (1). The mask evaluates to 1\n" +
	  "if (1) evaluates to a positive value, and 0 if (1) evaluates to 0 or negative.",
          ""
          )

newDocumentationItem( "clamp", SF, [SF,FL,FL], 
          "Generates a scalar field that evaluates to that of (1), but clamped between (2) and (3).\n",
          ""
          )

newDocumentationItem( "pow", SF, [SF,FL], 
          "Generates a scalar field that evaluates to that of (1) raised to the power of (2).\n",
          ""
          )

newDocumentationItem( "pow", SF, [SF,SF], 
          "Generates a scalar field that evaluates to that of (1) raised to the power of\n" +
	  "the evaluation of (2).",
          ""
          )

newDocumentationItem( "pow", CF, [CF,FL], 
          "Generates a color field that evaluates to that of (1) raised to the power of (2)\n" +
	  "on a component-by-component basis.",
          ""
          )

newDocumentationItem( "pow", CF, [CF,SF], 
          "Generates a scalar field that evaluates to that of (1) raised to the power of\n" +
	  "the evaluation of (2) on a component-by-component basis.",
          ""
          )

newDocumentationItem( "BlinnBlend", SF, [SF,SF,FL], 
          "Generates a scalar field that is the implicit function Blinn blend of the two\n" +
	  "scalar fields (1) and (2), with threshold (3), i.e.\n" +
	  "\n" +
	  "    BlinnBlend( s1, s2, a ) = exp(s1) + exp(s2) - a",
          "exp, constant"
          )

newDocumentationItem( "Union", SF, [SF,SF], 
          "Generates a scalar field that is the Constructive Solid Geometry union of\n" +
	  "scalar fields (1) and (2).\n",
          "intersection, cutout"
          )

newDocumentationItem( "intersection", SF, [SF,SF], 
          "Generates a scalar field that is the Constructive Solid Geometry intersection of\n" +
	  "scalar fields (1) and (2).\n",
          "Union, cutout"
          )

newDocumentationItem( "cutout", SF, [SF,SF], 
          "Generates a scalar field that is the Constructive Solid Geometry cutout of\n" +
	  "scalar fields (1) and (2).\n",
          "Union, intersection"
          )

newDocumentationItem( "outer", MF, [VF,VF], 
          "Generates a matrix field that is the outer product of two vector fields (1) and (2).",
          ""
          )

newDocumentationItem( "inverse", MF, [MF], 
          "Generates a matrix field that is the matrix inverse of the matrix field (1).",
          ""
          )

newDocumentationItem( "Pyroclast", SF, [VE,FL,FL,FL,FL,FL,VE,FL,FL], 
          "Generates a scalar field of a sphere displaced pyroclastically. The displacement\n" +
	  "uses fractal summed Perlin noise. The parameters are:\n\n" +
	  "   (1)  Location of the center of the displaced sphere.\n" +
	  "   (2)  Radius of the sphere.\n" +
	  "   (3)  Amplitude of the displacements.\n" +
	  "   (4)  Number of octaves in the fractal sum.\n" +
	  "   (5)  Scaling of the spatial frequency of the noise.\n" +
	  "   (6)  Roughness of the fractal sum.\n" +
	  "   (7)  Translation vector of the Perlin noise.\n" +
	  "   (8)  Time parameter for evolving noise in time.\n" +
	  "   (9)  Power exponent to apply to the displacement noise.",
          "RadialPyroclast"
          )


newDocumentationItem( "RadialPyroclast", SF, [VE,FL,FL,FL,FL,FL,FL,FL,FL], 
          "Generates a scalar field of a sphere displaced pyroclastically. The displacement\n" +
	  "uses fractal summed Perlin noise. The parameters are:\n\n" +
	  "   (1)  Location of the center of the displaced sphere.\n" +
	  "   (2)  Radius of the sphere.\n" +
	  "   (3)  Amplitude of the displacements.\n" +
	  "   (4)  Number of octaves in the fractal sum.\n" +
	  "   (5)  Scaling of the spatial frequency of the noise.\n" +
	  "   (6)  Roughness of the fractal sum.\n" +
	  "   (7)  Radial translation amount of the Perlin noise.\n" +
	  "   (8)  Time parameter for evolving noise in time.\n" +
	  "   (9)  Power exponent to apply to the displacement noise.",
          "Pyroclast"
          )



newDocumentationItem( "SFFFTNoise", SF, [FL,FL,FL,FL,IN], 
          "A scalar field that evaluates to Gaussian random values with spatial correlations.\n" +
	  "Uses Fast Fourier Transforms to create spectral coloration of the random values.\n" +
	  "The parameters are:\n\n" +
	  "   (1)  Power law from scale fall off.\n" +
	  "   (2)  Smallest spatial scale of the noise.\n" +
	  "   (3)  Largest spatial scale of the noise.\n" +
	  "   (4)  Length scale of the peak of the spectrum.\n" +
	  "   (5)  Number of grid points for the fft, applied to all three dimentions.",
          ""
          )




newDocumentationItem( "gridded", SF, [SG], 
          "Generates a scalar field from interpolating gridded data (1).",
          ""
          )

newDocumentationItem( "gridded", VF, [VG], 
          "Generates a vector field from interpolating gridded data (1).",
          ""
          )

newDocumentationItem( "gridded", CF, [CG], 
          "Generates a color field from interpolating gridded data (1).",
          ""
          )

newDocumentationItem( "gridded", MF, [MG], 
          "Generates a matrix field from interpolating gridded data (1).",
          ""
          )

newDocumentationItem( "advect", SF, [SF,VF,FL], 
          "Advect scalar field (1) via Semi-Lagrangian advection using the\n" +
	  "velocity field (2) and time step (3).",
          "warp"
          )
newDocumentationItem( "advect", VF, [VF,VF,FL], 
          "Advect vector field (1) via Semi-Lagrangian advection using the\n" +
	  "velocity field (2) and time step (3).",
          "warp"
          )

newDocumentationItem( "advect", CF, [CF,VF,FL], 
          "Advect color field (1) via Semi-Lagrangian advection using the\n" +
	  "velocity field (2) and time step (3).",
          "warp"
          )

newDocumentationItem( "advect", MF, [MF,VF,FL], 
          "Advect matrix field (1) via Semi-Lagrangian advection using the\n" +
	  "velocity field (2) and time step (3).",
          "warp"
          )

newDocumentationItem( "advect", FF, [FF,VF,FL], 
          "Advect form field (1) via Semi-Lagrangian advection using the\n" +
	  "velocity field (2) and time step (3).",
          "warp"
          )

newDocumentationItem( "warp", SF, [SF,VF], 
          "Remap the scalar field (1) using the mapping field (2).",
          "advect"
          )

newDocumentationItem( "warp", VF, [VF,VF], 
          "Remap the vector field (1) using the mapping field (2).",
          "advect"
          )

newDocumentationItem( "warp", CF, [CF,VF], 
          "Remap the color field (1) using the mapping field (2).",
          "advect"
          )

newDocumentationItem( "warp", MF, [MF,VF], 
          "Remap the matrix field (1) using the mapping field (2).",
          "advect"
          )

newDocumentationItem( "warp", FF, [FF,VF], 
          "Remap the form field (1) using the mapping field (2).",
          "advect"
          )


newDocumentationItem( "Periodic", SF, [SF,VE,VE], 
          "Make scalar field (1) periodic, centered on the point (2), with\n" +
	  "period length (3) in each direction.",
          ""
          )

newDocumentationItem( "Periodic", VF, [VF,VE,VE], 
          "Make vector field (1) periodic, centered on the point (2), with\n" +
	  "period length (3) in each direction.",
          ""
          )

newDocumentationItem( "Periodic", CF, [CF,VE,VE], 
          "Make color field (1) periodic, centered on the point (2), with\n" +
	  "period length (3) in each direction.",
          ""
          )

newDocumentationItem( "Periodic", FF, [FF,VE,VE], 
          "Make form field (1) periodic, centered on the point (2), with\n" +
	  "period length (3) in each direction.",
          ""
          )


newDocumentationItem( "wedge", FF, [FF,FF], 
	  "Wedge (outer) product of two form fields.",
          ""
          )


newDocumentationItem( "start", FF, [FF], 
	  "Hodge star of a form field.",
          ""
          )

newDocumentationItem( "contraction", FF, [VE,FF], 
	  "Contraction of a form field with a vector field.",
          ""
          )

newDocumentationItem( "dot", SF, [VF,VF], 
	  "Inner product of two vector fields.",
          ""
          )

newDocumentationItem( "unitvector", VF, [VF], 
	  "Generates a unit vector field from the vector field (1).",
          ""
          )

newDocumentationItem( "identity", VF, [], 
	  "Generates a vector field that evaluates to the evaluation position.",
          ""
          )

newDocumentationItem( "ImplicitSurfacePoint", VF, [SF,FL,IN], 
	  "Vector field that returns the corresponding position on the implicit surface (1)\n" +
	  "using (3) reversals and maximum step size of (2). ",
          ""
          )

newDocumentationItem( "cross", VF, [VF,VF], 
	  "Generates a vector field that evaluates to the cross product of the input vector fields.",
          ""
          )

newDocumentationItem( "div", SF, [VF], 
	  "Generates a scalar field that evaluates to the divergence of the input vector field.",
          ""
          )

newDocumentationItem( "curl", VF, [VF], 
	  "Generates a vector field that evaluates to the curl of the input vector field.",
          ""
          )

newDocumentationItem( "ContinuedFractionDisplacement", VF, [VF,IN], 
	  "Generates a vector field that evaluates to an iterative approximation of the inverse\n" +
	  "displacement. (1) is the displacement, and (2) is the number of iterations to solve for\n" +
	  "the inverse map. The displacement is frequently used to evaluate a field at a warped location,\n" +
	  "i.e. the warp function using a map X is\n\n" +
	  "   warp( f, X )  =  f( X )\n\n" +
	  "For some applications the inverse Y of the map X is needed, defined as\n\n" +
	  "   Y(X) = identity()\n\n" +
	  "ContinuedFractionDisplacement is an algorithm for constructing Y from the displacement field\n" +
	  "dX = X - x.",
          "warp,advect"
          )

newDocumentationItem( "XYZ", VF, [SF,SF,SF], 
	  "Generates a vector field whose components are the scalar fields (1), (2), and (3).\n" +
	  "Identical to the function \"component\"",
          "component"
          )

newDocumentationItem( "Chroma", CF, [CF], 
	  "Generates a color field consisting of the chroma of the input color field (1).",
          ""
          )

newDocumentationItem( "Blackbody", CF, [SF], 
	  "Generates a color field that returns the Planckian Locus blackbody color in XIE color space.\n" +
	  "The input scalar field (1) serves as the temperature in K.",
          ""
          )

newDocumentationItem( "RGB", CF, [SF,SF,SF], 
	  "Generates a color field whose components are the scalar fields (1), (2), and (3).",
          ""
          )

newDocumentationItem( "DetGrad", SF, [VF], 
	  "Generates a scalar field that evaluates to the determinant of the gradient of the\n" +
	  "vector field (1).",
          ""
          )

newDocumentationItem( "xComponent", SF, [VF], 
	  "Generates a scalar field that evaluates to the x component of the vector field (1).",
          ""
          )

newDocumentationItem( "yComponent", SF, [VF], 
	  "Generates a scalar field that evaluates to the y component of the vector field (1).",
          ""
          )

newDocumentationItem( "zComponent", SF, [VF], 
	  "Generates a scalar field that evaluates to the z component of the vector field (1).",
          ""
          )

newDocumentationItem( "component", VF, [SF,SF,SF], 
	  "Generates a vector field whose components are the scalar fields (1), (2), and (3).\n" +
	  "Identical to the function \"XYZ\"",
          "XYZ"
          )

newDocumentationItem( "zeroComponent", SF, [FF], 
	  "Generates a scalar field that evaluates to the 0-form part of the form field (1).",
          "oneComponent,twoComponent,threeComponent,component"
          )

newDocumentationItem( "oneComponent", VF, [FF], 
	  "Generates a vector field that evaluates to the 1-form part of the form field (1).",
          "zeroComponent,twoComponent,threeComponent,component"
          )

newDocumentationItem( "twoComponent", VF, [FF], 
	  "Generates a vector field that evaluates to the 2-form part of the form field (1).",
          "zeroComponent,oneComponent,threeComponent,component"
          )

newDocumentationItem( "threeComponent", SF, [FF], 
	  "Generates a scalar field that evaluates to the 3-form part of the form field (1).",
          "zeroComponent,oneComponent,twoComponent,component"
          )

newDocumentationItem( "component", FF, [SF,VF,VF,SF], 
	  "Generates a form field using (1) as the 0-form, (2) as the 1-form, (3) as the 2-form,\n" +
	  "and (4) as the 3-form.",
          "zeroComponent,oneComponent,twoComponent,threeComponent"
          )

newDocumentationItem( "lie", FF, [VF,FF], 
	  "Generates a form field of the Lie Derivative of (2) using the vector field (1).",
          ""
          )

newDocumentationItem( "makeGrid",SG, [GB,FL], 
	  "Initializes a scalar grid with dimensions and resolution (1), and default value (2).",
          ""
          )

newDocumentationItem( "makeGrid",VG, [GB,VE], 
	  "Initializes a vector grid with dimensions and resolution (1), and default value (2).",
          ""
          )

newDocumentationItem( "makeGrid",CG, [GB,CO], 
	  "Initializes a color grid with dimensions and resolution (1), and default value (2).",
	  ""
	  )

newDocumentationItem( "makeGrid",MG, [GB,MA], 
	  "Initializes a matrix grid with dimensions and resolution (1), and default value (2).",
	  ""
	  )



newDocumentationItem( "Blur",VO, [SG], 
	  "Simple blur of scalar grid (1) using averaging over all of the nearest neighbors.",
	  ""
	  )

newDocumentationItem( "Blur",VO, [VG], 
	  "Simple blur of vector grid (1) using averaging over all of the nearest neighbors.",
	  ""
	  )

newDocumentationItem( "Blur",VO, [CG], 
	  "Simple blur of color grid (1) using averaging over all of the nearest neighbors.",
	  ""
	  )

newDocumentationItem( "writeGrid",VO, [SG,ST], 
	  "Write a scalar grid (1) into a file with name (2).",
	  ""
	  )

newDocumentationItem( "writeGrid",VO, [VG,ST], 
	  "Write a vector grid (1) into a file with name (2).",
	  ""
	  )

newDocumentationItem( "writeGrid",VO, [CG,ST], 
	  "Write a color grid (1) into a file with name (2).",
	  ""
	  )

newDocumentationItem( "stampBlurredWisps",VO, [SG,VE,FL,VE,VE,FL,IN], 
	  "Stamped an antialiased streak into scalar grid (1), beginning at position (2),\n" +
	  "with a streak length from the velocity (4), acceleration (5), and time (3).\n" +
	  "The scalar value of the wisp is (5), and (6) is a seed for a random number sampling.",
	  ""
	  )

