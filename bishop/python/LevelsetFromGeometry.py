# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_LevelsetFromGeometry', [dirname(__file__)])
        except ImportError:
            import _LevelsetFromGeometry
            return _LevelsetFromGeometry
        if fp is not None:
            try:
                _mod = imp.load_module('_LevelsetFromGeometry', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _LevelsetFromGeometry = swig_import_helper()
    del swig_import_helper
else:
    import _LevelsetFromGeometry
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _LevelsetFromGeometry.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _LevelsetFromGeometry.SwigPyIterator_value(self)
    def incr(self, n = 1): return _LevelsetFromGeometry.SwigPyIterator_incr(self, n)
    def decr(self, n = 1): return _LevelsetFromGeometry.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _LevelsetFromGeometry.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _LevelsetFromGeometry.SwigPyIterator_equal(self, *args)
    def copy(self): return _LevelsetFromGeometry.SwigPyIterator_copy(self)
    def next(self): return _LevelsetFromGeometry.SwigPyIterator_next(self)
    def __next__(self): return _LevelsetFromGeometry.SwigPyIterator___next__(self)
    def previous(self): return _LevelsetFromGeometry.SwigPyIterator_previous(self)
    def advance(self, *args): return _LevelsetFromGeometry.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _LevelsetFromGeometry.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _LevelsetFromGeometry.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class Vector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Vector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Vector, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _LevelsetFromGeometry.new_Vector(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_Vector
    __del__ = lambda self : None;
    def set(self, *args): return _LevelsetFromGeometry.Vector_set(self, *args)
    def __add__(self, *args): return _LevelsetFromGeometry.Vector___add__(self, *args)
    def __sub__(self, *args): return _LevelsetFromGeometry.Vector___sub__(self, *args)
    def __div__(self, *args): return _LevelsetFromGeometry.Vector___div__(self, *args)
    def __mul__(self, *args): return _LevelsetFromGeometry.Vector___mul__(self, *args)
    def __xor__(self, *args): return _LevelsetFromGeometry.Vector___xor__(self, *args)
    def __iadd__(self, *args): return _LevelsetFromGeometry.Vector___iadd__(self, *args)
    def __isub__(self, *args): return _LevelsetFromGeometry.Vector___isub__(self, *args)
    def __imul__(self, *args): return _LevelsetFromGeometry.Vector___imul__(self, *args)
    def __idiv__(self, *args): return _LevelsetFromGeometry.Vector___idiv__(self, *args)
    def __call__(self, *args): return _LevelsetFromGeometry.Vector___call__(self, *args)
    def X(self): return _LevelsetFromGeometry.Vector_X(self)
    def Y(self): return _LevelsetFromGeometry.Vector_Y(self)
    def Z(self): return _LevelsetFromGeometry.Vector_Z(self)
    def magnitude(self): return _LevelsetFromGeometry.Vector_magnitude(self)
    def unitvector(self): return _LevelsetFromGeometry.Vector_unitvector(self)
    def normalize(self): return _LevelsetFromGeometry.Vector_normalize(self)
    def __eq__(self, *args): return _LevelsetFromGeometry.Vector___eq__(self, *args)
    def __ne__(self, *args): return _LevelsetFromGeometry.Vector___ne__(self, *args)
    def __lt__(self, *args): return _LevelsetFromGeometry.Vector___lt__(self, *args)
    def __le__(self, *args): return _LevelsetFromGeometry.Vector___le__(self, *args)
    def __gt__(self, *args): return _LevelsetFromGeometry.Vector___gt__(self, *args)
    def __ge__(self, *args): return _LevelsetFromGeometry.Vector___ge__(self, *args)
Vector_swigregister = _LevelsetFromGeometry.Vector_swigregister
Vector_swigregister(Vector)

class Color(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Color, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Color, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _LevelsetFromGeometry.new_Color(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_Color
    __del__ = lambda self : None;
    def set(self, *args): return _LevelsetFromGeometry.Color_set(self, *args)
    def __add__(self, *args): return _LevelsetFromGeometry.Color___add__(self, *args)
    def __sub__(self, *args): return _LevelsetFromGeometry.Color___sub__(self, *args)
    def __div__(self, *args): return _LevelsetFromGeometry.Color___div__(self, *args)
    def __mul__(self, *args): return _LevelsetFromGeometry.Color___mul__(self, *args)
    def __iadd__(self, *args): return _LevelsetFromGeometry.Color___iadd__(self, *args)
    def __isub__(self, *args): return _LevelsetFromGeometry.Color___isub__(self, *args)
    def __imul__(self, *args): return _LevelsetFromGeometry.Color___imul__(self, *args)
    def __idiv__(self, *args): return _LevelsetFromGeometry.Color___idiv__(self, *args)
    def __call__(self, *args): return _LevelsetFromGeometry.Color___call__(self, *args)
    def X(self): return _LevelsetFromGeometry.Color_X(self)
    def Y(self): return _LevelsetFromGeometry.Color_Y(self)
    def Z(self): return _LevelsetFromGeometry.Color_Z(self)
    def W(self): return _LevelsetFromGeometry.Color_W(self)
    def __eq__(self, *args): return _LevelsetFromGeometry.Color___eq__(self, *args)
    def __ne__(self, *args): return _LevelsetFromGeometry.Color___ne__(self, *args)
Color_swigregister = _LevelsetFromGeometry.Color_swigregister
Color_swigregister(Color)

class VolumeBase(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VolumeBase, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VolumeBase, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_VolumeBase()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_VolumeBase
    __del__ = lambda self : None;
    def transform(self, *args): return _LevelsetFromGeometry.VolumeBase_transform(self, *args)
    def addAttribute(self, *args): return _LevelsetFromGeometry.VolumeBase_addAttribute(self, *args)
    def attribute(self, *args): return _LevelsetFromGeometry.VolumeBase_attribute(self, *args)
VolumeBase_swigregister = _LevelsetFromGeometry.VolumeBase_swigregister
VolumeBase_swigregister(VolumeBase)


def SF(*args):
  return _LevelsetFromGeometry.SF(*args)
SF = _LevelsetFromGeometry.SF

def VF(*args):
  return _LevelsetFromGeometry.VF(*args)
VF = _LevelsetFromGeometry.VF

def CF(*args):
  return _LevelsetFromGeometry.CF(*args)
CF = _LevelsetFromGeometry.CF

def WriteFloatVolumeGrid(*args):
  return _LevelsetFromGeometry.WriteFloatVolumeGrid(*args)
WriteFloatVolumeGrid = _LevelsetFromGeometry.WriteFloatVolumeGrid

def ReadFloatVolumeGrid(*args):
  return _LevelsetFromGeometry.ReadFloatVolumeGrid(*args)
ReadFloatVolumeGrid = _LevelsetFromGeometry.ReadFloatVolumeGrid

def WriteVectorVolumeGrid(*args):
  return _LevelsetFromGeometry.WriteVectorVolumeGrid(*args)
WriteVectorVolumeGrid = _LevelsetFromGeometry.WriteVectorVolumeGrid

def ReadVectorVolumeGrid(*args):
  return _LevelsetFromGeometry.ReadVectorVolumeGrid(*args)
ReadVectorVolumeGrid = _LevelsetFromGeometry.ReadVectorVolumeGrid

def WriteColorVolumeGrid(*args):
  return _LevelsetFromGeometry.WriteColorVolumeGrid(*args)
WriteColorVolumeGrid = _LevelsetFromGeometry.WriteColorVolumeGrid

def ReadColorVolumeGrid(*args):
  return _LevelsetFromGeometry.ReadColorVolumeGrid(*args)
ReadColorVolumeGrid = _LevelsetFromGeometry.ReadColorVolumeGrid
SHARED_PTR_DISOWN = _LevelsetFromGeometry.SHARED_PTR_DISOWN
class ScalarVolume(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScalarVolume, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScalarVolume, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ScalarVolume()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ScalarVolume
    __del__ = lambda self : None;
    def eval(self, *args): return _LevelsetFromGeometry.ScalarVolume_eval(self, *args)
    def grad(self, *args): return _LevelsetFromGeometry.ScalarVolume_grad(self, *args)
ScalarVolume_swigregister = _LevelsetFromGeometry.ScalarVolume_swigregister
ScalarVolume_swigregister(ScalarVolume)

class VectorVolume(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VectorVolume, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VectorVolume, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_VectorVolume()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_VectorVolume
    __del__ = lambda self : None;
    def eval(self, *args): return _LevelsetFromGeometry.VectorVolume_eval(self, *args)
    def grad(self, *args): return _LevelsetFromGeometry.VectorVolume_grad(self, *args)
VectorVolume_swigregister = _LevelsetFromGeometry.VectorVolume_swigregister
VectorVolume_swigregister(VectorVolume)

class ColorVolume(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ColorVolume, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ColorVolume, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ColorVolume()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ColorVolume
    __del__ = lambda self : None;
    def eval(self, *args): return _LevelsetFromGeometry.ColorVolume_eval(self, *args)
    def grad(self, *args): return _LevelsetFromGeometry.ColorVolume_grad(self, *args)
ColorVolume_swigregister = _LevelsetFromGeometry.ColorVolume_swigregister
ColorVolume_swigregister(ColorVolume)

class ScalarVolumeGrid(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScalarVolumeGrid, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScalarVolumeGrid, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ScalarVolumeGrid()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ScalarVolumeGrid
    __del__ = lambda self : None;
    def init(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_init(self, *args)
    def setClearValue(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_setClearValue(self, *args)
    def getClearValue(self): return _LevelsetFromGeometry.ScalarVolumeGrid_getClearValue(self)
    def setOutsideValue(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_setOutsideValue(self, *args)
    def getOutsideValue(self): return _LevelsetFromGeometry.ScalarVolumeGrid_getOutsideValue(self)
    def rawPtr(self): return _LevelsetFromGeometry.ScalarVolumeGrid_rawPtr(self)
    def value(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_value(self, *args)
    def set(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_set(self, *args)
    def eval(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_eval(self, *args)
    def normalize(self, *args): return _LevelsetFromGeometry.ScalarVolumeGrid_normalize(self, *args)
ScalarVolumeGrid_swigregister = _LevelsetFromGeometry.ScalarVolumeGrid_swigregister
ScalarVolumeGrid_swigregister(ScalarVolumeGrid)

class VectorVolumeGrid(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VectorVolumeGrid, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VectorVolumeGrid, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_VectorVolumeGrid()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_VectorVolumeGrid
    __del__ = lambda self : None;
    def init(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_init(self, *args)
    def setClearValue(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_setClearValue(self, *args)
    def getClearValue(self): return _LevelsetFromGeometry.VectorVolumeGrid_getClearValue(self)
    def setOutsideValue(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_setOutsideValue(self, *args)
    def getOutsideValue(self): return _LevelsetFromGeometry.VectorVolumeGrid_getOutsideValue(self)
    def rawPtr(self): return _LevelsetFromGeometry.VectorVolumeGrid_rawPtr(self)
    def value(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_value(self, *args)
    def set(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_set(self, *args)
    def eval(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_eval(self, *args)
    def normalize(self, *args): return _LevelsetFromGeometry.VectorVolumeGrid_normalize(self, *args)
VectorVolumeGrid_swigregister = _LevelsetFromGeometry.VectorVolumeGrid_swigregister
VectorVolumeGrid_swigregister(VectorVolumeGrid)

class ColorVolumeGrid(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ColorVolumeGrid, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ColorVolumeGrid, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ColorVolumeGrid()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ColorVolumeGrid
    __del__ = lambda self : None;
    def init(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_init(self, *args)
    def setClearValue(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_setClearValue(self, *args)
    def getClearValue(self): return _LevelsetFromGeometry.ColorVolumeGrid_getClearValue(self)
    def setOutsideValue(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_setOutsideValue(self, *args)
    def getOutsideValue(self): return _LevelsetFromGeometry.ColorVolumeGrid_getOutsideValue(self)
    def rawPtr(self): return _LevelsetFromGeometry.ColorVolumeGrid_rawPtr(self)
    def value(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_value(self, *args)
    def set(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_set(self, *args)
    def eval(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_eval(self, *args)
    def normalize(self, *args): return _LevelsetFromGeometry.ColorVolumeGrid_normalize(self, *args)
ColorVolumeGrid_swigregister = _LevelsetFromGeometry.ColorVolumeGrid_swigregister
ColorVolumeGrid_swigregister(ColorVolumeGrid)

class ScalarField(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ScalarField, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ScalarField, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ScalarField()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ScalarField
    __del__ = lambda self : None;
ScalarField_swigregister = _LevelsetFromGeometry.ScalarField_swigregister
ScalarField_swigregister(ScalarField)

class VectorField(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, VectorField, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, VectorField, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_VectorField()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_VectorField
    __del__ = lambda self : None;
VectorField_swigregister = _LevelsetFromGeometry.VectorField_swigregister
VectorField_swigregister(VectorField)

class ColorField(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ColorField, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ColorField, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _LevelsetFromGeometry.new_ColorField()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_ColorField
    __del__ = lambda self : None;
ColorField_swigregister = _LevelsetFromGeometry.ColorField_swigregister
ColorField_swigregister(ColorField)

class LevelsetFromGeometry(ScalarVolume):
    __swig_setmethods__ = {}
    for _s in [ScalarVolume]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, LevelsetFromGeometry, name, value)
    __swig_getmethods__ = {}
    for _s in [ScalarVolume]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, LevelsetFromGeometry, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _LevelsetFromGeometry.new_LevelsetFromGeometry(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _LevelsetFromGeometry.delete_LevelsetFromGeometry
    __del__ = lambda self : None;
    def getGrid(self): return _LevelsetFromGeometry.LevelsetFromGeometry_getGrid(self)
    def eval(self, *args): return _LevelsetFromGeometry.LevelsetFromGeometry_eval(self, *args)
LevelsetFromGeometry_swigregister = _LevelsetFromGeometry.LevelsetFromGeometry_swigregister
LevelsetFromGeometry_swigregister(LevelsetFromGeometry)


def geom2ls(*args):
  return _LevelsetFromGeometry.geom2ls(*args)
geom2ls = _LevelsetFromGeometry.geom2ls
# This file is compatible with both classic and new-style classes.

