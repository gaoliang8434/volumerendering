
  // Types.h
  // typedef a whole bunch of templated objects from openvdb
  // we want to include the minimal number of header files
  // for compilation speed
 
#ifndef __OPENVDB_TYPEDEF_CODE_H__
#define __OPENVDB_TYPEDEF_CODE_H__
 
#include <openvdb/Grid.h>
#include <openvdb/math/Coord.h>
#include <openvdb/math/Transform.h>
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Vec4.h>
#include <openvdb/tree/Tree.h>
#include <openvdb/util/Util.h>
 
  // typedef a FloatTree and FloatGrid
  typedef openvdb::tree::Tree4<float, 5, 4, 3>::Type FloatTree;
  typedef openvdb::Grid<FloatTree> FloatGrid;
 
  // typedef a Vec3 grid
  typedef openvdb::tree::Tree4<openvdb::math::Vec3d, 5, 4, 3>::Type Vec3dTree;
  typedef openvdb::Grid<Vec3dTree> Vec3dGrid;
 
  // typedef a Vec4 grid (for color and other useful things)
  //typedef openvdb::tree::Tree4<openvdb::math::Vec4f, 5, 4, 3>::Type Vec4fTree;
  typedef openvdb::tree::Tree4<openvdb::Vec4f, 5, 4, 3>::Type Vec4fTree;
  typedef openvdb::Grid<Vec4fTree> Vec4fGrid;

  // typedef a Vec4 grid (for color and other useful things)
  typedef openvdb::tree::Tree4<openvdb::math::Vec4s, 5, 4, 3>::Type Vec4sTree;
  typedef openvdb::Grid<Vec4sTree> Vec4sGrid;
 
  // typedef math objects
  typedef openvdb::math::Coord Coord;
  typedef openvdb::math::Vec3s Vec3s;
  typedef openvdb::math::Vec3s Vec3f;
  typedef openvdb::math::Vec3d Vec3d;
  typedef openvdb::math::Vec3<uint32_t> Vec3I;
  typedef openvdb::math::Vec4<uint32_t> Vec4I;
#define INVALID_IDX openvdb::util::INVALID_IDX
 
  // grid stuff
  typedef openvdb::math::Transform Transform;
  typedef openvdb::math::UniformScaleMap UniformScaleMap;
 
#endif
