

#ifndef __VSHADER_H__
#define __VSHADER_H__


namespace lux
{

class RenderData;

struct ShaderData
{
   Vector *P;
   Vector *D;
   float *density;
   Color *litColor;
   Color *ambientColor;

   const RenderData* renderData;

   Color *Cf;
   float *T;

};

class VShader
{
  public:

    VShader(){}
    virtual ~VShader(){}

    virtual void eval( ShaderData& data ){}

};






}






#endif
