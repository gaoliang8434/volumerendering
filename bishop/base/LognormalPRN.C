
#include "LognormalPRN.h"

using namespace lux;

LognormalPRN::LognormalPRN() : mean (1.0), devFactor (-0.5) {}

LognormalPRN::~LognormalPRN(){}

void LognormalPRN::setParameters( const Noise_t& n ) 
{ 
       mean = n.lognormalmean;
       devFactor = -0.5*n.gaussianstandarddeviation*n.gaussianstandarddeviation;
       devScale = n.gaussianstandarddeviation;
       generator.setParameters( n ); 
}

const float LognormalPRN::eval()
{
       float x = devFactor + generator.eval() * devScale;
       return std::exp(x) * mean;
}


