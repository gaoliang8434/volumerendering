
%module bishop
%{
#include "SparseGrid.h"
%}

%include "SparseGrid.h"

%template(ScalarSGrid)  lux::SGrid<float>;
%template(VectorSGrid)  lux::SGrid<lux::Vector>;
%template(ColorSGrid)  lux::SGrid<lux::Color>;
%template(MatrixSGrid)  lux::SGrid<lux::Matrix>;
%template(ScalarGridBase)  std::shared_ptr<lux::SGrid<float> >;
%template(VectorGridBase)  std::shared_ptr<lux::SGrid<lux::Vector> >;
%template(ColorGridBase)  std::shared_ptr<lux::SGrid<lux::Color> >;

%template(ScalarFSGrid)  lux::FSGrid<float>;
%template(VectorFSGrid)  lux::FSGrid<lux::Vector>;
%template(ColorFSGrid)  lux::FSGrid<lux::Color>;
%template(MatrixFSGrid)  lux::FSGrid<lux::Matrix>;
%template(ScalarFrustumGridBase)  std::shared_ptr<lux::FSGrid<float> >;
%template(VectorFrustumGridBase)  std::shared_ptr<lux::FSGrid<lux::Vector> >;
%template(ColorFrustumGridBase)  std::shared_ptr<lux::FSGrid<lux::Color> >;
