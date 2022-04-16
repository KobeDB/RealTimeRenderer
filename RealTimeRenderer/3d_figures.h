//
// Created by Kobe De Broeck on 3-4-2022.
//

#ifndef ENGINE_3D_FIGURES_H
#define ENGINE_3D_FIGURES_H

#include "render_facilities.h"

#include "ini_configuration.h"


namespace KDBRenderFacilities {

    img::EasyImage drawWireFrame(const ini::Configuration& configuration);
    img::EasyImage drawZBufferedWireFrame(const ini::Configuration& configuration);

    img::EasyImage drawWireFrame(const std::vector<Figure>& figures, const Vector3D& eye, int size, const NormColor& backgroundColor);

}

namespace KDBRenderFacilities {
    Figure createCube(const NormColor& color);
    Figure createTetrahedron(const NormColor& color);
    Figure createOctahedron(const NormColor& color);
    Figure createIcosahedron(const NormColor& color);
    Figure createDodecahedron(const NormColor& color);
    Figure createSphere(double radius, int n, const NormColor& color);
    Figure createCone(int n, double h, const NormColor& color);
    Figure createCylinder(int n, double h, const NormColor& c);
    Figure createTorus(double r, double R, int n, int m, const NormColor& color);

    Figure generate3DLSystemFigure(const std::string& inputFile);
}

#endif //ENGINE_3D_FIGURES_H
