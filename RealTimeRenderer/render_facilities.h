//
// Created by Kobe De Broeck on 11-3-2022.
//

#ifndef ENGINE_RENDER_FACILITIES_H
#define ENGINE_RENDER_FACILITIES_H


#include "render_facilities.h"
#include "easy_image.h"
#include "vector3d.h"
#include <math.h>

#include <limits>

namespace KDBRenderFacilities {

    /*
     *      NormColor
     */

    /// Color with normalized range(0,1) rgb values
    struct NormColor {
        double r,g,b;
    };

    constexpr NormColor magenta {1,0,1};
    constexpr NormColor red {1,0,0};
    constexpr NormColor green {0,1,0};
    constexpr NormColor blue {0,0,1};

    img::Color toEasyImgColor(NormColor c);
    NormColor toNormColor(const std::vector<double>& vec);

    /*
     *      Point2D
     */

    struct Point2D {
        double x, y;
    };

    constexpr Point2D operator+(Point2D p0, Point2D p1)
    {
        return {p0.x + p1.x, p0.y + p1.y};
    }

    constexpr Point2D operator*(Point2D p, double d)
    {
        return {p.x * d, p.y * d};
    }

    constexpr Point2D operator*(double d, Point2D p)
    {
        return p * d;
    }

    struct ZBuffer {
    private:
        std::vector<double> buffer; // values stored in row major
        unsigned width, height;
    public:
        ZBuffer(unsigned width, unsigned height)
                : buffer(width*height, std::numeric_limits<double>::infinity()), width{width}, height{height} {}

        double& operator()(unsigned row, unsigned col) {return buffer[row * width + col];}
        double operator()(unsigned row, unsigned col) const {return buffer[row * width + col];}

        int getWidth() const {return width;}
        int getHeight() const {return height;}
    };

    /*
     *      Line2D
     */

    struct Line2D {
        Point2D p0, p1;
        NormColor c;

        double z1;
        double z2;

        Line2D(Point2D p0, Point2D p1, NormColor c, double z1, double z2) : p0{p0}, p1{p1}, c{c}, z1{z1}, z2{z2} {}
    };

    img::EasyImage draw2DLines(const std::vector<Line2D>& lines, unsigned size, NormColor backgroundcolor);
    img::EasyImage draw2DLines(ZBuffer& zBuffer, const std::vector<Line2D>& lines, unsigned size, NormColor backgroundcolor);

    /*
     *      3D stuff
     */

    struct PolarCoord {
        double theta{0};
        double phi{0};
        double r{0};
    };

    inline PolarCoord toPolar(const Vector3D& point)
    {
        PolarCoord p;
        p.r = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
        p.theta = atan2(point.y, point.x);
        p.phi = acos(point.z/p.r);
        return p;
    }

    struct Face {
        std::vector<int> point_indexes{};
        Face() = default;
        //Face(std::initializer_list<int> point_indexes) : point_indexes{point_indexes} {}
    };

    struct Figure {
        std::vector<Vector3D> vertices{};
        std::vector<Face> faces{};
        NormColor color{};
        Figure() = default;
        Figure(const std::vector<Vector3D>& vertices, std::vector<Face>& faces, NormColor c)
            : vertices{vertices}, faces{faces}, color{c} {}

        void triangulate();
    };

    inline Vector3D midPoint(const Vector3D& p0, const Vector3D& p1)
    {
        return Vector3D::point((p0.x+p1.x)/2, (p0.y+p1.y)/2, (p0.z+p1.z)/2);
    }

    inline Vector3D toVector3D(const std::vector<double>& v)
    {
        if(v.size() != 3) std::cerr << "ERROR: toVector3D: input vector is not of size 3\n";
        return Vector3D::point(v[0],v[1],v[2]);
    }


    Matrix scaling(double scale);

    Matrix rotationX(double angle);

    Matrix rotationY(double angle);

    Matrix rotationZ(double angle);

    Matrix translation(const Vector3D& vector);

    inline Matrix eyePointTransformation(PolarCoord p)
    {
        Matrix m;
        m(1,1) = -sin(p.theta);
        m(1,2) = -cos(p.theta)*cos(p.phi);
        m(1,3) = cos(p.theta)*sin(p.phi);
        m(2,1) = cos(p.theta);
        m(2,2) = -sin(p.theta)*cos(p.phi);
        m(2,3) = sin(p.theta)*sin(p.phi);
        m(3,2) = sin(p.phi);
        m(3,3) = cos(p.phi);
        m(4,3) = -p.r;
        m(4,4) = 1;
        return m;
    }

    void applyTransformation(Figure& fig, const Matrix& m);

    std::vector<Line2D> doProjection(const std::vector<Figure>& figures);

    std::vector<Line2D> toPerspectiveProjectedLines(const std::vector<Figure>& figures, PolarCoord eye);


    /*
     *  Angle conversions
     */

    constexpr double KDB_PI = 3.141592653;
    constexpr double rad(double angleDeg) { return angleDeg * KDB_PI / 180; }


}





#endif //ENGINE_RENDER_FACILITIES_H
