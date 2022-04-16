//
// Created by Kobe De Broeck on 11-3-2022.
//

#include "render_facilities.h"

#include <cmath>
#include <vector>
#include <limits>
#include "easy_image.h"

namespace KDBRenderFacilities {

    using namespace std;

    img::Color toEasyImgColor(NormColor c)
    {
        unsigned char red = c.r * 255;
        unsigned char green = c.g * 255;
        unsigned char blue = c.b * 255;
        return img::Color{red, green, blue};
    }

    NormColor toNormColor(const std::vector<double>& vec)
    {
        if(vec.size() != 3) throw std::runtime_error("ERROR toNormColor: Tried to convert vector of wrong size to color");
        return {vec[0], vec[1], vec[2]};
    }

    void updateIfBigger(double maybeBigger, double& max )
    {
        if(maybeBigger > max) {
            max = maybeBigger;
        }
    }

    void updateIfSmaller(double maybeSmaller, double& min)
    {
        if(maybeSmaller < min) {
            min = maybeSmaller;
        }
    }

    void Figure::triangulate() {
        vector<Face> triangle_faces;
        triangle_faces.reserve(faces.size());
        for(const Face& f : faces) {
            if(f.point_indexes.size() <= 3) {
                triangle_faces.push_back(f);
                continue;
            }
            for(vector<Face>::size_type i = 1; i < f.point_indexes.size()-1; i++) {
                triangle_faces.push_back(Face{{f.point_indexes[0], f.point_indexes[i], f.point_indexes[i+1]}});
            }
        }
        faces = triangle_faces;
    }

    void draw2DLine(img::EasyImage& image, Point2D left, Point2D right, img::Color c)
    {
        if( false )
        {
            cerr << "ERROR: Tried to draw a line whose point(s) lay out of bounds!\n";
            return;
        }

        if(left.x > right.x) {
            auto temp = left;
            left = right;
            right = temp;
        }

        int startX = ceil(left.x);
        int startY = ceil(left.y);
        int endX = ceil(right.x);
        int endY = ceil(right.y);

        // !!!! Test this before calculating m to prevent division by 0 !!!!
        if ( startX == endX ) {
            int StartDrawY = left.y < right.y ? ceil(left.y) : ceil(right.y);
            auto amtStepsY =  abs(startY - endY);
            for ( int i = 0; i < amtStepsY; i++ ) {
                image(startX, StartDrawY + i) = c;
            }
        }

        double m = (double)(endY - startY) / (endX - startX); // cast to double needed to perform floating point division

        if( startY == endY) {
            for( int i = 0; i < endX - startX; i++) {
                image(startX + i, startY) = c;
            }
        }

        if( (0 < m && m <= 1) || (-1 <= m && m < 0) ) {
            for( int i = 0; i < endX - startX; i++) {
                image(startX + i, round(startY + m * i)) = c;
            }
        }

        if(m > 1) {
            for( int i = 0; i < endY - startY; i++) {
                image(round(startX + i/m), startY + i) = c;
            }
        }

        if( m < -1) {
            for( int i = 0; i < startY - endY; i++) {
                image(round(startX - i/m), startY - i) = c;
            }
        }
    }

    void draw2DZBufLine(ZBuffer& zBuffer, img::EasyImage& image, Point2D point0, double z0, Point2D point1, double z1, img::Color c)
    {
        Point2D left = point0;
        Point2D right = point1;
        if(left.x > right.x) {
//            auto temp = left;
//            left = right;
//            right = temp;
            std::swap(left, right);
            std::swap(z0,z1); // Very Important, as z0 belongs to point0, and z1 belongs to point1!! (I think)
        }

        int startX = ceil(left.x);
        int startY = ceil(left.y);
        int endX = ceil(right.x);
        int endY = ceil(right.y);

        int dx = endX - startX;
        int dy = endY - startY;

        // !!!! Test this before calculating m to prevent division by 0 !!!!
        if ( startX == endX ) {
            int StartDrawY = left.y < right.y ? ceil(left.y) : ceil(right.y);
            auto amtStepsY =  abs(startY - endY);
            for ( int i = 0; i < amtStepsY; i++ ) {
                double p = (double)(dy - i) / dx;
                double reciprocal_z_I = p / z0 + ( 1 - p )/z1;
                unsigned pix_X = startX;
                unsigned pix_Y = StartDrawY + i;
                if(reciprocal_z_I < zBuffer(pix_X, pix_Y)) {
                    image(pix_X, pix_Y) = c;
                    zBuffer(pix_X, pix_Y) = reciprocal_z_I;
                }
            }
        }

        double m = (double)(endY - startY) / (endX - startX); // cast to double needed to perform floating point division

        if( startY == endY) {
            for( int i = 0; i < endX - startX; i++) {
                double p = (double)(dx - i) / dx;
                double reciprocal_z_I = p / z0 + ( 1 - p )/z1;
                unsigned pix_X = startX + i;
                unsigned pix_Y = startY;
                if(reciprocal_z_I < zBuffer(pix_X, pix_Y)) {
                    image(pix_X, pix_Y) = c;
                    zBuffer(pix_X, pix_Y) = reciprocal_z_I;
                }
            }
        }

        if( (0 < m && m <= 1) || (-1 <= m && m < 0) ) {
            for( int i = 0; i < endX - startX; i++)
            {
                double p = (double)(dx - i) / dx;
                double reciprocal_z_I = p / z0 + ( 1 - p )/z1;
                unsigned pix_X = startX + i;
                unsigned pix_Y = round(startY + m * i);
                if(reciprocal_z_I < zBuffer(pix_X, pix_Y)) {
                    image(pix_X, pix_Y) = c;
                    zBuffer(pix_X, pix_Y) = reciprocal_z_I;
                }
            }
        }

        if(m > 1) {
            for( int i = 0; i < endY - startY; i++) {
                double p = (double)(dy - i) / dy;
                double reciprocal_z_I = p / z0 + ( 1 - p )/z1;
                unsigned pix_X = round(startX + i/m);
                unsigned pix_Y = startY + i;
                if(reciprocal_z_I < zBuffer(pix_X, pix_Y)) {
                    image(pix_X, pix_Y) = c;
                    zBuffer(pix_X, pix_Y) = reciprocal_z_I;
                }
            }
        }

        if( m < -1) {
            for( int i = 0; i < startY - endY; i++) {
                double p = (double)(dy - i) / dy;
                double reciprocal_z_I = p / z0 + ( 1 - p )/z1;
                unsigned pix_X = round(startX - i/m);
                unsigned pix_Y = startY - i;
                if(reciprocal_z_I < zBuffer(pix_X, pix_Y)) {
                    image(pix_X, pix_Y) = c;
                    zBuffer(pix_X, pix_Y) = reciprocal_z_I;
                }
            }
        }
    }

    void calcPointTransformationsFor2DLines(const std::vector<Line2D>& lines, unsigned size,
                                            double& d, double& dx, double& dy, unsigned& imageX, unsigned& imageY)
    {
        constexpr double infinity{std::numeric_limits<int>::max()};
        double xmax{-infinity}, xmin{infinity}, ymax{-infinity}, ymin{infinity};
        for(const auto& line : lines) {
            updateIfBigger(line.p0.x, xmax);
            updateIfBigger(line.p1.x, xmax);
            updateIfBigger(line.p0.y, ymax);
            updateIfBigger(line.p1.y, ymax);
            updateIfSmaller(line.p0.x, xmin);
            updateIfSmaller(line.p1.x, xmin);
            updateIfSmaller(line.p0.y, ymin);
            updateIfSmaller(line.p1.y, ymin);
        }

        double xrange = xmax - xmin;
        double yrange = ymax - ymin;

        imageX = size * xrange / std::max(xrange, yrange);
        imageY = size * yrange / std::max(xrange, yrange);

        d = 0.95 * imageX/xrange;

        double DCx = d * (xmax + xmin) / 2;
        double DCy = d * (ymax + ymin) / 2;
        dx = ((double)imageX/2) - DCx;
        dy = ((double)imageY/2) - DCy;
    }

    img::EasyImage draw2DLines(const std::vector<Line2D>& lines, unsigned size, NormColor backgroundcolor)
    {
        double d, dx, dy;
        unsigned imageX, imageY;
        calcPointTransformationsFor2DLines(lines, size, d, dx, dy, imageX, imageY);
        img::EasyImage image{imageX, imageY, toEasyImgColor(backgroundcolor)};
        for(const auto& line : lines) {
            Point2D first = line.p0 * d + Point2D{dx, dy};
            Point2D second = line.p1 * d + Point2D{dx, dy};
            draw2DLine(image, first, second, toEasyImgColor(line.c));
            //draw2DZBufLine(zBuffer, image, first, line.z1, second, line.z2, toEasyImgColor(line.c));
        }

        return image;
    }

    img::EasyImage draw2DLines(ZBuffer& zBuffer, const std::vector<Line2D>& lines, unsigned size, NormColor backgroundcolor)
    {
        double d, dx, dy;
        unsigned imageX, imageY;
        calcPointTransformationsFor2DLines(lines, size, d, dx, dy, imageX, imageY);
        img::EasyImage image{imageX, imageY, toEasyImgColor(backgroundcolor)};
        for(const auto& line : lines) {
            Point2D first = line.p0 * d + Point2D{dx, dy};
            Point2D second = line.p1 * d + Point2D{dx, dy};
            draw2DZBufLine(zBuffer, image, first, line.z1, second, line.z2, toEasyImgColor(line.c));
        }

        return image;
    }

    Matrix scaling(double scale)
    {
        Matrix m;
        m(1,1) = scale;
        m(2,2) = scale;
        m(3,3) = scale;
        return m;
    }

    Matrix rotationX(double angle)
    {
        angle = rad(angle);
        Matrix m;
        m(2,2) = cos(angle);
        m(3,2) = -sin(angle);
        m(2,3) = sin(angle);
        m(3,3) = cos(angle);
        return m;
    }

    Matrix rotationY(double angle)
    {
        angle = rad(angle);

        Matrix m;
        m(1,1) = cos(angle);
        m(1,3) = -sin(angle);
        m(3,1) = sin(angle);
        m(3,3) = cos(angle);
        return m;
    }

    Matrix rotationZ(double angle)
    {
        angle = rad(angle);

        Matrix m;
        m(1,1) = cos(angle);
        m(1,2) = sin(angle);
        m(2,1) = -sin(angle);
        m(2,2) = cos(angle);
        return m;
    }

    Matrix translation(const Vector3D& vector)
    {
        Matrix m;
        m(4,1) = vector.x;
        m(4,2) = vector.y;
        m(4,3) = vector.z;
        return m;
    }

    void applyTransformation(Figure& fig, const Matrix& m)
    {
        for(Vector3D& vert : fig.vertices) {
            vert = vert * m;
        }
    }

    void doProjection(Figure& fig, double d = 1)
    {
        for(auto& v : fig.vertices) {
            v.x = d * v.x / -v.z;
            v.y = d * v.y / -v.z;
        }
    }

    std::vector<Line2D> toPerspectiveProjectedLines(const std::vector<Figure>& pFigures, PolarCoord eye)
    {
        vector<Figure> figures = pFigures; // make a local copy of pFigures

        vector<Line2D> lines;

        for(auto& fig : figures) {
            applyTransformation(fig, eyePointTransformation(eye));
            doProjection(fig);
            for(const auto& face : fig.faces) {
                for(int i = 0; i < face.point_indexes.size(); i++) {
                    int p0_i = face.point_indexes[i];
                    int p1_i = face.point_indexes[(i+1)%face.point_indexes.size()]; //This loops around to 0, to draw a line between the first and the last vertex
                    const Vector3D& p0 = fig.vertices[p0_i];
                    const Vector3D& p1 = fig.vertices[p1_i];
                    Point2D projP0 = {p0.x, p0.y};
                    Point2D projP1 = {p1.x, p1.y};
                    lines.push_back(Line2D{projP0, projP1, fig.color, p0.z, p1.z});
                }
            }
        }
        return lines;
    }

}
