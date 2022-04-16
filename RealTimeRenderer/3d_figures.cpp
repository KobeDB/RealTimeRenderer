//
// Created by Kobe De Broeck on 3-4-2022.
//

#include "3d_figures.h"
#include "l_parser.h"

#include <stack>
#include <fstream>
#include <vector>

using namespace std;

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

namespace KDBRenderFacilities {

    // Unused
    img::EasyImage drawWireFrame(const vector<Figure>& figures, const Vector3D& eye, int size, const NormColor& backgroundColor)
    {
        PolarCoord polarEye = toPolar(Vector3D::point(eye.x,eye.y,eye.z));
        auto lines = toPerspectiveProjectedLines(figures, polarEye);

        return draw2DLines(lines, size, backgroundColor);
    }

    struct WireFrameImageInfo{
        unsigned imageSize{};
        vector<Line2D> lines{};
        NormColor backgroundColor{};
    };

    WireFrameImageInfo getWireFrameImageInfo(const ini::Configuration& configuration)
    {
        unsigned size = configuration["General"]["size"].as_int_or_default(0);
        auto backgroundColorVec = configuration["General"]["backgroundcolor"].as_double_tuple_or_default({0,0,0});
        unsigned nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
        auto eye = configuration["General"]["eye"].as_double_tuple_or_default({0,0,0});

        vector<Figure> figures;

        for(int i = 0; i < nrFigures; i++) {
            string figName = "Figure" + to_string(i);
            const auto &figConf = configuration[figName];

            auto c = figConf["color"];
            NormColor figureColor = toNormColor(c);

            Figure fig {};

            string type = figConf["type"];

            if (type == "LineDrawing") {

                int nrPoints = figConf["nrPoints"];
                int nrLines = figConf["nrLines"];

                vector<Vector3D> vertices;
                vector<Face> faces;

                for (int j = 0; j < nrPoints; j++) {

                    auto point = figConf["point" + to_string(j)].as_double_tuple_or_die();
                    Vector3D v = toVector3D(point);
                    vertices.push_back(v);
                }

                for (int j = 0; j < nrLines; j++) {
                    auto line = figConf["line" + to_string(j)].as_int_tuple_or_die();
                    if (line.size() != 2) cerr << "ERROR: line vector not of size 2\n";
                    faces.push_back({{line[0], line[1]}});
                }

                fig = {vertices, faces, figureColor};
            }

            if(type == "Cube") {
                fig = createCube(figureColor);
            }

            if(type == "Tetrahedron") {
                fig = createTetrahedron(figureColor);
            }

            if(type == "Octahedron") {
                fig = createOctahedron(figureColor);
            }

            if(type == "Icosahedron") {
                fig = createIcosahedron(figureColor);
            }

            if(type == "Dodecahedron") {
                fig = createDodecahedron(figureColor);
            }

            if(type == "Sphere") {
                int n = figConf["n"].as_int_or_default(0);
                fig = createSphere(1,n,figureColor);
            }

            if(type == "Cone") {
                int n = figConf["n"].as_int_or_default(3);
                double height = figConf["height"].as_double_or_default(1);
                fig = createCone(n, height, figureColor);
            }

            if(type == "Cylinder") {
                int n = figConf["n"].as_int_or_default(3);
                double height = figConf["height"].as_double_or_default(1);
                fig = createCylinder(n, height, figureColor);
            }

            if(type == "Torus") {
                double r = figConf["r"].as_double_or_default(1);
                double R = figConf["R"].as_double_or_default(1);
                int n = figConf["n"].as_int_or_default(3);
                int m = figConf["m"].as_int_or_default(3);
                fig = createTorus(r, R, n, m, figureColor);
            }

            if(type == "3DLSystem") {
                string inputFile = figConf["inputfile"].as_string_or_die();
                fig = generate3DLSystemFigure(inputFile);
                fig.color = figureColor;
            }

            fig.triangulate();

            double scale = figConf["scale"];
            applyTransformation(fig, scaling(scale));

            double rotateX = figConf["rotateX"];
            applyTransformation(fig, rotationX(rotateX));
            double rotateY = figConf["rotateY"];
            applyTransformation(fig, rotationY(rotateY));
            double rotateZ = figConf["rotateZ"];
            applyTransformation(fig, rotationZ(rotateZ));

            const auto& figCenter = figConf["center"].as_double_tuple_or_default({0,0,0});
            applyTransformation(fig, translation(toVector3D(figCenter)));

            figures.push_back(fig);

        }
        PolarCoord polarEye = toPolar(Vector3D::point(eye[0],eye[1],eye[2]));
        auto lines = toPerspectiveProjectedLines(figures, polarEye);
        return {size, lines, toNormColor(backgroundColorVec)};
    }

    img::EasyImage drawWireFrame(const ini::Configuration& configuration)
    {
        if(configuration["General"]["type"].as_string_or_default("") != "Wireframe") {
            cerr << "drawWireFrame didn't get a Wireframe config!\n";
            return {};
        }
        WireFrameImageInfo info =  getWireFrameImageInfo(configuration);
        return draw2DLines(info.lines, info.imageSize, info.backgroundColor);
    }

    img::EasyImage drawZBufferedWireFrame(const ini::Configuration& configuration)
    {
        if(configuration["General"]["type"].as_string_or_default("") != "ZBufferedWireframe") {
            cerr << "drawWireFrame didn't get a ZBufferedWireframe config!\n";
            return {};
        }
        WireFrameImageInfo info =  getWireFrameImageInfo(configuration);
        ZBuffer zBuffer(info.imageSize, info.imageSize);
        return draw2DLines(zBuffer, info.lines, info.imageSize, info.backgroundColor);
    }

    struct PenState {
        Vector3D pos;

        Vector3D H; // look at
        Vector3D L; // left
        Vector3D U; // up

        double angleDelta{};

        PenState() : pos{Vector3D::point(0,0,0)}, H{Vector3D::point(1,0,0)}, L{Vector3D::point(0,1,0)}, U{Vector3D::point(0,0,1)} {}
    };



    void generate3DLSystemFigureRecursive(
            int nrIterations,
            const std::string& penInstructions,
            Figure& fig,
            const LParser::LSystem3D& lSystem,
            PenState& curPenState,
            stack<PenState>& penStateStack
            )
    {
        for(char c : penInstructions) {
            if(lSystem.get_alphabet().find(c) != lSystem.get_alphabet().end()) {
                if(nrIterations <= 0) {
                    if(lSystem.draw(c)) {
                        auto next = curPenState.pos + curPenState.H;
                        fig.vertices.push_back(curPenState.pos); int prevPos_i = fig.vertices.size()-1;
                        fig.vertices.push_back(next); int next_i = fig.vertices.size()-1;
                        fig.faces.push_back(Face{{prevPos_i, next_i}});

                        curPenState.pos = next;
                    }
                    continue;
                }
                generate3DLSystemFigureRecursive(nrIterations-1, lSystem.get_replacement(c), fig, lSystem, curPenState, penStateStack);
                continue;
            }
            if( c == '(') {
                penStateStack.push(curPenState);
                continue;
            }
            if(c == ')') {
                curPenState = penStateStack.top();
                penStateStack.pop();
                continue;
            }
            double angle = curPenState.angleDelta;
            if(c == '-' || c == '&' || c == '/') {
                angle = -angle;
            }
            if(c == '+' || c == '-') {
                auto H_new = curPenState.H * cos(angle) + curPenState.L * sin(angle);
                auto L_new = -curPenState.H * sin(angle) + curPenState.L * cos(angle);

                curPenState.H = H_new;
                curPenState.L = L_new;
                continue;
            }
            if(c == '^' || c == '&') {
                auto H_new = (curPenState.H * cos(angle)) + (curPenState.U * sin(angle));
                auto U_new = (-curPenState.H * sin(angle)) + (curPenState.U * cos(angle));

                curPenState.H = H_new;
                curPenState.U = U_new;
                continue;
            }
            if(c == '\\' || c == '/') {
                auto L_new = curPenState.L * cos(angle) - curPenState.U * sin(angle);
                auto U_new = curPenState.L * sin(angle) + curPenState.U * cos(angle);

                curPenState.L = L_new;
                curPenState.U = U_new;
                continue;
            }
            if(c == '|') {
                angle = KDB_PI;
                curPenState.H = -curPenState.H;
                curPenState.L = -curPenState.L;
                continue;
            }
            cerr << "Unhandled symbol in 3d lsystem\n";
        }
    }

    Figure generate3DLSystemFigure(const std::string& inputFile)
    {
        std::stack<PenState> penStateStack{};
        PenState curPenState{};

        LParser::LSystem3D lSystem;
        std::ifstream ifs{inputFile};
        ifs >> lSystem;
        if(!ifs) std::cerr << "Couldn't read l-system\n";

        curPenState.angleDelta = rad(lSystem.get_angle());

        Figure fig;
        generate3DLSystemFigureRecursive(lSystem.get_nr_iterations(), lSystem.get_initiator(), fig, lSystem, curPenState, penStateStack);

        return fig;
    }

    /// in: face point-indexes that start from 1 --> returns: indexes that start from 0
    vector<Face> toFacesStartingWithIndex0(const vector<Face>& in)
    {
        vector<Face> result;
        for(const auto& f : in ) {
            vector<int> newIndexes;
            for(auto index : f.point_indexes) {
                newIndexes.push_back(index-1);
            }
            result.push_back(Face{newIndexes});
        }
        return result;
    }

    Figure createSphereRecursive(int n, const Figure& sphereToRefine)
    {
        if(n <= 0) return sphereToRefine;

        vector<Vector3D> vertices;
        vector<Face> faces;
        for(const auto& f : sphereToRefine.faces) {
            if(f.point_indexes.size() != 3) cerr << "intermediate sphere's face doesn't have exactly 3 points!\n";
            const Vector3D& p0 = sphereToRefine.vertices[f.point_indexes[0]];
            const Vector3D& p1 = sphereToRefine.vertices[f.point_indexes[1]];
            const Vector3D& p2 = sphereToRefine.vertices[f.point_indexes[2]];

            Vector3D p0p1 = midPoint(p0,p1);
            Vector3D p1p2 = midPoint(p1,p2);
            Vector3D p2p0 = midPoint(p2,p0);

            vertices.push_back(p0);     int p0_i = vertices.size()-1; // Store the index of each point in the vertices vector
            vertices.push_back(p1);     int p1_i = vertices.size()-1;
            vertices.push_back(p2);     int p2_i = vertices.size()-1;
            vertices.push_back(p0p1);   int p0p1_i = vertices.size()-1;
            vertices.push_back(p1p2);   int p1p2_i = vertices.size()-1;
            vertices.push_back(p2p0);   int p2p0_i = vertices.size()-1;

            faces.push_back(Face{{p0_i, p0p1_i, p2p0_i}});
            faces.push_back(Face{{p0p1_i, p1_i, p1p2_i}});
            faces.push_back(Face{{p2p0_i, p1p2_i, p2_i}});
            faces.push_back(Face{{p2p0_i, p0p1_i, p1p2_i}}); // Middle triangle also has its own face apparently
        }

        return createSphereRecursive(n-1, Figure{vertices, faces, magenta});
    }

    Figure createSphere(double radius, int n, const NormColor& color)
    {
        Figure icosahedron = createIcosahedron(color);

        Figure sphere = createSphereRecursive(n, icosahedron);

        for(auto& p : sphere.vertices) {
            double r = p.length();
            p = Vector3D::point(p.x/r, p.y/r, p.z/r);
        }

        sphere.color = color;

        return sphere;
    }

    Figure createDodecahedron(const NormColor& color)
    {
        Figure icosahedron = createIcosahedron(color);

        vector<Vector3D> dodecaVertices;
        for(const auto& f : icosahedron.faces) {
            if(f.point_indexes.size() != 3) cerr << "icosahedron's face doesn't have exactly 3 points!\n";
            const Vector3D& p0 = icosahedron.vertices[f.point_indexes[0]];
            const Vector3D& p1 = icosahedron.vertices[f.point_indexes[1]];
            const Vector3D& p2 = icosahedron.vertices[f.point_indexes[2]];

            double dodecaX = (p0.x + p1.x + p2.x) / 3;
            double dodecaY = (p0.y + p1.y + p2.y) / 3;
            double dodecaZ = (p0.z + p1.z + p2.z) / 3;

            dodecaVertices.push_back(Vector3D::point(dodecaX, dodecaY, dodecaZ));
        }

        vector<Vector3D> vertices = dodecaVertices;
        vector<Face> faces = toFacesStartingWithIndex0({
                                                               Face{{1,2,3,4,5}},
                                                               {{1,6,7,8,2}},
                                                               {{2,8,9,10,3}},
                                                               {{3,10,11,12,4}},
                                                               {{4,12,13,14,5}},
                                                               {{5,14,15,6,1}},
                                                               {{20,19,18,17,16}},
                                                               {{20,15,14,13,19}},
                                                               {{19,13,12,11,18}},
                                                               {{18,11,10,9,17}},
                                                               {{17,9,8,7,16}},
                                                               {{16,7,6,15,20}}
                                                       });

        return {vertices, faces, color};
    }

    Figure createIcosahedron(const NormColor& color)
    {
        vector<Vector3D> vertices{12};
        vertices[1 - 1] = Vector3D::point(0,0,sqrt(5)/2);
        for(int i = 2; i <= 6; i++) {
            double magic = (i-2)*(2*KDB_PI/5);
            vertices[i-1] = Vector3D::point(cos(magic), sin(magic), 0.5);
        }
        for(int i = 7; i <= 11; i++) {
            double magic = (KDB_PI/5)+(i-7)*(2*KDB_PI/5);
            vertices[i-1] = Vector3D::point(cos(magic), sin(magic), -0.5);
        }
        vertices[12 - 1] = Vector3D::point(0,0,-sqrt(5)/2);

        vector<Face> faces = toFacesStartingWithIndex0({Face {{1,2,3}},
                                                        {{1,3,4}},
                                                        {{1,4,5}},
                                                        {{1,5,6}},
                                                        {{1,6,2}},
                                                        {{2,7,3}},
                                                        {{3,7,8}},
                                                        {{3,8,4}},
                                                        {{4,8,9}},
                                                        {{4,9,5}},
                                                        {{5,9,10}},
                                                        {{5,10,6}},
                                                        {{6,10,11}},
                                                        {{6,11,2}},
                                                        {{2,11,7}},
                                                        {{12,8,7}},
                                                        {{12,9,8}},
                                                        {{12,10,9}},
                                                        {{12,11,10}},
                                                        {{12,7,11}}});

        return Figure{vertices, faces, color};
    }

    Figure createOctahedron(const NormColor& color)
    {
        vector<Vector3D> vertices = { Vector3D::point(1,0,0),
                                      Vector3D::point(0,1,0),
                                      Vector3D::point(-1,0,0),
                                      Vector3D::point(0,-1,0),
                                      Vector3D::point(0,0,-1),
                                      Vector3D::point(0,0,1)};

        vector<Face> faces = toFacesStartingWithIndex0(
                {Face{{1,2,6}},{{2,3,6}},{{3,4,6}}, {{4,1,6}}, {{2,1,5}}, {{3,2,5}}, {{4,3,5}}, {{1,4,5}}});

        return {vertices, faces, color};
    }

    Figure createTetrahedron(const NormColor& color)
    {
        vector<Vector3D> vertices = {Vector3D::point(1,-1,-1),
                                     Vector3D::point(-1,1,-1),
                                     Vector3D::point(1,1,1),
                                     Vector3D::point(-1,-1,1)};

        vector<Face> faces = toFacesStartingWithIndex0({Face{{1,2,3}}, {{2,4,3}}, {{1,4,2}}, {{1,3,4}}});
        return {vertices, faces, color};
    }

    Figure createCube(const NormColor& color)
    {
        vector<Vector3D> vertices = {
                Vector3D::point(1, -1, -1),
                Vector3D::point(-1, 1, -1),
                Vector3D::point(1, 1, 1),
                Vector3D::point(-1, -1, 1),
                Vector3D::point(1, 1, -1),
                Vector3D::point(-1, -1, -1),
                Vector3D::point(1, -1, 1),
                Vector3D::point(-1, 1,1)};

        vector<Face> faces = toFacesStartingWithIndex0({Face{{1,5,3,7}}, {{5,2,8,3}}, {{2,6,4,8}}, {{6,1,7,4}}, {{7,3,8,4}}, {{1,6,2,5}}});
        return {vertices, faces, color};
    }

    Figure createCone(int n, double h, const NormColor& color)
    {
        if(n < 3) cerr << "Trying to create a cone with less than 3 faces huh? You're a funny lad\n";

        double r = 1;

        vector<Vector3D> vertices;

        Vector3D top = Vector3D::point(0,0,h);
        vertices.push_back(top);

        for(double a = 0; a < 2*KDB_PI; a += 2*KDB_PI/n) {
            Vector3D circlePoint = Vector3D::point(cos(a) * r, sin(a) * r, 0);
            vertices.push_back(circlePoint);
        }

        vector<Face> faces;
        for(int i = 2; i < vertices.size(); i++) { // index 2 -> p1 on circle (second point)
            faces.push_back(Face{{0,i-1,i}});
        }
        faces.push_back(Face{{0, (int)vertices.size()-1, 1}});

        return {vertices, faces, color};
    }

    Figure createCylinder(int n, double h, const NormColor& c)
    {
        if(n < 3) cerr << "Trying to create a cylinder with less than 3 faces huh? You're a funny lad\n";

        double r = 1;

        vector<Vector3D> vertices;

        // Bottom circle
        for(double a = 0; a < 2*KDB_PI; a += 2*KDB_PI/n) {
            Vector3D circlePoint = Vector3D::point(cos(a) * r, sin(a) * r, 0);
            vertices.push_back(circlePoint);
        }
        // Top circle
        for(double a = 0; a < 2*KDB_PI; a += 2*KDB_PI/n) {
            Vector3D circlePoint = Vector3D::point(cos(a) * r, sin(a) * r, h);
            vertices.push_back(circlePoint);
        }

        vector<Face> faces;
        int indexFirstPointTopCircle = (int)vertices.size()/2;
        for(int i = 1; i < indexFirstPointTopCircle; i++) {
            faces.push_back(Face{{i, i + indexFirstPointTopCircle, i-1 + indexFirstPointTopCircle, i-1}});
        }
        //last face
        faces.push_back(Face{{0, 0 + indexFirstPointTopCircle, (int)vertices.size()-1, indexFirstPointTopCircle-1 }});

        return{vertices, faces, c};
    }

    Figure createTorus(double r, double R, int n, int m, const NormColor& color)
    {
        vector<Vector3D> vertices;

        for(int u_i = 0; u_i < n; u_i++) {

            double u = u_i * (2 * KDB_PI / n);

            for(int v_i = 0; v_i < m; v_i++) {

                double v = v_i * (2 * KDB_PI / m);

                Vector3D p = Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v));
                vertices.push_back(p);
            }
        }
        // n == aantal cirkels waarin de buis wordt opgedeeld
        // m == aantal punten op elke cirkel
        vector<Face> faces;

        for(int u_i = 0; u_i < n; u_i++) {
            for(int v_i = 0; v_i < m; v_i++) {
                faces.push_back(Face{{
                                             u_i * m         + v_i,
                                             (u_i+1)%n * m   + v_i,
                                             (u_i+1)%n * m   + (v_i+1)%m,
                                             u_i * m         + (v_i+1)%m,
                                     }});
            }
        }

        return {vertices, faces, color};
    }
}