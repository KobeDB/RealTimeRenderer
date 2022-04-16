#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

#include <iostream>

#include "Camera.h"

#include "render_facilities.h"
#include "3d_figures.h"

using namespace std;
using namespace KDBRenderFacilities;

struct Line {
    glm::vec4 a, b;
};



// Override base class with your custom functionality
class Example : public olc::PixelGameEngine
{
public:
	Example()
	{
		// Name your application
		sAppName = "Example";
	}

public:

    Figure fig;

    glm::mat4 perspectiveProjectionTransform;

    glm::mat4 perspectiveTransform;
    glm::mat4 orthographicProjectionTransform;
    glm::mat4 viewportTransform;

	bool OnUserCreate() override
	{
        double nx = ScreenWidth();
        double ny = ScreenHeight();

        float viewportTransformArr[]{
            nx / 2,   0,      0,  (nx - 1) / 2,
            0,      ny / 2,   0,  (ny - 1) / 2,
            0,      0,      1,  0,
            0,      0,      0,  1,
        };
        viewportTransform = glm::transpose(glm::make_mat4(viewportTransformArr));

        double fov{ 120 };
        double n{ -1 }, f{ -5 };
        double t{ tan(fov / 2 * abs(n)) }, b{ -t };
        double r{ t * (nx / ny) }, l{ -r };

        float orthographicProjectionTransformArr[]{
            2 / (r - l),   0,     0,          -(r + l) / (r - l),
            0,      2 / (t - b),  0,          -(t + b) / (t - b),
            0,      0,          2 / (n - f),    -(n + f) / (n - f),
            0,      0,          0,          1,
        };
        orthographicProjectionTransform = glm::transpose(glm::make_mat4(orthographicProjectionTransformArr));

        float perspectiveTransformArr[]{
            n,      0,      0,      0,
            0,      n,      0,      0,
            0,      0,      n + f,    -n * f,
            0,      0,      1,      0,
        };
        perspectiveTransform = glm::transpose(glm::make_mat4(perspectiveTransformArr));

        perspectiveProjectionTransform = viewportTransform * orthographicProjectionTransform * perspectiveTransform;

        //fig = createSphere(1, 2, magenta);
        fig = KDBRenderFacilities::createTorus(1, 5, 30, 30, magenta);


		return true;
	}

    Camera camera;

    void drawFigure(const Figure& fig)
    {
        for (const auto& face : fig.faces) {
            for (int i = 0; i < face.point_indexes.size(); i++) {
                int p0_i = face.point_indexes[i];
                int p1_i = face.point_indexes[(i + 1) % face.point_indexes.size()]; //This loops around to 0, to draw a line between the first and the last vertex
                const Vector3D& p0 = fig.vertices[p0_i];
                const Vector3D& p1 = fig.vertices[p1_i];
                glm::vec4 homo_p0 = { p0.x, p0.y, p0.z, 1 };
                glm::vec4 homo_p1 = { p1.x, p1.y, p1.z, 1 };
                glm::mat4 M = perspectiveProjectionTransform * camera.GetCameraTransform();
                glm::vec4 p = M * homo_p0;
                glm::vec4 q = M * homo_p1;
                DrawLine(p.x / p.w, p.y / p.w, q.x / q.w, q.y / q.w);
            }
        }
    }

	bool OnUserUpdate(float fElapsedTime) override
	{
		Clear(olc::BLACK);

        if (GetKey(olc::Key::A).bHeld) camera.Move({ camera.movementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::D).bHeld) camera.Move({ -camera.movementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::W).bHeld) camera.Move({ 0, 0, camera.movementSpeed * fElapsedTime });
        if (GetKey(olc::Key::S).bHeld) camera.Move({ 0, 0, -camera.movementSpeed * fElapsedTime });

        if (GetKey(olc::Key::UP).bHeld) camera.Move({ 0, camera.movementSpeed * fElapsedTime , 0});
        if (GetKey(olc::Key::DOWN).bHeld) camera.Move({ 0, -camera.movementSpeed * fElapsedTime , 0 });

        drawFigure(fig);

        /*for (const Line& l : lines) {
            glm::mat4 M = viewportTransform * orthographicProjectionTransform * perspectiveTransform * camera.GetCameraTransform();
            glm::vec4 p =  M * l.a;
            glm::vec4 q = M * l.b;
            DrawLine(p.x / p.w, p.y / p.w, q.x / q.w, q.y / q.w);
        }*/

		return true;
	}
};

int main()
{
	Example demo;
	if (demo.Construct(400, 200, 3, 3))
		demo.Start();
	return 0;
}