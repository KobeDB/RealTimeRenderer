#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

#include <iostream>

#include "Camera.h"

using namespace std;

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
	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		return true;
	}

    vector<float> vertices 
    {
    -0.5f, -0.5f, -0.5f,
     0.5f, -0.5f, -0.5f,
     0.5f,  0.5f, -0.5f,
     0.5f,  0.5f, -0.5f,
    -0.5f,  0.5f, -0.5f,
    -0.5f, -0.5f, -0.5f,

    -0.5f, -0.5f,  0.5f,
     0.5f, -0.5f,  0.5f,
     0.5f,  0.5f,  0.5f,
     0.5f,  0.5f,  0.5f,
    -0.5f,  0.5f,  0.5f,
    -0.5f, -0.5f,  0.5f,

    -0.5f,  0.5f,  0.5f,
    -0.5f,  0.5f, -0.5f,
    -0.5f, -0.5f, -0.5f,
    -0.5f, -0.5f, -0.5f,
    -0.5f, -0.5f,  0.5f,
    -0.5f,  0.5f,  0.5f,

     0.5f,  0.5f,  0.5f,
     0.5f,  0.5f, -0.5f,
     0.5f, -0.5f, -0.5f,
     0.5f, -0.5f, -0.5f,
     0.5f, -0.5f,  0.5f,
     0.5f,  0.5f,  0.5f,

    -0.5f, -0.5f, -0.5f,
     0.5f, -0.5f, -0.5f,
     0.5f, -0.5f,  0.5f,
     0.5f, -0.5f,  0.5f,
    -0.5f, -0.5f,  0.5f,
    -0.5f, -0.5f, -0.5f,

    -0.5f,  0.5f, -0.5f,
     0.5f,  0.5f, -0.5f,
     0.5f,  0.5f,  0.5f,
     0.5f,  0.5f,  0.5f,
    -0.5f,  0.5f,  0.5f,
    -0.5f,  0.5f, -0.5f,
    };

    Camera camera;

	bool OnUserUpdate(float fElapsedTime) override
	{
		Clear(olc::BLACK);

        if (GetKey(olc::Key::A).bHeld) camera.Move({ camera.movementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::D).bHeld) camera.Move({ -camera.movementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::W).bHeld) camera.Move({ 0, 0, camera.movementSpeed * fElapsedTime });
        if (GetKey(olc::Key::S).bHeld) camera.Move({ 0, 0, -camera.movementSpeed * fElapsedTime });

        if (GetKey(olc::Key::UP).bHeld) camera.Move({ 0, camera.movementSpeed * fElapsedTime , 0});
        if (GetKey(olc::Key::DOWN).bHeld) camera.Move({ 0, -camera.movementSpeed * fElapsedTime , 0 });

        float s = 10;
        float scaleArr[]{
            s,0,0,0,
            0,s,0,0,
            0,0,s,0,
            0,0,0,1,
        };
        glm::mat4 scale = glm::make_mat4(scaleArr);

        vector<float> projectedVertices;
        projectedVertices.reserve(vertices.size());

        for (int i = 0; i <= vertices.size() - 3; i += 3) {
            glm::vec4 v = { vertices[i], vertices[i + 1] ,vertices[i + 2], 1 };
            v = camera.GetCameraTransform() * v;
            glm::vec2 projected = { -v.x / v.z, -v.y /v.z };
            projectedVertices.push_back(projected.x);
            projectedVertices.push_back(projected.y);
        }
        
        for (int i = 0; i <= projectedVertices.size() - 4; i += 2) {
            glm::vec2 v0 = { projectedVertices[i], projectedVertices[i+1]};
            glm::vec2 v1 = { projectedVertices[i+2], projectedVertices[i + 3] };
            float aspectRatio = ScreenWidth() / ScreenHeight();
            float virtualScreenHeight = 2;
            float virtualScreenWidth = virtualScreenHeight * aspectRatio;
            glm::vec2 offset = { virtualScreenWidth / 2, -1 };
            v0 = v0 + offset;
            v1 = v1 + offset;
            v0 = {v0.x, v0.y + virtualScreenHeight };
            v1 = {v1.x, v1.y + virtualScreenHeight };
            
            v0 = { v0.x * ScreenWidth()/virtualScreenWidth, v0.y * ScreenHeight() / virtualScreenHeight };
            v1 = { v1.x * ScreenWidth() / virtualScreenWidth, v1.y * ScreenHeight() / virtualScreenHeight };

            DrawLine({ int(v0.x), int(v0.y) }, { int(v1.x), int(v1.y) });
        }

		return true;
	}
};

int main()
{
	Example demo;
	if (demo.Construct(400, 300, 2, 2))
		demo.Start();
	return 0;
}