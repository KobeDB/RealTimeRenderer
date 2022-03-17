#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

#include <iostream>

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

    glm::vec3 eye{ 0,0, -10 };
    glm::vec3 center{ 0,0,-1 };
    glm::vec3 up{0,1,0};
    float eyeMovementSpeed = 4.0f;

    void MoveCamera(glm::vec3 dir)
    {
        eye.x += dir.x;
        eye.y += dir.y;
        eye.z += dir.z;
        center.x += dir.x;
        center.y += dir.y;
        center.z += dir.z;
    }

	bool OnUserUpdate(float fElapsedTime) override
	{
		Clear(olc::BLACK);

        if (GetKey(olc::Key::A).bHeld) MoveCamera({ eyeMovementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::D).bHeld) MoveCamera({ -eyeMovementSpeed * fElapsedTime, 0, 0 });
        if (GetKey(olc::Key::W).bHeld) MoveCamera({ 0, 0, eyeMovementSpeed * fElapsedTime });
        if (GetKey(olc::Key::S).bHeld) MoveCamera({ 0, 0, -eyeMovementSpeed * fElapsedTime });

        if (GetKey(olc::Key::UP).bHeld)MoveCamera({ 0, eyeMovementSpeed * fElapsedTime , 0});
        if (GetKey(olc::Key::DOWN).bHeld) MoveCamera({ 0, -eyeMovementSpeed * fElapsedTime , 0 });

        if (GetKey(olc::Key::J).bHeld) eye.y += eyeMovementSpeed * fElapsedTime;
        if (GetKey(olc::Key::L).bHeld) eye.y -= eyeMovementSpeed * fElapsedTime;

        float s = 10;
        float scaleArr[]{
            s,0,0,0,
            0,s,0,0,
            0,0,s,0,
            0,0,0,1,
        };
        glm::mat4 scale = glm::make_mat4(scaleArr);

        float eyeTranslationArr[]{
            1,  0,  0,  -eye.x,
            0,  1,  0,  -eye.y,
            0,  0,  1,  -eye.z,
            0,  0,  0,  1
        };
        // transpose because my source matrix has consecutive row elements,
        // make_mat4 instead interprets this as consecutive column elements!
        glm::mat4 eyeTranslation = glm::transpose(glm::make_mat4(eyeTranslationArr));

        glm::vec3 zc = glm::normalize(eye - center);
        glm::vec3 xc = glm::normalize(glm::cross(up, zc));
        glm::vec3 yc = glm::cross(zc, xc);
        
        float eyeRotationArr[]{
            xc.x,  xc.y,  xc.z,  0,
            yc.x,  yc.y,  yc.z,  0,
            zc.x,  zc.y,  zc.z,  0,
            0,  0,  0,           1,
        };
        // transpose because my source matrix has consecutive row elements,
        // make_mat4 instead interprets this as consecutive column elements!
        glm::mat4 eyeRotation = glm::transpose(glm::make_mat4(eyeRotationArr)); 
        

        vector<float> projectedVertices;
        projectedVertices.reserve(vertices.size());

        for (int i = 0; i <= vertices.size() - 3; i += 3) {
            glm::vec4 v = { vertices[i], vertices[i + 1] ,vertices[i + 2], 1 };
            v = eyeRotation * eyeTranslation * v;
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