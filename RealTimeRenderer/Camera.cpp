#include "Camera.h"

void Camera::Move(glm::vec3 dir)
{
    eye.x += dir.x;
    eye.y += dir.y;
    eye.z += dir.z;
    center.x += dir.x;
    center.y += dir.y;
    center.z += dir.z;
}

glm::mat4 Camera::GetCameraTransform() const
{
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

    return eyeRotation * eyeTranslation;
}