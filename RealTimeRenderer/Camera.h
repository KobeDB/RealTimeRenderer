#pragma once

#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

struct Camera {
    glm::vec3 eye{ 0, 0, -10 };
    glm::vec3 center{ 0, 0, -1 };
    glm::vec3 up{ 0, 1, 0 };
    float movementSpeed = 4.0f;

    void Move(glm::vec3 dir);

    glm::mat4 GetCameraTransform() const;
};
