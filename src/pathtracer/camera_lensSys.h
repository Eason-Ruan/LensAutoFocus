//
// Created by 阮立心 on 25-4-19.
//

#ifndef CAMERA_LENSSYS_H
#define CAMERA_LENSSYS_H

#include "camera.h"

namespace CGL {
    class CameraLensSys : public Camera {
        public:
        Ray generate_ray(double x, double y) const;
    };
}




#endif //CAMERA_LENSSYS_H
