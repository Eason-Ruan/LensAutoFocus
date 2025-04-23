//
// Created by 阮立心 on 25-4-20.
//

#ifndef SPECTRAL_RAY_H
#define SPECTRAL_RAY_H

#include "CGL/type.h"
#include "pathtracer/ray.h"

namespace CGL {
    class SpectralRay : Ray {
        Real lambda = 0;
    };
}
#endif //SPECTRAL_RAY_H
