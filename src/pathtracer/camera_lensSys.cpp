//
// Created by 阮立心 on 25-4-19.
//
#include "camera_lensSys.h"

#include "CGL/spectrum.h"
#include "lens-system/lens-system.h"

namespace CGL {
    Ray CameraLensSys::generate_ray(double x, double y) const {
        double u = (2.0 * (samplePoint.x + sampler->getNext()) - samplerBuffer.w) / (double) sampleBuffer.w;
        double v = (2.0 * (samplePoint.y + sampler->getNext()) - samplerBuffer.h) / (double) sampleBuffer.h;

        //采样波长
        double lambda_pdf = 1.0 / (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN);
        double lambda = SPD::LAMBDA_MIN + (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN) * sampler->getNext();

        //采样光线
        Ray ray_out;
        double ray_pdf;
        if (!lensSystem->sampleRay(u, v, lambda, *sampler, ray, ray_pdf, false)) {
            continue; // 采样失败跳过
        }
        return ray_out;
    }

}
