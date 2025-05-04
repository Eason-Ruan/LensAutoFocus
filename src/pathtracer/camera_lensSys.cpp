//
// Created by 阮立心 on 25-4-19.
//
#include "camera_lensSys.h"

#include "CGL/spectrum.h"
#include "lens-system/lens-system.h"

namespace CGL {
    Ray CameraLensSys::generate_ray(const double x, const double y) const {
        constexpr Real lambda_pdf = 1.0f / (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN);
        // 确保生成的波长在有效范围内，避免溢出
        Real lambda = random_sampler->getNext() * (SPD::LAMBDA_MAX - SPD::LAMBDA_MIN - 0.001) + SPD::LAMBDA_MIN;
        
        // 添加额外边界检查，确保波长在 LAMBDA_MIN 到 LAMBDA_MAX 之间
        if (lambda >= SPD::LAMBDA_MAX) {
            lambda = SPD::LAMBDA_MAX - 0.001;
        }
        if (lambda < SPD::LAMBDA_MIN) {
            lambda = SPD::LAMBDA_MIN;
        }
        
        //采样光线
        Ray ray_out;
        ray_out.o = Vector3D(0, 0, 0);
        ray_out.d = Vector3D(0, 0, 0);
        Real ray_pdf;
        const bool result = lensSys->sample_ray(x, y, lambda, *random_sampler, ray_out, ray_pdf, false);
        if (!result) {
            ray_out.d = Vector3D(0, 0, 0);
            return ray_out;
        }
        if (ray_out.d == Vector3D(0, 0, 0)) {
            std::cerr << "Error sampling ray" << std::endl;
        }
        ray_out.o /= 1000.0;
        ray_out.o =  c2w * ray_out.o + pos;
        ray_out.d = c2w * ray_out.d;
        ray_out.d.normalize();
        ray_out.lambda = lambda;
        ray_out.ray_pdf = ray_pdf;
        ray_out.lambda_pdf = lambda_pdf;
        return ray_out;
    }

    void CameraLensSys::focus_delta(const double delta_distance) {
        lensSys->focus_delta(delta_distance);
        lensSys->compute_exit_pupil_bounds();
        fprintf(stdout, "[CameraLensSys] Focused lens system by object delta distance %.2f mm \n", delta_distance); fflush(stdout);
    }

    void CameraLensSys::focus_mechanical(const double delta_distance) {
        lensSys->focus_mechanical_delta(delta_distance);
        lensSys->compute_exit_pupil_bounds();
        fprintf(stdout, "[CameraLensSys] Focused lens system by mechanical delta distance %.2f mm \n", delta_distance); fflush(stdout);
    }
}
