//
// Created by 阮立心 on 25-4-19.
//

#ifndef CAMERA_LENSSYS_H
#define CAMERA_LENSSYS_H

#include "camera.h"
#include "lens-sampler/random.h"
#include "lens-system/lens-system.h"
namespace CGL {
    class CameraLensSys final : public Camera {
    public:
        // TODO: Respecify the len system resolution
        CameraLensSys(const std::string& lensFile, const double width, const double height) : Camera(), lensSys(nullptr)  {
            // 设置相机参数，使其与镜头系统匹配
            lensSys  = new LensSystem(lensFile, 36, 27);
        }
        Ray generate_ray(double x, double y) const override;
        LensSystem* lensSys;
        Prl2::Sampler* random_sampler;
        ~CameraLensSys() override {
            delete lensSys;
        }
    };
}




#endif //CAMERA_LENSSYS_H
