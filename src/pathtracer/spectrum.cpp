#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>

#include "core/spectrum.h"

namespace Prl2 {
/* lens */
//均匀波长采样
inline Real sample_wavelength(Sampler& sampler, Real lambda) {
  // 采样的波长范围
  static const Real LAMBDA_MIN = SPD::LAMBDA_MIN;
  static const Real LAMBDA_MAX = SPD::LAMBDA_MAX;

  // 计算采样的波长
  return sampler.getNext() * (LAMBDA_MAX - LAMBDA_MIN) + LAMBDA_MIN;
}

// 从非等间隔的光谱功率分布(SPD)构造等间隔的SPD
// 通过线性插值计算非等间隔的SPD
SPD::SPD(const std::vector<Real>& _lambda, const std::vector<Real>& _phi) {
  assert(_lambda.size() == _phi.size());  // 确保波长和功率数组大小一致

  // 如果只有一个样本
  if (_lambda.size() == 1) {
    // 初始化其他波长的功率为0
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      phi[i] = 0;
    }
    addPhi(_lambda[0], _phi[0]);  // 添加单个样本的功率
  }
  // 如果有多个样本
  else {
    for (size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
      // 计算等间隔波长
      const Real lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;

      // 如果波长超出非等间隔SPD的范围，功率设为0
      if (lambda_value < _lambda.front() || lambda_value > _lambda.back()) {
        phi[i] = 0;
      } else if (lambda_value == _lambda.front()) {
        phi[i] = _phi.front();  // 如果波长等于最小值，直接取对应功率
      } else if (lambda_value == _lambda.back()) {
        phi[i] = _phi.back();  // 如果波长等于最大值，直接取对应功率
      } else {
        // 线性插值计算非等间隔SPD的功率
        const size_t lambda1_index = static_cast<size_t>(
            std::lower_bound(_lambda.begin(), _lambda.end(), lambda_value) -
            _lambda.begin());  // 找到大于等于当前波长的索引
        const size_t lambda0_index = lambda1_index - 1;  // 前一个索引
        assert(lambda0_index != lambda1_index);

        // 计算插值因子
        const Real t = (lambda_value - _lambda[lambda0_index]) /
                       (_lambda[lambda1_index] - _lambda[lambda0_index]);
        assert(t >= 0 && t <= 1);

        // 插值计算功率
        const Real interpolated_phi =
            (1.0f - t) * _phi[lambda0_index] + t * _phi[lambda1_index];
        phi[i] = interpolated_phi;
      }
    }
  }
}

// 添加功率到指定波长
void SPD::addPhi(const Real& _lambda, const Real& _phi) {
  // 如果波长超出范围，不进行加算
  if (_lambda < LAMBDA_MIN || _lambda >= LAMBDA_MAX) {
    return;
  }

  // 计算对应波长的索引
  const size_t lambda_index = (_lambda - LAMBDA_MIN) / LAMBDA_INTERVAL;

  // 将功率分配到两侧波长
  const Real lambda0 =
      LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;  // 左侧波长
  const Real t =
      (_lambda - lambda0) / LAMBDA_INTERVAL;  // 波长在区间中的位置
  phi[lambda_index] += (1.0f - t) * _phi;  // 分配到左侧波长
  phi[lambda_index + 1] += t * _phi;      // 分配到右侧波长
}

// 采样指定波长的功率
Real SPD::sample(const Real& l) const {
  assert(l >= LAMBDA_MIN && l < LAMBDA_MAX);  // 确保波长在范围内

  // 计算对应波长的索引
  const size_t lambda_index = (l - LAMBDA_MIN) / LAMBDA_INTERVAL;

  // 线性插值计算功率
  if (lambda_index == 0 || lambda_index == LAMBDA_SAMPLES - 1) {
    return phi[lambda_index];  // 如果在边界，直接返回功率
  } else {
    const Real lambda_nearest = LAMBDA_MIN + lambda_index * LAMBDA_INTERVAL;
    const Real t = (l - lambda_nearest) / LAMBDA_INTERVAL;  // 插值因子
    assert(t >= 0 && t <= 1);
    return (1.0f - t) * phi[lambda_index] + t * phi[lambda_index + 1];
  }
}

// 将SPD转换为XYZ颜色空间
XYZ SPD::toXYZ() const {
  XYZ xyz;

  for (std::size_t i = 0; i < LAMBDA_SAMPLES; ++i) {
    const Real lambda_value = LAMBDA_MIN + LAMBDA_INTERVAL * i;  // 当前波长
    const Real phi_value = phi[i];  // 当前波长的功率

    // 如果功率为0，跳过计算
    if (phi_value == 0) continue;

    // 计算对应等色函数的索引
    const int index = (lambda_value - 380) / 5;
    assert(index >= 0 && index < color_matching_func_samples);

    // 线性插值等色函数
    Real cmf_x, cmf_y, cmf_z;
    if (index != color_matching_func_samples - 1) {
      const Real cmf_lambda = 5 * index + 380;
      const Real t = (lambda_value - cmf_lambda) / 5;
      assert(t >= 0 && t <= 1);

      cmf_x = (1.0f - t) * color_matching_func_x[index] +
              t * color_matching_func_x[index + 1];
      cmf_y = (1.0f - t) * color_matching_func_y[index] +
              t * color_matching_func_y[index + 1];
      cmf_z = (1.0f - t) * color_matching_func_z[index] +
              t * color_matching_func_z[index + 1];
    } else {
      cmf_x = color_matching_func_x[index];
      cmf_y = color_matching_func_y[index];
      cmf_z = color_matching_func_z[index];
    }

    // 使用短册近似计算XYZ值
    xyz[0] += cmf_x * phi_value;
    xyz[1] += cmf_y * phi_value;
    xyz[2] += cmf_z * phi_value;
  }

  return xyz;
}

// 将RGB颜色转换为光谱表示
SPD RGB2Spectrum(const RGB& rgb) {
  // 采样的波长
  static const std::vector<Real> sampled_lambda = {
      380, 417.7, 455.55, 493.33, 531.11, 568.88, 606.66, 644.44, 682.22, 720};

  // 定义一些基础光谱
  static const SPD white_spectrum =
      SPD(sampled_lambda, {1, 1, 0.9999, 0.9993, 0.9992, 0.9998, 1, 1, 1, 1});
  static const SPD cyan_spectrum =
      SPD(sampled_lambda,
          {0.9710, 0.9426, 1.0007, 1.0007, 1.0007, 1.0007, 0.1564, 0, 0, 0});
  static const SPD magenta_spectrum = SPD(
      sampled_lambda, {1, 1, 0.9685, 0.2229, 0, 0.0458, 0.8369, 1, 1, 0.9959});
  static const SPD yellow_spectrum =
      SPD(sampled_lambda,
          {0.0001, 0, 0.1088, 0.6651, 1, 1, 0.9996, 0.9586, 0.9685, 0.9840});
  static const SPD red_spectrum =
      SPD(sampled_lambda,
          {0.1012, 0.0515, 0, 0, 0, 0, 0.8325, 1.0149, 1.0149, 1.0149});
  static const SPD green_spectrum = SPD(
      sampled_lambda, {0, 0, 0.0273, 0.7937, 1, 0.9418, 0.1719, 0, 0, 0.0025});
  static const SPD blue_spectrum =
      SPD(sampled_lambda,
          {1, 1, 0.8916, 0.3323, 0, 0, 0.0003, 0.0369, 0.0483, 0.0496});

  SPD ret;
  const Real red = rgb.x();   // 提取红色分量
  const Real green = rgb.y(); // 提取绿色分量
  const Real blue = rgb.z();  // 提取蓝色分量

  // 根据RGB分量的大小关系，组合基础光谱
  if (red <= green && red <= blue) {
    ret += red * white_spectrum;
    if (green <= blue) {
      ret += (green - red) * cyan_spectrum;
      ret += (blue - green) * blue_spectrum;
    } else {
      ret += (blue - red) * cyan_spectrum;
      ret += (green - blue) * green_spectrum;
    }
  } else if (green <= red && green <= blue) {
    ret += green * white_spectrum;
    if (red <= blue) {
      ret += (red - green) * magenta_spectrum;
      ret += (blue - red) * blue_spectrum;
    } else {
      ret += (blue - green) * magenta_spectrum;
      ret += (red - blue) * red_spectrum;
    }
  } else {
    ret += blue * white_spectrum;
    if (red <= green) {
      ret += (red - blue) * yellow_spectrum;
      ret += (green - red) * green_spectrum;
    } else {
      ret += (green - blue) * yellow_spectrum;
      ret += (red - green) * red_spectrum;
    }
  }

  return ret;
}

}  // namespace Prl2