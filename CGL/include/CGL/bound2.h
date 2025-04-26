//
// Created by Lixin Ruan on 25-4-15.
//

#ifndef BOUND2_H
#define BOUND2_H
#include "CGL/vector2D.h"
#include "CGL/type.h"

#include "util/random_util.h"

namespace CGL {
    class Bounds2 {
    public:
        Vector2D p0;
        Vector2D p1;

        Bounds2()
            : p0(Vector2D(std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max())),
              p1(Vector2D(std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest())) {}
        Bounds2(const Vector2D& _p0, const Vector2D& _p1) : p0(_p0), p1(_p1) {}

        Real area() const { return (p1.x - p0.x) * (p1.y - p0.y); }

        bool isValid() const { return p0.x < p1.x && p0.y < p1.y; }

        Vector2D samplePoint(Real& pdf) const {
            pdf = 1.0 / area();
            const double u = random_uniform();
            return p0 + u * (p1 - p0);
        }
    };

    inline Bounds2 extendBounds(const Bounds2& b, const Vector2D& p) {
        const Real p0x = std::min(b.p0.x, p.x);
        const Real p0y = std::min(b.p0.y, p.y);
        const Real p1x = std::max(b.p1.x, p.x);
        const Real p1y = std::max(b.p1.y, p.y);
        return {Vector2D(p0x, p0y), Vector2D(p1x, p1y)};
    }

    inline std::ostream& operator<<(std::ostream& stream, const Bounds2& bounds) {
        stream << bounds.p0 << ", " << bounds.p1;
        return stream;
    }
}
#endif //BOUND2_H
