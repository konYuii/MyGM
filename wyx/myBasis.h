#ifndef MYBASIS_H
#define MYBASIS_H
#include "glm/glm.hpp"
#include <vector>
#include <tuple>
#include <vector>
#include <functional>
#include "../include/tinynurbs/tinynurbs.h"

using namespace tinynurbs;

namespace myBasis
{
    template <typename T>
    int findSpan(const std::vector<T> &knots, int degree, T u);

    template <typename T>
    std::vector<T> bsplineBasis(unsigned int deg, int span, const std::vector<T> &knots, T u);

    template <typename T>
    std::vector<T> noZeroBsplineBasis(unsigned int deg, const std::vector<T> &knots, T u);

    template <typename T>
    array2<T> bsplineDerBasis(unsigned int deg, int span, const std::vector<T> &knots, T u, int num_ders);

    template <int dim, typename T>
    glm::vec<dim, T> curvePoint(unsigned int degree, const std::vector<T> &knots,
                                const std::vector<glm::vec<dim, T>> &control_points, T u);

    template <int dim, typename T>
    std::vector<glm::vec<dim, T>> curveDerivatives(unsigned int degree, const std::vector<T> &knots,
                                                   const std::vector<glm::vec<dim, T>> &control_points,
                                                   int num_ders, T u);

    template <typename T>
    glm::vec<3, T> curvePoint(const RationalCurve<T> &crv, T u);

    template <typename T>
    std::vector<glm::vec<3, T>> curveDerivatives(const RationalCurve<T> &crv, int num_ders, T u);

    template <typename T>
    glm::vec<3, T> curveTangent(const RationalCurve<T> &crv, T u);
}

#endif // MYBASIS_H