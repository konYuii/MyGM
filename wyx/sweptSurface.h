#ifndef SWEPTSURFACE_H
#define SWEPTSURFACE_H
#include "eigen/Eigen/Core"
#include "eigen/Eigen/Dense"
#include <iostream>
#include "myBasis.h"
#include "frame.h"

namespace sweptsurface
{
    template <typename T>
    void SurfMeshParams(int n, int m, std::vector<std::vector<glm::vec<3, T>>> &surf_points,
                        std::vector<T> &u_params, std::vector<T> &v_params);

    template <typename T>
    void GenerateKnotVector(std::vector<T> &knotvector, std::vector<T> &params, int p, int n);

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface0(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames);

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface1(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames);

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface2(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller, 
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames);

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface3(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames);

    template <typename T>
    void GenerateProfileCurves(std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                               const tinynurbs::RationalCurve<T> &contour,
                               const tinynurbs::RationalCurve<T> &trace,
                               std::function<glm::mat4(T)> shape_controller,
                               const std::vector<T> &v_bar,
                               const std::vector<std::vector<glm::vec3>> &frames);

    static glm::mat4 transform0(float t)
    {
        return glm::mat4(1.0f);
    }

    static glm::mat4 transform1(float t)
    {
        glm::mat4 shape_control((float)std::cos(M_PI * t), (float)-std::sin(M_PI * t), 0.0f, 0.0f,
                                (float)std::sin(M_PI * t), (float)std::cos(M_PI * t), 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f);

        return shape_control;
    }

    static glm::mat4 transform2(float t)
    {
        glm::mat4 shape_control(2.0f * (float)std::cos(M_PI * t), (float)-std::sin(M_PI * t), 0.0f, 0.0f,
                                (float)std::sin(M_PI * t), 3.0f * (float)std::cos(M_PI * t), 0.0f, 0.0f,
                                0.0f, 0.0f, 1.0f, 0.0f,
                                0.0f, 0.0f, 0.0f, 1.0f);
        
        return shape_control;
    }
}

#endif // SWEPTSURFACE_H