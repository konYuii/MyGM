#ifndef FRAME_H
#define FRAME_H

#include "myBasis.h"

namespace frame
{
    static std::vector<std::pair<std::vector<glm::vec3>, float>> frames_saved;

    glm::mat4 getFrameMatrix(const std::vector<glm::vec3> &frame);

    glm::mat4 getInverseFrameMatrix(const std::vector<glm::vec3> &frame);

    template <typename T>
    void computeFrames1(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params);

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame1(const tinynurbs::RationalCurve<T> &curve, const T &u_param);

    template <typename T>
    void computeFrames2(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params);

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame2(const tinynurbs::RationalCurve<T> &curve, const T &u_param);

    template <typename T>
    void computeFrames3(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params);

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame3(const tinynurbs::RationalCurve<T> &curve, const T &u_param);
}

#endif