#include "frame.h"

namespace frame {
    glm::mat4 getFrameMatrix(const std::vector<glm::vec3> &frame)
    {
        glm::vec3 X = frame[0];
        glm::vec3 Y = frame[1];
        glm::vec3 Z = frame[2];
        glm::mat4 frameMat = glm::mat4(glm::vec4(X, 0.0f), glm::vec4(Y, 0.0f), glm::vec4(Z, 0.0f), glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
        return glm::transpose(frameMat);
    }

    glm::mat4 getInverseFrameMatrix(const std::vector<glm::vec3> &frame)
    {
        return glm::transpose(getFrameMatrix(frame));
    }

    template <typename T>
    void computeFrames1(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params)
    {
        for(int i=0;i<u_params.size();i++)
        {
            std::vector<glm::vec3> frame(3);
            float u = u_params[i];
            std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u);
            glm::vec3 Bi = glm::normalize(glm::cross(ts[1], ts[2]));
            glm::vec3 Ta = glm::normalize(ts[1]);
            glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
            frame[0] = No;
            frame[1] = Bi;
            frame[2] = Ta;
            frames.push_back(frame);
            frames_saved.push_back(std::make_pair(frame, u));
        }
    }

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame1(const tinynurbs::RationalCurve<T> &curve, const T &u_param)
    {
        std::vector<glm::vec3> frame(3);
        std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u_param);
        glm::vec3 Bi = glm::normalize(glm::cross(ts[1], ts[2]));
        glm::vec3 Ta = glm::normalize(ts[1]);
        glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
        frame[0] = No;
        frame[1] = Bi;
        frame[2] = Ta;
        return frame;
    }

    template <typename T>
    void computeFrames2(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params)
    {
        for(int i=0;i<u_params.size();i++)
        {
            std::vector<glm::vec3> frame(3);
            float u = u_params[i];
            if(!i)
            {
                std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u);
                glm::vec3 Bi = glm::normalize(glm::cross(ts[1], ts[2]));
                glm::vec3 Ta = glm::normalize(ts[1]);
                glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
                frame[0] = No;
                frame[1] = Bi;
                frame[2] = Ta;
                frames.push_back(frame);
            }
            else
            {
                std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u);
                glm::vec3 Ta = glm::normalize(ts[1]);
                glm::vec3 bi = frames[i-1][1] - (glm::dot(frames[i-1][1], Ta) * Ta);
                glm::vec3 Bi = glm::normalize(bi);
                glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
                frame[0] = No;
                frame[1] = Bi;
                frame[2] = Ta;
                frames.push_back(frame);
            }
            frames_saved.push_back(std::make_pair(frame, u));
        }
        // 处理环状曲线Bi
    }

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame2(const tinynurbs::RationalCurve<T> &curve, const T &u_param)
    {
        std::vector<glm::vec3> frame(3);
        for(int i=0;i<frames_saved.size();i++)
        {
            if(frames_saved[i].second >= u_param)
            {
                if(i)
                {
                    std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 1, u_param);
                    glm::vec3 Ta = glm::normalize(ts[1]);
                    glm::vec3 bi = frames_saved[i-1].first[1] - (glm::dot(frames_saved[i-1].first[1], Ta) * Ta);
                    glm::vec3 Bi = glm::normalize(bi);
                    glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
                    frame[0] = No;
                    frame[1] = Bi;
                    frame[2] = Ta;
                    return frame;
                }
                else
                {
                    return frames_saved[0].first;
                }
            }            
        }
        return frames_saved[frames_saved.size()-1].first;
    }

    //显示模板化函数
    template void computeFrames1(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<float> &curve, const std::vector<float> &u_params);
    template void computeFrames2(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<float> &curve, const std::vector<float> &u_params);
    template std::vector<glm::vec3> computeSingleFrame1(const tinynurbs::RationalCurve<float> &curve, const float &u_param);
    template std::vector<glm::vec3> computeSingleFrame2(const tinynurbs::RationalCurve<float> &curve, const float &u_param);

}