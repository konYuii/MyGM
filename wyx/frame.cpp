#include "frame.h"
#include <iostream>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

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

    template <typename T>
    void computeFrames3(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<T> &curve, const std::vector<T> &u_params)
    {
        float y0[3];
        for(int i=0;i<u_params.size();i++)
        {
            std::vector<glm::vec3> frame(3);
            std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u_params[i]);
            if(!i)
            {
                glm::vec3 Bi = glm::normalize(glm::cross(ts[1], ts[2]));
                glm::vec3 Ta = glm::normalize(ts[1]);
                glm::vec3 No = glm::normalize(glm::cross(Bi, Ta));
                frame[0] = No;
                frame[1] = Bi;
                frame[2] = Ta;
                frames.push_back(frame);

                y0[0] = No.x;
                y0[1] = No.y;
                y0[2] = No.z;
                //std::cout<<y0[0]<<'\t'<<y0[1]<<'\t'<<y0[2]<<std::endl;
            }
            else
            {
                double t = 0.0;   // 初始时间
                double t1 = u_params[i]; // 结束时间
                double h = 0.01;   // 时间步长

                const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
                gsl_odeiv_step *s = gsl_odeiv_step_alloc(type, 3);
                gsl_odeiv_control *c = gsl_odeiv_control_standard_new(1e-6, 0.0, 0.0, 0.0);
                gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(3);

                double y[3] = {y0[0], y0[1], y0[2]};

                auto func = [](double t, const double y[], double df[], void *curve) -> int {
                    tinynurbs::RationalCurve<float> *c = (tinynurbs::RationalCurve<float> *)curve;
                    std::vector<glm::vec3> ts = myBasis::curveDerivatives(*c, 2, (float)t);
                    df[0] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].x  / glm::dot(ts[1], ts[1]);
                    df[1] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].y  / glm::dot(ts[1], ts[1]);
                    df[2] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].z  / glm::dot(ts[1], ts[1]);
                    return GSL_SUCCESS;
                };
                
                gsl_odeiv_system sys = {func, nullptr, 3, (void *)&curve};
                //std::cout << "break point1" << std::endl;
                // 进行积分
                while (t < t1)
                {
                    int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
                    if (status != GSL_SUCCESS)
                    {
                        std::cerr << "Integration error, status = " << status << std::endl;
                        break;
                    }

                    // 输出结果
                    // std::cout << t << '\t' << y[0] << '\t' << y[1] << '\t' << y[2] << std::endl;
                }

                // 释放资源
                gsl_odeiv_evolve_free(e);
                gsl_odeiv_control_free(c);
                gsl_odeiv_step_free(s);

                glm::vec3 Ta = glm::normalize(ts[1]);
                glm::vec3 No = glm::normalize(glm::vec3(y[0], y[1], y[2]));
                glm::vec3 Bi = glm::cross(Ta, No);
                frame[0] = No;
                frame[1] = Bi;
                frame[2] = Ta;
                frames.push_back(frame);
            }
            frames_saved.push_back(std::make_pair(frame, u_params[i]));
        }
    }

    template <typename T>
    std::vector<glm::vec3> computeSingleFrame3(const tinynurbs::RationalCurve<T> &curve, const T &u_param)
    {
        if(u_param < 1e-10)
            return frames_saved[0].first;

        std::vector<glm::vec3> frame(3);
        std::vector<glm::vec3> ts = myBasis::curveDerivatives(curve, 2, u_param);

        float y0[3] = {frames_saved[0].first[0].x, frames_saved[0].first[0].y, frames_saved[0].first[0].z};

        double t = 0.0;          // 初始时间
        double t1 = u_param;     // 结束时间
        double h = 0.01;         // 时间步长

        const gsl_odeiv_step_type *type = gsl_odeiv_step_rk8pd;
        gsl_odeiv_step *s = gsl_odeiv_step_alloc(type, 3);
        gsl_odeiv_control *c = gsl_odeiv_control_standard_new(1e-6, 0.0, 0.0, 0.0);
        gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(3);

        double y[3] = {y0[0], y0[1], y0[2]};

        auto func = [](double t, const double y[], double df[], void *curve) -> int
        {
            tinynurbs::RationalCurve<float> *c = (tinynurbs::RationalCurve<float> *)curve;
            std::vector<glm::vec3> ts = myBasis::curveDerivatives(*c, 2, (float)t);
            df[0] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].x / glm::dot(ts[1], ts[1]);
            df[1] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].y / glm::dot(ts[1], ts[1]);
            df[2] = -(ts[2].x * y[0] + ts[2].y * y[1] + ts[2].z * y[2]) * ts[1].z / glm::dot(ts[1], ts[1]);
            return GSL_SUCCESS;
        };

        gsl_odeiv_system sys = {func, nullptr, 3, (void *)&curve};
        // std::cout << "break point1" << std::endl;
        //  进行积分
        while (t < t1)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
            if (status != GSL_SUCCESS)
            {
                std::cerr << "Integration error, status = " << status << std::endl;
                break;
            }

            // 输出结果
            // std::cout << t << '\t' << y[0] << '\t' << y[1] << '\t' << y[2] << std::endl;
        }

        // 释放资源
        gsl_odeiv_evolve_free(e);
        gsl_odeiv_control_free(c);
        gsl_odeiv_step_free(s);

        glm::vec3 Ta = glm::normalize(ts[1]);
        glm::vec3 No = glm::normalize(glm::vec3(y[0], y[1], y[2]));
        glm::vec3 Bi = glm::cross(Ta, No);
        frame[0] = No;
        frame[1] = Bi;
        frame[2] = Ta;

        return frame;
    }

    //显示模板化函数
    template void computeFrames1(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<float> &curve, const std::vector<float> &u_params);
    template void computeFrames2(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<float> &curve, const std::vector<float> &u_params);
    template void computeFrames3(std::vector<std::vector<glm::vec3>> &frames, const tinynurbs::RationalCurve<float> &curve, const std::vector<float> &u_params);
    template std::vector<glm::vec3> computeSingleFrame1(const tinynurbs::RationalCurve<float> &curve, const float &u_param);
    template std::vector<glm::vec3> computeSingleFrame2(const tinynurbs::RationalCurve<float> &curve, const float &u_param);
    template std::vector<glm::vec3> computeSingleFrame3(const tinynurbs::RationalCurve<float> &curve, const float &u_param);

}