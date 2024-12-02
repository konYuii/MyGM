#include "sweptSurface.h"
#include <functional>

namespace sweptsurface
{
    template <typename T>
    void SurfMeshParams(int n, int m, std::vector<std::vector<glm::vec<3, T>>> &surf_points,
                        std::vector<T> &u_params, std::vector<T> &v_params)
    {
        int num = m + 1;
        u_params[0] = 0.0f;
        u_params[n] = 1.0f;
        std::vector<T> cds(n + 1);
        for (int l = 0; l <= m; l++)
        {
            float total = 0.0f;
            for (int k = 1; k <= n; k++)
            {
                cds[k] = glm::distance(surf_points[k - 1][l], surf_points[k][l]);
                total += cds[k];
            }
            if (total == 0.0f)
                num -= 1;
            else
            {
                float d = 0.0f;
                for (int k = 1; k < n; k++)
                {
                    d += cds[k];
                    u_params[k] += d / total;
                }
            }
        }
        if (num == 0)
        {
            std::cout << "Error: Invalid surface mesh parameters." << std::endl;
            exit(1);
        }
        for (int k = 1; k < n; k++)
            u_params[k] /= num;

        num = n + 1;
        v_params[0] = 0.0f;
        v_params[m] = 1.0f;
        cds.clear();
        cds.resize(m + 1);
        for (int k = 0; k <= n; k++)
        {
            float total = 0.0f;
            for (int l = 1; l <= m; l++)
            {
                cds[l] = glm::distance(surf_points[k][l - 1], surf_points[k][l]);
                total += cds[l];
            }
            if (total == 0.0f)
                num -= 1;
            else
            {
                float d = 0.0f;
                for (int l = 1; l < m; l++)
                {
                    d += cds[l];
                    v_params[l] += d / total;
                }
            }
        }
        if (num == 0)
        {
            std::cout << "Error: Invalid surface mesh parameters." << std::endl;
            exit(1);
        }
        for (int l = 1; l < m; l++)
            v_params[l] /= num;
    }

    template <typename T>
    void GenerateKnotVector(std::vector<T> &knotvector, std::vector<T> &params, int p, int n)
    {
        int m = knotvector.size() - 1;
        for (int i = 0; i <= p; i++)
            knotvector[i] = 0.0f;
        for (int i = m - p; i <= m; i++)
            knotvector[i] = 1.0f;
        for (int j = 1; j <= n - p; j++)
        {
            float total = 0.0f;
            for (int i = j; i <= j + p - 1; i++)
            {
                total += params[i];
            }
            knotvector[j + p] = total / p;
        }
    }

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface0(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames)
    {
        int n = contour.control_points.size() - 1;
        int m = trace.control_points.size() - 1;
        tinynurbs::RationalSurface<T> surf;
        surf.degree_u = contour.degree;
        surf.degree_v = trace.degree;
        surf.knots_u = contour.knots;
        surf.knots_v = trace.knots;
        surf.control_points.resize(n + 1, m + 1, glm::vec<3, T>(0.0f));
        surf.weights.resize(n + 1, m + 1, 1.0f);
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= m; j++)
            {
                surf.control_points(i, j) = contour.control_points[i] + trace.control_points[j];
                surf.weights(i, j) = contour.weights[i] * trace.weights[j];
            }
        }
        return surf;
    }

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface1(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames)
    {
        int n = contour.control_points.size() - 1;
        int m = trace.control_points.size() - 1;

        v_bar.resize(m + 1);
        for(int j=0;j<=m;j++)
            v_bar[j] = float(j) / float(m);     
        frame::computeFrames3(frames, trace, v_bar);
        GenerateProfileCurves(profile_curves, contour, trace, shape_controller, v_bar, frames);

        std::vector<std::vector<glm::vec<3, T>>> surf_points(n + 1, std::vector<glm::vec<3, T>>(m + 1));
        float u_step = 1.0f / (n);
        float v_step = 1.0f / (m);
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= m; j++)
            {
                std::vector<glm::vec3> frame = frame::computeSingleFrame3(trace, v_step * j);               
                std::vector<glm::vec3> frame0 = frames[0];
                glm::vec<4, T> tmp_p = frame::getInverseFrameMatrix(frame) *
                                       shape_controller(v_step * j) *
                                       frame::getFrameMatrix(frame0) *
                                       glm::vec4(myBasis::curvePoint(contour, u_step * i), 1.0f);
                glm::vec<3, T> p = glm::vec<3, T>(tmp_p/tmp_p.w) + myBasis::curvePoint(trace, v_step * j);
                surf_points[i][j] = p;
            }

        }
        // 利用数据点插值曲面
        tinynurbs::RationalSurface<T> surf;
        surf.degree_u = contour.degree;
        surf.degree_v = trace.degree;
        // 确定参数值 9.6
        std::vector<T> u_params(n + 1), v_params(m + 1);
        SurfMeshParams(n, m, surf_points, u_params, v_params);
        // for (int i = 0; i < u_params.size(); i++)
        //     std::cout << u_params[i] << " ";
        // std::cout << std::endl;
        // for (int i = 0; i < v_params.size(); i++)
        //     std::cout << v_params[i] << " ";
        // std::cout << std::endl;

        // 确定节点向量 9.8
        std::vector<T> u_knotvector(n + 2 + surf.degree_u), v_knotvector(m + 2 + surf.degree_v);
        GenerateKnotVector(u_knotvector, u_params, surf.degree_u, n);
        GenerateKnotVector(v_knotvector, v_params, surf.degree_v, m);
        // for (int i = 0; i < u_knotvector.size(); i++)
        //     std::cout << u_knotvector[i] << " ";
        // std::cout << std::endl;
        // for (int i = 0; i < v_knotvector.size(); i++)
        //     std::cout << v_knotvector[i] << " ";
        // std::cout << std::endl;
        surf.knots_u = u_knotvector;
        surf.knots_v = v_knotvector;

        // 确定控制点 A9.1
        surf.control_points.resize(n + 1, m + 1, glm::vec<3, T>(0.0f));
        // std::vector<std::vector<T>> A(m + 1, std::vector<T>(m + 1, 0.0f));
        Eigen::MatrixXf A(n + 1, n + 1);
        for (int k = 0; k <= n; k++)
        {
            int span = myBasis::findSpan(u_knotvector, surf.degree_u, u_params[k]);
            std::vector<T> N = myBasis::bsplineBasis(surf.degree_u, span, u_knotvector, u_params[k]);
            // for(int i=0;i<N.size();i++)
            //     std::cout << N[i] << " ";
            // std::cout << std::endl;
            A.row(k) = Eigen::Map<const Eigen::VectorXf>(N.data(), N.size());
        }
        std::vector<std::vector<glm::vec<3, T>>> R(n + 1, std::vector<glm::vec<3, T>>(m + 1));
        for (int l = 0; l <= m; l++)
        {
            for (int i = 0; i < 3; i++)
            {
                std::vector<T> rhs(n + 1, 0.0f);
                for (int j = 0; j <= n; j++)
                    rhs[j] = surf_points[j][l][i];
                Eigen::Map<Eigen::VectorXf> b(rhs.data(), rhs.size());
                Eigen::VectorXf x = A.colPivHouseholderQr().solve(b);
                for (int j = 0; j <= n; j++)
                    R[j][l][i] = x[j];
            }
        }
        Eigen::MatrixXf C(m + 1, m + 1);
        for (int l = 0; l <= m; l++)
        {
            int span = myBasis::findSpan(v_knotvector, surf.degree_v, v_params[l]);
            std::vector<T> N = myBasis::bsplineBasis(surf.degree_v, span, v_knotvector, v_params[l]);
            C.row(l) = Eigen::Map<const Eigen::VectorXf>(N.data(), N.size());
        }
        for (int k = 0; k <= n; k++)
        {
            for (int i = 0; i < 3; i++)
            {
                std::vector<T> rhs(m + 1, 0.0f);
                for (int j = 0; j <= m; j++)
                    rhs[j] = R[k][j][i];
                Eigen::Map<Eigen::VectorXf> b(rhs.data(), rhs.size());
                Eigen::VectorXf x = C.colPivHouseholderQr().solve(b);
                for (int j = 0; j <= m; j++)
                    surf.control_points(k, j)[i] = x[j];
            }
        }
        // for(int i=0;i<=n;i++)
        // {
        //     for(int j=0;j<=m;j++)
        //     {
        //         std::cout << surf.control_points(i, j)[0] << " " << surf.control_points(i, j)[1] << " " << surf.control_points(i, j)[2] << std::endl;
        //     }
        // }

        // 计算控制点权重
        surf.weights.resize(n + 1, m + 1, 1.0f);
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= m; j++)
            {
                // glm::vec<3, T> p = surf_points[i][j];
                // glm::vec<3, T> c = surf.control_points(i, j);
                // T w = glm::dot(p - c, p - c);
                // surf.weights(i, j) = w;
                surf.weights(i, j) = contour.weights[i] * trace.weights[j];
            }
        }
        return surf;
    }

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface2(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames)
    {
        int n = contour.control_points.size() - 1;
        int m = trace.control_points.size() - 1;
        tinynurbs::RationalSurface<T> surf;
        surf.degree_u = contour.degree;
        surf.degree_v = trace.degree;
        surf.knots_u = contour.knots;
        surf.knots_v = trace.knots;
        surf.control_points.resize(n + 1, m + 1, glm::vec<3, T>(0.0f));
        surf.weights.resize(n + 1, m + 1, 1.0f);
        
        v_bar.resize(m + 1);
        for(int i=0;i<=m;i++)
            v_bar[i] = (float)i / m;
        frame::computeFrames3(frames, trace, v_bar);
        GenerateProfileCurves(profile_curves, contour, trace, shape_controller, v_bar, frames);
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= m; j++)
            {
                std::vector<glm::vec3> frame = frame::computeSingleFrame3(trace, (float)j / m);
                std::vector<glm::vec3> frame0 = frames[0];
                std::vector<glm::vec3> ts = myBasis::curveDerivatives(trace, 2, (float)j / m);
                std::vector<glm::vec3> frenet_frame = frame::computeSingleFrame1(trace, (float)j / m);
                glm::vec<3, T> V1 = glm::vec3(frame::getInverseFrameMatrix(frame) *
                                              shape_controller((float)j / m) *
                                              frame::getFrameMatrix(frame0) *
                                              glm::vec4(contour.control_points[i], 1.0f));
                // std::cout<<"V1: "<<V1[0]<<" "<<V1[1]<<" "<<V1[2]<<std::endl;
                T vn = glm::dot(frenet_frame[0], V1);
                T vb = glm::dot(frenet_frame[1], V1);
                T vt = glm::dot(frenet_frame[2], V1);
                // std::cout<<"vn: "<<vn<<" vb: "<<vb<<" vt: "<<vt<<std::endl;
                glm::vec<3, T> S = myBasis::curvePoint(trace, (float)j / m);
                float ks = glm::length(glm::cross(ts[2], ts[1])) / std::pow(glm::length(ts[1]), 3);
                // std::cout<<ks<<std::endl;       
                glm::vec<3, T> N_1 = frenet_frame[0] - ks * (trace.control_points[j] - S);
                // glm::vec<3, T> N_1 = frenet_frame[0];
                T lamda = glm::length(N_1);
                glm::vec<3, T> N_2 = N_1 / lamda;
                glm::vec<3, T> T_2 = glm::cross(N_2, frenet_frame[1]);

                glm::vec<3, T> E1 = trace.control_points[j] + lamda * (vn * N_2 + vt * T_2) + vb * frenet_frame[1];
                surf.control_points(i, j) = E1;
                
                surf.weights(i, j) = contour.weights[i] * trace.weights[j];
            }
        }
        return surf;
    }

    template <typename T>
    tinynurbs::RationalSurface<T> GenerateSweptSurface3(const tinynurbs::RationalCurve<T> &contour,
                                                        const tinynurbs::RationalCurve<T> &trace,
                                                        std::function<glm::mat4(T)> shape_controller,
                                                        std::vector<T> &v_bar,
                                                        std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                                                        std::vector<std::vector<glm::vec3>> &frames)
    {
        tinynurbs::RationalSurface<T> surf;

        int q = trace.degree;
        int ktv = trace.knots.size();
        int nsect = ktv - q - 1;
        //std::cout<<"nsect: "<<nsect<<std::endl;
        v_bar.resize(nsect);
        v_bar[0] = 0.0f;
        v_bar[nsect - 1] = 1.0f;
        for(int k=1;k<nsect-1;k++)
        {
            v_bar[k] = 0.0f;
            for(int i=k+1;i<=k+q;i++)
                v_bar[k] += trace.knots[i];
            v_bar[k] /= q;            
        }
        
        // frame::computeFrames2(frames, trace, v_bar);
        frame::computeFrames3(frames, trace, v_bar);
        GenerateProfileCurves(profile_curves, contour, trace, shape_controller, v_bar, frames);
        for(int k=0;k<profile_curves.size();k++)
        {
            for(int i=0;i<profile_curves[k].weights.size();i++)
                profile_curves[k].weights[i] *= trace.weights[k];
        }

        int n = contour.control_points.size() - 1;
        surf.degree_u = contour.degree;
        surf.degree_v = trace.degree;
        surf.knots_u = contour.knots;
        surf.knots_v = trace.knots;
        surf.control_points.resize(n + 1, nsect, glm::vec<3, T>(0.0f));
        surf.weights.resize(n + 1, nsect, 1.0f);
        Eigen::MatrixXf C(nsect, nsect);
        for (int l = 0; l < nsect; l++)
        {
            int span = myBasis::findSpan(trace.knots, surf.degree_v, v_bar[l]);
            std::vector<T> N = myBasis::bsplineBasis(surf.degree_v, span, trace.knots, v_bar[l]);
            C.row(l) = Eigen::Map<const Eigen::VectorXf>(N.data(), N.size());
        }
        for (int k = 0; k <= n; k++)
        {
            for (int i = 0; i < 4; i++)
            {
                if(i == 3)
                {
                    std::vector<T> rhs(nsect, 0.0f);
                    for (int j = 0; j < nsect; j++)
                        rhs[j] = profile_curves[j].weights[k];
                    Eigen::Map<Eigen::VectorXf> b(rhs.data(), rhs.size());
                    Eigen::VectorXf x = C.colPivHouseholderQr().solve(b);
                    for (int j = 0; j < nsect; j++)
                        surf.weights(k, j) = x[j];
                    
                }
                else
                {
                    std::vector<T> rhs(nsect, 0.0f);
                    for (int j = 0; j < nsect; j++)
                        rhs[j] = profile_curves[j].control_points[k][i];
                    Eigen::Map<Eigen::VectorXf> b(rhs.data(), rhs.size());
                    Eigen::VectorXf x = C.colPivHouseholderQr().solve(b);
                    for (int j = 0; j < nsect; j++)
                        surf.control_points(k, j)[i] = x[j];
                }

            }
        }
        return surf;

    }

    

    template <typename T>
    void GenerateProfileCurves(std::vector<tinynurbs::RationalCurve<T>> &profile_curves,
                               const tinynurbs::RationalCurve<T> &contour,
                               const tinynurbs::RationalCurve<T> &trace,
                               std::function<glm::mat4(T)> shape_controller,
                               const std::vector<T> &v_bar,
                               const std::vector<std::vector<glm::vec3>> &frames)
    {
        int n = contour.control_points.size() - 1;
        int m = trace.control_points.size() - 1;

        for(int i=0;i<v_bar.size();i++)
        {
            tinynurbs::RationalCurve<T> profile_curve;
            profile_curve.degree = contour.degree;
            profile_curve.knots = contour.knots;
            profile_curve.weights = contour.weights;
            profile_curve.control_points.resize(n + 1, glm::vec<3, T>(0.0f));
            for(int j=0;j<=n;j++)
            {
                profile_curve.control_points[j] = myBasis::curvePoint(trace, v_bar[i]) +
                                                  glm::vec<3, float>(frame::getInverseFrameMatrix(frames[i]) *
                                                                     shape_controller(v_bar[i]) *
                                                                     frame::getFrameMatrix(frames[0]) *
                                                                     glm::vec4(contour.control_points[j], 1.0f));
            }

            profile_curves.push_back(profile_curve);
        }

    }

    //显示模板化函数
    template void SurfMeshParams(int n, int m, std::vector<std::vector<glm::vec<3, float>>> &surf_points, std::vector<float> &u_params, std::vector<float> &v_params);

    template void GenerateKnotVector(std::vector<float> &knotvector, std::vector<float> &params, int p, int n);

    template tinynurbs::RationalSurface<float> GenerateSweptSurface0(const tinynurbs::RationalCurve<float> &contour,
                                                                     const tinynurbs::RationalCurve<float> &trace,
                                                                     std::function<glm::mat4(float)> shape_controller,
                                                                     std::vector<float> &v_bar,
                                                                     std::vector<tinynurbs::RationalCurve<float>> &profile_curves,
                                                                     std::vector<std::vector<glm::vec3>> &frames);

    template tinynurbs::RationalSurface<float> GenerateSweptSurface1(const tinynurbs::RationalCurve<float> &contour,
                                                                     const tinynurbs::RationalCurve<float> &trace,
                                                                     std::function<glm::mat4(float)> shape_controller,
                                                                     std::vector<float> &v_bar,
                                                                     std::vector<tinynurbs::RationalCurve<float>> &profile_curves,
                                                                     std::vector<std::vector<glm::vec3>> &frames);

    template tinynurbs::RationalSurface<float> GenerateSweptSurface2(const tinynurbs::RationalCurve<float> &contour,
                                                                     const tinynurbs::RationalCurve<float> &trace,
                                                                     std::function<glm::mat4(float)> shape_controller,
                                                                     std::vector<float> &v_bar,
                                                                     std::vector<tinynurbs::RationalCurve<float>> &profile_curves,
                                                                     std::vector<std::vector<glm::vec3>> &frames);

    template tinynurbs::RationalSurface<float> GenerateSweptSurface3(const tinynurbs::RationalCurve<float> &contour,
                                                                     const tinynurbs::RationalCurve<float> &trace,
                                                                     std::function<glm::mat4(float)> shape_controller,
                                                                     std::vector<float> &v_bar,
                                                                     std::vector<tinynurbs::RationalCurve<float>> &profile_curves,
                                                                     std::vector<std::vector<glm::vec3>> &frames);

    template void GenerateProfileCurves(std::vector<tinynurbs::RationalCurve<float>> &profile_curves,
                                        const tinynurbs::RationalCurve<float> &contour,
                                        const tinynurbs::RationalCurve<float> &trace,
                                        std::function<glm::mat4(float)> shape_controller,
                                        const std::vector<float> &v_bar,
                                        const std::vector<std::vector<glm::vec3>> &frames);
}