#include "myBasis.h"
#include <iostream>

namespace myBasis
{
    template <typename T>
    int findSpan(const std::vector<T> &knots, int degree, T u)
    {
		// index of last control point
		int n = static_cast<int>(knots.size()) - degree - 2;
		assert(n >= 0);
		/*
			// For u that is equal to last knot value
			if (util::close(u, knots[n + 1])) {
				return n;
			}
		*/
		// For values of u that lies outside the domain
		if(u > (knots[n + 1] - std::numeric_limits<T>::epsilon())) {
			return n;
		}
		if(u < (knots[degree] + std::numeric_limits<T>::epsilon())) {
			return degree;
		}

		int low = degree;
		int high = n + 1;
		int mid = (low + high) >> 1;
		while(low < high) {
			mid = (low + high) >> 1;
			if(knots[mid] >= u) {
				high = mid;
			} else
				low = mid + 1;
		}
		return low - 1;
    }

    template <typename T>
    std::vector<T> bsplineBasis(unsigned int deg, int span, const std::vector<T> &knots, T u)
    {
        std::vector<T> N;
        N.resize(deg + 1, T(0));
        std::vector<T> left, right;
        left.resize(deg + 1, static_cast<T>(0.0));
        right.resize(deg + 1, static_cast<T>(0.0));
        T saved = 0.0, temp = 0.0;

        N[0] = 1.0;

        for (int j = 1; j <= static_cast<int>(deg); j++)
        {
            left[j] = (u - knots[span + 1 - j]);
            right[j] = knots[span + j] - u;
            saved = 0.0;
            for (int r = 0; r < j; r++)
            {
                temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }

        std::vector<T> basis;
        basis.resize(knots.size() - deg - 1, T(0.0));
        for (int j = 0; j <= deg; j++)
        {
            basis[span - deg + j] = N[j];
        }
        return basis;
    }

    template <typename T>
    std::vector<T> noZeroBsplineBasis(unsigned int deg, int span, const std::vector<T> &knots, T u)
    {
        std::vector<T> N;
        N.resize(deg + 1, T(0));
        std::vector<T> left, right;
        left.resize(deg + 1, static_cast<T>(0.0));
        right.resize(deg + 1, static_cast<T>(0.0));
        T saved = 0.0, temp = 0.0;

        N[0] = 1.0;

        for (int j = 1; j <= static_cast<int>(deg); j++)
        {
            left[j] = (u - knots[span + 1 - j]);
            right[j] = knots[span + j] - u;
            saved = 0.0;
            for (int r = 0; r < j; r++)
            {
                temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }

        return N;
    }

    template <typename T>
    array2<T> bsplineDerBasis(unsigned int deg, int span, const std::vector<T> &knots, T u,
                              int num_ders)
    {
        std::vector<T> left, right;
        left.resize(deg + 1, 0.0);
        right.resize(deg + 1, 0.0);
        T saved = 0.0, temp = 0.0;

        array2<T> ndu(deg + 1, deg + 1);
        ndu(0, 0) = 1.0;

        for (int j = 1; j <= static_cast<int>(deg); j++)
        {
            left[j] = u - knots[span + 1 - j];
            right[j] = knots[span + j] - u;
            saved = 0.0;

            for (int r = 0; r < j; r++)
            {
                // Lower triangle
                ndu(j, r) = right[r + 1] + left[j - r];
                temp = ndu(r, j - 1) / ndu(j, r);
                // Upper triangle
                ndu(r, j) = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }

            ndu(j, j) = saved;
        }

        array2<T> ders(num_ders + 1, deg + 1, T(0));

        for (int j = 0; j <= static_cast<int>(deg); j++)
        {
            ders(0, j) = ndu(j, deg);
        }

        array2<T> a(2, deg + 1);

        for (int r = 0; r <= static_cast<int>(deg); r++)
        {
            int s1 = 0;
            int s2 = 1;
            a(0, 0) = 1.0;

            for (int k = 1; k <= num_ders; k++)
            {
                T d = 0.0;
                int rk = r - k;
                int pk = deg - k;
                int j1 = 0;
                int j2 = 0;

                if (r >= k)
                {
                    a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                    d = a(s2, 0) * ndu(rk, pk);
                }

                if (rk >= -1)
                {
                    j1 = 1;
                }
                else
                {
                    j1 = -rk;
                }

                if (r - 1 <= pk)
                {
                    j2 = k - 1;
                }
                else
                {
                    j2 = deg - r;
                }

                for (int j = j1; j <= j2; j++)
                {
                    a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                    d += a(s2, j) * ndu(rk + j, pk);
                }

                if (r <= pk)
                {
                    a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                    d += a(s2, k) * ndu(r, pk);
                }

                ders(k, r) = d;

                int temp = s1;
                s1 = s2;
                s2 = temp;
            }
        }

        T fac = static_cast<T>(deg);
        for (int k = 1; k <= num_ders; k++)
        {
            for (int j = 0; j <= static_cast<int>(deg); j++)
            {
                ders(k, j) *= fac;
            }
            fac *= static_cast<T>(deg - k);
        }

        return ders;
    }

    template <int dim, typename T>
    glm::vec<dim, T> curvePoint(unsigned int degree, const std::vector<T> &knots,
                                const std::vector<glm::vec<dim, T>> &control_points, T u)
    {
        glm::vec<dim, T> point(T(0));

        int span = myBasis::findSpan(knots, degree, u);
        std::vector<T> N = myBasis::noZeroBsplineBasis(degree, span, knots, u);

        for (unsigned int j = 0; j <= degree; j++)
        {
            point += static_cast<T>(N[j]) * control_points[span - degree + j];
        }
        return point;
    }

    template <int dim, typename T>
    std::vector<glm::vec<dim, T>> curveDerivatives(unsigned int degree, const std::vector<T> &knots,
                                                   const std::vector<glm::vec<dim, T>> &control_points,
                                                   int num_ders, T u)
    {

        typedef glm::vec<dim, T> tvecn;
        using std::vector;

        std::vector<glm::vec<dim, T>> curve_ders;
        curve_ders.resize(num_ders + 1);

        // Assign higher order derivatives to zero
        for (int k = degree + 1; k <= num_ders; k++)
        {
            curve_ders[k] = tvecn(0.0);
        }

        // Find the span and corresponding non-zero basis functions & derivatives
        int span = myBasis::findSpan(knots, degree, u);
        array2<T> ders = myBasis::bsplineDerBasis<T>(degree, span, knots, u, num_ders);

        // Compute first num_ders derivatives
        int du = num_ders < static_cast<int>(degree) ? num_ders : static_cast<int>(degree);
        for (int k = 0; k <= du; k++)
        {
            curve_ders[k] = tvecn(0.0);
            for (int j = 0; j <= static_cast<int>(degree); j++)
            {
                curve_ders[k] += static_cast<T>(ders(k, j)) * control_points[span - degree + j];
            }
        }
        return curve_ders;
    }

    template <typename T>
    glm::vec<3, T> curvePoint(const RationalCurve<T> &crv, T u)
    {

        typedef glm::vec<4, T> tvecnp1;

        // Compute homogenous coordinates of control points
        std::vector<tvecnp1> Cw;
        Cw.reserve(crv.control_points.size());
        for (size_t i = 0; i < crv.control_points.size(); i++)
        {
            Cw.push_back(tvecnp1(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i])));
        }

        // Compute point using homogenous coordinates
        tvecnp1 pointw = myBasis::curvePoint(crv.degree, crv.knots, Cw, u);

        // Convert back to cartesian coordinates
        return util::homogenousToCartesian(pointw);
    }

    template <typename T>
    std::vector<glm::vec<3, T>> curveDerivatives(const RationalCurve<T> &crv, int num_ders, T u)
    {

        typedef glm::vec<3, T> tvecn;
        typedef glm::vec<4, T> tvecnp1;

        std::vector<tvecn> curve_ders;
        curve_ders.reserve(num_ders + 1);

        // Compute homogenous coordinates of control points
        std::vector<tvecnp1> Cw;
        Cw.reserve(crv.control_points.size());
        for (size_t i = 0; i < crv.control_points.size(); i++)
        {
            Cw.push_back(util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
        }

        // Derivatives of Cw
        std::vector<tvecnp1> Cwders =
            myBasis::curveDerivatives(crv.degree, crv.knots, Cw, num_ders, u);

        // Split Cwders into coordinates and weights
        std::vector<tvecn> Aders;
        std::vector<T> wders;
        for (const auto &val : Cwders)
        {
            Aders.push_back(util::truncateHomogenous(val));
            wders.push_back(val.w);
        }

        // Compute rational derivatives
        for (int k = 0; k <= num_ders; k++)
        {
            tvecn v = Aders[k];
            for (int i = 1; i <= k; i++)
            {
                v -= static_cast<T>(util::binomial(k, i)) * wders[i] * curve_ders[k - i];
            }
            curve_ders.push_back(v / wders[0]);
        }
        return curve_ders;
    }

    template <typename T>
    glm::vec<3, T> curveTangent(const RationalCurve<T> &crv, T u)
    {
        std::vector<glm::vec<3, T>> ders = myBasis::curveDerivatives(crv, 1, u);
        glm::vec<3, T> du = ders[1];
        T du_len = glm::length(du);
        if (!util::close(du_len, T(0)))
        {
            du /= du_len;
        }
        return du;
    }


    //显示模板化函数
    template int findSpan(const std::vector<float> &knotvector, int degree, float u);
    template std::vector<float> bsplineBasis(unsigned int deg, int span, const std::vector<float> &knots, float u);
    template std::vector<float> noZeroBsplineBasis(unsigned int deg, int span, const std::vector<float> &knots, float u);
    template array2<float> bsplineDerBasis(unsigned int deg, int span, const std::vector<float> &knots, float u, int num_ders);
    template glm::vec<3, float> curvePoint(unsigned int degree, const std::vector<float> &knots, const std::vector<glm::vec<3, float>> &control_points, float u);
    template glm::vec<4, float> curvePoint(unsigned int degree, const std::vector<float> &knots, const std::vector<glm::vec<4, float>> &control_points, float u);
    template std::vector<glm::vec<3, float>> curveDerivatives(unsigned int degree, const std::vector<float> &knots, const std::vector<glm::vec<3, float>> &control_points, int num_ders, float u);
    template std::vector<glm::vec<4, float>> curveDerivatives(unsigned int degree, const std::vector<float> &knots, const std::vector<glm::vec<4, float>> &control_points, int num_ders, float u);
    template glm::vec<3, float> curveTangent(const RationalCurve<float> &crv, float u);
    template glm::vec<3, float> curvePoint(const RationalCurve<float> &crv, float u);
    template std::vector<glm::vec<3, float>> curveDerivatives(const RationalCurve<float> &crv, int num_ders, float u);

}