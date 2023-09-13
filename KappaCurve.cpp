
#include "KappaCurve.h"

KappaCurve::KappaCurve(const PointVector& cpts, double u1, double u2) : u1(u1), u2(u2) {
    cps = cpts;
    NumCps = (int)cps.size() - 2;
    lambdai.resize(NumCps, defaultValue);
    ti.resize(NumCps, defaultValue);
    ci0.resize(NumCps);
    ci1.resize(NumCps);
    ci2.resize(NumCps);
    ci0[0] = cps[0];
    ci2[NumCps - 1] = cps[NumCps + 1];

    for (int i = 0; i < NumCps; ++i)
    {
        ci1[i] = cps[i + 1];
    }

    Calculate_ci02();

    for (int iter = 0; iter < ITER_NUM; ++iter)
    {
        PointVector prevci1;
        prevci1.resize(NumCps);
        prevci1 = ci1;
        Calculate_ci1();
        Calculate_ci02();
        Calculate_ti();
        Calculate_lambdai();
        if (CheckCvrg() || CheckItr(prevci1, NumCps))
        {
            //std::cout << "converge" << std::endl;
            break;
        }
    }

    C1 = BSCurve({ ci0[0], ci1[0], ci2[0] });
    C2 = BSCurve({ ci0[1], ci1[1], ci2[1] });
    l1 = C1.arcLength(u1, 1);
    l2 = C2.arcLength(0, u2);
    ratio = l1 / (l1 + l2);
}

KappaCurve::~KappaCurve() {
}

Point3D KappaCurve::eval(double u) const {
    auto [first, t] = parameter(u);
    if (first)
        return C1.eval(t);
    return C2.eval(t);
}

Point3D KappaCurve::eval(double u, size_t nr_der, VectorVector& der) const {
    auto [first, t] = parameter(u);
    if (first)
        return C1.eval(t, nr_der, der);
    return C2.eval(t, nr_der, der);
}

void KappaCurve::reverse() {
    C1.reverse(); C1.normalize();
    C2.reverse(); C2.normalize();
    std::swap(C1, C2);
    std::swap(u1, u2);
    u1 = 1 - u1;
    u2 = 1 - u2;
    std::swap(l1, l2);
    ratio = 1 - ratio;
}

double KappaCurve::arcLength(double from, double to) const {
    if (from >= to)
        return 0.0;
    auto [first1, t1] = parameter(from);
    auto [first2, t2] = parameter(to);
    if (first2)
        return C1.arcLength(t1, t2);
    if (!first1)
        return C2.arcLength(t1, t2);
    return C1.arcLength(t1, 1) + C2.arcLength(0, t2);
}

void KappaCurve::Calculate_ci02()
{
    for (int i = 0; i < (NumCps - 1); ++i)
    {
        ci2[i] = ci0[i + 1] = ci1[i].operator*(1 - lambdai[i]) + ci1[i + 1].operator*(lambdai[i]);
    }
}

void KappaCurve::Calculate_ci1()
{
    MatrixXd A(NumCps, NumCps);
    VectorXd b1(NumCps), b2(NumCps), b3(NumCps);
    A.setZero();
    if (NumCps == 1)
    {
        ci1[0] = (cps[1] - cps[0].operator*(pow(1 - ti[0], 2)) - cps[2].operator*(pow(ti[0], 2))) / (2 * ti[0] * (1 - ti[0]));
        return;
    }
    for (int i = 0; i < NumCps; ++i)
    {
        if (i == 0)
        {
            A(0, 0) = 2.0 * ti[0] * (1.0 - ti[0]) + ti[0] * ti[0] * (1.0 - lambdai[0]);
            A(0, 1) = ti[0] * ti[0] * lambdai[0];
            b1[0] = cps[1][0] - pow(1.0 - ti[0], 2) * cps[0][0];
            b2[0] = cps[1][1] - pow(1.0 - ti[0], 2) * cps[0][1];
            b3[0] = cps[1][2] - pow(1.0 - ti[0], 2) * cps[0][2];
        }
        else if (i == NumCps - 1)
        {
            A(NumCps - 1, NumCps - 2) = pow(1.0 - ti[NumCps - 1], 2) * (1.0 - lambdai[NumCps - 2]);
            A(NumCps - 1, NumCps - 1) = pow(1.0 - ti[NumCps - 1], 2) * lambdai[NumCps - 2] + 2.0 * ti[NumCps - 1] * (1.0 - ti[NumCps - 1]);
            b1[NumCps - 1] = cps[NumCps][0] - ti[NumCps - 1] * ti[NumCps - 1] * cps[NumCps + 1][0];
            b2[NumCps - 1] = cps[NumCps][1] - ti[NumCps - 1] * ti[NumCps - 1] * cps[NumCps + 1][1];
            b3[NumCps - 1] = cps[NumCps][2] - ti[NumCps - 1] * ti[NumCps - 1] * cps[NumCps + 1][2];
        }
        else
        {
            A(i, i - 1) = (1.0 - lambdai[i - 1]) * pow(1.0 - ti[i], 2);
            A(i, i) = lambdai[i - 1] * pow(1.0 - ti[i], 2) + (2.0 - (1.0 + lambdai[i]) * ti[i]) * ti[i];
            A(i, i + 1) = lambdai[i] * ti[i] * ti[i];
            b1[i] = cps[i + 1][0];
            b2[i] = cps[i + 1][1];
            b3[i] = cps[i + 1][2];
        }
    }
    VectorXd x1 = A.lu().solve(b1);
    VectorXd x2 = A.lu().solve(b2);
    VectorXd x3 = A.lu().solve(b3);

    for (int i = 0; i < NumCps; ++i)
    {
        ci1[i][0] = x1[i];
        ci1[i][1] = x2[i];
        ci1[i][2] = x3[i];
    }
}

void KappaCurve::Calculate_ti()
{
    for (int i = 0; i < NumCps; ++i)
    {
        Point3D ci2_ci0 = ci2[i] - ci0[i];
        Point3D ci0_piplus1 = ci0[i] - cps[i + 1];
        double a = ci2_ci0.normSqr();
        double b = 3 * ci2_ci0.operator*(ci0_piplus1);
        double c = (ci0[i].operator*(3) - cps[i + 1].operator*(2) - ci2[i]).operator*(ci0_piplus1);
        double d = -ci0_piplus1.normSqr();
        ti[i] = (a == 0) ? 0.5 : getRealSolutionOfCubicFunc(b / a, c / a, d / a);
    }
}

double KappaCurve::getRealSolutionOfCubicFunc(const double a, const double b, const double c)
{
    const double p = b - a * a / 3;
    const double q = 2 * a * a * a / 27 - a * b / 3 + c;
    const double tmp = q * q / 4 + p * p * p / 27;
    if (tmp < 0) return 0;
    const double mTmp = -q / 2 + sqrt(tmp);
    const double nTmp = -q / 2 - sqrt(tmp);
    const double m = (mTmp < 0) ? -pow(-mTmp, 1 / 3.0) : pow(mTmp, 1 / 3.0);
    const double n = (nTmp < 0) ? -pow(-nTmp, 1 / 3.0) : pow(nTmp, 1 / 3.0);
    return -a / 3 + m + n;
}

void KappaCurve::Calculate_lambdai()
{
    for (int i = 0; i < NumCps - 1; i++)
    {
        double A = sqrt(TrgArea(ci0[i], ci1[i], ci1[i + 1]));
        double B = sqrt(TrgArea(ci1[i], ci1[i + 1], ci2[i + 1]));
        lambdai[i] = A / (A + B);

        if (lambdai[i] < 0.0)
        {
            lambdai[i] = 0.0;
        }
        else if (lambdai[i] > 1.0)
        {
            lambdai[i] = 1.0;
        }
    }
}

bool KappaCurve::CheckCvrg()
{
    Point3D cti;

    double dis0;
    double dis1 = 0;

    for (int i = 0; i < NumCps; ++i)
    {
        cti = ci0[i].operator*((1 - ti[i]) * (1 - ti[i])) + ci1[i].operator*(2 * (1 - ti[i]) * ti[i]) + ci2[i].operator*(ti[i] * ti[i]);

        dis0 = (cti - cps[i + 1]).norm();

        dis1 += dis0;
    }

    //cout << "dis1:" << dis1 << endl;

    if (dis1 > 1.0e-8)
        return false;

    return true;
}

bool KappaCurve::CheckItr(const PointVector& prevci1, int n)
{
    double s = 0;
    for (int i = 0; i < n; ++i)
    {
        s += (ci1[i] - prevci1[i]).norm();
    }
    //cout << "s=" << s << endl;
    return s < 1.e-4;
}

double KappaCurve::TrgArea(const Point3D& p0, const Point3D& p1, const Point3D& p2)
{
    Point3D v1 = p1 - p0;
    Point3D v2 = p2 - p0;
    Point3D v3 = v1.operator^(v2);
    double Area = 0.5 * v3.norm();
    return Area;
}

std::pair<bool, double> KappaCurve::parameter(double u) const {
    if (u < ratio)
        return { true, u / ratio * (1 - u1) + u1 };
    return { false, (u - ratio) / (1 - ratio) * u2 };
}