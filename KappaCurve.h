#include <fstream>
#include <numeric>
#include <string>

#include <transfinite/surface-corner-based.hh>
#include <transfinite/surface-composite-ribbon.hh>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/StdVector>

using namespace Transfinite;
using namespace Eigen;

class KappaCurve : public Curve {

public:
	KappaCurve(const PointVector& cpts, double u1, double u2);
	virtual ~KappaCurve();
	Point3D eval(double u) const override;
	Point3D eval(double u, size_t nr_der, VectorVector& der) const override;
	void reverse() override;
	double arcLength(double from, double to) const override;
	void Calculate_ci02();
	void Calculate_ci1();
	void Calculate_ti();
	double getRealSolutionOfCubicFunc(const double a, const double b, const double c);
	void Calculate_lambdai();
	double TrgArea(const Point3D& p0, const Point3D& p1, const Point3D& p2);
	bool CheckCvrg();
	bool CheckItr(const PointVector& prevci1, int n);

private:
	std::pair<bool, double> parameter(double u) const;
	BSCurve C1, C2;
	double u1, u2, l1, l2, ratio;
	size_t degree[2];
	PointVector ci0, ci1, ci2;
	PointVector cps;
	int NumCps;
	DoubleVector lambdai, ti;
	double defaultValue = 0.5;
	int ITER_NUM = 400;
};
