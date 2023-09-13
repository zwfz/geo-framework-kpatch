
#include <fstream>
#include <numeric>
#include <string>

#include "object.hh"
#include <transfinite/surface-corner-based.hh>
#include <transfinite/surface-composite-ribbon.hh>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/StdVector>

using namespace Transfinite;
using namespace Eigen;

class KPatch : public Object {
public:
	KPatch(std::string filename);
	virtual ~KPatch();
	virtual void draw(const Visualization& vis) const override;
	virtual void drawWithNames(const Visualization& vis) const override;
	virtual Vector postSelection(int selected) override;
	virtual void movement(int selected, const Vector& pos) override;
	virtual void updateBaseMesh() override;
	virtual bool reload() override;
private:
	size_t degree[2];
	std::vector<Vector> control_points;
	CurveVector curves;
	PointVector PVs;

};