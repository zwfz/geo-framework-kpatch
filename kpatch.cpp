#include <fstream>

#include "KPatch.h"
#include "KappaCurve.h"

KPatch::KPatch(std::string filename) : Object(filename) {
	reload();
}

KPatch::~KPatch() {
}

void KPatch::draw(const Visualization& vis) const {
    Object::draw(vis);
    if (vis.show_control_points) {
        glDisable(GL_LIGHTING);
        glLineWidth(3.0);
        glColor3d(0.3, 0.3, 1.0);

        //Draw control points lines
        size_t index = 0;
        for (size_t i = 0, index = 0; i <= degree[0]; ++i) {
            glBegin(GL_LINE_STRIP);
            for (size_t j = 0; j <= degree[1]; ++j, ++index) {
                const auto& p = control_points[index];
                glVertex3dv(p.data());
            }
            glEnd();
            
        }

        // Draw curves
        glColor3d(1, 0.753, 0.796);
        for (size_t i = 0; i < 4; ++i) {
            glBegin(GL_LINE_STRIP);
            for (size_t j = 0; j <= 100; ++j)
            {
                const auto& p = curves[i]->eval(j / 100.0);
                glVertex3dv(p.data());
            }  
            glEnd();
        }

        //Draw control points
        glLineWidth(1.0);
        glPointSize(8.0);
        glColor3d(1.0, 0.0, 1.0);
        glBegin(GL_POINTS);
        for (const auto& p : control_points)
            glVertex3dv(p.data());
        glEnd();
        glPointSize(1.0);
        glEnable(GL_LIGHTING);
    }


}

void KPatch::drawWithNames(const Visualization& vis) const {
    if (!vis.show_control_points)
        return;
    for (size_t i = 0; i < control_points.size(); ++i) {
        const auto& p = control_points[i];
        glPushName(i);
        glRasterPos3dv(p.data());
        glPopName();
    }
}

Vector KPatch::postSelection(int selected) {
    return control_points[selected];
}

void KPatch::movement(int selected, const Vector& pos) {
    control_points[selected] = pos;
}

static void bernstein(size_t n, double u, std::vector<double>& coeff) {
    coeff.clear(); coeff.reserve(n + 1);
    coeff.push_back(1.0);
    double u1 = 1.0 - u;
    for (size_t j = 1; j <= n; ++j) {
        double saved = 0.0;
        for (size_t k = 0; k < j; ++k) {
            double tmp = coeff[k];
            coeff[k] = saved + tmp * u1;
            saved = tmp * u;
        }
        coeff.push_back(saved);
    }
}

void KPatch::updateBaseMesh() {
    
    size_t index = 0;
    for (size_t i = 0, index = 0; i <= degree[0]; ++i) {
        PVs.clear();
        for (size_t j = 0; j <= degree[1]; ++j, ++index) {
            PVs.push_back({ control_points[index][0], control_points[index][1], control_points[index][2] });
        }
        curves.push_back(std::make_shared<KappaCurve>( PVs, 0, 1));
    }
    size_t resolution = 50;

    //SurfaceCornerBased surf;
    SurfaceCompositeRibbon surf;
    surf.setCurves(curves);
    surf.setupLoop();
    surf.update();
    PointVector points;
    std::list<std::array<size_t, 3>> triangles;
    TriMesh mesh1;
    mesh1 = surf.eval(resolution);
    points= mesh1.points();
    triangles= mesh1.triangles();

    mesh.clear();
    std::vector<BaseMesh::VertexHandle> handles, tri;
    size_t n = degree[0], m = degree[1];

    for (const auto& p : points)
    {
        Vector pp(0.0, 0.0, 0.0);
        pp[0] = p[0];
        pp[1] = p[1];
        pp[2] = p[2];
        handles.push_back(mesh.add_vertex(pp));
    }

   
    for (const auto& t : triangles)
    {
        tri.clear();
        tri.push_back(handles[t[0]]);
        tri.push_back(handles[t[1]]);
        tri.push_back(handles[t[2]]);
        mesh.add_face(tri);
    }

    Object::updateBaseMesh(false, false);
}

bool KPatch::reload() {
    
    size_t n, m;
    try {
        std::ifstream f(filename.c_str());
        f.exceptions(std::ios::failbit | std::ios::badbit);
        f >> n >> m;
        degree[0] = n++; degree[1] = m++;
        control_points.resize(n * m);
        for (size_t i = 0, index = 0; i < n; ++i)
            for (size_t j = 0; j < m; ++j, ++index)
                f >> control_points[index][0] >> control_points[index][1] >> control_points[index][2];
    }
    catch (std::ifstream::failure&) {
        return false;
    }
    updateBaseMesh();
    return true;
}