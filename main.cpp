#include <memory>
#include <iostream>
#include <ogl_materials.h>
#include <gl3d_model.h>
#include <svector.h>
#include <color_map.h>

#define M_PI 3.14159265358
// --------------------- Experiment ---------------------

struct Coil {
    Coil(svector::float4 pos, float angle, int nloops, float radius, float current, float loop_gap)
        : pos(pos), angle(angle), n(nloops), r(radius), j(current), d(loop_gap)
    {}
    svector::float4 pos;  // vector position of the center
    float angle = 0;      // magnet angle on about the y axis
    int     n = 0;        // number of loops
    float   r = 0;        // radius meters
    float   j = 0;        // current
    float   d = 0;        // Distance between loops in meters
};

struct Plate {
    Plate(svector::float4 pos, int axis, float width, float height, float depth)
        : pos(pos), axis(axis), w(width), h(height), d(depth)
    {}
    svector::float4 pos;  // vector position of the center
    int     axis = 0;     // x = 0, y = 1, z = 2. Axis normal to the plate plane
    float   w = 0;        // width meters
    float   h = 0;        // height meters
    float   d = 0;        // depth meters
};

struct Field {
    Field(const svector::float4& pos, int axis, int cols, int rows, float width, float height)
        : pos(pos), axis(axis), cols(cols), rows(rows), width(width), height(height), field(size_t(rows) * cols) {}
    svector::float4 pos;  // position of the center of the field
    int axis;             // x = 0, y = 1, z = 2. Axis normal to the field plane
    int rows;
    int cols;
    float width;          // meters
    float height;         // meters
    std::vector<svector::float4> field;
};

using namespace svector;
class Experiment {
public:    
    Experiment()
        : m_coil(svector::float4(0, 0, 0), 0, 200, 0.025f, 1, 0.040f / (200 - 1))
        , m_plate(svector::float4(0.050f, 0, 0), 0, 0.150f, 0.150f, 0.003f)
        , m_plate_field(svector::float4(0.050f, 0, 0), 0, 400, 400, 0.150f, 0.150f)
        , m_coil_field(svector::float4(-0.050f, 0, 0), 1, 400, 400, 0.200f, 0.200f)
    {}

    void calculateMagneticField(Field& field, Coil& l) {
        std::cout << "starting calculation" << std::endl;
        float offsetx = field.width / 2;  // center of the field x
        float offsety = field.height / 2; // center of the field y
        float cx = field.width / field.cols;
        float cy = field.height / field.rows;
        // U0 / 4 pi * j constant
        // U0 = 4 pi x 10-7
        float k = 1E-7 * l.j;
        float alpha = -l.angle / 180 * M_PI; // counter clockwise
        float cosa = cos(alpha);
        float sina = sin(alpha);
        float h = (m_coil.n - 1) * m_coil.d;

        // Number of segments in the integration of the ring
        int N = 100;
        float ds = 2 * M_PI * l.r / N;
        int x, y;
        std::vector<svector::float4>& f = field.field;
        // Iterate through every point in space
#pragma omp parallel for schedule(static) private(x, y) shared(f)
        for (y = 0; y < field.rows; ++y) {
            for (x = 0; x < field.cols; ++x) {
                // p0 is the point on the grid for which the magnetic field 
                // due to every coild will be calculated to below
                // multiplying by d converts the x and y values to mm 
                // subtracting s moves it to the center of the grid
                float4 p0;
                if (field.axis == 0) p0 = float4(field.pos.x, field.pos.y + x * cx - offsetx, field.pos.z + y * cy - offsety);
                if (field.axis == 1) p0 = float4(field.pos.x + x * cx - offsetx, field.pos.y, field.pos.z + y * cy - offsety);

                // The integral value of the magnetic field at p0
                float4 sum = 0;

                // Loop through every loop in the coil
                for (int j = 0; j < l.n; ++j) {
                    // Loop through all segments of the loop
                    for (int i = 0; i < N; ++i) {
                        //point on circle in x, y, z
                        // The x component depends on j since each loop moves to the left with l.s
                        float theta = float(i) * 2.0f * M_PI / N;
                        float4 p1(h / 2 - j * l.d, l.r * sin(theta), l.r * cos(theta));
                        // rotate p1 along the y axis
                        float x1 = p1.x * cosa - p1.z * sina;
                        float z1 = p1.x * sina + p1.z * cosa;
                        p1.x = x1 - h / 2;
                        p1.z = z1;
                        p1 += l.pos;

                        // find unit vector tangent to point on circle in the direction of the current
                        // given by <0, -cos(theta), sin(theta)>
                        float4 vt(0, -cos(theta), sin(theta));
                        // rotate vt along the y axis
                        x1 = vt.x * cosa - vt.z * sina;
                        z1 = vt.x * sina + vt.z * cosa;
                        vt.x = x1;
                        vt.z = z1;

                        // find r vector between the two points
                        float4 vr = p0 - p1;

                        // find r distance between them 
                        float r = vr.norm();

                        // find cross product ds x r and plug into db equation
                        float4 db = k * cross3d(vt, vr) / (r * r * r) * ds;

                        sum += db;
                    }
                }
                // Assign calculated integral magnetic field value to segment
                f[size_t(x) + size_t(y) * field.cols] = sum;
            }
        }
        std::cout << "ended calculation" << std::endl;
    }

    float4 calculateMagneticField(float4& p0) {
        // U0 / 4 pi * j constant
        // U0 = 4 pi x 10-7
        float k = 1E-7 * m_coil.j;
        float alpha = -m_coil.angle / 180 * M_PI; // counter clockwise
        float cosa = cos(alpha);
        float sina = sin(alpha);
        float h = (m_coil.n - 1) * m_coil.d;

        // Number of segments in the integration of the ring
        int N = 100;
        float ds = 2 * M_PI * m_coil.r / N;
        int x, y;
        // The integral value of the magnetic field at p0
        float4 sum = 0;

        // Loop through every loop in the coil
        for (int j = 0; j < m_coil.n; ++j) {
            // Loop through all segments of the loop
            for (int i = 0; i < N; ++i) {
                //point on circle in x, y, z
                // The x component depends on j since each loop moves to the left with l.s
                float theta = float(i) * 2.0f * M_PI / N;
                // increaset p1.x by h/2 for rotation about the center
                float4 p1(h / 2 - j * m_coil.d, m_coil.r * sin(theta), m_coil.r * cos(theta));
                // rotate p1 along the y axis
                float x1 = p1.x * cosa - p1.z * sina;
                float z1 = p1.x * sina + p1.z * cosa;
                // substract h/2 to get back to the right position
                p1.x = x1 - h / 2;
                p1.z = z1;
                // add coil position after rotation, otherwise it will rotate the whole thing
                p1 += m_coil.pos;

                // find unit vector tangent to point on circle
                // given by <0, cos(theta), - sin(theta)>
                float4 vt(0, -cos(theta), sin(theta));
                // rotate vt along the y axis
                x1 = vt.x * cosa - vt.z * sina;
                z1 = vt.x * sina + vt.z * cosa;
                vt.x = x1;
                vt.z = z1;

                // find r vector between the two points
                float4 vr = p0 - p1;

                // find r distance between them 
                float r = vr.norm();

                // find cross product ds x r and plug into db equation
                float4 db = k * cross3d(vt, vr) / (r * r * r) * ds;
                sum += db;
            }
        }
        return sum;
    }

    void computePlateField() {
        calculateMagneticField(m_plate_field, m_coil);
    }
    void computeCoilField() {
        calculateMagneticField(m_coil_field, m_coil);
    }

    Coil m_coil;
    Plate m_plate;
    Field m_plate_field;
    Field m_coil_field;
};

std::shared_ptr<Experiment> experiment = std::make_shared<Experiment>();


// ---------------------- Graphics --------------------------------------------
class GLModel : public GL3DModel {
public:
    GLModel(std::shared_ptr<Experiment> experiment, float scale, bool drawAxes = false, bool enableLighting = false)
    : GL3DModel(scale, drawAxes, enableLighting), m_experiment(experiment), m_coil(experiment->m_coil), m_plate(experiment->m_plate)
    {}

    void drawMagnet() {
        glPushMatrix();
        float h = (m_coil.n - 1) * m_coil.d * 1000;
        glTranslatef(m_coil.pos.x * 1000 - h / 2, m_coil.pos.y * 1000, m_coil.pos.z * 1000);
        glRotatef(m_coil.angle + 90, 0, 1, 0);
        
        chrome.active();
        glColor3f(1, 1, 1);
        drawCylinder(h, m_coil.r * 1000);
        glPopMatrix();
    }

    void drawPlate() {
        glPushMatrix();
        glTranslatef(m_plate.pos.x * 1000+ m_plate.d * 1000 / 2 + 0.15, m_plate.pos.y * 1000, m_plate.pos.z * 1000);
        glRotatef(-90, 0, 1, 0);
        aluminum.active();
        glColor3f(1, 1, 1);
        drawCube(m_plate.w * 1000, m_plate.h * 1000, m_plate.d * 1000);
        glPopMatrix();
    }

    void drawFieldPlane(Field& f) {
        glPushMatrix();
        glDisable(GL_LIGHTING);
        glTranslatef(f.pos.x * 1000, f.pos.y * 1000, f.pos.z * 1000);
        if (f.axis == 0) { glRotatef(-90, 0, 0, 1); glRotatef(90, 1, 0, 0);}
        if (f.axis == 1) glRotatef(90, 1, 0, 0);
        float sx = f.width / f.cols * 1000;
        float sy = f.height / f.rows * 1000;
        float offsetx = f.width / 2 * 1000;
        float offsety = f.height / 2 * 1000;

        glColor3f(0.5, 0.5, 0.5);
        // itterates through every point on the grid, distance d away from each other
        for (int y = 0; y < f.rows; ++y) {
            for (int x = 0; x < f.cols; ++x) {
                float norm = f.field[size_t(x) + size_t(y) * f.cols].norm();
                float m = log(norm + 1) * 400;
                m = m < 0 ? 0 : m > 1 ? 1 : m;
                svector::float4 color = heat_map(m);
                glColor3f(color.r, color.g, color.b);
                glBegin(GL_QUADS);
                glVertex3f(x * sx - offsetx      , y * sy - offsety, 0);
                glVertex3f((x + 1) * sx - offsetx, y * sy - offsety, 0);
                glVertex3f((x + 1) * sx - offsetx, (y + 1) * sy - offsety, 0);
                glVertex3f(x * sx - offsetx      , (y + 1) * sy - offsety, 0);
                glEnd();
            }
        }
        glEnable(GL_LIGHTING);
        glPopMatrix();
    }

    void render3D(float scale) {
        drawMagnet();
        drawPlate();
        drawFieldPlane(m_experiment->m_plate_field);
        drawFieldPlane(m_experiment->m_coil_field);
    }
private:
    std::shared_ptr<Experiment> m_experiment;
    Coil m_coil;
    Plate m_plate;
};

GLModel glExperiment(experiment, 1, true, true);

void mouse_wheel(int wheel, int direction, int x, int y) { glExperiment.mouse_wheel(wheel, direction, x, y); }
void mouse_button(int button, int state, int x, int y) { glExperiment.mouse_button(button, state, x, y); }
void mouse_active_motion(int x, int y) { glExperiment.mouse_active_motion(x, y); }
void mouse_passive_motion(int x, int y) { glExperiment.mouse_passive_motion(x, y); }

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int width = viewport[2];
    int height = viewport[3];
    glExperiment.render(width, height);
    glutSwapBuffers();
}

void reshape(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    display();
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 27: exit(0);
        break;
    case 32: {
        experiment->computeCoilField();
        experiment->computePlateField();
        break;
    }
    case 'a': {
        break;
    }
    case 'd': {
        break;
    }
    }
    glutPostRedisplay();
}

//------------------------------------ Main ------------------------------------

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowPosition(500, 300);
    glutInitWindowSize(1000, 1000);
    glutCreateWindow("Simulation");

    glutDisplayFunc(display);
    glutIdleFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    // This function is only supported by freeglut
    // Linux will be emulated by implementation of GLModel3D
    // MacOS remove this and the implementation will emulate mouse wheel
    // Not tested on MacOS/Linux and I don't even know how the trackpad will see this events
    glutMouseWheelFunc(mouse_wheel);

    glutMouseFunc(mouse_button);
    glutMotionFunc(mouse_active_motion);
    glutPassiveMotionFunc(mouse_passive_motion);
    glutMouseWheelFunc(mouse_wheel);

    glClearColor(0, 0, 0, 1);

    glutMainLoop();
    return 0;
}