#pragma once

#include <GL/freeglut.h>
#include <GL/glut.h>
#include <svector.h>
#include <ogl_utilities.h>

// Most stupid hack to use mouse wheel on different platforms
#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP    3
#  define GLUT_WHEEL_DOWN  4
#  define GLUT_WHEEL_LEFT  5
#  define GLUT_WHEEL_RIGHT 6
#endif

#define HIGH_DPI_MOUSE    1       // enable if mouse is high dpi
#define MOUSE_ZOOM        2.0     // mouse zoom increment/decrement
#define MOUSE_ZOOM_NEAR   400     // mouse minimum near zoom 
#define MOUSE_ZOOM_FAR    -200    // mouse maximum far zoom
#define MOUSE_FACTOR      8000.0f // mouse_dx mouse_dy rotation factor DPI dependent 
// lights
static bool ShowLights = true;    // when enablelighting show ligths position for debugging
static bool Light0Enabled = true;
static bool Light1Enabled = true;
static bool Light2Enabled = true;
const svector::float4 Light0(0, 40, 70);
const svector::float4 Light1(-100, 100, 70);
const svector::float4 Light2(0, -100, 0);

class GL3DModel {
public:
    GL3DModel(float scale, bool drawAxes = false, bool enableLighting = false)
        : m_camera(0, 0, 0, 1), m_mouse_vz(-scale), m_scale(scale), m_drawAxes(drawAxes)
        , m_enableLighting(enableLighting) {
    }

    virtual void render3D(float scale) = 0;

    void render(size_t width, size_t height) {
        glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(60.0, float(width) / height, 0.2, 1000);
        gluLookAt(0, 0, 200, 0, 0, 0, 0.0, 1.0, 0.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glPushMatrix();
        glTranslatef(m_mouse_vx, -m_mouse_vy, m_mouse_vz);

        glEnable(GL_LINE_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_POINT_SMOOTH);
        glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


        float pi = 3.14159265358;
        svector::float4 qNew(0, 0, 0, 1);
        qNew.euler(m_mouse_dy / MOUSE_FACTOR, m_mouse_dx / MOUSE_FACTOR, 0.0);
        m_camera.quaternion_mult(qNew);
        svector::float4 qAxis = m_camera.axis();
        glRotatef(qAxis.w / pi * 180, qAxis.x, qAxis.y, qAxis.z);
        glRotatef(-90, 1, 0, 0);
        //glRotatef(-90, 0, 0, 1);
        glScalef(1, -1, 1);  // y = -y to have a RHS coordinate system

        glDisable(GL_LIGHTING);
        if (m_drawAxes)
            draw_axes_arrow(200, 200, 200, 0.3);
        if (m_enableLighting)
            enableLight();
        
        render3D(m_scale);
        
        glPopMatrix();
    }

    void mouse_wheel(int button, int dir, int x, int y) {
        // if this glutMouseWheelFunc is supported then we call mouse_button with right values
        // button is always 0 and dir = -1, 1 
        if (dir < 0) mouse_button(GLUT_WHEEL_DOWN, GLUT_DOWN, x, y);
        else mouse_button(GLUT_WHEEL_UP, GLUT_DOWN, x, y);
    }

    void mouse_button(int button, int status, int x, int y) {
        m_left_button_status = GLUT_UP;
        m_right_button_status = GLUT_UP;
        if ((button == GLUT_WHEEL_UP) || (button == GLUT_WHEEL_DOWN)) {
            if (status == GLUT_DOWN) {
                m_mouse_vz += button == GLUT_WHEEL_UP ? MOUSE_ZOOM : -MOUSE_ZOOM;
                if (m_mouse_vz > MOUSE_ZOOM_NEAR) m_mouse_vz = MOUSE_ZOOM_NEAR;     // limit minimum distance
                if (m_mouse_vz < MOUSE_ZOOM_FAR) m_mouse_vz = MOUSE_ZOOM_FAR; // limit maximum distance
            }
        }
        else {
            if (button == GLUT_LEFT_BUTTON) {
                if (status == GLUT_DOWN) {
                    m_left_button_status = GLUT_DOWN;
                    m_left_button_down_x = x;
                    m_left_button_down_y = y;
                }
            }
            if (button == GLUT_RIGHT_BUTTON) {
                if (status == GLUT_DOWN) {
                    m_right_button_status = GLUT_DOWN;
                    m_left_button_down_x = x;
                    m_left_button_down_y = y;
                }
            }
        }
    }
    void mouse_active_motion(int x, int y) {
        if (m_left_button_status == GLUT_DOWN) {
            m_mouse_dx = (x - m_left_button_down_x);
            m_mouse_dy = (y - m_left_button_down_y);
#if !HIGH_DPI_MOUSE
            if (fabs(m_mouse_dx) < 2) m_mouse_dx = 0;
            if (fabs(m_mouse_dy) < 2) m_mouse_dy = 0;
#endif
        }
        if (m_right_button_status == GLUT_DOWN) {
            float dx = (x - m_left_button_down_x) / float(50);
            float dy = (y - m_left_button_down_y) / float(50);
            m_mouse_vy += dy;
            m_mouse_vx += dx;
        }
        m_left_button_down_x = x;
        m_left_button_down_y = y;
    }

    void mouse_passive_motion(int x, int y) {
        m_left_button_down_y = y;
        m_left_button_down_x = x;
    }

    void enableLight() {
        glShadeModel(GL_SMOOTH);
        GLfloat lightingmodelAmbient[] = { 0.8, 0.8, 0.8, 1.0 };
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightingmodelAmbient);

        auto enableLight = [](int light, const svector::float4& pos, bool show) {
            static float position[] = { 0, 0, 0, 1.0f };
            static float colorWhite[] = { 1.0, 1.0, 1.0, 1.0 };
            glPushMatrix();
            glTranslatef(pos.x, pos.y, pos.z);
            if (show) {
                glDisable(GL_LIGHTING);
                glColor3f(1, 1, 0);
                glutSolidSphere(1, 20, 20);
                glEnable(GL_LIGHTING);
            }
            glLightfv(light, GL_POSITION, position);
            glLightfv(light, GL_DIFFUSE, colorWhite);
            glLightfv(light, GL_SPECULAR, colorWhite);
            glEnable(light);
            glPopMatrix();
        };
        
        if (Light0Enabled) enableLight(GL_LIGHT0, Light0, ShowLights);
        if (Light1Enabled) enableLight(GL_LIGHT1, Light1, ShowLights);
        if (Light2Enabled) enableLight(GL_LIGHT2, Light2, ShowLights);

        glColor3f(1, 1, 1);
        
    }
private:
    svector::float4 m_camera;
    int m_left_button_status = 0;
    int m_right_button_status = 0;
    int m_left_button_down_x = 0;
    int m_left_button_down_y = 0;
    float m_scale;
    float m_mouse_dx = 0;
    float m_mouse_dy = 0;
    float m_mouse_vx = 0;
    float m_mouse_vy = 0;
    float m_mouse_vz;
    bool m_drawAxes;
    bool m_enableLighting;
};


