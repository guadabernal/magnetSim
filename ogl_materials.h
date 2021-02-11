#pragma once
#include <GL/glut.h>
#include <GL/glu.h>
#include <vector>

class material {
public:
    material(std::vector<float>&& ambient, std::vector<float>&& diffuse, std::vector<float>&& specular, float shiness)
    : m_ambient(ambient), m_diffuse(diffuse), m_specular(specular), m_shininess(shiness)
    {}
    void active() {
        glMaterialfv(GL_FRONT, GL_AMBIENT, m_ambient.data());
        glMaterialfv(GL_FRONT, GL_DIFFUSE, m_diffuse.data());
        glMaterialfv(GL_FRONT, GL_SPECULAR, m_specular.data());
        glMaterialf(GL_FRONT, GL_SHININESS, m_shininess);
    }
private:
    std::vector<float> m_ambient;
    std::vector<float> m_diffuse;
    std::vector<float> m_specular;
    float m_shininess;
};


static material gold{ {0.24725f, 0.1995f, 0.0745f, 1.0f}, { 0.75164f, 0.60648f, 0.22648f, 1.0f }, { 0.628281f, 0.555802f, 0.366065f, 1.0f }, 51.2f};
static material silver{ {0.23125f, 0.23125f, 0.23125f, 1.0f}, {0.2775f, 0.2775f, 0.2775f, 1.0f}, {0.773911f, 0.773911f, 0.773911f, 1.0f}, 89.6f};
static material chrome{ {0.25f, 0.25f, 0.25f, 1.0f}, {0.4f, 0.4f, 0.4f, 1.0f}, {0.85, 0.85, 0.85, 1.0f}, 76.8f };
static material aluminum{ {0.3f, 0.3f, 0.3f, 1.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, {0.3f, 0.3f, 0.3f, 1.0f}, 10.8f };
