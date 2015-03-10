//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.

struct Sphere{
    string m_name;
    mat4 m_transform;
    mat4 m_invTransform;
    float m_rgb[3];
    float m_Ka; float m_Kd; float m_Ks; float m_Kr;
    float m_n;
    Sphere(string name, mat4 transform, mat4 invTransform, float r, float g, float b,
           float Ka, float Kd, float Ks, float Kr, float n):
           m_name(name), m_transform(transform), m_invTransform(invTransform),
           m_Ka(Ka), m_Kd(Kd), m_Ks(Ks), m_Kr(Kr), m_n(n)
    {
        m_rgb[0] = r;
        m_rgb[1] = g;
        m_rgb[2] = b;
    }
};

struct Light{
    string m_name;
    float m_x; float m_y; float m_z;
    float m_rgb[3];
    Light(string name, float x, float y, float z, float r, float g, float b):
          m_name(name), m_x(x), m_y(y), m_z(z)
    {
        m_rgb[0] = r;
        m_rgb[1] = g;
        m_rgb[2] = b;
    }
};

vector<vec4> g_colors;
vector<Sphere> spheres;
vector<Light> lights;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
    if (vs[0] == "NEAR") {
        g_near = toFloat(vs[1]);
    }
    if (vs[0] == "LEFT") {
        g_left = toFloat(vs[1]);
    }
    if (vs[0] == "RIGHT") {
        g_right = toFloat(vs[1]);
    }
    if (vs[0] == "BOTTOM") {
        g_bottom = toFloat(vs[1]);
    }
    if (vs[0] == "TOP") {
        g_top = toFloat(vs[1]);
    }
    if (vs[0] == "SPHERE") {
        string name = vs[1];
        mat4 transform = Scale(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7])) * Translate(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]));
        mat4 invTransform;
        InvertMatrix(transform, invTransform);
        Sphere sphere(name, transform, invTransform, toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10]), toFloat(vs[11]), toFloat(vs[12]), toFloat(vs[13]), toFloat(vs[14]), toFloat(vs[15]));
        spheres.push_back(sphere);
    }
    if (vs[0] == "LIGHT") {
        Light light(vs[1], toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
        lights.push_back(light);
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    return vec4(0.0f, 0.0f, 0.0f, 1.0f);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    dir = vec4(0.0f, 0.0f, -1.0f, 0.0f);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, "output.ppm", buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

