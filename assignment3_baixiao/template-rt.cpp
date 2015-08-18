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
#include <ctime>
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
    float m_t;
    vec4 intersection;
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

float g_ambiant[3];
float g_backgrond[3];

string g_outputFile;

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
        InvertMatrix( Translate(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4])) * Scale(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7])), invTransform);
        
        Sphere sphere(name, transform, invTransform, toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10]), toFloat(vs[11]), toFloat(vs[12]), toFloat(vs[13]), toFloat(vs[14]), toFloat(vs[15]));
        spheres.push_back(sphere);
    }
    if (vs[0] == "LIGHT") {
        Light light(vs[1], toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]), toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
        lights.push_back(light);
    }
    if (vs[0] == "BACK") {
        g_backgrond[0] = toFloat(vs[1]);
        g_backgrond[1] = toFloat(vs[2]);
        g_backgrond[2] = toFloat(vs[3]);
    }
    if (vs[0] == "AMBIENT") {
        g_ambiant[0] = toFloat(vs[1]);
        g_ambiant[1] = toFloat(vs[2]);
        g_ambiant[2] = toFloat(vs[3]);
    }
    if (vs[0] == "OUTPUT") {
        g_outputFile = vs[1];
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

vec4 trace(const Ray& ray, int depth=-1, int mySphere=-1)
{
    // TODO: implement your ray tracing routine here.
    if (depth == 0) {
        return vec4(0.0f, 0.0f, 0.0f, 0.0f);
    }
    
    vector<vec4> t_origin;
    vector<vec4> t_dir;
    vector<float> intersectionZ;
    vec4 rgb;
    if (depth == -1){
        rgb.x = g_backgrond[0]; rgb.y = g_backgrond[1]; rgb.z = g_backgrond[2]; rgb.w = 1.0;
    }
    else{
        rgb.x = 0; rgb.y = 0; rgb.z = 0; rgb.w = 1.0;
    }
    //inverse transform the ray
    for (int i=0; i<spheres.size(); i++) {
        t_origin.push_back(spheres[i].m_invTransform * ray.origin);
        t_dir.push_back(spheres[i].m_invTransform * ray.dir);
    }
    //find intersection of ray for each sphere
    for (int i=0; i<spheres.size(); i++) {
        float A = length(t_dir[i]) * length(t_dir[i]);
        float B = (t_origin[i][0]*t_dir[i][0] + t_origin[i][1]*t_dir[i][1] + t_origin[i][2]*t_dir[i][2]) * 2 / t_origin[i][3];
        float C = t_origin[i][0]*t_origin[i][0] + t_origin[i][1]*t_origin[i][1] + t_origin[i][2]*t_origin[i][2] - 1;
        float solnFactor = B*B - 4*A*C;
        if (solnFactor >= 0) {
            float t1 = (-B - sqrt(solnFactor))/(2*A);
            float t2 = (-B + sqrt(solnFactor))/(2*A);
            vec4 intersectionT1 = ray.origin + ray.dir * t1;
            vec4 intersectionT2 = ray.origin + ray.dir * t2;
            
            if (intersectionT1.z < -g_near && t1 > 0) {
                intersectionZ.push_back(intersectionT1.z);
                spheres[i].m_t = t1;
                spheres[i].intersection = intersectionT1;
            } else if (intersectionT2.z < -g_near && t2 > 0){
                intersectionZ.push_back(intersectionT2.z);
                spheres[i].m_t = t2;
                spheres[i].intersection = intersectionT2;
            } else{
                // no valid intersection
                intersectionZ.push_back(0);
                spheres[i].m_t = 0;
            }
        }else if (solnFactor < 0){
            //no intersection for this sphere
            intersectionZ.push_back(0);
            spheres[i].m_t = 0;
        }
    }
    
    //get the front most sphere
    float sphereNum = -1;
    float minDistance = 1000;
    for (int i=0; i<spheres.size(); i++) {
        if ((0 - intersectionZ[i]) < minDistance && intersectionZ[i] != 0) {
            minDistance = 0 - intersectionZ[i];
            sphereNum = i;
        }
    }
    //we can find an object
    if (sphereNum != mySphere && sphereNum != -1){
        //get surface normal vector
        vec4 normal = t_origin[sphereNum] + t_dir[sphereNum] * spheres[sphereNum].m_t;
        normal.w = 0;
        normal = spheres[sphereNum].m_invTransform * normal;
        normal = normal / length(normal);
        vec4 eye = -ray.dir;
        eye = eye / length(eye);
        
        rgb.x = 0; rgb.y = 0; rgb.z = 0;
        for (int l=0; l<lights.size(); l++){
            vec4 lightVec;
            vec4 refVec;
            lightVec.x = lights[l].m_x-spheres[sphereNum].intersection.x;
            lightVec.y = lights[l].m_y-spheres[sphereNum].intersection.y;
            lightVec.z = lights[l].m_z-spheres[sphereNum].intersection.z;
            lightVec.w = 0;
            vec4 u_lightVec = lightVec;
            lightVec = lightVec / length(lightVec);
            
            refVec = 2*normal*dot(normal, lightVec) - lightVec;
            
            //check for shadow
            bool obstructed = false;
            for (int i=0; i<spheres.size(); i++) { // loop through spheres
                //if (i == sphereNum) continue; //checking against itself
                vec4 trans_origin = spheres[i].m_invTransform * spheres[sphereNum].intersection;
                vec4 trans_dir = spheres[i].m_invTransform * lightVec;
                float A = length(trans_dir) * length(trans_dir);
                float B = (trans_origin[0]*trans_dir[0] + trans_origin[1]*trans_dir[1] + trans_origin[2]*trans_dir[2]) * 2 / trans_origin[3];
                float C = trans_origin[0]*trans_origin[0] + trans_origin[1]*trans_origin[1] + trans_origin[2]*trans_origin[2] - 1;
                float solnFactor = B*B - 4*A*C;
                if (solnFactor >= 0) {
                    float t1 = (-B - sqrt(solnFactor))/(2*A);
                    float t2 = (-B + sqrt(solnFactor))/(2*A);
                    if (t1 > 0.01 || t2 > 0.01){ // make sure t is not close to zero
                        obstructed = true;
                        break;
                    }
                }
            }
            if (obstructed) {
                continue;
            }
        
            
            float dotProduct = dot(normal, lightVec);
            if (dotProduct > 1 || dotProduct < 0) {
                dotProduct = 0;
            }
            //add diffusion
            rgb.x += spheres[sphereNum].m_rgb[0]*lights[l].m_rgb[0]*spheres[sphereNum].m_Kd*dotProduct;
            rgb.y += spheres[sphereNum].m_rgb[1]*lights[l].m_rgb[1]*spheres[sphereNum].m_Kd*dotProduct;
            rgb.z += spheres[sphereNum].m_rgb[2]*lights[l].m_rgb[2]*spheres[sphereNum].m_Kd*dotProduct;
            
            //add specular component
            float specularDotproduct = dot(refVec, eye);
            float result = 1;
            if (specularDotproduct > 1 || specularDotproduct < 0) {
                specularDotproduct = 0;
            }
            for (int i=0; i<spheres[sphereNum].m_n; i++) {
                result *= specularDotproduct;
            }
            rgb.x += lights[l].m_rgb[0]*spheres[sphereNum].m_Ks*result;
            rgb.y += lights[l].m_rgb[1]*spheres[sphereNum].m_Ks*result;
            rgb.z += lights[l].m_rgb[2]*spheres[sphereNum].m_Ks*result;
            
        }
        //ambiant lights
        if (spheres[sphereNum].m_Ka != 0) {
            rgb.x += spheres[sphereNum].m_rgb[0]*g_ambiant[0] * spheres[sphereNum].m_Ka;
            rgb.y += spheres[sphereNum].m_rgb[1]*g_ambiant[1] * spheres[sphereNum].m_Ka;
            rgb.z += spheres[sphereNum].m_rgb[2]*g_ambiant[2] * spheres[sphereNum].m_Ka;
        }
        
        //reflection ray
        vec4 nextRay_origin = spheres[sphereNum].intersection;
        vec4 nextRay_dir = 2 * normal * dot(normal, eye) - eye;
        nextRay_dir = nextRay_dir / length(nextRay_dir);
        Ray nextRay; nextRay.origin = nextRay_origin; nextRay.dir = nextRay_dir;
        if (depth == -1) {
            vec4 reflection = trace(nextRay, 3, sphereNum);
            rgb += spheres[sphereNum].m_Kr * reflection;
        }else{
            vec4 reflection = trace(nextRay, depth-1, sphereNum);
            rgb += spheres[sphereNum].m_Kr * reflection;
        }
    }
    //clamp the values
    rgb[0] = rgb[0] <= 1.0f ? rgb[0] : 1.0;
    rgb[1] = rgb[1] <= 1.0f ? rgb[1] : 1.0;
    rgb[2] = rgb[2] <= 1.0f ? rgb[2] : 1.0;
    return rgb;
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    float x = g_left + (ix/(float)g_width)*(g_right-g_left);
    float y = g_bottom + (iy/(float)g_height)*(g_top-g_bottom);
    dir = vec4(x, y, -1.0f, 0.0f);
    //dir = dir / length(dir);
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
    char * filename = new char[g_outputFile.size() + 1];
    copy(g_outputFile.begin(), g_outputFile.end(), filename);
    filename[g_outputFile.size()] = '\0';
    savePPM(g_width, g_height, filename, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    //float start = clock();
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
    //printf("program lasts %f\n", (clock()-start)/(float)CLOCKS_PER_SEC);
	return 0;
}

