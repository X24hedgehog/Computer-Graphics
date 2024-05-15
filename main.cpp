#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <math.h>
#include <ctime>
#include <iostream>
#define M_PI 3.14159265358979323846

#include <string>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <omp.h>

double epsilon = 1e-9;
int num_ray_max = 2;
int num_reflect_max = 20; // directlighting: 1

double double_rand()
{
    static thread_local std::mt19937 *generator = nullptr;
    if (!generator)
        generator = new std::mt19937(clock() + omp_get_thread_num());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(*generator);
}

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const
    {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const
    {
        return sqrt(norm2());
    }
    void normalize()
    {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double &operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector gamma_correction(const Vector &a)
{
    return Vector(std::min<double>(255, pow(a[0], 1 / 2.2)), std::min<double>(255, pow(a[1], 1 / 2.2)), std::min<double>(255, pow(a[2], 1 / 2.2)));
}
Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

class Ray
{
public:
    Vector O;
    Vector u;
    Ray(const Vector &Origin, const Vector &direction)
    {
        O = Origin;
        u = direction;
    }
};

class Geometry
{
public:
    Vector Alb;
    bool m;
    Geometry(Vector Albedo, bool IsMirror)
    {
        Alb = Albedo;
        m = IsMirror;
    }
    virtual double check_intersect(const Ray &ray, Vector &N) = 0;
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class BoundingBox
{
public:
    BoundingBox()
    {
    }

    Vector verticemin = Vector(INT_MAX, INT_MAX, INT_MAX);
    Vector verticemax = -1 * Vector(INT_MAX, INT_MAX, INT_MAX);

    void change_bound(std::vector<Vector> &indice_list)
    {
        for (int i = 0; i < indice_list.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (indice_list[i][j] > verticemax[j])
                {
                    verticemax[j] = indice_list[i][j];
                }
                else if (indice_list[i][j] < verticemin[j])
                {
                    verticemin[j] = indice_list[i][j];
                }
            }
        }
    }

    bool check_intersect(const Ray &ray)
    {
        Vector left_bound, right_bound;
        for (int i = 0; i < 3; i++)
        {
            left_bound[i] = (verticemin[i] - ray.O[i]) / ray.u[i];
            right_bound[i] = (verticemax[i] - ray.O[i]) / ray.u[i];
            if (left_bound[i] > right_bound[i])
            {
                double temp = right_bound[i];
                right_bound[i] = left_bound[i];
                left_bound[i] = temp;
            }
        }
        double min_bound = std::max<double>(std::max<double>(left_bound[0], left_bound[1]), left_bound[2]);
        double max_bound = std::min<double>(std::min<double>(right_bound[0], right_bound[1]), right_bound[2]);
        return (min_bound < max_bound && max_bound > 0);
    }
};

class TriangleMesh : public Geometry
{
public:
    TriangleMesh(Vector Albedo, bool IsMirror = false) : Geometry(Albedo, IsMirror)
    {
    }

    double check_intersect_triangle(const Ray &ray, TriangleIndices &indice, Vector &final_N)
    {
        Vector origin = ray.O;
        Vector direction = ray.u;
        Vector A = vertices[indice.vtxi];
        Vector B = vertices[indice.vtxj];
        Vector C = vertices[indice.vtxk];
        Vector e1 = B - A;
        Vector e2 = C - A;
        Vector N = cross(e1, e2);
        double beta = dot(e2, cross(A - origin, direction)) / dot(direction, N);
        double gamma = -1 * dot(e1, cross(A - origin, direction)) / dot(direction, N);
        double alpha = 1 - beta - gamma;

        if (beta < 0 || gamma < 0 || alpha < 0)
        {
            return -1;
        }

        double t = dot(A - origin, N) / dot(direction, N);
        if (t < 0)
        {
            return -1;
        }
        N.normalize();
        final_N = N;
        return t;
    }

    double check_intersect(const Ray &ray, Vector &N) override
    {
        // if (box.check_intersect(ray))
        // {
        double t = -1;
        double temp;
        TriangleIndices object;
        Vector N_temp;
        for (int i = 0; i < indices.size(); i++)
        {
            temp = check_intersect_triangle(ray, indices[i], N_temp);
            if ((t < 0) || (temp > 0 && temp < t))
            {
                t = temp;
                N = N_temp;
            }
        }
        return t;
        // }
        // return -1;
    }

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vec = vec * 0.8 + Vector(0, -10, 0);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
        box.change_bound(vertices);
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BoundingBox box;
};

class Sphere : public Geometry
{
public:
    double R;
    Vector C;
    Sphere(double Radius, Vector Center, Vector Albedo, bool IsMirror = false) : Geometry(Albedo, IsMirror)
    {
        R = Radius;
        C = Center;
    }

    double check_intersect(const Ray &ray, Vector &N) override
    {
        double delta = pow(dot(ray.u, ray.O - C), 2) - (ray.O - C).norm2() + pow(R, 2);
        if (delta >= 0)
        {
            double t1 = dot(ray.u, C - ray.O) + sqrt(delta);
            double t2 = dot(ray.u, C - ray.O) - sqrt(delta);
            if (t2 >= 0)
            {
                Vector P = ray.O + t2 * ray.u;
                N = P - C;
                N.normalize();
                return t2;
            }
            else
            {
                if (t1 >= 0)
                {
                    Vector P = ray.O + t1 * ray.u;
                    N = P - C;
                    N.normalize();
                    return t1;
                }
            }
        }
        return -1.;
    }
};

Vector randomCos(const Vector &N)
{
    double r1 = double_rand();
    double r2 = double_rand();
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    int iMin = 0;
    double min = abs(N[0]);
    for (int i = 1; i < 3; i++)
    {
        if (abs(N[i]) < min)
        {
            min = abs(N[i]);
            iMin = i;
        }
    }

    Vector T1;
    if (iMin == 0)
        T1 = Vector(0, N[2], -N[1]);
    else if (iMin == 1)
        T1 = Vector(N[2], 0, -N[0]);
    else if (iMin == 2)
        T1 = Vector(N[1], -N[0], 0);
    T1.normalize();

    Vector T2 = cross(T1, N);

    return T1 * x + T2 * y + N * z;
}

class Scene
{
public:
    std::vector<Geometry *> G_list;
    Vector cam_pos;
    int H;
    int W;
    double alpha;
    Vector l_pos;
    double l_int;
    double l_rad;
    Scene(const Vector &camera_position, int height, int width, double angle, Vector light_position, double light_intensity, double light_radius)
    {
        cam_pos = camera_position;
        H = height;
        W = width;
        alpha = angle;
        l_pos = light_position;
        l_int = light_intensity;
        l_rad = light_radius;
    }
    void add_geometry(Geometry *object_pointer)
    {
        G_list.push_back(object_pointer);
    }

    Vector get_pixel_color(const Ray &ray, int num_reflex = 0)
    {
        if (num_reflex >= num_reflect_max)
        {
            return Vector(0., 0., 0.);
        }
        Vector Lo(0., 0., 0.);
        double t = -1.;
        Vector N_temp, N;
        int object_index;
        for (int i = 0; i < G_list.size(); i++)
        {
            double temp = G_list[i]->check_intersect(ray, N_temp);
            if ((t < 0) || (temp > 0 && temp < t))
            {
                t = temp;
                object_index = i;
                N = N_temp;
            }
        }
        if (t < 0)
        {
            return Vector(4000, 4000, 4000);
        }
        else
        {
            Vector P = ray.O + t * ray.u;
            // Vector N = P - S_list[sphere_index].C;
            N.normalize();
            P = P + epsilon * N;
            if (G_list[object_index]->m == true)
            {
                Vector reflection_direction = ray.u - 2 * (dot(ray.u, N)) * N;
                return get_pixel_color(Ray(P, reflection_direction), num_reflex + 1);
            }
            else
            {
                double visibility = 1.;
                for (int i = 0; i < G_list.size(); i++)
                {
                    Vector point_to_light = l_pos - P;
                    point_to_light.normalize();
                    double point_to_light_intersect = G_list[i]->check_intersect(Ray(P, point_to_light), N_temp);
                    if ((point_to_light_intersect >= 0) && (point_to_light_intersect < (l_pos - P).norm()))
                    {
                        // std::cout << i;
                        visibility = 0.;
                        break;
                    }
                }
                // std::cout << '1';
                // Lo = (visibility) ? Vector(255., 0., 0.) : Vector(0., 255., 0.);
                double I = l_int;
                Vector l_dir = l_pos - P;
                l_dir.normalize();
                double component1 = I / (4 * M_PI * (l_pos - P).norm2());
                Vector component2 = G_list[object_index]->Alb / M_PI;
                double component3 = dot(N, l_dir);
                Lo = gamma_correction(component1 * component2 * component3) * visibility;

                // double I = l_int;
                // double R = l_rad;
                // Vector x = G_list[object_index]->C;
                // Vector x_to_light = x - l_pos;
                // x_to_light.normalize();
                // Vector xprime = R * randomCos(x_to_light) + l_pos;
                // Vector Nprime = xprime - l_pos;
                // Nprime.normalize();
                // double d = (xprime - P).norm();
                // Vector omega = xprime - P;
                // omega.normalize();
                // double pdf = dot(Nprime, x_to_light) / (M_PI * R * R);
                // Vector rho = S_list[sphere_index].Alb;
                // Lo = I / pow(2 * M_PI * R, 2) * rho / M_PI * visibility * std::max(dot(N, omega), 0.) * std::max(dot(Nprime, (-1.) * omega), 0.) / (pow((xprime - P).norm(), 2) * pdf);

                Ray myRandomRay = Ray(P, randomCos(N));
                return Lo + G_list[object_index]->Alb * get_pixel_color(myRandomRay, num_reflex + 1);
            }
        }
    }

    void render(std::vector<unsigned char> &image)
    {
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < H; i++)
        {
            for (int j = 0; j < W; j++)
            {
                Vector total_color(0., 0., 0.);
                for (int num_ray = 0; num_ray < num_ray_max; num_ray++)
                {
                    double x_ray = j - double(W) / 2. + 0.5;
                    double y_ray = double(H) / 2. - i - 0.5;
                    double z_ray = -W / (2. * tan(alpha / 2.));

                    Vector direction_to_pixel = Vector(x_ray, y_ray, z_ray);
                    direction_to_pixel.normalize();
                    Vector output_color = get_pixel_color(Ray(cam_pos, direction_to_pixel));
                    total_color = total_color + output_color;
                }
                // std::cout << total_color[0] << ' ' << total_color[1] << ' ' << total_color[2] << std::endl;
                Vector average_color = total_color / num_ray_max;
                // Vector gamma_color = gamma_correction(average_color);

                image[(i * W + j) * 3 + 0] = std::min(255., average_color[0]);
                image[(i * W + j) * 3 + 1] = std::min(255., average_color[1]);
                image[(i * W + j) * 3 + 2] = std::min(255., average_color[2]);
            }
        }
    }
};

int main()
{
    int W = 512;
    int H = 512;
    time_t start, end;
    std::time(&start);

    // std::vector<unsigned char> image(W * H * 3, 0);
    // for (int i = 0; i < H; i++)
    // {
    //     for (int j = 0; j < W; j++)
    //     {

    //         image[(i * W + j) * 3 + 0] = 255;
    //         image[(i * W + j) * 3 + 1] = 0;
    //         image[(i * W + j) * 3 + 2] = 0;
    //     }
    // }

    // return 0;
    Scene scene = Scene(Vector(0, 0, 55), W, H, M_PI / 3., Vector(-10, 20, 40), 1e10, 3.);
    std::vector<Sphere> object_list;
    TriangleMesh Cat = TriangleMesh(Vector(0.3, 0.3, 0.3), false);
    Cat.readOBJ("cat.obj");
    object_list.push_back(Sphere(10, Vector(-20, 0, 0), Vector(1., 1., 1.)));
    object_list.push_back(Sphere(10, Vector(0, 0, 0), Vector(1., 1., 1.), true));
    object_list.push_back(Sphere(10, Vector(20, 0, 0), Vector(1., 1., 1.)));
    object_list.push_back(Sphere(940, Vector(0, 0, 1000), Vector(1., 0., 1.)));
    object_list.push_back(Sphere(940, Vector(0, 1000, 0), Vector(1., 0., 0.)));
    object_list.push_back(Sphere(940., Vector(0, 0, -1000), Vector(0., 1., 0.)));
    object_list.push_back(Sphere(990., Vector(0, -1000, 0), Vector(0., 0., 1.)));
    object_list.push_back(Sphere(940., Vector(-1000, 0, 0), Vector(0., 1., 1.)));
    object_list.push_back(Sphere(940., Vector(1000, 0, 0), Vector(1., 1., 0.)));

    for (int i = 0; i < object_list.size(); i++)
    {
        scene.add_geometry(&object_list[i]);
    }
    // scene.add_geometry(&Cat);

    std::vector<unsigned char> image(W * H * 3, 0);

    scene.render(image);
    stbi_write_png("image2_1.png", W, H, 3, &image[0], 0);
    std::time(&end);
    std::cout << end - start << std::endl;
}