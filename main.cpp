#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <math.h>
#define M_PI 3.14159265358979323846

double epsilon = 1e-9;

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

class Sphere
{
public:
    double R;
    Vector C;
    Vector Alb;
    bool m;
    Sphere(double Radius, Vector Center, Vector Albedo, bool IsMirror = false)
    {
        R = Radius;
        C = Center;
        Alb = Albedo;
        m = IsMirror;
    }

    double check_intersect(const Ray &ray)
    {
        double t = -1.;
        double delta = pow(dot(ray.u, ray.O - C), 2) - (ray.O - C).norm2() + pow(R, 2);
        if (delta >= 0)
        {
            double t1 = dot(ray.u, C - ray.O) + sqrt(delta);
            double t2 = dot(ray.u, C - ray.O) - sqrt(delta);
            if (t2 >= 0)
            {
                return t2;
            }
            else
            {
                if (t1 >= 0)
                {
                    return t1;
                }
            }
        }
        return t;
    }
};

class Scene
{
public:
    std::vector<Sphere> S_list;
    Vector cam_pos;
    int H;
    int W;
    double alpha;
    Vector l_pos;
    double l_int;
    Scene(const Vector &camera_position, int height, int width, double angle, Vector light_position, double light_intensity)
    {
        cam_pos = camera_position;
        H = height;
        W = width;
        alpha = angle;
        l_pos = light_position;
        l_int = light_intensity;
    }
    void add_sphere(const Sphere &sphere)
    {
        S_list.push_back(sphere);
    }

    Vector get_pixel_color(const Ray &ray, int num_reflex = 0)
    {
        double t = -1.;
        int sphere_index;
        for (int i = 0; i < S_list.size(); i++)
        {
            double temp = S_list[i].check_intersect(ray);
            if ((t < 0) || (temp > 0 && temp < t))
            {
                t = temp;
                sphere_index = i;
            }
        }
        if (t < 0)
        {
            return Vector(4000, 4000, 4000);
        }
        else
        {
            Vector P = ray.O + t * ray.u;
            Vector N = P - S_list[sphere_index].C;
            N.normalize();
            P = P + epsilon * N;
            if (S_list[sphere_index].m == true)
            {
                if (num_reflex < 2)
                {
                    Vector reflection_direction = ray.u - 2 * (dot(ray.u, N)) * N;
                    return get_pixel_color(Ray(P, reflection_direction), num_reflex + 1);
                }
                else
                {
                    return (Vector(0., 0., 0.));
                }
            }
            bool intersect = false;
            for (int i = 0; i < S_list.size(); i++)
            {
                Vector point_to_light = l_pos - P;
                point_to_light.normalize();
                double point_to_light_intersect = S_list[i].check_intersect(Ray(P, point_to_light));
                if ((point_to_light_intersect >= 0) && (point_to_light_intersect < (l_pos - P).norm()))
                {
                    intersect = true;
                    break;
                }
            }
            if (intersect == false)
            {

                double color1 = l_int / (4 * M_PI * (l_pos - P).norm2());
                Vector color2 = S_list[sphere_index].Alb / M_PI;
                double color3 = dot(N, (l_pos - P) / (l_pos - P).norm());
                Vector final_color = color1 * color2 * color3;
                return final_color;
            }
            else
            {
                return Vector(0, 0, 0);
            }
        }
    }
    void render(std::vector<unsigned char> &image)
    {
        for (int i = 0; i < H; i++)
        {
            for (int j = 0; j < W; j++)
            {
                double x_ray = j - double(W) / 2. + 0.5;
                double y_ray = double(H) / 2. - i - 0.5;
                double z_ray = -W / (2. * tan(alpha / 2.));

                Vector direction_to_pixel = Vector(x_ray, y_ray, z_ray);
                direction_to_pixel.normalize();
                Vector output_color = get_pixel_color(Ray(cam_pos, direction_to_pixel));
                Vector gamma_color = gamma_correction(output_color);

                image[(i * W + j) * 3 + 0] = gamma_color[0];
                image[(i * W + j) * 3 + 1] = gamma_color[1];
                image[(i * W + j) * 3 + 2] = gamma_color[2];
            }
        }
    }
};

int main()
{
    int W = 512;
    int H = 512;

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
    Scene scene = Scene(Vector(0, 0, 55), 512, 512, M_PI / 3., Vector(-10, 20, 40), 1e10);
    scene.add_sphere(Sphere(10, Vector(-20, 0, 0), Vector(1., 1., 1.)));
    scene.add_sphere(Sphere(10, Vector(0, 0, 0), Vector(1., 1., 1.), true));
    scene.add_sphere(Sphere(10, Vector(20, 0, 0), Vector(1., 1., 1.)));
    scene.add_sphere(Sphere(940, Vector(0, 0, 1000), Vector(1., 0., 1.)));
    scene.add_sphere(Sphere(940, Vector(0, 1000, 0), Vector(1., 0., 0.)));
    scene.add_sphere(Sphere(940., Vector(0, 0, -1000), Vector(0., 1., 0.)));
    scene.add_sphere(Sphere(990., Vector(0, -1000, 0), Vector(0., 0., 1.)));
    scene.add_sphere(Sphere(940., Vector(-1000, 0, 0), Vector(0., 1., 1.)));
    scene.add_sphere(Sphere(940., Vector(1000, 0, 0), Vector(1., 1., 0.)));

    std::vector<unsigned char> image(W * H * 3, 0);

    scene.render(image);
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
}