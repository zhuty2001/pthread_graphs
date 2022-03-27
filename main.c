#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "image.h"
#include "scene.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

/**
 * @brief You never need to read this function.
 * 
 * @details Randomly generate a ray from camera.
 */
ray generate_ray(const ray *camera, unsigned short xi[3],
                 const int image_width, const int image_height,
                 const int x, const int y,
                 const vec3 *cx, const vec3 *cy,
                 const int sx, const int sy) {
    /* For the current subpixel, randomly generate direction of the ray */
    double r1 = 2 * erand48(xi);
    double r2 = 2 * erand48(xi);
    double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
    double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
    vec3 camera_ray_dir = vec3_unit(vec3_add(vec3_mulf(*cx, (((sx + .5 + dx) / 2 + x) / image_width - .5)),
                                             vec3_add(vec3_mulf(*cy, (((sy + .5 + dy) / 2 + y) / image_height - .5)), camera->d)));
    /* Shoot the ray from camera */
    return ray_init(vec3_add(camera->o, vec3_mulf(camera_ray_dir, camera_offset)), camera_ray_dir); /* Offset = 140 is added to place the ray deep into the scene because the box-like walls are actually made of large spheres */
}

struct Ppthread
{
    int image_height,image_width,samples_per_pixel;
    vec3 *image;
    ray* camera;
    vec3* cx;
    vec3* cy;
    int bef;
};



/**
 * @brief You never need to read this function.
 * 
 * @details Ray-scene intersection test, return false if no intersection, 
 *          otherwise record the nearest hit time into `t` and index of the corresponding sphere into `id`.
 */
bool intersect(const ray *r, double *t, int *id) {
    const int n = sizeof(spheres) / sizeof(sphere);
    const double inf = 1e20;
    double d;
    int i;
    *t = inf;
    for (i = n; i--; )
        if ((d = sphere_intersect(r, &spheres[i])) && d < *t) {
            *t = d;
            *id = i;
        }
    return *t < inf;
}



/**
 * @brief You never need to read this function.
 * 
 * @details Called recursively, trace the ray to spheres 
 *          and return the color got by the ray.
 */
vec3 radiance(const ray r, int depth, unsigned short *xi) {
    double t; /* Hit time */
    int id = 0; /* ID of intersected object */
    const sphere *obj;
    vec3 x, n, nl, f;
    double p;
    if (depth > 10 || !intersect(&r, &t, &id)) /* If miss, return black */
        return vec3_init(0, 0, 0);
    obj = &spheres[id]; /* The hit object */
    x = vec3_add(r.o, vec3_mulf(r.d, t));
    n = vec3_unit(vec3_sub(x, obj->position));
    nl = vec3_dot(n, r.d) < 0 ? n : vec3_mulf(n, -1);
    f = obj->color;
    p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; /* Max extent as the absorption rate */
    if (++depth > 5) { /* Russian roulette */
        if (erand48(xi) < p)
            f = vec3_mulf(f, 1 / p);
        else
            return obj->emission;
    }
    if (obj->material == DIFFUSE) { /* Ideal diffuse reflection */
        double r1 = 2 * PI * erand48(xi), r2 = erand48(xi), r2s = sqrt(r2);
        vec3 w = nl;
        vec3 u = vec3_unit(vec3_cross(fabs(w.x) > .1 ? vec3_init(0, 1, 0) : vec3_init(1, 0, 0), w));
        vec3 v = vec3_cross(w, u);
        vec3 d = vec3_unit(vec3_add(vec3_mulf(u, cos(r1) * r2s), vec3_add(vec3_mulf(v, sin(r1) * r2s), vec3_mulf(w, sqrt(1 - r2)))));
        return vec3_add(obj->emission, vec3_mul(f, radiance(ray_init(x, d), depth, xi)));
    }
    else if (obj->material == SPECULAR) { /* Ideal specular reflection */
        return vec3_add(obj->emission, vec3_mul(f, radiance(ray_init(x, vec3_sub(r.d, vec3_mulf(n, (2 * vec3_dot(n, r.d))))), depth, xi)));
    }
    else if (obj->material == REFRACTIVE) { /* Ideal dielectric refraction */
        ray r_reflect = ray_init(x, vec3_sub(r.d, vec3_mulf(n, 2 * vec3_dot(n, r.d))));
        vec3 tdir;
        ray r_transmit;
        bool into = vec3_dot(n, nl) > 0; /* Whether ray from outside going in */
        double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = vec3_dot(r.d, nl), cos2t;
        double a, b, R0, c;
        double Re, Tr, P, RP, TP;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) /* Total internal reflection */
            return vec3_add(obj->emission, vec3_mul(f, radiance(r_reflect, depth, xi)));
        tdir = vec3_unit(vec3_sub(vec3_mulf(r.d, nnt), vec3_mulf(n, ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))));
        r_transmit = ray_init(x, tdir);
        a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : vec3_dot(tdir, n));
        Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
        return vec3_add(obj->emission, vec3_mul(f, (depth > 2 ? (erand48(xi) < P ? /* Russian roulette */
                                                                 vec3_mulf(radiance(r_reflect, depth, xi), RP) : vec3_mulf(radiance(r_transmit, depth, xi), TP)) :
                                                    vec3_add(vec3_mulf(radiance(r_reflect, depth, xi), Re), vec3_mulf(radiance(r_transmit, depth, xi), Tr)))));
    }
    else
        return vec3_init(0, 0, 0);
}



/**
 * @brief You never need to understand the ray tracing algorithm, 
 *        but you can still make it parallel by using pthread.
 * 
 * @details Given the image settings, camera settings and scene settings, 
 *          trace rays from camera into the scene and output an amazing image of the scene.
 */

void* the_pthread(void* p)
{
    int i;
    int y, x, sy, sx, s;
    int image_width = 1024;
    int image_height = 768;
    struct Ppthread *temp_p=(struct Ppthread*)p;
    for (y = temp_p->bef; y < temp_p->bef + image_height / 4; ++y) {
        unsigned short xi[3] = { 0, 0 };
        xi[2] = y * y * y;
        for (x = 0; x < image_width; ++x) {
            i = (image_height - y - 1) * image_width + x;
            for (sy = 0; sy < 2; ++sy) {
                for (sx = 0; sx < 2; ++sx) {
                    vec3 accumulated_color = vec3_init(0, 0, 0);
                    for (s = 0; s < temp_p->samples_per_pixel; s++) {
                        /* Randomly generate a ray from camera */
                        ray camera_ray = generate_ray(temp_p->camera, xi, temp_p->image_width, temp_p->image_height,
                                                      x, y, temp_p->cx, temp_p->cy, sx, sy);
                        /* Trace the ray and compute the color */
                        vec3 camera_ray_color = radiance(camera_ray, 0, xi);
                        /* Add the color to `accumulated_color` for this subpixel */
                        accumulated_color = vec3_add(accumulated_color, vec3_mulf(camera_ray_color, (1. / temp_p->samples_per_pixel)));
                    }
                    /* Add `accumulated_color` for this subpixel to the color for the image pixel, which can be imagined as 1 = 0.25 * 4 */
                    temp_p->image[i] = vec3_add(temp_p->image[i], vec3_mulf(vec3_init(clamp(accumulated_color.x), clamp(accumulated_color.y), clamp(accumulated_color.z)), .25));
                }
            }
        }
    }
    return NULL;
}
int main(int argc, char *argv[]) {
    /* Image settings */
    int image_width = 1024;
    int image_height = 768;
    int samples_per_pixel = argc == 2 ? atoi(argv[1]) / 4 : 1;
    vec3 *image = malloc(sizeof(vec3) * image_width * image_height);

    /* Camera settings */
    ray camera = ray_init(vec3_init(50, 52, 295.6), vec3_unit(vec3_init(0, -0.042612, -1)));
    vec3 cx = vec3_init(image_width * .5135 / image_height, 0, 0);
    vec3 cy = vec3_mulf(vec3_unit(vec3_cross(cx, camera.d)), .5135);

    /* Ray tracing */
    pthread_t tid1;
    pthread_t tid2;
    pthread_t tid3;
    pthread_t tid4;
    struct Ppthread *ptr1 = (struct Ppthread*)malloc(sizeof(struct Ppthread));
    struct Ppthread *ptr2 = (struct Ppthread*)malloc(sizeof(struct Ppthread));
    struct Ppthread *ptr3 = (struct Ppthread*)malloc(sizeof(struct Ppthread));
    struct Ppthread *ptr4 = (struct Ppthread*)malloc(sizeof(struct Ppthread));
    ptr1->samples_per_pixel = samples_per_pixel;
    ptr2->samples_per_pixel = samples_per_pixel;
    ptr3->samples_per_pixel = samples_per_pixel;
    ptr4->samples_per_pixel = samples_per_pixel;
    ptr1->image=image;
    ptr2->image=image;
    ptr3->image=image;
    ptr4->image=image;
    ptr1->camera=&camera;
    ptr2->camera=&camera;
    ptr3->camera=&camera;
    ptr4->camera=&camera;
    ptr1->image_height=image_height;
    ptr2->image_height=image_height;
    ptr3->image_height=image_height;
    ptr4->image_height=image_height;
    ptr1->image_width=image_width;
    ptr2->image_width=image_width;
    ptr3->image_width=image_width;
    ptr4->image_width=image_width;
    ptr1->cx = &cx;
    ptr1->cy = &cy;
    ptr2->cx = &cx;
    ptr2->cy = &cy;
    ptr3->cx = &cx;
    ptr3->cy = &cy;
    ptr4->cx = &cx;
    ptr4->cy = &cy;
    ptr1->bef=(0);
    ptr2->bef=(image_height*1/4);
    ptr3->bef=(image_height*2/4);
    ptr4->bef=(image_height*3/4);
    pthread_create(&tid1,NULL,the_pthread,(void*)ptr1);
    pthread_create(&tid2,NULL,the_pthread,(void*)ptr2);
    pthread_create(&tid3,NULL,the_pthread,(void*)ptr3);
    pthread_create(&tid4,NULL,the_pthread,(void*)ptr4);
    pthread_join(tid1,NULL);
    pthread_join(tid2,NULL);
    pthread_join(tid3,NULL);
    pthread_join(tid4,NULL);
    /* Write the image from memory to image.ppm */
    write_image("image.ppm", image, image_width, image_height);
    free(image);
    free(ptr1);free(ptr2);free(ptr3);free(ptr4);

    return 0;
}