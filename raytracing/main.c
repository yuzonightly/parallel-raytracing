/**
 * Copyright (c) 2020, Evgeniy Morozov
 * All rights reserved.
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 * mpicc test.c -o output && mpirun -np 4 output
 */
// 44922 160 177 194
#include <stdio.h>

#include "hittables/rt_hittable_list.h"
#include "rt_camera.h"
#include "rt_skybox_simple.h"
#include <errno.h>
#include <string.h>
#include <scenes/rt_scenes.h>
#include <assert.h>
#include <mpi.h>

static void show_usage(const char *program_name, int err);

static colour_t ray_colour(const ray_t *ray, const rt_hittable_list_t *list, rt_skybox_t *skybox, int child_rays)
{
    if (child_rays <= 0)
    {
        return colour(0, 0, 0);
    }

    rt_hit_record_t record;
    if (!rt_hittable_list_hit_test(list, ray, 0.001, INFINITY, &record))
    {
        return rt_skybox_value(skybox, ray);
    }

    ray_t scattered;
    colour_t attenuation;
    colour_t emitted = rt_material_emit(record.material, record.u, record.v, &record.p);
    if (rt_material_scatter(record.material, ray, &record, &attenuation, &scattered))
    {
        return vec3_sum(emitted, vec3_multiply(attenuation, ray_colour(&scattered, list, skybox, child_rays - 1)));
    }
    return emitted;
}

int main(int argc, char *argv[])
{

    // MPI
    int PROCESS_RANK, SIZE_OF_CLUSTER;
    int ROOT = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &SIZE_OF_CLUSTER);
    MPI_Comm_rank(MPI_COMM_WORLD, &PROCESS_RANK);

    // Constants
    const double ASPECT_RATIO = 3.0 / 2.0;
    const int IMAGE_WIDTH = 300;
    const int IMAGE_HEIGHT = (int)(IMAGE_WIDTH / ASPECT_RATIO);
    const int CHILD_RAYS = 50;

    // ROOT parameters
    FILE *out_file = stdout;

    // Parameters
    long number_of_samples;
    // World
    rt_hittable_list_t *world = NULL;
    rt_skybox_t *skybox = NULL;

    const char *number_of_samples_str = NULL;
    const char *scene_id_str = NULL;
    const char *file_name = NULL;
    bool verbose = false;

    // Parse console arguments
    for (int i = 1; i < argc; ++i)
    {
        if (0 == strcmp(argv[i], "-s") || 0 == strcmp(argv[i], "--samples"))
        {
            if (i + 1 >= argc)
            {
                fprintf(stderr, "Fatal error: Argument '%s' doesn't have a value\n", argv[i]);
                show_usage(argv[0], EXIT_FAILURE);
            }
            number_of_samples_str = argv[++i];
            continue;
        }
        else if (0 == strcmp(argv[i], "--scene"))
        {
            if (i + 1 >= argc)
            {
                fprintf(stderr, "Fatal error: Argument '%s' doesn't have a value\n", argv[i]);
                show_usage(argv[0], EXIT_FAILURE);
            }
            scene_id_str = argv[++i];
            continue;
        }
        else if (0 == strcmp(argv[i], "-v") || 0 == strcmp(argv[i], "--verbose"))
        {
            verbose = true;
        }
        else if (0 == strcmp(argv[i], "-h"))
        {
            show_usage(argv[0], EXIT_SUCCESS);
        }
        else if ('-' == *argv[i])
        {
            fprintf(stderr, "Fatal error: Unknown argument '%s'\n", argv[i]);
            show_usage(argv[0], EXIT_FAILURE);
        }
        else if (i + 1 < argc)
        {
            fprintf(stderr, "Fatal error: Too many positional arguments (1 expected)\n");
            show_usage(argv[0], EXIT_FAILURE);
        }
        file_name = argv[i];
    }

    if (verbose)
    {
        fprintf(stderr, "Non-parsed parameters:\n");
        fprintf(stderr, "\t- number of samples: %s\n", number_of_samples_str);
        fprintf(stderr, "\t- scene ID:          %s\n", scene_id_str);
        fprintf(stderr, "\t- file_name:         %s\n", file_name);
    }

    // Parse resulting parameters
    number_of_samples = 1000;
    if (NULL != number_of_samples_str)
    {
        char *end_ptr = NULL;
        number_of_samples = strtol(number_of_samples_str, &end_ptr, 10);
        if (*end_ptr != '\0')
        {
            fprintf(stderr, "Fatal error: Value of 'samples' is not a correct number\n");
            show_usage(argv[0], EXIT_FAILURE);
        }
    }
    rt_scene_id_t scene_id = RT_SCENE_SHOWCASE;
    if (NULL != scene_id_str)
    {
        scene_id = rt_scene_get_id_by_name(scene_id_str);
        if (RT_SCENE_NONE == scene_id)
        {
            fprintf(stderr, "Fatal error: Invalid scene identifier\n");
            show_usage(argv[0], EXIT_FAILURE);
        }
    }

    if (verbose)
    {
        fprintf(stderr, "Parsed parameters:\n");
        fprintf(stderr, "\t- number of samples: %ld\n", number_of_samples);
        fprintf(stderr, "\t- scene ID:          %d\n", scene_id);
        fprintf(stderr, "\t- file_name:         %s\n", file_name);
    }

    // Declare Camera parameters
    point3_t look_from, look_at;
    vec3_t up = point3(0, 1, 0);
    double focus_distance = 10.0, aperture = 0.0, vertical_fov = 40.0;

    // Select a scene from a pre-defined one
    switch (scene_id)
    {
        case RT_SCENE_RANDOM:
            look_from = point3(13, 2, 3);
            look_at = point3(0, 0, 0);
            aperture = 0.1;
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_random();
            break;

        case RT_SCENE_TWO_SPHERES:
            look_from = point3(13, 2, 3);
            look_at = point3(0, 0, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_two_spheres();
            break;

        case RT_SCENE_TWO_PERLIN_SPHERES:
            look_from = point3(13, 2, 3);
            look_at = point3(0, 0, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_two_perlin_spheres();
            break;

        case RT_SCENE_EARTH:
            look_from = point3(13, 2, 3);
            look_at = point3(0, 0, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_earth();
            break;

        case RT_SCENE_LIGHT_SAMPLE:
            look_from = point3(26, 3, 6);
            look_at = point3(0, 2, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_background(colour(0, 0, 0));
            world = rt_scene_light_sample();
            break;

        case RT_SCENE_CORNELL_BOX:
            look_from = point3(278, 278, -800);
            look_at = point3(278, 278, 0);

            skybox = rt_skybox_new_background(colour(0, 0, 0));
            world = rt_scene_cornell_box();
            break;

        case RT_SCENE_INSTANCE_TEST:
            look_from = point3(0, 5, -20);
            look_at = point3(0, 0, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_instance_test();
            break;

        case RT_SCENE_CORNELL_SMOKE:
            look_from = point3(278, 278, -800);
            look_at = point3(278, 278, 0);

            skybox = rt_skybox_new_background(colour(0, 0, 0));
            world = rt_scene_cornell_box_smoke_boxes();
            break;

        case RT_SCENE_SHOWCASE:
            look_from = point3(478, 278, -600);
            look_at = point3(278, 278, 0);

            skybox = rt_skybox_new_background(colour(0, 0, 0));
            world = rt_scene_showcase();
            break;

        case RT_SCENE_METAL_TEST:
            look_from = point3(0, 5, -10);
            look_at = point3(0, 2, 0);
            vertical_fov = 20.0;

            skybox = rt_skybox_new_gradient(colour(1, 1, 1), colour(0.5, 0.7, 1));
            world = rt_scene_metal_test();
            break;
        case RT_SCENE_NONE:
            fprintf(stderr, "Fatal error: scene id is undefined after parsing the parameters\n");
            return EXIT_FAILURE;
    }

    rt_camera_t *camera =
        rt_camera_new(look_from, look_at, up, vertical_fov, ASPECT_RATIO, aperture, focus_distance, 0.0, 1.0);

    int *range, *slices;
    int total_image_size = IMAGE_HEIGHT * IMAGE_WIDTH;
    // Used in MPI_Gatherv
    int *rcvcounts;
    rcvcounts = (int *)(malloc(SIZE_OF_CLUSTER * sizeof(int)));
    // Used in MPI_Gatherv
    int *displays;
    displays = (int *)(malloc(SIZE_OF_CLUSTER * sizeof(int)));
    if (PROCESS_RANK == ROOT)
    {
        int slice = total_image_size / SIZE_OF_CLUSTER;
        // We use remainder to better distribute the ranges
        int remainder = total_image_size % SIZE_OF_CLUSTER;

        // Create array with ranges.
        slices = (int *)malloc((2 * SIZE_OF_CLUSTER) * sizeof(int));
        int current = 0;
        for (int i = 0; i < SIZE_OF_CLUSTER; i++)
        {
            displays[i] = current * 3;
            if (i != SIZE_OF_CLUSTER - 1)
            {
                slices[i * 2 + 0] = current;
                // This !(!(remainder ^ 0)) compares remainder and 0
                // if != returns 1, if not 0
                current = current + slice + !(!(remainder ^ 0));
                slices[i * 2 + 1] = current;
            }
            else
            {
                slices[i * 2 + 0] = current;
                slices[i * 2 + 1] = total_image_size;
            }
            if (remainder > 0)
            {
                remainder--;
            }
            // Multiply by 3 -> RGB.
            rcvcounts[i] = (slices[i * 2 + 1] - slices[i * 2 + 0]) * 3;
        }
    }

    // range is a tuple containing [initial_index, final_index]
    range = (int *)malloc(2 * sizeof(int));
    // Scatter the ranges
    MPI_Scatter(slices, 2, MPI_INT, range, 2, MPI_INT, ROOT, MPI_COMM_WORLD);
    printf("\n%d, %d\n", range[0], range[1]);

    // Since we have triples (RGB), multiply local_range * 3
    int *partial_buffer;
    int local_range = range[1] - range[0];
    int partial_buffer_size = 3 * local_range;
    partial_buffer = (int *)(malloc(partial_buffer_size * sizeof(int)));
    // 0 to 59999
    int buffer_offset = 0;
    for (int r = range[0]; r < range[1]; r++)
    {
        // Extract i and j from r
        int j = (int)(r / IMAGE_WIDTH);
        int i = r % IMAGE_WIDTH;

        colour_t pixel = colour(0, 0, 0);
        for (int s = 0; s < number_of_samples; ++s)
        {
            double u = (double)(i + rt_random_double(0, 1)) / (IMAGE_WIDTH - 1);
            double v = (double)(j + rt_random_double(0, 1)) / (IMAGE_HEIGHT - 1);

            ray_t ray = rt_camera_get_ray(camera, u, v);
            vec3_add(&pixel, ray_colour(&ray, world, skybox, CHILD_RAYS));
        }
        // Add result to partial_buffer
        rt_write_colour(partial_buffer, buffer_offset, pixel, number_of_samples);
        buffer_offset = buffer_offset + 1;
    }

    int complete_buffer_size = total_image_size * 3;
    int *complete_buffer = (int *)(malloc(complete_buffer_size * sizeof(int)));
    // Gather partial_buffer from all processes
    // printf("\n%d\n", rcvcounts[0]);
    // printf("\n%d\n", rcvcounts[1]);
    // printf("\n%d\n", displays[0]);
    // printf("\n%d\n", displays[1]);
    MPI_Gatherv(partial_buffer, partial_buffer_size, MPI_INT, complete_buffer, rcvcounts, displays, MPI_INT, ROOT,
                MPI_COMM_WORLD);

    if (PROCESS_RANK == ROOT)
    {
        if (NULL != file_name)
        {
            out_file = fopen(file_name, "w");
            if (NULL == out_file)
            {
                fprintf(stderr, "Fatal error: Unable to open file %s: %s", file_name, strerror(errno));
                goto cleanup;
            }
        }

        // Write final results to file
        fprintf(out_file, "P3\n%d %d\n255\n", IMAGE_WIDTH, IMAGE_HEIGHT);
        // Since we execute differently, we changed the way the result is written
        for (int i = IMAGE_HEIGHT - 1; i >= 0; i--)
        {
            // get row
            int idx = i * IMAGE_WIDTH;
            for (int j = 0; j < IMAGE_WIDTH; j++)
            {
                // get equivalent index in complete_buffer
                // Store row by row starting from last
                fprintf(out_file, "%d %d %d\n", complete_buffer[(idx + j) * 3 + 0], complete_buffer[(idx + j) * 3 + 1],
                        complete_buffer[(idx + j) * 3 + 2]);
            }
        }
    }

    fprintf(stderr, "\nDone\n");
cleanup:
    // Cleanup
    rt_hittable_list_deinit(world);
    rt_camera_delete(camera);
    rt_skybox_delete(skybox);

    MPI_Finalize();
    return EXIT_SUCCESS;
}

static void show_usage(const char *program_name, int err)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "%s [-s|--samples N] [--scene SCENE] [-v|--verbose] [output_file_name]\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-s | --samples      <int>       Number of rays to cast for each pixel\n");
    fprintf(
        stderr,
        "\t--scene             <string>    ID of the scene to render. List of available scenes is printed below.\n");
    fprintf(stderr, "\t-v | --verbose                  Enable verbose output\n");
    fprintf(stderr, "\t-h                              Show this message and exit\n");
    fprintf(stderr, "Positional arguments:\n");
    fprintf(stderr,
            "\toutput_file_name                Name of the output file. Outputs image to console if not specified.\n");
    fprintf(stderr, "Available scenes:\n");
    rt_scene_print_scenes_info(stderr);

    exit(err);
}
