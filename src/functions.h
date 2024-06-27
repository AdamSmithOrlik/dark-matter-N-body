#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <functional>
#include <cmath>
#include <random>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include "particle.h"
#include <iostream>
#include <gsl/gsl_integration.h>

// Constants
const double PI = 3.14159265358979323846;

// Maxwell-Boltzmann PDF
double maxwell_boltzmann_pdf(double v, double mass = 1.0, double temperature = 1.0)
{
    const double kB = 1.0; // Boltzmann constant in suitable units
    double factor = std::sqrt((2.0 * mass) / (PI * kB * temperature));
    return factor * v * v * std::exp(-mass * v * v / (2.0 * kB * temperature));
}

// Integration function using GSL
double integrate(std::function<double(double)> func, double a, double b, size_t limit = 1000)
{
    gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(limit);
    double result, error;

    auto gsl_func = [](double x, void *params) -> double
    {
        return static_cast<std::function<double(double)> *>(params)->operator()(x);
    };

    gsl_function F;
    F.function = gsl_func;
    F.params = &func;

    gsl_integration_qags(&F, a, b, 0, 1e-7, limit, workspace, &result, &error);

    gsl_integration_workspace_free(workspace);
    return result;
}

// Compute the CDF for Maxwell-Boltzmann distribution
std::pair<std::vector<double>, std::vector<double>> compute_maxwell_boltzmann_cdf(double mass = 1.0, double temperature = 1.0, double v_max = 10.0, int num_points = 1000)
{
    std::vector<double> v(num_points);
    std::vector<double> cdf(num_points, 0.0);

    double dv = v_max / (num_points - 1);
    for (int i = 0; i < num_points; ++i)
    {
        v[i] = i * dv;
    }

    for (int i = 1; i < num_points; ++i)
    {
        cdf[i] = integrate([=](double x)
                           { return maxwell_boltzmann_pdf(x, mass, temperature); }, 0.0, v[i]);
    }

    // Normalize the CDF
    for (double &val : cdf)
    {
        val /= cdf[num_points - 1];
    }

    return std::make_pair(v, cdf);
}

// Invert the CDF using interpolation
std::function<double(double)> invert_cdf(const std::vector<double> &v, const std::vector<double> &cdf)
{
    if (v.size() != cdf.size())
    {
        throw std::invalid_argument("v and cdf must be the same size");
    }

    return [v, cdf](double p)
    {
        if (p <= 0.0)
        {
            return v.front();
        }
        if (p >= 1.0)
        {
            return v.back();
        }

        auto it = std::lower_bound(cdf.begin(), cdf.end(), p);
        size_t idx = std::distance(cdf.begin(), it);
        if (idx == 0)
        {
            return v[0];
        }
        double t = (p - cdf[idx - 1]) / (cdf[idx] - cdf[idx - 1]);
        return v[idx - 1] + t * (v[idx] - v[idx - 1]);
    };
}

// Maxwell-Boltzmann inverse CDF function
std::function<double(double)> maxwell_boltzmann(double mass = 1.0, double temperature = 1.0)
{
    auto [v, cdf] = compute_maxwell_boltzmann_cdf(mass, temperature);
    return invert_cdf(v, cdf);
}

// NFW Density function
double nfw_density(double r, double rho_0 = 1.0, double rs = 1.0)
{
    double x = r / rs;
    double result = rho_0 / (x * (1 + x) * (1 + x));
    // std::cout << "nfw_density(" << r << ") = " << result << std::endl; // Debug print
    return result;
}

// Compute the CDF for a given density function
std::pair<std::vector<double>, std::vector<double>> compute_cdf(std::function<double(double)> density_func, double rmax, int num_points = 1000)
{
    std::vector<double> r(num_points);
    std::vector<double> cdf(num_points, 0.0);

    double dr = rmax / (num_points - 1);
    for (int i = 0; i < num_points; ++i)
    {
        r[i] = i * dr;
    }

    for (int i = 1; i < num_points; ++i)
    {
        cdf[i] = integrate([=](double x)
                           { return 4 * PI * x * x * density_func(x); }, 0.0, r[i]);
    }

    // Normalize the CDF
    for (double &val : cdf)
    {
        val /= cdf[num_points - 1];
    }

    // // Debug print
    // for (int i = 0; i < num_points; ++i)
    // {
    //     std::cout << "r[" << i << "] = " << r[i] << ", cdf[" << i << "] = " << cdf[i] << std::endl;
    // }

    return std::make_pair(r, cdf);
}

// Generate starting positions based on a given density function
std::vector<std::array<double, 2>> starting_positions(std::function<double(double)> density_func, double rmax = 100, int num_particles = 100)
{
    auto [r, cdf] = compute_cdf(density_func, rmax);
    auto inverse_cdf = invert_cdf(r, cdf);

    // Sample uniform random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> uniform_samples(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        uniform_samples[i] = dis(gen);
    }

    // Use the inverse CDF to get the particle positions
    std::vector<double> radii(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        radii[i] = inverse_cdf(uniform_samples[i]);
        // std::cout << "radii[" << i << "] = " << radii[i] << std::endl;
    }

    // Convert spherical coordinates to Cartesian coordinates
    std::uniform_real_distribution<> angle_dis(0.0, 2 * PI);
    std::vector<std::array<double, 2>> positions(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        double theta = angle_dis(gen);
        positions[i] = {radii[i] * std::cos(theta), radii[i] * std::sin(theta)};
        // std::cout << "positions[" << i << "] = (" << positions[i][0] << ", " << positions[i][1] << ")" << std::endl;
    }

    return positions;
}

// Generate starting velocities based on a given velocity distribution function
std::vector<std::array<double, 2>> starting_velocities(std::function<double(double)> velocity_func, int num_particles = 100)
{
    // Sample uniform random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> uniform_samples(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        uniform_samples[i] = dis(gen);
    }

    // Use the inverse CDF to get the particle speeds
    std::vector<double> speeds(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        speeds[i] = velocity_func(uniform_samples[i]);
    }

    // Assign random directions
    std::uniform_real_distribution<> angle_dis(0.0, 2 * PI);
    std::vector<std::array<double, 2>> velocities(num_particles);
    for (int i = 0; i < num_particles; ++i)
    {
        double angle = angle_dis(gen);
        velocities[i] = {speeds[i] * std::cos(angle), speeds[i] * std::sin(angle)};
    }

    return velocities;
}

// Initialize particles
std::vector<Particle> initialize_particles(int num_particles = 100)
{
    auto density_func = [](double r)
    { return nfw_density(r); };
    auto positions = starting_positions(density_func, 100, num_particles);
    auto maxwell_func = maxwell_boltzmann(1.0, 1.0);
    auto velocities = starting_velocities(maxwell_func, num_particles);

    std::vector<Particle> particles;
    for (int i = 0; i < num_particles; ++i)
    {
        particles.emplace_back(1.0, positions[i], velocities[i]);
    }

    return particles;
}

#endif // FUNCTIONS_H