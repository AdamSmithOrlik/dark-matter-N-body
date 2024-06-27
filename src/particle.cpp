#include "particle.h"
#include <algorithm> // For std::fill
#include <iostream>

Particle::Particle(double mass, array<double, 2> position, array<double, 2> velocity)
    : mass(mass), position(position), velocity(velocity), cross_section(10.0)
{
    fill(force.begin(), force.end(), 0.0);
    fill(acceleration.begin(), acceleration.end(), 0.0);
}

double Particle::scattering_probability(const Particle &other) const
{
    double cross_section = 10;
    return cross_section;
}

void Particle::print() const
{
    std::cout << "\nPosition: (" << position[0] << ", " << position[1] << ")"
              << "\nVelocity: (" << velocity[0] << ", " << velocity[1] << ")"
              << std::endl;
}