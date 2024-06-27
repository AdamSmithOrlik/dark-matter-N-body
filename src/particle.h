#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>
using namespace std;

class Particle
{
public:
    Particle(double mass, array<double, 2> position, array<double, 2> velocity);

    double mass;
    array<double, 2> position;
    array<double, 2> velocity;
    array<double, 2> force;
    array<double, 2> acceleration;

    double scattering_probability(const Particle &other) const;

    void print() const;

private:
    double cross_section;
};

#endif // PARTICLE_H