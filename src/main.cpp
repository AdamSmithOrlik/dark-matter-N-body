#include <iostream>
#include "functions.h"
#include "particle.h"
#include "constants.h"

using namespace std;

int main()
{
    // // Testing the Particle class
    // array<double, 2> position = {1.0, 2.0};
    // array<double, 2> velocity = {0.5, -0.5};
    // Particle p1(1.0, position, velocity);

    // array<double, 2> other_position = {-1.0, -2.0};
    // array<double, 2> other_velocity = {-0.5, 0.5};
    // Particle p2(1.0, other_position, other_velocity);

    // double probability = p1.scattering_probability(p2);
    // cout << "Scattering probability: " << probability << endl;
    // // Print Particle 1 properties
    // cout << "Particle 1 properties:" << endl;
    // cout << "Mass: " << p1.mass << endl;
    // cout << "Position: (" << p1.position[0] << ", " << p1.position[1] << ")" << endl;
    // cout << "Velocity: (" << p1.velocity[0] << ", " << p1.velocity[1] << ")" << endl;
    // cout << "Force: (" << p1.force[0] << ", " << p1.force[1] << ")" << endl;
    // cout << "Acceleration: (" << p1.acceleration[0] << ", " << p1.acceleration[1] << ")" << endl;
    // cout << "Cross Section: " << p1.scattering_probability(p2) << endl;
    // cout << endl;

    // // Print Particle 2 properties
    // cout << "Particle 2 properties:" << endl;
    // cout << "Mass: " << p2.mass << endl;
    // cout << "Position: (" << p2.position[0] << ", " << p2.position[1] << ")" << endl;
    // cout << "Velocity: (" << p2.velocity[0] << ", " << p2.velocity[1] << ")" << endl;
    // cout << "Force: (" << p2.force[0] << ", " << p2.force[1] << ")" << endl;
    // cout << "Acceleration: (" << p2.acceleration[0] << ", " << p2.acceleration[1] << ")" << endl;
    // cout << "Cross Section: " << p2.scattering_probability(p1) << endl;
    // cout << endl;
    int num_particles = 1000000;
    auto particles = initialize_particles(num_particles);

    // for (int i; i < num_particles; i++)
    // {
    //     auto particle = particles[i];
    //     // only print every 100th element
    //     if (i % 100 == 0)
    //     {
    //         particle.print();
    //     }
    // }

    for (const auto &particle : particles)
    {
        // particle.print();
    }

    return 0;
}