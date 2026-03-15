# Mount Charles — A Particle Simulator

This project is a custom GEANT4-style particle propagation simulation.

# Conventions
Mirroring GEANT4, the simulation uses millimetres (mm), nanoseconds (ns) and mega-electronvolts (MeV) as the base units for length, time and energy respectively.

Natural units are used for all the electronvolt-related quantities, so the units of mass, energy and momentum are all MeV. 

# Structure

## World
The top-level object in the simulation is the world. This stores:
1. Global time: the absolute time of the simulation, which starts at 0ns when the particle(s) is produced.
2. Time step: the time step used in the simulation. Each simulation step will advance the global time by this time step.
3. List of particles: a list of all the particles in the event. These are kept in the list even when they are no longer being simulated (due to exiting the volume, decaying, etc.).
4. Volume: the simulation volume in which the particles are contained. The particles are killed upon exiting this volume.
5. RNG: the random number generator instance used for generating pseudo-random numbers in all functions that need them (Klein-Nishina sampling, MCS deflection angle, etc.). This single top-level instance is used so that simulations are reproducible, determined by a single seed value.

## Volume
The simulation volume is a cube centred on the origin and characterized by a single `size` parameter: this corresponds to the edge length of the simulation cube. Particle [interaction](#interactions) and [propagation](#propagation) is only calculated inside this volume. For simplicity, the volume is taken to be made of liquid water.

## Particle
A particle is made of these components:
1. Particle type: currently one of electron ($e^-$), muon ($\mu^-$) or gamma ($\gamma$).
2. Particle state: this describes the particle's properties, namely
   1. Position: the 3D position of the particle (mm).
   2. Momentum: the momentum of the particle (MeV).
   3. Mass: the mass of the particle (MeV).
   4. Alive: whether the particle is considered "alive" or not, i.e. whether it is still being simulated.
3. Particle track: this stores the position, time and energy deposited (can be None) for each point in the particle's track.
4. Distance to next interaction: this is calculated at every interaction, after which the particle will be propagated but not interacted until the next time this variable drops back to 0. This is used for regulating the discrete nature of gamma interactions. For electrons and muons, this is always set to 0.
   
# Physics processes

## Propagation
For every time step, every active particle in the simulation is propagated with the small position update

$$\Delta\vec{r} = \hat{p}\cdot\beta c\cdot \Delta t$$

where $\vec{r}$ is the particle position, $\hat{p}$ its normalized momentum vector, $\beta$ the speed parameter of the particle ($v/c$), $c$ the speed of light (in appropriate units, mm/ns) and $\Delta t$ the time step (in ns).

## Interactions
As this is a very simple particle simulator, only the dominant interaction for each particle type is implemented.

### Electron
Electron interactions are treated in two parts: energy loss and scattering. These are essentially the manifestations of the change in momentum at each time step: energy loss corresponds to the change in the momentum magnitude, while scattering corresponds to the change in the momentum direction.

#### Energy loss
Energy loss for ionizing particles is governed by the stopping power of the target material for that particle, notated dE/dx, which varies as a function of particle momentum. For electrons in this simulation, the [NIST ESTAR database](https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html) is used, with the material set to liquid water. This yields a list of discrete kinetic energy values with their corresponding dE/dx values. In order to obtain an analytical form that can be computed quickly during the simulation, a function (see [Table of coefficients](#table-of-coefficients)) is fit to the database values. Then, at every time step, the energy of the particle is decreased by

$$\Delta E = \frac{dE}{dx}\cdot dx = \frac{dE}{dx}\cdot\beta c\cdot\Delta t$$

Once the electron momentum drops below 0.103 MeV, it is out of the tabulated range, so a constant 8 MeV/cm energy loss is assumed.

#### Scattering
The dominant process in electron scattering (at MeV energies) is multiple Coulomb scattering (MCS). This process can be approximated using the Highland formula (given in this [PDG review](https://pdg.lbl.gov/2019/reviews/rpp2018-rev-passage-particles-matter.pdf) for example), which gives the RMS scattering angle for the electron for a small step $dx$:

$$\theta_0 = \frac{13.6\text{ MeV}}{\beta p}\sqrt{\frac{dx}{X_0}}\left[1 + 0.038\ln{\left(\frac{dx}{X_0}\right)}\right]$$

where $p$ is the momentum of the electron, $\beta$ the speed parameter and $X_0$ the radiation length of the material. In this formula, the $\frac{z^2}{\beta^2}$ term has been omitted from the logarithm, as it is considered to be negligible. In this simulation, the value for water of $X_0=36.08\text{ cm}$ is used.

### Muon
Muons are treated very similarly to electrons. The two components of the interaction are described below.

#### Energy loss
Energy loss for muons is separated into two cases: momentum above and below 50MeV ("high" and "low" momentum).

The dE/dx curve for the high momentum range is computed from the tabulated values in the full version of [this paper](https://pdg.lbl.gov/2023/AtomicNuclearProperties/adndt.pdf). These values are again fit to a degree 8 "log polynomial" (see [Table of coefficients](#table-of-coefficients)) in order to interpolate between the discrete given values.

For the low momentum range, the Bethe-Bloch equation is used:

$$-\frac{dE}{dx} = K\frac{Z_{\text{eff}}}{A}\frac{1}{\beta^2}\left[\frac{1}{2}\ln{\left(\frac{2m_e \beta^2\gamma^2 T_{\text{max}}}{I^2}\right)}-\beta^2-\frac{C}{Z}-\delta\right]$$

where $Z_{\text{eff}}$ and $A$ are the effective atomic number and mass of the target, $T_{\text{max}}$ the maximum kinetic energy transfer to an electron, $I$ the mean excitation energy, $\frac{C}{Z}$ the shell correction and $\delta$ the density effect correction which we take to be negligible. Further, $T_{\text{max}}$ is given by (in natural units)

$$T_{\text{max}} = \frac{2m_e\beta^2\gamma^2}{1 + 2\gamma m_e/m_\mu + (m_e/m_\mu)^2}$$

whereas the shell correction is given by (from [this paper](https://pdg.lbl.gov/2023/AtomicNuclearProperties/adndt.pdf))

$$\frac{C}{Z} \approx 0.42237\cdot\left(\frac{1}{\beta\gamma}\right)^2 + 0.0304\cdot\left(\frac{1}{\beta\gamma}\right)^4$$

and the values of the other parameters for liquid water are shown in this table:

| Parameter        | Value    |
| :--------------- | -------: |
| $Z_{\text{eff}}$ | 10       |
| $A$              | 18 g/mol |
| $I$              | 75 eV    |

As for electrons, once the muon momentum drops below 8.9 MeV, it is out of the tabulated range, so a constant 8 MeV/cm energy loss is assumed.

#### Scattering
Muons are scattering identically to electrons. Generally, their higher momentum will mean that their tracks are "straighter".

### Gamma
Gammas are fundamentally different in that they interact in discrete events. The two processes considered for gammas in this simulation are photoelectric absorption — which fully absorbs the gamma — and Compton scattering — also called incoherent scattering, which scatters gammas while reducing their energy. Note the absence of pair production as one of the simulated processes.

#### Interactions
The gamma interactions are computed with the following steps:
1. Obtain the total gamma mass attenuation coefficients, by adding individual interaction coefficients given in the [NIST XCOM database](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html) (fit with log polynomial and power law functions, see [Table of coefficients](#table-of-coefficients)):

   $$\mu_{\text{total}}(E)=\mu_{\text{photo}}(E) + \mu_{\text{Compton}}(E)$$

   Convert this to a total macroscopic attenuation coefficient by multiplying by the material density (for liquid water, this is a factor of 1).

2. Calculate the mean free path:

   $$\lambda(E)=\frac{1}{\mu_{\text{total}}(E)}$$

   This corresponds to the mean distance between each interaction, the distribution of which is given by

   $$P(x)=\frac{1}{\lambda}e^{-\lambda x}$$

3. After each interaction, calculate the distance to the next interaction using

   $$d=-\lambda\ln{(\alpha)}$$

   where $\alpha$ is a random number sampled from $(0, 1]$.
   
4. To decide which interaction happens at the interaction point, the probability for each can be calculated simply with 

   $$P_{\text{photo}}=\frac{\mu_{\text{photo}}}{\mu_{\text{total}}}\hspace{2cm}P_{\text{Compton}}=\frac{\mu_{\text{Compton}}}{\mu_{\text{total}}}$$

#### Energy loss
To compute the energy loss of the gamma, we treat the interactions separately. Photoelectric absorption is trivial as the gamma simply loses all of its energy. For Compton scattering the energy loss is computed via the scattering angle, which is done with Klein-Nishina sampling using a rejection algorithm. For an incoming gamma energy $E$:
1. Sample $\mu=\cos{(\theta)}$ uniformly between -1 and 1.
2. Calculate the energy ratio

   $$\frac{E'}{E}=\frac{1}{1 + \alpha(1-\mu)}$$
   
   with $\alpha=\frac{E}{m_e}$ where $m_e$ is the mass of the electron.
   
3. Compute the Klein-Nishina weight

   $$f(\mu)=\left(\frac{E'}{E}\right)^2 \left[\frac{E'}{E}+\frac{E}{E'}-(1-\mu^2)\right]$$
   
4. Accept the candidate $\mu$ with probability $P=\frac{f(\mu)}{f_{\text{max}}}$, where $f_{\text{max}}$ is set to 2.0 for this simulation. This variable essentially decides how stringent the rejection algorithm is: too low values provide less precise results, too high values mean more iterations before a candidate is accepted. If the candidate is rejected, repeat the process with a new $\mu$ candidate. Otherwise, extract the angle with $\theta=\arccos{(\mu)}$.

Once the angle is found for the scatter, the energy loss is calculated with

$$\Delta E = E-E' = E\left(1 - \frac{1}{1 + \frac{E}{m_e}(1-\cos{(\theta)})}\right)$$

### Table of coefficients
The function that is used to recreate the dE/dx curves for electrons and muons, as well as the Compton scattering cross-section is the so-called "log polynomial" of degree $D$, given by

$$\frac{dE}{dx}(p)=\sum_{n=0}^D c_n\ln{(p)}^n$$

For this simulation, the log polynomial of degree 8 is used. The optimal coefficients for each particle are shown in the table below.

| Coefficient | Electrons | Muons (<50 MeV) | Muons (>50 MeV) | Gammas (Compton) |
| :---------: | --------: | --------------: | --------------: | ----------------: |
| $c_0$ | $1.97185875\cdot 10^0$     | $-2.21192313\cdot 10^5$ | $1.13754387\cdot 10^3$     | $7.05838611\cdot 10^{-2}$  | 
| $c_1$ | $-4.90322067\cdot 10^{-1}$ | $4.16349323\cdot 10^5$  | $-1.13642381\cdot 10^3$    | $-3.55852266\cdot 10^{-2}$ |
| $c_2$ | $5.67984147\cdot 10^{-1}$  | $-3.02334049\cdot 10^5$ | $4.96588219\cdot 10^2$     | $4.99571928\cdot 10^{-3}$  |
| $c_3$ | $-3.78515229\cdot 10^{-1}$ | $9.22330794\cdot 10^4$  | $-1.23563655\cdot 10^2$    | $7.02986758\cdot 10^{-4}$  |
| $c_4$ | $1.96937857\cdot 10^{-1}$  | $1.78846389\cdot 10^3$  | $1.91190645\cdot 10^1$     | $-2.81938080\cdot 10^{-4}$ |
| $c_5$ | $-6.69875048\cdot 10^{-2}$ | $-9.81957228\cdot 10^3$ | $-1.88126582\cdot 10^0$    | $1.60699821\cdot 10^{-5}$  |
| $c_6$ | $1.30714285\cdot 10^{-2}$  | $2.97223872\cdot 10^3$  | $1.14850292\cdot 10^{-1}$  | $2.97942065\cdot 10^{-6}$  |
| $c_7$ | $-1.31646064\cdot 10^{-3}$ | $-3.90203242\cdot 10^2$ | $-3.97495919\cdot 10^{-3}$ | $-3.27163146\cdot 10^{-7}$ |
| $c_8$ | $5.29555090\cdot 10^{-5}$  | $1.99344973\cdot 10^1$  | $5.96940644\cdot 10^{-5}$  | $3.12212653\cdot 10^{-9}$  |

To recreate the photoelectric absorption cross-section for gammas, the sum of $N$ power laws is used:

$$\sigma(E)=\sum_{i=1}^N A_i E^{-p_i}$$

For this simulation, $N=2$ is used. The optimal coefficients are shown below.

| Coefficient | Gammas (photoelectric) |
| :---------: | ----------------------: |
| $A_1$ | $1.99979566\cdot 10^{-6}$ |
| $A_2$ | $1.52497758\cdot 10^{-6}$ |
| $p_1$ | $3.16002041$ |
| $p_2$ | $1.03879142$ |

# Usage
A small example code for running the simulation is shown below:

``` rust
mod utils;
mod particle;
mod geometry;
mod sim;

use std::io::{Write, stdout};

use sim::world::World;
use geometry::volume::Volume;
use particle::particle::{Particle, ParticleType};
use utils::vec3::Vec3;

fn main() {
    // Define detector volume
    let volume = Volume::new(
        1000.0,  // size (1m cube)
        360.8    // radiation length (water)
    );

    // Define particles
    let particles = vec![
        Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(5.0,   0.0,  0.0), ParticleType::Electron),
        Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 500.0,  0.0), ParticleType::Muon),
        Particle::new(Vec3(0.0, 0.0, 0.0), Vec3(0.0,   0.0, 10.0), ParticleType::Gamma),
    ];

    // Define world
    let mut world = World::new(
        particles,
        volume,
        0.01,      // dt
        42         // random seed
    );

    // Simulate
    let mut step = 0;
    while world.has_alive_particles() {
        // Step the simulation
        world.step();
        step += 1;

        print!("\rStep {}", step);
        stdout().flush().unwrap();

        if step > 1000 {
            println!("Max steps reached!");
            break;
        }
    }
    println!("\nFinished after {} steps\n\n", step);

    for track in world.tracks().iter() {
        let particle_name = match track.particle_type{
            ParticleType::Electron => "e-",
            ParticleType::Muon     => "mu-",
            ParticleType::Gamma    => "gamma",
        };
        println!("Particle {}:", particle_name);
        for point in track.points.iter() {
            println!("   Position: ({:.2}, {:.2}, {:.2})     Time: {:.2}     Energy deposited: {:.2?}",
                     point.r.0,
                     point.r.1,
                     point.r.2,
                     point.t,
                     point.E
            );
        }
        println!("");
    }
}
```

