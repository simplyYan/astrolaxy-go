package astrolaxy

import (
	"fmt"
	"math"
)

//Use the _help function to get an explanation of all formulas. It explains the basics of how to use them and what they are for.
func _help() {
	fmt.Println(`
	Drake equation:
         N = R* × fp × ne × fl × fi × fc × L
         Where:
             N is the number of civilizations in our galaxy with which we could communicate.
             R* is the rate of star formation in the galaxy.
             fp is the fraction of stars that have planets.
             ne is the average number of planets that can potentially support life per star that has planets.
             fl is the fraction of these planets that actually develop life.
             fi is the fraction of planets that develop intelligent life.
             fc is the fraction of civilizations that develop technology capable of being detected.
             L is the average time these civilizations emit detectable signals into space.

     Newton's Law of Universal Gravitation:
         F = G × (m₁ × m₂) / r²
         Where:
             F is the gravitational force between two masses.
             G is the gravitational constant.
             m₁ and m₂ are the masses of the two bodies.
             r is the distance between the centers of masses.

     Coulomb's Law:
         F = k × (|q₁ × q₂|) / r²
         Where:
             F is the electrical force between two charges.
             k is the electrostatic constant.
             q₁ and q₂ are the charges of the two particles.
             r is the distance between the charges.

     Newton's Second Law (Law of Force or Fundamental Law of Dynamics):
         F = m × a
         Where:
             F is the applied force.
             m is the mass of the object.
             a is the acceleration produced by the force.

     Snell-Descartes Law:
         n₁ × sin(θ₁) = n₂ × sin(θ₂)
         Where:
             n₁ and n₂ are the refractive indices of media 1 and 2, respectively.
             θ₁ and θ₂ are the angles of incidence and refraction, respectively.

     Schrödinger equation (time-independent version):
         ĤΨ = EΨ
         Where:
             Ĥ is the Hamiltonian operator.
             Ψ is the wave function of the system.
             E is the total energy of the system.

     Stefan-Boltzmann Law:
         P = σ × A × T⁴
         Where:
             P is the total radiant power emitted by a black body.
             σ is the Stefan-Boltzmann constant.
             A is the surface area of the blackbody.
             T is the blackbody temperature in kelvin.

     Ampere-Maxwell Law:
         ∇ × B = μ₀ × J + μ₀ × ε₀ × ∂E/∂t
         Where:
             B is the magnetic field.
             J is the electrical current density.
             And it's the electric field.
             μ₀ is the magnetic permeability of the vacuum.
             ε₀ is the electrical permittivity of the vacuum.

     Bernoulli equation:
         P₁ + (1/2)ρv₁² + ρgh₁ = P₂ + (1/2)ρv₂² + ρgh₂
         Where:
             P is the pressure.
             ρ is the density of the fluid.
             v is the fluid velocity.
             h is the height in relation to a chosen reference.
             g is the acceleration due to gravity.

     Equation of E=mc² (Mass-energy equivalence):
         E = m × c²
         Where:
             And it's energy.
             m is the mass.
             c is the speed of light in a vacuum.
            
     Faraday's Law of Electromagnetic Induction:
         ε = -dΦ/dt
         Where:
             ε is the induced electromotive force.
             Φ is the magnetic flux through a surface.
             dt is the time variation.

     Ohm's Law:
         V = I × R
         Where:
             V is the potential difference (voltage).
             I is the electric current.
             R is the electrical resistance.

     Maxwell's Equation for the Propagation of Electromagnetic Waves:
         ∇²E = με(∂²E/∂t²)
         Where:
             And it's the electric field.
             μ is the magnetic permeability.
             ε is the electrical permittivity.
             ∇² is the Laplacian operator.
             ∂²E/∂t² is the second time derivative of the electric field.

     Planck-Einstein equation:
         E = h × f
         Where:
             And it's the energy of the photon.
             h is Planck's constant.
             f is the frequency of electromagnetic radiation.

     Einstein's equation for Special Relativity (energy-quantum of motion):
         E² = (pc)² + (m₀c²)²
         Where:
             And it's the total energy.
             p is the quantity of movement.
             c is the speed of light in a vacuum.
             m₀ is the rest mass of the object.

     Navier-Stokes equation (fundamental equation of fluid dynamics):
         ρ(∂v/∂t) + ρ(v ⋅ ∇)v = -∇P + μ∇²v + ρg
         Where:
             ρ is the density of the fluid.
             v is the velocity vector.
             t is time.
             P is the pressure.
             μ is the dynamic viscosity.
             g is the acceleration due to gravity.

     Schrödinger equation (time-dependent version):
         ĤΨ = iħ(∂Ψ/∂t)
         Where:
             Ĥ is the Hamiltonian operator.
             Ψ is the wave function of the system.
             t is time.
             ħ is the reduced Planck constant.

     Equation of the Kinetic Theory of Gases (equation of state):
         PV = nRT
         Where:
             P is the pressure.
             V is the volume.
             n is the number of moles.
             R is the ideal gas constant.
             T is the temperature in Kelvin.

     Heisenberg Equation of Uncertainty:
         Δx × Δp ≥ ħ/2
         Where:
             Δx is the uncertainty in position.
             Δp is the uncertainty in linear momentum.
             ħ is the reduced Planck constant.

     Energy Conservation Law:
         ΔE = Q - W
         Where:
             ΔE is the variation in the internal energy of the system.
             Q is the heat added to the system.
             W is the work done by the system.
            
            
    Coulomb's Law Equation:
         F = k × (|q₁ × q₂|) / r²
         Where:
             F is the electrostatic force between two charges.
             k is the electrostatic constant.
             q₁ and q₂ are the magnitudes of the charges involved.
             r is the distance between the charges.

     Snell-Descartes Law for Refraction:
         n₁ × sin(θ₁) = n₂ × sin(θ₂)
         Where:
             n₁ and n₂ are the refractive indices of media 1 and 2, respectively.
             θ₁ and θ₂ are the angles of incidence and refraction, respectively.

     Hooke's Law Equation:
         F = -k × x
         Where:
             F is the applied force.
             k is the spring constant.
             x is the deformation of the spring from its equilibrium length.

     Bernoulli equation (simplified form for incompressible flow without viscosity):
         P + (1/2)ρv² + ρgh = constant
         Where:
             P is the pressure.
             ρ is the density of the fluid.
             v is the fluid velocity.
             h is the height in relation to a chosen reference.
             g is the acceleration due to gravity.

     Poisson's equation:
         ∇²φ = -ρ/ε₀
         Where:
             ∇²φ is the Laplacian operator of φ, a scalar field.
             ρ is the charge density.
             ε₀ is the electrical permittivity of the vacuum.

     Schrödinger Equation (Wave Equation):
         iħ ∂Ψ/∂t = ĤΨ
         Where:
             Ĥ is the Hamiltonian operator.
             Ψ is the wave function of the system.
             t is time.
             ħ is the reduced Planck constant.

     Continuity Equation (for fluid flow):
         A₁v₁ = A₂v₂
         Where:
             A is the cross-sectional area of the flow.
             v is the fluid velocity.

     Ampere's law (integral form):
         ∮ B ⋅ dl = μ₀I + μ₀ε₀ dΦE/dt
         Where:
             B is the magnetic field.
             dl is a length element along the closed path of integration.
             I is the electric current that crosses the surface delimited by the closed path.
             ΦE is the electric flux through the surface bounded by the closed path.

     Maxwell's Equation for Electric Field Divergence:
         ∇ ⋅ E = ρ/ε₀
         Where:
             ∇ ⋅ E is the divergence of the electric field.
             ρ is the charge density.
             ε₀ is the electrical permittivity of the vacuum.

     Linear Momentum Conservation Equation (integral form):
         ∫ F ⋅ dA = d/dt ∫ p ⋅ dA
         Where:
             F is the external force acting on a system.
             p is the linear momentum of the system.
             A is the surface that delimits the volume of the system.
            
    
    Kepler's Law: The three laws of planetary motion formulated by Johannes Kepler, the first of which is: The orbit of a planet around the Sun is an ellipse with the Sun at one of the foci.

Newton's Law of Universal Gravitation: The law formulated by Isaac Newton that describes the gravitational attraction between two massive bodies. The formula is: F = G × (m1 × m2) / r^2, where F is the gravitational force, G is the gravitational constant, m1 and m2 are the masses of the two bodies and r is the distance between them.

Escape Velocity Equation: This equation is used to calculate the speed required for an object to escape the gravity of a celestial body. The formula is: v = sqrt(2 × G × M / r), where v is the escape velocity, G is the gravitational constant, M is the mass of the celestial body and r is the radius of the celestial body.

Wien's Law: This law describes the relationship between the temperature of an object and the maximum intensity of its blackbody radiation. The formula is: λ_max = b / T, where λ_max is the peak radiation wavelength, b is the Wien constant and T is the temperature of the object in kelvin.

Stefan-Boltzmann Law: This law describes the total amount of radiant energy emitted by a black body as a function of its temperature. The formula is: L = σ × A × T^4, where L is the luminosity, σ is the Stefan-Boltzmann constant, A is the surface area of the object and T is the temperature in kelvin.

Schwarzschild equation: This is a solution to the spacetime metric around a massive spherical body, such as a star or black hole. The formula is used in Einstein's theory of general relativity.

Hubble's Law: This law describes the expansion of the universe, stating that the speed of recession of a galaxy is proportional to its distance. The formula is: v = H0 × d, where v is the recession velocity, H0 is the Hubble constant and d is the distance to the galaxy.

Roche's Formula: This formula is used to determine the limit at which a celestial body, such as a moon, can hold together due to gravitational force relative to another body, such as a planet. This is important for understanding the formation of planetary rings. The formula is derived from gravitational and centrifugal equilibrium equations.

Radiation Transfer Equations: These are a series of partial differential equations that describe how radiation propagates and interacts with matter in space, and are fundamental for modeling stars, planets and other astronomical objects.

Friedmann equation: This is one of the fundamental equations of cosmology that describes the evolution of the universe on a large scale. It is a solution of Einstein's field equations in general relativity.

     Law of Conservation of Angular Momentum: This law describes the conservation of angular momentum in astronomical systems such as stars, galaxies and planetary systems. Mathematically, it can be expressed as L = Iω, where L is the angular momentum, I is the moment of inertia and ω is the angular velocity.

     Newton's Law of Motion: This is one of the fundamental laws of physics that describes how a body responds to a force applied to it. Newton's second law is often expressed as F = ma, where F is the applied force, m is the mass of the body, and a is the resulting acceleration.

     Polythermal Equation of State: This equation is used in stellar astrophysics to describe the behavior of a polythermal gas (a gas that does not follow the ideal gas law) in terms of pressure, density and temperature.

     Titius-Bode Law: Although not a true physical law, this is an empirical rule that suggests a relationship between the average distances of planets around the Sun. Although it was initially considered a coincidence, it has historically been used to predict the existence of planets beyond Mars.

     Lane-Emden equation: This is a differential equation that describes the internal structure of a star in hydrostatic equilibrium, assuming that the star is spherically symmetric and composed of a polythermal gas.

     Chandrasekhar equation: This is an equation that defines the Chandrasekhar limit, which is the maximum mass a white dwarf can have before collapsing into a supernova or becoming a black hole, due to electron degeneracy pressure.

     Parallax Formula: This is a formula used to measure the distance of nearby astronomical objects, such as stars, using the parallax technique. The formula relates the distance to the observed parallax angle.

     Laplace-Poisson equation: This is a partial differential equation used in astrophysics to describe the gravitational field generated by mass distributions.

     Law of the Special Case of Relativity: This is one of the fundamental equations of Einstein's special theory of relativity, which relates the energy of a particle to its mass and velocity. The famous equation is E = mc^2, where E is the energy, m is the mass and c is the speed of light in a vacuum.
    
Inverse Square Law of Distance: This law describes how the intensity of light or radiation decreases as the distance increases. Mathematically, it can be expressed as I = P / (4πd^2), where I is the intensity, P is the radiated power and d is the distance.

     Law of Radioactive Half-Life: This law describes the rate at which a radioactive isotope decays over time. It is expressed by the formula N(t) = N0 × e^(-λt), where N(t) is the amount of isotope remaining at time t, N0 is the initial amount, λ is the decay constant and e is the base of the natural logarithm.

     Luminous Flux Formula: This formula is used to calculate the luminous flux of a light source in a given area. Mathematically, it is expressed as F = L / (4πd^2), where F is the luminous flux, L is the luminosity of the source and d is the distance to the source.

     Law of Conservation of Energy: This fundamental law of physics states that the total energy of an isolated system remains constant over time. It is expressed by the equation E = K + U, where E is the total energy, K is the kinetic energy and U is the potential energy.

     Law of Cosmological Redshift: This law describes how light from distant objects is redshifted due to the expansion of the universe. It is expressed by the formula z = (λ_obs - λ_em) / λ_em, where z is the redshift, λ_obs is the observed wavelength and λ_em is the emitted wavelength.

     Absolute Magnitude Formula: This formula is used to calculate the absolute magnitude of an astronomical object, which is a measure of its intrinsic luminosity. Mathematically, it is expressed as M = m - 5(log(d) - 1), where M is the absolute magnitude, m is the apparent magnitude and d is the distance in parsecs.

     Solar Constant Equation: This equation describes the rate at which solar energy is radiated at the Earth's surface. The solar constant is approximately 1361 watts per square meter.

     Critical Density Formula: This formula is used in cosmology to determine the critical density of the universe, which is the density necessary for the expansion of the universe to eventually stop. It is expressed as ρ_c = 3H^2 / (8πG), where ρ_c is the critical density, H is the Hubble constant and G is the gravitational constant.

     Age of the Universe Formula: This formula is used to calculate the approximate age of the universe based on cosmological observations. It is expressed as t = 1 / H0, where t is the age of the universe and H0 is the Hubble constant.

     Maxwell-Boltzmann Equilibrium Law: This law describes the distribution of particle velocities in an ideal gas. It is expressed by the equation f(v) = (m / 2πkT)^(3/2) * 4πv^2 * e^(-mv^2 / 2kT), where f(v) is the velocity distribution function, m is the mass of the particle, v is the velocity, k is the Boltzmann constant and T is the temperature.


Orbital Eccentricity Formula: This formula is used to describe the shape of the orbit of one celestial body around another. It is expressed as e = (r_max - r_min) / (r_max + r_min), where e is the orbital eccentricity, r_max is the radius of the orbit at the furthest point and r_min is the radius at the closest point.

     Tully-Fisher Relation Formula: This formula is used to estimate the luminosity of a spiral galaxy based on its spectral line width of 21 cm. It is expressed as L ∝ v^4, where L is the luminosity and v is the width of the spectral line.

     Radiation Pressure Formula: This formula is used to calculate the pressure exerted by electromagnetic radiation on a body. It is expressed as P = (4σ/c) × T^4, where P is the pressure, σ is the Stefan-Boltzmann constant, c is the speed of light in a vacuum and T is the temperature.

     Law of Conservation of Linear Momentum: This law states that the total linear momentum of an isolated system remains constant as long as no external force acts on it. It is expressed by the equation ∑p = 0, where ∑p is the sum of the linear moments of all particles in the system.

     Roche Period Formula: This formula is used to calculate the rotation period of an orbiting moon due to the gravitational pull of the main body. It is expressed as T = 2π × √(R^3 / (G × M)), where T is the period of rotation, R is the radius of the moon, G is the gravitational constant and M is the mass of the main body.

     Flux Density Formula: This formula is used to calculate the flux density of radiation or energy passing through a specific area. It is expressed as F = Φ / A, where F is the flux density, Φ is the total flow and A is the area.

     Cosmological Radial Distance Formula: This formula is used to calculate the radial distance of a celestial object based on the cosmological redshift. It is expressed as d = c × z / H0, where d is the radial distance, c is the speed of light, z is the redshift, and H0 is the Hubble constant.

     Orbital Speed Formula: This formula is used to calculate the orbital speed of an object in orbit around another body. It is expressed as v = √(G × M / r), where v is the orbital velocity, G is the gravitational constant, M is the mass of the central body, and r is the distance of the object from the center of the central body.

     Effective Temperature Formula: This formula is used to calculate the effective temperature of an astronomical object, such as a star, based on its luminosity and radius. It is expressed as Teff = (L / 4πσR^2)^1/4, where Teff is the effective temperature, L is the luminosity, σ is the Stefan-Boltzmann constant and R is the radius.

     Star Burst Size Ratio Formula: This formula is used to determine the relationship between a star's burst size and its intrinsic luminosity. It is expressed as L ∝ R^2 × T^4, where L is the luminosity, R is the radius of the star and T is the effective temperature.
    
    
    Radial Velocity Equation: This equation is used to calculate the radial velocity of a celestial object relative to the Earth along the line of sight. It is expressed as v_rad = c × (Δλ / λ), where v_rad is the radial velocity, c is the speed of light, Δλ is the change in wavelength of the observed light, and λ is the original wavelength.

     Blueshift Formula: This formula is used to calculate the blueshift of light from a celestial object moving towards Earth. It is expressed as z = Δλ / λ, where z is the blueshift, Δλ is the change in wavelength of the observed light, and λ is the original wavelength.

     Luminosity-Mass Relationship Formula: This formula is used to relate the luminosity of a star to its mass. It is expressed as L ∝ M^3.5, where L is the luminosity and M is the mass of the star.

     Period-Luminosity Relationship Formula: This formula is used to relate the period of variation in brightness of a variable star to its intrinsic luminosity. It is especially useful for Cepheid stars. It is expressed as L ∝ P^3.5, where L is the luminosity and P is the period.

     Generalized Lane-Emden Equation: This is a generalized form of the Lane-Emden equation, used to describe the internal structure of a star in hydrostatic equilibrium, considering different equations of state for different types of stars.

     Apparent Magnitude Formula: This formula is used to calculate the apparent magnitude of an astronomical object, which is a measure of its brightness as observed from Earth. It is expressed as m = -2.5 × log(F), where m is the apparent magnitude and F is the energy flux received on Earth.

     White Dwarf Formula: This formula is used to calculate the radius of a white dwarf, which is a final step in stellar evolution for stars of similar mass to the Sun. It is expressed as R = (ℏ^2 / (2πme^2) ) × (3 / (8πGρ))^(1/3), where R is the radius, ℏ is the reduced Planck constant, me is the mass of the electron, G is the gravitational constant, and ρ is the density of the white dwarf .

     Orbital Period Formula: This formula is used to calculate the orbital period of a celestial object in orbit around another body. It is expressed as T = 2π × √(a^3 / (G(M1 + M2))), where T is the orbital period, a is the semi-major axis of the orbit, G is the gravitational constant, and M1 and M2 are the masses of the orbiting bodies.

     Cosmological Comoable Distance Formula: This formula is used to calculate the comoving distance between celestial objects in the expanding universe. It is expressed as D = c × ∫(1/a) dt, where D is the comoving distance, c is the speed of light and a is the cosmological scale factor.

     Dark Energy Density Formula: This formula is used to calculate the density of dark energy, a hypothetical form of energy believed to be responsible for the accelerated expansion of the universe. It is expressed as ρ_Λ = Λc^2 / (8πG), where ρ_Λ is the dark energy density, Λ is the cosmological constant, c is the speed of light and G is the gravitational constant.
    
    Comovable Distance Formula in Cosmology: This formula is used to calculate the comoving distance between two objects in an expanding universe. It is expressed as D_c = c × ∫(dt/a), where D_c is the comoving distance, c is the speed of light, dt is the time element, and a is the cosmological scale factor.

     Brightness Temperature Formula: This formula is used to calculate the brightness temperature of an astronomical object based on the flux density received on Earth. It is expressed as T_b = (F / σ)^1/4, where T_b is the brightness temperature and σ is the Stefan-Boltzmann constant.

     Angular Distance Formula: This formula is used to calculate the angular distance between two celestial objects as seen from Earth. It is expressed as θ = d / D, where θ is the angular distance, d is the linear distance between objects and D is the distance to objects.

     Age of a Star Formula: This formula is used to estimate the age of a star based on its mass and composition. It is expressed as t = (1 / τ) × ln(M_i / M_f), where t is the age of the star, τ is the characteristic time constant of the evolutionary stage, M_i is the initial mass of the star and M_f is the final mass of the star.

     Mass-Luminosity Relationship Formula: This formula is used to relate the mass of a star to its intrinsic luminosity. It is expressed as L ∝ M^α, where L is the luminosity and M is the mass of the star. The value of α depends on the age and composition of the star.

     Schwarzschild Ratio Formula: This formula is used to calculate the Schwarzschild radius, which defines the event horizon of a black hole. It is expressed as r_s = 2GM / c^2, where r_s is the Schwarzschild radius, G is the gravitational constant, M is the mass of the black hole and c is the speed of light.

     Tidality Ratio Formula: This formula is used to calculate the tidality ratio between two celestial bodies in a gravitationally interacting orbit. It is expressed as η = (M1 / M2) × (R / r)^3, where η is the tidal ratio, M1 and M2 are the masses of the bodies, R is the radius of the larger body and r is the distance between the bodies .

     Eddington's Luminosity Formula: This formula is used to calculate the maximum luminosity a star can reach before the radiation pressure exceeds the gravitational force, leading to instability. It is expressed as L_E = (4πGMm_p / σ_Tc) ≈ 3.86 × 10^26 W, where L_E is the Eddington luminosity, G is the gravitational constant, m_p is the mass of the proton, σ_T is the Thompson cross section and c is the velocity from light.

     Stefan-Boltzmann Law Formula for Effective Temperature: This formula is used to calculate the effective temperature of an astronomical object based on its luminosity and radius. It is expressed as T_eff = (L / (4πσR^2))^1/4, where T_eff is the effective temperature, L is the luminosity, σ is the Stefan-Boltzmann constant and R is the radius of the object.

     Rayleigh-Jeans Distribution Formula: This is one of two classical spectral distribution formulas that describe the energy distribution of a blackbody at high frequencies. It is expressed as B(λ, T) = (8πkT / λ^4), where B is the spectral irradiance, k is the Boltzmann constant, T is the temperature and λ is the wavelength.
    
    Virial Theorem Formula: This formula is a relationship between the kinetic energy and the potential energy of a stellar or galactic system in equilibrium. For an isolated system, it is expressed as 2K + U = 0, where K is the kinetic energy and U is the potential energy.

     Habitable Zone Ratio Formula: This formula is used to calculate the width of the habitable zone around a star, where conditions suitable for the existence of liquid water can be found. It is expressed as d_HZ = √(L / L_☉), where d_HZ is the width of the habitable zone, L is the luminosity of the star and L_☉ is the solar luminosity.

     Roche Ratio Formula for Binary Stars: This formula is used to calculate the Roche radius of a component of a binary star. It is expressed as r_R = 0.49 × (m / (M + m))^(1/3) × d, where r_R is the Roche radius, m is the mass of the smaller component, M is the mass of the larger component, and d is the distance between the centers of the components.

     Solar Mass Formula: This formula is used to calculate the mass of a star in terms of the solar mass. It is expressed as M = (L / L_☉) × (R / R_☉)^2, where M is the mass of the star, L is the luminosity, R is the radius, and L_☉ and R_☉ are the luminosity and radius solar, respectively.

     Distance Modulus Formula: This formula is used to calculate the distance modulus, which is a measurement of the distance between the Earth and an astronomical object based on its apparent and absolute magnitude. It is expressed as m - M = 5 × log(d / 10), where m is the apparent magnitude, M is the absolute magnitude and d is the distance in parsecs.

     Scale Distance Formula: This formula is used to calculate the moving distance between celestial objects in an expanding universe. It is expressed as d = c / H_0, where d is the scaling distance, c is the speed of light and H_0 is the Hubble constant.

     Warp Speed Formula: This formula is used to calculate the warp speed, which is the relative speed between two celestial objects due to the expansion of the universe. It is expressed as v = H_0 × d, where v is the distortion velocity, H_0 is the Hubble constant and d is the comoving distance between the objects.

     Decay Parameter Formula: This formula is used to calculate the decay parameter of an unstable subatomic particle. It is expressed as λ = ln(2) / τ, where λ is the decay parameter and τ is the decay time constant.

     Energy Difference Formula: This formula is used to calculate the energy difference between two energy levels in an atom or molecule. It is expressed as ΔE = hf, where ΔE is the energy difference, h is Planck's constant and f is the frequency of the absorbed or emitted radiation.

     Bode-Titius Law Formula Revisited: This formula is a modified version of the Bode-Titius Law, which attempts to predict the approximate distances of the planets in the Solar System from the Sun. It is expressed as a = (n + 4) / 10 , where a is the approximate distance in astronomical units (AU) and n is the number of the planet in the sequence.
    
    Black Hole Escape Velocity Formula: This formula is used to calculate the escape velocity from a black hole. It is expressed as v_esc = c × √(2GM / rc^2), where v_esc is the escape velocity, c is the speed of light, G is the gravitational constant, M is the mass of the black hole and r is the radius of the horizon of events.

     Cosmological Decomprehension Time Formula: This formula is used to calculate the time from the Big Bang to the present, also known as the age of the universe. It is expressed as t = 1 / H_0, where t is the decompression time and H_0 is the Hubble constant.

     Gravitational Shift Formula: This formula is used to calculate the shift of a ray of light passing near a massive body, such as a star or a galaxy. It is expressed as δθ = 4GM / (c^2b), where δθ is the angular deviation, G is the gravitational constant, M is the mass of the massive body, c is the speed of light and b is the impact parameter.

     Goldreich-Julian Ratio Formula: This formula is used to calculate the charge current density in a rapidly rotating particle beam, such as in pulsars. It is expressed as J = (ρB/e), where J is the current density, ρ is the charge density, B is the magnetic field and e is the electron charge.

     Critical Density Formula of the Universe: This formula is used to calculate the density of matter and energy required for the universe to be flat. It is expressed as ρ_c = 3H_0^2 / (8πG), where ρ_c is the critical density, H_0 is the Hubble constant and G is the gravitational constant.

     Spatial Signature Formula: This formula is used to characterize the curvature of spacetime around a massive object. It is expressed as s = (ds^2 / c^2) = g_ab dx^a dx^b, where s is the spatial signature, ds is the interval element, c is the speed of light, g_ab is the metric tensor and dx^a and dx^b are displacement elements.

     Tidal Force Formula: This formula is used to calculate the tidal force exerted on an object due to the gravitational field gradient of another object. It is expressed as F_tidal = (d^2Φ / dr^2) × Δm, where F_tidal is the tidal force, Φ is the gravitational potential, r is the radial distance and Δm is the mass difference.

     Space Curvature Formula: This formula is used to calculate the intrinsic curvature of space in a specific region. It is expressed as R_abcd = g_ad,g_bc - g_ac,g_bd + g_ac,g_bd, where R_abcd is the curvature tensor, g_ab is the metric tensor and commas indicate partial derivatives.

     Light Energy Flow Formula: This formula is used to calculate the energy flow that passes through a given area in a given time interval. It is expressed as F = dE / (Adt), where F is the energy flow, dE is the energy change, A is the area and dt is the time interval.

     Average Life of a Particle Formula: This formula is used to calculate the average life of an unstable particle. It is expressed as τ = h / (Γ), where τ is the half-life, h is Planck's constant and Γ is the decay width.
    
    Probability Density Function Formula: This formula is used to describe the probability of finding a particle at a certain position in a quantum system. It is expressed as ψ(x) = |Ψ(x)|^2, where ψ(x) is the probability density function and Ψ(x) is the wave function.

     Energy Flow Formula in General Relativity: This formula is used to calculate the energy flow in a region of curved spacetime according to the theory of general relativity. It is expressed as dE/dt = -T^ab × u_b, where dE/dt is the rate of change of energy, T^ab is the energy-momentum tensor and u_b is the four-velocity vector.

     Friedmann-Lemaître-Robertson-Walker Equation Formula: This formula describes the evolution of the universe on a large scale according to the theory of general relativity. It is expressed as (da/dt)^2 = (8πG/3)ρ - k/a^2, where a is the scale factor, t is the cosmic time, G is the gravitational constant, ρ is the density of matter -energy and k is the spatial curvature.

     Nuclear Binding Energy Formula: This formula is used to calculate the energy released or absorbed during nuclear reactions. It is expressed as ΔE = (Δm)c^2, where ΔE is the energy released or absorbed, Δm is the change in mass, and c is the speed of light.

     Gravitational Potential Formula: This formula is used to calculate the gravitational potential at a specific point around a massive object. It is expressed as Φ = -GM/r, where Φ is the gravitational potential, G is the gravitational constant, M is the mass of the object and r is the distance from the object to the point.

     Radioactive Decay Rate Formula: This formula is used to calculate the rate at which a sample of radioactive material decays. It is expressed as N(t) = N0 × e^(-λt), where N(t) is the number of nuclei remaining at time t, N0 is the initial number of nuclei, λ is the decay constant, and e is the basis of the natural logarithm.

     Chandrasekhar Mass Formula: This formula is used to calculate the maximum mass of a stable white dwarf. It is expressed as M = (ℏ^3 / (Gm_e^2)) × (3/2π)^(2/3) × (1 / (μ_e^2))^(1/3), where M is the mass from Chandrasekhar, ℏ is the reduced Planck constant, G is the gravitational constant, m_e is the electron mass, and μ_e is the average atomic mass per electron.

     Universe Expansion Rate Formula: This formula is used to calculate the expansion rate of the universe with respect to cosmic time. It is expressed as H(t) = H0 × √(Ω_r(1 + z)^4 + Ω_m(1 + z)^3 + Ω_k(1 + z)^2 + Ω_Λ), where H(t) is the rate of expansion, H0 is the Hubble constant, Ω_r, Ω_m, Ω_k and Ω_Λ are the densities of radiation, matter, curvature and dark energy, respectively, and z is the redshift.

     Tolman-Oppenheimer-Volkoff Limit Formula: This formula is used to calculate the maximum mass limit of a stable neutron star. It is expressed as M_max ≈ 0.7 × (ℏc / G) × (1 / μ_n^2) × (1 / m_p^2), where M_max is the maximum Tolman-Oppenheimer-Volkoff mass, ℏ is the reduced Planck constant, c is the speed of light, G is the gravitational constant, μ_n is the average atomic mass per neutron and m_p is the mass of the proton.

     Schwarzschild Rotation Rate Formula: This formula is used to calculate the rotation rate required for an object in orbit around a black hole to remain in stable orbit. It is expressed as ω = √(GM / r^3), where ω is the rotation rate, G is the gravitational constant, M is the mass of the black hole, and r is the distance of the object from the center of the black hole.
    
    Schwarzschild Radius Formula: This formula calculates the Schwarzschild radius, which defines the event horizon of a black hole. It is expressed as rs=2GMc2rs​=c22GM​, where rsrs​ is the Schwarzschild radius, GG is the gravitational constant, MM is the mass of the black hole and cc is the speed of light in a vacuum.

Event Horizon Area Formula: This formula calculates the area of the spherical surface that constitutes the event horizon of a black hole. It is expressed as A=4πrs2A=4πrs2​, where AA is the area of the event horizon and rsrs​ is the Schwarzschild radius.

Hawking Temperature Formula: This formula calculates the temperature of Hawking radiation emitted by a black hole. It is expressed as T=ℏc38πkGMT=8πkGMℏc3​, where TT is the Hawking temperature, ℏℏ is the reduced Planck constant, cc is the speed of light, kk is the Boltzmann constant, GG is the gravitational constant and MM is the mass of the Black Hole.

Stable Circular Orbital Period Formula: This formula calculates the orbital period of a particle in stable circular orbit around a black hole. It is expressed as T=2πr3GMT=2πGMr3​

​, where TT is the orbital period, rr is the orbital radius and MM is the mass of the black hole.

Gravitational Field Energy Density Formula: This formula calculates the energy density of the gravitational field around a black hole. It is expressed as ρ=3c632πG3M2ρ=32πG3M23c6​, where ρρ is the energy density, cc is the speed of light, GG is the gravitational constant and MM is the mass of the black hole.

Exotic Matter Energy Density Formula: This formula calculates the energy density of exotic matter needed to keep a wormhole open. It is expressed as ρ=−c48πG(g00+g11)ρ=−8πGc4​(g00​+g11​), where ρρ is the energy density, cc is the speed of light, GG is the gravitational constant and g00g00​ and g11g11 ​ are components of the Morris-Thorne metric.

     Spacetime Stress Flow Formula: This formula calculates the stress flow (gravitational stress) required to keep a wormhole open. It is a complex expression based on the Morris-Thorne metric.

     Dynamic Stability Formula: This formula evaluates the dynamic stability of a wormhole, determining whether it can remain open or will collapse. It is based on detailed analyzes of Einstein's field equations and specific boundary conditions for wormholes.

     Geodesic Extent Formula: This formula describes the geodesics (trajectories of particle motion) around a wormhole, providing information about how particles would travel through the wormhole. It is derived from the equations of particle motion in curved spacetime.
    
    Minkowski Metric Formula: The Minkowski metric describes the flat spacetime of the universe in the absence of mass or energy. It is expressed as ds2=−c2dt2+dx2+dy2+dz2ds2=−c2dt2+dx2+dy2+dz2, where dsds is the space-time interval element, cc is the speed of light, dtdt is the time difference, and dxdx , dydy, dzdz are the spatial differences in the xx, yy and zz directions.

Schwarzschild Metric Formula: The Schwarzschild metric describes the spacetime around a spherical, non-rotating black hole. It is expressed as ds2=−(1−rsr)c2dt2+dr21−rsr+r2(dθ2+sin⁡2θdϕ2)ds2=−(1−rrs​​)c2dt2+1−rrs​​dr2​+r2(dθ2+sin2θdϕ2 ), where rsrs​ is the Schwarzschild radius, rr is the radial distance, and θθ, ϕϕ are the spherical angles.

Kerr Metric Formula: The Kerr metric describes the spacetime around a rotating black hole. It is a generalization of the Schwarzschild metric and includes the rotation parameter aa. The Kerr metric is quite complex and includes additional terms that describe the effects of rotation.

Friedmann-Lemaître-Robertson-Walker (FLRW) Metric Formula: The FLRW metric describes the spacetime of a large-scale homogeneous and isotropic universe, as predicted by the standard cosmological model. It is expressed as ds2=−c2dt2+a(t)2(dr21−kr2+r2(dθ2+sin⁡2θdϕ2))ds2=−c2dt2+a(t)2(1−kr2dr2​+r2(dθ2+sin2θdϕ2)) , where a(t)a(t) is the scale factor of the universe and kk is the spatial curvature.

Time Contraction Formula (Time Dilation): Time contraction describes the phenomenon in which time passes more slowly for a moving observer compared to a stationary observer. It is expressed as dt′=dt1−v2/c2dt′=1−v2/c2

​dt​, where dt′dt′ is the time interval measured by the moving observer, dtdt is the time interval measured by the stationary observer, vv is the relative speed and cc is the speed of light.

Length Dilation Formula: Length dilation describes the phenomenon in which the length of a moving object is measured to be greater than that when it is at rest. It is expressed as L′=L1−v2c2L′=L1−c2v2​

​, where L′L′ is the length measured by the moving observer, LL is the proper length of the object, vv is the relative speed and cc is the speed of light.

Spatial Curvature Formula: Spatial curvature describes the intrinsic curvature of space in a given region. In an FLRW metric, spatial curvature is described by the parameter kk. It can be positive (closed space), negative (open space) or zero (flat space).

Metric Tensor Formula: The metric tensor describes the geometry of spacetime in a given region. In Cartesian coordinates, the metric tensor for the Minkowski metric is a diagonal matrix with elements −1,1,1,1−1,1,1,1 along the diagonal.

Lorentz Contraction Formula: Lorentz contraction describes the reduction in length of a moving object in the direction of motion, according to the theory of special relativity. It is expressed as L′=LγL′=γL​, where L′L′ is the length measured in the moving frame, LL is the proper length of the object and γγ is the Lorentz factor, given by γ=11−v2c2γ=1 −c2v2​

     ​1​.

     Gravitational Contraction Formula: Gravitational contraction is an effect predicted by the theory of general relativity in which space-time is curved in the presence of mass and energy. The Schwarzschild metric is an example of how the presence of mass causes spacetime to warp around a massive body.
	
	`)

func DrakeEquation(N, fp, ne, fl, fi, fc, L float64) float64 {
	return N * fp * ne * fl * fi * fc * L
}

func KeplersThirdLaw(a, T float64) float64 {
	const G = 6.67430e-11
	return math.Pow((math.Pow(a, 3)/(4*math.Pow(math.Pi, 2))), 1/3) * math.Sqrt((4*math.Pow(math.Pi, 2))/(G*T*T))
}

func UniversalGravitationForce(m1, m2, r float64) float64 {
	const G = 6.67430e-11
	return (G * m1 * m2) / (r * r)
}

func EscapeVelocity(m, r float64) float64 {
	const G = 6.67430e-11
	return math.Sqrt((2 * G * m) / r)
}

func WienDisplacementLaw(T float64) float64 {
	const b = 2.8977729e-3
	return b / T
}

func StefanBoltzmannLaw(T float64) float64 {
	const sigma = 5.670374419e-8
	return sigma * math.Pow(T, 4)
}

func SchwarzschildRadius(m float64) float64 {
	const G = 6.67430e-11
	const c = 299792458
	return (2 * G * m) / (c * c)
}

func HubbleLaw(v, d float64) float64 {
	const H0 = 70
	return v / d
}

func RocheLimit(density1, density2, radius1, radius2 float64) float64 {
	return 2.44 * (radius1 + radius2) * math.Pow(density1/density2, 1/3)
}

func FriedmannEquation(H, rho, k float64) float64 {
	const G = 6.67430e-11
	return math.Pow(H, 2) - (8*math.Pi*G*rho)/3 + k
}

func ConservationOfAngularMomentum(m, r, v float64) float64 {
	return m * r * v
}

func NewtonsSecondLaw(F, m float64) float64 {
	return F / m
}

func PolytropicEquationOfState(P, K, rho, n float64) float64 {
	return P - K*math.Pow(rho, n)
}

func TitiusBodeLaw(n float64) float64 {
	if n == 0 {
		return 0
	}
	return 0.4 + (0.3 * n)
}

func ChandrasekharEquation(m float64) float64 {
	const k = 1.8751
	const mp = 1.6726219e-27
	const h = 6.62607015e-34
	const c = 299792458
	const G = 6.67430e-11

	return math.Sqrt((math.Pow(k, 3) * math.Pow(mp, 5)) / (math.Pow(h, 3) * math.Pow(c, 3) * G * m))
}

func ParallaxFormula(d, θ float64) float64 {
	return d / math.Tan(θ)
}

func SpecialRelativityEquation(E, m, c float64) float64 {
	return E / (m * c * c)
}

func InverseSquareLaw(P, d float64) float64 {
	return P / (4 * math.Pi * math.Pow(d, 2))
}

func RadioactiveHalfLife(N0, lambda, t float64) float64 {
	return N0 * math.Exp(-lambda*t)
}

func LuminousFlux(L, d float64) float64 {
	return L / (4 * math.Pi * math.Pow(d, 2))
}

func ConservationOfEnergy(K, U float64) float64 {
	return K + U
}

func cosmologicalRedshift(lambda_obs, lambda_em float64) float64 {
	return (lambda_obs - lambda_em) / lambda_em
}

func AbsoluteMagnitude(m, d float64) float64 {
	return m - 5*(math.Log10(d)-1)
}

const solarConstant = 1361

func CriticalDensity(H float64) float64 {
	const G = 6.67430e-11
	return 3 * math.Pow(H, 2) / (8 * math.Pi * G)
}

func AgeOfUniverse(H0 float64) float64 {
	return 1 / H0
}

func MaxwellBoltzmannEquilibrium(m, v, T float64) float64 {
	const k = 1.380649e-23
	return math.Pow(m/(2*math.Pi*k*T), 3/2) * 4 * math.Pi * math.Pow(v, 2) * math.Exp(-m*v*v/(2*k*T))
}

func OrbitalEccentricity(a, b float64) float64 {
	return math.Sqrt(1 - math.Pow(b, 2)/math.Pow(a, 2))
}

func TullyFisherRelation(v, a, b float64) float64 {
	return math.Log10(v) - a*math.Log10(b)
}

func RadiationPressure(I, c float64) float64 {
	return I / c
}

func ConservationOfLinearMomentum(m1, v1, m2, v2 float64) float64 {
	return m1*v1 + m2*v2
}

func FluxDensity(F, A float64) float64 {
	return F / A
}

func OrbitalVelocity(G, M, r float64) float64 {
	return math.Sqrt((G * M) / r)
}

func StellarBurstSizeRatio(L1, L2 float64) float64 {
	return math.Sqrt(L1 / L2)
}

func RadialVelocity(c, deltaLambda, lambda float64) float64 {
	return c * (deltaLambda / lambda)
}

func BlueShift(deltaLambda, lambda float64) float64 {
	return deltaLambda / lambda
}

func LuminosityMassRelation(M float64) float64 {
	return math.Pow(M, 3.5)
}

func PeriodLuminosityRelation(P float64) float64 {
	return math.Pow(P, 3.5)
}

func ApparentMagnitude(F float64) float64 {
	return -2.5 * math.Log10(F)
}

func WhiteDwarfRadius(hbar, me, G, rho float64) float64 {
	return math.Pow((hbar*hbar/(2*math.Pi*me*me))*(3/(8*math.Pi*G*rho)), 1/3)
}

func OrbitalPeriod(a, G, M1, M2 float64) float64 {
	return 2 * math.Pi * math.Sqrt(math.Pow(a, 3)/(G*(M1+M2)))
}

func ComovingCosmologicalDistance(c, a, dt float64) float64 {
	return c * integral(1/a, dt)
}

func DarkEnergyDensity(Lambda, c, G float64) float64 {
	return (Lambda * c * c) / (8 * math.Pi * G)
}

func StarAge(tau, Mi, Mf float64) float64 {
	return (1 / tau) * math.Log(Mi/Mf)
}

func TidalRatio(M1, M2, R, r float64) float64 {
	return (M1 / M2) * math.Pow(R/r, 3)
}

func EddingtonLuminosity(G, M, mp, sigma_T, c float64) float64 {
	return (4 * math.Pi * G * M * mp) / (sigma_T * c)
}

func StefanBoltzmannEffectiveTemperature(L, R float64) float64 {
	const sigma = 5.670374419e-8
	return math.Pow((L / (4 * math.Pi * sigma * math.Pow(R, 2))), 1/4)
}

func RayleighJeansDistribution(k, T, lambda float64) float64 {
	return (8 * math.Pi * k * T) / math.Pow(lambda, 4)
}

func VirialTheorem(K, U float64) float64 {
	return 2*K + U
}

func HabitableZoneRatio(L, Lsun float64) float64 {
	return math.Sqrt(L / Lsun)
}

func RocheLimitRatio(m, M, d float64) float64 {
	return 0.49 * math.Pow(m/(M+m), 1/3) * d
}

func SolarMass(L, Lsun, R, Rsun float64) float64 {
	return (L / Lsun) * math.Pow(R/Rsun, 2)
}

func DistanceModulus(m, M, d float64) float64 {
	return m - M - 5*math.Log10(d/10)
}

func ScaleDistance(c, H0 float64) float64 {
	return c / H0
}

func DistortionVelocity(H0, d float64) float64 {
	return H0 * d
}

func DecayParameter(tau float64) float64 {
	return math.Log(2) / tau
}

func EnergyDifference(h, f float64) float64 {
	return h * f
}

func BodeTitiusLaw(n float64) float64 {
	return (n + 4) / 10
}

func BlackHoleEscapeVelocity(c, G, M, r float64) float64 {
	return c * math.Sqrt(2*G*M/(r*c*c))
}

func CosmologicalDecompressionTime(H0 float64) float64 {
	return 1 / H0
}

func GravitationalLensing(G, M, c, b float64) float64 {
	return 4 * G * M / (c * c * b)
}

func GoldreichJulianRatio(rho, B, e float64) float64 {
	return rho * B / e
}

func CriticalDensityUniverse(H0, G float64) float64 {
	return 3 * math.Pow(H0, 2) / (8 * math.Pi * G)
}

func SpatialSignature(ds, c float64) float64 {
	return ds / (c * c)
}

func TidalForce(d2Phi_dr2, delta_m float64) float64 {
	return d2Phi_dr2 * delta_m
}

func SpatialCurvature(g_ad_g_bc, g_ac_g_bd, g_ac_g_bd float64) float64 {
	return g_ad_g_bc - g_ac_g_bd + g_ac_g_bd
}

func LuminousEnergyFlux(dE, A, dt float64) float64 {
	return dE / (A * dt)
}

func ParticleLifetime(h, Gamma float64) float64 {
	return h / Gamma
}

func ProbabilityDensityFunction(psi float64) float64 {
	return math.Pow(math.Abs(psi), 2)
}

func SchwarzschildRadius(G, M, c float64) float64 {
	return (2 * G * M) / (c * c)
}

func EventHorizonArea(rs float64) float64 {
	return 4 * math.Pi * math.Pow(rs, 2)
}

func HawkingTemperature(G, M, h, c, k float64) float64 {
	return (h * c * c * c) / (8 * math.Pi * G * M * k)
}

func StableOrbitalPeriod(G, M, r float64) float64 {
	return 2 * math.Pi * math.Sqrt(math.Pow(r, 3)/(G*M))
}

func GravitationalFieldEnergyDensity(G, M, c float64) float64 {
	return (3 * c * c) / (32 * math.Pi * G * math.Pow(M, 2))
}

func ExoticMatterEnergyDensity(c, G, g00, g11 float64) float64 {
	return -(c * c) / (8 * math.Pi * G * (g00 + g11))
}

func MinkowskiMetric(ds, dt, dx, dy, dz, c float64) float64 {
	return -c*c*dt*dt + dx*dx + dy*dy + dz*dz
}

func TimeDilation(dt, v, c float64) float64 {
	return dt * math.Sqrt(1-(v*v)/(c*c))
}

func LengthContraction(L, v, c float64) float64 {
	return L * math.Sqrt(1-(v*v)/(c*c))
}

func LorentzContraction(L, gamma float64) float64 {
	return L / gamma
}

func AmpereMaxwellLaw(B, J, E, mu0, epsilon0, dE_dt float64) float64 {
	return mu0*J + mu0*epsilon0*dE_dt
}

func BernoulliEquation(P1, rho, v1, h1, g, P2, v2, h2 float64) float64 {
	return P1 + (0.5 * rho * v1 * v1) + (rho * g * h1) - P2 - (0.5 * rho * v2 * v2) - (rho * g * h2)
}

func EnergyMassEquivalence(m, c float64) float64 {
	return m * c * c
}

func CoulombLaw(k, q1, q2, r float64) float64 {
	return k * math.Abs(q1*q2) / (r * r)
}

func SnellDescartesLaw(n1, theta1, n2 float64) float64 {
	return (n1 / n2) * math.Sin(theta1)
}

func KineticGasEquation(n, R, T float64) float64 {
	return n * R * T
}

func UncertaintyHeisenberg(h float64) float64 {
	return h / 2
}

func PlanckEinsteinEquation(h, f float64) float64 {
	return h * f
}

func EinsteinEquation(p, c, m0 float64) float64 {
	return math.Sqrt(math.Pow((p*c), 2) + math.Pow((m0*math.Pow(c, 2)), 2))
}

func FaradayLaw(dPhi_dt float64) float64 {
	return -dPhi_dt
}

func OhmsLaw(I, R float64) float64 {
	return I * R
}

func HookeLaw(k, x float64) float64 {
	return -k * x
}

func AmpereLaw(B_integral, μ0, I, ε0, dΦE_dt float64) float64 {
	return B_integral - μ0*I - μ0*ε0*dΦE_dt
}

func MaxwellDivergence(E_divergence, p, ε0 float64) float64 {
	return E_divergence - p/ε0
}


}
