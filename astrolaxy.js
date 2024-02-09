function drakeEquation(N, fp, ne, fl, fi, fc, L) {
    return N * fp * ne * fl * fi * fc * L;
}

function keplersThirdLaw(a, T) {
    const G = 6.67430e-11; 
    return Math.pow((Math.pow(a, 3) / (4 * Math.pow(Math.PI, 2))), 1 / 3) * Math.sqrt((4 * Math.pow(Math.PI, 2)) / (G * T * T));
}

function universalGravitationForce(m1, m2, r) {
    const G = 6.67430e-11; 
    return (G * m1 * m2) / (r * r);
}

function escapeVelocity(m, r) {
    const G = 6.67430e-11; 
    return Math.sqrt((2 * G * m) / r);
}

function wienDisplacementLaw(T) {
    const b = 2.8977729e-3; 
    return b / T;
}

function stefanBoltzmannLaw(T) {
    const sigma = 5.670374419e-8; 
    return sigma * Math.pow(T, 4);
}

function schwarzschildRadius(m) {
    const G = 6.67430e-11; 
    const c = 299792458; 
    return (2 * G * m) / (c * c);
}

function hubbleLaw(v, d) {
    const H0 = 70; 
    return v / d;
}

function rocheLimit(density1, density2, radius1, radius2) {
    return 2.44 * (radius1 + radius2) * Math.pow(density1 / density2, 1 / 3);
}

function friedmannEquation(H, rho, k) {
    return Math.pow(H, 2) === (8 * Math.PI * G * rho) / 3 - k;
}


function conservationOfAngularMomentum(m, r, v) {
    return m * r * v;
}

function newtonsSecondLaw(F, m) {
    return F / m;
}

function polytropicEquationOfState(P, K, rho, n) {
    return P = K * Math.pow(rho, n);
}

function titiusBodeLaw(n) {
    if (n === 0) return 0;
    return 0.4 + (0.3 * n);
}

function chandrasekharEquation(m) {
    const k = 1.8751; 
    const mp = 1.6726219e-27; 
    const h = 6.62607015e-34; 
    const c = 299792458; 
    const G = 6.67430e-11; 

    return Math.sqrt((Math.pow(k, 3) * Math.pow(mp, 5)) / (Math.pow(h, 3) * Math.pow(c, 3) * G * m));
}

function parallaxFormula(d, θ) {
    return d / Math.tan(θ);
}

function specialRelativityEquation(E, m, c) {
    return E / (m * c * c);
}

function inverseSquareLaw(P, d) {
    return P / (4 * Math.PI * Math.pow(d, 2));
}

function radioactiveHalfLife(N0, lambda, t) {
    return N0 * Math.exp(-lambda * t);
}

function luminousFlux(L, d) {
    return L / (4 * Math.PI * Math.pow(d, 2));
}

function conservationOfEnergy(K, U) {
    return K + U;
}

function cosmologicalRedshift(lambda_obs, lambda_em) {
    return (lambda_obs - lambda_em) / lambda_em;
}

function absoluteMagnitude(m, d) {
    return m - 5 * (Math.log10(d) - 1);
}

const solarConstant = 1361; 

function criticalDensity(H) {
    const G = 6.67430e-11; 
    return 3 * Math.pow(H, 2) / (8 * Math.PI * G);
}

function ageOfUniverse(H0) {
    return 1 / H0;
}

function maxwellBoltzmannEquilibrium(m, v, T) {
    const k = 1.380649e-23; 
    return Math.pow(m / (2 * Math.PI * k * T), 3 / 2) * 4 * Math.PI * Math.pow(v, 2) * Math.exp(-m * v * v / (2 * k * T));
}

function orbitalEccentricity(a, b) {
    return Math.sqrt(1 - Math.pow(b, 2) / Math.pow(a, 2));
}

function tullyFisherRelation(v, a, b) {
    return Math.log10(v) - a * Math.log10(b);
}

function radiationPressure(I, c) {
    return I / (c);
}

function conservationOfLinearMomentum(m1, v1, m2, v2) {
    return m1 * v1 + m2 * v2;
}

function fluxDensity(F, A) {
    return F / A;
}

function orbitalVelocity(G, M, r) {
    return Math.sqrt((G * M) / r);
}

function stellarBurstSizeRatio(L1, L2) {
    return Math.sqrt(L1 / L2);
}

function radialVelocity(c, deltaLambda, lambda) {
    return c * (deltaLambda / lambda);
}

function blueShift(deltaLambda, lambda) {
    return deltaLambda / lambda;
}

function luminosityMassRelation(M) {
    return Math.pow(M, 3.5);
}

function periodLuminosityRelation(P) {
    return Math.pow(P, 3.5);
}

function generalizedLaneEmdenEquation() {

}

function apparentMagnitude(F) {
    return -2.5 * Math.log10(F);
}

function whiteDwarfRadius(hbar, me, G, rho) {
    return Math.pow((hbar * hbar / (2 * Math.PI * me * me)) * (3 / (8 * Math.PI * G * rho)), 1 / 3);
}

function orbitalPeriod(a, G, M1, M2) {
    return 2 * Math.PI * Math.sqrt(Math.pow(a, 3) / (G * (M1 + M2)));
}

function comovingCosmologicalDistance(c, a, dt) {
    return c * integral(1 / a, dt);
}

function darkEnergyDensity(Lambda, c, G) {
    return (Lambda * c * c) / (8 * Math.PI * G);
}

function starAge(tau, Mi, Mf) {
    return (1 / tau) * Math.log(Mi / Mf);
}

function tidalRatio(M1, M2, R, r) {
    return (M1 / M2) * Math.pow(R / r, 3);
}

function eddingtonLuminosity(G, M, mp, sigma_T, c) {
    return (4 * Math.PI * G * M * mp) / (sigma_T * c);
}

function stefanBoltzmannEffectiveTemperature(L, R) {
    const sigma = 5.670374419e-8; 
    return Math.pow((L / (4 * Math.PI * sigma * Math.pow(R, 2))), 1/4);
}

function rayleighJeansDistribution(k, T, lambda) {
    return (8 * Math.PI * k * T) / Math.pow(lambda, 4);
}

function virialTheorem(K, U) {
    return 2 * K + U;
}

function habitableZoneRatio(L, Lsun) {
    return Math.sqrt(L / Lsun);
}

function rocheLimitRatio(m, M, d) {
    return 0.49 * Math.pow(m / (M + m), 1 / 3) * d;
}

function solarMass(L, Lsun, R, Rsun) {
    return (L / Lsun) * Math.pow(R / Rsun, 2);
}

function distanceModulus(m, M, d) {
    return m - M - 5 * Math.log10(d / 10);
}

function scaleDistance(c, H0) {
    return c / H0;
}

function distortionVelocity(H0, d) {
    return H0 * d;
}

function decayParameter(tau) {
    return Math.log(2) / tau;
}

function energyDifference(h, f) {
    return h * f;
}

function bodeTitiusLaw(n) {
    return (n + 4) / 10;
}

function blackHoleEscapeVelocity(c, G, M, r) {
    return c * Math.sqrt(2 * G * M / (r * c * c));
}

function cosmologicalDecompressionTime(H0) {
    return 1 / H0;
}

function gravitationalLensing(G, M, c, b) {
    return 4 * G * M / (c * c * b);
}

function goldreichJulianRatio(rho, B, e) {
    return rho * B / e;
}

function criticalDensityUniverse(H0, G) {
    return 3 * Math.pow(H0, 2) / (8 * Math.PI * G);
}

function spatialSignature(ds, c) {
    return ds / (c * c);
}

function tidalForce(d2Phi_dr2, delta_m) {
    return d2Phi_dr2 * delta_m;
}

function spatialCurvature(g_ad_g_bc, g_ac_g_bd, g_ac_g_bd) {
    return g_ad_g_bc - g_ac_g_bd + g_ac_g_bd;
}

function luminousEnergyFlux(dE, A, dt) {
    return dE / (A * dt);
}

function particleLifetime(h, Gamma) {
    return h / Gamma;
}

function probabilityDensityFunction(psi) {
    return Math.pow(Math.abs(psi), 2);
}

function schwarzschildRadius(G, M, c) {
    return (2 * G * M) / (c * c);
}

function eventHorizonArea(rs) {
    return 4 * Math.PI * Math.pow(rs, 2);
}

function hawkingTemperature(G, M, h, c, k) {
    return (h * c * c * c) / (8 * Math.PI * G * M * k);
}

function stableOrbitalPeriod(G, M, r) {
    return 2 * Math.PI * Math.sqrt(Math.pow(r, 3) / (G * M));
}

function gravitationalFieldEnergyDensity(G, M, c) {
    return (3 * c * c) / (32 * Math.PI * G * Math.pow(M, 2));
}

function exoticMatterEnergyDensity(c, G, g00, g11) {
    return -(c * c) / (8 * Math.PI * G * (g00 + g11));
}

function minkowskiMetric(ds, dt, dx, dy, dz, c) {
    return -c * c * dt * dt + dx * dx + dy * dy + dz * dz;
}

function timeDilation(dt, v, c) {
    return dt * Math.sqrt(1 - (v * v) / (c * c));
}

function lengthContraction(L, v, c) {
    return L * Math.sqrt(1 - (v * v) / (c * c));
}

function lorentzContraction(L, gamma) {
    return L / gamma;
}

function ampereMaxwellLaw(B, J, E, mu0, epsilon0, dE_dt) {
    return mu0 * J + mu0 * epsilon0 * dE_dt;
}

function bernoulliEquation(P1, rho, v1, h1, g, P2, v2, h2) {
    return P1 + (0.5 * rho * v1 * v1) + (rho * g * h1) - P2 - (0.5 * rho * v2 * v2) - (rho * g * h2);
}

function energyMassEquivalence(m, c) {
    return m * c * c;
}

function coulombLaw(k, q1, q2, r) {
    return k * Math.abs(q1 * q2) / (r * r);
}

function snellDescartesLaw(n1, theta1, n2) {
    return (n1 / n2) * Math.sin(theta1);
}

function kineticGasEquation(n, R, T) {
    return n * R * T;
}

function uncertaintyHeisenberg(h) {
    return h / 2;
}

function planckEinsteinEquation(h, f) {
    return h * f;
}

function einsteinEquation(p, c, m0) {
    return Math.sqrt(Math.pow((p * c), 2) + Math.pow((m0 * Math.pow(c, 2)), 2));
}


function faradayLaw(dPhi_dt) {
    return -dPhi_dt;
}

function ohmsLaw(I, R) {
    return I * R;
}

function hookeLaw(k, x) {
    return -k * x;
}

function ampereLaw(B_integral, μ0, I, ε0, dΦE_dt) {
    return B_integral - μ0 * I - μ0 * ε0 * dΦE_dt;
}

function maxwellDivergence(E_divergence, p, ε0) {
    return E_divergence - p / ε0;
}

