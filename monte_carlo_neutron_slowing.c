#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

// Physical constants
#define K_BOLTZMANN 1.380649e-23
#define M_NEUTRON 1.674927e-27
#define PI 3.14159265358979323846

// Simulation parameters
#define NUM_HISTORIES 1000000
#define NUM_INTERVALS 512
#define MAX_SPEED_RATIO 32.0
#define SOURCE_SPEED_RATIO 10.0
#define MAX_COLLISIONS 10000

typedef struct {
    double A;          // Moderator mass ratio (M/m_neutron)
    double K;          // Absorption parameter
    double T;          // Temperature (K)
    double v_T;        // Thermal speed
    
    // Cross sections
    double sigma_s;    // Scattering cross-section
    double sigma_a0;   // Absorption cross-section normalization
    double alpha_s;    // Alpha parameter
    
    // Tallies
    int scattering_tally[NUM_INTERVALS];
    int absorption_tally[NUM_INTERVALS];
    double speed_bins[NUM_INTERVALS + 1];
    
    // Statistics
    int total_histories;
    int total_scatterings;
    int total_absorptions;
    
    // RNG state
    uint64_t rng_state;
} Simulation;

// Function declarations
void init_simulation(Simulation* sim, double A, double K, double T);
void run_simulation(Simulation* sim);
void run_history(Simulation* sim);
double calculate_flux(Simulation* sim, double* y, double* flux, int* count);
void print_results(Simulation* sim);
void save_results(Simulation* sim, const char* filename);
double random_uniform(Simulation* sim);
void elastic_collision(Simulation* sim, double* v, double V, double mu);
double erf_approximation(double x);
double sample_maxwellian_speed(Simulation* sim, double temperature);
double maxwell_flux(double y);
double sample_cosine(Simulation* sim);
double calculate_effective_sigma_s(Simulation* sim, double v);

// Random number generator (xorshift64*)
double random_uniform(Simulation* sim) {
    sim->rng_state ^= sim->rng_state >> 12;
    sim->rng_state ^= sim->rng_state << 25;
    sim->rng_state ^= sim->rng_state >> 27;
    return (double)(sim->rng_state * 0x2545F4914F6CDD1DULL) / (double)UINT64_MAX;
}

// Initialize simulation parameters
void init_simulation(Simulation* sim, double A, double K, double T) {
    sim->A = A;
    sim->K = K;
    sim->T = T;
    
    // Calculate thermal speed (most probable speed for Maxwell-Boltzmann)
    sim->v_T = sqrt(2.0 * K_BOLTZMANN * T / M_NEUTRON);
    
    // Cross-section parameters
    sim->sigma_s = 1.0;
    sim->sigma_a0 = K * sim->sigma_s * sim->v_T;
    sim->alpha_s = sim->A * M_NEUTRON / (2.0 * K_BOLTZMANN * T);
    
    // Initialize speed bins
    double max_speed = MAX_SPEED_RATIO * sim->v_T;
    for (int i = 0; i <= NUM_INTERVALS; i++) {
        sim->speed_bins[i] = i * max_speed / NUM_INTERVALS;
    }
    
    // Initialize tallies
    for (int i = 0; i < NUM_INTERVALS; i++) {
        sim->scattering_tally[i] = 0;
        sim->absorption_tally[i] = 0;
    }
    
    sim->total_histories = 0;
    sim->total_scatterings = 0;
    sim->total_absorptions = 0;
    
    // Initialize RNG
    sim->rng_state = (uint64_t)time(NULL) + (uint64_t)(uintptr_t)sim;
}

// Error function approximation
double erf_approximation(double x) {
    // Abramowitz and Stegun approximation
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    int sign = (x < 0) ? -1 : 1;
    x = fabs(x);
    
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    
    return sign * y;
}

// Calculate effective scattering cross-section (paper's equation 4)
double calculate_effective_sigma_s(Simulation* sim, double v) {
    double alpha_moderator = (sim->A * M_NEUTRON) / (2.0 * K_BOLTZMANN * sim->T);
    
    double x = sqrt(alpha_moderator) * v;
    double term1 = exp(-alpha_moderator * v * v) / (sqrt(PI * alpha_moderator) * v);
    double term2 = (1.0 + 1.0/(2.0 * alpha_moderator * v * v)) * erf_approximation(x);
    
    return sim->sigma_s * (term1 + term2);
}

// Sample from Maxwell-Boltzmann distribution using Box-Muller
double sample_maxwellian_speed(Simulation* sim, double temperature) {
    double alpha = M_NEUTRON / (2.0 * K_BOLTZMANN * temperature);
    
    // Generate three independent normal random variables using Box-Muller
    double u1, u2, u3, u4, u5, u6;
    
    u1 = random_uniform(sim);
    u2 = random_uniform(sim);
    u3 = random_uniform(sim);
    u4 = random_uniform(sim);
    u5 = random_uniform(sim);
    u6 = random_uniform(sim);
    
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
    double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * PI * u2);
    double z2 = sqrt(-2.0 * log(u3)) * cos(2.0 * PI * u4);
    double z3 = sqrt(-2.0 * log(u3)) * sin(2.0 * PI * u4);
    
    // Magnitude of 3D velocity vector (speed)
    double speed = sqrt(z0*z0 + z1*z1 + z2*z2) / sqrt(2.0 * alpha);
    
    return speed;
}

// Sample cosine uniformly for isotropic scattering
double sample_cosine(Simulation* sim) {
    return 2.0 * random_uniform(sim) - 1.0;
}

// TARGET VELOCITY SAMPLING
double sample_target_velocity(Simulation* sim, double v) {
    double A = sim->A;
    double alpha_moderator = (A * M_NEUTRON) / (2.0 * K_BOLTZMANN * sim->T);
    
    // Calculate D(v) from paper Eq (6a-c)
    double sigma_s_eff = calculate_effective_sigma_s(sim, v);
    double sigma_a = sim->sigma_a0 / v;
    
    double D = v * (sigma_a + sigma_s_eff);
    
    // P1, P2, P3 from paper Eq (6a-c)
    double P1 = v * sigma_a / D;           // Absorption
    double P2 = sim->sigma_s * v / D;      // Scattering type 1
    double P3 = (2.0 * sim->sigma_s / sqrt(PI * alpha_moderator)) / D;  // Scattering type 2
    
    double V, mu, v_rel;
    int accepted = 0;
    
    while (!accepted) {
        // Choose which distribution to sample from
        double eta = random_uniform(sim);
        
        if (eta < P2 / (P2 + P3)) {
            // Sample from C2 * V² * exp(-αV²)
            // C2 = 4α^(3/2)/√π
            double C2 = 4.0 * pow(alpha_moderator, 1.5) / sqrt(PI);
            
            // Acceptance-rejection for V² exp(-αV²)
            double V_max = 5.0 / sqrt(alpha_moderator);  // 5σ cutoff
            do {
                V = V_max * random_uniform(sim);
                double f = C2 * V * V * exp(-alpha_moderator * V * V);
                double f_max = C2 * pow(2.0/alpha_moderator, 1.0) * exp(-2.0);
                if (random_uniform(sim) * f_max < f) break;
            } while (1);
        } else {
            // Sample from C3 * V³ * exp(-αV²)
            // C3 = 2α²
            double C3 = 2.0 * alpha_moderator * alpha_moderator;
            
            // Acceptance-rejection for V³ exp(-αV²)
            double V_max = 6.0 / sqrt(alpha_moderator);
            do {
                V = V_max * random_uniform(sim);
                double f = C3 * V * V * V * exp(-alpha_moderator * V * V);
                double f_max = C3 * pow(3.0/alpha_moderator, 1.5) * exp(-3.0);
                if (random_uniform(sim) * f_max < f) break;
            } while (1);
        }
        
        // Sample μ uniformly [-1, 1]
        mu = 2.0 * random_uniform(sim) - 1.0;
        
        // Calculate relative speed
        v_rel = sqrt(v*v + V*V - 2.0*v*V*mu);
        
        // PAPER'S REJECTION: accept with probability |v-V|/(v+V)
        if (random_uniform(sim) < v_rel / (v + V)) {
            accepted = 1;
        }
    }
    
    return V;
}

//collision physics 3D treatment
void elastic_collision(Simulation* sim, double* v, double V, double mu) {
    double A = sim->A;
    
    // Neutron velocity: assume along z-axis for convenience
    // v = (0, 0, v) in spherical coordinates
    double v_mag = *v;
    
    // Target velocity in 3D: V with direction (θ, φ)
    // μ = cosθ between neutron and target
    double theta = acos(mu);  // θ = arccos(μ)
    double phi = 2.0 * PI * random_uniform(sim);  // Random φ
    
    // Convert target velocity to Cartesian
    double Vx = V * sin(theta) * cos(phi);
    double Vy = V * sin(theta) * sin(phi);
    double Vz = V * mu;  // V * cosθ
    
    // Neutron velocity (along z-axis)
    double vz = v_mag;
    
    // Center of mass velocity (PAPER: g = (v + A*V)/(A+1))
    double gx = (0.0 + A * Vx) / (A + 1.0);
    double gy = (0.0 + A * Vy) / (A + 1.0);
    double gz = (vz + A * Vz) / (A + 1.0);
    
    // Relative velocity: v_rel = |v - V|
    double rel_x = 0.0 - Vx;
    double rel_y = 0.0 - Vy;
    double rel_z = vz - Vz;
    double v_rel = sqrt(rel_x*rel_x + rel_y*rel_y + rel_z*rel_z);
    
    // Neutron speed in CM frame: β|v_rel|, where β = A/(A+1)
    double beta = A / (A + 1.0);
    double v_cm_mag = v_rel * beta;
    
    // ISOTROPIC SCATTERING IN CM FRAME (3D sphere)
    // Uniform on unit sphere: cosθ_cm ~ Uniform[-1,1], φ_cm ~ Uniform[0,2π]
    double cos_theta_cm = 2.0 * random_uniform(sim) - 1.0;
    double sin_theta_cm = sqrt(1.0 - cos_theta_cm * cos_theta_cm);
    double phi_cm = 2.0 * PI * random_uniform(sim);
    
    // New neutron velocity in CM frame
    double v_cm_x = v_cm_mag * sin_theta_cm * cos(phi_cm);
    double v_cm_y = v_cm_mag * sin_theta_cm * sin(phi_cm);
    double v_cm_z = v_cm_mag * cos_theta_cm;
    
    // Transform back to lab frame: v' = v_cm + g
    double vx_new = v_cm_x + gx;
    double vy_new = v_cm_y + gy;
    double vz_new = v_cm_z + gz;
    
    // New neutron speed
    *v = sqrt(vx_new*vx_new + vy_new*vy_new + vz_new*vz_new);
}

// Run one neutron history 
void run_history(Simulation* sim) {
    if (sim->total_histories < 5) {
        printf("DEBUG: A=%.1f, alpha_s=%.6f\n", sim->A, sim->alpha_s);
    }

    double v = SOURCE_SPEED_RATIO * sim->v_T;
    int collision_count = 0;
    
    while (collision_count < MAX_COLLISIONS) {
        collision_count++;
        
        // Find speed bin index
        int bin_idx = -1;
        for (int i = 0; i < NUM_INTERVALS; i++) {
            if (v >= sim->speed_bins[i] && v < sim->speed_bins[i + 1]) {
                bin_idx = i;
                break;
            }
        }
        
        // Break conditions
        if (bin_idx < 0 || bin_idx >= NUM_INTERVALS) {
            break;
        }
        
        // Calculate speed-dependent absorption probability using paper's method
        double y = v / sim->v_T;
        double sigma_s_eff = calculate_effective_sigma_s(sim, v);
        double sigma_a = sim->sigma_a0 / v;  // 1/v absorption
        
        double absorption_prob = sigma_a / (sigma_a + sigma_s_eff);
        
        if (random_uniform(sim) < absorption_prob) {
            // Absorption event
            sim->absorption_tally[bin_idx]++;
            sim->total_absorptions++;
            break;
        } else {
            // Scattering event
            sim->scattering_tally[bin_idx]++;
            sim->total_scatterings++;
            
            // Sample target velocity using paper's rejection method
            double V = sample_target_velocity(sim, v);
            
            // Sample collision cosine
            double mu = sample_cosine(sim);
            
            // Use 3D collision physics
            elastic_collision(sim, &v, V, mu);
            
            // Additional break condition for very low speeds
            if (v < 0.01 * sim->v_T) {
                break;
            }
        }
    }
    
    sim->total_histories++;
}

// Run complete simulation
void run_simulation(Simulation* sim) {
    printf("Running %d neutron histories...\n", NUM_HISTORIES);
    printf("Parameters: A=%.2f, K=%.2f, T=%.0f K\n", sim->A, sim->K, sim->T);
    printf("Thermal speed: %.2e m/s\n", sim->v_T);
    
    clock_t start = clock();
    
    for (int i = 0; i < NUM_HISTORIES; i++) {
        if ((i + 1) % 10000 == 0) {
            printf("Completed %d histories\n", i + 1);
        }
        run_history(sim);
    }
    
    clock_t end = clock();
    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Simulation completed in %.2f seconds\n", cpu_time);
    printf("Total scatterings: %d\n", sim->total_scatterings);
    printf("Total absorptions: %d\n", sim->total_absorptions);
    printf("Average collisions per history: %.2f\n", 
           (double)(sim->total_scatterings + sim->total_absorptions) / NUM_HISTORIES);
}

// Calculate flux from tallies
double calculate_flux(Simulation* sim, double* y, double* flux, int* count) {
    double bin_width = sim->speed_bins[1] - sim->speed_bins[0];
    int valid_points = 0;
    
    for (int i = 0; i < NUM_INTERVALS; i++) {
        double bin_center = (sim->speed_bins[i] + sim->speed_bins[i + 1]) / 2.0;
        y[i] = bin_center / sim->v_T;  // Dimensionless speed
        
        if (sim->scattering_tally[i] > 0) {
            // Calculate effective scattering sigma at this speed
            double sigma_es = calculate_effective_sigma_s(sim, bin_center);
            
            // Correct Formula: Flux = Tally / (Sigma * Width * N_histories)
            flux[i] = (double)sim->scattering_tally[i] / 
                     (NUM_HISTORIES * sigma_es * bin_width);
            valid_points++;
        } else {
            flux[i] = 0.0;
        }
    }
    
    *count = valid_points;
    
    // Normalize to thermal peak (approximate range y=0.8-1.2)
    // BUT EXCLUDE SOURCE REGION (y > 8.0)
    double max_flux = 0.0;
    for (int i = 0; i < NUM_INTERVALS; i++) {
        if (y[i] > 0.8 && y[i] < 1.2 && y[i] < 8.0 && flux[i] > max_flux) {
            max_flux = flux[i];
        }
    }
    
    if (max_flux > 0.0) {
        for (int i = 0; i < NUM_INTERVALS; i++) {
            flux[i] /= max_flux;
        }
    }
    
    return max_flux;
}
 

// Maxwell-Boltzmann flux for comparison
double maxwell_flux(double y) {
    return y * y * exp(-y * y);  // Speed distribution (most probable speed normalization)
}

// Print results
void print_results(Simulation* sim) {
    double y[NUM_INTERVALS];
    double flux[NUM_INTERVALS];
    int valid_points;
    
    calculate_flux(sim, y, flux, &valid_points);
    
    printf("\n=== RESULTS ===\n");
    printf("Valid data points: %d\n", valid_points);
    printf("y\tFlux\tMB Flux\n");
    printf("---\t-----\t--------\n");
    
    for (int i = 0; i < NUM_INTERVALS; i += NUM_INTERVALS / 16) {
        if (flux[i] > 0) {
            double mb_flux = maxwell_flux(y[i]);
            printf("%.3f\t%.6f\t%.6f\n", y[i], flux[i], mb_flux);
        }
    }
}

// Save results to file
void save_results(Simulation* sim, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error opening file %s\n", filename);
        return;
    }
    
    double y[NUM_INTERVALS];
    double flux[NUM_INTERVALS];
    int valid_points;
    
    calculate_flux(sim, y, flux, &valid_points);
    
    fprintf(file, "y,flux,mb_flux\n");
    for (int i = 0; i < NUM_INTERVALS; i++) {
        if (flux[i] > 0) {
            double mb_flux = maxwell_flux(y[i]);
            fprintf(file, "%.6f,%.6f,%.6f\n", y[i], flux[i], mb_flux);
        }
    }
    
    fclose(file);
    printf("Results saved to %s\n", filename);
}

void test_carbon_energy_loss() {
    printf("\n=== TESTING CARBON ENERGY LOSS ===\n");
    
    Simulation sim;
    init_simulation(&sim, 12.0, 0.1, 293.0);
    
    // Test with stationary target
    printf("\n1. Stationary target (V=0):\n");
    double v = 10.0 * sim.v_T;
    double V = 0.0;
    double mu = 1.0;  // Head-on
    
    double v_initial = v;
    elastic_collision(&sim, &v, V, mu);
    printf("   Head-on: v_final/v_initial = %.3f (expected: %.3f)\n", 
           v/v_initial, (12.0-1.0)/(12.0+1.0));
    
    // Test average energy loss
    printf("\n2. Average over many collisions (V = v_T):\n");
    v = 10.0 * sim.v_T;
    V = sim.v_T;
    
    double total_energy_ratio = 0.0;
    int n = 10000;
    for (int i = 0; i < n; i++) {
        double v_test = v;
        mu = sample_cosine(&sim);
        elastic_collision(&sim, &v_test, V, mu);
        total_energy_ratio += (v_test*v_test) / (v*v);
    }
    printf("   Average E_final/E_initial = %.4f\n", total_energy_ratio/n);
    printf("   Theoretical (3D): %.4f\n", (144.0 + 1.0)/(13.0*13.0));
}

int main() {
    printf("=== Monte Carlo Neutron Transport Simulation ===\n");
    
    // Hydrogen moderator (A=1), K=0.18
    printf("Case 1: Hydrogen moderator (A=1.0, K=0.18)\n");
    Simulation sim;
    init_simulation(&sim, 1.0, 0.18, 293.0);
    run_simulation(&sim);
    print_results(&sim);
    save_results(&sim, "neutron_flux_hydrogen.csv");
    
    // Carbon moderator (A=12)
    printf("\nCase 2: Carbon moderator (A=12.0, K=0.1)\n");
    test_carbon_energy_loss();
    Simulation sim_carbon;
    init_simulation(&sim_carbon, 12.0, 0.1, 293.0);
    run_simulation(&sim_carbon);
    print_results(&sim_carbon);
    save_results(&sim_carbon, "neutron_flux_carbon.csv");
    
    // Higher absorption case (A=1), K=0.36
    printf("\nCase 3: High absorption (A=1.0, K=0.36)\n");
    Simulation sim_highK;
    init_simulation(&sim_highK, 1.0, 0.36, 293.0);
    run_simulation(&sim_highK);
    print_results(&sim_highK);
    save_results(&sim_highK, "neutron_flux_highK.csv");
    
    printf("\n=== All simulations completed ===\n");
    printf("Output files created:\n");
    printf("  - neutron_flux_hydrogen.csv\n");
    printf("  - neutron_flux_carbon.csv\n");
    printf("  - neutron_flux_highK.csv\n");
    
    return 0;
}
