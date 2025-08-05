/*
    ============================================================================
    pmx_chi.cpp | Susceptibility tensor calculation functions for PicoMax
    ============================================================================

    This file implements the calculation of the susceptibility tensor χ(q,ω) using 
    the empirical pseudopotential method (EPM). The susceptibility tensor describes 
    how a material responds to external electromagnetic fields.

    OVERVIEW OF CALCULATION METHODOLOGY:
    ===================================
    
    1. PHYSICS BACKGROUND:
       - χ(q,ω) describes linear response: P(q,ω) = χ(q,ω) E(q,ω)
       - Calculated using Fermi's Golden Rule for electronic transitions
       - Inter-band transitions: valence → conduction bands
       - Energy conservation: E_c(k) - E_v(k+q) = ħω
    
    2. CALCULATION FLOW:
       a) Setup: Choose broadening function, allocate memory
       b) Electronic structure: Solve H(k)|ψ_n(k)⟩ = E_n(k)|ψ_n(k)⟩
       c) Overlap integrals: ⟨k,c|ĵ e^{i(q+G)·r}|k+q,v⟩
       d) Sum over transitions with energy conservation
       e) Apply physical scaling factors and broadening
    
    3. TENSOR COMPONENTS:
       - LL: Longitudinal-longitudinal (charge density response)
       - TT: Transverse-transverse (current density response)  
       - LT: Longitudinal-transverse (mixed coupling)
       - xyz: Cartesian components for specific directions
    
    4. SPIN-ORBIT COUPLING (SOC) HANDLING:
       - Hamiltonian size: NPW → 2×NPW (spinor basis)
       - Eigenvectors: Include both spin-up/down components
       - Overlap integrals: Sum over spin degrees of freedom
       - Memory allocation: Account for doubled basis size

    Authors: Jungho Mun, Pronoy Das
    ============================================================================
*/


#include "pmx_chi.hpp"      // Function declarations for susceptibility calculations
#include "pmx_epm.hpp"      // Empirical pseudopotential method Hamiltonian
#include "pmx_basis.hpp"    // Basis set and G-vector utilities

/*
    FUNCTION DEFINITIONS
    ===================
*/
namespace pmx{

/*
    ============================================================================
    SECTION 1: BROADENING FUNCTIONS (DIRAC DELTA APPROXIMATIONS)
    ============================================================================
    
    These functions provide different approximations to the Dirac delta function
    δ(E_c - E_v - ħω) used in the susceptibility calculation. The broadening
    accounts for:
    - Finite electronic lifetimes
    - Phonon interactions  
    - Instrumental resolution
    - Computational discretization
    
    The function pointer allows runtime selection of broadening scheme.
*/

/*
    Function pointer for dirac delta function selection
    Allows dynamic choice of broadening method at runtime
*/
inline double (*diracdelta) (double value);

/*
    Lorentzian broadening function
    
    Mathematical form: δ(x) = ε/π/(ε² + x²)
    
    Physical interpretation:
    - Natural line shape for radiative transitions
    - Long tails represent collision broadening
    - FWHM = 2ε (full width at half maximum)
    - Commonly used in optical spectroscopy
    
    Parameters:
    - x: Energy difference (E_c - E_v - ħω)
    - EPSILON: Global smearing parameter (controls line width)
*/
inline double diracdelta_lorentzian(double x){
    return EPSILON/pi/(pow(EPSILON,2)+pow(x,2));
}

/*
    Gaussian broadening function
    
    Mathematical form: δ(x) = (1/(ε√(2π))) × exp(-x²/(2ε²))
    
    Physical interpretation:
    - Models thermal broadening from phonon interactions
    - Instrument resolution function
    - Standard deviation _sum = ε
    - No long tails (compared to Lorentzian)
    - Computationally well-behaved
    
    Parameters:
    - x: Energy difference 
    - EPSILON: Standard deviation of Gaussian distribution
*/
inline double diracdelta_gaussian(double x){
    return exp(-0.5*x*x/(EPSILON*EPSILON))/(EPSILON*sqrt(2.0*pi));
}

/*
    Triangular broadening function
    
    Mathematical form: 
    δ(x) = (4/ε) × (0.5 - |x|/ε)  for |x| < ε/2
         = 0                      otherwise
    
    Physical interpretation:
    - Simple piecewise linear approximation
    - Finite support (goes to zero at boundaries)
    - Computationally efficient
    - Less realistic than Gaussian/Lorentzian
    - Useful for testing and debugging
    
    Parameters:
    - x: Energy difference
    - EPSILON: Half-width of triangle base
*/
inline double diracdelta_triangle(double x){
    if (abs(x)<0.5*EPSILON){
        return 4.0/EPSILON*(0.5-abs(x)/EPSILON);
    }else{
        return 0;
    }
}

/*
    Rectangular (box) broadening function
    
    Mathematical form:
    δ(x) = 1/ε  for |x| < ε/2
         = 0    otherwise
    
    Physical interpretation:
    - Uniform weight within energy window
    - Sharp cutoff at boundaries
    - Simplest possible broadening
    - Can cause artificial discontinuities in results
    - Mainly used for testing purposes
    
    Parameters:
    - x: Energy difference
    - EPSILON: Half-width of rectangular window
*/
inline double diracdelta_rect(double x){
    if (abs(x)<0.5*EPSILON){
        return 1.0/EPSILON;
    }else{
        return 0;
    }
}

/*
    ============================================================================
    SECTION 2: KRAMERS-KRONIG TRANSFORMATION
    ============================================================================
    
    Implements the Hilbert transform to connect real and imaginary parts of
    the susceptibility tensor through causality principles.
    
    PHYSICS BACKGROUND:
    - Kramers-Kronig relations arise from causality constraints
    - Connect Re[χ(ω)] and Im[χ(ω)] through integral transforms
    - Essential for obtaining complete frequency response
    - Allow calculation of one component from the other
    
    MATHEMATICAL FORMS:
    Re[χ(ω)] = (2/π) P∫ ω'Im[χ(ω')]/(ω'²-ω²) dω'
    Im[χ(ω)] = -(2/π) P∫ ω'Re[χ(ω')]/(ω'²-ω²) dω'
    
    where P denotes principal value integral
    
    IMPLEMENTATION:
    - Discrete version for finite frequency grid
    - Careful handling of pole at ω'=ω (principal value)
    - Preserves original data while adding transformed components

    References:
    [1] Landau & Lifshitz, "Electrodynamics of Continuous Media"
    [2] Phys. Rev. B 74, 035101 (2006)
    
    Parameters:
    - ReX: Real parts [input/output array]
    - ImX: Imaginary parts [input/output array]
    - w: Frequency grid points
    - dw: Frequency spacing
*/


inline Eigen::Vector3cd socCurrent(Eigen::Vector3d Ki,      // bra  G‑vector
                                   Eigen::Vector3d Kj,      // ket  G‑vector
                                   pmx::mater      mat_params,     // material params (λ_S, λ_A maps)
                                   std::vector<Eigen::Vector3d> atomic_pos){
    using namespace pmx;
    Eigen::Vector3d q  = Kj  - Ki;               //   G_j − G_i

    /* λ_S(q²), λ_A(q²) — same lookup as in HamiltonianEPM() */
    double q_sq_norm = q.squaredNorm();
    double lambda_A_Ry = pmx::get_soc_form_factor(mat_params.Ua_SO_Ry, q_sq_norm);
    double lambda_S_Ry = pmx::get_soc_form_factor(mat_params.Us_SO_Ry, q_sq_norm);
    double lambda_S_eV = lambda_S_Ry * Ry2eV; // Ry2eV is global
    double lambda_A_eV = lambda_A_Ry * Ry2eV;  // there is a scale factor for SOC in eV for easy parameter tuning
    /* structure factor  _sum_s e^{ i q·τ_s } /N */
    // Structure factor part for SOC
    double cos_sum = 0.0;
    double sin_sum = 0.0;
    if (NATOM > 0) { // NATOM is global
        for (const auto& atom_pos_a_units : atomic_pos) { // atomic_pos are in units of [a]
            // q_vec_soc_dimless is G_j - G_i (dimless multiples of 2pi/a)
            // atom_pos_a_units is tau_s (dimless multiples of a)
            // dot product is dimensionless
            cos_sum += std::cos(2.0 * pi * q.dot(atom_pos_a_units));
            sin_sum += std::sin(2.0 * pi * q.dot(atom_pos_a_units));
        }
    }

    if (NATOM > 0) { // NATOM is global
        cos_sum /=NATOM;
        sin_sum /=NATOM;
    }
    
    // Scalar part of VSO from Mathematica: (Lambda_S * cos_qT + Lambda_A * sin_qT)
    // double VSO_scalar_struct_part = lambda_S_eV * cos_sum + im * lambda_A_eV * sin_sum;
    std::complex<double> VSO= lambda_S_eV * cos_sum + im * lambda_A_eV * sin_sum;

    // Compute the SOC matrix element: ⟨G_i| im * VSO * q × σ |G_j⟩
    // 
    // The cross product q × σ gives a vector operator:
    // (q × σ)_x = q_y * σ_z - q_z * σ_y
    // (q × σ)_y = q_z * σ_x - q_x * σ_z  
    // (q × σ)_z = q_x * σ_y - q_y * σ_x
    //
    // For the matrix elements between plane wave states |G_i⟩ and |G_j⟩,
    // the spatial part gives δ(G_j - G_i), so this is only non-zero when G_i = G_j
    // The spin matrix elements for the cross product components are:
    
    Eigen::Vector3cd soc_cross_product;
    
    // For a spinor basis |G,↑⟩, |G,↓⟩, the matrix elements are:
    // x-component: ⟨↑| (q_y σ_z - q_z σ_y) |↓⟩ = q_y * 1 - q_z * (-im) = q_y + im * q_z
    //              ⟨↓| (q_y σ_z - q_z σ_y) |↑⟩ = q_y * (-1) - q_z * (im) = -q_y - im * q_z
    soc_cross_product(0) = im * VSO * (q(1) + im * q(2));  // x-component: im*VSO*(q_y + i*q_z)
    
    // y-component: ⟨↑| (q_z σ_x - q_x σ_z) |↓⟩ = q_z * 1 - q_x * 1 = q_z - q_x
    //              ⟨↓| (q_z σ_x - q_x σ_z) |↑⟩ = q_z * 1 - q_x * (-1) = q_z + q_x
    soc_cross_product(1) = im * VSO * (q(2) - q(0));       // y-component: im*VSO*(q_z - q_x)
    
    // z-component: ⟨↑| (q_x σ_y - q_y σ_x) |↓⟩ = q_x * (-im) - q_y * 1 = -im * q_x - q_y
    //              ⟨↓| (q_x σ_y - q_y σ_x) |↑⟩ = q_x * (im) - q_y * 1 = im * q_x - q_y
    soc_cross_product(2) = im * VSO * (-im * q(0) - q(1)); // z-component: im*VSO*(-i*q_x - q_y)

    return soc_cross_product;
}





void kramerskronigtransform(double *ReX, double *ImX, double *w, double dw){

    // STEP 1: Preserve original real part for later use
    // We need the original values because we modify ReX in the first loop
    double ReX1 [NFREQ];
    for (int f=0; f<NFREQ; f++)
        ReX1[f] = ReX[f];

    // Integration prefactor from Kramers-Kronig relation: (2/π) × dω
    double ds = dw*2.0/pi;

    // STEP 2: Generate real part from imaginary part
    // Re[χ(ω)] = (2/π) P∫ ω' Im[χ(ω')] / (ω'² - ω²) dω'
    for (int f=0; f<NFREQ; f++){
        double ReXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            // Skip pole at ω'=ω (principal value integral)
            if (s!=f)
                ReXf += w[s]*ImX[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ReX[f] += ds*ReXf;  // Add transformed component to original
    }

    // STEP 3: Generate imaginary part from original real part  
    // Im[χ(ω)] = -(2/π) P∫ ω' Re[χ(ω')] / (ω'² - ω²) dω'
    // Note: We use ReX1 (original) not the modified ReX from step 2
    for (int f=0; f<NFREQ; f++){
        double ImXf = 0.0;
        for (int s=0; s<NFREQ; s++){
            // Skip pole at ω'=ω (principal value integral)
            if (s!=f)
                ImXf += w[s]*ReX1[s]/(w[s]*w[s]-w[f]*w[f]);
        }
        ImX[f] += ds*ImXf;  // Add transformed component to original
    }

    return;
}

/*
    calculate the susceptibility tensor matrix using the current-current correlation function

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
    [3] 
*/
void chi_tensor(env &dat){
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;

    // use indirect function calls for diracdelta function
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;
            break;
        case 2:
            diracdelta = &diracdelta_triangle;
            break;
        case 3:
            diracdelta = &diracdelta_rect;
            break;
    }
    NTSR = 3;

    // scalefactor
    double SCALEFACTOR = 16.0*dat.lat.bz_volume*pow(pi,5)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    double SF_SOC = 1; // scale factor for spin-orbit coupling susceptibility tensor
    

    // initialize susceptibility tensor matrix elements
    dat.ReXij = new double *****[NQ];
    dat.ImXij = new double *****[NQ];
    for (int q=0; q<NQ; q++){
        dat.ReXij[q] = new double ****[NTSR];
        dat.ImXij[q] = new double ****[NTSR];
        for (int i=0; i<NTSR; i++){
            dat.ReXij[q][i] = new double ***[NTSR];
            dat.ImXij[q][i] = new double ***[NTSR];
            for (int j=0; j<NTSR; j++){
                dat.ReXij[q][i][j] = new double **[NEPS];
                dat.ImXij[q][i][j] = new double **[NEPS];
                for (int m=0; m<NEPS; m++){
                    dat.ReXij[q][i][j][m] = new double *[NEPS];
                    dat.ImXij[q][i][j][m] = new double *[NEPS];
                    for (int n=0; n<=m; n++){
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    int NBAND_C [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC x NPW]
    C_k = new std::complex<double> **[NKPT];
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // vector{k}
        Eigen::Vector3d K = dat.lat.K[k];

        // Hamiltonian at k
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // compute eigenvalue problem at k
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        E_k[k] = new double [NBAND];
        for (int c=0; c<NBAND; c++){
            E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
            if (E_k[k][c]<=0) {
                NBAND_C[k] = c;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_C[k]];
        for (int c=0; c<NBAND_C[k]; c++){
            const int NPW_SOC = 2 * NPW;
            // C_k[k][c] = new std::complex<double> [NPW_SOC];
            C_k_up[k][c] = new std::complex<double> [NPW];
            C_k_down[k][c] = new std::complex<double> [NPW];
            // for (int p=0; p<NPW_SOC; p++){
            //     C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
            // }// loop over p
            for (int p=0; p<NPW; p++){
                C_k_up[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(2*p);
                C_k_down[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(2*p+1);
            }// loop over p
        }// loop over c
    }
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    // std::complex<double> *****oint;
    std::complex<double> *****ointup;
    std::complex<double> *****ointdown;
    std::complex<double> *****ointupdown;
    std::complex<double> *****ointdownup;
    ointup = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    ointdown = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    ointupdown = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    ointdownup = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        ointup[k] = new std::complex<double> ***[3];
        ointdown[k] = new std::complex<double> ***[3];
        ointupdown[k] = new std::complex<double> ***[3];
        ointdownup[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            ointup[k][i] = new std::complex<double> **[NEPS];
            ointdown[k][i] = new std::complex<double> **[NEPS];
            ointupdown[k][i] = new std::complex<double> **[NEPS];
            ointdownup[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                ointup[k][i][m] = new std::complex<double> *[NBAND_C[k]];
                ointdown[k][i][m] = new std::complex<double> *[NBAND_C[k]];
                ointupdown[k][i][m] = new std::complex<double> *[NBAND_C[k]];
                ointdownup[k][i][m] = new std::complex<double> *[NBAND_C[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements..." << std::endl;
    int NBAND_V [NKPT]; // number of valence bands
    std::complex<double> ***C_kq; // eigenvector at k+q
    double **E_kq; // eigenvalue at k
    E_kq = new double *[NKPT];
    C_kq = new std::complex<double> **[NKPT];
    for (int k=0; k<NKPT; k++){
        E_kq[k] = new double [NBAND];
    }
    // #pragma omp parallel
    for (int q=0; q<NQ; q++){
        std::cout << "  (" << (q+1) << "/" << NQ << ")..." << std::endl; // (q/NQ)
        Eigen::Vector3d Q = dat.lat.Q[q];

        std::cout << "    Solving for overlap integrals..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for
        for (int k=0; k<NKPT; k++){
            
            // Hamiltonian at k+q
            Eigen::Vector3d K = dat.lat.K[k];
            Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K+Q,dat.lat.atomic,dat.mat);   

            // compute eigenvalue problem at k+q
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_kq(H);
            for (int v=0; v<NBAND; v++){
                E_kq[k][v] = double(eigsolver_kq.eigenvalues()(v))-dat.energyoffset;
                if (E_kq[k][v]>0) {
                    NBAND_V[k] = v;
                    break;
                }
            }
            C_kq[k] = new std::complex<double> *[NBAND_V[k]];
            for (int v=0; v<NBAND_V[k]; v++){
                const int NPW_SOC = 2 * NPW;
                // C_kq[k][v] = new std::complex<double> [NPW_SOC];
                C_kq_up[k][v] = new std::complex<double> [NPW];
                C_kq_down[k][v] = new std::complex<double> [NPW];
                // for (int p=0; p<NPW_SOC; p++){
                //     C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
                // }
                for (int p=0; p<NPW; p++){
                    C_kq_up[k][v][p] = eigsolver_kq.eigenvectors().col(v)(2*p);
                    C_kq_down[k][v][p] = eigsolver_kq.eigenvectors().col(v)(2*p+1);
                }
            }
            
            // overlap integrals
            for (int m=0; m<NEPS; m++){
                std::vector<Eigen::Vector3cd> uvec_m = uvec_LT(Q+dat.lat.G[m]);
                
                // Pre-compute active G indices for this m ONCE per m
                std::vector<int> active_G_indices;
                for (int p=0; p<NPW; p++){
                    if (dat.lat.loci[m][p] != -1){
                        active_G_indices.push_back(p);
                    }
                }
                
                // Pre-calculate SOC contribution once per (k,m) pair for active indices only
                std::vector<Eigen::Vector3cd> v_soc_cache(active_G_indices.size());
                for (size_t i_active = 0; i_active < active_G_indices.size(); i_active++){
                    int p = active_G_indices[i_active];
                    int loci_p = dat.lat.loci[m][p];
                    v_soc_cache[i_active] = socCurrent(
                        dat.lat.G[p], dat.lat.G[loci_p], dat.mat, dat.lat.atomic);
                }


                for (int i=0; i<NTSR; i++){
                    for (int c=0; c<NBAND_C[k]; c++){
                        ointup[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        ointdown[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        ointupdown[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        ointdownup[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        for (int v=0; v<NBAND_V[k]; v++){
                            // u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r}j_0|k+q,c>
                            ointup[k][i][m][c][v] = 0;
                            ointupdown[k][i][m][c][v] = 0;
                            ointdownup[k][i][m][c][v] = 0;
                            ointdown[k][i][m][c][v] = 0;
                            for (size_t i_active = 0; i_active < active_G_indices.size(); i_active++){
                                int p = active_G_indices[i_active];
                                int loci_p = dat.lat.loci[m][p];                          
                                Eigen::Vector3cd v_orb = (dat.lat.G[p] + K + Q/2 + dat.lat.G[m]/2).cast<std::complex<double>>();
                                std::complex<double> soc_contribution = 0.0 * uvec_m[i].dot(v_soc_cache[i_active]);
                                std::complex<double> orb_contribution = uvec_m[i].dot(v_orb);
                                ointup[k][i][m][c][v] += conj(C_k_up[k][c][p]) * C_kq_up[k][v][loci_p] * orb_contribution;
                                ointdown[k][i][m][c][v] -= conj(C_k_down[k][c][p]) * C_kq_down[k][v][loci_p] * orb_contribution;

                                // ointupdown[k][i][m][c][v] += conj(C_k[k][c][2*p]) * C_kq[k][v][2*loci_p+1] * soc_contribution;
                                // ointdownup[k][i][m][c][v] -= conj(C_k[k][c][2*p+1]) * C_kq[k][v][2*loci_p] * soc_contribution;
                            }
                        }//loop over v
                    }//loop over c
                }//loop over i
            }//loop over m
        }//loop over k
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        
        std::cout << "    Solving for Xijmn..." << std::endl;
        time_0 = std::chrono::system_clock::now();
        #pragma omp parallel for collapse(4)
        for (int i=0; i<NTSR; i++){for (int j=0; j<NTSR; j++){for (int m=0; m<NEPS; m++){for (int n=0; n<=m; n++){
            for (int f=0; f<NFREQ; f++){
                double tmp_real_1 = 0;
                double tmp_imag_1 = 0;
                double tmp_real_2 = 0;
                double tmp_imag_2 = 0;
                // sum over k,c,v
                std::complex<double> Oij;
                for (int k=0; k<NKPT; k++){
                    for (int c=0; c<NBAND_C[k]; c++){
                        for (int v=0; v<NBAND_V[k]; v++){
                            int c_spin = c % 2;  // 0=up, 1=down
                            int v_spin = v % 2;
                            double dE = E_k[k][c]-E_kq[k][v];
                            Oij = ointup[k][i][m][c][v] * conj(ointup[k][j][n][c][v]) + ointdown[k][i][m][c][v] * conj(ointdown[k][j][n][c][v]);
                            // if (c_spin == v_spin) {
                            //     if (c_spin == 0) {
                            //         // Both spin-up
                            //         Oij = ointup[k][i][m][c][v] * conj(ointup[k][j][n][c][v]);
                            //     } else {
                            //         // Both spin-down
                            //         Oij = ointdown[k][i][m][c][v] * conj(ointdown[k][j][n][c][v]);
                            //     }
                            // } 
                            // else {
                            //     if (c_spin == 0) {
                            //         // up-down spin
                            //         Oij = ointupdown[k][i][m][c][v] * conj(ointupdown[k][j][n][c][v]);
                            //     } else {
                            //         // down-up spin
                            //         Oij = ointdownup[k][i][m][c][v] * conj(ointdownup[k][j][n][c][v]);
                            //     }
                            // }

                            tmp_imag_1 += dat.lat.KW[k] * Oij.real() * (*diracdelta)(dE-dat.freq[f]);
                            tmp_real_1 -= dat.lat.KW[k] * Oij.imag() * (*diracdelta)(dE-dat.freq[f]);
                        }
                    }   
                }

                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    tmp_imag_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                    tmp_real_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                }

                dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (tmp_imag_1);
                dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (tmp_real_1);
            }
        }}}}
        #pragma omp barrier
        time_1 = std::chrono::system_clock::now();
        elapsed_time = time_1-time_0;
        std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;


        
        // kramer-kronig transform
        if(dat.kk){
            std::cout << "    Kramers-Kronig transform..." << std::endl;
            time_0 = std::chrono::system_clock::now();
            #pragma omp parallel for collapse(4)
            for (int i=0; i<NTSR; i++){ for (int j=0; j<NTSR; j++){ for (int m=0; m<NEPS; m++){ for (int n=0; n<=m; n++){
                kramerskronigtransform(dat.ReXij[q][i][j][m][n],dat.ImXij[q][i][j][m][n],dat.freq,dat.dfreq);
            }}}}
            #pragma omp barrier
            time_1 = std::chrono::system_clock::now();
            elapsed_time = time_1-time_0;
            std::cout << "    Elapsed time: " << elapsed_time.count() << " s" << std::endl;
        }
        

        // // free dynamically allocated memory
        // for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int v=0; v<NBAND_V[k]; v++){
        //     delete [] oint[k][i][m][v];
        // }}}}
        // for (int k=0; k<NKPT; k++){
        //     for (int c=0; c<NBAND_C[k]; c++){
        //         delete [] C_kq[k][c];
        //     }
        //     delete [] C_kq[k];
        // }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] ointup[k][i][m];
                delete [] ointdown[k][i][m];
                delete [] ointupdown[k][i][m];
                delete [] ointdownup[k][i][m];
            }
            delete [] ointup[k][i];
            delete [] ointdown[k][i];
            delete [] ointdownup[k][i];
            delete [] ointupdown[k][i];
        }
        delete [] ointup[k];
        delete [] ointdown[k];
        delete [] ointdownup[k];
        delete [] ointupdown[k];
    }
    delete [] ointup;
    delete [] ointdown;
    delete [] ointupdown;
    delete [] ointdownup;

    for (int k=0; k<NKPT; k++){
        for (int v=0; v<NBAND_V[k]; v++){
            delete [] C_k[k][v];
        }
        delete [] C_k[k];
        delete [] E_k[k];
        delete [] E_kq[k];
    }
    delete [] E_k;
    delete [] E_kq;
    delete [] C_k;
    delete [] C_kq;

    return;
}

} /* namespace pmx */
