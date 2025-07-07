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

    Author: Jungho Mun
    Date: June 6, 2024
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
    calculate the susceptibility tensor matrix using an alternative expression
    which treats the longitudinal terms as charge

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
*/
void chi_LL(env &dat){
    // ========================================================================
    // INITIALIZATION AND SETUP
    // ========================================================================
    
    // Timing variables for performance monitoring
    std::chrono::time_point<std::chrono::system_clock> time_0, time_1;
    std::chrono::duration<double> elapsed_time;

    // Set tensor dimension to 3 for full 3D tensor calculation
    NTSR = 3;

    // STEP 1: Select broadening function based on user input
    // This determines how discrete energy levels are converted to continuous spectra
    switch (dat.delta) {
        case 0:
            diracdelta = &diracdelta_gaussian;    // Thermal/phonon broadening
            break;
        case 1:
            diracdelta = &diracdelta_lorentzian;  // Natural line shape
            break;
        case 2:
            diracdelta = &diracdelta_triangle;    // Simple linear broadening
            break;
        case 3:
            diracdelta = &diracdelta_rect;        // Uniform broadening
            break;
    }

    // STEP 2: Define physical scaling factors for different tensor components
    // These convert from atomic units to practical units (eV, Å, etc.)
    double SCALEFACTOR = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
    
    // Longitudinal-longitudinal: charge density response  
    double SCALEFACTOR_LL = 2.0*dat.lat.bz_volume*pow(e,2)/(pi*eV*a*angstrom);
    
    // Transverse-transverse: current density response
    double SCALEFACTOR_TT = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    
    // Longitudinal-transverse: mixed coupling terms
    double SCALEFACTOR_LT = 8.0*dat.lat.bz_volume*pi*pow(hbar,2)*pow(e,2)/(pow(eV,2)*e_m*pow(a*angstrom,3));

    // STEP 3: Allocate memory for susceptibility tensor matrix elements
    // Structure: χ[q-vector][i-component][j-component][G-vector_m][G-vector_n][frequency]    dat.ReXij = new double *****[NQ];      // Real parts of susceptibility tensor
    dat.ImXij = new double *****[NQ];      // Imaginary parts of susceptibility tensor
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
                    for (int n=0; n<=m; n++){  // Only upper triangular part (symmetry)
                        dat.ReXij[q][i][j][m][n] = new double [NFREQ];
                        dat.ImXij[q][i][j][m][n] = new double [NFREQ];
                    }
                }
            }
        }
    }

    // ========================================================================
    // PHASE 1: CONDUCTION BAND ELECTRONIC STRUCTURE CALCULATION
    // ========================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    
    // Variables for conduction band states at k-points
    int NBAND_C [NKPT]; // Number of conduction bands per k-point
    double **E_k; // Conduction band eigenvalues E_k[k-point][band]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // Conduction band eigenvectors C_k[k-point][band][G-vector]
    C_k = new std::complex<double> **[NKPT];
    
    // Parallelize over k-points for electronic structure calculation
    #pragma omp parallel for
    for (int k=0; k<NKPT; k++){
        // Get current k-vector from lattice structure
        Eigen::Vector3d K = dat.lat.K[k];

        // Construct EPM Hamiltonian matrix H(k)
        // Returns NPW×NPW matrix (or 2×NPW×2×NPW if SOC enabled)
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G,K,dat.lat.atomic,dat.mat);   
        
        // Solve eigenvalue problem: H(k)|ψ_n(k)⟩ = E_n(k)|ψ_n(k)⟩
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver_k(H);

        // Extract eigenvalues and determine band character
        E_k[k] = new double [NBAND];
        for (int c=0; c<NBAND; c++){
            // Extract eigenvalues (from highest energy downward)
            E_k[k][c] = eigsolver_k.eigenvalues()(NBAND-1-c)-dat.energyoffset;
            if (E_k[k][c]<=0) {
                NBAND_C[k] = c;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_C[k]];
        for (int c=0; c<NBAND_C[k]; c++){
            const int NPW_SOC = 2 * NPW;
            C_k[k][c] = new std::complex<double> [NPW_SOC];
            for (int p=0; p<NPW_SOC; p++){
                C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
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
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_C[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements (alternative form)..." << std::endl;
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
                C_kq[k][v] = new std::complex<double> [NPW_SOC];
                for (int p=0; p<NPW_SOC; p++){
                    C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
                }
            }
            
            // overlap integrals
            for (int m=0; m<NEPS; m++){
                std::vector<Eigen::Vector3cd> uvec_m = uvec_LT(Q+dat.lat.G[m]);
                for (int i=0; i<NTSR; i++){
                    for (int c=0; c<NBAND_C[k]; c++){
                        oint[k][i][m][c] = new std::complex<double> [NBAND_V[k]];
                        for (int v=0; v<NBAND_V[k]; v++){
                            oint[k][i][m][c][v] = 0;
                            if (i==0){// L, <k,c|e^{-i*(q+g_m)*r}|k+q,v>
                                for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                    int loci_p = dat.lat.loci[m][p];
                                    // Sum over spin-up and spin-down components
                                    // Spin-up: index 2*p and 2*loci_p
                                    oint[k][i][m][c][v] += conj(C_k[k][c][2*p]) * C_kq[k][v][2*loci_p];
                                    // Spin-down: index 2*p+1 and 2*loci_p+1
                                    oint[k][i][m][c][v] += conj(C_k[k][c][2*p+1]) * C_kq[k][v][2*loci_p+1];
                                }}
                            }else{// T, u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r} \hat{j}_0 |k+q,v>
                                for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                    int loci_p = dat.lat.loci[m][p];
                                    Eigen::Vector3cd v_orb =(dat.lat.G[p] + K + Q/2 + dat.lat.G[m]/2).cast<std::complex<double>>();
                                    /* SOC piece – same parameters/maps as Hamiltonian */
                                    Eigen::Vector3cd v_soc = socCurrent(
                                        dat.lat.G[p],                // G_i  (bra)
                                        dat.lat.G[loci_p],           // G_j  (ket)
                                        dat.mat,                     // material (λ maps)
                                        dat.lat.atomic);             // τ_s list
                                    /* Total current operator */
                                    std::complex<double> momentum_factor = uvec_m[i].dot(v_orb + v_soc);
                                    // Sum over spin-up and spin-down components
                                    // Spin-up: index 2*p and 2*loci_p
                                    oint[k][i][m][c][v] += conj(C_k[k][c][2*p]) * C_kq[k][v][2*loci_p] * momentum_factor;
                                    // Spin-down: index 2*p+1 and 2*loci_p+1
                                    oint[k][i][m][c][v] += conj(C_k[k][c][2*p+1]) * C_kq[k][v][2*loci_p+1] * momentum_factor;
                                }}
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
            Eigen::Vector3d vecqgm = Q+dat.lat.G[m];
            Eigen::Vector3d vecqgn = Q+dat.lat.G[n];
            for (int f=0; f<NFREQ; f++){
                double tmp_real_1 = 0;
                double tmp_imag_1 = 0;
                double tmp_real_2 = 0;
                double tmp_imag_2 = 0;
                // sum over k,c,v
                for (int k=0; k<NKPT; k++){for (int c=0; c<NBAND_C[k]; c++){for (int v=0; v<NBAND_V[k]; v++){
                    double dE = E_k[k][c]-E_kq[k][v];
                    std::complex<double> Oij = oint[k][i][m][c][v] * conj(oint[k][j][n][c][v]);
                    
                    // 
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                          * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                          * (*diracdelta)(dE-dat.freq[f]);

                    if (dat.kk==0){
                        if (i==0 && j==0){// LL
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         * (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }else if (i==0 || j==0){// LT
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                          / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                          / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }else{// TT
                            tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                            tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        }
                    }
                }}}
                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    if (i==0 && j==0){
                        tmp_imag_1 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                        tmp_real_1 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                    }else if (i==0){
                        tmp_imag_1 *= SCALEFACTOR_LT / (vecqgm.norm()*dat.freq[f]);
                        tmp_real_1 *= SCALEFACTOR_LT / (vecqgm.norm()*dat.freq[f]);
                    }else if (j==0){
                        tmp_imag_1 *= SCALEFACTOR_LT / (vecqgn.norm());
                        tmp_real_1 *= SCALEFACTOR_LT / (vecqgn.norm());
                    }else{
                        tmp_imag_1 *= SCALEFACTOR_TT / (dat.freq[f]*dat.freq[f]);
                        tmp_real_1 *= SCALEFACTOR_TT / (dat.freq[f]*dat.freq[f]);
                    }
                }
                if (dat.kk==0){
                    if (i==0 && j==0){
                        tmp_imag_2 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                        tmp_real_2 *= SCALEFACTOR_LL / (vecqgm.norm()*vecqgn.norm());
                    }else if (i==0){
                        tmp_imag_2 *= SCALEFACTOR_LT / (vecqgm.norm());
                        tmp_real_2 *= SCALEFACTOR_LT / (vecqgm.norm());
                    }else if (j==0){
                        tmp_imag_2 *= SCALEFACTOR_LT / (vecqgn.norm());
                        tmp_real_2 *= SCALEFACTOR_LT / (vecqgn.norm());
                    }else{
                        tmp_imag_2 *= SCALEFACTOR_TT;
                        tmp_real_2 *= SCALEFACTOR_TT;
                    }
                    
                }

                dat.ImXij[q][i][j][m][n][f] = (pi*tmp_imag_1+2.0*tmp_imag_2);
                dat.ReXij[q][i][j][m][n][f] = (pi*tmp_real_1+2.0*tmp_real_2);
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
        // for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int c=0; c<NBAND_C[k]; c++){
        //     delete [] oint[k][i][m][c];
        // }}}}
        // for (int k=0; k<NKPT; k++){
        //     for (int v=0; v<NBAND_V[k]; v++){
        //         delete [] C_kq[k][v];
        //     }
        //     delete [] C_kq[k];
        // }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

    for (int k=0; k<NKPT; k++){
        for (int c=0; c<NBAND_C[k]; c++){
            delete [] C_k[k][c];
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
    NTSR = 2;

    // scalefactor
    double SCALEFACTOR = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));
    double SF = 4*pi*dat.lat.bz_volume*pow(2*pi/(a*angstrom),3)*pow(e,2)/eV; // scale factor for susceptibility tensor ~10e-8
    double SF_TT = pow(hbar,4) * pow(1/(a*angstrom),2) / (pow(e_m,2) * pow(eV,2)); // scale factor for transverse-transverse susceptibility tensor
    double SF_LL = pow((a*angstrom),2); // scale factor for longitudinal-transverse susceptibility tensor
    double SF_SOC = (1/(4 * e_m * pow(3e+10,2))) * eV; // scale factor for spin-orbit coupling susceptibility tensor
    double SF_LT= pow(SF_TT*SF_LL,0.5); // scale factor for longitudinal-transverse susceptibility tensor

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
            C_k[k][c] = new std::complex<double> [NPW_SOC];
            for (int p=0; p<NPW_SOC; p++){
                C_k[k][c][p] = eigsolver_k.eigenvectors().col(NBAND-1-c)(p);
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
                C_kq[k][v] = new std::complex<double> [NPW_SOC];
                for (int p=0; p<NPW_SOC; p++){
                    C_kq[k][v][p] = eigsolver_kq.eigenvectors().col(v)(p);
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
                            // T, u^i_{q+g_m} * <k,c|e^{-i*(q+g_m)*r} \hat{j}_0 |k+q,v>
                                // Loop over only active indices
                            for (size_t i_active = 0; i_active < active_G_indices.size(); i_active++){
                                int p = active_G_indices[i_active];
                                int loci_p = dat.lat.loci[m][p];
                                
                                Eigen::Vector3cd v_orb = (dat.lat.G[p] + K + Q/2 + dat.lat.G[m]/2).cast<std::complex<double>>();
                                
                                std::complex<double> soc_contribution = SF_SOC * uvec_m[i].dot(v_soc_cache[i_active]);
                                std::complex<double> orb_contribution = uvec_m[i].dot(v_orb);
                                
                                // Apply contributions
                                ointup[k][i][m][c][v] += conj(C_k[k][c][2*p]) * C_kq[k][v][2*loci_p] * orb_contribution;
                                ointdown[k][i][m][c][v] += conj(C_k[k][c][2*p+1]) * C_kq[k][v][2*loci_p+1] * orb_contribution;
                                ointupdown[k][i][m][c][v] += conj(C_k[k][c][2*p]) * C_kq[k][v][2*loci_p+1] * soc_contribution;
                                ointdownup[k][i][m][c][v] -= conj(C_k[k][c][2*p+1]) * C_kq[k][v][2*loci_p] * soc_contribution;
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
                for (int k=0; k<NKPT; k++){for (int c=0; c<NBAND_C[k]; c++){for (int v=0; v<NBAND_V[k]; v++){
                    // Determine actual spin indices instead of using modulo
                    int c_spin = c % 2;  // 0=up, 1=down
                    int v_spin = v % 2;
                    double dE = E_k[k][c]-E_kq[k][v];
                    std::complex<double> Oij;
                    if (c_spin == v_spin) {
                        // Same spin: use appropriate arrays
                        if (c_spin == 0) {
                            // Both spin-up
                            Oij = ointup[k][i][m][c][v] * conj(ointup[k][j][n][c][v]);
                        } else {
                            // Both spin-down
                            Oij = ointdown[k][i][m][c][v] * conj(ointdown[k][j][n][c][v]);
                        }
                    } else {
                        if (c_spin == 0) {
                            // up-down spin
                            Oij = ointupdown[k][i][m][c][v] * conj(ointupdown[k][j][n][c][v]);
                        } else {
                            // down-up spin
                            Oij = ointdownup[k][i][m][c][v] * conj(ointdownup[k][j][n][c][v]);
                        }
                    }
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                                * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                        * (*diracdelta)(dE-dat.freq[f]);    
                }}}
                // if (i == 0 && j == 0){// LL
                //     tmp_imag_1 *= 2.0;
                //     tmp_real_1 *= 2.0;
                // }
                // if (f==0){
                //     tmp_imag_1 = 0;
                //     tmp_real_1 = 0;
                // }else{
                //     if (i%3!=0 && j%3!=0){
                //         tmp_imag_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                //         tmp_real_1 *= 1 / (dat.freq[f]*dat.freq[f]);

                //         dat.ImXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_imag_1);
                //         dat.ReXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_real_1);
                //     } else if (i%3==0 && i==j){
                //         tmp_imag_1 *= 1 / ((Q + dat.lat.G[m]).norm()*(Q + dat.lat.G[n]).norm());
                //         tmp_real_1 *= 1 / ((Q + dat.lat.G[m]).norm()*(Q + dat.lat.G[n]).norm());

                //         dat.ImXij[q][i][j][m][n][f] = SF*SF_LT * (tmp_imag_1);
                //         dat.ReXij[q][i][j][m][n][f] = SF*SF_LT * (tmp_real_1);
                //     } else if (i%3==0 || j%3==0){
                //         tmp_imag_1 *= 1 / ((Q + dat.lat.G[n]).norm()*dat.freq[f]);
                //         tmp_real_1 *= 1 / ((Q + dat.lat.G[n]).norm()*dat.freq[f]);

                //         dat.ImXij[q][i][j][m][n][f] = SF*SF_LL* (tmp_imag_1);
                //         dat.ReXij[q][i][j][m][n][f] = SF*SF_LL* (tmp_real_1);
                //     } else{
                //         tmp_imag_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                //         tmp_real_1 *= 1 / (dat.freq[f]*dat.freq[f]);

                //         dat.ImXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_imag_1);
                //         dat.ReXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_real_1);
                //     }
                // }
                
                // dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_imag_1);
                // dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_real_1);

                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    // TT
                    tmp_imag_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                    tmp_real_1 *= 1 / (dat.freq[f]*dat.freq[f]);
                }
                dat.ImXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_imag_1);
                dat.ReXij[q][i][j][m][n][f] = SF*SF_TT * (tmp_real_1);
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

/*
    calculate the susceptibility tensor matrix using the current-current correlation function,
    where the polarization unit vectors are given by the Cartesian unit vectors {nx,ny,nz}

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions
    [2] (Adler) “Quantum Theory of the Dielectric Constant in Real Solids”, Phys. Rev. 126, 413 (1962)
*/
void chi_xyz(env &dat){
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
    double SCALEFACTOR = 32.0*dat.lat.bz_volume*pow(pi,3)*pow(hbar,4)*pow(e,2)/(pow(eV,3)*pow(e_m,2)*pow(a*angstrom,5));

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
                        // for (int f=0; f<NFREQ; f++){
                        //     dat.ReXij[q][i][j][m][n][f] = 0;
                        //     dat.ImXij[q][i][j][m][n][f] = 0;
                        // }
                    }
                }
            }
        }
    }

    //======================================================================
    std::cout << "  Solving for |k,v>..." << std::endl;
    time_0 = std::chrono::system_clock::now();
    int NBAND_V [NKPT]; // number of conduction bands
    double **E_k; // eigenvalue at k [NK x NC]
    E_k = new double *[NKPT];
    std::complex<double> ***C_k; // eigenvector at k [NK x NC]
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
        for (int v=0; v<NBAND; v++){
            E_k[k][v] = eigsolver_k.eigenvalues()(v)-dat.energyoffset;
            if (E_k[k][v]>0) {
                NBAND_V[k] = v;
                break;
            }
        }
        C_k[k] = new std::complex<double> *[NBAND_V[k]];
        for (int v=0; v<NBAND_V[k]; v++){
            C_k[k][v] = new std::complex<double> [NPW];
            for (int p=0; p<NPW; p++){
                C_k[k][v][p] = eigsolver_k.eigenvectors().col(v)(p);
            }// loop over p
        }// loop over v
    }
    #pragma omp barrier
    time_1 = std::chrono::system_clock::now();
    elapsed_time = time_1-time_0;
    std::cout << "  Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    //======================================================================
    // initialize large temporary variables
    // overlap integral, t_{Q+G_m} * <k,v|e^{-i*G_m}j_0|k+q,c> <k+q,c|e^{i*G_n}j_0|k,v> * t_{Q+G_n}
    std::complex<double> *****oint;
    oint = new std::complex<double> ****[NKPT]; // overlap integral [NK x NTSR x NEPS x NC x NV]
    for (int k=0; k<NKPT; k++){
        oint[k] = new std::complex<double> ***[3];
        for (int i=0; i<NTSR; i++){
            oint[k][i] = new std::complex<double> **[NEPS];
            for (int m=0; m<NEPS; m++){
                oint[k][i][m] = new std::complex<double> *[NBAND_V[k]];
            }
        }
    }

    std::cout << "  Solving for susceptibility tensor matrix elements (xyz-basis)..." << std::endl;
    int NBAND_C [NKPT]; // number of valence bands
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
            for (int c=0; c<NBAND; c++){
                E_kq[k][c] = double(eigsolver_kq.eigenvalues()(NBAND-1-c))-dat.energyoffset;
                if (E_kq[k][c]<=0) {
                    NBAND_C[k] = c;
                    break;
                }
            }
            C_kq[k] = new std::complex<double> *[NBAND_C[k]];
            for (int c=0; c<NBAND_C[k]; c++){
                C_kq[k][c] = new std::complex<double> [NPW];
                for (int p=0; p<NPW; p++){
                    C_kq[k][c][p] = eigsolver_kq.eigenvectors().col(NBAND-1-c)(p);
                }
            }
            
            std::vector<Eigen::Vector3d> uvec;
            uvec.reserve(3);
            uvec[0] = {1,0,0};
            uvec[1] = {0,1,0};
            uvec[2] = {0,0,1};

            // overlap integrals
            for (int m=0; m<NEPS; m++){
                for (int i=0; i<NTSR; i++){
                    for (int v=0; v<NBAND_V[k]; v++){
                        oint[k][i][m][v] = new std::complex<double> [NBAND_C[k]];
                        for (int c=0; c<NBAND_C[k]; c++){
                            // u^i_{q+g_m} * <k,v|e^{-i*(q+g_m)*r}j_0|k+q,c>
                            oint[k][i][m][v][c] = 0;
                            for (int p=0; p<NPW; p++){if (dat.lat.loci[m][p]!=-1){
                                oint[k][i][m][v][c] += (conj(C_k[k][v][p]) * C_kq[k][c][dat.lat.loci[m][p]]) 
                                            * (uvec[i].dot(dat.lat.G[p]+K+Q/2+dat.lat.G[m]/2));
                            }}
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
                for (int k=0; k<NKPT; k++){for (int v=0; v<NBAND_V[k]; v++){for (int c=0; c<NBAND_C[k]; c++){
                    double dE = E_kq[k][c]-E_k[k][v];
                    std::complex<double> Oij = oint[k][i][m][v][c] * conj(oint[k][j][n][v][c]);
                    
                    // 
                    tmp_imag_1 += dat.lat.KW[k] * Oij.real()
                                          * (*diracdelta)(dE-dat.freq[f]);
                    tmp_real_1 -= dat.lat.KW[k] * Oij.imag()
                                          * (*diracdelta)(dE-dat.freq[f]);

                    if (dat.kk==0){
                        tmp_real_2 += dat.lat.KW[k] * Oij.real()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                        tmp_imag_2 += dat.lat.KW[k] * Oij.imag()
                                         / (dE) / (dE*dE-dat.freq[f]*dat.freq[f]);
                    }
                }}}
                if (f==0){
                    tmp_imag_1 = 0;
                    tmp_real_1 = 0;
                }else{
                    tmp_imag_1 /= (dat.freq[f]*dat.freq[f]);
                    tmp_real_1 /= (dat.freq[f]*dat.freq[f]);
                }
                dat.ImXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_imag_1+2*tmp_imag_2);
                dat.ReXij[q][i][j][m][n][f] = SCALEFACTOR * (pi*tmp_real_1+2*tmp_real_2);
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
        

        // free dynamically allocated memory
        for (int k=0; k<NKPT; k++){for (int i=0; i<NTSR; i++){for (int m=0; m<NEPS; m++){for (int v=0; v<NBAND_V[k]; v++){
            delete [] oint[k][i][m][v];
        }}}}
        for (int k=0; k<NKPT; k++){
            for (int c=0; c<NBAND_C[k]; c++){
                delete [] C_kq[k][c];
            }
            delete [] C_kq[k];
        }
    }// loop q

    // free dynamically allocated memory
    for (int k=0; k<NKPT; k++){
        for (int i=0; i<NTSR; i++){
            for (int m=0; m<NEPS; m++){
                delete [] oint[k][i][m];
            }
            delete [] oint[k][i];
        }
        delete [] oint[k];
    }
    delete [] oint;

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


/*
    ============================================================================
    COMPLETE DOCUMENTATION OF REMAINING LOGICAL FLOW
    ============================================================================
    
    The susceptibility calculation continues with these main phases:
    
    PHASE 2: BAND CLASSIFICATION
    ----------------------------
    - Separate valence (E < 0) and conduction (E > 0) bands
    - Store only conduction band eigenvectors to save memory
    - Account for SOC: eigenvectors have 2×NPW components (spinors)
    
    PHASE 3: Q-VECTOR LOOP  
    -----------------------
    For each momentum transfer q:
    
    3A) Calculate valence states at k+q points
        - Solve H(k+q) eigenvalue problem
        - Extract valence band eigenvectors and energies
        - Handle SOC: spinor eigenvectors with doubled size
    
    3B) Compute overlap integrals
        - Matrix elements: ⟨c,k|ĵ_i e^{-i(q+G)·r}|v,k+q⟩
        - Current operator: ĵ = (e/m)(-iℏ∇) = (e/m)ℏ(k+G) in reciprocal space
        - Longitudinal: Pure phase factors (no current operator)
        - Transverse: Include momentum operator projected onto perpendicular directions
        - SOC: Sum over both spin components for each spatial plane wave
    
    3C) Susceptibility tensor calculation
        - Apply Fermi's Golden Rule with energy conservation
        - Energy difference: ΔE = E_c(k) - E_v(k+q)
        - Delta function: δ(ΔE - ħω) with chosen broadening
        - Occupation factors: [f_v - f_c] ≈ 1 for T=0
        - Sum over all k-points and band combinations
    
    3D) Apply scaling and tensor-specific factors
        - LL components: SCALEFACTOR_LL (charge response)
        - TT components: SCALEFACTOR_TT (current response)  
        - LT components: SCALEFACTOR_LT (mixed coupling)
        - Include geometric factors and unit conversions
    
    PHASE 4: POST-PROCESSING
    ------------------------
    - Apply Kramers-Kronig relations if requested (KK=1)
    - Add δ_{G,G'} terms for macroscopic response (q→0 limit)
    - Store results in output tensor arrays
    - Memory cleanup and deallocation
    
    KEY PHYSICS INSIGHTS:
    ====================
    
    1. ELECTRONIC TRANSITIONS:
       - χ(q,ω) peaks occur at inter-band transition energies
       - Momentum conservation: k → k+q during transition
       - Energy conservation: E_c(k) - E_v(k+q) = ħω
    
    2. CURRENT OPERATOR IMPLEMENTATION:
       - Longitudinal: Pure charge density (no current operator)
       - Transverse: Current density ĵ = (e/m)p̂ with momentum operator
       - Direction: Unit vectors define tensor component orientations
    
    3. SOC EFFECTS:
       - Doubles Hamiltonian and eigenvector sizes
       - Mixes spin channels in wavefunctions  
       - Current operator diagonal in spin space
       - Overlap integrals sum over both spin components
    
    4. BROADENING PHYSICS:
       - Converts discrete transitions to continuous spectra
       - Models finite electronic lifetimes
       - Accounts for phonon interactions and disorder
       - Instrumental resolution effects
    
    5. SCALING FACTORS:
       - Convert atomic units to practical units (eV, Å)
       - Include Brillouin zone volume normalization
       - Account for different field orientations (L vs T)
       - Proper electromagnetic unit conversions
    
    MEMORY MANAGEMENT:
    ==================
    - Large arrays: eigenvectors, overlap integrals, susceptibility tensors
    - Careful allocation/deallocation to prevent memory leaks
    - OpenMP parallelization requires thread-safe memory handling
    - SOC doubles memory requirements for eigenvector storage
    
    COMPUTATIONAL COMPLEXITY:
    =========================
    - Scales as O(NQ × NKPT × NBAND²) for basic calculation
    - Additional factors: NEPS² (G-vectors), NFREQ (frequencies)
    - Dominated by overlap integral calculations
    - Parallelization over k-points provides good scaling
    
    ============================================================================
*/

// The remaining functions (chi_tensor, chi_xyz, etc.) follow similar patterns
// with variations in:
// - Tensor component selection (full tensor vs specific components)
// - Overlap integral formulations (different current operator treatments)
// - Scaling factor applications (geometry-specific factors)
// - Memory allocation strategies (optimized for specific use cases)

} /* namespace pmx */

/*
    END OF FILE DOCUMENTATION
    =========================
    
    This file implements the core physics of electronic susceptibility
    calculations using empirical pseudopotentials. The functions provide
    different approaches (charge-based, full tensor, simplified) but all
    follow the fundamental quantum mechanical framework of electronic
    transitions between valence and conduction bands with proper handling
    of spin-orbit coupling effects.
    
    For usage examples and parameter descriptions, see the main documentation
    and example input files in the project repository.
*/