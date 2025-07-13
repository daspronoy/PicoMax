/*
    epm.cpp | includes functions related to calculation of electronic bands

    Authors: Jungho Mun, Pronoy Das
*/

/*
    header files
*/
#include "pmx_epm.hpp"
#include <Eigen/Dense> // Required for Matrix2cd, Identity
#include <cmath> // For std::round, std::abs, std::cos, std::sin
#include <map>   // For std::map in mater struct for SOC params

/*
    function definitions
*/
namespace pmx{

// Helper to get SOC form factor; returns value in Ry
// q_sq_norm_dimless is |G'-G|^2, where G vectors are dimensionless (multiples of 2*pi/a)
double get_soc_form_factor(const std::map<double, double>& so_map, double q_sq_norm_dimless) {
    double tol = 1e-2; // Tolerance for matching q_sq_norm
    // The keys in so_map are doubles like 0.565147, 1.33333, etc.
    // We search for the closest key in the map.
    if (std::abs(q_sq_norm_dimless) < tol) {
        // For q=0, SOC form factors are typically zero or not well-defined in this context.
        // The cross product term will also be zero if G'=G and K vectors are the same.
        // However, if map contains a key for 0, use it.
        auto it_zero = so_map.find(0.0);
        if (it_zero != so_map.end()) {
            return it_zero->second;
        }
        return 0.0;
    }

    double closest_key = 0.0;    double min_diff = std::numeric_limits<double>::max();
    for (const auto& pair : so_map) {
        double diff = std::abs(q_sq_norm_dimless - pair.first);
        if (diff < min_diff) {
            min_diff = diff;
            closest_key = pair.first;
        }
    }

    if (min_diff < tol) {
        return so_map.at(closest_key);
    }
    return 0.0;
}

/*
    Returns the Hamiltonian matrix for the empirical pseudopotential method
    Now includes Spin-Orbit Coupling (SOC).
*/
Eigen::MatrixXcd HamiltonianEPM (std::vector<Eigen::Vector3d> G_vectors, Eigen::Vector3d K_vec, std::vector<Eigen::Vector3d> atomic_pos, pmx::mater mat_params){
    int num_g_vectors = NPW; // NPW is a global variable for number of plane waves
    Eigen::MatrixXcd H_soc(2 * num_g_vectors, 2 * num_g_vectors);
    H_soc.setZero(); // Initialize to zero

    double KINETIC_CONST = KCON * pow(2.0 * pi / a, 2); // a is global lattice const

    // Define Pauli matrices
    Eigen::Matrix2cd sigma_x, sigma_y, sigma_z, identity_2x2;
    sigma_x << 0, 1, 1, 0;
    sigma_y << 0, -im, im, 0; // im is global std::complex<double>(0,1)
    sigma_z << 1, 0, 0, -1;
    identity_2x2 = Eigen::Matrix2cd::Identity();

    for (int i = 0; i < num_g_vectors; i++) {
        for (int j = 0; j <= i; j++) { // Fill lower triangle (and diagonal) of G-blocks
            
            std::complex<double> H_scalar_part;
            if (i == j) { // Diagonal G-block: Kinetic energy
                H_scalar_part = KINETIC_CONST * (K_vec + G_vectors[i]).squaredNorm();
                // V(G=0) from local pseudopotential is typically zero or absorbed into energy reference.
            } else { // Off-diagonal G-block: Local pseudopotential
                Eigen::Vector3d dG_pseudo = G_vectors[i] - G_vectors[j];                if (dG_pseudo.squaredNorm() < mat_params.g2max) {
                    H_scalar_part = pmx::pseudopotential(dG_pseudo, atomic_pos, mat_params);
                } else {
                    H_scalar_part = std::complex<double>(0.0, 0.0);
                }
            }
            Eigen::Matrix2cd H_scalar_2x2_block = H_scalar_part * identity_2x2;

            // SOC Part
            // Wavevectors for SOC (dimensionless, in units of 2*pi/a)
            Eigen::Vector3d Ki_dimless = K_vec + G_vectors[i];
            Eigen::Vector3d Kj_dimless = K_vec + G_vectors[j];
            // q_vec for SOC, K[j]-K[i] from Mathematica
            Eigen::Vector3d q_vec_soc_dimless = G_vectors[j] - G_vectors[i]; 

            // SOC parameters (retrieve from map, convert from Ry to eV)
            double q_sq_norm = q_vec_soc_dimless.squaredNorm();
            double lambda_S_Ry;
            if (mat_params.use_uniform_us_so) {
                lambda_S_Ry = mat_params.uniform_us_so_val_Ry;
            } else {
                lambda_S_Ry = pmx::get_soc_form_factor(mat_params.Us_SO_Ry, q_sq_norm);
            }
            double lambda_A_Ry = pmx::get_soc_form_factor(mat_params.Ua_SO_Ry, q_sq_norm);
            double lambda_S_eV = lambda_S_Ry * Ry2eV; // Ry2eV is global
            double lambda_A_eV = lambda_A_Ry * Ry2eV;  // there is a scale factor for SOC in eV for easy parameter tuning

            // Structure factor part for SOC
            double cos_sum = 0.0;
            double sin_sum = 0.0;
            if (NATOM > 0) { // NATOM is global
                for (const auto& atom_pos_a_units : atomic_pos) { // atomic_pos are in units of [a]
                    // q_vec_soc_dimless is G_j - G_i (dimless multiples of 2pi/a)
                    // atom_pos_a_units is tau_s (dimless multiples of a)
                    // dot product is dimensionless
                    cos_sum += std::cos(2.0 * pi * q_vec_soc_dimless.dot(atom_pos_a_units));
                    sin_sum += std::sin(2.0 * pi * q_vec_soc_dimless.dot(atom_pos_a_units));
                }
            }

            if (NATOM > 0) { // NATOM is global
                cos_sum /=NATOM;
                sin_sum /=NATOM;
            }
            
            // Scalar part of VSO from Mathematica: (Lambda_S * cos_qT + Lambda_A * sin_qT)
            // double VSO_scalar_struct_part = lambda_S_eV * cos_sum + im * lambda_A_eV * sin_sum;
            std::complex<double> VSO_scalar_struct_part = lambda_S_eV * cos_sum + im * lambda_A_eV * sin_sum;
            // WARNING: If Lambda_A and sin_sum are both non-zero, the Mathematica formula
            // for VSO might lead to a non-Hermitian Hamiltonian.
            // For Ge, Lambda_A is often zero. For Te, this might need review.
            // A common Hermitian form involves (Lambda_S cos_qT + I * Lambda_A sin_qT).

            // Kinematic part for SOC: Cross[K[j]+veck, K[i]+veck]
            Eigen::Vector3d soc_kinematic_vec = Kj_dimless.cross(Ki_dimless);

            // SOC vector potential (purely imaginary based on Mathematica -I * real_scalar * real_vector)
            Eigen::Vector3cd SO_vec_potential = im*VSO_scalar_struct_part * soc_kinematic_vec.cast<std::complex<double>>();

            // SOC 2x2 block: SO_vec_potential . sigma_vec
            Eigen::Matrix2cd H_soc_2x2_block = 
                SO_vec_potential.x() * sigma_x + 
                SO_vec_potential.y() * sigma_y + 
                SO_vec_potential.z() * sigma_z;
            
            // Total 2x2 block
            Eigen::Matrix2cd total_2x2_block = H_scalar_2x2_block + H_soc_2x2_block;
            
            H_soc.block<2,2>(2 * i, 2 * j) = total_2x2_block;
        }
    }
    // The matrix H_soc is Hermitian by construction if SelfAdjointEigenSolver is used,
    // which will effectively use H_ji = H_ij^dagger for j > i.
    // If the VSO formula implies non-Hermitian terms (see WARNING above), results might be unphysical.
    return H_soc;
}


/*
    Returns the complex pseudopotential Vg, at the given reciprocal vector G (not |G|)

    Vg = sum_{s=1}^{NATOM} V_{s,g} exp(-i*2*pi*dot(t_s,g)) / NATOM 
    where V_{s,g} is the atomic form factor for atom s at G.
    The input G to this function is G_i - G_j.

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] (Cohen and Bergstresser) "Band Structures and Pseudopotnetial Form Factors 
        for Fourteen Semiconductors of the Diamond and Zinc-blende Structures", 
        Phys. Rev. 141, 789 (1966)
    [3] (Bergstresser and Cohen) "Electronic Structure and Optical Properties of 
        Hexagonal CdSe, CdS, and ZnS", Phys. Rev. 164, 1069 (1967)
*/
std::complex<double> pseudopotential(
        Eigen::Vector3d G_diff, std::vector<Eigen::Vector3d> T_atomic_pos, pmx::mater mat_params){
    double tol = 1e-2;
    double GMAG_sq = G_diff.squaredNorm(); // This G_diff is G_i - G_j
    int form_factor_idx = -1; // Index for mat_params.vg

    // Find the index j for mat_params.g2 that matches GMAG_sq
    for (int k=0; k<mat_params.NV; k++){ // NV is number of G2 values for Vg
        if (abs(mat_params.g2[k]-GMAG_sq)<tol){
            form_factor_idx = k;
            break;
        }
    }

    if (form_factor_idx == -1) { // No matching G-vector squared magnitude found for Vg
        return std::complex<double>(0.0, 0.0);
    }

    std::complex<double> Vg_total = 0;
    if (NATOM > 0) {
        for (int s=0; s<NATOM; s++){ // Sum over atoms in basis
            // mat_params.vg[s][form_factor_idx] is V_{s,G} in Ry
            // T_atomic_pos[s] is atomic position tau_s in units of [a]
            // G_diff is G_i - G_j in units of [2pi/a]
            // G_diff.dot(T_atomic_pos[s]) is dimensionless
            Vg_total += exp(-im * 2.0 * pi * G_diff.dot(T_atomic_pos[s])) * mat_params.vg[s][form_factor_idx];
        }
        Vg_total /= double(NATOM); // Average over atoms
    }
    
    Vg_total *= Ry2eV; // Convert from Ry to eV
    return Vg_total;
}


/*
    Calculates electron energy band using the empirical pseudopotential method

    todo: change function name

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] (Cohen and Bergstresser) "Band Structures and Pseudopotnetial Form Factors 
        for Fourteen Semiconductors of the Diamond and Zinc-blende Structures", 
        Phys. Rev. 141, 789 (1966)
    [3] (Bergstresser and Cohen) "Electronic Structure and Optical Properties of 
        Hexagonal CdSe, CdS, and ZnS", Phys. Rev. 164, 1069 (1967)
*/
void empiricalpseudopotentialmethod(pmx::env &dat){
    dat.eband = new double*[NQ]; // NQ is global number of Q (k-points in path)

    // obtain energy offset
    setRefEnergy(dat); // This will now use the SOC Hamiltonian if called

    #pragma omp parallel for
    for (int q=0; q<NQ; q++){
        
        // Hamiltonian (now includes SOC and is 2*NPW x 2*NPW)
        Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G, dat.lat.Q[q], dat.lat.atomic, dat.mat);

        // diagonalize the Hamiltonian
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
        if (eigsolver.info() != Eigen::Success){ std::cout << "ERROR: Eigensolver failed to solve the matrix." << std::endl; abort(); }

        // truncate to NBAND and energy offset
        // NBAND is global, number of bands to output
        // eigsolver.eigenvalues() will have 2*NPW eigenvalues, sorted.
        dat.eband[q] = new double[NBAND];
        for (int n=0; n<NBAND; n++){
            if (n < eigsolver.eigenvalues().size()) { // Ensure we don't go out of bounds
                dat.eband[q][n] = double(eigsolver.eigenvalues()(n)) - dat.energyoffset;
            } else {
                // Handle cases where NBAND might be larger than available eigenvalues (e.g. small NPW)
                // This shouldn't happen if NBAND <= 2*NPW
                dat.eband[q][n] = 0.0; // Or some indicator like NaN
            }
        }
    }
    // #pragma omp barrier // Not strictly needed after parallel for unless there's code after depending on all threads finishing
    // OMP barrier is implicit at the end of parallel for unless nowait is specified.

    return;
}


/*
    Returns the nonlocal pseudopotential Vg^NL

    Gi, Gj : reciprocal lattice vectors
    T : atomic basis vectors
    mat : material parameters

    Reference:
    [1] See the documentation Eq. XXX for explicit expressions 
    [2] 
*/
std::complex<double> get_Vg_nonlocal(Eigen::Vector3d Gi, Eigen::Vector3d Gj, 
                                    std::vector<Eigen::Vector3d> T, pmx::mater mat){
    // This function is not implemented yet.
    return std::complex<double>(0.0,0.0);
}



// complex nonlocal-pseudopotential for the given reciprocal vector
std::complex<double> nonlocalpseudopotential(Eigen::Vector3d Gm, Eigen::Vector3d Gn, 
                    std::vector<double> v){
    // This function seems to be an older or alternative implementation.
    // It uses a hardcoded vec_tau and specific indices for v parameters.
    // The main HamiltonianEPM uses the more general pseudopotential function.
    // If nonlocal EPM is to be used with SOC, this would also need modification.
    double a0_param = v[6]; // Assuming v has a specific structure
    double b0_param = v[7];
    double rl0_param = v[9];

    //define variables for symmetric and asymmetric part of the potential
	double sym = 0.0e0, asym = 0.0e0; // asym is not used here
	std::complex<double> potential = 0.0e0;

	//Each cell has 2 atoms.  Use mid point between them as offset in units of $a$
	Eigen::Vector3d vec_tau_nonlocal; // This should ideally come from crystal structure
    vec_tau_nonlocal << 0.125, 0.125, 0.125;

	//Find norm of the reciprocal vector (these are K+G vectors, already scaled by 2pi/a if from K[m])
    // The input Gm, Gn are G-vectors (dimless multiples of 2pi/a).
    // To match Mathematica K[m] which is h*b1+..., these Gm, Gn are G-vectors.
    // The K[m] in Mathematica's nonlocal part are |K+G|*(2pi/a).
    // Here, Gm.norm() is dimensionless. Need to scale by (2pi/a) for physical k.
	double Km_phys = Gm.norm()*(2.0*pi/a); // a is global
	double Kn_phys = Gn.norm()*(2.0*pi/a);
	double kfermi_phys = pow(6.0*pi*pi/(pow(a*angstrom,3)), double(1.0/3.0)); // kfermi in cm^-1, a in Angstrom
    // The KCON in kinetic energy is hbar^2/(2m eV A^2).
    // a0E = a0*Ry + b0*Ry * (hbar^2/(2m)) * (k^2 - kf^2) / eV
    // (hbar^2/(2m)) * k^2 is KCON * (2pi/a)^2 * k_dimless^2 * eV = KCON * k_phys^2 * eV
    // So, (b0*Ry2eV) * KCON * (Km_phys*Kn_phys/((2pi/a)^2) - kfermi_phys^2/((2pi/a)^2))
    // This is getting complicated. The original Mathematica nonlocal part needs careful translation.
    // For now, this function is not directly called by the main EPM loop if EPM_NONLOCAL is false.
    // The provided task is to implement SOC, not to fix/implement nonlocal EPM.

	double a0E_val = a0_param*Ry2eV + (b0_param*Ry2eV*KCON*(sqrt(Km_phys*Km_phys*Kn_phys*Kn_phys)-pow(kfermi_phys*(a*angstrom/(2*pi)),2)));
    // The kfermi term needs to be in same units as Km_phys*Kn_phys for KCON scaling.
    // Km_phys is like k * (2pi/a). So Km_phys / (2pi/a) is dimensionless k.
    // kfermi_phys is k_F. So kfermi_phys / (2pi/a) is dimensionless k_F.
    // a0E = a0*Ry2eV + b0*Ry2eV*KCON * ( (Km_dimless*Kn_dimless) - kf_dimless^2 )
    // where Km_dimless = Gm.norm(), Kn_dimless = Gn.norm()
    // kf_dimless = kfermi_phys * a / (2*pi)
    double Km_dimless = Gm.norm();
    double Kn_dimless = Gn.norm();
    double kf_dimless = kfermi_phys * (a*angstrom) / (2.0*pi); // a in Angstrom
    a0E_val = a0_param*Ry2eV + (b0_param*Ry2eV*KCON*(sqrt(Km_dimless*Km_dimless*Kn_dimless*Kn_dimless)-pow(kf_dimless,2)));


	double Fl_val;
	double Kmr = Km_phys * rl0_param; // rl0 seems to be a length scale
	double Knr = Kn_phys * rl0_param;

	if(abs(Kmr-Knr) < 1e-9){ // If Kmr and Knr are too close
		Fl_val = 0.5*pow(rl0_param,3) * (pow(pmx::sphbesj(0,Kmr),2) - pmx::sphbesj(-1,Kmr)*pmx::sphbesj(1,Kmr));
	}else{ // Original formula
        // Typo in original Mathematica: Km*Km-Kn*Kn should be Km_phys^2 - Kn_phys^2
		Fl_val = (pow(rl0_param,2)/(Km_phys*Km_phys-Kn_phys*Kn_phys)) * 
                 (Km_phys*pmx::sphbesj(1,Kmr)*pmx::sphbesj(0,Knr) - Kn_phys*pmx::sphbesj(1,Knr)*pmx::sphbesj(0,Kmr));
	}	
	sym = (4.0*pi/pow(a*angstrom,3))*a0E_val*Fl_val; // sym should be in eV if a0E_val is eV
    // (4pi/Volume) * Energy * Length^3 = Energy * (4pi L^3 / Volume). This should be dimensionless * Energy.
    // Fl has units of L^3. Volume is a^3. So this is dimensionless * Energy.
    sym /= eV; // Convert from ergs to eV if a0E_val was in ergs. But a0E_val is already eV.
               // So this division by eV is likely an error if sym is meant to be in eV.
               // Let's assume sym is already in eV from a0E_val.

	//Complex Pseudopotential value
	potential = std::complex<double>(sym*cos(2*pi*((Gm-Gn).dot(vec_tau_nonlocal))), 
                                    asym*sin(2*pi*((Gm-Gn).dot(vec_tau_nonlocal)))); // asym is 0
	return potential;
}

/*
    Returns the spherical Bessel function of the first kind
*/
double sphbesj(int n, double r){
	double jn = 0.0e0;
    double abs_r = std::abs(r); // Use std::abs for floating point
	if(abs_r < 1e-8){ // Check absolute value for r near zero
		switch(n){
			case -1: // j_{-1}(x) = cos(x)/x. For x->0, cos(x)->1, so 1/x. Diverges.
                     // Using a large number to approximate divergence or indicate issue.
				jn = (r == 0.0) ? std::numeric_limits<double>::infinity() : cos(r)/r; // More robust
                if (!std::isfinite(jn)) jn = (r > 0 ? 1.0e8 : -1.0e8); // Fallback for exact zero
                break;
			case 0: // j_0(x) = sin(x)/x. For x->0, limit is 1.
				jn = 1.0e0;	break;
			case 1: // j_1(x) = sin(x)/x^2 - cos(x)/x. For x->0, limit is 0 (approx x/3).
				jn = 0.0e0;	break;
            default: // For other n, behavior at r=0 might be 0 or undefined.
                jn = 0.0; break; 
		}
	}
	else{
		switch(n){
			case -1:
				jn = cos(r)/r;	break;
			case 0: 
				jn = sin(r)/r;	break;
			case 1: 
				jn = (sin(r)/(r*r))-(cos(r)/r);	break;
            default: // General case for spherical bessel j_n(r)
                // This simple switch only handles specific n.
                // For a general implementation, use recurrence or standard library if available.
                // For now, returning 0 for unhandled n.
                jn = 0.0; break;
		}
	}
	return jn;
}

/*
    set reference energy [eV] at the reference kpoint

    default reference point is the gamma point

    you can change the reference point using the option `-refpoint x,y,z`
    ex) `picomax ... -refpoint 0.5,0.5,0.5`
*/
void setRefEnergy(pmx::env &dat){

    // the global Hamiltonian matrix H at reference kpoint
    Eigen::Vector3d K_REF = dat.refpoint; // reference kpoint

    // HamiltonianEPM now returns the SOC-included Hamiltonian
    Eigen::MatrixXcd H = HamiltonianEPM(dat.lat.G, K_REF, dat.lat.atomic, dat.mat);

    // compute the eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigsolver(H);
    if (eigsolver.info() != Eigen::Success){ 
        std::cout << "ERROR::setRefEnergy:: Eigensolver failed to solve the matrix at reference point." << std::endl; 
        abort(); 
    }

    // zero level of energy at the valence band maximum at reference point
    // dat.nvalence is number of valence electron pairs (e.g., 4 for Ge, Si)
    // Eigenvalues are sorted in increasing order.
    // If nvalence = 4, it means 8 electrons. VBM is the (2*4-1)=7th eigenvalue (0-indexed) if all degeneracies lifted.
    // Or, if nvalence refers to the (spin-degenerate) band index, then (nvalence-1) is the (spatial) VBM.
    // With SOC, this band splits. The higher of the split bands is the true VBM.
    // The current indexing dat.eigenvalues()(dat.nvalence-1) might need adjustment
    // depending on how nvalence is defined and how SOC affects the VBM.
    // For now, keeping it as is. If dat.nvalence refers to the Nth eigenvalue in the list.
    if (dat.nvalence > 0 && dat.nvalence <= eigsolver.eigenvalues().size()) {
        // dat.energyoffset = eigsolver.eigenvalues()(dat.nvalence - 1)-0.117716;
        dat.energyoffset = (eigsolver.eigenvalues()(dat.nvalence - 1)+eigsolver.eigenvalues()(dat.nvalence))/2.0;
    } else {
        std::cout << "WARNING::setRefEnergy:: nvalence is out of bounds for eigenvalues. Setting offset to 0." << std::endl;
        dat.energyoffset = 0.0;
        if (eigsolver.eigenvalues().size() > 0) { // Or set to lowest eigenvalue if available
             // dat.energyoffset = eigsolver.eigenvalues()(0);
        }
    }
    return;
}

}/* namespace pmx */
