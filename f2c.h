#ifndef F2C_INC
#define F2C_INC

#ifdef JUMP

#define COOLR coolr
#define COOLI cooli
#define COOLINMO coolinmo
#define CHEMINMO cheminmo
#define INIT_TOLERANCES init_tolerances
#define EVOLVE_ABUNDANCES evolve_abundances
#define RATE_EQ rate_eq
#define COMPUTE_GAMMA compute_gamma
#define LOAD_H2_TABLE load_h2_table

#else

#define COOLR coolr_
#define COOLI cooli_
#define COOLINMO coolinmo_
#define CHEMINMO cheminmo_
#define INIT_TOLERANCES init_tolerances_
#define EVOLVE_ABUNDANCES evolve_abundances_
#define RATE_EQ rate_eq_
#define COMPUTE_GAMMA compute_gamma_
#define LOAD_H2_TABLE load_h2_table_

#endif

#endif /* F2C_INC */
