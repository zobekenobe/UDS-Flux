#include"udf.h"
#include"sg.h"

// Constant values for simulation
real L 	 	= 0.25;			// channel length (m)
real W 	 	= 0.040;		// channel width (m)
real h      = 0.0008;		// channel height (m)
real m_A0   = 0.005837;		// solute mass fraction at inlet (kg/ kg) 0.6% for NaCl
real Re     = 278.0; 		// feed Reynolds number dimensionless
real del_P  = 5e+5;   		// Pressure Drop across membrane
real Lp     = 1.80556e-11;	// permeability of membrane (conductance (m/Pa.s))
real Ps	 	= 5.44782e-5; 	// overall solute permeability (m/Pa.s)
real sigma0 = 0.957308625;	// initial reflection coefficient 0.9797
real J_v0   = 2.77321e-6;  	// J_v0 (m/s) initial guess v_vel through membrane 
real Dh 	= 4.0*(30.0/31.0);  // Hydraulic Diameter
real R_prime= 0.0;

// user defined memory C_UDMI
int x_nondim = 0; 		// x/h dimensionaless distance from the wall
int y_nondim = 1; 		// y/h dimensionaless distance fron the wall
int gamma    = 2;		// concentration polarization layer thickness (m)
int delta 	 = 3;		// concentration polarization layer thickness (eq 11)
int k_delta	 = 4; 		// mass transfer coefficient (eq 13)

int m_Ap   	 = 5; 		// Solute mass fraction at the wall apparently not needed

int J_v 	 = 6;		// Flux at the wall
int R 		 = 7;		// Dynamic Membrane Resistance
int sigma	 = 8;		// reflection coefficient 
int schmidt	 = 9;		// Schmidt Number

// user defined scalar transport C_UDSI
int m_A = 0;

// Initialize the solute concentration, flux, membrane resistance, sigma
DEFINE_INIT(initial_setup, domain_pointer)
{
	Thread* thread_pointer;
	cell_t cellt;

	thread_loop_c(thread_pointer, domain_pointer)
	{
		begin_c_loop(cellt, thread_pointer)
		{
			C_UDSI(cellt, thread_pointer, m_A) 		= m_A0;
			C_UDMI(cellt, thread_pointer, J_v) 		= J_v0;
			C_UDMI(cellt, thread_pointer, R)   		= 0.98;
			C_UDMI(cellt, thread_pointer, m_Ap)		= 0.0;
			C_UDMI(cellt, thread_pointer, sigma)	= sigma0;
			R_prime = 1.0 - C_UDMI(cellt, thread_pointer, R);
		}
		end_c_loop(cellt, thread_pointer)
	}
}


DEFINE_PROPERTY(Nacl_density, cellt,  thread_pointer)
{
	return 997.1 + 694.*C_UDSI(cellt, thread_pointer, m_A);
}

DEFINE_DIFFUSIVITY(Nacl_diffusivity, cellt, thread_pointer, position_index)
{
	 if (C_UDSI(cellt, thread_pointer, m_A) > 0.006) 
	 	return C_R(cellt, thread_pointer)*1.45e-09;
	 else
	 	return C_R(cellt, thread_pointer)*1.61e-09 * (1 - 14. * C_UDSI(cellt, thread_pointer, m_A));
}

/* defining the fluid viscosity in the solvent*/
DEFINE_PROPERTY(NaCl_viscosity, cellt, thread_pointer)
{
	return 0.89e-03 * ( 1 + 1.63 * C_UDSI(cellt, thread_pointer, m_A));
}

/* Initial solute concentration */
/* Initial solute concentration through the side inlet */
DEFINE_PROFILE(Nacl_mA_inlet, thread_pointer, position_index)
{
	face_t facet;

	begin_f_loop(facet, thread_pointer)
	{
		F_PROFILE(facet, thread_pointer, position_index) = m_A0;
	}
	end_f_loop(facet, thread_pointer)
}

/*Initial solute concentration through the bottom wall */

DEFINE_PROFILE(Nacl_mA_bottom, thread_pointer, position_index)
{
	Thread *t0 = thread_pointer->t0;
	cell_t c0;
	cell_t cellt;
	face_t facet;

	real A[ND_ND];
	real ds; 			// distance between cell and the face centroids
	real es[ND_ND];
	real A_by_es;
	real dr0[ND_ND];
	real source;

	real dphi[ND_ND];
	real k;

	real temp1;
	real temp2; 
	real m_Aw; 		   // m_A at the wall

	real osmotic_m_Aw;
	real osmotic_m_Ap;
	real nothing;


	begin_f_loop(facet, thread_pointer)
	{
		BOUNDARY_FACE_GEOMETRY(facet, thread_pointer, A, ds, es, A_by_es, dr0);

		c0 = F_C0(facet, thread_pointer);
		k = C_UDSI_DIFF(c0, t0, m_A)/C_R(c0, t0);

		if(NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(m_A)))) source = 0.0;
			
		else
			BOUNDARY_SECONDARY_GRADIENT_SOURCE(source, SV_UDSI_G(m_A), dphi, es, A_by_es, k);

		temp1 = k * A_by_es/ds;
		temp2 = C_UDMI(c0,t0,J_v)* C_R(c0, t0) * NV_MAG(A);
		m_Aw  = (temp1 * C_UDSI(c0, t0, m_A) - source + temp2 * C_UDMI(c0, t0,m_Ap))/(temp1 + temp2);
		//nothing = C_UDSI(cellt, thread_pointer, m_A);

		// ======================================================================== Calculating the Flux =======================================================================
		osmotic_m_Aw = 805.1e+05 * F_UDSI(facet,thread_pointer,m_A);
		osmotic_m_Ap = 805.1e+05 * R_prime * F_UDSI(facet,thread_pointer,m_A);

		C_UDMI(c0,t0,J_v) = Lp * (del_P - C_UDMI(c0,t0,sigma) * (osmotic_m_Aw - osmotic_m_Ap));

		// ======================================================================== Calculating the Resistance =================================================================
		C_UDMI(c0,t0,R) = 1.-(1.-C_UDMI(c0,t0,sigma))/(1.-C_UDMI(c0,t0,sigma) * exp(-1.0 * (1.-C_UDMI(c0,t0,sigma))*C_UDMI(c0,t0,J_v)/Ps));
		R_prime = 1.-C_UDMI(c0,t0,R);
		C_UDMI(c0,t0,m_Ap) = R_prime * m_Aw;*/
		// ======================================================================== Calculating the solute concentration at the wall ===========================================
		F_PROFILE(facet, thread_pointer, position_index) = m_Aw;
	}
	end_f_loop(facet,  thread_pointer)
	Message("Flux: %f \n", nothing);
}


// Using DEFINE_PROFILE to define the flux at the bottom wall
DEFINE_PROFILE(Nacl_yvel_bottom, thread_pointer, position_index)
{
	face_t facet;
	cell_t c0;

	Thread *t0 = thread_pointer->t0;

	begin_f_loop(facet, thread_pointer)
	{
		c0 = F_C0(facet, thread_pointer);
		F_PROFILE(facet, thread_pointer, position_index) = -1.0 * C_UDMI(c0, t0, J_v);
	}
	end_f_loop(facet, thread_pointer)
}

