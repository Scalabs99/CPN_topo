#ifndef CPN_MEAS_C
#define CPN_MEAS_C

#include "../include/cpn_conf.h"
#include <time.h>

// measure local observables on the periodic copy
void perform_measures_localobs(CPN_Conf *conf, CPN_Conf *flow_temp, Geometry const *const geo,
							   CPN_Param const *const param, FILE *datafilep, FILE *topofilep, FILE *topogradfilep, CPN_Conf *aux_conf)
{
	int cool_step = 0, i, grad_step = 0;
	double magn_susc[2], Q[3], chi_p[3];
	double energy;
	cmplx Pol_loop; // Polyakov loop;

	energy = energy_density(conf, geo, param);
	magnetic_susceptibility(conf, param, magn_susc);
	Pol_loop = compute_Polyakov(conf, geo, param);

	fprintf(datafilep, "%ld %.16lf %.16lf %.16lf %.16lf %.16lf\n", conf->update_index, energy, magn_susc[0], magn_susc[1], creal(Pol_loop), cimag(Pol_loop));
	fflush(datafilep);

	// compute topological observables of hot configuration
	for (i = 0; i < 3; i++) // 1=geometric charge U, 2=geometric charge z, 3=non-geometric plaquette charge
	{
		Q[i] = topo_charge(conf, geo, param, i);   // compute topological charge using i^th discretization
		chi_p[i] = chi_prime(conf, geo, param, i); // compute chi' using i^th discretization
	}
	// print topological observable of hot configuration (and energy)
	fprintf(topofilep, "%ld %d", conf->update_index, cool_step);
	for (i = 0; i < 3; i++)
		fprintf(topofilep, " %.16lf", Q[i]);
	for (i = 0; i < 3; i++)
		fprintf(topofilep, " %.16lf", chi_p[i]);
	fprintf(topofilep, " %.16lf", energy);
	fprintf(topofilep, "\n");

	// print topological observable of hot configuration (and energy) in the topograd file
	fprintf(topogradfilep, "%ld %d", conf->update_index, grad_step);
	for (i = 0; i < 3; i++)
		fprintf(topogradfilep, " %.16lf", Q[i]);
	for (i = 0; i < 3; i++)
		fprintf(topogradfilep, " %.16lf", chi_p[i]);
	fprintf(topogradfilep, " %.16lf", energy);
	fprintf(topogradfilep, "\n");

	// refresh stored topo charge of periodic configuration (used only for multicanonic)
	conf->stored_topo_charge = Q[0];

	// aux_conf = conf (used for cooling)
	copyconf(conf, param, aux_conf);

	// perform cooling on aux_conf, using the improved cooling function
	for (cool_step = 1; cool_step < (param->d_coolsteps + 1); cool_step++)
	{
		cooling_improved(aux_conf, geo, param);	   // perform 1 cooling step
		if ((cool_step % param->d_coolevery) == 0) // perform measures on cooled conf
		{
			// compute topological observables of cooled configuration
			for (i = 0; i < 3; i++)
			{
				Q[i] = topo_charge(aux_conf, geo, param, i);
				chi_p[i] = chi_prime(aux_conf, geo, param, i);
			}
			energy = energy_density(aux_conf, geo, param);
			// print topological observable of cooled configuration
			fprintf(topofilep, "%ld %d", conf->update_index, cool_step);
			for (i = 0; i < 3; i++)
				fprintf(topofilep, " %.16lf", Q[i]);
			for (i = 0; i < 3; i++)
				fprintf(topofilep, " %.16lf", chi_p[i]);
			fprintf(topofilep, " %.16lf", energy);
			fprintf(topofilep, "\n");
		}
	}

	// aux_conf = conf (gradient flow)
	copyconf(conf, param, aux_conf);

	// perform the gradient flow on aux_conf, using the gradient flow function
	for (grad_step = 1; grad_step < (param->d_grad_steps + 1); grad_step++)
	{
		gradient_flow(aux_conf, flow_temp, geo, param); // perform 1 cooling step
		if ((grad_step % param->d_gradevery) == 0)		// perform measures on cooled conf
		{
			// compute topological observables of cooled configuration
			for (i = 0; i < 3; i++)
			{
				Q[i] = topo_charge(aux_conf, geo, param, i);
				chi_p[i] = chi_prime(aux_conf, geo, param, i);
			}
			energy = energy_density(aux_conf, geo, param);
			// print topological observable of cooled configuration
			fprintf(topogradfilep, "%ld %d", conf->update_index, grad_step);
			for (i = 0; i < 3; i++)
				fprintf(topogradfilep, " %.16lf", Q[i]);
			for (i = 0; i < 3; i++)
				fprintf(topogradfilep, " %.16lf", chi_p[i]);
			fprintf(topogradfilep, " %.16lf", energy);
			fprintf(topogradfilep, "\n");
		}
	}

	fflush(topofilep);
	fflush(topogradfilep);
}

void perform_measure_gradient_flow(CPN_Conf const *const conf, Geometry const *const geo,
								   CPN_Param const *const param, FILE *gradfilep, FILE *argPfilep, CPN_Conf *flow_temp, CPN_Conf *aux_conf)
{
	int i, more_steps = 3e4;
	double energy, energy_out, ftheta_mean, fz_mean; // energy_in;
	double Q[3];
	long Lx = param->d_size[1];
	long j;
	double arg_P;

	// open the gradient flow force file
	FILE *f_force_grad = fopen("forces_grad.dat", "w");
	if (f_force_grad != NULL)
	{
		fprintf(f_force_grad, "# |F_z_tg|^2 \t |F_theta|^2\n");
		fflush(f_force_grad);
	}

	// aux_conf = conf ( we work on the aux conf and not on the conf )
	copyconf(conf, param, aux_conf);

	// compute and print the energy and the topological charge of the conf before the gradient flow
	energy = energy_density(aux_conf, geo, param);

	for (i = 0; i < 3; i++)
	{
		Q[i] = topo_charge(aux_conf, geo, param, i);
	}

	fprintf(gradfilep, "%.16lf", energy);
	for (i = 0; i < 3; i++)
		fprintf(gradfilep, " %.16lf", Q[i]);
	fprintf(gradfilep, "\n");

	fflush(gradfilep);

	// Compute the value of fz_mean and fu_mean on the starting configuration
	fz_mean = mean_force_z_tang(aux_conf, param, geo);
	ftheta_mean = mean_force_theta(aux_conf, param, geo);

	// print these values on the forces file
	if (f_force_grad != NULL)
	{
		fprintf(f_force_grad, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
		fflush(f_force_grad);
	}
	// Initialize energy_out with the value of energy
	energy_out = energy;

	// perform gradient flow
	do
	{

		// compute the energy before the integration step
		// energy_in = energy_out;

		// perform the integration step
		gradient_flow_tg(aux_conf, flow_temp, geo, param);

		// compute the energy after the integration step
		energy_out = energy_density(aux_conf, geo, param);

		// compute the topological charge of the configuration after the integration step
		for (i = 0; i < 3; i++)
		{
			Q[i] = topo_charge(aux_conf, geo, param, i);
		}

		// print the energy and the topological charge of the out configuration
		fprintf(gradfilep, "%.16lf", energy_out);
		for (i = 0; i < 3; i++)
			fprintf(gradfilep, " %.16lf", Q[i]);
		fprintf(gradfilep, "\n");
		fflush(gradfilep);

		// Compute the lattice mean of the forces
		fz_mean = mean_force_z_tang(aux_conf, param, geo);
		ftheta_mean = mean_force_theta(aux_conf, param, geo);

		// Print them on the file
		if (f_force_grad != NULL)
		{
			fprintf(f_force_grad, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
			fflush(f_force_grad);
		}

	} while (max(fz_mean, ftheta_mean) > 1e-9); // max(fz_mean, ftheta_mean) > 1e-9 , fabs(energy_out - energy_in) > (param->d_tolerance)

	for (j = 0; j < more_steps; j++)
	{
		// compute the energy before the integration step
		// energy_in = energy_out;

		// perform the integration step
		gradient_flow_tg(aux_conf, flow_temp, geo, param);

		// compute the energy after the integration step
		energy_out = energy_density(aux_conf, geo, param);

		// compute the topological charge of the configuration after the integration step
		for (i = 0; i < 3; i++)
		{
			Q[i] = topo_charge(aux_conf, geo, param, i);
		}

		// print the energy and the topological charge of the out configuration
		fprintf(gradfilep, "%.16lf", energy_out);
		for (i = 0; i < 3; i++)
			fprintf(gradfilep, " %.16lf", Q[i]);
		fprintf(gradfilep, "\n");
		fflush(gradfilep);

		// Compute the lattice mean of the forces
		fz_mean = mean_force_z_tang(aux_conf, param, geo);
		ftheta_mean = mean_force_theta(aux_conf, param, geo);

		// Print them on the file
		if (f_force_grad != NULL)
		{
			fprintf(f_force_grad, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
			fflush(f_force_grad);
		}
	}

	// compute the argP(n_x) on the final configuration
	for (j = 0; j < Lx; j++)
	{
		arg_P = compute_arg_Pol(aux_conf, geo, param, j);
		fprintf(argPfilep, "%ld %.16lf\n", j, arg_P);
		fflush(argPfilep);
	}

	// FIX: Chiudi il file delle forze del Gradient Flow
	if (f_force_grad != NULL)
	{
		fclose(f_force_grad);
	}
}

void measure_RK23(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param, CPN_Conf *flow1, CPN_Conf *flow2, CPN_Conf *flow3, CPN_Conf *aux_conf)
{
	int more_steps = 3e4, result;
	double energy, energy_out, ftheta_mean, fz_mean, step; // energy_in;
	long j;

	// initialize the step variable 
	step = param->d_int_step; 

	// open the gradient flow force file
	FILE *f_force_grad23 = fopen("forces_grad23.dat", "w");
	FILE *f_ener_grad23 = fopen("energy_grad23.dat", "w");
	FILE *f_step_grad23 = fopen("step_grad23.dat", "w");
	if (f_force_grad23 != NULL)
	{
		fprintf(f_force_grad23, "# |F_z_tg|^2 \t |F_theta|^2\n");
		fflush(f_force_grad23);
	}

	if (f_ener_grad23 != NULL)
	{
		fprintf(f_ener_grad23, "# Energy\n");
		fflush(f_ener_grad23);
	}

	if (f_step_grad23 != NULL)
	{
		fprintf(f_step_grad23, "# Step\n");
		fflush(f_step_grad23);
	}

	// aux_conf = conf ( we work on the aux conf and not on the conf )
	copyconf(conf, param, aux_conf);

	// compute and print the energy and the topological charge of the conf before the gradient flow
	energy = energy_density(aux_conf, geo, param);

	// print this value on the energy file
	if (f_ener_grad23 != NULL)
	{
		fprintf(f_ener_grad23, " %.16le\n", energy);
		fflush(f_ener_grad23);
	}

	// Compute the value of fz_mean and fu_mean on the starting configuration
	fz_mean = mean_force_z_tang(aux_conf, param, geo);
	ftheta_mean = mean_force_theta(aux_conf, param, geo);

	// print these values on the forces file
	if (f_force_grad23 != NULL)
	{
		fprintf(f_force_grad23, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
		fflush(f_force_grad23);
	}

	if (f_step_grad23 != NULL)
	{
		fprintf(f_step_grad23, "%.16le\n", step);
		fflush(f_step_grad23);
	}
	// Initialize energy_out with the value of energy
	energy_out = energy;

	// perform gradient flow
	do
	{

		// compute the energy before the integration step
		// energy_in = energy_out;

		// perform the integration step
		do
		{
			result = adaptive_step_RK23(aux_conf, flow1, flow2, flow3, geo, param, step);

		} while (result == 0);

		// compute the energy after the integration step
		energy_out = energy_density(aux_conf, geo, param);

		// Compute the lattice mean of the forces
		fz_mean = mean_force_z_tang(aux_conf, param, geo);
		ftheta_mean = mean_force_theta(aux_conf, param, geo);

		// Print them on the file
		if (f_force_grad23 != NULL)
		{
			fprintf(f_force_grad23, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
			fflush(f_force_grad23);
		}

		// print this value on the energy file
		if (f_ener_grad23 != NULL)
		{
			fprintf(f_ener_grad23, " %.16le\n", energy_out);
			fflush(f_ener_grad23);
		}

		// print the modified integration step on file
		if (f_step_grad23 != NULL)
		{
			fprintf(f_step_grad23, " %.16le\n", step);
			fflush(f_step_grad23);
		}

	} while (max(fz_mean, ftheta_mean) > 1e-9); // max(fz_mean, ftheta_mean) > 1e-9 , fabs(energy_out - energy_in) > (param->d_tollerance)

	for (j = 0; j < more_steps; j++)
	{
		// compute the energy before the integration step
		// energy_in = energy_out;

		do
		{
			// perform the integration step
			result = adaptive_step_RK23(aux_conf, flow1, flow2, flow3, geo, param, step);
		} while (result == 0);

		// compute the energy after the integration step
		energy_out = energy_density(aux_conf, geo, param);

		// print this value on the energy file
		if (f_ener_grad23 != NULL)
		{
			fprintf(f_ener_grad23, " %.16le\n", energy_out);
			fflush(f_ener_grad23);
		}

		// Compute the lattice mean of the forces
		fz_mean = mean_force_z_tang(aux_conf, param, geo);
		ftheta_mean = mean_force_theta(aux_conf, param, geo);

		// Print them on the file
		if (f_force_grad23 != NULL)
		{
			fprintf(f_force_grad23, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
			fflush(f_force_grad23);
		}

		// print the modified integration step on file
		if (f_step_grad23 != NULL)
		{
			fprintf(f_step_grad23, " %.16le\n", step);
			fflush(f_step_grad23);
		}


	}

	// FIX: Chiudi il file delle forze, dell'energia e dello step del Gradient Flow
	if (f_force_grad23 != NULL)
	{
		fclose(f_force_grad23);
	}
	if (f_ener_grad23 != NULL)
	{
		fclose(f_ener_grad23); 
	}
	if (f_step_grad23 != NULL)
	{
		fclose(f_step_grad23); 
	}
}

void perform_measure_cooling(CPN_Conf const *const conf, Geometry const *const geo,
							 CPN_Param const *const param, FILE *coolfilep, FILE *argPfilep, CPN_Conf *aux_conf)
{
	int i;
	double energy, energy_in, energy_out, fz_mean, ftheta_mean;
	double Q[3];
	long Lx = param->d_size[1];
	long j;
	double arg_P;

	// Open the cooling force file
	FILE *f_force_cool = fopen("forces_cool.dat", "w");
	if (f_force_cool != NULL)
	{
		fprintf(f_force_cool, "# |F_z_tang|^2 \t |F_theta|^2\n");
		fflush(f_force_cool);
	}

	// aux_conf = conf we work on the aux_conf
	copyconf(conf, param, aux_conf);

	// compute the energy
	energy = energy_density(aux_conf, geo, param);

	for (i = 0; i < 3; i++)
	{
		Q[i] = topo_charge(aux_conf, geo, param, i);
	}

	fprintf(coolfilep, "%.16lf", energy);
	for (i = 0; i < 3; i++)
		fprintf(coolfilep, " %.16lf", Q[i]);
	fprintf(coolfilep, "\n");

	// Forces the writing of the first value on the file
	fflush(coolfilep);

	// Compute the value of fz_mean and fu_mean on the starting configuration
	fz_mean = mean_force_z_tang(aux_conf, param, geo);
	ftheta_mean = mean_force_theta(aux_conf, param, geo);

	// Print the values on the file
	if (f_force_cool != NULL)
	{
		fprintf(f_force_cool, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
		fflush(f_force_cool);
	}

	// Initialize energy_out with the value of energy
	energy_out = energy;

	do
	{

		// compute the energy before the integration step
		energy_in = energy_out;

		// perform the cooling step
		cooling_improved(aux_conf, geo, param);

		// compute the energy after the integration step
		energy_out = energy_density(aux_conf, geo, param);

		// compute the topological charge of the configuration after the cooling step
		for (i = 0; i < 3; i++)
		{
			Q[i] = topo_charge(aux_conf, geo, param, i);
		}

		// print the energy and the topological charge of the out configuration
		fprintf(coolfilep, "%.16lf", energy_out);
		for (i = 0; i < 3; i++)
			fprintf(coolfilep, " %.16lf", Q[i]);
		fprintf(coolfilep, "\n");
		fflush(coolfilep);

		// Compute the lattice mean of the forces
		fz_mean = mean_force_z_tang(aux_conf, param, geo);
		ftheta_mean = mean_force_theta(aux_conf, param, geo);

		// Print the values on the file
		if (f_force_cool != NULL)
		{
			fprintf(f_force_cool, "%.16le \t %.16le\n", fz_mean, ftheta_mean);
			fflush(f_force_cool);
		}

	} while (fabs(energy_out - energy_in) > (param->d_tolerance));

	// compute the argP(n_x) on the final configuration
	for (j = 0; j < Lx; j++)
	{
		arg_P = compute_arg_Pol(aux_conf, geo, param, j);
		fprintf(argPfilep, "%ld %.16lf\n", j, arg_P);
		fflush(argPfilep);
	}

	// FIX: Chiudi il file delle forze del Gradient Flow
	if (f_force_cool != NULL)
	{
		fclose(f_force_cool);
	}
}

// compute plaquette Pi_{mu nu}(i) on site i and plane (mu,nu)
// Pi_{mu nu}(i) = U(i)_mu U(i+mu)_nu conj( U(i+nu) )_mu conj( U(i) )_nu
cmplx plaquette(CPN_Conf const *const conf, Geometry const *const geo, long const i, int const mu, int const nu)
{
	return (conf->U[i][mu] * conf->U[geo->up[i][mu]][nu] * conj(conf->U[geo->up[i][nu]][mu]) * conj(conf->U[i][nu]));
}

// compute the energy density for single link E = S_{Symanzik}(theta=0) / (2 V N beta) with TBCs
double energy_density(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param)
{
	int i, mu;
	double e = 0.0;
	cmplx e1 = 0.0 + I * 0.0, e2 = 0.0 + I * 0.0;
	cmplx aux1, aux2;

	for (i = 0; i < param->d_volume; i++)
	{
		for (mu = 0; mu < 2; mu++)
		{
			aux1 = vector_scalar_product_matrix(conf->z[geo->up[i][mu]], conf->z[i], conf->M1[i][mu]); // conj( z(i+mu) )* conj( M1(i)_mu ) * z(i)
			e1 += conj(conf->U[i][mu]) * aux1;

			aux2 = vector_scalar_product_matrix(conf->z[geo->up[geo->up[i][mu]][mu]], conf->z[i], conf->M2[i][mu]); // conj( z(i+2mu) ) * conj( M2(i)_mu ) * z(i)
			e2 += conj(conf->U[geo->up[i][mu]][mu] * conf->U[i][mu]) * aux2;
		}
	}

	e = 2.0 * (c1 * creal(e1) + c2 * creal(e2));
	e /= (double)(param->d_volume);
	e = 2.0 * (c1 + c2) - e / 2.0;
	return e;
}

// compute the geometric topological charge density expressed in terms of the scalar field z on site i
double geo_topo_charge_z_density(CPN_Conf const *const conf, Geometry const *const geo, long const i)
{
	cmplx aux_1, aux_2, aux_3, p1, p2;
	double q;
	int mu = 0;
	int nu = 1 - mu;

	aux_1 = vector_scalar_product_matrix(conf->z[geo->up[geo->up[i][mu]][nu]], conf->z[i], conf->M1[i][mu]);
	aux_2 = vector_scalar_product_two_matrices(conf->z[geo->up[i][mu]], conf->z[geo->up[geo->up[i][mu]][nu]], conf->M1[i][mu], conf->M1[i][mu]);
	aux_3 = vector_scalar_product_matrix_inverted(conf->z[i], conf->z[geo->up[i][mu]], conf->M1[i][mu]);
	p1 = aux_1 * aux_2 * aux_3;

	aux_1 = vector_scalar_product(conf->z[geo->up[i][nu]], conf->z[i]);
	aux_2 = vector_scalar_product_matrix(conf->z[geo->up[geo->up[i][mu]][nu]], conf->z[geo->up[i][nu]], conf->M1[i][mu]);
	aux_3 = vector_scalar_product_matrix_inverted(conf->z[i], conf->z[geo->up[geo->up[i][mu]][nu]], conf->M1[i][mu]);
	p2 = aux_1 * aux_2 * aux_3;

	q = -1.0 * (arg(p1) + arg(p2)) / (2.0 * pi);
	return q;
}

// compute the geometric topological charge expressed in terms of the gauge field U on site i
double geo_topo_charge_U_density(CPN_Conf const *const conf, Geometry const *const geo, long const i)
{
	int mu = 0;
	int nu = 1 - mu;
	cmplx plaq;
	double q;

	plaq = plaquette(conf, geo, i, mu, nu); // Plaq = plaquette(i)_{01}
	q = arg(plaq) / (2.0 * pi);				// {Im log(plaq) } / (2 pi)
	return q;
}

// compute the non-geometric topological charge expressed in terms of the plaquettes on site i
double plaq_topo_charge_density(CPN_Conf const *const conf, Geometry const *const geo, long const i)
{
	int mu = 0;
	int nu = 1 - mu;
	double q;

	q = cimag(plaquette(conf, geo, i, mu, nu)); // Im{ plaquette(i)_{01} }
	q /= (2.0 * pi);
	return q;
}

// compute the total topological charge
double topo_charge(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param, int const which_charge)
{
	// pointer to desired discretization of the topological charge density
	double (*q_ptr)(CPN_Conf const *const, Geometry const *const, long const);
	if (which_charge == 0)
	{
		q_ptr = &geo_topo_charge_U_density;
	}
	if (which_charge == 1)
	{
		q_ptr = &geo_topo_charge_z_density;
	}
	if (which_charge == 2)
	{
		q_ptr = &plaq_topo_charge_density;
	}

	double Q = 0.0;
	long i;
	for (i = 0; i < param->d_volume; i++)
	{
		Q += ((*q_ptr)(conf, geo, i));
	} // Q = sum_i { q(i) }
	return Q;
}

// compute quantity G = (1/4) sum_{x} |x|^2 q(x)q(0) needed to compute chi'
// the mean value of this sum in the continuum limit is chi' = (1/4) \int d^2x |x|^2 <q(x)q(0)> and is related to the first derivative with respect to k=p^2 of
// the Fourier transform FG of the topological charge density correlator: chi' = - (1/4) lim_{k->0} d FG(k) / dk
// where FG(p^2) = \int d^2 x e^(ipx) <q(x)q(0)>
double chi_prime(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param, int const which_charge)
{
	// pointer to desired discretization of the topological charge density
	double (*q_ptr)(CPN_Conf const *const, Geometry const *const, long const);
	if (which_charge == 0)
	{
		q_ptr = &geo_topo_charge_U_density;
	}
	if (which_charge == 1)
	{
		q_ptr = &geo_topo_charge_z_density;
	}
	if (which_charge == 2)
	{
		q_ptr = &plaq_topo_charge_density;
	}

	double d2, G = 0.0;
	long i;
	for (i = 0; i < param->d_volume; i++)
	{
		d2 = square_distance(i, 0, param);	// i->x ==> d2 = d(x,0)^2 = |x|^2 = |i|^2 where distance is computed on the torus
		G += ((*q_ptr)(conf, geo, i)) * d2; // sum_i { q(i) |i|^2 }
	}
	G *= ((*q_ptr)(conf, geo, 0)); // sum_i { q(i)q(0) |i|^2 }
	G /= 4.0;					   // G = (1/4) sum_i { q(i)q(0) |i|^2 }
	return G;
}

// compute quantities needed for the "magnetic susceptibility" at momentum p=0 (stored in magn_susc[0]) and at momentum p = q = 2pi/L (stored in magn_susc[1])
// q is the smallest momentum possible on a lattice with L sites, where L is the smallest size of the lattice
// the magnetic susceptibility can be expressed in terms of G(p) = Fourier transform of G(x) = P_(ij)(x) P_(ij)(0) - 1/N, where P_(ij)(x) = z*_i(x) z_j(x)
// magnetic susceptibility at p=0: chi_m(0) = <G(p=0)> (this quantity is trivially real)
// magnetic susceptibility at p=q: chi_m(q) = <G(p=q)> (here I take just the real part, as the imaginary part averages to zero over the ensamble)
// mang_susc[0] = G(p=0), magn_susc[1]=G(p=q), the correlation length of the system can be expressed in terms of chi_m(q)/chi_m(0)
void magnetic_susceptibility(CPN_Conf const *const conf, CPN_Param const *const param, double *magn_susc)
{
	int i;
	double L, sum1 = 0.0, q;
	cmplx a;
	cmplx sum2 = 0.0 + I * 0.0;
	long x[2];

	// L = min(Lx,Lt)
	if (param->d_size[0] > param->d_size[1])
		L = ((double)param->d_size[1]);
	else
		L = ((double)param->d_size[0]);
	q = 2.0 * pi / L;

	for (i = 0; i < param->d_volume; i++)
	{
		a = vector_scalar_product(conf->z[i], conf->z[0]);
		sum1 += cmplx_norm(a);
		si_to_cart(x, i, param); // i->x
		sum2 += (cmplx_norm(a) - 1.0 / ((double)N)) * cexp(I * (q * ((double)x[0])));
	}
	sum1 -= ((double)(param->d_volume)) / ((double)N);

	magn_susc[0] = sum1;
	magn_susc[1] = creal(sum2);
}

// perform a single cooling step on the given conf (this function should receive an aux conf)
// cooling = each field is aligned to its local force: z = F_z / |F_z|, U = F_U / |F_U|
// forces are determined using the non-improved lattice action S_L \propto sum_{i, mu} U(i)_mu conj(z)(i+mu) z(i)
void cooling(CPN_Conf *conf, Geometry const *const geo, CPN_Param const *const param)
{
	int i, mu;
	cmplx F_U;
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx aux[N] __attribute__((aligned(DOUBLE_ALIGN)));

	// cool U fields
	for (i = 0; i < param->d_volume; i++)
	{
		for (mu = 0; mu < 2; mu++)
		{
			F_U = vector_scalar_product(conf->z[geo->up[i][mu]], conf->z[i]); // F_U = conj(z(i+mu)) z(i)
			conf->U[i][mu] = F_U / cmplx_abs(F_U);							  // align U along the local force on link (i,mu): U = F_U/|F_U|
		}
	}

	// cool z fields
	for (i = 0; i < param->d_volume; i++)
	{
		vector_zero(F_z); // F_z = 0
		for (mu = 0; mu < 2; mu++)
		{
			// aux = U(i)_mu * z(i+mu) + conj( U(i-mu)_mu ) * z(i-mu)
			vector_linear_combination_cmplx_coeff(aux, conf->z[geo->up[i][mu]], conf->z[geo->dn[i][mu]], conf->U[i][mu], conj(conf->U[geo->dn[i][mu]][mu]));
			vector_sum(F_z, aux); // F_z += aux;
		}
		vector_normalization(F_z);	   // F_z -> F_z/|F_z|
		vector_equal(conf->z[i], F_z); // align z along the local force on site i: z = F_z/|F_z|
	}
}

// perform a single cooling step on the given conf (this function should receive an aux conf)
// cooling = each field is aligned to its local force: z = F_z / |F_z|, U = F_U / |F_U|
// forces are determined using the improved lattice action ( Symanzik ) with TBC
void cooling_improved(CPN_Conf *conf, Geometry const *const geo, CPN_Param const *const param)
{
	int i, mu;
	cmplx F_U;
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));

	// cool U fields
	for (i = 0; i < param->d_volume; i++)
	{
		for (mu = 0; mu < 2; mu++)
		{
			F_U = force_U(conf, geo, param, i, mu); // compute the force using the improved formula;
			conf->U[i][mu] = F_U / cmplx_abs(F_U);	// align U along the local force on link (i,mu): U = F_U/|F_U|
		}
	}

	// cool z fields
	for (i = 0; i < param->d_volume; i++)
	{
		force_z(conf, geo, i, F_z);	   // compute the force using the improved formula;
		vector_normalization(F_z);	   // F_z -> F_z/|F_z|
		vector_equal(conf->z[i], F_z); // align z along the local force on site i: z = F_z/|F_z|
	}

	fix_gauge_conf(conf, param, geo);
}

// Perform a single integration step for the gradient flow equations
// The integration scheme is the simple Euler scheme
void gradient_flow(CPN_Conf *conf, CPN_Conf *flow_temp, Geometry const *const geo, CPN_Param const *const param)
{
	int i, mu;
	double c = 2.0 * (param->d_beta) * N * (param->d_int_step);
	cmplx F_U;
	cmplx F_z[N] __attribute__((aligned(DOUBLE_ALIGN)));

	for (i = 0; i < param->d_volume; i++)
	{
		// Compute the force F_U and the Euler evolution step
		for (mu = 0; mu < 2; mu++)
		{
			F_U = force_U(conf, geo, param, i, mu);
			flow_temp->U[i][mu] = conf->U[i][mu] + c * F_U; // save the updated link variable U
			flow_temp->U[i][mu] = flow_temp->U[i][mu] / cmplx_abs(flow_temp->U[i][mu]);
		}

		// Compute the force F_z and the Euler evolution step
		force_z(conf, geo, i, F_z);
		vector_times_real_const(F_z, c);

		// Assign the value of conf->z[i] to flow_temp->z[i]
		vector_equal(flow_temp->z[i], conf->z[i]);

		// Euler step
		vector_sum(flow_temp->z[i], F_z);

		// normalize z[i] 
		vector_normalization(flow_temp->z[i]);
	}

	for (i = 0; i < param->d_volume; i++)
	{
		vector_equal(conf->z[i], flow_temp->z[i]);
		conf->U[i][0] = flow_temp->U[i][0];
		conf->U[i][1] = flow_temp->U[i][1];
	}
}

// Perform a single integration step for the gradient flow equations
// The integration scheme is the simple Euler scheme
void gradient_flow_tg(CPN_Conf *conf, CPN_Conf *flow_temp, Geometry const *const geo, CPN_Param const *const param)
{
	int i, mu;
	double c_theta = 2.0 * (param->d_beta) * N * (param->d_int_step);
	double c_z = 2.0 * (param->d_beta) * N * (param->d_int_step);
	double f_theta;
	cmplx F_z_tg[N] __attribute__((aligned(DOUBLE_ALIGN)));

	for (i = 0; i < param->d_volume; i++)
	{

		// Compute the force F_theta and the Euler evolution step
		for (mu = 0; mu < 2; mu++)
		{
			f_theta = F_theta(conf, geo, i, mu);
			flow_temp->U[i][mu] = conf->U[i][mu] * cexp(I * c_theta * f_theta); // use the phase to compute the auxiliary conf
		}

		// Compute the force F_z and the Euler evolution step
		force_z_tangent(conf, geo, i, F_z_tg);
		vector_times_real_const(F_z_tg, c_z);

		// Assign the value of conf->z[i] to flow_temp->z[i]
		vector_equal(flow_temp->z[i], conf->z[i]);

		// Euler step
		vector_sum(flow_temp->z[i], F_z_tg);
		// normalize z[i]
		vector_normalization(flow_temp->z[i]);
	}

	for (i = 0; i < param->d_volume; i++)
	{
		vector_equal(conf->z[i], flow_temp->z[i]);
		conf->U[i][0] = flow_temp->U[i][0];
		conf->U[i][1] = flow_temp->U[i][1];
	}
}

int adaptive_step_RK23(CPN_Conf *conf, CPN_Conf *flow_temp1, CPN_Conf *flow_temp2, CPN_Conf *flow_temp3, Geometry const *const geo, CPN_Param const *const param, double step)
{
	int i, mu, j;
	double c_theta = 2.0 * (param->d_beta) * N * (step);
	double c_z = 2.0 * (param->d_beta) * N * (step);
	double error_abs_z, error_absU, error_abs_in, error_abs = 0.0;
	double errorU_0, errorU_1;
	double q, p = 2.0;
	cmplx F_z_tg1[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx F_z_tg2[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx F_z_tg3[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx error[N] __attribute__((aligned(DOUBLE_ALIGN)));
	double f_theta1[2] __attribute__((aligned(DOUBLE_ALIGN)));
	double f_theta2[2] __attribute__((aligned(DOUBLE_ALIGN)));
	double f_theta3[2] __attribute__((aligned(DOUBLE_ALIGN)));

	cmplx z_ord2[N] __attribute__((aligned(DOUBLE_ALIGN)));
	cmplx U_ord2[2] __attribute__((aligned(DOUBLE_ALIGN)));



	for (i = 0; i < param->d_volume; i++)
	{
		// Compute the force F_theta and the Euler evolution step
		for (mu = 0; mu < 2; mu++)
		{
			f_theta1[mu] = F_theta(conf, geo, i, mu);
			flow_temp1->U[i][mu] = conf->U[i][mu] * cexp(I * c_theta * (2.0 / 3) * f_theta1[mu]); // this is csi 2
		}
		// Now repeat the same steps with the z field
		// Compute csi2
		force_z_tangent(conf, geo, i, F_z_tg1);
		vector_times_real_const(F_z_tg1, (2.0 / 3) * c_z);
		// Assign the value of conf->z[i] to flow_temp->z[i]
		vector_equal(flow_temp1->z[i], conf->z[i]);
		// Euler step
		vector_sum(flow_temp1->z[i], F_z_tg1);
		// normalize z[i]
		vector_normalization(flow_temp1->z[i]);
	}

	for (i = 0; i < param->d_volume; i++)
	{
		// Compute the force F_theta and the Euler evolution step
		for (mu = 0; mu < 2; mu++)
		{
			f_theta2[mu] = F_theta(flow_temp1, geo, i, mu);
			flow_temp2->U[i][mu] = conf->U[i][mu] * cexp(I * c_theta * (2.0 / 3) * f_theta2[mu]); // this is csi 3
		}
		// Now repeat the same steps with the z field
		// Compute csi3
		force_z_tangent(flow_temp1, geo, i, F_z_tg2);
		vector_times_real_const(F_z_tg2, (2.0 / 3) * c_z);
		// Assign the value of conf->z[i] to flow_temp->z[i]
		vector_equal(flow_temp2->z[i], conf->z[i]);
		// Euler step
		vector_sum(flow_temp2->z[i], F_z_tg2);
		// normalize z[i]
		vector_normalization(flow_temp2->z[i]);
	}

	for (i = 0; i < param->d_volume; i++)
	{

		// Compute f(csi3), f(csi2), f(csi1)
		force_z_tangent(conf, geo, i, F_z_tg1);
		force_z_tangent(flow_temp1, geo, i, F_z_tg2);
		force_z_tangent(flow_temp2, geo, i, F_z_tg3);

		for (mu = 0; mu < 2; mu++)
		{
			f_theta1[mu] = F_theta(conf, geo, i, mu);
			f_theta2[mu] = F_theta(flow_temp1, geo, i, mu);
			f_theta3[mu] = F_theta(flow_temp2, geo, i, mu);

			z_ord2[mu] = conf->U[i][mu] * cexp(I * c_theta * (0.25 * f_theta1[mu] + 0.75 * f_theta2[mu]));
			flow_temp3->U[i][mu] = conf->U[i][mu] * cexp(I * c_theta * (0.25 * f_theta1[mu] + 0.375 * (f_theta2[mu] + f_theta3[mu])));
		}
		// Compute the candidate solution (two stages)
		vector_equal(z_ord2, conf->z[i]);
		vec_lin_comb_3_real_coeff(z_ord2, F_z_tg1, F_z_tg2, 1.0, 0.25 * c_z, 0.75 * c_z);
		vector_normalization(z_ord2);

		// Compute the control solution
		vector_equal(flow_temp3->z[i], conf->z[i]);
		vec_lin_comb_4_real_coeff(flow_temp3->z[i], F_z_tg1, F_z_tg2, F_z_tg3, 1.0, 0.25 * c_z, 0.375 * c_z, 0.375 * c_z);
		vector_normalization(flow_temp3->z[i]);

		// Compute the error for z;
		for (j = 0; j < N; j++)
			error[j] = flow_temp3->z[i][j] - z_ord2[j];

		error_abs_z = vector_abs(error);

		// Do the same for U[i][mu]
		errorU_0 = cmplx_abs((flow_temp3->U[i][0] - U_ord2[0]));
		errorU_1 = cmplx_abs((flow_temp3->U[i][1] - U_ord2[1]));
		error_absU = max(errorU_0, errorU_1);

		// Take the max between the two
		error_abs_in = max(error_abs_z, error_absU);
		error_abs = max(error_abs, error_abs_in);
	}

	// Compute the factor q = (safety_factor) * (epsilon/error_abs)^(1/p+1);
	double safety = 0.9;
	q = safety * pow((param->d_epsilon / (error_abs + 1e-15)), 1.0 / (p + 1.0)); // add 1e-15 to avoid dividing by zero


	// Limits to avoid a variation too big or too small of the step
	if (q > 2.0)
		q = 2.0;
	if (q < 0.2)
		q = 0.2;

	step = step * q;

	// Accept/reject step
	if (error_abs <= param->d_epsilon)
	{
		// STEP ACCEPTED: update using the better solution (flow_temp2, order 3)
		for (i = 0; i < param->d_volume; i++)
		{
			vector_equal(conf->z[i], flow_temp3->z[i]);
			conf->U[i][0] = flow_temp3->U[i][0];
			conf->U[i][1] = flow_temp3->U[i][1];
		}
		return 1; // success
	}

	// STEP REJECTED: the error is to big
	// the starting conf remains the same
	return 0;
}

// compute the arg(P(n_x)) for the cooled ( either with GF or cooling ) configuration
double compute_arg_Pol(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param, long const i)
{
	double arg_pol = 0.0;
	long j, r;
	int Lt = param->d_size[0];
	long cart_coord[2] = {0, i};
	r = cart_to_si(cart_coord, param);
	cmplx temp = 1.0;

	for (j = 0; j < Lt; j++)
	{
		temp *= conf->U[r][0];
		r = geo->up[r][0]; // jump to the next site ( in the 0 direction )
	}

	arg_pol = arg(temp);
	return arg_pol;
}

// compute the spatially averaged Polyakov loop;
cmplx compute_Polyakov(CPN_Conf const *const conf, Geometry const *const geo, CPN_Param const *const param)
{

	cmplx Pol = 0.0; // spatially averaged Polyakov loop
	int Lx = param->d_size[1];
	int Lt = param->d_size[0];

	long i, j, r;
	for (i = 0; i < Lx; i++)
	{
		cmplx temp = 1.0; // initialize temp = 1 + 0i;
		long cart_coord[2] = {0, i};
		r = cart_to_si(cart_coord, param);

		for (j = 0; j < Lt; j++)
		{
			temp *= conf->U[r][0];
			r = geo->up[r][0]; // jump to the next site ( in the 0 direction )
		}

		Pol += temp;
	}

	Pol = Pol / Lx;
	return Pol;
}

// Compute the mean over the lattice of the squared modulus of the forces F_z
double mean_force_z_tang(CPN_Conf const *const conf, CPN_Param const *const param, Geometry const *const geo)
{
	double fz_sq_mean = 0.0;
	int i;
	for (i = 0; i < param->d_volume; i++)
	{
		double fz_sq = 0.0;

		// Compute z force on site i
		cmplx f_z_array[N] __attribute__((aligned(DOUBLE_ALIGN)));
		force_z_tangent(conf, geo, i, f_z_array);

		fz_sq = vector_norm(f_z_array); // squared modulus of F_z

		fz_sq_mean += fz_sq;
	}

	fz_sq_mean = fz_sq_mean / (double)param->d_volume;

	return fz_sq_mean;
}

// Compute the mean over the lattice of the squared modulus of the forces F_U
double mean_force_U(CPN_Conf const *const conf, CPN_Param const *const param, Geometry const *const geo)
{
	double fu_sq_mean = 0.0;
	int i;
	for (i = 0; i < param->d_volume; i++)
	{

		double fu_sq = 0.0;
		int mu;

		// Compute the U force on site i, and sum over the links (mu=0 and mu=1)
		for (mu = 0; mu < 2; mu++)
		{
			cmplx fu_val = force_U(conf, geo, param, i, mu);
			fu_sq += cmplx_norm(fu_val); // sum the squared modulus over mu;
		}

		fu_sq_mean += fu_sq;
	}

	fu_sq_mean = fu_sq_mean / (double)param->d_volume;

	return fu_sq_mean;
}

// Compute the mean over the lattice of the squared modulus of the forces F_theta
double mean_force_theta(CPN_Conf const *const conf, CPN_Param const *const param, Geometry const *const geo)
{
	double ftheta_sq_mean = 0.0;
	int i;
	for (i = 0; i < param->d_volume; i++)
	{

		double ftheta_sq = 0.0;
		int mu;

		// Compute the U force on site i, and sum over the links (mu=0 and mu=1)
		for (mu = 0; mu < 2; mu++)
		{
			cmplx ftheta_val = F_theta(conf, geo, i, mu);
			ftheta_sq += cmplx_norm(ftheta_val); // sum the squared modulus over mu;
		}

		ftheta_sq_mean += ftheta_sq;
	}

	ftheta_sq_mean = ftheta_sq_mean / (double)param->d_volume;

	return ftheta_sq_mean;
}

#endif
