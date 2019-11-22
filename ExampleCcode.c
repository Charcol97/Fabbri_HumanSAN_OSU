//==============================================================================
// CellML file:   E:\humanSAN\HumanSAN_Fabbri_Fantini_Wilders_Severi_2017.cellml
// CellML model:  Human_SAN_Fabbri_Fantini_Wilders_Severi_2017
// Date and time: 04/10/2017 at 09:52:16
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C was done using COR (0.9.31.1409)
//    Copyright 2002-2017 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define X   10
#define Y   9
#define Z   4 	
#define D1S   1.2		
#define D2S   0.6		
#define DDS   1.08		
#define dx  0.03		// mm

double  Y0[X+1][Y+1][Z+1];
double  Y1[X+1][Y+1][Z+1];
double  Y2[X+1][Y+1][Z+1];
double  Y3[X+1][Y+1][Z+1];
double  Y4[X+1][Y+1][Z+1];
double  Y5[X+1][Y+1][Z+1];
double  Y6[X+1][Y+1][Z+1];
double  Y7[X+1][Y+1][Z+1];
double  Y8[X+1][Y+1][Z+1];
double  Y9[X+1][Y+1][Z+1];
double Y10[X+1][Y+1][Z+1];
double Y11[X+1][Y+1][Z+1];
double Y12[X+1][Y+1][Z+1];
double Y13[X+1][Y+1][Z+1];
double Y14[X+1][Y+1][Z+1];
double Y15[X+1][Y+1][Z+1];
double Y16[X+1][Y+1][Z+1];
double Y17[X+1][Y+1][Z+1];
double Y18[X+1][Y+1][Z+1];
double Y19[X+1][Y+1][Z+1];
double Y20[X+1][Y+1][Z+1];
double Y21[X+1][Y+1][Z+1];
double Y22[X+1][Y+1][Z+1];
double Y23[X+1][Y+1][Z+1];
double Y24[X+1][Y+1][Z+1];
double Y25[X+1][Y+1][Z+1];
double Y26[X+1][Y+1][Z+1];
double Y27[X+1][Y+1][Z+1];
double Y28[X+1][Y+1][Z+1];
double Y29[X+1][Y+1][Z+1];
double Y30[X+1][Y+1][Z+1];
double Y31[X+1][Y+1][Z+1];
double Y32[X+1][Y+1][Z+1];
double dY0;
double dY1;
double dY2;
double dY3;
double dY4;
double dY5;
double dY6;
double dY7;
double dY8;
double dY9;
double dY10;
double dY11;
double dY12;
double dY13;
double dY14;
double dY15;
double dY16;
double dY17;
double dY18;
double dY19;
double dY20;
double dY21;
double dY22;
double dY23;
double dY24;
double dY25;
double dY26;
double dY27;
double dY28;
double dY29;
double dY30;
double dY31;
double dY32;

double EC50_SR;   // millimolar (in Ca_SR_release)
double HSR;   // dimensionless (in Ca_SR_release)
double MaxSR;   // dimensionless (in Ca_SR_release)
double MinSR;   // dimensionless (in Ca_SR_release)
double kiCa;   // per_millimolar_second (in Ca_SR_release)
double kim;   // per_second (in Ca_SR_release)
double koCa;   // per_millimolar2_second (in Ca_SR_release)
double kom;   // per_second (in Ca_SR_release)
double ks;   // per_second (in Ca_SR_release)
double CM_tot;   // millimolar (in Ca_buffering)
double CQ_tot;   // millimolar (in Ca_buffering)
double Mgi;   // millimolar (in Ca_buffering)
double TC_tot;   // millimolar (in Ca_buffering)
double TMC_tot;   // millimolar (in Ca_buffering)
double kb_CM;   // per_second (in Ca_buffering)
double kb_CQ;   // per_second (in Ca_buffering)
double kb_TC;   // per_second (in Ca_buffering)
double kb_TMC;   // per_second (in Ca_buffering)
double kb_TMM;   // per_second (in Ca_buffering)
double kf_CM;   // per_millimolar_second (in Ca_buffering)
double kf_CQ;   // per_millimolar_second (in Ca_buffering)
double kf_TC;   // per_millimolar_second (in Ca_buffering)
double kf_TMC;   // per_millimolar_second (in Ca_buffering)
double kf_TMM;   // per_millimolar_second (in Ca_buffering)
double K_up;   // millimolar (in Ca_intracellular_fluxes)
double P_up_basal;   // millimolar_per_second (in Ca_intracellular_fluxes)
double slope_up;   // millimolar (in Ca_intracellular_fluxes)
double tau_dif_Ca;   // second (in Ca_intracellular_fluxes)
double tau_tr;   // second (in Ca_intracellular_fluxes)
double L_cell;   // micrometre (in Cell_parameters)
double L_sub;   // micrometre (in Cell_parameters)
double R_cell;   // micrometre (in Cell_parameters)
double V_i_part;   // dimensionless (in Cell_parameters)
double V_jsr_part;   // dimensionless (in Cell_parameters)
double V_nsr_part;   // dimensionless (in Cell_parameters)
double Cao;   // millimolar (in Ionic_values)
double Ki;   // millimolar (in Ionic_values)
double Ko;   // millimolar (in Ionic_values)
double Nao;   // millimolar (in Ionic_values)
double C;   // microF (in Membrane)
double F;   // coulomb_per_mole (in Membrane)
double R2;   // joule_per_kilomole_kelvin (R in Membrane)
double T;   // kelvin (in Membrane)
double clamp_mode;   // dimensionless (in Membrane)
double Nai_clamp;   // dimensionless (in Nai_concentration)
double ACh;   // millimolar (in Rate_modulation_experiments)
double Iso_1_uM;   // dimensionless (in Rate_modulation_experiments)
double V_holding;   // millivolt (in Voltage_clamp)
double V_test;   // millivolt (in Voltage_clamp)
double t_holding;   // second (in Voltage_clamp)
double t_test;   // second (in Voltage_clamp)
double V_dL;   // millivolt (in i_CaL_dL_gate)
double k_dL;   // millivolt (in i_CaL_dL_gate)
double Km_fCa;   // millimolar (in i_CaL_fCa_gate)
double alpha_fCa;   // per_second (in i_CaL_fCa_gate)
double k_fL;   // millivolt (in i_CaL_fL_gate)
double shift_fL;   // millivolt (in i_CaL_fL_gate)
double P_CaL;   // nanoA_per_millimolar (in i_CaL)
double offset_fT;   // second (in i_CaT_fT_gate)
double P_CaT;   // nanoA_per_millimolar (in i_CaT)
double ACh_on;   // dimensionless (in i_KACh)
double g_KACh;   // microS (in i_KACh)
double g_Kr;   // microS (in i_Kr)
double g_Ks_;   // microS (in i_Ks)
double g_Kur;   // microS (in i_Kur)
double K1ni;   // millimolar (in i_NaCa)
double K1no;   // millimolar (in i_NaCa)
double K2ni;   // millimolar (in i_NaCa)
double K2no;   // millimolar (in i_NaCa)
double K3ni;   // millimolar (in i_NaCa)
double K3no;   // millimolar (in i_NaCa)
double K_NaCa;   // nanoA (in i_NaCa)
double Kci;   // millimolar (in i_NaCa)
double Kcni;   // millimolar (in i_NaCa)
double Kco;   // millimolar (in i_NaCa)
double Qci;   // dimensionless (in i_NaCa)
double Qco;   // dimensionless (in i_NaCa)
double Qn;   // dimensionless (in i_NaCa)
double blockade_NaCa;   // dimensionless (in i_NaCa)
double Km_Kp;   // millimolar (in i_NaK)
double Km_Nap;   // millimolar (in i_NaK)
double i_NaK_max;   // nanoA (in i_NaK)
double delta_m;   // millivolt (in i_Na_m_gate)
double g_Na;   // microS (in i_Na)
double g_Na_L;   // microS (in i_Na)
double y_shift;   // millivolt (in i_f_y_gate)
double Km_f;   // millimolar (in i_f)
double alpha;   // dimensionless (in i_f)
double blockade;   // dimensionless (in i_f)
double g_f;   // microS (in i_f)
double g_to;   // microS (in i_to)
double j_SRCarel;   // millimolar_per_second (in Ca_SR_release)
double kCaSR;   // dimensionless (in Ca_SR_release)
double kiSRCa;   // per_millimolar_second (in Ca_SR_release)
double koSRCa;   // per_millimolar2_second (in Ca_SR_release)
double delta_fCMi;   // per_second (in Ca_buffering)
double delta_fCMs;   // per_second (in Ca_buffering)
double delta_fCQ;   // per_second (in Ca_buffering)
double delta_fTC;   // per_second (in Ca_buffering)
double delta_fTMC;   // per_second (in Ca_buffering)
double delta_fTMM;   // per_second (in Ca_buffering)
double P_up;   // millimolar_per_second (in Ca_intracellular_fluxes)
double b_up;   // dimensionless (in Ca_intracellular_fluxes)
double j_Ca_dif;   // millimolar_per_second (in Ca_intracellular_fluxes)
double j_tr;   // millimolar_per_second (in Ca_intracellular_fluxes)
double j_up;   // millimolar_per_second (in Ca_intracellular_fluxes)
double V_cell;   // millimetre3 (in Cell_parameters)
double V_i;   // millimetre3 (in Cell_parameters)
double V_jsr;   // millimetre3 (in Cell_parameters)
double V_nsr;   // millimetre3 (in Cell_parameters)
double V_sub;   // millimetre3 (in Cell_parameters)
double E_Ca;   // millivolt (in Ionic_values)
double E_K;   // millivolt (in Ionic_values)
double E_Na;   // millivolt (in Ionic_values)
double RTonF;   // millivolt (in Membrane)
double V;   // millivolt (in Membrane)
double i_tot;   // nanoA (in Membrane)
double Nai;   // millimolar (in Nai_concentration)
double V_clamp;   // millivolt (in Voltage_clamp)
double Iso_shift_dL;   // millivolt (in i_CaL_dL_gate)
double Iso_slope_dL;   // dimensionless (in i_CaL_dL_gate)
double adVm;   // millivolt (in i_CaL_dL_gate)
double alpha_dL;   // per_second (in i_CaL_dL_gate)
double bdVm;   // millivolt (in i_CaL_dL_gate)
double beta_dL;   // per_second (in i_CaL_dL_gate)
double dL_infinity;   // dimensionless (in i_CaL_dL_gate)
double tau_dL;   // second (in i_CaL_dL_gate)
double fCa_infinity;   // dimensionless (in i_CaL_fCa_gate)
double tau_fCa;   // second (in i_CaL_fCa_gate)
double fL_infinity;   // dimensionless (in i_CaL_fL_gate)
double tau_fL;   // second (in i_CaL_fL_gate)
double ACh_block;   // dimensionless (in i_CaL)
double Iso_increase_1;   // dimensionless (Iso_increase in i_CaL)
double i_CaL;   // nanoA (in i_CaL)
double i_siCa;   // nanoA (in i_CaL)
double i_siK;   // nanoA (in i_CaL)
double i_siNa;   // nanoA (in i_CaL)
double dT_infinity;   // dimensionless (in i_CaT_dT_gate)
double tau_dT;   // second (in i_CaT_dT_gate)
double fT_infinity;   // dimensionless (in i_CaT_fT_gate)
double tau_fT;   // second (in i_CaT_fT_gate)
double i_CaT;   // nanoA (in i_CaT)
double b_infinity;   // dimensionless (in i_KACh_a_gate)
double alpha_b;   // per_second (in i_KACh_a_gate)
double beta_b;   // per_second (in i_KACh_a_gate)
double tau_b;   // second (in i_KACh_a_gate)
double i_KACh;   // nanoA (in i_KACh)
double alfapaF;   // per_second (in i_Kr_pa_gate)
double betapaF;   // per_second (in i_Kr_pa_gate)
double pa_infinity;   // dimensionless (in i_Kr_pa_gate)
double tau_paF;   // second (in i_Kr_pa_gate)
double tau_paS;   // second (in i_Kr_pa_gate)
double pi_infinity;   // dimensionless (in i_Kr_pi_gate)
double tau_pi;   // second (in i_Kr_pi_gate)
double i_Kr;   // nanoA (in i_Kr)
double Iso_shift_1;   // millivolt (Iso_shift in i_Ks_n_gate)
double alpha_n;   // per_second (in i_Ks_n_gate)
double beta_n;   // per_second (in i_Ks_n_gate)
double n_infinity;   // dimensionless (in i_Ks_n_gate)
double tau_n;   // second (in i_Ks_n_gate)
double E_Ks;   // millivolt (in i_Ks)
double g_Ks;   // microS (in i_Ks)
double i_Ks;   // nanoA (in i_Ks)
double r_Kur_infinity;   // dimensionless (in i_Kur_rKur_gate)
double tau_r_Kur;   // second (in i_Kur_rKur_gate)
double s_Kur_infinity;   // dimensionless (in i_Kur_sKur_gate)
double tau_s_Kur;   // second (in i_Kur_sKur_gate)
double i_Kur;   // nanoA (in i_Kur)
double di;   // dimensionless (in i_NaCa)
double dodo;   // dimensionless (in i_NaCa)
double i_NaCa;   // nanoA (in i_NaCa)
double k12;   // dimensionless (in i_NaCa)
double k14;   // dimensionless (in i_NaCa)
double k21;   // dimensionless (in i_NaCa)
double k23;   // dimensionless (in i_NaCa)
double k32;   // dimensionless (in i_NaCa)
double k34;   // dimensionless (in i_NaCa)
double k41;   // dimensionless (in i_NaCa)
double k43;   // dimensionless (in i_NaCa)
double x1;   // dimensionless (in i_NaCa)
double x2;   // dimensionless (in i_NaCa)
double x3;   // dimensionless (in i_NaCa)
double x4;   // dimensionless (in i_NaCa)
double Iso_increase_2;   // dimensionless (Iso_increase in i_NaK)
double i_NaK;   // nanoA (in i_NaK)
double alpha_h;   // per_second (in i_Na_h_gate)
double beta_h;   // per_second (in i_Na_h_gate)
double h_infinity;   // dimensionless (in i_Na_h_gate)
double tau_h;   // second (in i_Na_h_gate)
double E0_m;   // millivolt (in i_Na_m_gate)
double alpha_m;   // per_second (in i_Na_m_gate)
double beta_m;   // per_second (in i_Na_m_gate)
double m_infinity;   // dimensionless (in i_Na_m_gate)
double tau_m;   // second (in i_Na_m_gate)
double E_mh;   // millivolt (in i_Na)
double i_Na;   // nanoA (in i_Na)
double i_Na_;   // nanoA (in i_Na)
double i_Na_L;   // nanoA (in i_Na)
double ACh_shift;   // millivolt (in i_f_y_gate)
double Iso_shift_2;   // millivolt (Iso_shift in i_f_y_gate)
double tau_y;   // second (in i_f_y_gate)
double y_infinity;   // dimensionless (in i_f_y_gate)
double G_f;   // microS (in i_f)
double G_f_K;   // microS (in i_f)
double G_f_Na;   // microS (in i_f)
double g_f_K;   // microS (in i_f)
double g_f_Na;   // microS (in i_f)
double i_f;   // nanoA (in i_f)
double i_fK;   // nanoA (in i_f)
double i_fNa;   // nanoA (in i_f)
double q_infinity;   // dimensionless (in i_to_q_gate)
double tau_q;   // second (in i_to_q_gate)
double r_infinity;   // dimensionless (in i_to_r_gate)
double tau_r;   // second (in i_to_r_gate)
double i_to;   // nanoA (in i_to)
int g[X+1][Y+1][Z+1];
int h[X+1][Y+1][Z+1];
double xx[X+1][Y+1][Z+1];
double yy[X+1][Y+1][Z+1];
double zz[X+1][Y+1][Z+1];
double dc[X+1][Y+1][Z+1][10];
double df[X+1][Y+1][Z+1][10];
double D1;
double D2;
double DD;

int main (int argc, char **argv)
{
	int x,y,z, i,j,k;
	int num;
	int isbound;
	int gg[27];
	float root2 = sqrt(2.0);
	float root3 = sqrt(3.0);
	float ic, ir, il, imax;
	float tflt;
	double du;
	double dudx, dudy, dudz;
	double dudx2, dudy2, dudz2;
	double dudxdy,dudxdz, dudydz;
	char *str1;
	FILE *out1;
	FILE *in;
	double t, udt;    
	double end_time;
	double steps; 
	long int increment;

	t = 0.0; /* Time (s) */
	udt = 0.0000025;/* Time step (s) */
	end_time = 1.5; //s
	steps = end_time/udt;
	
	char c;
	double gx = 0.0, gy = 0.0, gz = 0.0;
	for (z = 0; z < Z; z++) {
		for (y = 0; y < Y; y++) {
			for (x = 0; x < X; x++) {     
				g[x][y][z] = 0;
				xx[x][y][z] = 0.0;
				yy[x][y][z] = 0.0;
				zz[x][y][z] = 0.0;
			}
		}
	}
	for(x = 2; x <= X-2; x++){
		for(y = 2; y <= Y-2; y++){
			for(z = 2; z <= Z-2; z++){
				g[x][y][z] = 2;
			}
		}
	}

	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++)
				if (g[x][y][z] > 0)
					h[x][y][z] = 1;

	num = 1;      
	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++)
				if (h[x][y][z] > 0){
					h[x][y][z] = num;
					num++;
				}

	for (x = 1; x < X; x++) 
		for (y = 1; y < Y; y++) 
			for (z = 1; z < Z; z++){
				gg[1] = h[x - 1][y - 1][z - 1];
				gg[2] = h[x - 1][y - 1][z];
				gg[3] = h[x - 1][y - 1][z + 1];
				gg[4] = h[x - 1][y][z - 1];
				gg[5] = h[x - 1][y][z];
				gg[6] = h[x - 1][y][z + 1];
				gg[7] = h[x - 1][y + 1][z - 1];
				gg[8] = h[x - 1][y + 1][z];
				gg[9] = h[x - 1][y + 1][z + 1];

				gg[10] = h[x][y - 1][z - 1];
				gg[11] = h[x][y - 1][z];
				gg[12] = h[x][y - 1][z + 1];
				gg[13] = h[x][y][z - 1];
				gg[14] = h[x][y][z + 1];
				gg[15] = h[x][y + 1][z - 1];
				gg[16] = h[x][y + 1][z];
				gg[17] = h[x][y + 1][z + 1];

				gg[18] = h[x + 1][y - 1][z - 1];
				gg[19] = h[x + 1][y - 1][z];
				gg[20] = h[x + 1][y - 1][z + 1];
				gg[21] = h[x + 1][y][z - 1];
				gg[22] = h[x + 1][y][z];
				gg[23] = h[x + 1][y][z + 1];
				gg[24] = h[x + 1][y + 1][z - 1];
				gg[25] = h[x + 1][y + 1][z];
				gg[26] = h[x + 1][y + 1][z + 1];

				isbound = 0;
				for(i = 1; i <= 26; i++){ 
					if (gg[i] > 0){
						gg[i] = 1;
						isbound++;
					} 
					else
						gg[i] = 0;
				}

				if (h[x][y][z] == 0 && isbound > 0) {
					ic = (gg[3]/root3) - (gg[1]/root3) + (gg[6]/root2) +
						(gg[9]/root3) - (gg[7]/root3) - (gg[4]/root2) +
						(gg[12]/root2) - (gg[10]/root2) + gg[14] +
						(gg[17]/root2) - (gg[15]/root2) - gg[13] +
						(gg[20]/root3) - (gg[18]/root3) + (gg[23]/root2) +
						(gg[26]/root3) - (gg[24]/root3) - (gg[21]/root2);

					ir = (gg[9]/root3) - (gg[2]/root2) - (gg[3]/root3) - 
						(gg[1]/root3) + (gg[8]/root2) + (gg[7]/root3) +
						(gg[17]/root2) - gg[11] - (gg[12]/root2) -
						(gg[10]/root2) + gg[16] + (gg[15]/root2) +
						(gg[26]/root3) - (gg[19]/root2) - (gg[20]/root3) -
						(gg[18]/root3) + (gg[25]/root2) + (gg[24]/root3);

					il = (gg[18]/root3) + (gg[19]/root2) + (gg[20]/root3) +
						(gg[21]/root2) + gg[22] + (gg[23]/root2) +
						(gg[24]/root3) + (gg[25]/root2) + (gg[26]/root3) -
						(gg[1]/root3) - (gg[2]/root2) - (gg[3]/root3) -
						(gg[4]/root2) - gg[5] - (gg[6]/root2) - 
						(gg[7]/root3) - (gg[8]/root2) - (gg[9]/root3);

					imax = fabs(ic);
					if (fabs(ir) > imax)
						imax = fabs(ir);
					if (fabs(il) > imax)
						imax = fabs(il);

					i = 0;
					j = 0; 
					k = 0;

					tflt = ir / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5) 
						i = 0;
					else if (tflt > 0.5) 
						i = 1;
					else if (tflt < -0.5) 
						i = -1;

					tflt = ic / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5) 
						j = 0;
					else if (tflt > 0.5) 
						j = 1;
					else if (tflt < -0.5) 
						j = -1;

					tflt = il / fabs(imax);
					if (tflt <= 0.5 && tflt >= -0.5)
						k = 0;
					else if (tflt > 0.5)
						k = 1;
					else if (tflt < -0.5)
						k = -1;

					if (imax == 0){
						i = 0;
						j = 0;
						k = 0;
					}
					if (h[x + k][y + i][z + j] > 0)    
						h[x][y][z] = -1 * h[x + k][y + i][z + j];   
					else
						h[x][y][z] = h[x + k][y + i][z + j];
				}
			}

	for (z = 0; z < Z; z++) {
		for (y = 0; y < Y; y++) {
			for (x = 0; x < X; x++) {
				if(g[x][y][z] > 0){
					Y0[x][y][z]  = 1.09961e-9;   
					Y1[x][y][z]  = 4.12507e-9;  
					Y2[x][y][z]  = 0.789534;   
					Y3[x][y][z]  = 0.210465; 
					Y4[x][y][z]  = 0.303076;   
					Y5[x][y][z]  = 0.162453; 
					Y6[x][y][z]  = 0.096923;   
					Y7[x][y][z]  = 0.0278476;   
					Y8[x][y][z]  = 0.35637;   
					Y9[x][y][z]  = 0.568582;   
					Y10[x][y][z] = 0.273853;   
					Y11[x][y][z] = 0.409962;  
					Y12[x][y][z] = 0.0000640697;   
					Y13[x][y][z] = 0.000143098;   
					Y14[x][y][z] = -58.8898;   
					Y15[x][y][z] = 5.0;   
					Y16[x][y][z] = 0.0000564044;   
					Y17[x][y][z] = 0.643182;   
					Y18[x][y][z] = 0.536436;  
					Y19[x][y][z] = 0.0231233;   
					Y20[x][y][z] = 0.405591;  
					Y21[x][y][z] = 0.00203218; 
					Y22[x][y][z] = 0.0111612;   
					Y23[x][y][z] = 0.704326;  
					Y24[x][y][z] = 0.854558;   
					Y25[x][y][z] = 0.164813;   
					Y26[x][y][z] = 0.00359723;  
					Y27[x][y][z] = 0.833405;   
					Y28[x][y][z] = 0.0461467;  
					Y29[x][y][z] = 0.116474;  
					Y30[x][y][z] = 0.00349197;   
					Y31[x][y][z] = 0.480295; 
					Y32[x][y][z] = 0.00544038;  
				}
			}
		}
	}

	EC50_SR = 0.45;   // millimolar (in Ca_SR_release)
	HSR = 2.5;   // dimensionless (in Ca_SR_release), gets hardcoded for performance
	MaxSR = 15.0;   // dimensionless (in Ca_SR_release)
	MinSR = 1.0;   // dimensionless (in Ca_SR_release)
	kiCa = 500.0;   // per_millimolar_second (in Ca_SR_release)
	kim = 5.0;   // per_second (in Ca_SR_release)
	koCa = 10000.0;   // per_millimolar2_second (in Ca_SR_release)
	kom = 660.0;   // per_second (in Ca_SR_release)
	ks = 148041085.1;   // per_second (in Ca_SR_release)
	CM_tot = 0.045;   // millimolar (in Ca_buffering)
	CQ_tot = 10.0;   // millimolar (in Ca_buffering)
	Mgi = 2.5;   // millimolar (in Ca_buffering)
	TC_tot = 0.031;   // millimolar (in Ca_buffering)
	TMC_tot = 0.062;   // millimolar (in Ca_buffering)
	kb_CM = 542.0;   // per_second (in Ca_buffering)
	kb_CQ = 445.0;   // per_second (in Ca_buffering)
	kb_TC = 446.0;   // per_second (in Ca_buffering)
	kb_TMC = 7.51;   // per_second (in Ca_buffering)
	kb_TMM = 751.0;   // per_second (in Ca_buffering)
	kf_CM = 1.642e6;   // per_millimolar_second (in Ca_buffering)
	kf_CQ = 175.4;   // per_millimolar_second (in Ca_buffering)
	kf_TC = 88800.0;   // per_millimolar_second (in Ca_buffering)
	kf_TMC = 227700.0;   // per_millimolar_second (in Ca_buffering)
	kf_TMM = 2277.0;   // per_millimolar_second (in Ca_buffering)
	K_up = 0.000286113;   // millimolar (in Ca_intracellular_fluxes)
	P_up_basal = 5.0;   // millimolar_per_second (in Ca_intracellular_fluxes)
	slope_up = 5.0e-5;   // millimolar (in Ca_intracellular_fluxes)
	tau_dif_Ca = 5.469e-5;   // second (in Ca_intracellular_fluxes)
	tau_tr = 0.04;   // second (in Ca_intracellular_fluxes)
	L_cell = 67.0;   // micrometre (in Cell_parameters)
	L_sub = 0.02;   // micrometre (in Cell_parameters)
	R_cell = 3.9;   // micrometre (in Cell_parameters)
	V_i_part = 0.46;   // dimensionless (in Cell_parameters)
	V_jsr_part = 0.0012;   // dimensionless (in Cell_parameters)
	V_nsr_part = 0.0116;   // dimensionless (in Cell_parameters)
	Cao = 1.8;   // millimolar (in Ionic_values)
	
	Ko = 5.4;   // millimolar (in Ionic_values)
	Nao = 140.0;   // millimolar (in Ionic_values)
	C = 5.7e-5;   // microF (in Membrane)
	F = 96485.0;//96485.3415;   // coulomb_per_mole (in Membrane)
	R2 = 8314.0;//8314.472;   // joule_per_kilomole_kelvin (R in Membrane)
	T = 310.0;   // kelvin (in Membrane)
	clamp_mode = 0.0;   // dimensionless (in Membrane)
	Nai_clamp = 1.0;   // dimensionless (in Nai_concentration)
	ACh = 0.0;//2.5e-6*1;//10.0e-6;   // millimolar (in Rate_modulation_experiments)
	Iso_1_uM = 0.0;   // dimensionless (in Rate_modulation_experiments)
	V_holding = -45.0;   // millivolt (in Voltage_clamp)
	V_test = -35.0;   // millivolt (in Voltage_clamp)
	t_holding = 0.5;   // second (in Voltage_clamp)
	t_test = 0.5;   // second (in Voltage_clamp)
	V_dL = -16.4508;   // millivolt (in i_CaL_dL_gate)
	k_dL = 4.3371;   // millivolt (in i_CaL_dL_gate)
	Km_fCa = 0.000338;   // millimolar (in i_CaL_fCa_gate)
	alpha_fCa = 0.0075;   // per_second (in i_CaL_fCa_gate)
	k_fL = 0.0;   // millivolt (in i_CaL_fL_gate)
	shift_fL = 0.0;   // millivolt (in i_CaL_fL_gate)
	P_CaL = 0.4578;   // nanoA_per_millimolar (in i_CaL)
	offset_fT = 0.0;   // second (in i_CaT_fT_gate)
	P_CaT = 0.04132;   // nanoA_per_millimolar (in i_CaT)
	ACh_on = 1.0;   // dimensionless (in i_KACh)
	g_KACh = 0.00345;   // microS (in i_KACh)
	g_Kr = 0.00424;   // microS (in i_Kr)
	g_Ks_ = 0.00065;   // microS (in i_Ks)
	g_Kur = 0.1539e-3;   // microS (in i_Kur)
	K1ni = 395.3;   // millimolar (in i_NaCa)
	K1no = 1628.0;   // millimolar (in i_NaCa)
	K2ni = 2.289;   // millimolar (in i_NaCa)
	K2no = 561.4;   // millimolar (in i_NaCa)
	K3ni = 26.44;   // millimolar (in i_NaCa)
	K3no = 4.663;   // millimolar (in i_NaCa)
	K_NaCa = 3.343;   // nanoA (in i_NaCa)
	Kci = 0.0207;   // millimolar (in i_NaCa)
	Kcni = 26.44;   // millimolar (in i_NaCa)
	Kco = 3.663;   // millimolar (in i_NaCa)
	Qci = 0.1369;   // dimensionless (in i_NaCa)
	Qco = 0.0;   // dimensionless (in i_NaCa)
	Qn = 0.4315;   // dimensionless (in i_NaCa)
	blockade_NaCa = 0.0;   // dimensionless (in i_NaCa)
	Km_Kp = 1.4;   // millimolar (in i_NaK)
	Km_Nap = 14.0;   // millimolar (in i_NaK)
	i_NaK_max = 0.08105;   // nanoA (in i_NaK)
	delta_m = 1.0e-5;   // millivolt (in i_Na_m_gate)
	g_Na = 0.0223;   // microS (in i_Na)
	g_Na_L = 0.0;   // microS (in i_Na)
	y_shift = 0.0;   // millivolt (in i_f_y_gate)
	Km_f = 45.0;   // millimolar (in i_f)
	alpha = 0.5927;   // dimensionless (in i_f)
	blockade = 0.0;   // dimensionless (in i_f)
	g_f = 0.00427;   // microS (in i_f)
	g_to = 3.5e-3;   // microS (in i_to)

	if (Iso_1_uM > 0.0)
		b_up = -0.25;
	else if (ACh > 0.0)
		b_up = 0.7*ACh/(0.00009+ACh);
	else
		b_up = 0.0;
	P_up = P_up_basal*(1.0-b_up);
	V_cell = 1e-9*3.14159265358979*pow(R_cell, 2.0)*L_cell; //mm^3
	V_sub = 1e-9*2.0*3.14159265358979*L_sub*(R_cell-L_sub/2.0)*L_cell; //mm^3
	V_nsr = V_nsr_part*V_cell; //mm^3
	V_i = V_i_part*V_cell-V_sub;
	V_jsr = V_jsr_part*V_cell;
	RTonF = R2*T/F;
	k34 = Nao/(K3no+Nao);
	G_f = g_f/(Ko/(Ko+Km_f));
	G_f_K = G_f/(alpha+1.0);
	G_f_Na = alpha*G_f_K;
	g_f_Na = G_f_Na*Ko/(Ko+Km_f);
	g_f_K = G_f_K*Ko/(Ko+Km_f);

	ACh_block = 0.31*ACh/(ACh+0.00009);

	if (Iso_1_uM > 0.0){
		g_Ks = 1.2*g_Ks_;
		Iso_increase_2 = 1.2;
		Iso_increase_1 = 1.23;
		Iso_shift_dL = -8.0;
		Iso_slope_dL = -27.0;
		Iso_shift_1 = -14.0;
		Iso_shift_2 = 7.5;
	}
	else{
		g_Ks = g_Ks_;
		Iso_increase_2 = 1.0;
		Iso_increase_1 = 1.0;
		Iso_shift_dL = 0.0;
		Iso_slope_dL = 0.0;
		Iso_shift_1 = 0.0;
		Iso_shift_2 = 0.0;
	}

	alpha_b = (3.5988-0.025641)/(1.0+0.0000012155/pow(1.0*ACh, 1.6951))+0.025641;
	if (ACh > 0.0)
		ACh_shift = -1.0-9.898*pow(1.0*ACh, 0.618)/(pow(1.0*ACh, 0.618)+0.00122423);
	else
		ACh_shift = 0.0;
	
	D1 = D1S;
	D2 = D2S;
	DD = DDS;
	for (x = 1; x < X; x++)
    	for (y = 1; y < Y; y++)
          	for (z = 1; z < Z; z++)
      			if (g[x][y][z] > 0) {
					dc[x][y][z][1] = D2 + (DD * xx[x][y][z] * xx[x][y][z]);
					dc[x][y][z][2] = DD * xx[x][y][z] * yy[x][y][z];
					dc[x][y][z][3] = DD * xx[x][y][z] * zz[x][y][z];
					dc[x][y][z][4] = DD * yy[x][y][z] * xx[x][y][z];
					dc[x][y][z][5] = D2 + (DD * yy[x][y][z] * yy[x][y][z]);
					dc[x][y][z][6] = DD * yy[x][y][z] * zz[x][y][z];
					dc[x][y][z][7] = DD * zz[x][y][z] * xx[x][y][z];
					dc[x][y][z][8] = DD * zz[x][y][z] * yy[x][y][z];
					dc[x][y][z][9] = D2 + (DD * zz[x][y][z] * zz[x][y][z]);

					df[x][y][z][1] = (xx[x + 1][y][z] - xx[x - 1][y][z]) / (2*dx);
					df[x][y][z][2] = (xx[x][y + 1][z] - xx[x][y - 1][z]) / (2*dx);
					df[x][y][z][3] = (xx[x][y][z + 1] - xx[x][y][z - 1]) / (2*dx);
					df[x][y][z][4] = (yy[x + 1][y][z] - yy[x - 1][y][z]) / (2*dx);
					df[x][y][z][5] = (yy[x][y + 1][z] - yy[x][y - 1][z]) / (2*dx);
					df[x][y][z][6] = (yy[x][y][z + 1] - yy[x][y][z - 1]) / (2*dx);
					df[x][y][z][7] = (zz[x + 1][y][z] - zz[x - 1][y][z]) / (2*dx);
					df[x][y][z][8] = (zz[x][y + 1][z] - zz[x][y - 1][z]) / (2*dx);
					df[x][y][z][9] = (zz[x][y][z + 1] - zz[x][y][z - 1]) / (2*dx); 
				}

	for (increment = 0; increment < steps; increment++){     
		for (x = 1; x < X; x++) {
			for (y = 1; y < Y; y++){
				for (z = 1; z < Z; z++){ 
					if (h[x][y][z] < 0){
						for (i = -1; i <= 1; i++){
							for (j = -1; j <= 1; j++){
								for (k = -1; k <= 1; k++){
									if (h[x][y][z] == -h[x + i][y + j][z + k]){
										Y14[x][y][z] = Y14[x + i][y + j][z + k];
									}
								}
							}
						}
					}
				}
			}
		}
		for (x = 1; x < X; x++) {
			for (y = 1; y < Y; y++){
				for (z = 1; z < Z; z++){
					if (g[x][y][z] == 2){
						Ki = 140.0;   // millimolar (in Ionic_values)
						E_K = RTonF*log(Ko/Ki);
						j_SRCarel = ks*Y1[x][y][z]*(Y10[x][y][z]-Y12[x][y][z]);
						kCaSR = MaxSR-(MaxSR-MinSR)/(1.0+pow(EC50_SR/Y10[x][y][z], HSR));
						koSRCa = koCa/kCaSR;
						kiSRCa = kiCa*kCaSR;
						dY2 = kim*Y3[x][y][z]-kiSRCa*Y12[x][y][z]*Y2[x][y][z]-(koSRCa*pow(Y12[x][y][z], 2.0)*Y2[x][y][z]-kom*Y1[x][y][z]);						
						dY1 = koSRCa*pow(Y12[x][y][z], 2.0)*Y2[x][y][z]-kom*Y1[x][y][z]-(kiSRCa*Y12[x][y][z]*Y1[x][y][z]-kim*Y0[x][y][z]);						
						dY0 = kiSRCa*Y12[x][y][z]*Y1[x][y][z]-kim*Y0[x][y][z]-(kom*Y0[x][y][z]-koSRCa*pow(Y12[x][y][z], 2.0)*Y3[x][y][z]);						
						dY3 = kom*Y0[x][y][z]-koSRCa*pow(Y12[x][y][z], 2.0)*Y3[x][y][z]-(kim*Y3[x][y][z]-kiSRCa*Y12[x][y][z]*Y2[x][y][z]);						
						delta_fTC = kf_TC*Y13[x][y][z]*(1.0-Y7[x][y][z])-kb_TC*Y7[x][y][z];
						dY7 = delta_fTC;
						delta_fTMC = kf_TMC*Y13[x][y][z]*(1.0-(Y8[x][y][z]+Y9[x][y][z]))-kb_TMC*Y8[x][y][z];
						dY8 = delta_fTMC;
						delta_fTMM = kf_TMM*Mgi*(1.0-(Y8[x][y][z]+Y9[x][y][z]))-kb_TMM*Y9[x][y][z];
						dY9 = delta_fTMM;
						delta_fCMi = kf_CM*Y13[x][y][z]*(1.0-Y4[x][y][z])-kb_CM*Y4[x][y][z];
						dY4 = delta_fCMi;
						delta_fCMs = kf_CM*Y12[x][y][z]*(1.0-Y5[x][y][z])-kb_CM*Y5[x][y][z];
						dY5 = delta_fCMs;
						delta_fCQ = kf_CQ*Y10[x][y][z]*(1.0-Y6[x][y][z])-kb_CQ*Y6[x][y][z];
						dY6 = delta_fCQ;
						j_Ca_dif = (Y12[x][y][z]-Y13[x][y][z])/tau_dif_Ca;
						j_up = P_up/(1.0+exp((-Y13[x][y][z]+K_up)/slope_up));
						dY13 = 1.0*(j_Ca_dif*V_sub-j_up*V_nsr)/V_i-(CM_tot*delta_fCMi+TC_tot*delta_fTC+TMC_tot*delta_fTMC);
						if ((udt > t_holding) && (udt < t_holding+t_test))
							V_clamp = V_test;
						else
							V_clamp = V_holding;

						if (clamp_mode >= 1.0)
							V = V_clamp;
						else
							V = Y14[x][y][z];
						i_siCa = 2.0*P_CaL*(V-0.0)/(RTonF*(1.0-exp(-1.0*(V-0.0)*2.0/RTonF)))*(Y12[x][y][z]-Cao*exp(-2.0*(V-0.0)/RTonF))*Y16[x][y][z]*Y18[x][y][z]*Y17[x][y][z];
						i_CaT = 2.0*P_CaT*V/(RTonF*(1.0-exp(-1.0*V*2.0/RTonF)))*(Y12[x][y][z]-Cao*exp(-2.0*V/RTonF))*Y19[x][y][z]*Y20[x][y][z];
						k32 = exp(Qn*V/(2.0*RTonF));
						Nai = Y15[x][y][z];
						k43 = Nai/(K3ni+Nai);
						di = 1.0+Y12[x][y][z]/Kci*(1.0+exp(-Qci*V/RTonF)+Nai/Kcni)+Nai/K1ni*(1.0+Nai/K2ni*(1.0+Nai/K3ni));
						k14 = Nai/K1ni*Nai/K2ni*(1.0+Nai/K3ni)*exp(Qn*V/(2.0*RTonF))/di;
						k12 = Y12[x][y][z]/Kci*exp(-Qci*V/RTonF)/di;
						k41 = exp(-Qn*V/(2.0*RTonF));
						x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
						dodo = 1.0+Cao/Kco*(1.0+exp(Qco*V/RTonF))+Nao/K1no*(1.0+Nao/K2no*(1.0+Nao/K3no));
						k21 = Cao/Kco*exp(Qco*V/RTonF)/dodo;
						k23 = Nao/K1no*Nao/K2no*(1.0+Nao/K3no)*exp(-Qn*V/(2.0*RTonF))/dodo;
						x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
						x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
						x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);
						i_NaCa = (1.0-blockade_NaCa)*K_NaCa*(x2*k21-x1*k12)/(x1+x2+x3+x4);
						dY12 = j_SRCarel*V_jsr/V_sub - ((i_siCa+i_CaT-2.0*i_NaCa)/(2.0*F*V_sub) + j_Ca_dif+CM_tot*delta_fCMs);
						j_tr = (Y11[x][y][z]-Y10[x][y][z])/tau_tr;
						dY11 = j_up - j_tr*V_jsr/V_nsr;
						dY10 = j_tr-(j_SRCarel+CQ_tot*delta_fCQ);
						E_Na = RTonF*log(Nao/Nai);
						E_Ca = 0.5*RTonF*log(Cao/Y12[x][y][z]);
						i_fNa = Y30[x][y][z]*g_f_Na*(V-E_Na)*(1.0-blockade);
						i_fK = Y30[x][y][z]*g_f_K*(V-E_K)*(1.0-blockade);
						i_f = i_fNa+i_fK;
						i_Kr = g_Kr*(V-E_K)*(0.9*Y22[x][y][z]+0.1*Y23[x][y][z])*Y24[x][y][z];
						E_Ks = RTonF*log((Ko+0.12*Nao)/(Ki+0.12*Nai));
						i_Ks = g_Ks*(V-E_Ks)*pow(Y25[x][y][z], 2.0);					
						i_to = g_to*(V-E_K)*Y31[x][y][z]*Y32[x][y][z];
						i_NaK = Iso_increase_2*i_NaK_max*pow(1.0+pow(Km_Kp/Ko, 1.2), -1.0)*pow(1.0+pow(Km_Nap/Nai, 1.3), -1.0)*pow(1.0+exp(-(V-E_Na+110.0)/20.0), -1.0);						
						E_mh = RTonF*log((Nao+0.12*Ko)/(Nai+0.12*Ki));
						i_Na_ = g_Na*pow(Y29[x][y][z], 3.0)*Y28[x][y][z]*(V-E_mh);						
						i_Na_L = g_Na_L*pow(Y29[x][y][z], 3.0)*(V-E_mh);						
						i_Na = i_Na_+i_Na_L;
						i_siK = 0.000365*P_CaL*(V-0.0)/(RTonF*(1.0-exp(-1.0*(V-0.0)/RTonF)))*(Ki-Ko*exp(-1.0*(V-0.0)/RTonF))*Y16[x][y][z]*Y18[x][y][z]*Y17[x][y][z];
						i_siNa = 0.0000185*P_CaL*(V-0.0)/(RTonF*(1.0-exp(-1.0*(V-0.0)/RTonF)))*(Nai-Nao*exp(-1.0*(V-0.0)/RTonF))*Y16[x][y][z]*Y18[x][y][z]*Y17[x][y][z];
						i_CaL = (i_siCa+i_siK+i_siNa)*(1.0-ACh_block)*1.0*Iso_increase_1;
						if (ACh > 0.0){
							i_KACh = ACh_on*g_KACh*(V-E_K)*(1.0+exp((V+20.0)/20.0))*Y21[x][y][z];
						}
						else{
							i_KACh = 0.0;
						}

						i_Kur = g_Kur*Y26[x][y][z]*Y27[x][y][z]*(V-E_K);
						i_tot = i_f+i_Kr+i_Ks+i_to+i_NaK+i_NaCa+i_Na+i_CaL+i_CaT+i_KACh+i_Kur;
						dY14 = -i_tot/C;
						dY15 = (1.0-Nai_clamp)*-1.0*(i_Na+i_fNa+i_siNa+3.0*i_NaK+3.0*i_NaCa)/(1.0*(V_i+V_sub)*F);
						dL_infinity = 1.0/(1.0+exp(-(V-V_dL-Iso_shift_dL)/(k_dL*(1.0+Iso_slope_dL/100.0))));
						if (V == -41.8)
							adVm = -41.80001;
						else if (V == 0.0)
							adVm = 0.0;
						else if (V == -6.8)
							adVm = -6.80001;
						else
							adVm = V;
						alpha_dL = -0.02839*(adVm+41.8)/(exp(-(adVm+41.8)/2.5)-1.0)-0.0849*(adVm+6.8)/(exp(-(adVm+6.8)/4.8)-1.0);
						if (V == -1.8)
							bdVm = -1.80001;
						else
							bdVm = V;
						beta_dL = 0.01143*(bdVm+1.8)/(exp((bdVm+1.8)/2.5)-1.0);
						tau_dL = 0.001/(alpha_dL+beta_dL);
						dY16 = (dL_infinity-Y16[x][y][z])/tau_dL;
						fCa_infinity = Km_fCa/(Km_fCa+Y12[x][y][z]);
						tau_fCa = 0.001*fCa_infinity/alpha_fCa;
						dY17 = (fCa_infinity-Y17[x][y][z])/tau_fCa;
						fL_infinity = 1.0/(1.0+exp((V+37.4+shift_fL)/(5.3+k_fL)));
						tau_fL = 0.001*(44.3+230.0*exp(-pow((V+36.0)/10.0, 2.0)));
						dY18 = (fL_infinity-Y18[x][y][z])/tau_fL;
						dT_infinity = 1.0/(1.0+exp(-(V+38.3)/5.5));
						tau_dT = 0.001/(1.068*exp((V+38.3)/30.0)+1.068*exp(-(V+38.3)/30.0));
						dY19 = (dT_infinity-Y19[x][y][z])/tau_dT;
						fT_infinity = 1.0/(1.0+exp((V+58.7)/3.8));
						tau_fT = 1.0/(16.67*exp(-(V+75.0)/83.3)+16.67*exp((V+75.0)/15.38))+offset_fT;
						dY20 = (fT_infinity-Y20[x][y][z])/tau_fT;
						beta_b = 10.0*exp(0.0133*(V+40.0));
						b_infinity = alpha_b/(alpha_b+beta_b);
						tau_b = 1.0/(alpha_b+beta_b);
						dY21 = (b_infinity-Y21[x][y][z])/tau_b;
						alfapaF = 1.0/(1.0+exp(-(V+23.2)/6.6))/(0.84655354/(37.2*exp(V/11.9)+0.96*exp(-V/18.5)));
						betapaF = 4.0*((37.2*exp(V/15.9)+0.96*exp(-V/22.5))/0.84655354-1.0/(1.0+exp(-(V+23.2)/10.6))/(0.84655354/(37.2*exp(V/15.9)+0.96*exp(-V/22.5))));
						pa_infinity = 1.0/(1.0+exp(-(V+10.0144)/7.6607));
						tau_paS = 0.84655354/(4.2*exp(V/17.0)+0.15*exp(-V/21.6));
						tau_paF = 1.0/(30.0*exp(V/10.0)+exp(-V/12.0));
						dY23 = (pa_infinity-Y23[x][y][z])/tau_paS;
						dY22 = (pa_infinity-Y22[x][y][z])/tau_paF;
						tau_pi = 1.0/(100.0*exp(-V/54.645)+656.0*exp(V/106.157));
						pi_infinity = 1.0/(1.0+exp((V+28.6)/17.1));
						dY24 = (pi_infinity-Y24[x][y][z])/tau_pi;
						n_infinity = sqrt(1.0/(1.0+exp(-(V+0.6383-Iso_shift_1)/10.7071)));
						alpha_n = 28.0/(1.0+exp(-(V-40.0-Iso_shift_1)/3.0));
						beta_n = 1.0*exp(-(V-Iso_shift_1-5.0)/25.0);
						tau_n = 1.0/(alpha_n+beta_n);
						dY25 = (n_infinity-Y25[x][y][z])/tau_n;
						r_Kur_infinity = 1.0/(1.0+exp((V+6.0)/-8.6));
						tau_r_Kur = 0.009/(1.0+exp((V+5.0)/12.0))+0.0005;
						dY26 = (r_Kur_infinity-Y26[x][y][z])/tau_r_Kur;
						s_Kur_infinity = 1.0/(1.0+exp((V+7.5)/10.0));
						tau_s_Kur = 0.59/(1.0+exp((V+60.0)/10.0))+3.05;
						dY27 = (s_Kur_infinity-Y27[x][y][z])/tau_s_Kur;
						h_infinity = 1.0/(1.0+exp((V+69.804)/4.4565));
						alpha_h = 20.0*exp(-0.125*(V+75.0));
						beta_h = 2000.0/(320.0*exp(-0.1*(V+75.0))+1.0);
						tau_h = 1.0/(alpha_h+beta_h);
						dY28 = (h_infinity-Y28[x][y][z])/tau_h;
						m_infinity = 1.0/(1.0+exp(-(V+42.0504)/8.3106));
						E0_m = V+41.0;
						if (fabs(E0_m) < delta_m)
							alpha_m = 2000.0;
						else
							alpha_m = 200.0*E0_m/(1.0-exp(-0.1*E0_m));
						beta_m = 8000.0*exp(-0.056*(V+66.0));
						tau_m = 1.0/(alpha_m+beta_m);
						dY29 = (m_infinity-Y29[x][y][z])/tau_m;
						tau_y = 1.0/(0.36*(V+148.8-ACh_shift-Iso_shift_2)/(exp(0.066*(V+148.8-ACh_shift-Iso_shift_2))-1.0)+0.1*(V+87.3-ACh_shift-Iso_shift_2)/(1.0-exp(-0.2*(V+87.3-ACh_shift-Iso_shift_2))))-0.054;
						if (V < -(80.0-ACh_shift-Iso_shift_2-y_shift))
							y_infinity = 0.01329+0.99921/(1.0+exp((V+97.134-ACh_shift-Iso_shift_2-y_shift)/8.1752));
						else
							y_infinity = 0.0002501*exp(-(V-ACh_shift-Iso_shift_2-y_shift)/12.861);
						dY30 = (y_infinity-Y30[x][y][z])/tau_y;
						q_infinity = 1.0/(1.0+exp((V+49.0)/13.0));
						tau_q = 0.001*0.6*(65.17/(0.57*exp(-0.08*(V+44.0))+0.065*exp(0.1*(V+45.93)))+10.1);
						dY31 = (q_infinity-Y31[x][y][z])/tau_q;
						r_infinity = 1.0/(1.0+exp(-(V-19.3)/15.0));
						tau_r = 0.001*0.66*1.4*(15.59/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84)))+2.98);
						dY32 = (r_infinity-Y32[x][y][z])/tau_r;					

						dudx2  = (Y14[x - 1][y][z] + Y14[x + 1][y][z] - 2 * Y14[x][y][z]) / (dx*dx);          
						dudy2  = (Y14[x][y - 1][z] + Y14[x][y + 1][z] - 2 * Y14[x][y][z]) / (dx*dx);
						dudz2  = (Y14[x][y][z - 1] + Y14[x][y][z + 1] - 2 * Y14[x][y][z]) / (dx*dx);
						dudxdy = (Y14[x + 1][y + 1][z] + Y14[x - 1][y - 1][z] - Y14[x + 1][y - 1][z] - Y14[x - 1][y + 1][z])/(4*dx*dx);  
						dudxdz = (Y14[x + 1][y][z + 1] + Y14[x - 1][y][z - 1] - Y14[x + 1][y][z - 1] - Y14[x - 1][y][z + 1])/(4*dx*dx);
						dudydz = (Y14[x][y + 1][z + 1] + Y14[x][y - 1][z - 1] - Y14[x][y + 1][z - 1] - Y14[x][y - 1][z + 1])/(4*dx*dx);
						dudx   = (Y14[x + 1][y][z] - Y14[x - 1][y][z])/(2*dx);  
						dudy   = (Y14[x][y + 1][z] - Y14[x][y - 1][z])/(2*dx);
						dudz   = (Y14[x][y][z + 1] - Y14[x][y][z - 1])/(2*dx);
						du= dc[x][y][z][1]*dudx2  + dudx*DD*(xx[x][y][z]*df[x][y][z][1] + xx[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][2]*dudxdy + dudy*DD*(xx[x][y][z]*df[x][y][z][4] + yy[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][3]*dudxdz + dudz*DD*(xx[x][y][z]*df[x][y][z][7] + zz[x][y][z]*df[x][y][z][1]) +
							dc[x][y][z][4]*dudxdy + dudx*DD*(yy[x][y][z]*df[x][y][z][2] + xx[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][5]*dudy2  + dudy*DD*(yy[x][y][z]*df[x][y][z][5] + yy[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][6]*dudydz + dudz*DD*(yy[x][y][z]*df[x][y][z][8] + zz[x][y][z]*df[x][y][z][5]) +
							dc[x][y][z][7]*dudxdz + dudx*DD*(zz[x][y][z]*df[x][y][z][3] + xx[x][y][z]*df[x][y][z][9]) +
							dc[x][y][z][8]*dudydz + dudy*DD*(zz[x][y][z]*df[x][y][z][6] + yy[x][y][z]*df[x][y][z][9]) +
							dc[x][y][z][9]*dudz2  + dudz*DD*(zz[x][y][z]*df[x][y][z][9] + zz[x][y][z]*df[x][y][z][9]);

						Y0[x][y][z] = Y0[x][y][z] + udt*dY0;
						Y1[x][y][z] = Y1[x][y][z] + udt*dY1;
						Y2[x][y][z] = Y2[x][y][z] + udt*dY2;
						Y3[x][y][z] = Y3[x][y][z] + udt*dY3;
						Y4[x][y][z] = Y4[x][y][z] + udt*dY4;
						Y5[x][y][z] = Y5[x][y][z] + udt*dY5;
						Y6[x][y][z] = Y6[x][y][z] + udt*dY6;
						Y7[x][y][z] = Y7[x][y][z] + udt*dY7;
						Y8[x][y][z] = Y8[x][y][z] + udt*dY8;
						Y9[x][y][z] = Y9[x][y][z] + udt*dY9;							
						Y10[x][y][z] = Y10[x][y][z] + udt*dY10;
						Y11[x][y][z] = Y11[x][y][z] + udt*dY11;
						Y12[x][y][z] = Y12[x][y][z] + udt*dY12;
						Y13[x][y][z] = Y13[x][y][z] + udt*dY13;
						Y14[x][y][z] = Y14[x][y][z] + udt*(du + dY14);				
						Y15[x][y][z] = Y15[x][y][z] + udt*dY15;
						Y16[x][y][z] = Y16[x][y][z] + udt*dY16;
						Y17[x][y][z] = Y17[x][y][z] + udt*dY17;
						Y18[x][y][z] = Y18[x][y][z] + udt*dY18;
						Y19[x][y][z] = Y19[x][y][z] + udt*dY19;
						Y20[x][y][z] = Y20[x][y][z] + udt*dY20;
						Y21[x][y][z] = Y21[x][y][z] + udt*dY21;						
						Y22[x][y][z] = Y22[x][y][z] + udt*dY22;
						Y23[x][y][z] = Y23[x][y][z] + udt*dY23;
						Y24[x][y][z] = Y24[x][y][z] + udt*dY24;
						Y25[x][y][z] = Y25[x][y][z] + udt*dY25;
						Y26[x][y][z] = Y26[x][y][z] + udt*dY26;
						Y27[x][y][z] = Y27[x][y][z] + udt*dY27;
						Y28[x][y][z] = Y28[x][y][z] + udt*dY28;
						Y29[x][y][z] = Y29[x][y][z] + udt*dY29;
						Y30[x][y][z] = Y30[x][y][z] + udt*dY30;
						Y31[x][y][z] = Y31[x][y][z] + udt*dY31;
						Y32[x][y][z] = Y32[x][y][z] + udt*dY32;
					}
				}
			}
		}
		t = t+udt;
		if (increment % 400 == 0) {
			str1 = malloc(100 * sizeof(char));
			sprintf(str1, "test/aa%.5ld.vtk", increment / 400);
			out1 = fopen(str1, "wt");
			for (z = 0; z < Z; z++) {
				for (y = 0; y < Y; y++) {
					for (x = 0; x < X; x++) {
						fprintf(out1, "%2.3f ", Y14[x][y][z]);
					}
					fprintf(out1, "\n");
				}
				fprintf(out1, "\n");
			}
			free(str1);
			fclose(out1);
		}
		
	}
	return 0;
}
