/*
* Copyright (C) 2004 Boston University
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of the
*  License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include <virtual-dendrite.h>
#include <math.h>

#define sm 0.08
#define sh -0.06
#define sn 0.07
#define sq -0.08
#define vm -39.0
#define vh -45.0
#define vn -29.8
#define vq -50.0
#define tm 1.0
#define th 400.0
#define tn 40.0
#define tq 10000.0
#define rm 0.0
#define rh 0.0
#define rn 0.0
#define rq 0.0

#define V (y[0])
#define m (y[1])
#define h (y[2])
#define n (y[3])
#define q (y[4])
#define dV (dy[0])
#define dm (dy[1])
#define dh (dy[2])
#define dn (dy[3])
#define dq (dy[4])

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new vdend();
}

static DefaultGUIModel::variable_t vars[] =
{
	{ "Vright", "", DefaultGUIModel::INPUT, },
	{ "Vleft", "", DefaultGUIModel::INPUT, },
	{ "Vcent", "", DefaultGUIModel::INPUT, },
	{ "Iapp", "", DefaultGUIModel::INPUT, },
	{ "Vold", "", DefaultGUIModel::OUTPUT, },
	{ "Iout", "", DefaultGUIModel::OUTPUT, },
	{ "Gna_max", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Vna", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Gk_max", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Vk", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Gh_max", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Vh", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Gl_max", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Vl", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "rho_side", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "rho_cent", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "c_dend", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "twoway", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Ioffset", "", DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "m", "", DefaultGUIModel::STATE, },
	{ "h", "", DefaultGUIModel::STATE, },
	{ "n", "", DefaultGUIModel::STATE, },
	{ "q", "", DefaultGUIModel::STATE, }, 
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

vdend::vdend(void) : DefaultGUIModel("Virtual Dendrite", ::vars, ::num_vars), dt(RT::System::getInstance()->getPeriod() * 1e-6), Gna_max(0.01), Gk_max(100.0), Gh_max(0.0), Gl_max(0.63), Vna(50.0), Vk(-80.0), Vh(-35.0), Vl(-51.0), rho_side(.034), rho_cent(.034), c_dend(34), Ioffset(0.0), twoway(0.0) {

	DefaultGUIModel::createGUI(vars, num_vars);
	update(INIT);
	refresh();
	QTimer::singleShot(0, this, SLOT(resizeMe()));
}

vdend::~vdend(void) {}

void vdend::solve(double *y) {
		
	/*  static double tau_m, tau_h, tau_n, tau_q;
	static double alpha_m, alpha_h, alpha_n, alpha_q;
	static double V_inf, m_inf, h_inf, n_inf, q_inf;
	static double Gz_right, Gz_left, Gz_cent;
	static double Gm, Gna, Gk, Gh, Gl;
	static double Vright, Vleft, Vcent, Iapp, Iout;
	static double toR, toL, toC;
	static int twoway;
	*/
	Vright = input(0) * 1e3;
	Vleft = input(1) * 1e3;
	Vcent = input(2) * 1e3;
	Iapp = input(3) * 1e12 + Ioffset;
	
	toR = (int) fabs(Vright) > 1e-9;
	toL = (int) fabs(Vleft) > 1e-9;
	toC = (int) fabs(Vcent) > 1e-9;
	
	Vright = Vright + V * (1 - toR);
	Vleft = Vleft + V * (1 - toL);
	Vcent = Vcent + V * (1 - toC);
	
	alpha_m = exp((V - vm) * (2 * sm - rm)) / (2 * tm);
	alpha_h = exp((V - vh) * (2 * sh - rh)) / (2 * th);
	alpha_n = exp((V - vn) * (2 * sn - rn)) / (2 * tn);
	alpha_q = exp((V - vq) * (2 * sq - rq)) / (2 * tq);
	
	beta_m = exp((-V + vm) * (2 * sm + rm)) / (2 * tm);
	beta_h = exp((-V + vh) * (2 * sh + rh)) / (2 * th);
	beta_n = exp((-V + vn) * (2 * sn + rn)) / (2 * tn);
	beta_q = exp((-V + vq) * (2 * sq + rq)) / (2 * tq);
	
	tau_m = 1.0 / (alpha_m + beta_m);
	tau_h = 1.0 / (alpha_h + beta_h);
	tau_n = 1.0 / (alpha_n + beta_n);
	tau_q = 1.0 / (alpha_q + beta_q);
	
	m_inf = alpha_m * tau_m;
	h_inf = alpha_h * tau_h;
	n_inf = alpha_n * tau_n;
	q_inf = alpha_q * tau_q;
	
	m = m_inf + (m - m_inf) * exp(-dt / tau_m);
	h = h_inf + (h - h_inf) * exp(-dt / tau_h);
	n = n_inf + (n - n_inf) * exp(-dt / tau_n);
	q = q_inf + (q - q_inf) * exp(-dt / tau_q);
	
	if (m > 1) m = 1;
	else if (m < 0) m = 0;

	if (h > 1)h = 1;
	else if (h < 0) h = 0;

	if (n > 1) n = 1;
	else if (n < 0) n = 0;

	if (q > 1) q = 1;
	else if (q < 0) q = 0;
	
	Gz_right = (1 / rho_side) * toR;
	Gz_left = (1 / rho_side) * toL;
	Gz_cent = (1 / rho_cent) * toC;
	
	Gna = Gna_max * m;
	Gk = Gk_max * h * n;
	Gh = Gh_max * q;
	Gl = Gl_max;
	
	Gm = Gz_right + Gz_left + Gz_cent + Gna + Gk + Gh + Gl;
	
	V_inf = (Gl * Vl + Gna * Vna + Gh * Vh + Gk * Vk + Gz_left * Vleft + Gz_right * Vright + twoway * Gz_cent * Vcent
	+ Iapp) / Gm;
	
	V = V_inf - (V_inf - V) * exp(-dt * Gm / c_dend);
	
	Iout = (V - Vcent) * Gz_cent;
	output(1) = Iout * 1e-12;
}

void vdend::execute(void) {
	output(0) = V * 1e-3;
	solve(y);
}

void vdend::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
		case INIT:
			setState("m", m);
			setState("h", h);
			setState("n", n);
			setState("q", V);
		
			setParameter("Gna_max", Gna_max);
			setParameter("Vna", Vna);
			setParameter("Gk_max", Gk_max);
			setParameter("Vk", Vk);
			setParameter("Gh_max", Gh_max);
			setParameter("Vh", Vh);
			setParameter("Gl_max", Gl_max);
			setParameter("Vl", Vl);
			setParameter("rho_side", rho_side * 1e3);
			setParameter("rho_cent", rho_cent * 1e3);
			setParameter("c_dend", c_dend);
			setParameter("twoway", twoway);
			setParameter("Ioffset", Ioffset);
			break;
	
		case MODIFY:
			Gna_max = getParameter("Gna_max").toDouble();
			Vna = getParameter("Vna").toDouble();
			Gk_max = getParameter("Gk_max").toDouble();
			Vk = getParameter("Vk").toDouble();
			Gh_max = getParameter("Gh_max").toDouble();
			Vh = getParameter("Vh").toDouble();
			Gl_max = getParameter("Gl_max").toDouble();
			Vl = getParameter("Vl").toDouble();
			rho_side = getParameter("rho_side").toDouble() * 1e-3;
			rho_cent = getParameter("rho_cent").toDouble() * 1e-3;
			c_dend = getParameter("c_dend").toDouble();
			twoway = getParameter("twoway").toDouble();
			Ioffset = getParameter("Ioffset").toDouble();
			break;
	
		case PERIOD:
			dt = RT::System::getInstance()->getPeriod() * 1e-6;
			break;
		
		default:
			break;
	}
}
