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

#include <default_gui_model.h>

class vdend : public DefaultGUIModel {
	
	public:
		vdend(void);
		virtual ~vdend(void);
		virtual void execute(void);
	
	protected:
		virtual void update(DefaultGUIModel::update_flags_t);
	
	private:
		void solve(double *);
	
		double y[5];
	
		double dt;
		double Gna_max, Gk_max, Gh_max, Gl_max;
		double Vna, Vk, Vh, Vl;
		double rho_side, rho_cent;
		double c_dend, Ioffset;
		double twoway;
	
		double Vright, Vleft, Vcent, Iapp, Iout;
	
		int toR, toL, toC;
	
		double tau_m, tau_h, tau_n, tau_q;
		double alpha_m, alpha_h, alpha_n, alpha_q;
		double beta_m, beta_h, beta_n, beta_q;
		double V_inf, m_inf, h_inf, n_inf, q_inf;
		double Gz_right, Gz_left, Gz_cent;
		double Gm, Gna, Gk, Gh, Gl;
};
