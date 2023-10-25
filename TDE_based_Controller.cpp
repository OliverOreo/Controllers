// Flexible ureteroscopy robot
// Nankai Univeristy
// Xiangyu Wang, 2023/6/18

#include "TDE_based_Controller.h"
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <math.h>
#include <queue>

using namespace cv;
using namespace std;

// Time-delay control gain-adapting laws
void Dynamic::Adaptive_law(double & tracking_error)
{
	// Update the buffer of function memory
	if(Func.size() >= N_memory){
		Func1.pop();
		Func1.push(pow(Sign(tracking_error), alpha1));
		Func2.pop();
		Func2.push(pow(Sign(tracking_error), alpha2));
	}
	else{
		Func1.push(pow(Sign(tracking_error), alpha1));
		Func2.push(pow(Sign(tracking_error), alpha2));
	}
	
	h_a1 = h^alpha1;
	h_a2 = h^alpha2;
	
	for(int j = 0; j < N_memory; j++){
		GL_sum_s += lambda1*(-1)^j * tgamma(r1 + 1) / (tgamma(j+1) * tgamma(r1 -j + 1)) * Func1(N_memory - j) / h_a1 + lambda2*(-1)^j * tgamma(r2 - 1 + 1) / (tgamma(j+1) * tgamma(r2 - 1 -j + 1)) * Func2(N_memory - j) / h_a2;
		GL_sum_phi = lambda1*(-1)^j * tgamma(r1 + 1 + 1) / (tgamma(j+1) * tgamma(r1 + 1 -j + 1)) * Func1(N_memory - j) / h_a1 + lambda2*(-1)^j * tgamma(r2 + 1) / (tgamma(j+1) * tgamma(r2 -j + 1)) * Func2(N_memory - j) / h_a2;
	}
	// Manifold
	s = (tracking_error - last_tracking_error) / Ts + GL_sum_s;
   		
	// Gain-adapting Law
	if (xi >= xi_max || xi < 0){
		d_xi = mu1*fabs(s)*Sign(xi_max / 2 - xi);
	}
	else{
		if (fabs(s) >= Delta){
			eta = tanh(kappa*s / lambda_c);
			d_xi = mu2*(kappa*k_eta*k_eta / lambda_c)*(1 - eta*eta) / (1 - k_eta*k_eta*eta*eta)*eta / (1 - k_eta*k_eta*eta*eta);
		}
		else{
			d_xi = -mu1 / (fabs(s) + delta);
		}
	}
	xi = xi + d_xi*Ts;
	M = M0*(1 + kappa1*xi);
	k = k0*(1 + kappa2*xi);	
	
	// Temporary variable
	phi = RL_sum_phi / (1 + kappa1*xi);
	last_tracking_error = tracking_error;
	// Reaching Law
	s = s + Ts*(-w1*s - w2*pow(Sign(s), 1));
}

// Calculation of control input
double Dynamic::Controller(double & tracking_error)
{	
	ATDC_adaptive_law(tracking_error);
	u1 = ddot_ref + phi + w1*s + w2*pow(fabs(s),p)*Sign(s);
	u2 = k*fabs(s);
	u_input = u1 + u2;
	tau = M*u_input + tau - M0*ddot_beta;
	return tau;
}
