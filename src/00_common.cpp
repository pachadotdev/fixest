#include "01_convergence.h"
#include "02_demeaning.h"

inline bool continue_criterion(double a, double b, double diffMax)
{
	// continuing criterion of the algorithm
	double diff = fabs(a - b);
	return ((diff > diffMax) && (diff / (0.1 + fabs(a)) > diffMax));
}

inline bool stopping_criterion(double a, double b, double diffMax)
{
	// stopping criterion of the algorithm
	double diff = fabs(a - b);
	return ((diff < diffMax) || (diff / (0.1 + fabs(a)) < diffMax));
}

// TODO: interruption in cpp11
// for interruption
// int pending_interrupt() {
//     if (check_user_interrupt()) {
//         return 1;
//     } else {
//         return 0;
//     }
// }

// IT update + returns numerical convergence indicator
// nb_coef_no_KQ = K or Q
bool update_X_IronsTuck(int nb_coef_no_KQ, vector<double> &X,
						const vector<double> &GX, const vector<double> &GGX,
						vector<double> &delta_GX, vector<double> &delta2_X)
{

	for (int i = 0; i < nb_coef_no_KQ; ++i)
	{
		double GX_tmp = GX[i];
		delta_GX[i] = GGX[i] - GX_tmp;
		delta2_X[i] = delta_GX[i] - GX_tmp + X[i];
		// delta_GX[i] = GGX[i] - GX[i];
		// delta2_X[i] = delta_GX[i] - GX[i] + X[i];
	}

	// delta_GX %*% delta2_X and crossprod(delta2_X)
	double vprod = 0, ssq = 0;
	for (int i = 0; i < nb_coef_no_KQ; ++i)
	{
		double delta2_X_tmp = delta2_X[i];
		vprod += delta_GX[i] * delta2_X_tmp;
		ssq += delta2_X_tmp * delta2_X_tmp;
		// vprod += delta_GX[i] * delta2_X[i];
		// ssq += delta2_X[i] * delta2_X[i];
	}

	bool res = false;

	if (ssq == 0)
	{
		res = true;
	}
	else
	{
		double coef = vprod / ssq;

		// update of X:
		for (int i = 0; i < nb_coef_no_KQ; ++i)
		{
			X[i] = GGX[i] - coef * delta_GX[i];
		}
	}

	return (res);
}
