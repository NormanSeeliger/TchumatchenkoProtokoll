/********************************************
 *** Useful generic functions and classes ***
 ***    (c) Jannik Luboeinski 2015/2016	***
 ********************************************/

// namespace std is required

#include <vector> // therefore, add -std=c++0x to the command line
#include <chrono> // for dateStr() and getClockSeed()
#include <limits> // for maximum/minimum magnitude of integers
#include <cmath>
#include <cstring>
//#include <sys/ioctl.h> // for getLinebreak()
#include <errno.h>
#include <memory> // for unique pointers
//#include <algorithm>  // for max_element, min_element

#if defined(NETWORKSIMULATION) || defined(NEURONSIMULATION) // GSL functions for fitting - therefore, link with "-lgsl -lgslcblas"
/*#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h> */
#endif

/*** Class idpair ***
 * a pair of an integer and a double / co-ordinate values */		
class idpair
{
public:
	int x;
	double y;

	bool equals(idpair idp)
	{
		if (idp.x == x && idp.y == y)
			return true;
		else
			return false;
	}
	idpair(int a, double b)
	{
		x = a;
		y = b;
	}
	idpair() // necessary for creating an array of idpair objects
	{
		x = 0;
		y = 0.0;
	}
};

#if defined(NETWORKSIMULATION) || defined(NEURONSIMULATION)
/*** Structure fitdata ***
 * contains arrays for the data points (x, y, yerror) *
 * is used for GSL fits (e.g. in findFit()) */	
struct fitdata_pairs {
	size_t n; // number of values
	double *x; // data x values
	double *y; // data y values
	double *sigma; // data y error estimates
};

/*** Structure fitdata_idpairs ***
 * contains arrays for the data points (x, y, yerror) [where x is an integer]*
 * is used for GSL fits (e.g. in findFit()) */	
struct fitdata_idpairs {
	size_t n; // number of values
	idpair *val; // data (x,y) values
	double *sigma; // data y error estimates
};

/*** findFit ***
 * Fits an arbitrary nonlinear function on a given set of data points using GSL; returns the fitted parameters *
 * in par and their errors in err (the constant parameters are left untouched) *
 * - func: the gsl fit function, which has to be set up beforehand *
 * - par: all necessary parameters (both fit parameters and constant parameters), has to contain the initial values upon calling *
 * - err: parameter error estimates (are set in the end of the fit procedure) *
 * - chisq: if pointer is not NULL, chi squared is written to it after the fit was successful *
 * - return: number of iterations if the fit converged, 0 if not */
/*int findFit(gsl_multifit_function_fdf &func, double *par, double *err, double *chisq)
{
	/*const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; // Levenberg-Marquardt derivative solver
	const double abs_err = 1e-4; // absolute error tolerance for exit condition
	const double rel_err = 1e-4; // relative error tolerance for exit condition
	int status; // status of GSL fitting methods (error output)
	int iter = 0; // iteration number
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc(T, func.n, func.p); // solver of type T for n data points and p parameters
	gsl_matrix *covar = gsl_matrix_alloc(func.p, func.p); // covariance matrix
	gsl_vector_view x = gsl_vector_view_array(par, func.p); // vector view on array of initial values

	gsl_multifit_fdfsolver_set(s, &func, &x.vector); // initialize the solver using the function and derivative and the initial guess x

	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate(s);

		if (status)
			break;

		status = gsl_multifit_test_delta(s->dx, s->x, abs_err, rel_err); // check for convergence depending on absolute and relative error

	} while (status == GSL_CONTINUE && iter < 500);

	if (iter >= 500)
	{
		std::cout << "Fit error: fit has not converged after 500 iterations." << std::endl;
		return 0;
	}
	else if (status != GSL_SUCCESS)
	{
		std::cout << "Fit error: " << gsl_strerror(status) << std::endl;
		return 0;
	}

	gsl_multifit_covar(s->J, 0.0, covar); // computes the covariance of best-fit parameters

	for (int i=0; i<func.p; i++) // writes the resulting values to the given arrays
	{
		par[i] = gsl_vector_get(s->x, i);
		err[i] = sqrt(gsl_matrix_get(covar,i,i));
	}
	
	if (chisq)
		*chisq = pow(gsl_blas_dnrm2(s->f), 2); // computes chi squared (gsl_blas_dnrm2 returns the Euclidian norm of s->f)

	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	return 0;
}*/
#endif

/*** concat ***
 * Concatenates two character arrays/strings and returns a string *
 * - str1: the first character array
 * - str2: the second character array
 * - return: the concatenated character array in a new string variable */
string concat(const char *str1, const char *str2)
{
	string s = string(str1) + string(str2);
	return s;
}
string concat(string str1, const char *str2)
{
	string s = str1 + string(str2);
	return s;
}
string concat(const char *str1, string str2)
{
	string s = string(str1) + str2;
	return s;
}

/*** dtos ***
 * Converts a double variable to a string, rounded to the given number of internal decimal places, *
 * using the C method sprintf *
 * - num: the double floating point number to be converted *
 * - n: the number of internal decimal places (has to be greater than or equal to 0) *
 * - dont_allow_zero: specifies if zeros shall be replaced by "nan" values *
 * - return: the string containing the rounded floating point number */
string dtos(double num, int n, bool dont_allow_zero = false)
{
	string s; // the string to be returned
	char *buf; // the buffer the converted number is written to
	char *format; // the format string for sprintf
	int prepoint; // number of places before point, has to be at least 1
	int dn = (n > 0 ? dn = (int) floor(log10(n)) + 1 : dn = 1); // decimal places of the variable n (needed to generate the format string)

	if (num <= 0.0 || floor(log10(num)) < 0.0)
		prepoint = 1;
	else
		prepoint = (int) floor(log10(num)) + 1;

	buf = new char[prepoint + n + 2]; // decimal places plus internal decimal places of num plus 2 for point and zero-terminator
	format = new char[dn + 4]; // decimal places of n plus 4 for the characters "%.f\0"

	sprintf(format, "%%.%df", n); // create format string
	sprintf(buf, (const char*) format, num); // create character array from double
	s.assign(buf); // create buffer

	if (num == 0. && dont_allow_zero)
		s.assign("nan");	

	delete[] buf;
	delete[] format;
	return s;
}

/*** dateStr ***
 * Adds a string containing a time stamp of 18 char's *
 * to the start of another string *
 * - str: the character array to which the time stamp will be added  *
 * - fixdate [optional]: true if this time shall be fixed, false by default *
 * - return: the date string plus the specified string, if stated */
string dateStr(string str, bool fixdate = false) 
{
	const int len = 19; // time stamp length: 18 characters (plus terminating \0)
	static string datestr;

	if (fixdate == true) // if this date is to be fixed, update date
	{
		char buf[len];
		time_t t = time(NULL);
   	struct tm *now  = localtime(&t);
		strftime(buf, len, "%y-%m-%d_%H-%M-%S", now); 
		datestr.assign(buf);
	}

	return datestr + string(str);
}

/*** stod_acc ***
 * String to double conversion alternative to stod *
 * - str: a string containing only a floating point number, must contain a point
 * - return: the string converted into a double */
double stod_acc(const string &str) 
{
	int place = str.find(".")-1; // decimal power
	double d = 0.0;
	
	for (int i=0; i<str.length(); i++)
	{
		double base = pow(10, place);

		if (str[i] == '1')
			d += 1.*base;
		else if (str[i] == '2')
			d += 2.*base;
		else if (str[i] == '3')
			d += 3.*base;
		else if (str[i] == '4')
			d += 4.*base;
		else if (str[i] == '5')
			d += 5.*base;
		else if (str[i] == '6')
			d += 6.*base;
		else if (str[i] == '7')
			d += 7.*base;
		else if (str[i] == '8')
			d += 8.*base;
		else if (str[i] == '9')
			d += 9.*base;
		else if (str[i] == '.')
			continue;
		place--;
	}
	return d;
}

/*** min ***
 * Finds the minimum value of a double array/a vector/an idpair array/two numbers *
 * - obj: an array of type double or a vector of type idpair or an array of type idpair (considering y)
 * - len: the number of array elements
 * - return: the minimum value */
double min (double *obj, int len)
{
	double min_n = (len > 0) ? obj[0] : std::numeric_limits<double>::lowest();
	for (int i=1; i < len; i++)
	{
		if (obj[i] < min_n)
			min_n = obj[i];
	}
	return min_n;
}
idpair min (vector<idpair> obj)
{
	double min_n = (obj.size() > 0) ? obj.front().y : std::numeric_limits<double>::lowest();
	int x = 0;
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i).y < min_n)
		{
			min_n = obj.at(i).y;
			x = obj.at(i).x;
		}
	}
	return idpair(x, min_n);
	//return *min_element(&obj.front().y, &obj.back().y); // SEEMS to work
}
idpair min (idpair* obj, int len)
{
	double max_n = (len > 0) ? obj[0].y : std::numeric_limits<double>::lowest();
	int x = 0;
	for (int i=1; i < len; i++)
	{
		if (obj[i].y < max_n)
		{
			max_n = obj[i].y;
			x = i;
		}
	}
	return idpair(x, max_n);
}

/*** max ***
 * Finds the maximum value of a double array *
 * - obj: an array of type double or a vector of type idpair or an array of type idpair (considering y)
 * - len: the number of array elements
 * - return: the maximum value */
double max (double *obj, int len)
{
	double max_n = (len > 0) ? obj[0] : std::numeric_limits<double>::max();
	for (int i=1; i < len; i++)
	{
		if (obj[i] > max_n)
			max_n = obj[i];
	}
	return max_n;
}
double max (vector<double> obj)
{
	double max_n = (obj.size() > 0) ? obj.front() : std::numeric_limits<double>::max();
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i) > max_n)
		{
			max_n = obj.at(i);
		}
	}
	return max_n;
}
idpair max (vector<idpair> obj)
{
	double max_n = (obj.size() > 0) ? obj.front().y : std::numeric_limits<double>::max();
	int x = 0;
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i).y > max_n)
		{
			max_n = obj.at(i).y;
			x = obj.at(i).x;
		}
	}
	return idpair(x, max_n);
	//return *max_element(&obj.front().y, &obj.back().y);
}
idpair max (idpair* obj, int len)
{
	double max_n = (len > 0) ? obj[0].y : std::numeric_limits<double>::max();
	int x = 0;
	for (int i=1; i < len; i++)
	{
		if (obj[i].y > max_n)
		{
			max_n = obj[i].y;
			x = i;
		}
	}
	return idpair(x, max_n);
}

/*** secondmax ***
 * Finds the second maximum value of a vector, works only for psoitive values *
 * - obj: an array of type double or a vector of type idpair (considering y then)
 * - len: the number of array elements
 * - return: the maximum value */
idpair secondmax (vector<idpair> obj)
{
	double max_n = (obj.size() > 0) ? obj.front().y : std::numeric_limits<double>::max();
	double max_n2 = -1;
	int x = 0, x2 = 0;
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i).y > max_n)
		{
			max_n = obj.at(i).y;
			x = obj.at(i).x;
		}
	}
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i).y > max_n2 && obj.at(i).y < max_n)
		{
			max_n2 = obj.at(i).y;
			x2 = obj.at(i).x;
		}
	}
	if (max_n2 == -1) 
		max_n2 = max_n;
	return idpair(x2, max_n2);
}

/*** firstIdpair ***
 * Returns the first/second/third... idpair from an idpair vector after a certain x offset *
 * - offset: the x value from which to start
 * - obj: a vector of type idpair
 * - n: specifies if to take the first (1), second (2), third (3)... idpair after the offset
 * - return: the first idpair after offset */
idpair firstIdpair (const int offset, const vector<idpair> obj, int n)
{
	int s=1;
	for (int i=0; i < obj.size(); i++)
	{
		if (obj.at(i).x > offset)
		{
			if (s < n)
				s++;
			else
				return obj.at(i);
		}
	}
	return idpair();
}

/*** isCloser ***
 * Determines if a is closer to b than c is to b *
 * - a: a floating point number
 * - b: a floating point number
 * - return: true, if a is closer to b, false otherwise */
bool isCloser(double a, double b, double c)
{
	if ( abs(a - b) < abs(b - c) )
		return true;
	else
		return false;
}

/*** timeMeasure ***
 * Starts or stops a time measurement (accuracy is one second) *
 * - start: boolean which is true when intending to start and false when intending to stop measurement
 * - return: the seconds elapsed */
int timeMeasure (bool start)
{
	static int start_time;
	if (start == true) 
	{
		start_time = time(NULL);
	}
	else
	{
		return (time(NULL) - start_time);
	}
}

/*** getMean ***
 * Calculates the mean value of a part of an array of double values *
 * - array: array of double values *
 * - start: starting point in array *
 * - len: number of array elements to be considered *
 * - return: the mean of all considered numbers */
double getMean(const double *array, const int start, const int len)
{
	double m = 0.0;
	for(int i=start; i<start+len; i++)
	{
		m += array[i];
	}
	return (m / double(len));
}

/*** corrFunc (in construction) ***
 * Calculates the correlation function between two datasets and returns the values in an array *
 * - q1: array of values of the first quantity
 * - q2: array of values of the second quantity
 * - n: number of values
 * - dur: number of time steps (total time)
 * - dt: time step size
 * - return: array of values of the correlation function */
/*double* corrFunc(double *q1, double *q2, int n, int dur, double dt)
{
	double *c = new double[dur*2+1];
	
	for(int tau=-dur; tau<=dur; tau++)
	{
		c[tau+dur] = 0.0; // [tau+dur] since we are starting the loop from -dur
		for(int t=-dur; t<=dur; t++)
		{
			c[tau+dur] += q1[t+dur] * q2[t+dur+tau] * dt;
		}
		q1[
	}
	
	return NULL;
}*/


/*** getSeparator ***
 * Determines the length of a line in console/terminal and creates a line separator *
 * - return: the separator string */
string getSeparator()
{
	/*const char sep_char = '-'; // character from which the separator is created
	string separator;
	
	struct winsize w;
   ioctl(0, TIOCGWINSZ, &w);
	
	for (int i=0; i<w.ws_col; i++)
		separator += sep_char;

	return separator;*/
	return string("");
}

/*** showChDirErrMessage ***
 * Prints the last errno message out of the possible chdir() errors via cout */
void showChDirErrMessage()
{
	cout << "Error while changing directory: ";
	switch(errno)
	{
		case EACCES:
			cout << "EACCES" << endl;
			break;
		case ENAMETOOLONG:
			cout << "ENAMETOOLONG" << endl;
			break;
		case ENOENT:
			cout << "ENOENT" << endl;
			break;
		case ENOTDIR:
			cout << "ENOTDIR" << endl;
			break;
		case ELOOP:
			cout << "ELOOP" << endl;
			break;
	}
}

/*** getClockSeed ***
 * Gets a random generator seed from the computer clock, guarantees not to return *
 * the same seed in two subsequent calls (very important!) *
 * - return: clock counts since epoch of the clock */
static unsigned int getClockSeed()
{
	static int last_seed;
	int seed;

	while ( (seed = chrono::system_clock::now().time_since_epoch().count()) == last_seed ) {}
	last_seed = seed;

	return seed;
}



#ifdef NETWORKSIMULATION
/*========================================================================================================================================
  = Below are model functions that are used to fit data =*/

/*** gauss_f ***
 * Gaussian fit function of the form A * exp(-pow(r, 2) / (2*var)) + b *
 * This function does not specify which parameters will be fitted */
/*int gauss_f (const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n = ((fitdata_pairs*)data)->n;
	double *r = ((fitdata_pairs*)data)->x;
	double *y = ((fitdata_pairs*)data)->y;
	double *sigma = ((fitdata_pairs*) data)->sigma;

	double A = gsl_vector_get(x, 0);
	double var = gsl_vector_get(x, 1);
	double b = gsl_vector_get(x, 2);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double Yi = A * exp(-pow(r[i], 2) / (2*var)) + b; // equation of the Gaussian
		gsl_vector_set(f, i, (Yi - y[i])/sigma[i]);
	}

	return GSL_SUCCESS;
} */

/*** gauss_df_A_var ***
 * Jacobian for the Gaussian fit function of the form A * exp(-r^2 / (2*var)) + b *
 * with fit parameters A, var */
/*int gauss_df_A_var (const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n = ((fitdata_pairs*)data)->n;
	double *r = ((fitdata_pairs*) data)->x;
	double *y = ((fitdata_pairs*) data)->y;
	double *sigma = ((fitdata_pairs*) data)->sigma;

	double A = gsl_vector_get(x, 0);
	double var = gsl_vector_get(x, 1);

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-r^2 / (2*var))  */
		/* and the xj are the parameters (A,var) */
		/*
		double e = exp(-pow(r[i], 2) / (2*var));
		gsl_matrix_set(J, i, 0, e/sigma[i]); // derivative wrt A
		gsl_matrix_set(J, i, 1, pow(r[i], 2) / (2*pow(var,2)) * A * e/sigma[i]); // derivative wrt var
	}
	return GSL_SUCCESS;
} */

/*** gauss_fdf_A_var ***
 * Calls Gaussian function gauss_f and Jacobian function gauss_df_A_var */
/*int gauss_fdf_A_var (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	gauss_f(x, data, f);
	gauss_df_A_var(x, data, J);

	return GSL_SUCCESS;
} */
#endif

#ifdef NEURONSIMULATION
/*** exp_f ***
 * Exponential fit function of the form A * exp(-(t-t_0) / tau) + b *
 * This function does not specify which parameters will be fitted */
/*int exp_f (const gsl_vector *x, void *data, gsl_vector *f)
{
	size_t n = ((fitdata_idpairs*) data)->n;
	idpair *val = ((fitdata_idpairs*) data)->val;
	double *sigma = ((fitdata_idpairs*) data)->sigma;

	double A = gsl_vector_get(x, 0);
	double tau = gsl_vector_get(x, 1);
	double b = gsl_vector_get(x, 2);
	double t_0 = gsl_vector_get(x, 3); // constant, not a fit parameter

	size_t i;

	for (i = 0; i < n; i++)
	{
		double t = val[i].x;
		double Yi = A * exp(-(t-t_0) / tau) + b;
		gsl_vector_set(f, i, (Yi - val[i].y)/sigma[i]);
	}

	return GSL_SUCCESS;
}

/*** exp_df_A_tau_b ***
 * Jacobian for the exponential fit function of the form A * exp(-(t-t_0) / tau) + b *
 * with fit parameters A, tau */
/*int exp_df_A_tau (const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n = ((fitdata_idpairs*) data)->n;
	idpair *val = ((fitdata_idpairs*) data)->val;
	double *sigma = ((fitdata_idpairs*) data)->sigma;

	double A = gsl_vector_get(x, 0);
	double tau = gsl_vector_get(x, 1);
	double t_0 = gsl_vector_get(x, 3); // constant, not a fit parameter

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-(t-t_0) / tau) + b  */
		/* and the xj are the parameters (A,tau,b) */
		/*double t = val[i].x;
		double s = sigma[i];
		double e = exp(-(t-t_0) / tau);
		gsl_matrix_set(J, i, 0, e/s); // derivative wrt A
		gsl_matrix_set(J, i, 1, (t-t_0) / pow(tau, 2) * A * e/s); // derivative wrt tau
	}
	return GSL_SUCCESS;
}

/*** exp_df_A_tau_b ***
 * Jacobian for the exponential fit function of the form A * exp(-(t-t_0) / tau) + b *
 * with fit parameters A, tau, b */
/*int exp_df_A_tau_b (const gsl_vector *x, void *data, gsl_matrix *J)
{
	size_t n = ((fitdata_idpairs*) data)->n;
	idpair *val = ((fitdata_idpairs*) data)->val;
	double *sigma = ((fitdata_idpairs*) data)->sigma;

	double A = gsl_vector_get(x, 0);
	double tau = gsl_vector_get(x, 1);
	double t_0 = gsl_vector_get(x, 3); // constant, not a fit parameter

	size_t i;

	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-(t-t_0) / tau) + b  */
		/* and the xj are the parameters (A,tau,b) */
		/*double t = val[i].x;
		double s = sigma[i];
		double e = exp(-(t-t_0) / tau);
		gsl_matrix_set(J, i, 0, e/s); // derivative wrt A
		gsl_matrix_set(J, i, 1, (t-t_0) / pow(tau, 2) * A * e/s); // derivative wrt tau
		gsl_matrix_set(J, i, 2, 1/s); // derivative wrt b
	}
	return GSL_SUCCESS;
}

/*** exp_fdf_A_tau ***
 * Calls exponential function exp_f and Jacobian function exp_df_A_tau */
/*int exp_fdf_A_tau (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	exp_f(x, data, f);
	exp_df_A_tau(x, data, J);

	return GSL_SUCCESS;
}

/*** exp_fdf_A_tau_b ***
 * Calls exponential function exp_f and Jacobian function exp_df_A_tau_b */
/*int exp_fdf_A_tau_b (const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
	exp_f(x, data, f);
	exp_df_A_tau_b(x, data, J);

	return GSL_SUCCESS;
}*/
#endif

