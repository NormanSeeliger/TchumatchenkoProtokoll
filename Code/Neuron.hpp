/*******************************************************
 *** LIF model of a single neuron with ChR2 channels ***
 ***       (c) Jannik Luboeinski 2015			        ***
 *******************************************************/

#include <random>

using namespace std;

#include "ChR2_3state.hpp"

#define TYPE_EXC 1
#define TYPE_INH 2
//#define INWARD_RECT // specifies if voltage-dependent inward rectification shall be considered

/*** Neuron class ***
 * Represents one neuron with a stochastically modelled population of ChR2 */
class Neuron {	

private:

/*** Computational parameters ***/
double dt; // ms, one time step for numerical simulation

/*** State variables ***/
double V; // mV, the current membrane potential
double I_ChR2; // nA, the current evoked by the ChR2 channels
double I_cst; // nA, the externally applied stimulus current
double I_fluct; // nA, the fluctuating current evoked by external synaptic inputs (computed using an OU process with mean 0)
double I_syn; // nA, the synaptic input affecting this neuron
double refractory; // ms, time span until refractory period is over
vector<int> spike_history; // vector of all spike times (in units of time steps) in the process since the last reset
int inh_incoming; // number of incoming inhibitory connections in a network
int exc_incoming; // number of incoming excitatory connections in a network
int inh_outgoing; // number of outgoing inhibitory connections in a network
int exc_outgoing; // number of outgoing excitatory connections in a network
Stimulus cst; // current stimulus for this neuron
Stimulus lst; // light stimulus for this neuron
bool cst_set; // true if a current stimulus has been set since the last reset
bool lst_set; // true if a light stimulus has been set since the last reset
ChR2 ch; // the stochastic ChR2 channel instance
minstd_rand0 rg; // default uniform generator for random numbers to generate noise (seed is chosen in constructor)
//sse_random_gen rg; // SSE uniform generator for random numbers (seed is chosen in constructor)
normal_distribution<double> norm_dist; // normal distribution for random numbers to generate noise, constructed in Neuron class constructor

protected:

/*** Physical parameters ***/
double tau_mem; // ms, the membrane time constant
double tau_ref; // ms, refractory period - should be at least one time step!
double G_mem; // µS, conductance of the cell membrane
double V_rev; // mV, the reversal potential of the neuron
double V_reset; // mV, the reset potential of the neuron
double V_th; // mV, the threshold potential of the neuron
double V_spike; // mV, the height of an action potential
int N_ChR2; // expression level/number of ChR2 channels in this neuron
#ifdef INWARD_RECT
double V_inw1; // mV, constant for voltage-dependent inward rectification
double V_inw2; // mV, constant for voltage-dependent inward rectification
#endif
double tau_OU; // ms, correlation time of the Ornstein-Uhlenbeck process
double sigma_WN; // nA, standard deviation of the Gaussian white noise inside the OU process
double I_const; // nA, the constant current evoked by external synaptic inputs 
int type; // the type of this neuron (inhibitory/excitatory)


/*** normalRandomNumber ***
 * Returns a random number drawn from a normal distribution with standard deviation 1 and mean 0 *
 * - return: the random number of type double (technically in units of sqrt(s)) */
double normalRandomNumber()
{
	double nd = norm_dist(rg);

	/* Check for mean and standard deviation *
	static int i = 0;
	static double mean =0.0;
	static double var =0.0;
	const int av_steps = 1000000;
	i++;
	mean += nd;
	var += pow(nd, 2); // compute variance based on an assumed mean of 0.0
	if (i == av_steps)
	{
		cout << endl << "mean: " << mean / double(av_steps) << endl;
		cout << "stddev: " << sqrt(var / double(av_steps)) << endl;
		mean = 0.0;
		var = 0.0;
		i=0;
	} */
	return nd;
}

public:

/*** saveNeuronParams ***
 * Saves all the neuron parameters (including the channel parameters) to a given file */
void saveNeuronParams(ofstream *f) const
{
	*f << endl;
	*f << "Neuron parameters:" << endl;
	*f << "tau_mem = " << tau_mem << " ms" << endl;
	*f << "tau_ref = " << tau_ref << " ms" << endl;
	*f << "G_mem = " << G_mem << " µS" << endl;
	*f << "V_rev = " << V_rev << " mV" << endl;
	*f << "V_reset = " << V_reset << " mV" << endl;
	*f << "V_th = " << V_th << " mV" << endl;
	*f << "V_spike = " << V_spike << " mV" << endl;
	*f << "N_ChR2 = " << N_ChR2 << endl;
#ifdef INWARD_RECT
	*f << "V_inw1 = " << V_inw1 << " mV" << endl;
	*f << "V_inw2 = " << V_inw2 << " mV" << endl;
#endif
	*f << "tau_OU = " << tau_OU << " ms" << endl;
	*f << "sigma_WN = " << sigma_WN << " nA" << endl;
	*f << "I_const = " << I_const << " nA" << endl;

	ch.saveChannelParams(f);
}

/*** getNumberIncoming ***
 * Returns the number of either inhibitory or excitatory incoming connections to this neuron *
 * from other network neurons *
 * - int type: the type of incoming connections (inh./exc.)
 * - return: the number of incoming connections */
int getNumberIncoming(int type) const
{
	if (type == TYPE_INH)
		return inh_incoming;
	else if (type == TYPE_EXC)
		return exc_incoming;
}
int getNumberOutgoing(int type) const
{
	if (type == TYPE_INH)
		return inh_outgoing;
	else if (type == TYPE_EXC)
		return exc_outgoing;
}
/*** incNumberIncomingInh ***
 * Increases the number of incoming inhibitory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberIncomingInh()
{
	inh_incoming++;
}

/*** incNumberIncomingExc ***
 * Increases the number of incoming excitatory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberIncomingExc()
{
	exc_incoming++;
}


/*** incNumberOutgoingInh ***
 * Increases the number of incoming inhibitory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberOutgoingInh()
{
	inh_outgoing++;
}

/*** incNumberOutGoingExc ***
 * Increases the number of incoming excitatory connections to this neuron (only to be used while *
 * a network is being built) */
void incNumberOutgoingExc()
{
	exc_outgoing++;
}

/*** getVoltage ***
 * Returns the membrane potential of the neuron *
 * - return: the membrane potential in mV */
double getVoltage() const
{
	return V;
}

/*** getCurrent ***
 * Returns total current effecting the neuron *
 * - return: the instantaneous current in nA */
double getCurrent() const
{
	return I_ChR2+I_cst+I_const+I_fluct;
}

/*** getChR2Current ***
 * Returns current evoked by the ChR2 molecules *
 * - return: the instantaneous ChR2 channel current in nA */
double getChR2Current() const
{
	return I_ChR2;
}

/*** getStimulusCurrent ***
 * Returns current evoked by external stimulation *
 * - return: the instantaneous current stimulus in nA */
double getStimulusCurrent() const
{
	return I_cst;
}

/*** getFluctCurrent ***
 * Returns fluctuating current evoked by external synaptic inputs *
 * - return: the instantaneous fluctuating current in nA */
double getFluctCurrent() const
{
	return I_fluct;
}

/*** getConstCurrent ***
 * Returns the constant current elicited by the surrounding network *
 * - return: the constant current in nA */
double getConstCurrent() const
{
	return I_const;
}

/*** getSynapticCurrent ***
 * Returns the synaptic current that arrived in the previous time step *
 * - return: the synaptic current in nA */
double getSynapticCurrent() const
{
	return I_syn;
}

/*** getActivity ***
 * Returns true if the neuron is spiking in this instant of duration dt *
 * - return: whether neuron is firing or not */
bool getActivity() const
{
	if (V == V_spike)
		return true;
	else
		return false;
}

/*** getSpikeTime ***
 * Returns the spike time for a given spike number *
 * - int n: the number of the spike (in temporal order, starting with 1)
 * - return: the spike time in units of time bins for the n-th spike (or -1 if it does not exist) */
int getSpikeTime(int n) const
{
	if (n <= spike_history.size() && n >= 1)
		return spike_history.at(n-1);
	else
		return -1;
}

/*** hideSpikeTime ***
 * "Hides" the spike time for a given spike number, which means that *
 * it will be stored as a negative number, such that it will not anymore *
 * be taken into consideration for synaptic transmission *
 * - int n: the number of the spike (in temporal order, starting with 1)
 * - return: the spike time for the n-th spike (or -1.0 if it does not exist) */
void hideSpikeTime(int n)
{
	if (n <= spike_history.size() && n >= 1)
		spike_history.at(n-1) *= -1;
}

/*** getSpikeCount ***
 * Returns the number of spikes that have occurred since the last reset *
 * - return: the number of spikes */
int getSpikeCount() const
{
	return spike_history.size();
}

/*** getSteadyVoltage ***
 * Returns the steady-state voltage (V_ss) for a constant light stimulus of given irradiance *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - return: the steady-state voltage at this irradiance in mV */
double getSteadyVoltage(double irradiance) const
{
	double ph_inf, gamma_r;
	ph_inf = ch.getPhotonInflux(irradiance);
	gamma_r = ch.getGammaR();

	return (1.0/(14.0*ch.gamma_d0*(ph_inf + gamma_r)*G_mem))*(1250.0*ph_inf*gamma_r*(G_mem + ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(ph_inf + gamma_r)*
			 (7.0*I_const + G_mem*(760.0 + 7.0*V_rev)) - sqrt(50000.0*ph_inf*gamma_r*(76.0*ph_inf*ch.gamma_d0 + 125.0*ph_inf*gamma_r + 76.0*ch.gamma_d0*gamma_r)*
			 ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(ch.gamma_d0*gamma_r*(7.0*I_const + G_mem*(-760+7*V_rev)) + ph_inf*(-1250*gamma_r*(G_mem - 
			 ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(7.0*I_const + G_mem*(-760+7.0*V_rev))), 2))); // mind the 1e-9!
}

/*** getSteadyCurrent ***
 * Returns the steady-state current (I_ss) for a constant light stimulus of given irradiance *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - return: the steady-state current at this irradiance in nA */
double getSteadyCurrent(double irradiance) const
{
	double ph_inf, gamma_r;
	ph_inf = ch.getPhotonInflux(irradiance);
	gamma_r = ch.getGammaR();
	
	return -((1.0/(14.0*ch.gamma_d0*(ph_inf + gamma_r)))*(-1250.0*ph_inf*gamma_r*(G_mem + ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(ph_inf + gamma_r)*
			   (7.0*I_const + G_mem*(-760.0 + 7.0*V_rev)) + sqrt(50000.0*ph_inf*gamma_r*(76.0*ph_inf*ch.gamma_d0 + 125.0*ph_inf*gamma_r + 76.0*ch.gamma_d0*gamma_r)*
			   ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(ch.gamma_d0*gamma_r*(7.0*I_const + G_mem*(-760+7*V_rev)) + ph_inf*(-1250*gamma_r*(G_mem - 
			   ch.G_ChR2*1e-9*N_ChR2) + ch.gamma_d0*(7.0*I_const + G_mem*(-760+7.0*V_rev))), 2)))); // mind the 1e-9!
}

/*** getSaturationCurrent ***
 * Returns the saturation current for a constant light stimulus in the infinite irradiance limit *
 * - return: the saturation current for the set parameters in nA */
double getSaturationCurrent() const
{
	double gamma_r = ch.getGammaR();
	return -((1.0/(14.0*ch.gamma_d0))*(-760.0*ch.gamma_d0*G_mem - 1250.0*gamma_r*G_mem + 7.0*ch.gamma_d0*I_const - 1250.0*gamma_r*ch.G_ChR2*1e-9*N_ChR2 + 
   		   7.0*ch.gamma_d0*G_mem*V_rev + sqrt(3800000.0*ch.gamma_d0*gamma_r*ch.G_ChR2*1e-9*G_mem*N_ChR2 + 
     		   6250000.0*pow(gamma_r,2)*ch.G_ChR2*1e-9*G_mem*N_ChR2 + pow(-1250.0*gamma_r*(G_mem - ch.G_ChR2*1e-9*N_ChR2) + 
       	   ch.gamma_d0*(7.0*I_const + G_mem*(-760.0 + 7.0*V_rev)),2)))); // mind the 1e-9!
}

/*** getOpenProb ***
 * Returns the probability that a ChR2 channel is open *
 * - return: the open probability */
double getOpenProb() const
{
	return ch.getO();
}

/*** processTimeStep ***
 * Processes one time step (of duration dt) for the neuron * 
 * - int t_step: time step at which to evaluate stimulus (< 0 before stimulus onset) *
 * - double _I_syn [optional]: synaptic current evoked by other neurons, only of importance in a network */
void processTimeStep(int t_step, double _I_syn = 0.0)
{
	double sigma_OU = sigma_WN / sqrt(2.*tau_OU/1000.); // leads to the sigma of the OU process; times 1000 because required unit is [s]
#ifdef INWARD_RECT
	double G = (1-exp(-V/V_inw1))/(V/V_inw2); // from Grossman et al., 2011
#else
	double G = 1.;
#endif

	if (lst_set)
		ch.setIrradiance(lst.getStimulusAt(t_step)); // set the current light stimulus value as irradiance for the channel
	ch.processTimeStep(t_step, V); // process time step for channels (is considering the set light stimulus)
	if (cst_set)
		I_cst = cst.getStimulusAt(t_step); // get stimulus current in nA
	if (V < 0.) // no current during action potential
		I_ChR2 = - N_ChR2 * 1e-9 * ch.G_ChR2 * G * V * ch.getO(); // compute channel current in nA (therefore, multiply G_ChR2 by 10^-9)
	else
		I_ChR2 = - N_ChR2 * 1e-9 * ch.G_ChR2 * G * V_reset * ch.getO(); // use V_reset as voltage for time bin in which a spike occurred

#ifndef LIF_FR_TEST
#ifdef DELTA_SYNAPSES
	I_fluct = normalRandomNumber() * sigma_WN; // units wrong??? sqrt(s)*nA???
#else
	I_fluct = I_fluct * exp(-dt/tau_OU) + normalRandomNumber() * sqrt(1. - exp(-2.*dt/tau_OU)) * sigma_OU; // compute fluctuating synaptic input in nA
#endif
#endif
	I_syn = _I_syn;

	//V += dt/tau_mem * (- V + V_rev + (I_const + I_fluct + I_cst + I_ChR2 + I_syn) / G_mem); // compute mem. pot. in mV (Euler method)
	
	V = V * exp(-dt/tau_mem) + (V_rev + (I_const + I_fluct + I_cst + I_ChR2 + I_syn) / G_mem) * (1. - exp(-dt/tau_mem)); // compute mem. pot. in mV (analytical solution)

	if (refractory > 0.0) // in refractory period
	{
		V = V_reset;
		refractory -= dt;
	}
	else if (V >= V_th) // threshold crossing
	{
		V = V_spike;
		refractory = tau_ref;
		if (t_step >= 0) // count only those spikes that occur after the stimulus onset
			spike_history.push_back(t_step);
	}
}

/*** setCurrentStimulus ***
 * Sets a current stimulus for the neuron *
 * - Stimulus& _cst: shape of one stimulus period */
void setCurrentStimulus(const Stimulus& _cst)
{
	cst = Stimulus(_cst); // use copy constructor -> TODO: why use copy constructor here instead of assignment operator?
	cst_set = true;
}

/*** setLightStimulus ***
 * Sets a light stimulus for the ChR2 molecules in the neuron *
 * - Stimulus& _lst: shape of one stimulus period */
void setLightStimulus(const Stimulus& _lst)
{
	lst = _lst; // use assignment operator
	lst_set = true;
}

/*** getCurrentStimulusAt ***
 * Returns the current stimulus magnitude at a given time *
 * - t_step: time step at which to evaluate stimulus
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getCurrentStimulusAt(int t_step) const
{
	if (cst_set)
		return cst.getStimulusAt(t_step);
	else
		return 0.0;
}

/*** getLightStimulusAt ***
 * Returns the light stimulus magnitude at a given time *
 * - t_step: time step at which to evaluate stimulus
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getLightStimulusAt(int t_step) const
{
	if (lst_set)
		return lst.getStimulusAt(t_step);
	else
		return 0.0;
}

/*** multiplyLightStimulus ***
 * Multiplies the set light stimulus by a real number *
 * - double r: number to multiply */
void multiplyLightStimulus(double r)
{
	lst.multiplyBy(r);
}

/*** multiplyCurrentStimulus ***
 * Multiplies the set current stimulus by a real number *
 * - double r: number to multiply */
void multiplyCurrentStimulus(double r)
{
	cst.multiplyBy(r);
}

/*** setConstCurrent ***
 * Sets the constant current to a newly defined value *
 * - double _I_const: constant current in nA */
void setConstCurrent(double _I_const)
{
	I_const = _I_const;
}

/*** setTauOU ***
 * Sets the time constant of the Ornstein-Uhlenbeck process *
 * (synaptic time constant of assumed input synapses) *
 * - int _type: the neuron type */
void setTauOU(double _tau_OU)
{
	tau_OU = _tau_OU;
}

/*** setType ***
 * Sets the type of this neuron *
 * - int _type: the neuron type */
void setType(int _type)
{
	type = _type;
}

/*** getType ***
 * Returns the type of this neuron *
 * - return: the neuron type */
int getType()
{
	return type;
}

/*** reset ***
 * Resets neuron to initial state */
void reset()
{
	vector<int>().swap(spike_history); // additionally to clearing the vector, reallocates it
	inh_incoming = 0;
	exc_incoming = 0;
	inh_outgoing = 0;
	exc_outgoing = 0;
	V = V_rev;
	I_ChR2 = 0.0; // assume closed channel
	I_cst = 0.0;
	I_fluct = 0.0;
	I_syn = 0.0;
	refractory = 0.0;
	cst_set = false;
	lst_set = false;
	
	rg.seed(getClockSeed()); // set new seed by clock's epoch
	norm_dist.reset(); // reset the normal distribution for random numbers
	ch.reset();
}

/*** Empty constructor ***
 * Just for syntactic reasons - when creating an array of Neuron class */
/*Neuron() : dt(0.0), ch(0.0), norm_dist(0.0,1.0), rg(getClockSeed())
{
}*/

/*** Constructor ***
 * Sets all parameters on experimentally determined values */
Neuron(const double _dt) : 
	dt(_dt), ch(_dt), norm_dist(0.0,1.0), rg(getClockSeed())
{
	//tau_mem = 30.0; // from Bachelor's thesis
	tau_mem = 10.0; // estimated
	tau_ref = 3.0; // from Bachelor's thesis
	//G_mem = 0.2; // from Bachelor's thesis
	G_mem = 0.1; // estimated
	V_rev = -65.0; // from Bachelor's thesis
	V_reset = -70.0; // from Bachelor's thesis 
	//V_th = -50.0; // from Bachelor's thesis 
	V_th = -55.0; // estimated
	V_spike = 35.0; // estimated from sketch in Bachelor's thesis
	//N_ChR2 = 120000; // from Grossman et al., 2011 (cf. calculation)
	//N_ChR2 = 50000; // minimum from Michael's BSc thesis
	//N_ChR2 = 625000; // maximum from Michael's BSc thesis
	//N_ChR2 = 2e9; // from Nagel et al., 2003 for oocyte
	//N_ChR2 = 400000; // maximum from Michael's BSc thesis
	N_ChR2 = 60000; // low estimation
	//N_ChR2 = 300000; // high estimation
#ifdef INWARD_RECT
	V_inw1 = 40.; // from Grossman et al., 2011
	V_inw2 = 15.; // from Grossman et al., 2011
#endif
	//tau_OU = 7.1; // from Bachelor's thesis
	//sigma_WN = 0.003; // estimated, "small" sigma
	//sigma_WN = 0.095; // estimated, "large" sigma
	sigma_WN = 0.010; // estimated, "medium" sigma
	I_const = 0.116; // estimated (just before steady-state is above threshold), is overwritten by setConstCurrent()
	setTauOU(5.0); // estimated

	
	reset();
}

/*** Destructor ***
 * Frees the allocated memory */
~Neuron()
{
	//cout << "Destructor ~Neuron() was called!" << endl;
}


};


















