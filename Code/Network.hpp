/********************************************************************
 *** Model of a network, consisting of neurons with ChR2 channels ***
 ***       		  (c) Jannik Luboeinski 2015/2016			       	***
 ********************************************************************/

#include <random>

using namespace std;

//#define DELTA_SYNAPSES // if defined, delta synapses (synaptic time constant -> 0) are are used, else, exponential decay is used

#include "Neuron.hpp"
#include "GaussianNoise.cpp"

/*** Network class ***
 * Represents a network of neurons */
class Network {	

private:

/*** Computational parameters ***/
double  dt;				// ms, one time step for numerical simulation
//double I_neglect; 	// nA, threshold current below which PSCs for one neuron will be neglected
//double t_neglect_exc; // tau_syn, time after which PSCs will be neglected in exc. neurons (in units of tau_syn)
//double t_neglect_inh; // tau_syn, time after which PSCs will be neglected in inh. neurons (in units of tau_syn)
//double expDecaySyn; 	// pre-computed value for better performance
//int tb_max_exc; 		// negligence time (including tau_syn) for exc. neurons in time bins - for faster computation
//int tb_max_inh; 		// negligence time (including tau_syn) for inh. neurons in time bins - for faster computation
int 	N; 				// total number of excitatory plus inhibitory neurons
double 	sigma_in2ex; 	// Sigma for gauss distribution from inhibitory to  excitatory neurons
double 	mu_in2ex;   	// Mu for gauss distribution from inhibitory to 	excitatory neurons
double 	sigma_in2in;	// Sigma for gauss distribution from inhibitory to  inhibitory neurons
double 	mu_in2in;  		// Mu for gauss distribution from inhibitory to 	inhibitory neurons

double 	sigma_ex2ex; 	// Sigma for gauss distribution from excitatory to  excitatory neurons
double 	mu_ex2ex;   	// Mu for gauss distribution from excitatory to 	excitatory neurons
double 	sigma_ex2in;	// Sigma for gauss distribution from excitatory to  inhibitory neurons
double 	mu_ex2in;  		// Mu for gauss distribution from excitatory to 	inhibitory neurons

double amp_ex2ex;
double amp_ex2in;

/*** State variables ***/
vector<Neuron> neurons; // vector of all N neuron instances (first excitatory, then inhibitory)
bool**  w; 				// the connectivity matrix, the main diagonal is zero (because there is no self-coupling)
minstd_rand0 rg; 		// default uniform generator for random numbers to establish connections (seed is chosen in constructor)
uniform_real_distribution<double> u_dist;
//sse_random_gen rg; 	// SSE uniform generator for random numbers (seed is chosen in constructor)

protected:

/*** Physical parameters ***/
int  	Nl; 			// number of neurons in one line (row or column) of the exc. population (better choose an odd number, for there exists a "central" neuron)
int 	Nl_inh; 		// number of neurons in one line (row or column) of the inh. population (better choose an odd number, for there exists a "central" neuron)
double 	tau_syn; 		// ms, the synaptic time constant
double 	pc; 			// connection probability (prob. that a directed connection exists)
double 	J_ee; 			// pC, magnitude of excitatory PSP effecting an excitatory postsynaptic neuron
double 	J_ei; 			// pC, magnitude of excitatory PSP effecting an inhibitory postsynaptic neuron
double 	J_ie; 			// pC, magnitude of inhibitory PSP effecting an excitatory postsynaptic neuron
double 	J_ii; 			// pC, magnitude of inhibitory PSP effecting an inhibitory postsynaptic neuron
double 	sigma_light; 	// standard deviation of light distribution (dimensionless)

public:

/*** cNN (macro) ***
 * Returns a consecutive number for excitatory neuron (i|j) rather than a pair of numbers like (i|j), be  *
 * aware that it starts with zero, unlike i and j *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
#define cNN(i, j) (((i)-1)*Nl + ((j)-1))

/*** cInh (macro) ***
 * Returns a consecutive number for inhibitory neuron (i|j) rather than a pair of numbers like (i|j), be  *
 * aware that it starts with zero, unlike i and j *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
#define cInh(i, j) (((i)-1)*Nl_inh + ((j)-1))

#define cInhG(i, j) ((((i)-1)*Nl_inh + ((j)-1))+Nl*Nl)

/*** rowNN (macro) ***
 * Returns the row number for excitatory neuron n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
#define rowNN(n) ((((n) - ((n) % Nl)) / Nl) + 1)

/*** rowInh (macro) ***
 * Returns the row number for inhibitory neuron n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
// WAS_NN (15.04.2016) : #define rowInh(n) ((((n) - ((n) % Nl_inh)) / Nl_inh) + 1)
#define rowInh(n) (2*(((n-(Nl*Nl) - ((n-(Nl*Nl)) % Nl_inh)) / Nl_inh) + 1)-1)

/*** colNN (macro) ***
 * Returns the column number for excitatory neuron n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
#define colNN(n) (((n) % Nl) + 1) 

/*** colInh (macro) ***
 * Returns the column number for inhibitory neuron n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
 // WAS_NN (15.04.2016) : #define rowInh (((n) % N1_inh) + 1)
#define colInh(n) (2*(((n-(Nl*Nl)) % Nl_inh) + 1)-1) 

/*** shallBeConnected ***
 * Draws a uniformly distributed random number from the interval 0.0 to 1.0 and returns, *
 * depending on the connection probability, whether or not a connection shall be established *
 * - return: true if connection shall be established, false if not */
bool shallBeConnectedEx(int m,int n)
/* Excitatory to ecxitatory, medium range*/
{
	double distance = sqrt(pow(double(rowNN(m))-double(rowNN(n)),2)+pow(double(colNN(m))-double(colNN(n)),2));
	if (u_dist(rg) <= ((amp_ex2ex*exp((-1*(pow(distance-mu_ex2ex, 2)))/(2*pow(sigma_ex2ex,2))))))
		return true;
	else
		return false; 
} 

bool shallBeConnectedEx2(int m,int n)
/* Excitatory to inhibitory neurons, large range*/
{
	double distance = sqrt(pow(double(rowInh(n))-double(rowNN(m))+0.5,2)+pow(double(colInh(n))-double(colNN(m))+0.5,2));
	// Reversed computation - m is excitatory and gets subtracted
	if (u_dist(rg) <= ((amp_ex2in*exp((-1*(pow(distance-mu_ex2in, 2)))/(2*pow(sigma_ex2in,2))))))
		return true;
	else
		return false; 
} 

bool shallBeConnectedInh(int m, int n)
/* Inhibitory to excitatory, medium range*/
{ 
//	if (generateGaussianNoise(1/(sqrt(2*3.14159265*pow(3,2))),3) <= (1/(sqrt(2*3.14159265*pow(3,2))))*exp((-1*(pow(sqrt(pow(rowInh(m)-rowNN(n)+0.5,2)+pow(colInh(m)-colNN(n)+0.5,2)),2)))/(2*pow(3,2)))) // -(1/sqrt(2))
	double distance = sqrt(pow(double(rowInh(m))-double(rowNN(n))+0.5,2)+pow(double(colInh(m))-double(colNN(n))+0.5,2));
	if (u_dist(rg) <= ((1.0/1.0)*exp((-1*(pow(distance-mu_in2ex, 2)))/(2*pow(sigma_in2ex,2))))) // -(1/sqrt(2)) // u_dist(rg) // enerateGaussianNoise(1/sqrt(2),5)
		return true;
	else 
		return false; 
}
bool shallBeConnectedInh2(int m, int n)
/* Inhibitory to inhibitory, large range */
{
//	if (GaussianNoise(0.70,2) <= (1/sqrt(pow(rowInh(m)-rowInh(n),2)+pow(colInh(m)-colInh(n),2)))) 
	if (u_dist(rg) <= (10.0/10.0)*exp((-1*(pow(sqrt(pow(double(rowInh(m))-double(rowInh(n)),2)+pow(double(colInh(m))-double(colInh(n)),2))-mu_in2in,2)))/(2*pow(sigma_in2in,2)))) // -(1/sqrt(2))
		return true;
	else
		return false; 
}


/*** areConnected ***
 * Returns whether or not there is a synapse from neuron m to neuron n *
 * - int m: the number of the first neuron in consecutive order *
 * - int n: the number of the second neuron in consecutive order *
 * - return: true if connection from m to n exists, false if not */
bool areConnected(int m, int n) const
{
	if (w[m][n])
		return true;
	else
		return false;
}

/*** saveNetworkParams ***
 * Saves all the network parameters (including the neuron and channel parameters) to a given file */
void saveNetworkParams(ofstream *f) const
{
	*f << endl;
	*f << "Network parameters:" << endl;
	*f << "Ne = " << pow(Nl, 2) << " (" << Nl << " x " << Nl << ")" << endl;
	*f << "Ni = " << pow(Nl_inh, 2) << " (" << Nl_inh << " x " << Nl_inh << ")" << endl;
	*f << "tau_syn = "

#ifdef DELTA_SYNAPSES
		<< 0
#else
	 	<< tau_syn 
#endif
		<< " ms" << endl;
	*f << "pc = " << pc << endl;
	*f << "J_ee = " << J_ee << " pC" << endl;
	*f << "J_ei = " << J_ei << " pC" << endl;
	*f << "J_ie = " << J_ie << " pC" << endl;
	*f << "J_ii = " << J_ii << " pC" << endl;
	*f << "sigma_light = " << sigma_light << endl;
	*f << "sigma_in2ex = " << sigma_in2ex << endl;
	*f << "sigma_in2in = " << sigma_in2in << endl; 
	*f << "mu_in2ex = " << mu_in2ex << endl;
	*f << "mu_in2in = " << mu_in2in << endl;
	*f << "sigma_ex2ex = " << sigma_ex2ex << endl;
	*f << "sigma_ex2in = " << sigma_ex2in << endl; 
	*f << "mu_ex2ex = " << mu_ex2ex << endl;
	*f << "mu_ex2in = " << mu_ex2in << endl;
	*f << "amp_ex2ex = " << amp_ex2ex << endl;
	*f << "amp_ex2in = " << amp_ex2in << endl;
	// TODO: 04.05.2016 SIGMA and MU INTO CONSTRUCTOR + REPLACE
	//*f << "I_neglect = " << I_neglect << " nA" << endl;
	//*f << "t_neglect_exc = " << t_neglect_exc << " tau_syn" << endl;
	//*f << "t_neglect_inh = " << t_neglect_inh << " tau_syn" << endl;
	
	neurons[0].saveNeuronParams(f); // all neurons have the same parameters, take the first one
}


/*** processTimeStep ***
 * Processes one time step (of duration dt) for the network * 
 * - int tb_step: current time step (for evaluating stimulus and for computing spike contributions) *
 * - return: number of spikes that occurred within the considered time step in the whole network */
int processTimeStep(int tb_step)
{
	double I_syn;
	int spike_count = 0;

	for (int m=0; m<N; m++)
	{
		// compute synaptic current affecting neuron (i|j)
#ifdef DELTA_SYNAPSES
		I_syn = 0.;
#else
		I_syn = neurons[m].getSynapticCurrent() * exp(-dt/tau_syn); // exponential decay of previous synaptic current contributions
#endif

		for (int n=0; n<N; n++) // loop over neurons (in consecutive order)
		{
			if (w[n][m]) // check if there is a connection n->m
			{
				if (neurons[n].getSpikeTime(neurons[n].getSpikeCount()) == (tb_step-1) && tb_step > 0) // did neuron n spike in the previous time step? tb_step > 0 crucial
				{
#ifdef DELTA_SYNAPSES
					if (neurons[n].getType() == TYPE_EXC && neurons[m].getType() == TYPE_EXC) // E -> E
						I_syn += J_ee; // Heaviside function can be left out
					else if (neurons[n].getType() == TYPE_EXC && neurons[m].getType() == TYPE_INH) // E -> I
						I_syn += J_ei;
					else if (neurons[n].getType() == TYPE_INH && neurons[m].getType() == TYPE_EXC) // I -> E
						I_syn -= J_ie;
					else // I -> I
						I_syn -= J_ii;
#else
					if (neurons[n].getType() == TYPE_EXC && neurons[m].getType() == TYPE_EXC) // E -> E
						I_syn += J_ee/tau_syn; // Heaviside function can be left out
					else if (neurons[n].getType() == TYPE_EXC && neurons[m].getType() == TYPE_INH) // E -> I
						I_syn += J_ei/tau_syn;
					else if (neurons[n].getType() == TYPE_INH && neurons[m].getType() == TYPE_EXC) // I -> E
						I_syn -= J_ie/tau_syn;
					else // I -> I
						I_syn -= J_ii/tau_syn;
#endif
				}
			}
		}

		// process time step for neuron (i|j)
		neurons[m].processTimeStep(tb_step, I_syn);
		
		// count spikes in this time step
		spike_count += int(neurons[m].getActivity());
	}
	return spike_count;
}

/*** setGaussianLightStimulus ***
 * Sets a light stimulus with its irradiance amplitude centered on neuron (_i|_j) and with *
 * Gaussian decreasing irradiance for the surrounding neurons (only affecting exc. population) *
 * - Stimulus& _lst: shape of one stimulus period *
 * - double _x: the row where the center neuron is located (can be a non-integer number) *
 * - double _y: the column where the center neuron is located (can be a non-integer number) */
void setGaussianLightStimulus(Stimulus& _lst, double _x, double _y)
{
	for (int m=0; m<Nl*Nl; m++)
	{
		neurons[m].setLightStimulus(_lst); // set temporal course of irradiance E
		neurons[m].multiplyLightStimulus( exp( -(pow(rowNN(m)-_x, 2) + pow(colNN(m)-_y, 2)) / (2*pow(sigma_light, 2)) ) ); // set spatial course (Gaussian)
	}
}

/*** setUniformLightStimulus ***
 * Sets the same light stimulus for all excitatory neurons *
 * - Stimulus& _lst: shape of one stimulus period */
void setUniformLightStimulus(Stimulus& _lst)
{
	for (int m=0; m<Nl*Nl; m++)
		neurons[m].setLightStimulus(_lst); // set temporal course of irradiance E
}

/*** setGaussianCurrentStimulus ***
 * Sets a current stimulus with its amplitude centered on neuron (_i|_j) and with *
 * Gaussian decreasing irradiance for the surrounding neurons (only affecting exc. population) *
 * - Stimulus& _cst: shape of one stimulus period *
 * - double _x: the row where the center neuron is located (can be a non-integer number) *
 * - double _y: the column where the center neuron is located (can be a non-integer number) */
void setGaussianCurrentStimulus(Stimulus& _cst, double _x, double _y)
{
	for (int m=0; m<Nl*Nl; m++)
	{
		neurons[m].setCurrentStimulus(_cst); // set temporal course of current stimulus
		neurons[m].multiplyCurrentStimulus( exp( -(pow(rowNN(m)-_x, 2) + pow(colNN(m)-_y, 2)) / (2*pow(sigma_light, 2)) ) ); // set spatial course (Gaussian)
	}
}

/*** setUniformCurrentStimulus ***
 * Sets the same current stimulus for all excitatory neurons *
 * - Stimulus& _cst: shape of one stimulus period */
void setUniformCurrentStimulus(Stimulus& _cst)
{
	for (int m=0; m<Nl*Nl; m++)
		neurons[m].setCurrentStimulus(_cst); // set temporal course of current stimulus
}

/*** setConstCurrent ***
 * Sets the constant current for all neurons to a newly defined value
 * - double _I_const: constant current in nA */
void setConstCurrent(double _I_const)
{
	for (int m=0; m<N; m++)
	{
		neurons[m].setConstCurrent(_I_const);
	}
}

/*** setSynTimeConstant ***
 * Sets the synaptic time constant *
 * - double _tau_syn: synaptic time constant in ms */
void setSynTimeConstant(double _tau_syn)
{
	tau_syn = _tau_syn;
	for (int m=0; m<N; m++)
	{
		neurons[m].setTauOU(tau_syn);
	}
}

/*** setCouplingStrengths ***
 * Sets the synaptic coupling strengths *
 * - double _J_ee: coupling strength for exc. -> exc. connections in pC *
 * - double _J_ei: coupling strength for exc. -> inh. connections in pC *
 * - double _J_ie: coupling strength for inh. -> exc. connections in pC *
 * - double _J_ii: coupling strength for inh. -> inh. connections in pC */
void setCouplingStrengths(double _J_ee, double _J_ei, double _J_ie, double _J_ii)
{
	J_ee = _J_ee;
	J_ei = _J_ei;
	J_ie = _J_ie;
	J_ii = _J_ii;
}

/*** setPC ***
 * Sets the connection probability, requires to change all the synaptic connections *
 * - double _pc: connection probability */
void setPC(double _pc)
{
	// set connection probability
	pc = _pc; 

	// compute PSC negligence time according to "pN" approximation
	//t_neglect_exc = log(pc*(Nl*Nl*J_ee+Nl_inh*Nl_inh*J_ie)/(I_neglect*tau_syn)); // for excitatory neurons
	//t_neglect_inh = log(pc*(Nl*Nl*J_ei+Nl_inh*Nl_inh*J_ii)/(I_neglect*tau_syn)); // for inhibitory neurons

	//tb_max_exc = int(round(tau_syn * t_neglect_exc / dt)); // number of time bins until negligence
	//tb_max_inh = int(round(tau_syn * t_neglect_inh / dt)); // number of time bins until negligence

	// renew synaptic connections
	for (int m=0; m<N; m++)
	{
		for (int n=0; n<N; n++)
		{
			w[m][n] = false; // necessary for resetting the connections

			if (m != n) // if not on main diagsonal (which should remain zero)
			{
			// use random generator depending on connection probability
			// WAS_NN: if (w[m][n] = shallBeConnected()) 
				if (m < Nl*Nl) { // excitatory input
					if (n >= Nl*Nl) {
						if (w[m][n] = shallBeConnectedEx2(m,n)) {// is there a connection m->n?
							neurons[n].incNumberIncomingExc();
							neurons[m].incNumberOutgoingExc();
						}
					}
					else {
						if (w[m][n] = shallBeConnectedEx(m,n)) {// is there a connection m->n?
							neurons[n].incNumberIncomingExc();
							neurons[m].incNumberOutgoingExc();
						}	
					}
				}
				else { // inhibitory input
					if (n >= Nl*Nl) { // Inhibitory PV possess gap junctions as connections
						if (w[m][n] = shallBeConnectedInh2(m,n)){ // is there a connection m->n? /TODO 14.04
							neurons[n].incNumberIncomingInh();
							neurons[m].incNumberOutgoingInh();
						}
					}
					else {
						if (w[m][n] = shallBeConnectedInh(m,n)) {
							neurons[n].incNumberIncomingInh();
							neurons[m].incNumberOutgoingInh();
						}
					}
				}
			}		
		}
	}
}

/*** getGaussSigma ***
 * Returns the standard deviation underlying the the Gaussian distribution of stimuli */
double getGaussSigma()
{
	return sigma_light;
}

/*** reset ***
 * Resets the network and all neurons to initial state */
void reset()
{
	rg.seed(getClockSeed()); // set new seed by clock's epoch
	u_dist.reset(); // reset the uniform distribution for random numbers
	
	for (int m=0; m<N; m++)
	{
		neurons[m].reset();
	}
}

/*** setConnection
	Setting parameters for inhibitory connections according to a gauss distribution with variable sigma and mu*/
//void setConnection() {

//}
/*** Constructor ***
 * Sets all parameters on experimentally determined values, creates neurons and synapses *
 * -> it is required to call setPC, setSynTimeConstant and setCouplingStrengths immediately *
 *    after calling this constructor! *
 * - double _dt: the length of one time step in ms *
 * - int _Nl: the number of neurons in one line in excitatory population (row/column) *
 * - int _Nl_inh: the number of neurons in one line in inhibitory population (row/column) - line structure so that stimulation of inhib.
						population could be implemented more easily */
Network(const double _dt, const int _Nl, const int _Nl_inh) 
	: dt(_dt), Nl(_Nl), Nl_inh(_Nl_inh), u_dist(0.0,1.0), rg(getClockSeed())
{
	sigma_light = int(round(Nl/5)); // estimated
	//I_neglect = 1e-10; // estimated
	N = Nl*Nl + Nl_inh*Nl_inh; // total number of neurons
	
	//3,5,1/sqrt(2),2,12,12,0,0,0.6,0.6
	sigma_in2ex = 2;
	sigma_in2in = 2;
	mu_in2ex = 0;
	mu_in2in = 0;
	sigma_ex2ex = 4;
	sigma_ex2in = 4;
	mu_ex2ex = 0;
	mu_ex2in = 0;
	amp_ex2in = 0.41;
	amp_ex2ex = 0.41;	

	// Create neurons and connection matrix
	neurons = vector<Neuron> (N, Neuron(_dt));
	w = new bool* [N];

	for (int m=0; m<N; m++)
	{
		if (m < Nl*Nl) // first Nl*Nl neurons are excitatory
			neurons[m].setType(TYPE_EXC);
		else // remaining neurons are inhibitory
			neurons[m].setType(TYPE_INH);
		w[m] = new bool [N];
	}
}

/*** Destructor ***
 * Cleans up the allocated memory for arrays */
~Network()
{
	//cout << "Destructor ~Network() was called!" << endl;
	for(int i=0; i<N; i++)
	//{
		//cout << "Line " << i+1;
		//cout << " (" << w[i][0];	
		//if (i < Nl)
		//{
			//cout << " " << neurons[i][0].getVoltage();
			//delete[] neurons[i];
			//neurons[i] = NULL;
			//cout << ") - Neurons deleted";
		//}	
		delete[] w[i];
		//w[i] = NULL;
		//cout << ") - connection deleted" << endl;
	//}
	//cout << "For-loop done" << endl;
	//delete[] neurons;
	delete[] w;
}

/* =============================================================================================================================== */
/* ==== Functions redirecting to corresponding functions in Neuron class ========================================================= */

/*** getVoltage ***
 * Returns the membrane potential of neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the membrane potential in mV */
double getVoltage(int i, int j) const
{
	return neurons[cNN(i,j)].getVoltage();
}
double getVoltage(int m) const
{
	return neurons[m].getVoltage();
}

/*** getCurrent ***
 * Returns total current effecting neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous current in nA */
double getCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getCurrent();
}
double getCurrent(int m) const
{
	return neurons[m].getCurrent();
}

/*** getChR2Current ***
 * Returns current evoked by the ChR2 molecules in neuron (i|j)*
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous ChR2 channel current in nA */
double getChR2Current(int i, int j) const
{
	return neurons[cNN(i,j)].getChR2Current();
}
double getChR2Current(int m) const
{
	return neurons[m].getChR2Current();
}

/*** getStimulusCurrent ***
 * Returns current evoked by external stimulation in neuron (i|j)*
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous current stimulus in nA */
double getStimulusCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getStimulusCurrent();
}
double getStimulusCurrent(int m) const
{
	return neurons[m].getStimulusCurrent();
}

/*** getFluctCurrent ***
 * Returns fluctuating current evoked by external synapses in neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous fluctuating current in nA */
double getFluctCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getFluctCurrent();
}
double getFluctCurrent(int m) const
{
	return neurons[m].getFluctCurrent();
}

/*** getConstCurrent ***
 * Returns the constant current elicited by the surrounding network (not this network!) in neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the constant current in nA */
double getConstCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getConstCurrent();
}
double getConstCurrent(int m) const
{
	return neurons[m].getConstCurrent();
}

/*** getSynapticCurrent ***
 * Returns the synaptic current that arrived in the previous time step *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the synaptic current in nA */
double getSynapticCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getSynapticCurrent();
}
double getSynapticCurrent(int m) const
{
	return neurons[m].getSynapticCurrent();
}

/*** getActivity ***
 * Returns true if neuron (i|j) is spiking in this instant of duration dt *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: whether neuron is firing or not */
bool getActivity(int i, int j) const
{
	return neurons[cNN(i,j)].getActivity();
}
bool getActivity(int m) const
{
	return neurons[m].getActivity();
}

/*** getSpikeTime ***
 * Returns the spike time for a given spike number (in temporal order, starting with 1) of neuron (i|j) *
 * - int n: the number of the considered spike (in temporal order, starting with 1) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the spike time for the n-th spike (or -1.0 if there exists none) */
double getSpikeTime(int n, int i, int j) const
{
	return neurons[cNN(i,j)].getSpikeTime(n);
}
double getSpikeTime(int n, int m) const
{
	return neurons[m].getSpikeTime(n);
}

/*** getSpikeCount ***
 * Returns the number of spikes that have occurred since the last reset of neuron (i|j) *
 * - return: the number of spikes */
int getSpikeCount(int i, int j) const
{
	return neurons[cNN(i,j)].getSpikeCount();
}
int getSpikeCount(int m) const
{
	return neurons[m].getSpikeCount();
}

/*** getSteadyVoltage ***
 * Returns the steady-state voltage (V_ss) for a constant light stimulus of given irradiance for neuron (i|j); *
 * returns only different values for different neurons when the constant current is distinct for the neurons *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the steady-state voltage at this irradiance in mV */
double getSteadyVoltage(double irradiance, int i, int j) const
{
	return neurons[cNN(i,j)].getSteadyVoltage(irradiance);
}
double getSteadyVoltage(double irradiance, int m) const
{
	return neurons[m].getSteadyVoltage(irradiance);
}

/*** getSteadyCurrent ***
 * Returns the steady-state current (I_ss) for a constant light stimulus of given irradiance for neuron (i|j); *
 * returns only different values for different neurons when the constant current is distinct for the neurons *
 * - double irradiance: the irradiance value in mW/mm^2 *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the steady-state current at this irradiance in nA */
double getSteadyCurrent(double irradiance, int i, int j) const
{
	return neurons[cNN(i,j)].getSteadyCurrent(irradiance);
}
double getSteadyCurrent(double irradiance, int m) const
{
	return neurons[m].getSteadyCurrent(irradiance);
}

/*** getSaturationCurrent ***
 * Returns the saturation current for a constant light stimulus in the infinite irradiance limit *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the saturation current for the set parameters in nA */
double getSaturationCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getSaturationCurrent();
}
double getSaturationCurrent(int m) const
{
	return neurons[m].getSaturationCurrent();
}

/*** getOpenProb ***
 * Returns the probability that a ChR2 channel of neuron (i|j) is open *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the open probability */
double getOpenProb(int i, int j) const
{
	return neurons[cNN(i,j)].getOpenProb();
}
double getOpenProb(int m) const
{
	return neurons[m].getOpenProb();
}

/*** setCurrentStimulus ***
 * Sets a current stimulus for neuron (i|j)
 * - Stimulus& _cst: shape of one stimulus period *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
void setCurrentStimulus(Stimulus& _cst, int i, int j)
{
	neurons[cNN(i,j)].setCurrentStimulus(_cst);
}
void setCurrentStimulus(Stimulus& _cst, int m)
{
	neurons[m].setCurrentStimulus(_cst);
}

/*** setLightStimulus ***
 * Sets a light stimulus for the ChR2 molecules in neuron (i|j) *
 * - Stimulus& _lst: shape of one stimulus period *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
void setLightStimulus(Stimulus& _lst, int i, int j)
{
	neurons[cNN(i,j)].setLightStimulus(_lst);
}
void setLightStimulus(Stimulus& _lst, int m)
{
	neurons[m].setLightStimulus(_lst);
}

/*** getCurrentStimulusAt ***
 * Returns the current stimulus magnitude at a given time *
 * - t_step: time step at which to evaluate stimulus
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getCurrentStimulusAt(int t_step, int i, int j) const
{
	return neurons[cNN(i,j)].getCurrentStimulusAt(t_step);
}
double getCurrentStimulusAt(int t_step, int m) const
{
	return neurons[m].getCurrentStimulusAt(t_step);
}

/*** getLightStimulusAt ***
 * Returns the light stimulus magnitude at a given time *
 * - t_step: time step at which to evaluate stimulus
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: stimulus at given time (if stimulus is not set, 0.0) */
double getLightStimulusAt(int t_step, int i, int j) const
{
	return neurons[cNN(i,j)].getLightStimulusAt(t_step);
}
double getLightStimulusAt(int t_step, int m) const
{
	return neurons[m].getLightStimulusAt(t_step);
}


/*** getNumberIncoming ***
 * Returns the number of either inhibitory or excitatory incoming connections to this neuron *
 * from other network neurons *
 * - int type: the type of incoming connections (inh./exc.)
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the number of incoming connections */
double getNumberIncoming(int type, int i, int j) const
{
	return neurons[cNN(i,j)].getNumberIncoming(type);
}
double getNumberIncoming(int type, int m) const
{
	return neurons[m].getNumberIncoming(type);
}

double getNumberOutgoing(int type, int i, int j) const
{
	return neurons[cInhG(i,j)].getNumberOutgoing(type);
}
double getNumberOutgoing(int type, int m) const
{
	return neurons[m].getNumberOutgoing(type);
}
/*** saveNeuronParams ***
 * Saves all the neuron parameters (including the channel parameters) to a given file; *
 * all neurons have the same parameters, so the first one is taken */
void saveNeuronParams(ofstream *f) const
{
	neurons[0].saveNeuronParams(f);
}

/* =============================================================================================================================== */

};


















