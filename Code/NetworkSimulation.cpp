/*************************************************************************
 *** Simulation of a neuronal network with Channelrhodopsin-2 channels ***
 ***           	(c) Jannik Luboeinski 2015/2016			       		  ***
 *************************************************************************/

#define NETWORKSIMULATION // important for Neuron class
#include <iostream>
#include <fstream>
#include <iomanip>
#include <windows.h>
//#include <stdlib.h>
using namespace std;

//#define STIMULUS_COMPARISON // performs a comparison between pulsed light stimulus and pulsed current stimulus
//#define CONNECTION_PLOT // creates a gnuplot file for a plot showing the detailed interneuronal connections within the excitatory network
//#define CONNECTION_PLOT_CREATION // actually creates the plot of interneuronal connections within the excitatory network
#define SPIKES_PLOT // creates a plot of the number of spikes per ms in the whole network as a function of time (only for minimum stimulus amplitude!)

#include "Tools.hpp"
#include "Plots.hpp"
#include "Network.hpp"

/*** NetworkSimulation class ***
 * simulates a network of neurons, has instances of Neuron class */
class NetworkSimulation {	

private:

/*** Simulation parameters ***/
const double dt; // ms, one time step for numerical integration
const double dE; // mW/mm^2, irradiance accuracy (in case of STIMULUS_COMPARISON: also dimension nA)
const int Nl; // number of neurons in one line (row or column) of the network
const int Nl_inh; // number of neurons in one line (row or column) of the network
const double center; // the center position of the excitatory neuron distribution
const double center_inh; // the center position of the inhibitory neuron distribution
const int t_offset; // time until the stimulus onset, just for plotting purposes (set it to zero for no effect) - does not contribute to t_max
double pc; // connection probability for unidirectional neuron connections
double tau_syn; // ms, synaptic time constant
double t_max;  // ms, total duration of simulation
double E_max; // mW/mm^2, irradiance to go up to (in case of STIMULUS_COMPARISON: also dimension nA)
double t_pulse; // ms, duration of one stimulus pulse
double frequency; // Hz, frequency for stimulus pulses
double J; // synaptic coupling strength (scale factor for all four different coupling strengths)
Network net;  // the network
#ifdef STIMULUS_COMPARISON
Network net2; // instance of another network 
#endif

/*** Output parameters ***/
const vector<double> stimplot; // stimulus amplitudes for which data and plots will be created (only multiples of dE with precision up to .2 are processed)
ofstream *contour; // pointer to a data file for creating contour plot (includes freq., intensity, mean maxima, mean FWHM, tau_adapt)
#ifdef SEEK_I_CONST
double *seekic; // pointer to a variable to communicate with NetworkBatch class while seeking I_const
#endif
string purpose; // a small text describing the purpose of this simulation

/*** considerAmplitude ***
 * Checks if given stimulus amplitude must be considered (based on array stimplot) *
 * - E: discretized irradiance (in units of dE)
 * - return: true, if given irradiance must be considered */
bool considerAmplitude(int E) const
{
	for (int i=0; i<stimplot.size(); i++)
	{
		if (int(round(stimplot[i]/dE)) == E)
			return true;
	}
	return false;	
}


/*** saveParams ***
 * Saves the crucial parameters in a file *
 * - str: string containing additional information, like the elapsed time */
void saveParams(const char* str) 
{
	ofstream f (dateStr("_PARAMS.txt"));
	f << dateStr("") << endl; // time stamp
	f << endl;

	// Parameters
	f << "Simulation parameters:" << endl;
	f << "dt = " << dt << " ms" << endl;
	f << "dE = " << dE << " mW/mm^2" << endl;
	f << "t_max = " << t_max << " ms (" << int(ceil(t_max / dt)) << " steps)" << endl;
	f << "t_offset = " << t_offset << " ms (" << int(ceil(t_offset / dt)) << " steps)" << endl;
	f << "E_max = " << E_max << " mW/mm^2 (" << int(ceil(E_max / dE)) + 1 << " steps)" << endl;
	f << "t_pulse = " << t_pulse << " ms" << endl;
	f << "frequency = " << frequency << " Hz" << endl;
	f << "J = " << J << endl;
	net.saveNetworkParams(&f);
	f << endl;

	// Additional information
	f << "Constant stimulation beyond " << 1000.0/t_pulse << " Hz" << endl;
	f << "Purpose: " << purpose << endl;
	f << str << endl;
	f.close();
}

public:

/*** Simulation function ***
 * Runs the network simulation *
 * - working_dir: working directory *
 * - first_sim: tells if this is the first simulation run by the batch code */
int simulate(string working_dir, bool first_sim)
{
// ==============================================================================================================================
	// Initialize

	// Start time measurement
	timeMeasure(true);

	// Demand the purpose
	if (first_sim)
	{
		cout << "Purpose: ";
		getline(cin, purpose);
	}

	// Constants
	const int m = int(ceil(E_max / dE)) + 1; // number of irradiance steps (including zero)
	const int n = int(ceil(t_max / dt)); // number of time steps
	const int offset = int(ceil(t_offset / dt)); // number of offset time steps
	const int pulse_length = int(ceil(t_pulse / dt)); // number of time steps of one pulse
	const int period_steps = int(round((1000/frequency)/dt)); // number of time steps of one period
	const string path = working_dir + "\\f=" + dtos(frequency, 1) + ",pc=" + dtos(pc, 3) + ",tau_syn=" + dtos(tau_syn, 0) 
									+ ",J=" + dtos(J, 1); // path to working directory
	const string separator = getSeparator(); cout << separator << endl; // string of characters for a separator in command line

	// Try to change directory
	system(concat("mkdir \"", path+string("\"")).c_str());
	
	if (!SetCurrentDirectory(path.c_str())) {
		showChDirErrMessage();
		cout << separator << endl;
		return -1;
	}

	// Output with general information
	// if this is the first simulation of this batch (first_sim = true), use time stamp from NetworkBatch.cpp, else, set a new time stamp
	cout << "\x1b[33mSimulating network with Ne = " << Nl*Nl << ", Ni = " << Nl_inh*Nl_inh 
		  << " for " << t_max << " ms (" << dateStr("", !first_sim) << ")\x1b[0m" << endl;
	cout << "Stimulus: rectangle, \x1b[34mt_pulse = " << t_pulse << " ms, \x1b[35mfrequency = " << frequency << " Hz\x1b[0m" << endl;
	cout << "Connectivity: \x1b[32mpc = " << pc << ", \x1b[31mtau_syn = " 
#ifdef DELTA_SYNAPSES
		  << 0
#else
		  << tau_syn 
#endif
		  << " ms, \x1b[36mJ = " << J << "\x1b[0m" << endl;
	cout << "Other parameters: \x1b[35mI_const = " << net.getConstCurrent(1,1) << " nA\x1b[0m" << endl;

	// Declarations and initializations
	double max_firing_rate = 0.0; // contains the maximum firing rate over all exc. neurons and stimulus amplitudes
	double min_firing_rate = numeric_limits<double>::max(); // contains the minimum firing rate over all exc. neurons and stimulus amplitudes
	double max_firing_rate_inh = 0.0; // contains the maximum firing rate over all inh. neurons and stimulus amplitudes
	double min_firing_rate_inh = numeric_limits<double>::max(); // contains the minimum firing rate over all inh. neurons and stimulus amplitudes
	ofstream gpl_mf("mean_fr.gpl"); // gnuplot file for creating a plot of mean firing rate over irradiance/current amplitude in a PDF file
	ofstream gpl_sf("gauss_sigma_fr.gpl"); // gnuplot file for creating a plot of the firing rate std. dev. over irradiance/current amplitude in a PDF file
	ofstream txt_msf(dateStr("_fr.txt")); // data file for creating plots over stimulus amplitude
	ofstream txt_stim; // data for creating map plot of the stimulus for specific stimulus amplitude
	ofstream* txt_fr = new ofstream[stimplot.size()]; // data for creating exc. plots for specific stimulus amplitudes over time
	ofstream* txt_fr_inh = new ofstream[stimplot.size()]; // data for creating inh. plots for specific stimulus amplitudes over time
	ofstream gpl_light_stim; // gnuplot file for creating irradiance map plot for last specific stimulus amplitude (only one is necessary, for all the other
									 // amplitudes the only difference is the colorcode range, which always goes up to the particular amplitude)
	ofstream gpl_curr_stim; // gnuplot file for creating current stimulus map plot for last specific stimulus amplitude (only one is necessary, for all the other
									 // amplitudes the only difference is the colorcode range, which always goes up to the particular amplitude)
	ofstream* gpl_lfr = new ofstream[stimplot.size()]; // gnuplot files for creating exc. firing rate map plots for specific light stimulus amplitudes
	ofstream* gpl_lfr_inh = new ofstream[stimplot.size()]; // gnuplot files for creating inh. firing rate map plots for specific light stimulus amplitudes
	ofstream* gpl_fit_light = new ofstream[stimplot.size()]; // gnuplot files for Gauss sigma fit with light stimulus for specific light stimulus amplitudes
#ifdef STIMULUS_COMPARISON
	ofstream* gpl_cfr = new ofstream[stimplot.size()]; // gnuplot files for creating exc. firing rate map plots for specific current stimulus amplitudes
	ofstream* gpl_cfr_inh = new ofstream[stimplot.size()]; // gnuplot files for creating inh. firing rate map plots for specific current stimulus amplitudes
	ofstream* gpl_fit_current = new ofstream[stimplot.size()]; // gnuplot files for Gauss sigma fit with current stimulus for specific light stimulus amplitudes
#endif
#ifdef CONNECTION_PLOT
	ofstream gpl_cp("connections_arrows.gpl"); // gnuplot file for creating a plot of the interneuronal connections
#endif
	ofstream gpl_cn("inc_connection_map.gpl"); // gnuplot file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream txt_cn(dateStr("_inc_connection.txt")); // data file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream txt_cf_in(dateStr("_connection_function_in.txt")); // data file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream txt_cf_ex(dateStr("_connection_function_ex.txt")); // data file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream gplscript("gpl"); // shell script for calling all the gnuplot scripts (to re-generate the plots, e.g. for another gnuplot version)
#ifdef SPIKES_PLOT
	ofstream gpl_spikes("spike_num.gpl"); // gnuplot file for creating a plot of the number of spikes over time
	ofstream txt_spikes(dateStr("_spike_num.txt")); // data file for creating a plot of the number of spikes over time
#endif
	ofstream logf("gslfit.log"); //TODO remove?

	writePalParula(); // create palette file for gnuplot color plots

	Stimulus lst = Stimulus(period_steps); // light stimulus
#ifdef STIMULUS_COMPARISON
	Stimulus cst = Stimulus(period_steps); // current stimulus
#endif

	double p = 0.0; // percentage of process completeness
	int current_E = 0; // specifies current position (index) in stimplot[] array

// ==============================================================================================================================
	// Prepare files for writing

	if (stimplot.size() > 0)
	{
		// Prepare file for irradiance map plot for specific stimulus amplitude
		string filename;
		filename = "_irradiance_" + dtos(stimplot.back(), 2) + ".txt"; // choose last specific stimulus amplitude (stimplot.back()) for this plot
		txt_stim.open(dateStr(filename));
		filename = "irradiance_" + dtos(stimplot.back(), 2) + "_map.gpl";
		gpl_light_stim.open(filename);
#ifdef STIMULUS_COMPARISON
		filename = "currentstim_" + dtos(stimplot.back(), 2) + "_map.gpl";
		gpl_curr_stim.open(filename);
#endif
		if (!txt_stim.is_open() || !gpl_light_stim.is_open()
#ifdef STIMULUS_COMPARISON
			 || !gpl_curr_stim.is_open()
#endif
												 )
		{
			cout << "Unable to open file!" << endl << separator << endl;
			return -1;
		}

		// Prepare files for specific stimulus plots
		for (int i=0; i<stimplot.size(); i++)
		{
			filename = "_fr_exc_" + dtos(stimplot[i], 2) + ".txt";
			txt_fr[i].open(dateStr(filename)); // open data file
			filename = "fr_exc_" + dtos(stimplot[i], 2) + "_map_light.gpl";
			gpl_lfr[i].open(filename); // open gnuplot file for firing rate map (exc. population) with light stimulus

			filename = "_fr_inh_" + dtos(stimplot[i], 2) + ".txt";
			txt_fr_inh[i].open(dateStr(filename)); // open data file
			filename = "fr_inh_" + dtos(stimplot[i], 2) + "_map_light.gpl";
			gpl_lfr_inh[i].open(filename); // open gnuplot file for firing rate map (inh. population) with light stimulus

			filename = "gauss_sigma_fit_" + dtos(stimplot[i], 2) + "_light.gpl";
			gpl_fit_light[i].open(filename); // open gnuplot file for Gauss sigma fit with light stimulus
#ifdef STIMULUS_COMPARISON
			filename = "fr_exc_" + dtos(stimplot[i], 2) + "_map_current.gpl";
			gpl_cfr[i].open(filename); // open gnuplot file for firing rate map (exc. population) with current stimulus

			filename = "fr_inh_" + dtos(stimplot[i], 2) + "_map_current.gpl";
			gpl_cfr_inh[i].open(filename); // open gnuplot file for firing rate map (inh. population) with current stimulus

			filename = "gauss_sigma_fit_" + dtos(stimplot[i], 2) + "_current.gpl";
			gpl_fit_current[i].open(filename); // open gnuplot file for Gauss sigma fit with current stimulus
#endif
			if ( !txt_fr[i].is_open() || !gpl_lfr[i].is_open() || !txt_fr_inh[i].is_open() || !gpl_lfr_inh[i].is_open() || !gpl_fit_light[i].is_open()
#ifdef STIMULUS_COMPARISON
					|| !gpl_cfr[i].is_open() || !gpl_cfr_inh[i].is_open() || !gpl_fit_current[i].is_open()
#endif
				)
			{
				cout << "Unable to open file!" << endl << separator << endl;
				return -1;
			}
		}
	}
	
	// Check if files have been opened properly
  	if ( !gplscript.is_open() || !gpl_mf.is_open() || !txt_msf.is_open() || ((contour != NULL) && (!contour->is_open()))	)
	{
		cout << "Unable to open file!" << endl << separator << endl;
		return -1;
	}

	if (
#ifdef CONNECTION_PLOT
		!gpl_cp.is_open() || 
#endif
#ifdef SPIKES_PLOT
		!gpl_spikes.is_open() || 
		!txt_spikes.is_open() ||
#endif
		!gpl_cn.is_open() || !txt_cn.is_open())
	{
		cout << "Unable to open file!" << endl << separator << endl;
		return -1;
	}


// ==============================================================================================================================
	// Output script for connection plot

	int total_c_count_exc = 0; // total number of exc. connections
	int total_c_count_inh = 0; // total number of inh. connections
	//int max_c_count = 0; // the maximum number of incoming connections received by a neuron // TODO REMOVE?
#ifdef CONNECTION_PLOT
	int color_n = 0;
	char color[10];
	gpl_cp << "set term pdf enhanced font \"FrontPage, 12\" color lw 1" << endl;
	gpl_cp << "set output '" << dateStr("_connections.pdf'") << endl << endl;
	gpl_cp << "set xrange [1" << ":" << Nl << "]" << endl;
	gpl_cp << "set yrange [1" << ":" << Nl << "]" << endl;
	gpl_cp << "set xlabel \"Neuron\"" << endl;
	gpl_cp << "set ylabel \"Neuron\"" << endl;
#endif
	//net.setConnection(3,5,1/sqrt(2),2,12,12,0,0,0.6,0.6);
		/*
    // Single neuron in the middle - EXCITATORY
    int a = (Nl*Nl/2+Nl/2); // FOR MIDDLE NEURON
    //int a = (Nl*Nl+Nl_inh/8+(Nl_inh*Nl_inh)/8); // SIDE NEURON
    for (int k = 1; k <= 2*Nl; k ++) 
	{
		for (int l = 1; l <= 2*Nl; l++)
    	{
    		int test = rowNN(a);
    		int test2 = colNN(a);
    		if (k==2*rowNN(a) and l == 2*colNN(a)) {
    			// The neuron to be considered for distribution
    			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 10 << endl;
				continue;			
			}
    		if ((k % 2 != 0) && (l % 2 != 0)) // excitatory..
    		{	
    		    int n = cNN((k+1)/2,(l + 1)/2);
				if (net.areConnected(a,n)) {
					txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << -4 << endl;
					continue;			
				}
				txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << -1 << endl;
				continue;
			}
			else if (((k % 2 == 0) && (l % 2 == 0)) && ((k/2) % 2 != 0) && ((l / 2) % 2 != 0)) // Inhibitory
			{
    		    int n = cInhG(((k+2)/4),((l+2)/4)); // Inhibitory
				if (net.areConnected(a,n)) {
					txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 6 << endl;			
					continue;
				}
				txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 1 << endl;
				continue;

			}
			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 0 << endl;			
		}
		txt_cn << endl; 
	}
	txt_cn.close();
	createNetworkColorPlot(gpl_cn, 2*Nl, -1.0, 3, "inc_connection", "", true, "Connection distribution"); //, 0.0, double(max_c_count), 1);
	gplscript << "gnuplot inc_connection_map.gpl" << endl; */
    // Single neuron in the middle - INHIBITORY
    int a = (Nl*Nl+Nl_inh/2+(Nl_inh*Nl_inh)/2); // FOR MIDDLE NEURON
    //int a = (Nl*Nl+Nl_inh/8+(Nl_inh*Nl_inh)/8); // SIDE NEURON
    for (int k = 1; k <= 2*Nl; k ++) 
	{
		for (int l = 1; l <= 2*Nl; l++)
    	{
    		int test = rowInh(a);
    		int test2 = colInh(a);
    		if (k==2*rowInh(a) and l == 2*colInh(a)) {
    			// The neuron to be considered for distribution
    			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 10 << endl;
				continue;			
			}
    		if ((k % 2 != 0) && (l % 2 != 0)) // excitatory..
    		{	
    		    int n = cNN((k+1)/2,(l + 1)/2);
				if (net.areConnected(a,n)) {
					txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << -4 << endl;
					continue;			
				}
				txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << -1 << endl;
				continue;
			}
			else if (((k % 2 == 0) && (l % 2 == 0)) && ((k/2) % 2 != 0) && ((l / 2) % 2 != 0)) // Inhibitory
			{
    		    int n = cInhG(((k+2)/4),((l+2)/4)); // Inhibitory
				if (net.areConnected(a,n)) {
					txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 6 << endl;			
					continue;
				}
				txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 1 << endl;
				continue;

			}
			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << 0 << endl;			
		}
		txt_cn << endl; 
	}
	txt_cn.close();
	createNetworkColorPlot(gpl_cn, 2*Nl, -1.0, 3, "inc_connection", "", true, "Connection distribution"); //, 0.0, double(max_c_count), 1);
	gplscript << "gnuplot inc_connection_map.gpl" << endl; 
	//cout << "Number of I->E connections: " << net.getNumberOutgoing(TYPE_INH,a) << endl;

	//Network group [10] {{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30},{0.1,60,30}};
	
	// Gaussian distribution test
	int counter = 0;
	int counter_total = 0;
	int vec_x = 0;
	int vec_y = 0;
	int fixed = a;// (Nl*Nl+Nl_inh/2+(Nl_inh*Nl_inh)/2); 
	double gaussian_result = 0.0;
	
	for (int entry = 1; entry < 2; entry ++)
	{
		Network net_new = Network(0.1,60,30);
		//net_new.setConnection(3,5,1/sqrt(2),2,15,20,1,1/sqrt(2),0.4,0.4);
		//net_new.setConnection(3,5,1/sqrt(2),2,12,12,0,0,0.6,0.6);
		net_new.setPC(0.001);
		int counter_exex = 0; // / 3600
		int counter_exin = 0; // / 3600 
		int counter_inex = 0; // / 900
		int counter_inin = 0; // / 900
		//fixed = 2065; 
		for (int neuronnumber = 1; neuronnumber < 3600; neuronnumber ++)
		{
			counter_exex = counter_exex + net_new.getNumberIncoming(TYPE_EXC,neuronnumber);
			counter_exin = counter_exin + net_new.getNumberIncoming(TYPE_INH,neuronnumber);
		}
		cout << "Ex2Ex connections" << counter_exex << endl;
		cout << "Ex2Ex connections mean: " << counter_exex/3600 << endl;
		cout << "In2Ex connections" << counter_exin << endl;
		cout << "In2Ex connections mean: " << counter_exin/3600 << endl;
		for (int neuronnumber = 3600; neuronnumber < 4500; neuronnumber ++)
		{
			counter_inex = counter_inex + net_new.getNumberIncoming(TYPE_EXC,neuronnumber);
			counter_inin = counter_inin + net_new.getNumberIncoming(TYPE_INH,neuronnumber);
		}
		cout << "Ex2In connections" << counter_inex << endl;
		cout << "Ex2In connections mean: " << counter_inex/900 << endl;
		cout << "In2In connections" << counter_inin << endl;
		cout << "In2In connections mean: " << counter_inin/900 << endl;
			}
	// FOR EXCITATORY NEURONS	
	/*
	for (double x = 1.0; x <= Nl; x++){ // NL instead of Nl_inh, for inhibitory neurons being in a 60x60 matrix by only being 30x30 in number
		for (double y = 0.0; y <=x; y++){
			for (int z = 0; z < Nl * Nl; z ++ ){
				if (net_new.areConnected(fixed,z) && (sqrt(pow(rowNN(fixed)-rowNN(z),2)+pow(colNN(fixed)-colNN(z),2)) == sqrt(pow(x,2)+pow(y,2)))){
					counter += 1;
					counter_total += 1;
				}	
				else if (!(net_new.areConnected(fixed,z)) && (sqrt(pow(rowNN(fixed)-rowNN(z),2)+pow(colNN(fixed)-colNN(z),2)) == sqrt(pow(x,2)+pow(y,2))))
				{
					counter_total += 1;
				}		
			}
			//gaussian_result = (4.0/10.0)*exp((-1*(pow(sqrt(pow(x,2)+pow(y,2))-1.0,2)))/(2*pow(15.0,2)));
			//txt_cf_ex << std::fixed << sqrt(pow(x,2)+pow(y,2)) << "\t\t\t" << counter << "\t\t\t" << counter_total << "\t\t\t" ;
			//txt_cf_ex << std::fixed << setprecision(10) << gaussian_result << endl;
			//counter = 0;
			//counter_total = 0;
		}
	}}	
	for (int neuronnumber = 1; neuronnumber < 900
		for (double x = 1.0; x <= Nl; x++){
			for (double y = 0.0; y <=x; y++){
				for (int z = 1; z < Nl_inh * Nl_inh; z ++ ){
					double test_A = sqrt(pow(rowInh(z+Nl*Nl)-rowNN(fixed)+0.5,2)+pow(colInh(z+Nl*Nl)-colNN(fixed)+0.5,2));
					double test_B =  (sqrt(pow(x/2.0,2)+pow(y/2.0,2)));
					if (net_new.areConnected(fixed,z+Nl*Nl) && (sqrt(pow(rowInh(z+Nl*Nl)-rowNN(fixed)+0.5,2)+pow(colInh(z+Nl*Nl)-colNN(fixed)+0.5,2)) == (sqrt(pow(x/2.0,2)+pow(y/2.0,2))))){
						// + 1 because x-y is a 1:1 system while inhibitory normally are in between (0.5), so the distances are scaled twice
						counter += 1;
						counter_total += 1;
					}	
					else if (!(net_new.areConnected(fixed,z+Nl*Nl)) && (sqrt(pow(rowInh(z+Nl*Nl)-rowNN(fixed)+0.5,2)+pow(colInh(z+Nl*Nl)-colNN(fixed)+0.5,2)) == (sqrt(pow(x/2.0,2)+pow(y/2.0,2)))))
					{
						counter_total += 1;
					}		
				}
				gaussian_result = (4.0/10.0)*exp((-1*(pow(sqrt(pow(x/2.0,2)+pow(y/2.0,2))-(1/sqrt(2.0)),2)))/(2*pow(20,2)));
				txt_cf_in << std::fixed << sqrt(pow(x/2.0,2)+pow(y/2.0,2)) << "\t\t\t" << counter << "\t\t\t" << counter_total << "\t\t\t";
				txt_cf_in << std::fixed << setprecision(10) << gaussian_result << endl;
				counter = 0;
				counter_total = 0;
			}
		} 	  
	}	
		
		/* FOR INHIBITORY NEURONS
	for (double x = 1.0; x <= Nl/2; x++){ // NL instead of Nl_inh, for inhibitory neurons being in a 60x60 matrix by only being 30x30 in number
		for (double y = 0.0; y <=x; y++){
			for (int z = 0; z < Nl_inh * Nl_inh; z ++ ){
				int a = (sqrt(pow(rowInh(fixed)-rowInh(z+Nl*Nl),2)+pow(colInh(fixed)-colInh(Nl*Nl+z),2)))/2;
				int b = sqrt(pow(x,2)+pow(y,2));
				if (net_new.areConnected(fixed,z+Nl*Nl) && (sqrt(pow(rowInh(fixed)-rowInh(z+Nl*Nl),2)+pow(colInh(fixed)-colInh(z+Nl*Nl),2)) == sqrt(pow(x,2)+pow(y,2)))){
					counter += 1;
					counter_total += 1;
				}	
				else if (!(net_new.areConnected(fixed,z+Nl*Nl)) && (sqrt(pow(rowInh(fixed)-rowInh(z+Nl*Nl),2)+pow(colInh(fixed)-colInh(z+Nl*Nl),2)) == sqrt(pow(x,2)+pow(y,2))))
				{
					counter_total += 1;
				}		
			}
			gaussian_result = (10.0/10.0)*exp((-1*(pow(sqrt(pow(x,2)+pow(y,2))-2.0,2)))/(2*pow(5.0,2)));
			txt_cf_in << std::fixed << sqrt(pow(x,2)+pow(y,2)) << "\t\t\t" << counter << "\t\t\t" << counter_total << "\t\t\t" ;
			txt_cf_in << std::fixed << setprecision(10) << gaussian_result << endl;
			counter = 0;
			counter_total = 0;
		}
	}	
		for (double x = 1.0; x <= Nl/2; x++){
			for (double y = 0.0; y <=x; y++){
				for (int z = 0; z < Nl * Nl; z ++ ){
					if (net_new.areConnected(fixed,z) && (sqrt(pow(rowInh(fixed)-rowNN(z)+0.5,2)+pow(colInh(fixed)-colNN(z)+0.5,2)) == (sqrt(pow(x/2.0,2)+pow(y/2.0,2))))){
						// + 1 because x-y is a 1:1 system while inhibitory normally are in between (0.5), so the distances are scaled twice
						counter += 1;
						counter_total += 1;
					}	
					else if (!(net_new.areConnected(fixed,z)) && (sqrt(pow(rowInh(fixed)-rowNN(z)+0.5,2)+pow(colInh(fixed)-colNN(z)+0.5,2)) == (sqrt(pow(x/2.0,2)+pow(y/2.0,2)))))
					{
						if ((sqrt(pow(x/2.0,2)+pow(y/2.0,2))) == sqrt(0.5)) {
							cout << "ERROR AT " << z << endl;
						}
						counter_total += 1;
					}		
				}
				gaussian_result = (10.0/10.0)*exp((-1*(pow(sqrt(pow(x/2.0,2)+pow(y/2.0,2))-(1/sqrt(2.0)),2)))/(2*pow(3.0,2)));
				txt_cf_ex << std::fixed << sqrt(pow(x/2.0,2)+pow(y/2.0,2)) << "\t\t\t" << counter << "\t\t\t" << counter_total << "\t\t\t";
				txt_cf_ex << std::fixed << setprecision(10) << gaussian_result << endl;
				counter = 0;
				counter_total = 0;
			}
		} 	  
	}	*/
	//text_cf.close();
	/*
	ofstream Dingens("machwas.gpl");
	if (Dingens.is_open()){
		Dingens << "set term pdf enhanced font \"FrontPage, 20\" color lw 1" << endl;
		Dingens << "set output '" << dateStr("_connection_function.pdf") << "'" << endl;
//		Dingens << "plot '" << dateStr("_connection_function_in2ex.txt") << "' using 1:($2 / $3) notitle, 1 * exp( -(x -(1/sqrt(2)))**2 / (2*3**2) )" << endl; // GAUSS TODO 25.04 18:33
//		Dingens << "plot '" << dateStr("_connection_function_in2in.txt") << "' using 1:($2 / $3) notitle, 1 * exp( -((x-2)**2) / (2*5**2) )" << endl; // GAUSS TODO 25.04 18:33
		Dingens << "plot '" << dateStr("_connection_function_ex.txt") << "' using 1:($2 / $3) notitle, 0.4 * exp( -(x -1)**2 / (2*15**2) )" << endl; // GAUSS TODO 25.04 18:33
		Dingens << "plot '" << dateStr("_connection_function_in.txt") << "' using 1:($2 / $3) notitle, 0.4 * exp( -((x-1/sqrt(2))**2) / (2*20**2) )" << endl; // GAUSS TODO 25.04 18:33
	}
	else {
		return 0;
	} 
	Dingens.close(); */
	//system("gnuplot machwas.gpl"); 
	 
	 
	/*
	// loops over k and l replace one loop over m in consecutive numbering
	for (int k=1; k<=Nl_inh; k++) // row/column neuron numbering
	{
		for (int l=1; l<=Nl_inh; l++) // row/column neuron numbering
		{
			int m = cInhG(k,l); // achieve consecutive number
			//const int c_count_exc = net.getNumberIncoming(TYPE_EXC,k,l); // number of incoming excitatory connections to neuron m
			//const int c_count_inh = net.getNumberIncoming(TYPE_INH,k,l); // number of incoming inhibitory connections to neuron m
			const int c_count_inh = net.getNumberOutgoing(TYPE_INH,m);
#ifdef CONNECTION_PLOT
			// draw connections between excitatory neurons
			for (int n=0; n<Nl*Nl; n++) // consecutive neuron numbering
			{
				switch (color_n)
				{
					case 0:
						strcpy(color, "red");
						break;
					case 1:
						strcpy(color, "blue");
						break;
					case 2:
						strcpy(color, "green");
						break;
					case 3:
						strcpy(color, "purple");
						break;
					default:
						strcpy(color, "orange");
						color_n = -1;
						break;
				}
			
				if (net.areConnected(n, m)) // check if there is connection n->m
				{
					gpl_cp << "set arrow from " << net.rowNN(n) << "," << net.colNN(n) << " to " << k << "," << l
							 << " head size 0.2,20 filled back lc rgb '" << color << "'" << endl;

					if (net.areConnected(m, n) && (m > n)) // check if there is also a connection m->n (if yes, draw a yellow dashed second arrow),
																		// but not to do this twice, only consider lower triangle of connection matrix (m > n condition)
					{
						gpl_cp << "set arrow from " << k << "," << l << " to " << net.rowNN(n) << "," << net.colNN(n) 
								 << " head size 0.2,20 filled back lc rgb 'yellow' dt 3" << endl;
					}
					color_n++;
				}
			}

			// draw neurons
			if (c_count_exc > 0) // if neurons possesses incoming connections, draw filled circle
				gpl_cp << "set object circle at first " << k << "," << l
					 	 << " radius char 0.15 front fillcolor rgb 'black' fillstyle solid noborder" << endl;
			else  // if not, draw empty circle
				gpl_cp << "set object circle at first " << k << "," << l
						  << " radius char 0.15 front fillstyle empty border lc rgb 'black' lw 1" << endl;
#endif //CONNECTION_PLOT

			//if (max_c_count < c_count_exc)
			//	max_c_count = c_count_exc;
			//total_c_count_exc += c_count_exc;
			total_c_count_inh += c_count_inh;

			txt_cn << fixed << k << "\t\t\t" << l << "\t\t\t" << c_count_inh << endl;
		} // end of for(l)
		txt_cn << endl; // empty line for color plot
	} 
 */
#ifdef CONNECTION_PLOT
	gpl_cp << "plot 0 notitle lc rgb 'white'" << endl;
	gpl_cp.close();
#ifdef CONNECTION_PLOT_CREATION
	system("gnuplot connections_arrows.gpl"); // resulting PDF file needs several MB of disk space
#endif
#endif //CONNECTION_PLOT

	txt_cn.close();
	createNetworkColorPlot(gpl_cn, Nl, -1.0, 3, "inc_connection", "", true, "# of incoming connections"); //, 0.0, double(max_c_count), 1);
	gplscript << "gnuplot inc_connection_map.gpl" << endl;
	//
	//cout << "Number of E->E connections: " << total_c_count_exc << " (expected: " << pc * (pow(Nl*Nl, 2)-(Nl*Nl)) << ")" << endl;
//	cout << "Number of I->E connections: " << total_c_count_inh << " (expected: " << pc * pow(Nl, 2) * pow(Nl_inh, 2) << ")" << endl;

// ==============================================================================================================================
	// (Light/current) stimulation loop

	for(int i = 0; i < m; i++)
	{

		// Reset network
		net.reset();
#ifdef STIMULUS_COMPARISON
		net2.reset();
#endif

		// Clear stimuli
		lst.clear();
#ifdef STIMULUS_COMPARISON
		cst.clear();
#endif
#ifdef SPIKES_PLOT
		int spike_num = 0; // reset number of spikes in currently considered interval
#endif

		// Time loop
		for (int j = 0; j <= n+offset; j++)
		{
			// Update percentage
			double p_new = double(round(double(i*(n+offset) + j) / double(m*(n+offset)) * 1000.0)) / 10.0; // round to first decimal place
			if (p_new > p) 
			{
				p = p_new;
				printf("\rProgress: %.1f %% completed.", p);
				fflush(stdout);
			}

			// Stimuli start at offset time
			if (j == offset)
			{
				// Set light stimulus
				lst.addRectPulse(i*dE, 0, pulse_length);
				net.setGaussianLightStimulus(lst, center, center); // set amplitude of light stimulus in the center of the network
				
				// Set current stimulus
#ifdef STIMULUS_COMPARISON
				cst.addRectPulse(net2.getSteadyCurrent(i*dE, 0), 0, pulse_length); // TODO steady current also if f < constant?
				net2.setGaussianCurrentStimulus(cst, center, center); // set amplitude of light stimulus in the center of the network
#endif
			}
			
			// Calculate next step for Network (and implicitly for Neurons and Channelrhodopsin)
#ifdef SPIKES_PLOT
			spike_num +=
#endif
				net.processTimeStep(j-offset); // necessarily negative time steps during offset - those spikes will not be counted
#ifdef SPIKES_PLOT
			if ((i == 0) && (j-offset) % int(1./dt) == 0)
			{
				txt_spikes << (j-offset)*dt << "\t\t\t" << spike_num << endl; // time and number of spikes within one millisecond
				spike_num = 0; // reset spike count
			}
#endif
#ifdef STIMULUS_COMPARISON
			net2.processTimeStep(j-offset);
#endif
			
		} // end of for(j)

// ==============================================================================================================================
		// Compute mean firing rate, standard deviation for the firing rates and fit the Gaussian distribution

		double mfr_light = 0.; // mean firing rate of the whole excitatory population in the light-stimulated network
		double mfr_light_inh = 0.; // mean firing rate of the whole inhibitory population in the light-stimulated network
		double sdfr_light = 0.; // standard deviation of the firing rate of the whole excitatory population in the light-stimulated network
		double sdfr_light_inh = 0.; // standard deviation of the firing rate of the whole inhibitory population in the light-stimulated network
		int sdfr_light_part1 = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part2 = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part1_inh = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_light_part2_inh = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int denom_light = 4; // denominator for sum of slice results
		/* Sigma Gauss determination by average fitting
		int nu_max = net.getSpikeCount(center, center); // number of spikes of center neuron in light-stimulated network 
		int nu_rmax = 0; 
		int denom_light = Nl*Nl - 1; // denominator for sigma_gauss_light*/

		double sigma_gauss_light = 0.; // standard deviation of Gaussian firing rate distribution evoked by light stimulus (not to be confused with standard deviation 												 // of firing rates!)
		double sigma_gauss_err_light = 0.; // error of the standard deviation of Gaussian firing rate distribution evoked by light stimulus
		double sigma_gauss_light_inh = 0.; // standard deviation of Gaussian firing rate distribution in inhibitory population evoked by light stimulus (not to be
															// confused with standard deviation of firing rates!)
		double sigma_gauss_err_light_inh = 0.; // error of the standard deviation of Gaussian firing rate distribution in inhibitory population evoked by light stimulus
#ifdef STIMULUS_COMPARISON
		double mfr_current = 0.; // mean firing rate of the whole excitatory population in the current-stimulated network
		double mfr_current_inh = 0.; // mmean firing rate of the whole inhibitory population in the current-stimulated network
		double sdfr_current = 0.; // standard deviation of the firing rate of the whole excitatory population in the light-stimulated network
		double sdfr_current_inh = 0.; // standard deviation of the firing rate of the whole inhibitory population in the light-stimulated network
		int sdfr_current_part1 = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part2 = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part1_inh = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int sdfr_current_part2_inh = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
		int denom_current = 4; // denominator for sum of slice results
		/* Sigma Gauss determination by average fitting
		int nu_max2 = net2.getSpikeCount(center, center); // number of spikes of center neuron in current-stimulated network
		int denom_current = Nl*Nl - 1; // denominator for sigma_gauss_light */
		
		double sigma_gauss_current = 0.; // standard deviation of Gaussian firing rate distribution evoked by current stimulus (not to be confused with standard 													// deviation of firing rates!)
		double sigma_gauss_err_current = 0.; // error of the standard deviation of Gaussian firing rate distribution evoked by current stimulus
		double sigma_gauss_current_inh = 0.; // standard deviation of Gaussian firing rate distribution in inhibitory population evoked by current stimulus (not to be
															// confused with standard deviation of firing rates!)
		double sigma_gauss_err_current_inh = 0.; // error of the standard deviation of Gaussian firing rate distribution in inhibitory population evoked by current stimulus
		
#endif

		
		//double nu_slice[4][Nl]; // 2-dimensional array containing the firing rates within all four "slices" of firing rates that are used to fit
		//double sigma_gauss[4]; fill_n(sigma_gauss, 4, 0.); // array to contain the final sigma values of all 4 slices
		//double sigma_gauss_err[4]; fill_n(sigma_gauss_err, 4, 0.); // array to contain the final sigma error values of all 4 slices
		//double nu_err[Nl]; fill_n(nu_err, Nl, 0.2); // data error estimates for firing rates -- estimate by 0.2 

		double fit_radius[Nl*Nl]; // x-values for fitting (radius)
		double fit_radius_inh[Nl_inh*Nl_inh]; // x-values for fitting (radius) in inh. population
		double fit_fr[Nl*Nl]; // y-values for fitting (firing rate of light-stimulated network)
		double fit_fr_inh[Nl_inh*Nl_inh]; // y-values for fitting (firing rate of light-stimulated network) in inh. population
#ifdef STIMULUS_COMPARISON
		double fit_fr2[Nl*Nl]; // y-values for fitting (firing rate of current-stimulated network)
#endif
		double fit_fr_err[Nl*Nl]; fill_n(fit_fr_err, Nl*Nl, 0.2); // data error estimates for firing rates -- estimate by 0.2 // TODO OK?
		double fit_fr_err_inh[Nl*Nl]; fill_n(fit_fr_err_inh, Nl_inh*Nl_inh, 0.2); // data error estimates for firing rates -- estimate by 0.2 // TODO OK?

		double b = 0.; // baseline shift of the Gaussian - should be zero, just as in the Gaussian stimulus
		double var_init = pow(net.getGaussSigma(),2); // use variance of Gaussian stimulus distribution as initial value
		double par_light[3] = {0., var_init, b}; // set initial values for A & var and fixed value for b (for light stimulus)
		double par_light_inh[3] = {0., var_init, b}; // set initial values for A & var and fixed value for b (for light stimulus)
		double err_light[3]; // parameter error array (for light stimulus)
#ifdef STIMULUS_COMPARISON
		double par_curr[3] = {0., var_init, b}; // set initial values for A & var and fixed value for b (for current stimulus)
		double err_curr[3]; // parameter error array (for current stimulus)
#endif
		
		// OUTPUT (TODO REMOVE)
		//ofstream diagfit(dateStr("diagfit_") + dtos(i*dE,2) + ".txt");
		//if (!diagfit.is_open()) {
		//	cout << "Error while opening diagfit file!" << endl;
		//	return -3;
		//} //------

		for (int m=0; m<pow(Nl,2)+pow(Nl_inh,2); m++)
		{
			int nu = net.getSpikeCount(m); // get number of spikes of neuron m in light-stimulated network
#ifdef STIMULUS_COMPARISON
			int nu2 = net2.getSpikeCount(m); // get number of spikes of neuron m in current-stimulated network
#endif

			if (m < pow(Nl,2)) // neuron m is in excitatory population
			{
				mfr_light += double(nu);
				sdfr_light_part1 += pow(nu,2);
				sdfr_light_part2 += nu;

				fit_radius[m] = sqrt(pow(rowNN(m)-center,2) + pow(colNN(m)-center,2)); // store radius in fit data array
				fit_fr[m] = nu/(t_max/1000.);	// store firing rate in fit data array

#ifdef STIMULUS_COMPARISON
				mfr_current += double(nu2);
				sdfr_current_part1 += pow(nu2,2);
				sdfr_current_part2 += nu2;

				fit_fr2[m] = nu2/(t_max/1000.);	// store firing rate in fit data array
#endif
			}
			else // neuron m is in inhibitory population
			{
				int mprime = m - pow(Nl,2);

				mfr_light_inh += double(nu);
				sdfr_light_part1_inh += pow(nu,2);
				sdfr_light_part2_inh += nu;

				fit_radius_inh[mprime] = sqrt(pow(rowInh(mprime)-center_inh,2) + pow(colInh(mprime)-center_inh,2)); // store radius in fit data array
				fit_fr_inh[mprime] = nu/(t_max/1000.);	// store firing rate in fit data array

#ifdef STIMULUS_COMPARISON
				mfr_current_inh += double(nu2);
				sdfr_current_part1_inh += pow(nu2,2);
				sdfr_current_part2_inh += nu2;
#endif
			}
		}

		// OUTPUT (TODO REMOVE)
		//for (int a=0; a<Nl; a++)
		//	diagfit << a+1 << "\t\t\t" << nu_slice[2][a] << "\t\t\t" << nu_slice[3][a] << endl;
		//diagfit.close();
		//--------

		if (i > 0) // computation of Gaussian sigma only makes sense if there is a stimulus
		{
			// findFit for excitatory population
			/*gsl_multifit_function_fdf func; // fit function plus its Jacobian (fdf)
			fitdata_pairs fd = {Nl*Nl, fit_radius, fit_fr, fit_fr_err}; // structure of arrays for the data points (x, y, yerror)
	
			// set the properties of the fit function
			func.f = &gauss_f; // the fit function
			func.df = &gauss_df_A_var; // the Jacobian of the fit function
			func.fdf = &gauss_fdf_A_var; // function calling _f and _df functions
			func.n = Nl*Nl; // number of data points
			func.p = 3; // number of parameters (both fit parameters and constant parameters)
			func.params = &fd; // fitdata structure (contains arrays for (x, y, yerror))

			if (findFit(func, par_light, err_light, NULL)) // do the fit procedure
			{
				sigma_gauss_light = sqrt(abs(par_light[1])); // read out sigma_gauss
				sigma_gauss_err_light = sqrt(abs(err_light[1])); // read out the sigma_gauss error
					
				logf << "Irradiance " << double(i)*dE << ": sigma = " << sigma_gauss_light << " +- " << sigma_gauss_err_light
					  << " (b = " << par_light[2] << ", A = " << par_light[0] << ")" << endl;
			}
			else 
			{
				logf << "Irradiance " << double(i)*dE << ": error" << endl;
			}

			if (Nl_inh > 0)
			{
				// findFit for inhibitory population
				fitdata_pairs fd_inh = {Nl_inh*Nl_inh, fit_radius_inh, fit_fr_inh, fit_fr_err_inh}; // structure of arrays for the data points (x, y, yerror)
	
				// set the properties of the fit function
				func.f = &gauss_f; // the fit function
				func.df = &gauss_df_A_var; // the Jacobian of the fit function
				func.fdf = &gauss_fdf_A_var; // function calling _f and _df functions
				func.n = Nl_inh*Nl_inh; // number of data points
				func.p = 3; // number of parameters (both fit parameters and constant parameters)
				func.params = &fd_inh; // fitdata structure (contains arrays for (x, y, yerror))

				if (findFit(func, par_light_inh, err_light, NULL)) // do the fit procedure
				{
					sigma_gauss_light_inh = sqrt(abs(par_light_inh[1])); // read out sigma_gauss
					sigma_gauss_err_light_inh = sqrt(abs(err_light[1])); // read out the sigma_gauss error
					
					logf << "Irradiance " << double(i)*dE << ": sigma = " << sigma_gauss_light_inh << " +- " << sigma_gauss_err_light_inh
						  << " (b = " << par_light_inh[2] << ", A = " << par_light_inh[0] << ") (inh.)" << endl;
				}
				else 
				{
					logf << "Irradiance " << double(i)*dE << ": error (inh.)" << endl;
				}
			} */
	
			// TODO : the same for current stimulus
		}

/* Sigma Gauss determination by average fitting
		logf << nu_max << " / " << nu_rmax << " - " << denom_light << endl;*/
		mfr_light /= (double(Nl*Nl)*t_max / 1000.0); // compute the mean firing rate from the number of spikes
		if (Nl_inh > 0)
			mfr_light_inh /= (double(Nl_inh*Nl_inh)*t_max / 1000.0); // compute the mean firing rate from the number of spikes
		sdfr_light = sqrt( ( double(sdfr_light_part1) - double(pow(sdfr_light_part2, 2)) / pow(Nl, 2) ) / double(Nl*Nl-1) ) / (t_max / 1000.0); // compute standard 																														// deviation of the firing rates according to Steiner's translation theorem
		if (Nl_inh > 0)
			sdfr_light_inh = sqrt( ( double(sdfr_light_part1_inh) - double(pow(sdfr_light_part2_inh, 2)) / pow(Nl_inh, 2) ) / 
										double(Nl_inh*Nl_inh-1) ) / (t_max / 1000.0); // compute standard deviation of the firing rates according to Steiner's translation 																										 // theorem

		//sigma_gauss_light = //(sigma_gauss[0] + sigma_gauss[1] + sigma_gauss[2] + sigma_gauss[3]) / denom_light; // compute mean standard deviation
		//sigma_gauss_err_light = sqrt(pow(sigma_gauss_err[0],2) + pow(sigma_gauss_err[1],2) 
									  		  //+ pow(sigma_gauss_err[2],2) + pow(sigma_gauss_err[3],2)) / denom_light;  // compute mean error of the standard deviation

/* Sigma Gauss determination by average fitting
		sigma_gauss_light /= double(denom_light); // compute the Gaussian sigma - divide by number of considered neurons
#ifdef STIMULUS_COMPARISON
		sigma_gauss_current /= double(denom_current); 
#endif*/

#ifdef STIMULUS_COMPARISON
		mfr_current /= (double(Nl*Nl)*t_max / 1000.0);
		sdfr_current = sqrt( ( double(sdfr_current_part1) - double(pow(sdfr_current_part2, 2)) / pow(Nl, 2) ) / double(Nl*Nl-1) ) / (t_max / 1000.0);
		if (Nl_inh > 0)
			sdfr_current_inh = sqrt( ( double(sdfr_current_part1_inh) - double(pow(sdfr_current_part2_inh, 2)) / pow(Nl_inh, 2) ) 
										/ double(Nl_inh*Nl_inh-1) ) / (t_max / 1000.0);
		//sigma_gauss_current = (sigma_gauss[0] + sigma_gauss[1] + sigma_gauss[2] + sigma_gauss[3]) / denom_current; // TODO ??
		//sigma_gauss_current = (sigma_gauss[0] + sigma_gauss[1] + sigma_gauss[2] + sigma_gauss[3]) / denom_current;
#endif

// ==============================================================================================================================
		// Output in data files

		// Write mean firing rate, standard deviation of the firing rate and Gauss sigma of the firing rate (over irradiance/current) in data file
		txt_msf << fixed << double(i)*dE << "\t\t\t" //1 Light intensity
				  << dtos(mfr_light, 6, true) << "\t\t\t" //2 Mean firing rate
				  << dtos(sdfr_light, 6, true) << "\t\t\t" //3 Standard deviation of the mean firing rate
				  << dtos(sigma_gauss_light, 6, true) << "\t\t\t" //4 Gaussian sigma of the firing rate
				  << dtos(sigma_gauss_err_light, 6, true) << "\t\t\t" //5 Fit error for Gaussian sigma of the firing rate
				  << dtos(mfr_light_inh, 6, true) << "\t\t\t" //6 Mean firing rate (inh. population)
				  << dtos(sdfr_light_inh, 6, true) << "\t\t\t" //7 Standard deviation of the mean firing rate (inh. population)
				  << dtos(sigma_gauss_light_inh, 6, true) << "\t\t\t" //8 Gaussian sigma of the firing rate (inh. population)
				  << dtos(sigma_gauss_err_light_inh, 6, true); //9 Fit error for Gaussian sigma of the firing rate (inh. population)
#ifdef STIMULUS_COMPARISON
		txt_msf << "\t\t\t" << net2.getSteadyCurrent(i*dE, 0) << "\t\t\t" // TODO//10 Current amplitude
				  << dtos(mfr_current, 6, true) << "\t\t\t" //11 Mean firing rate
				  << dtos(sdfr_current, 6, true) << "\t\t\t" //12 Standard deviation of the mean firing rate
				  << dtos(sigma_gauss_current, 6, true) << "\t\t\t" //13 Gaussian sigma of the firing rate
				  << dtos(sigma_gauss_err_current, 6, true) << "\t\t\t" //14 Fit error for Gaussian sigma of the firing rate
				  << dtos(mfr_current_inh, 6, true) << "\t\t\t" //15 Mean firing rate (inh. population)
				  << dtos(sdfr_current_inh, 6, true) << "\t\t\t" //16 Standard deviation of the mean firing rate (inh. population)
				  << dtos(sigma_gauss_current_inh, 6, true) << "\t\t\t" //17 Gaussian sigma of the firing rate (inh. population)
				  << dtos(sigma_gauss_err_current_inh, 6, true); //18 Fit error for Gaussian sigma of the firing rate (inh. population)
#endif
		txt_msf << endl;


		// Write contour plot file
		if (contour != NULL)
		{
			*contour << fixed << frequency << "\t\t\t" << pc << "\t\t\t" << tau_syn << "\t\t\t" << J << "\t\t\t" << double(i)*dE << "\t\t\t"
						<< dtos(mfr_light, 6, true) << "\t\t\t"
						<< dtos(sigma_gauss_light, 6, true) << "\t\t\t"
						<< dtos(mfr_light_inh, 6, true) << "\t\t\t"
						<< dtos(sigma_gauss_light_inh, 6, true);
#ifdef STIMULUS_COMPARISON // TODO: current stimulus?
			*contour << "\t\t\t" << dtos(mfr_current, 6, true)
					   << "\t\t\t" << dtos(sigma_gauss_current, 6, true)
						<< "\t\t\t" << dtos(mfr_current_inh, 6, true)
					   << "\t\t\t" << dtos(sigma_gauss_current_inh, 6, true);
#endif
			*contour << endl;
		}

#ifdef SEEK_I_CONST
		*seekic = mfr_light; // shall be done only once because in "SEEK_I_CONST" mode only irradiance E=0 is used
#endif

// ==============================================================================================================================
		// Preparations for output scripts for firing rate map plots, Gauss sigma fit plots and irradiance map plot (for specific stimulus amplitudes)
		if (considerAmplitude(i)) 
		{
			for (int m=0; m < pow(Nl,2)+pow(Nl_inh,2); m++)
			{
				double fr = 1000. * double(net.getSpikeCount(m)) / t_max;
#ifdef STIMULUS_COMPARISON
				double fr2 = 1000. * double(net2.getSpikeCount(m)) / t_max;
#endif
				
				if (m < pow(Nl,2))	// neuron m is in excitatory population
				{
					if (fr > max_firing_rate)
						max_firing_rate = fr;
					if (fr < min_firing_rate)
						min_firing_rate = fr;
#ifdef STIMULUS_COMPARISON
					if (fr2 > max_firing_rate)
						max_firing_rate = fr2;
					if (fr2 < min_firing_rate)
						min_firing_rate = fr2;
#endif

					txt_fr[current_E] << fixed << rowNN(m) << "\t\t\t" << colNN(m) << "\t\t\t" << fr
#ifdef STIMULUS_COMPARISON
											<< "\t\t\t" << fr2
#endif
											<< endl;

					if (current_E == stimplot.size()-1) // data output for light stimulus plot (for last specific amplitude)
					{
						txt_stim << fixed << rowNN(m) << "\t\t\t" << colNN(m) << "\t\t\t" << net.getLightStimulusAt(offset, m) << endl;
						if ((m+1) % Nl == 0)
							txt_stim << endl; // another free line necessary for contour plot
					}	
					if ((m+1) % Nl == 0)
						txt_fr[current_E] << endl; // another free line necessary for contour plot				
				}
				else // neuron m is in inhibitory population
				{
					if (fr > max_firing_rate_inh)
						max_firing_rate_inh = fr;
					if (fr < min_firing_rate_inh)
						min_firing_rate_inh = fr;
#ifdef STIMULUS_COMPARISON
					if (fr2 > max_firing_rate_inh)
						max_firing_rate_inh = fr2;
					if (fr2 < min_firing_rate_inh)
						min_firing_rate_inh = fr2;
#endif
					int m_eff = m - pow(Nl, 2);
					txt_fr_inh[current_E] << fixed << rowInh(m_eff) << "\t\t\t" << colInh(m_eff) << "\t\t\t" << fr
#ifdef STIMULUS_COMPARISON
												 << "\t\t\t" << fr2
#endif
												 << endl;
					
					if ((m_eff+1) % Nl_inh == 0)
						txt_fr_inh[current_E] << endl; // another free line necessary for contour plot	
				}
			}

			// Close firing rate map data files
			txt_fr[current_E].close(); 
			txt_fr_inh[current_E].close();

			// Create fit plot (Gaussian FR distribution) of excitatory population for light stimulus
			createGaussFitPlot(gpl_fit_light[current_E], par_light[0], sigma_gauss_light, par_light[2], i*dE, center, "_light");
			gplscript << "gnuplot gauss_sigma_fit_" << dtos(i*dE,2) << "_light.gpl" << endl;

#ifdef STIMULUS_COMPARISON 
			// Create fit plot (Gaussian FR distribution) of excitatory population for current stimulus
			createGaussFitPlot(gpl_fit_light[current_E], par_curr[0], sigma_gauss_current, par_curr[2], i*dE, center, "_current");
			gplscript << "gnuplot gauss_sigma_fit_" << dtos(i*dE,2) << "_current.gpl" << endl;
#endif

			// Create irradiance plot for light stimulus
			if (current_E == stimplot.size()-1) 
			{
				txt_stim.close(); // close data file
				createNetworkColorPlot(gpl_light_stim, Nl, stimplot.back(), 3, "irradiance", "", true, "E / mW/mm²", 0., stimplot.back());
				gplscript << "gnuplot irradiance_" << dtos(stimplot.back(),2) << "_map.gpl" << endl;
#ifdef STIMULUS_COMPARISON 
				createNetworkColorPlot(gpl_curr_stim, Nl, stimplot.back(), 4, "currentstim", "", true, "E / mW/mm²", 0., stimplot.back());
				gplscript << "gnuplot currentstim_" << dtos(stimplot.back(),2) << "_map.gpl" << endl;
#endif
			}
			current_E++;
		} 
	} // end of for(i)
	cout << endl;

// ==============================================================================================================================
	// Create output scripts for firing rate map plots (do this here to have the final min_firing_rate and max_firing_rate)
	for (int e=0; e<stimplot.size(); e++)
	{
		// Create firing rate plot of excitatory population for light stimulus
		createNetworkColorPlot(gpl_lfr[e], Nl, stimplot[e], 3, "fr_exc", "_light", true, "{/Symbol n} / Hz", min_firing_rate, max_firing_rate);
		gplscript << "gnuplot fr_exc_" << dtos(stimplot[e],2) << "_map_light.gpl" << endl;

		// Create firing rate plot of inhibitory population for light stimulus
		createNetworkColorPlot(gpl_lfr_inh[e], Nl_inh, stimplot[e], 3, "fr_inh", "_light", true, "{/Symbol n} / Hz", min_firing_rate_inh, max_firing_rate_inh);
		gplscript << "gnuplot fr_inh_" << dtos(stimplot[e],2) << "_map_light.gpl" << endl;

#ifdef STIMULUS_COMPARISON 
		// Create firing rate plot of excitatory population for external current stimulus
		createNetworkColorPlot(gpl_cfr[e], Nl, stimplot[e], 4, "fr_exc", "_current", true, "{/Symbol n} / Hz", min_firing_rate, max_firing_rate);
		gplscript << "gnuplot fr_exc_" << dtos(stimplot[e],2) << "_map_current.gpl" << endl;
	
		// Create firing rate plot of inhibitory population for external current stimulus
		createNetworkColorPlot(gpl_cfr_inh[e], Nl_inh, stimplot[e], 4, "fr_inh", "_current", true, "{/Symbol n} / Hz", min_firing_rate_inh, max_firing_rate_inh);
		gplscript << "gnuplot fr_inh_" << dtos(stimplot[e],2) << "_map_current.gpl" << endl;
#endif

	}

// ==============================================================================================================================
	// Create plot for spikes per millisecond
#ifdef SPIKES_PLOT
	gpl_spikes << "set term pdf enhanced font \"FrontPage, 20\" color solid lw 1" << endl
				  << "set output '" << dateStr("_spike_num.pdf'") << endl
				  << "set xlabel \"time / ms\"" << endl
				  << "set ylabel \"spikes/ms\"" << endl
				  << "set tmargin at screen 0.95" << endl
				  << "set bmargin at screen 0.23" << endl << endl
				  << "plot [x=0:" << t_max << "] '" << dateStr("_spike_num.txt") << "' \\" << endl
				  << 		"\tusing 1:2 notitle with lines";
	txt_spikes.close();
	gpl_spikes.close();
	system("gnuplot spike_num.gpl");
	gplscript << "gnuplot spike_num.gpl" << endl;
#endif

// ==============================================================================================================================
	// Create plots for mean firing rate and Gauss sigma over irradiance/current stimulus
	createFROverStimPlot(gpl_mf, dE, E_max, PLOT_TYPE_MEAN);
	gplscript << "gnuplot mean_fr.gpl" << endl;
	createFROverStimPlot(gpl_sf, dE, E_max, PLOT_TYPE_GAUSS_SIGMA);
	gplscript << "gnuplot gauss_sigma_fr.gpl" << endl;

	// Close files that have remained open
	txt_msf.close();
	gplscript.close(); // script is not executed, it is only created in case it is needed later on

	// Call GNUplot to process the files
	system("gnuplot mean_fr.gpl");
	system("gnuplot gauss_sigma_fr.gpl");

// ==============================================================================================================================
	// Display final information and save parameters 
	int tsec = timeMeasure(false);
	char el_time [32];
	sprintf(el_time, "Elapsed time: %d min %d sec", int(floor(tsec / 60)), int(tsec % 60));
	cout << el_time << endl;
	saveParams(el_time); // additional informaton: elapsed time, threshold frequency for constant illumination
	if (contour != NULL)
		*contour << endl; // attach empty line to contour plot file
	if (system("cd ../..") == -1)
		showChDirErrMessage();
	cout << separator << endl;

// ==============================================================================================================================
	// Free allocated memory
	delete[] txt_fr;
	delete[] txt_fr_inh;
	delete[] gpl_lfr;
	delete[] gpl_lfr_inh;
	delete[] gpl_fit_light;
#ifdef STIMULUS_COMPARISON
	delete[] gpl_cfr;
	delete[] gpl_cfr_inh;
	delete[] gpl_fit_current;
#endif
logf.close();
	return 0;
}

#ifdef SEEK_I_CONST
/*** setSeekICVar ***
 * Sets the variable to communicate with NetworkBatch class for seeking I_const */
void setSeekICVar(double *_seekic) 
{
	seekic = _seekic;
}
#endif

/*** setParams ***
 * Sets the simulation parameters on given values and resets network(s) */
void setParams(double _t_pulse, double _frequency, double _const_current, double _tau_syn, double _J,
					double _J_ee, double _J_ei, double _J_ie, double _J_ii, double _pc) 
{
	t_pulse = _t_pulse;
	frequency = _frequency;
	pc = _pc;
	//net.setConnection(3,5,1/sqrt(2),2,15,20,1,1/sqrt(2),0.4,0.4);
	tau_syn = _tau_syn;
	J = _J;
	net.setConstCurrent(_const_current);
	net.setSynTimeConstant(_tau_syn);
	net.setCouplingStrengths(_J_ee*J, _J_ei*J, _J_ie*J, _J_ii*J);
	net.reset();
	net.setPC(_pc);
#ifdef STIMULUS_COMPARISON
	net2.setConstCurrent(_const_current);
	net2.setSynTimeConstant(_tau_syn);
	net2.setCouplingStrengths(_J_ee*J, _J_ei*J, _J_ie*J, _J_ii*J);
	net2.reset();
	net2.setPC(_pc);
#endif
}

/*** Constructor ***
 * Sets all parameters on given values and calls constructors for Neuron instances *
 * - _contour: ofstream of contour plot data file or NULL */
NetworkSimulation(const int _Nl, const int _Nl_inh, const double _dt, const double _dE, double _t_offset, double _t_max, 
						double _E_max, const vector<double> _stimplot, ofstream *_contour) 
	: Nl(_Nl), Nl_inh(_Nl_inh), net(_dt, _Nl, _Nl_inh), dt(_dt), dE(_dE), center(double(_Nl)/2.+0.5), center_inh(double(_Nl_inh)/2.+0.5), 
	  stimplot(_stimplot), t_offset(_t_offset)
#ifdef STIMULUS_COMPARISON
	  , net2(_dt, _Nl, _Nl_inh)
#endif
{
	contour = _contour;
	t_max = _t_max; 
	E_max = _E_max;
}	

};


















