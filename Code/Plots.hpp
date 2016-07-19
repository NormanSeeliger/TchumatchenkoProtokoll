/************************************************************
 *** Functions creating scripts for plotting with gnuplot ***
 ***       (c) Jannik Luboeinski 2015/2016	       		 ***
 ************************************************************/

/*** writePalParula ***
 * Writes a palette file for gnuplot similar to Matlab's parula colors */
void writePalParula()
{
	ofstream pal("palette.pal");
	pal << "# Matlab color map parula" << endl << endl;
	pal << "set palette defined (\\" << endl;
	pal << " 0    0.2081    0.1663    0.5292,\\" << endl;
	pal << " 1    0.2116    0.1898    0.5777,\\" << endl;
	pal << " 2    0.2123    0.2138    0.6270,\\" << endl;
	pal << " 3    0.2081    0.2386    0.6771,\\" << endl;
	pal << " 4    0.1959    0.2645    0.7279,\\" << endl;
	pal << " 5    0.1707    0.2919    0.7792,\\" << endl;
	pal << " 6    0.1253    0.3242    0.8303,\\" << endl;
	pal << " 7    0.0591    0.3598    0.8683,\\" << endl;
	pal << " 8    0.0117    0.3875    0.8820,\\" << endl;
	pal << " 9    0.0060    0.4086    0.8828,\\" << endl;
	pal << "10    0.0165    0.4266    0.8786,\\" << endl;
	pal << "11    0.0329    0.4430    0.8720,\\" << endl;
	pal << "12    0.0498    0.4586    0.8641,\\" << endl;
	pal << "13    0.0629    0.4737    0.8554,\\" << endl;
	pal << "14    0.0723    0.4887    0.8467,\\" << endl;
	pal << "15    0.0779    0.5040    0.8384,\\" << endl;
	pal << "16    0.0793    0.5200    0.8312,\\" << endl;
	pal << "17    0.0749    0.5375    0.8263,\\" << endl;
	pal << "18    0.0641    0.5570    0.8240,\\" << endl;
	pal << "19    0.0488    0.5772    0.8228,\\" << endl;
	pal << "20    0.0343    0.5966    0.8199,\\" << endl;
	pal << "21    0.0265    0.6137    0.8135,\\" << endl;
	pal << "22    0.0239    0.6287    0.8038,\\" << endl;
	pal << "23    0.0231    0.6418    0.7913,\\" << endl;
	pal << "24    0.0228    0.6535    0.7768,\\" << endl;
	pal << "25    0.0267    0.6642    0.7607,\\" << endl;
	pal << "26    0.0384    0.6743    0.7436,\\" << endl;
	pal << "27    0.0590    0.6838    0.7254,\\" << endl;
	pal << "28    0.0843    0.6928    0.7062,\\" << endl;
	pal << "29    0.1133    0.7015    0.6859,\\" << endl;
	pal << "30    0.1453    0.7098    0.6646,\\" << endl;
	pal << "31    0.1801    0.7177    0.6424,\\" << endl;
	pal << "32    0.2178    0.7250    0.6193,\\" << endl;
	pal << "33    0.2586    0.7317    0.5954,\\" << endl;
	pal << "34    0.3022    0.7376    0.5712,\\" << endl;
	pal << "35    0.3482    0.7424    0.5473,\\" << endl;
	pal << "36    0.3953    0.7459    0.5244,\\" << endl;
	pal << "37    0.4420    0.7481    0.5033,\\" << endl;
	pal << "38    0.4871    0.7491    0.4840,\\" << endl;
	pal << "39    0.5300    0.7491    0.4661,\\" << endl;
	pal << "40    0.5709    0.7485    0.4494,\\" << endl;
	pal << "41    0.6099    0.7473    0.4337,\\" << endl;
	pal << "42    0.6473    0.7456    0.4188,\\" << endl;
	pal << "43    0.6834    0.7435    0.4044,\\" << endl;
	pal << "44    0.7184    0.7411    0.3905,\\" << endl;
	pal << "45    0.7525    0.7384    0.3768,\\" << endl;
	pal << "46    0.7858    0.7356    0.3633,\\" << endl;
	pal << "47    0.8185    0.7327    0.3498,\\" << endl;
	pal << "48    0.8507    0.7299    0.3360,\\" << endl;
	pal << "49    0.8824    0.7274    0.3217,\\" << endl;
	pal << "50    0.9139    0.7258    0.3063,\\" << endl;
	pal << "51    0.9450    0.7261    0.2886,\\" << endl;
	pal << "52    0.9739    0.7314    0.2666,\\" << endl;
	pal << "53    0.9938    0.7455    0.2403,\\" << endl;
	pal << "54    0.9990    0.7653    0.2164,\\" << endl;
	pal << "55    0.9955    0.7861    0.1967,\\" << endl;
	pal << "56    0.9880    0.8066    0.1794,\\" << endl;
	pal << "57    0.9789    0.8271    0.1633,\\" << endl;
	pal << "58    0.9697    0.8481    0.1475,\\" << endl;
	pal << "59    0.9626    0.8705    0.1309,\\" << endl;
	pal << "60    0.9589    0.8949    0.1132,\\" << endl;
	pal << "61    0.9598    0.9218    0.0948,\\" << endl;
	pal << "62    0.9661    0.9514    0.0755,\\" << endl;
	pal << "63    0.9763    0.9831    0.0538)\\" << endl;
	pal.close();
}

/*** writePalRainbow ***
 * Writes a palette file for gnuplot with the color that were used for the first single channel *
 * color plots */
void writePalRainbow()
{
	ofstream pal("palette.pal");
	pal << "# Rainbow color map" << endl << endl;
	pal << "set palette defined (\\" << endl;
	pal << " 0    \"#0000FF\",\\" << endl;
	pal << " 1    \"#0000FF\",\\" << endl;
	pal << " 2    \"#0080FF\",\\" << endl;
	pal << " 3    \"#00FFFF\",\\" << endl;
	pal << " 4    \"#00FF40\",\\" << endl;
	pal << " 5    \"#80FF00\",\\" << endl;
	pal << " 6    \"#FFFF00\",\\" << endl;
	pal << " 7    \"#FFBF00\",\\" << endl;
	pal << " 8    \"#FF4000\",\\" << endl;
	pal << " 9    \"#FF0000\",\\" << endl;
	pal.close();
}

/*** createOverviewColorPlot ***
 * Creates a color plot of data that has to be specified, representing an overview over different parameter sets *
 * - dE: size of one irradiance step
 * - E_min: minimum irradiance
 * - E_max: maximum irradiance
 * - file: a partial file name for the gnuplot script and the resulting plot file
 * - ylabel: label of the y-axis
 * - ystart: lowest y-value
 * - yend: highest y-value
 * - ystep: one step of the y-axis values
 * - yformat: format (decimal places) of the y-axis tics
 * - cblabel: label of the colorbar
 * - xcolumn: column in the data file in which the x values are located
 * - ycolumn: column in the data file in which the y values are located
 * - cbcolumn: column in the data file in which the color values are located
 * - threedim [optional]: three-dimensional plot (true/false) */
void createOverviewColorPlot(double dE, double E_min, double E_max, const char* file, const char *ylabel, double ystart, double yend, double ystep, 
								     const char* yformat, const char *cblabel, int xcolumn, int ycolumn, int cbcolumn, bool threedim = false)
{
	ofstream gpl(concat(file, "_cplot.gpl"));
	gpl << "# Output configuration" << endl;
	gpl << "set term pdf enhanced font \"FrontPage, 20\" color solid lw 2.5 size 5,4" << endl;
	gpl << "set output '" << dateStr("_cplot_") << file << ".pdf'" << endl;
	gpl << "#set term png enhanced font \"FrontPage, 16\" size 1280,1024" << endl;
	gpl << "#set output '" << dateStr("_cplot_") << file << ".png'" << endl << endl;

	gpl << "# Color map plot configuration" << endl;
	if (!threedim)
	{
		gpl << "set view map	# use 2-dim. plot instead of a 3-dim. plot" << endl;
	}
	else 
	{
		gpl << "set view 70,75 # just for 3-dim. plot to rotate" << endl;
		gpl << "set ticslevel 0 # just for 3-dim. plot to shift z=0 down into the xy plane" << endl;
	}
	gpl << "unset key" << endl;
	gpl << "load 'palette.pal'" << endl << endl;

	gpl << "# Axis configuration" << endl;
	gpl << "xmin=" << E_min << endl;
	gpl << "xmax=" << E_max << endl;
	gpl << "ymin=" << ystart << endl;
	gpl << "ymax=" << yend << endl;
	gpl << "set xtics " << dE << endl;
	gpl << "set ytics " << 2.*ystep << endl;
	gpl << "set mytics 2" << endl;
	gpl << "set xrange [xmin-" << 0.5*dE << ":xmax+" << 0.5*dE << "]" << endl; // extend for half the step size for fields in matrix
	gpl << "set yrange [ymin-" << 0.5*ystep << ":ymax+" << 0.5*ystep << "]" << endl;
	gpl << "set xlabel \"E / mW/mm²\"" << endl;
	gpl << "set ylabel \"" << ylabel << "\" offset -1,0" << endl;
	gpl << "set cblabel \"" << cblabel << "\" offset 1" << endl;
	gpl << "set format x '%.1f'" << endl;
	gpl << "set format y " << yformat << endl;
	gpl << "set format z '%.0f'" << endl << endl;

	gpl << "# Insert grid" << endl;
	gpl << "set for [i=0:" << int(round((E_max-E_min)/dE)) << "] arrow from (xmin+" << 0.5*dE << ")+i*" << dE << ",ymin-" << 0.5*ystep
		 << 	" to (xmin+" << 0.5*dE << ")+i*" << dE << ",ymax+" << 0.5*ystep << " nohead front lw 1 lc rgb 'black' lw 1" << endl;
	gpl << "set for [i=0:" << int(round((yend-ystart)/ystep))-1 << "] arrow from xmin-" << 0.5*dE << ",(ymin+" << 0.5*ystep << ")+i*" << ystep
		 <<   " to xmax+" << 0.5*dE << ",(ymin+" << 0.5*ystep << ")+i*" << ystep << " nohead front lw 1 lc rgb 'black' lw 1" << endl << endl;

	
	gpl << "# Set margins" << endl;
	gpl << "set tmargin at screen 0.97" << endl;
	gpl << "set bmargin at screen 0.20" << endl;
	gpl << "set rmargin at screen 0.75" << endl;
	gpl << "set lmargin at screen 0.21" << endl << endl;

	gpl << "# Plotting" << endl;
	gpl << "set datafile missing \"nan\"" << endl;
	gpl << "splot '" << dateStr("_cplotdata.txt'") << " using " << xcolumn << ":" << ycolumn << ":" << cbcolumn << " notitle with image" << endl;
	gpl.close();

	//system(concat("gnuplot ", concat(file, "_cplot.gpl")).c_str());
}

/*** createNetworkColorPlot ***
 * Creates a network color plot of data that has to be specified *
 * - f: the gnuplot file to write to
 * - Nl: the number of neurons in one line
 * - irradiance: the irradiance (amplitude) used for the simulation (set this to a negative value if irrelevant)
 * - column: the column in the data file to be used for color-coded values
 * - prefix: the principal name of the plot
 * - postfix: if necessary, the type of stimulation (actually either "_light" or "_current") 
 * - matrix: if true, a non-interpolated matrix plot is created, else an interpolated color plot
 * - cblabel: the label for the colorbar 
 * - cbmin [optional]: the minimum for the colorbar
 * - cbmax [optional]: the maximum for the colorbar 
 * - cbtics [optional]: one increment for the colorbar axis */
void createNetworkColorPlot(ofstream &f, int Nl, double irradiance, int column, const char* prefix, const char* postfix, bool matrix, const char* cblabel, 
							double cbmin = -1.0, double cbmax = -1.0, double cbtics = -1.0)
{
	char gplcall[100];

	if (irradiance >= 0.0)
		sprintf(gplcall, "%s_%.2f", prefix, irradiance);
	else	
		sprintf(gplcall, "%s", prefix);

	f << "# Output configuration" << endl;
	f << "set term pdf enhanced font \"FrontPage, 20\" color solid lw 2.5 size 5,4.5" << endl;
	f << "set output '" << dateStr("_") << gplcall << "_map" << postfix << ".pdf'" << endl << endl;

	f << "# Contour plot configuration" << endl;
	if (!matrix)
	{
		f << "set pm3d" << endl;
		f << "unset surface" << endl;
	}
	f << "set view map	# 2-dim. plot" << endl;
	f << "#set view 50,75 # just for 3-dim. plot to rotate # 3-dim. plot" << endl;
	f << "#set ticslevel 0 # just for 3-dim. plot to shift z=0 down into the xy plane # 3-dim. plot" << endl;
	f << "unset key" << endl;
	if (!matrix)
		f << "set pm3d interpolate 100,100 # interpolate the color (0: gnuplot automatically decides the steps to use)" << endl;
	f << "load 'palette.pal' # load palette" << endl;
	f.precision(5); // set number of significant digits in output
	if (strcmp(prefix, "firingrate") == 0)
		f << "set title \"Firing rates across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "irradiance") == 0)
		f << "set title \"Light intensities across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "inc_connection") == 0)
		f << "set title \"Incoming connections per neuron" << endl;
	f << "set size square" << endl << endl;  // quadratic size

	f << "# Axis configuration" << endl;
	if (!matrix)
	{
		f << "set xrange [1" << ":" << Nl << "]" << endl;
		f << "set yrange [1" << ":" << Nl << "]" << endl;
	}
	else
	{
		f << "set xrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
		f << "set yrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
	}
	f << "set mxtics 2" << endl;
	f << "set mytics 2" << endl;
	if (cbmin >= 0.0 && cbmax >= 0.0)
	{
		if (cbmin == cbmax)
			f << "set cbrange [" << 0.95*cbmin << ":" << 1.05*cbmin << "]" << endl;
		else
			f << "set cbrange [" << cbmin << ":" << cbmax << "]" << endl;
	}
	if (cbtics > 0.0)
		f << "set cbtics " << cbtics << endl;
	f << "set xlabel \"Neuron column\"" << endl;
	f << "set ylabel \"Neuron row\" offset -1.5" << endl; // offset to shift the label to the left
	f << "set cblabel \"" << cblabel << "\"" << endl;
	f << "set format x '%1.0f'" << endl;
	f << "set format y '%1.0f'" << endl;
	f << "set format z '%.5f'" << endl << endl;

	f << "# Set margins" << endl;
	f << "set tmargin at screen 0.90" << endl;
	f << "set bmargin at screen 0.17" << endl;
	f << "set rmargin at screen 0.78" << endl;
	f << "set lmargin at screen 0.15" << endl << endl;

	f << "# Plotting" << endl;
	f << "set datafile missing \"nan\"" << endl;
	if (!matrix)
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " notitle with lines lt 1" << endl;
	else
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " with image" << endl;
	f.close();

	if (irradiance >= 0.0)
		sprintf(gplcall, "gnuplot %s_%.2f_map%s.gpl", prefix, irradiance, postfix);
	else
		sprintf(gplcall, "gnuplot %s_map%s.gpl", prefix, postfix);
	system(gplcall);
}


/*** createFROverStimPlot ***
 * Creates a simple plot of either the mean firing rate or the firing rate standard deviation over the stimulus quantity *
 * - f: the gnuplot file to write to
 * - dE: size of one irradiance step
 * - E_max: maximum irradiance
 * - type: type of plot (types defined below by column in data file) */
#define PLOT_TYPE_MEAN 2
#define PLOT_TYPE_GAUSS_SIGMA 4
void createFROverStimPlot(ofstream &f, const double dE, const double E_max, const int type)
{
	f << "set term pdf enhanced font \"FrontPage, 20\" color solid lw 2.5" << endl; // solid
	f << "set output '";
	if (type == PLOT_TYPE_MEAN)
		f << dateStr("_mean_fr.pdf'") << endl;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << dateStr("_gauss_sigma_fr.pdf'") << endl;
	//f << "set log x" << endl;
	f << "set xlabel \"Ê / mW/mm^2\"" << endl;
#ifdef STIMULUS_COMPARISON
	f << "set x2label \"Î_{cst} / nA\"" << endl;
	f << "set key top left" << endl;
	f << "set key samplen 2" << endl;
	//f << "set log x2" << endl;
	f << "set x2range [" << 0.1 << ":" << 1 << "]" << endl; //TODO
	f << "set xtics nomirror" << endl;
	f << "set x2tics" << endl;
#endif
	if (type == PLOT_TYPE_MEAN)
		f << "set ylabel \"⟨{/Symbol n}⟩ / Hz\"" << endl;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << "set ylabel \"{/Symbol s}_{FR}\"" << endl;
	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	f << "plot [x=";
	if (type == PLOT_TYPE_MEAN)
		f << -0.5*dE;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << +0.5*dE; // no value exists for E=0

	f << ":" << E_max+0.5*dE << "] '" << dateStr("_fr.txt'") 
	  //<< " using ($1 == 0 ? NaN : $1):" // leave out irradiance 0.00 (replace it with NaN)
	  << " using 1:" << type << ":" << type+1 << " lc 1 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 1:" << type << " lc 1 " // lines
	  << "title \"exc.\" with lines," << " \\" << endl
	  << "\t\t'" << dateStr("_fr.txt'") 
	  << " using 1:" << type+4 << ":" << type+5 << " lc 2 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 1:" << type+4 << " lc 2 " // lines
	  << "title \"inh.\" with lines";

#ifdef STIMULUS_COMPARISON
	f << "," << " \\" << endl
	  << "\t\t'" << dateStr("_fr.txt'") 
	  //<< " using ($6 == 0 ? NaN : $6):"
	  << " using 10:" << type+9 << ":" << type+10 << " lc 3 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'") 
	  << " using 10:" << type+9 << " lc 3 " // lines
	  << "title \"cst\" axes x2y1 with lines" << endl;
#endif
	
	f.close();
}


/*** createPeakPlot ***
 * Creates a simple plot of an average peak of either the instantaneous firing rate or the open probability over time *
 * - f: the gnuplot file to write to
 * - period_len: length of one stimulus period in ms
 * - E: irradiance for this peak plot in mW/mm²
 * - type: type of plot (types defined below by column in data file) 
 * - amp: absolute amplitude of the peak
 * - mu: the "mean" of the peak (the time where the amplitude is located)
 * - base: the baseline (lowest value within the peak)
 * - t_hm_left: the time in ms at which lies the left "half maximum"
 * - t_hm_right: the time in ms at which lies the right "half maximum"
 [not needed anymore (were necessary for Gauss fit):]
 * - freq: stimulus frequency for this peak plot in Hz
 * - cplotfile: file that contains data for color plot of several measurements */
#define PLOT_TYPE_FR 2
#define PLOT_TYPE_O 3
void createPeakPlot(ofstream &f, double period_len, double E, int type, double amp, double mu, double base, double t_hm_left, double t_hm_right/*, int freq, string cplotfile*/)
{
	f << "set term pdf enhanced font \"FrontPage, 20\" color lw 2.5" << endl;

	// Set output file
	f << "set output '";
	if (type == PLOT_TYPE_FR)
		f << dateStr(string("_fr_peak_") + dtos(E,2) + string(".pdf'")) << endl << endl;
	else if (type == PLOT_TYPE_O)
		f << dateStr(string("_O_peak_") + dtos(E,2) + string(".pdf'")) << endl << endl;

	
	// Set variables (or rather initial values in case of Gauss fit)
	f << "A = " << amp-base << endl 
	  << "mu = " << mu << endl 
	  << "b = " << base << endl
	  << "FWHM = " << t_hm_right-t_hm_left << endl << endl;
	
	// Set labels
	f << "set xlabel \"t / ms\"" << endl;
	if (type == PLOT_TYPE_O)
		f << "set ylabel \"O(t)\" offset 1.8" << endl << endl;
	else if (type == PLOT_TYPE_FR)
	{
		f << "set ylabel \"{/Symbol n}_{inst} / Hz\" offset 1.8" << endl << endl;
		/* Gauss fit to obtain FWHM (not used anymore because peaks are in general not symmetric)
		f << "# Fit function" << endl 
		  << "set fit quiet" << endl // to suppress fit output in the terminal (fit.log is still created)
		  << "set fit errorvariables" << endl // creates *_err type variables for fitting errors
		  << "sigma = 2" << endl
		  << "g(x) = A * exp(- (x - mu)**2 / (2*sigma**2)) + b" << endl
		  << "fit g(x) '" << dateStr("_peak_") << dtos(E,2) << ".txt' using 1:2 via A, mu, sigma, b" << endl
		  << "FWHM = 2*sqrt(2*log(2))*sigma" << endl
		  << "FWHM_err = 2*sqrt(2*log(2))*sigma_err" << endl << endl;

		f << "# Write FWHM result to colorplot data file" << endl
		  << "set print \"../" << cplotfile << "\" append" << endl
		  << "print " << dtos(freq,1) << ", " << dtos(E,2) << ", A+b, sqrt(A_err**2 + b_err**2), FWHM, FWHM_err" << endl << endl;*/
	}
	// Print parameters obtained from fitting to labels and to file
	f << "# Create labels and lines" << endl 
	  << "set label sprintf(\"Height = %.2f";
	if (type == PLOT_TYPE_FR)
		f << " Hz";
	f << "\", A+b) at screen 0.56,screen 0.88" << endl 
	  << "set label sprintf(\"Base = %.2f";
	if (type == PLOT_TYPE_FR)
		f << " Hz";
	f << "\", b) at screen 0.56,screen 0.81 " << endl
	  << "set label sprintf(\"FWHM = %.1f ms\", FWHM) at screen 0.56,screen 0.74 " << endl
	  << "set arrow from " << t_hm_left << ",b+A/2 to " << t_hm_right << ",b+A/2 nohead lc rgb 'blue' dt 2" << endl << endl;

	// Set margins
	f << "# Set margins" << endl ;
	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	// Plot command
	f << "plot [x=0" << ":" << period_len << "] '" << dateStr("_peak_") << dtos(E,2) << ".txt' "
	  << "using 1:" << type << " notitle lc rgb 'red' with lines";
	/* Gauss fit
	if (type == PLOT_TYPE_FR)
		f << ", \\" << endl << "\t\tg(x) notitle lc rgb 'green'" << endl;*/
	
	f.close();

	// Call gnuplot script
	if (type == PLOT_TYPE_FR)
		system(concat("gnuplot fr_peak_", dtos(E,2) + string(".gpl")).c_str());
	else if (type == PLOT_TYPE_O)
		system(concat("gnuplot O_peak_", dtos(E,2) + string(".gpl")).c_str());
}

/*** createGaussFitPlot ***
 * Creates a plot of data and fit curve of a half Gaussian fit (firing rate over radius from center), if GPL_FIT_COMPARISON is defined, *
 * a comparative gnuplot fit is also included in the script *
 * - f: the gnuplot file to write to
 * - A: the amplitude of the Gaussian obtained by a fit
 * - sigma: the standard deviation of the Gaussian obtained by a fit
 * - b: the offset of the Gaussian obtained by a fit
 * - E: the currently adjusted irradiance
 * - center: the center of the quadratic neuron population
 * - postfix: posterior term for the filename */
void createGaussFitPlot(ofstream &f, const double A, const double sigma, const double b, const double E, const double center, string postfix)
{
	f << "set term pdf enhanced font \"FrontPage, 20\" color lw 1" << endl;
	f << "set output '" << dateStr("_gauss_sigma_fit_") << dtos(E,2) << postfix << ".pdf'" << endl << endl;

	f << "set xlabel \"r\"" << endl;
	f << "set ylabel \"{/Symbol n} / Hz\"" << endl;
	f << "set key samplen 1" << endl << endl;

	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	f << "f(x) = " << A << " * exp( -x**2 / (2*" << sigma << "**2) ) + " << b << endl << endl;

#ifdef GPL_FIT_COMPARISON
	f << "g(x) = A * exp( -x**2 / (2*sigma**2) )" << endl;
	f << "A = " << A << endl;
	f << "sigma = " << sigma << endl;
	f << "fit g(x) '" << dateStr("_fr_exc_") << dtos(E,2) << ".txt' using (sqrt(($1-" << center << ")**2 + ($2-" << center << ")**2)):3 via A, sigma" << endl << endl;
#endif

	f << "plot '" << dateStr("_fr_exc_") << dtos(E,2) << ".txt' using (sqrt(($1-" << center << ")**2 + ($2-" << center << ")**2)):3 lc 1 title 'data', \\" << endl
		 << "\t\tf(x) lc 4 lw 2 title 'GSL fit'";
#ifdef GPL_FIT_COMPARISON
	f << ", g(x) lc 2 lw 2 dt 2 title 'GPL fit'";
#endif
	f.close();

	system(concat("gnuplot gauss_sigma_fit_" + dtos(E,2) + postfix, ".gpl").c_str());
}



