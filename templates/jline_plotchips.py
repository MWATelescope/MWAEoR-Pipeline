#!/usr/bin/env python

"""
Taken from /astro/mwaeor/jline/software/chips_2019/scripts/plotchips_all.py

placeholder - this does CHIPS plots

Tsys is something that is hard coded inside CHIPS. If sub-banding, this needs
to be accounted for, and this is currently **not** done in this script. If you
care about this you'll have to do some cosmologically based Tsys calculations
which will effect the thermal noise estimation (not the power, but the upper
limits as they are a combo of power and noise)

example:

```
salloc --nodes=1 --mem=350G --time=8:00:00 --clusters=garrawarla --partition=workq --account=mwaeor --ntasks=20 --tmp=500G
module load singularity
mkdir /dev/shm/deleteme
cd /dev/shm/deleteme
cp /astro/mwaeor/dev/nfdata/1094751504/ps_metrics-newhyp+qa/* .
export ext="grid_30l_src4k_8s_80kHz_eor0high"
singularity exec --cleanenv --home /astro/mwaeor/dev/mplhome /pawsey/mwa/singularity/ssins/ssins_latest.sif python \
    /pawsey/mwa/mwaeor/dev/MWAEoR-Pipeline/templates/jline_plotchips.py \
    --basedir=./ \
    --polarisation=xx \
    --chips_tag=grid_30l_src4k_8s_80kHz_eor0high \
    --min_power=1e3 \
    --max_power=1e15
cp chips2D_xx_grid_30l_src4k_8s_80kHz_eor0high_crosspower.png /astro/mwaeor/dev/nfdata/1094751504/ps_metrics-newhyp+qa/
```
"""

import numpy as np
from numpy import ma
import os
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib import ticker,cm
from matplotlib.colors import LogNorm, Normalize
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import warnings
from astropy.cosmology import LambdaCDM
from copy import deepcopy

##Speed of light in km/s
VELC_KMS = 2.995e5

##Speed of light in m/s
SPEED_LIGHT = 299792458.0

##21cm wavelength in m
WAVELENGTH_21CM = 0.21

##Boltzmann constant
BOLTZMANN = 1.38e-23

##Lower k perp bins that always get ignored from outputs
KPERP_START = 2

##Upper k_parallel bins that always get ignored from outputs
KPARRA_END = 2

def get_args(argv=None):
    """Parse command line arugments using argparse. Returns the args"""

    parser = argparse.ArgumentParser()
    parser.add_argument("--title", default="Crosspower")

    file_group = parser.add_argument_group('FILE LOCATIONS')
    file_group.add_argument("--basedir", default='/astro/mwaeor/MWA/output/',
        help="The directory where the CHIPS outputs are located")
    file_group.add_argument("--chips_tag",
        help="Use this when doing a 2D power spectrum (not a ratio). "\
        "This is the tag that was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag=real_cool_name")

    file_group.add_argument("--chips_tag_one",
        help="Use this data as the numerator in a 2D ratio/diff plot. " \
        "This is the string that was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag_one=real_cool_name")
    file_group.add_argument("--chips_tag_two",
        help="Use this data as the denominator in a 2D ratio/diff plot. " \
        "This is thethat was passed to the 'lssa_thermal' step, e.g. if you "\
        "have an output called crosspower_xx_0.iter.real_cool_name.dat then "\
        "you should enter --chips_tag_two=real_cool_name")
    file_group.add_argument("--outputdir", default='./',
        help="Directory in which to output resultant plots. Default = './'")

    plot_group = parser.add_argument_group('PLOTTING OPTIONS')
    plot_group.add_argument("--plot_type",
        help='Which type of plot to make. Options are: 1D, 2D, 2D_Ratio, '\
              '2D_diff, 1D_comp. Defaults to 2D',default='2D')
    plot_group.add_argument("--polarisation", default='both',
        help='Which polarisation (XX or YY) to plot. Options are: xx, yy, both. Defaults to both')
    plot_group.add_argument("--plot_mode", default='png',
        help="Either 'png' or 'screen'. 'png' saves to .png, 'screen' does plt.show(). " \
        "Default = png")
    plot_group.add_argument("--max_power", default=0.0, type=float,
        help="Maximum power value to show. Defaults to 1e12 for a 2D or 1D" \
             "plot, and 10 for a 2D_Ratio plot. See --max_neg_power if" \
             "doing a 2D difference plot or using --colourscale=pos_and_negs")
    plot_group.add_argument("--min_power", default=0.0, type=float,
        help="Minimum power value to show. Defaults to 1e3 for a 2D or 1D" \
             "plot, and 0.1 for a 2D_Ratio plot. See --min_neg_power if" \
             "doing a 2D difference plot or using --colourscale=pos_and_negs")

    plot_group_2D = parser.add_argument_group('2D PLOTTING OPTIONS')
    plot_group_2D.add_argument("--colourscale", default='all_positive',
        help='How to handle negative powers in the colour scale for 2D plot.' \
        "Choices are: negs_are_grey (anything negative is blocked out as grey), " \
        "all_positive (anything below --min_power is same colour), "\
        "pos_and_negs (positive values are orange, negative purple)" \
        "Defaults to --colourscale=all_positive")
    plot_group_2D.add_argument("--max_neg_power", default=False, type=float,
        help="When plotting with --colourscale=pos_and_negs, use this value" \
        "as the upper limit of negative colourbar. Should be larger than "\
        "--min_neg_power. Default = negative --min_power")
    plot_group_2D.add_argument("--min_neg_power", default=False, type=float,
        help="When plotting with --colourscale=pos_and_negs, use this value" \
        "as the upper limit of negative colourbar. Should be smaller than "\
        "--max_neg_power. Default = negative --max_power")

    plot_group_1D = parser.add_argument_group('1D PLOTTING OPTIONS')
    plot_group_1D.add_argument("--wedge_factor", default=False, type=float,
        help="The scaling factor between k_parrallel and k_perpendicular" \
             "to use during cut. Defaults to the horizon")
    plot_group_1D.add_argument("--plot_delta", default=False, action='store_true',
        help="Plot in unitless Delta rather than P(k)")

    plot_group_1D.add_argument("--ktot_bin_edges", default=False,
        help="Path to text file containing edge values of k-bins for a 1D " \
             "plot. Overrides --low_k_edge, --high_k_edge, --num_k_edges." \
             "Will grid the data to the centre of each pair of bin edges, " \
             "e.g. if you provide 21 bin edges, your data will be gridded " \
             "to 20 bin centres." )

    plot_group_1D.add_argument("--low_k_edge", type=float, default=1e-2,
        help="Lowest k-value to grid the data to (applies to both perp and parra)")
    plot_group_1D.add_argument("--high_k_edge", type=float, default=5,
        help="Highest k-value to grid the data to (applies to both perp and parra)")
    plot_group_1D.add_argument("--num_k_edges", type=int, default=21,
        help="Number of k-bins to grid to between --low_k_edge and --high_k_edge")

    plot_group_1D.add_argument("--kperp_max", default=20, type=float,
        help="Maximum k_perp to use in 1D averaging. Can use this to "\
        "chop off small spatial scales. Defaults to 20")
    plot_group_1D.add_argument("--kperp_min", default=0, type=float,
        help="Minimum k_perp to use in 1D averaging. Can use this to "\
        "chop off large spatial scales. Defaults to no cut")

    plot_group_1D.add_argument("--kparra_min", default=0, type=float,
        help="Minimum k_parra to use in 1D averaging. Can use this to "\
        "chop off things close to the wedge that are leaking")

    plot_group_1D.add_argument("--plot_wedge_cut_2D", default=False,
        action='store_true',
        help="TAdd to plot the 2D power spectra with and without the cuts " \
             "that are being applied before binning to a 1D power spectra")

    plot_group_1D.add_argument("--chips_tag_one_label", default=False,
        help="When plotting two 1D power spectra on same axes, use this to " \
             "label the the data from --chips_tag_one")
    plot_group_1D.add_argument("--chips_tag_two_label", default=False,
        help="When plotting two 1D power spectra on same axes, use this to " \
             "label the the data from --chips_tag_two")




    chips_group = parser.add_argument_group('CHIPS OPTIONS')
    chips_group.add_argument("--N_kperp",type=int, default=80,
        help="The number of kperp bins used in CHIPS 'fft_thermal' command. Default=80")
    chips_group.add_argument("--N_chan",type=int, default=384,
        help="The number of frequency channels used in CHIPS 'prepare_diff' and "\
        "'fft_thermal' commands. Default=384")
    ##TODO make options for low/high band that grab the frequency automagically for you
    chips_group.add_argument("--lowerfreq", default=167.035e6, type=float,
        help="Lowest frequency channel in data (Hz). Default is 167.035e+6")
    chips_group.add_argument("--chan_width", default=80e+3, type=float,
        help="Width of individual spectral channel (Hz). Default = 80e+3")
    # chips_group.add_argument("--deltat", default=8., type=float,
    #     help="Time resolution of data (s). Default = 8")
    chips_group.add_argument("--umax", default=300., type=float,
        help="Maximum u-value used in 'fft_thermal' stage (wavelengths). Default = 300")
    chips_group.add_argument("--density_correction", default=True, action='store_false',
        help="Density correction based on Barry et al. 2019a (Appendix A)")

    plot_group_1D.add_argument("--non_RTS_outputs", default=False,
        action='store_true',
        help="This script defaults to scaling the output power spectra " \
        "based on the RTS weighting scheme. If using data with a different " \
        "weighting scheme, add this flag to switch off this scaling")

    chips_group = parser.add_argument_group('COSMOLOGICAL CONSTANTS')
    chips_group.add_argument("--omega_matter", default=0.272, type=float,
        help="Matter density parameter. Default = 0.272")
    chips_group.add_argument("--omega_baryon", default=0.046, type=float,
        help="Baryon density parameter. Default = 0.046")
    chips_group.add_argument("--omega_lambda", default=0.7, type=float,
        help="Dark energy dentisty parameter. Default = 0.7")
    chips_group.add_argument("--hubble" ,default=70.4, type=float,
        help="Hubble param in km/s/Mpc, default=70.4")

    args = parser.parse_args(argv)

    args.Neta = int(args.N_chan/2)

    ##If people like the caps, let them eat cake
    if args.polarisation == 'XX': args.polarisation = 'xx'
    if args.polarisation == 'YY': args.polarisation = 'yy'

    return args

class FakeArgs(object):
    """When not using argparse, use this object instead which can be called
    from jupyter notebook"""

    def __init__(self):
        self.N_kperp = None
        self.N_chan = None
        self.lowerfreq  = None
        self.chan_width = None
        self.umax = None
        self.omega_matter = None
        self.omega_baryon = None
        self.omega_lambda = None
        self.hubble = None
        self.Neta = None
        self.wedge_factor = 0.0
        self.basedir = None
        self.non_RTS_outputs = False
        self.max_1D_kmode = 20
        self.kperp_max = None
        self.density_correction = False
        self.plot_wedge_cut_2D = False

class ChipsDataProducts(object):
    """Class to read in a bunch of cosmological constants and obserational
    parameters, read in CHIPS data products, and apply them to the CHIPS
    outputs to format them into a 1D or 2D array for plotting"""

    def __init__(self, parser_args):
        """Setup CHIPS cosmological coords and obserational setup based on
        the arguments provided by the parser from argparse"""

        self.parser_args = parser_args
        self.setup_chips_params(parser_args)

    def _create_k_coords(self, u_arr, v_arr, eta, z):
        """
        Convert u distance and eta into k_perp, k_parallel as per Morales et al. (2004).
        Constructs a LambdaCDM cosmology using the values from self.parser_args.

        Stores some cosmological constants for later calculations, being
        self.DM, self.Ez, self.hubble_distance.

        Parameters
        ----------
        u_arr : numpy array, float
            NDarray of u values. Should be in wavelengths.
        eta : numpy array, float
            1Darray of eta values.
        z : float
            Redshift at the central frequency of the band.

        Returns
        -------
        k_parra : numpy array, float
            NDarray of kparra values. Should be in units of h*Mpc^-1.
        k_perp : numpy array, float
            NDarray of kperp values. Should be in units of h*Mpc^-1.
        k_z : numpy array, float
            1Darray of kz values. Should be in units of h*Mpc^-1.
        """
        ##TODO allow the user to select a default astropy cosmology like Planck18
        # from astropy.cosmology import Planck18
        # cosmology = Planck18

        ##Create the cosmology and report the proprties
        cosmology = LambdaCDM(H0=self.parser_args.hubble,
                              Om0=self.parser_args.omega_matter,
                              Ode0=self.parser_args.omega_lambda,
                              Ob0=self.parser_args.omega_baryon)

        print('Cosmology being used has the following parameters:')
        print(f"{chr(11)}H0 = {cosmology.H0:.2f}")
        print(f"{chr(11)}Omega Matter = {cosmology.Om0:.4f}")
        print(f"{chr(11)}Omega Lambda = {cosmology.Ode0:.4f}")
        print(f"{chr(11)}Omega Baryon = {cosmology.Ob0:.4f}")

        nu_21 = (SPEED_LIGHT)/(WAVELENGTH_21CM) #[Hz]

        # Cosmological scaling parameter:
        h = cosmology.H(0).value/100 # Hubble parameter.

        E_z = cosmology.efunc(z) ## Scaling function, see (Hogg 2000)

        # Cosmological distances:
        Dm = cosmology.comoving_distance(z).value*h #[Mpc/h] Transverse co-moving distance.
        DH = 3000 # [Mpc/h] Hubble distance.

        # k parallel
        k_parra = eta * (2*np.pi*nu_21*E_z)/(DH*(1 + z)**2) # [Mpc^-1 h]

        # k perpendicular
        k_perp = u_arr * (2*np.pi/Dm) # [Mpc^-1 h]

        self.DM = Dm
        self.Ez = E_z
        self.h = h
        self.hubble_distance = DH

        return k_parra, k_perp

    def setup_chips_params(self, parser_args):
        """When converting data from output CHIPS files into a 1D or 2D power
        spectrum, do all the the things that are common to the two operations,
        like setting up coordinates and DMs and things.

        There are many cosmological things here that I do not understand"""

        self.central_freq = parser_args.lowerfreq + int(parser_args.N_chan / 2)*parser_args.chan_width

        ##Frequency bandwidth of the data
        bandwidth = float(parser_args.N_chan)*parser_args.chan_width

        ##The eta coords (FT pair of frequency)
        self.eta = np.zeros(int(parser_args.Neta))
        for i in range (0, int(parser_args.Neta)):
            self.eta[i] = (float(i)-0.5) / bandwidth
        self.eta[0] = self.eta[1]/2.

        #The bin length on the u,v plane the data were gridded to (in wavelengths)
        self.u_arr = np.zeros(parser_args.N_kperp)
        for i in range(0,parser_args.N_kperp):
            self.u_arr[i] = float(i)*parser_args.umax*1.1/float(parser_args.N_kperp)
        self.u_arr[0] = self.u_arr[1]/2.

        ##21cm radiation frequency in m/s
        f21 = SPEED_LIGHT / WAVELENGTH_21CM

        ##Redshift based on central frequency
        z = f21 / self.central_freq - 1.

        ##Calculate the k_parallel and k_perpendicular coords
        ##also sets:
        ##    self.DM ( Transverse co-moving distance, [Mpc/h])
        ##    self.Ez (scaling function from Hogg et al 2000)

        ## based on astropy cosmology
        k_parra, k_perp = self._create_k_coords(self.u_arr, self.u_arr, self.eta, z)

        ##There are parts of the 2D outputs that always get thrown away,
        ##so make the same cut on the output coords
        self.kpa = k_parra
        self.kper = k_perp

        ##There are parts of the 2D outputs that always get thrown away,
        ##so make the same cut on the output coords
        self.kpa = self.kpa[:-KPARRA_END]
        self.kper = self.kper[KPERP_START:]


        ##New way of doing it===================================================
        cent_wavelength = SPEED_LIGHT / self.central_freq
        beam_area_steradian = 0.07597
        ##Frequency bandwidth of the data
        bandwidth = float(parser_args.N_chan)*parser_args.chan_width

        norm_numer = cent_wavelength**4 * self.DM**2*self.hubble_distance * bandwidth * (1 + z)**2 * 1.e6
        norm_denom = parser_args.N_chan*(2*BOLTZMANN*1e26)**2*beam_area_steradian*f21*self.Ez

        self.normalisation = norm_numer / norm_denom

        ##Gridding causes decohence due to the grid points not being
        ##infite. Must multiply by this factor
        ##This is a density factor correction based on Barry et al 2019.
        self.decoherence_factor = 2.

        print("Params either set or calculated:")
        print(f"{chr(11)}Central wavelength (m) {cent_wavelength:.3f}")
        print(f"{chr(11)}DM {self.DM:.1f}", )
        print(f"{chr(11)}Hubble_distance {self.hubble_distance}")
        print(f"{chr(11)}Bandwidth {bandwidth:.5e}")
        print(f"{chr(11)}redshift {z:.3f}")
        print(f"{chr(11)}Num freq chan {parser_args.N_chan}")
        print(f"{chr(11)}Ez", self.Ez)

        print(f"OVERALL NORMALISATION APPLIED (includes decoherence etc):  {self.decoherence_factor*self.normalisation:.8e}")

        if parser_args.wedge_factor > 0:
            self.wedge_factor = parser_args.wedge_factor
        else:
            ##Used to do the wedge cut for 1D
            self.wedge_factor = self.DM * self.Ez / (self.hubble_distance * (z + 1))

    def _open_and_reshape(self, filename):
        """Read in a CHIPS binary output and reshape into a 2D array"""

        with open(filename, 'rb') as train_xf:
            # print(filename)
            ##Reads in as 1D
            data = np.fromfile(train_xf, dtype=np.float32)
            ##Make 2D - data were written out by spatial direction, then spectral,
            ##reshape into 2D using the parser arguments
            data = np.reshape(data, (self.parser_args.N_kperp,self.parser_args.N_chan))

            ##TODO through some useful error if the reshaping cannot be done

            ##python lists the 'y' axes first, 'x' second so swap the axes
            ##to play nicely with imshow and the like
            data = np.swapaxes((data),1,0)

            ##Limit data to only the positive spectral frequencies
            data = data[self.parser_args.Neta-1:,:]

            # print(f"pca Maximum {data.max():.3e}")

            ##There are some bins that always get ignored so throw them away now
            data = data[:-KPARRA_END, KPERP_START:]

        return data

    def _read_in_data_and_convert(self, polarisation, chips_tag=False, oneD=False):
        """Attempt to read in the data based on user provided paths and polarisation
        options. Will first check for files that have had kriging (have a zero
        in the title), if it can't find those, look for ones without kriging
        (have a one in the title).

        Convert to a 2D array for 2D plot by default, or if specified, a 1D array.
        1D array requires one extra CHIPS output file to have been downloaded"""

        ##For ratio plots, need to pick denominator or numerator
        if chips_tag:
            pass
        ##if not, should only need to use the chips_tag in self.parser_args.chips_tag
        else:
            chips_tag = self.parser_args.chips_tag

        filename = False
        kriging = -1

        for kriging_attempt in [0, 1, 20, 21, 22]:

            filename_attempt = f"{self.parser_args.basedir}crosspower_{polarisation}_{kriging_attempt}.iter.{chips_tag}.dat"

            if os.path.isfile(filename_attempt):
                filename = filename_attempt
                kriging = kriging_attempt

        if not filename:
            msg = f'Could not open crosspower files based in input params.{chr(10)}' \
                f'Searched for filenames based on:{chr(10)}' \
                f"{self.parser_args.basedir}crosspower_{polarisation}_*.iter.{chips_tag}.dat"
            sys.exit(msg)

        crosspower = self._open_and_reshape(filename)

        filename = f"{self.parser_args.basedir}outputweights_{polarisation}_{kriging}.iter.{chips_tag}.dat"

        if not os.path.isfile(filename):
            msg = f'Could not open the outputweights file with the same kriging{chr(10)}' \
            f'number as the crosspower. Searched for the following file:{chr(10)}' \
            f'{filename}'
            sys.exit(msg)

        weights = self._open_and_reshape(filename)

        crosspower = masked_array(crosspower, weights == 0)

        # print("Just read in", crosspower)
        crosspower = crosspower / weights

        ##32 is a CHIPS based number hard coded to make the weights sensible
        ##For anything else, stick to one??
        ##TODO used to be 36, pls why?
        ##TODO sort this max weight nightmare
        # max_weights = 2
        # weight_adjustment = 32 / 2

        if self.parser_args.non_RTS_outputs:
            weight_scheme = 1
        else:
            weight_scheme = 32

        self.weight_data = weights/(self.normalisation)**2*weight_scheme*np.sqrt(self.parser_args.N_chan) / self.decoherence_factor
        self.crosspower = crosspower*self.decoherence_factor*self.normalisation

        print(f"Max power in file {self.crosspower.max():.4e}")

        # ##There are parts of the 2D outputs that always get thrown away,
        # ##so make the same cut on the output coords
        # self.kpa = self.kpa[:-KPARRA_END]
        # self.kper = self.kper[KPERP_START:]

    def read_data_and_create_2Darray(self, polarisation, chips_tag=False):
        """Attempt to read in the data based on user provided paths and polarisation
        options. Will first check for files that have had kriging (have a zero
        in the title), if it can't find those, look for ones without kriging
        (have a one in the title).

        Convert data to a 2D array for a 2D plot"""

        self._read_in_data_and_convert(polarisation, chips_tag=chips_tag, oneD=False)

        ##make an 'extent' list for imshow, that details x,y coords for a 2D plot
        ##goes as [low x coord, high x coord, low y coord, high y coord]
        extent = [self.kper[0], self.kper[-1], self.kpa[0], self.kpa[-1]]

        twoD_ps_array = self.crosspower

        # ##meshgrid them to find the length in kspace of all bins
        # k_perp_mesh, k_parr_mesh = np.meshgrid(self.kper, self.kpa)
        # k_lengths_mesh = np.sqrt(k_perp_mesh**2 + k_parr_mesh**2)
        #
        # u_arr_mesh, eta_mesh = np.meshgrid(self.u_arr[KPERP_START:], self.eta[:-KPARRA_END])
        # np.savez("2D_coords_and_power.npz", u_arr_mesh=u_arr_mesh,
        #                                     eta_mesh=eta_mesh,
        #                                     k_perp_mesh=k_perp_mesh,
        #                                     k_parr_mesh=k_parr_mesh,
        #                                     twoD_power=self.crosspower)

        return twoD_ps_array, extent

    def _grid_2D_to_1D_k3_inside(self, k_perp, k_parra, twoD_data, twoD_weights, ktot_bin_edges, convert_to_delta=True):
        twoD_weights_sqrt = np.sqrt(twoD_weights)

        ##meshgrid them to find the length in kspace of all bins
        k_perp_mesh, k_parr_mesh = np.meshgrid(k_perp, k_parra)
        k_lengths_mesh = np.sqrt(k_perp_mesh**2 + k_parr_mesh**2)

        ##This does the wedge cut?
        wedge_cut = k_parr_mesh > k_perp_mesh*self.wedge_factor

        ##Cuts off in k perpendicular, avoids small spatial scales in the 1D
        k_perp_cut = k_perp_mesh <= self.parser_args.kperp_max

        ##Cuts off in k perpendicular, avoids large spatial scales in the 1D
        k_perp_cut_min = k_perp_mesh > self.parser_args.kperp_min

        ##Cuts off in k parallel, avoiding things close to the wedge
        kparra_cut = k_parr_mesh > self.parser_args.kparra_min

        nozero_per = k_perp_mesh > 0.0
        nozero_par = k_parr_mesh > 0.0

        ##Find the centre of all the bins as the gridding coords
        ktot_bins = (ktot_bin_edges[1:] + ktot_bin_edges[:-1])/2

        ##How many bins we have, and make a zero array for gridding 1D power
        num_ktot_bins = len(ktot_bins)
        oneD_power = np.zeros(int(num_ktot_bins))
        oneD_delta = np.zeros(int(num_ktot_bins))
        oneD_weights = np.zeros(int(num_ktot_bins))
        ##Keep track of the locations of the bins on the 2D PS if we want to
        ##plot the wedge cut
        binning_array = np.ones(twoD_data.shape)*-1.0

        for k_tot_ind in range(num_ktot_bins):
            ##This finds all bins that sit inside the current annulus
            above_min = k_lengths_mesh > ktot_bin_edges[k_tot_ind]
            below_max = k_lengths_mesh <= ktot_bin_edges[k_tot_ind + 1]

            # cut_inds = np.where(above_min & below_max & wedge_cut & k_perp_cut & nozero_par & nozero_per)
            cut_inds = np.where(above_min & below_max & wedge_cut & k_perp_cut & nozero_par & nozero_per & kparra_cut & k_perp_cut_min)

            ##Always get annoying warnnging that a Masked element has been set
            ##to NaN here, so ignore them
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning,
                    message="Warning: converting a masked element to nan.")

                oneD_delta[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights_sqrt[cut_inds]**2*k_lengths_mesh[cut_inds]**3)
                oneD_weights[k_tot_ind] = np.nansum(twoD_weights_sqrt[cut_inds]**2)

                oneD_power[k_tot_ind] = np.nansum(twoD_data[cut_inds]*twoD_weights_sqrt[cut_inds]**2)

            binning_array[cut_inds] = k_tot_ind + 1

        self.binning_array = binning_array

        oneD_power = oneD_power / oneD_weights
        oneD_delta = oneD_delta / oneD_weights

        ##Which sigma to report the noise to
        sigma = 2
        sqrt_weights = np.ones(len(oneD_weights))
        sqrt_weights[np.where(oneD_weights != 0)] = np.sqrt(oneD_weights[np.where(oneD_weights != 0)])
        oneD_noise = sigma / sqrt_weights

        ##Convert to delta
        oneD_delta = oneD_delta / (2*np.pi**2)
        oneD_noise = oneD_noise*ktot_bins**3 / (2*np.pi**2)

        return ktot_bins, oneD_noise, oneD_power, oneD_delta

    def _plot_wedge_cut(self):
        """Plot the 2D power, what was cut, retained, and which bins were used
        during the gridding to a 1D PS."""

        ##These indexes will represent where power was cut
        wedge_cut_2D = np.where(self.binning_array == -1.0)

        ##Make copies of the data to cut
        power_cut = deepcopy(self.crosspower)
        power_retained = deepcopy(self.crosspower)

        power_retained[wedge_cut_2D] = 0.0

        ##Just subtract off what was cut to see what was retained
        power_cut -= power_retained


        fig, axs = plt.subplots(2,2, figsize=(8,10))

        extent = [self.kper[0], self.kper[-1], self.kpa[0], self.kpa[-1]]

        ##Mask negative power, too much of a hassle otherwise
        twoD_ps_array = masked_array(self.crosspower, self.crosspower <= 0)
        power_cut = masked_array(power_cut, power_cut <= 0)
        power_retained = masked_array(power_retained, power_retained <= 0)

        ##log10 hates all the zeros and negatives etc so mute a tonne
        ##of warnings

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)

            im0 = axs[0,0].imshow(np.log10(twoD_ps_array), origin='lower',
                                  interpolation='none', extent=extent,
                                  aspect='auto')

            im1 = axs[0,1].imshow(np.log10(power_cut), origin='lower',
                                  interpolation='none', extent=extent,
                                  aspect='auto')

            im2 = axs[1,0].imshow(np.log10(power_retained), origin='lower',
                                  interpolation='none', extent=extent,
                                  aspect='auto')

        ##Mask out bins that aren't used in gridding
        binning_array = masked_array(self.binning_array, self.binning_array == -1)

        im3 = axs[1,1].imshow(binning_array, origin='lower', interpolation='none',
                      extent=extent, aspect='auto')

        axs[0,0].set_title('All power')
        axs[0,1].set_title('Power cut')
        axs[1,0].set_title('Power retained')
        axs[1,1].set_title('Bins applied')

        for ax,im in zip(axs.flatten(), [im0, im1, im2, im3]):

            ax.set_xscale("log")
            ax.set_yscale("log")

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)

            label = "log10(P(k))"
            if im == im3: label = 'Bin number'
            cb = fig.colorbar(im, cax=cax,extend='both', label=label) # format='%.0e'

        # print(np.nanmax(binning_array))

        axs[1,0].set_xlabel(f"{chr(36)}k_{chr(92)}bot{chr(36)}")
        axs[1,1].set_xlabel(f"{chr(36)}k_{chr(92)}bot{chr(36)}")
        axs[0,0].set_ylabel(f"{chr(36)}k_{chr(92)}parallel{chr(36)}")
        axs[0,1].set_ylabel(f"{chr(36)}k_{chr(92)}parallel{chr(36)}")

        plt.tight_layout()
        print("Saving wedgecut_2D.png")
        fig.savefig('wedgecut_2D.png',bbox_inches='tight')
        plt.close()

    def read_data_and_create_1Darray(self, polarisation, chips_tag=False):

        self._read_in_data_and_convert(polarisation, chips_tag=chips_tag, oneD=True)

        if self.parser_args.ktot_bin_edges:

            ktot_bin_edges = np.loadtxt(self.parser_args.ktot_bin_edges)

            # if os.path.isfile("./"+self.parser_args.ktot_bin_edges):
            #     ktot_bin_edges = np.loadtxt(self.parser_args.ktot_bin_edges)
            # else:
            #     exit(f"Cannot find --ktot_bin_edges={self.parser_args.ktot_bin_edges}, file doesn't exist")
        else:
            low_k_edge = self.parser_args.low_k_edge
            high_k_edge = self.parser_args.high_k_edge
            num_k_edges = self.parser_args.num_k_edges
            ktot_bin_edges = 10**np.linspace(np.log10(low_k_edge), np.log10(high_k_edge), num_k_edges)

        self.ktot_bin_edges = ktot_bin_edges

        # print("BOTTOM EDGE, UPPER EDGE", ktot_bin_edges[0],  ktot_bin_edges[-1])

        ktot_bins, oneD_noise, oneD_power, oneD_delta = self._grid_2D_to_1D_k3_inside(self.kper,
             self.kpa, self.crosspower, self.weight_data, ktot_bin_edges)

        ##If requested, make a 2D plot of the cuts applied
        if self.parser_args.plot_wedge_cut_2D:

            self._plot_wedge_cut()

        return ktot_bins, oneD_noise, oneD_power, oneD_delta

def create_positives_cmap(args, set_under_colour=False):
    """Sets up a colourmap (cmap) with a log10 normalisation based on the input
    arguments. Optionally, sets anything under zero to a specific colour"""
    ##Setup up a bespoke colour map
    cmap = mpl.cm.get_cmap("Spectral_r")

    upper = np.ceil(np.log10(args.max_power))
    lower = np.ceil(np.log10(args.min_power))

    if set_under_colour:
        cmap.set_under(set_under_colour)
        cmap.set_bad(set_under_colour)
        ##Set zero as the lower bound, so that everything under is the desired
        ##colour
        bounds = [0]
    else:
        bounds = []
        cmap.set_bad(cmap(0))

    for bound in np.linspace(lower, upper, cmap.N):
        # bounds.append(10**bound)
        bounds.append(bound)

    # ticks = [10**tick for tick in np.arange(lower, upper + 1, 2)]
    ticks = [tick for tick in np.arange(lower, upper + 1, 2)]
    labels = ["1e+{:d}".format(int(tick)) for tick in ticks]

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    return cmap, norm, ticks, labels

def do_2D_axes_labels(ax, title, polarisation,
                      hide_cbar_label=False, hide_k_perp_label=False):
    """Makes the 2D axis log scale, sticks labels on the axes, sets a title"""

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    if polarisation == 'xx':
        title_pol = 'XX'
    else:
        title_pol = 'YY'

    ax.set_title(f'{title_pol} {title}',size=16)
    ax.set_xlabel(f'k{chr(36)}_{chr(92)}bot{chr(36)} ({chr(36)}h{chr(36)}Mpc{chr(36)}^{-1}{chr(36)})',fontsize=18)
    if not hide_k_perp_label:
        ax.set_ylabel(f'k{chr(36)}_{chr(91)}parallel{chr(36)} ({chr(36)}h{chr(36)}Mpc{chr(36)}^{-1}{chr(36)})',fontsize=18)
    ax.set_xscale("log")
    ax.set_yscale("log")

def plot_2D_on_ax(twoD_ps_array, extent, ax, fig, polarisation,
                  args, hide_cbar_label=False, hide_k_perp_label=False):
    """Plot the 2D PS data in `twoD_ps_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`. Uses a single
    cmap"""

    if args.colourscale == "negs_are_grey":
        cmap, norm, ticks, ticklabels = create_positives_cmap(args, "Grey")

    else:
        cmap, norm, ticks, ticklabels = create_positives_cmap(args)

    ##Gotta mask out the negatives or matplotlib just wigs out and shows
    ##incorrect colours everywhere
    twoD_ps_array = masked_array(twoD_ps_array, twoD_ps_array <= 0)

    ##Stop warnings when plotting log10 of zero
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")


        im = ax.imshow(np.log10(twoD_ps_array), cmap=cmap, origin='lower',
                       norm=norm, aspect='auto', extent=extent,
                       interpolation='none')

    ##Append a smaller axis to plot the colourbar on
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)

    cb = fig.colorbar(im, cax=cax, format='%.0e',extend='both',
                      ticks=ticks)
    cb.ax.set_yticklabels(ticklabels)

    if not hide_cbar_label:
        cax.set_ylabel(f'P(k) mK{chr(36)}^2{chr(36)} {chr(36)}h^{-3}{chr(36)} Mpc{chr(36)}^3{chr(36)}',fontsize=14)
    cb.ax.tick_params(labelsize=11)

    if args.colourscale == "negs_are_grey":
        cax.text(1.13, 0.01, f'{chr(36)}{chr(92)}leq0{chr(36)}', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=11)

    do_2D_axes_labels(ax, args.title, polarisation, hide_cbar_label, hide_k_perp_label)

def plot_2D_on_ax_two_colour_bars(twoD_ps_array, extent, ax, fig, cax_pos,
                  cax_neg, polarisation, args, title, cmap='PurpOrang',
                  hide_cbar_label=False, hide_k_perp_label=False):
    """Plot the 2D PS data in `twoD_ps_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`. Uses two different
    cmaps, one for positive values, another for negative. Positive colourbar
    is plotted on axes `cax_pos`, negative on `cax_neg`."""

    pos_data = masked_array(twoD_ps_array, twoD_ps_array <= 0)
    neg_data = masked_array(twoD_ps_array, twoD_ps_array >= 0)

    ##==========================================================================
    ##Do the positive plot
    ##==========================================================================
    if cmap == 'PurpOrang':
        # cmap_pos = mpl.cm.get_cmap("Oranges").copy()
        cmap_pos = "Oranges"
    elif cmap == 'BlueRed':
        # cmap_pos = mpl.cm.get_cmap("Blues").copy()
        cmap_pos = "Blues"

    vmax = np.log10(args.max_power)
    vmin = np.log10(args.min_power)

    tick_high = np.floor(vmax)
    tick_low = np.ceil(vmin)

    inc = int(np.ceil((tick_high - tick_low) / 6))
    pos_ticks =  np.arange(tick_low, tick_high + 1, inc)

    # pos_bounds = np.linspace(lower, upper, 100)
    # pos_ticks =  np.arange(lower, upper + 1, 2)
    pos_labels = [f"1e+{int(tick)}" for tick in pos_ticks]

    ##Stop warnings when plotting log10 of zero
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")

        im_pos = ax.imshow(np.log10(pos_data), cmap=cmap_pos, origin='lower',
                           vmin=vmin, vmax=vmax,
                           aspect='auto', extent=extent, interpolation='none')

    cb_pos = fig.colorbar(im_pos, cax=cax_pos, extend='both', ticks=pos_ticks)
    cb_pos.ax.set_yticklabels(pos_labels)
    cb_pos.ax.tick_params(labelsize=11)

    ##==========================================================================
    ##Do the negative plot
    ##==========================================================================

    if args.max_neg_power:
        neg_upper = abs(args.min_neg_power)
    else:
        neg_upper = args.max_power

    if args.min_neg_power:
        neg_lower = abs(args.max_neg_power)
    else:
        neg_lower = args.min_power

    if cmap == 'PurpOrang':
        # cmap_neg = mpl.cm.get_cmap("Purples_r").copy()
        cmap_neg = "Purples_r"
    elif cmap == 'BlueRed':
        # cmap_neg = mpl.cm.get_cmap("Reds_r").copy()
        cmap_neg = "Reds_r"

    vmin = -np.log10(neg_upper)
    vmax = -np.log10(neg_lower)

    tick_low = np.ceil(np.log10(neg_lower))
    tick_high = np.floor(np.log10(neg_upper))

    inc = int(np.ceil((tick_high - tick_low) / 6))
    neg_ticks =  -np.arange(tick_low, tick_high + 1, inc)

    neg_labels = [f"-1e+{int(abs(tick))}" for tick in neg_ticks]

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="invalid value encountered in log10")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
            message="divide by zero encountered in log10")

        im_neg = ax.imshow(-np.log10(np.abs(neg_data)), cmap=cmap_neg, origin='lower',
                           vmin=vmin, vmax=vmax,
                           aspect='auto', extent=extent, interpolation='none')

    cb_neg = fig.colorbar(im_neg, cax=cax_neg, extend='both',
                          ticks=neg_ticks)
    cb_neg.ax.set_yticklabels(neg_labels)
    cb_neg.ax.tick_params(labelsize=11)

    if not hide_cbar_label:
        cax_pos.yaxis.set_label_coords(4.5,0.0)
        cax_pos.set_ylabel(f'P(k) mK{chr(36)}^2{chr(36)} {chr(36)}h^{-3}{chr(36)} Mpc{chr(36)}^3{chr(36)}',fontsize=14)

    do_2D_axes_labels(ax, title, polarisation, hide_cbar_label, hide_k_perp_label)

def setup_ax_and_cax_for_double_colourbar(args, fig, num_axes, polarisation):
    """When using two cmaps, sets up an axes and two colour bar axes (cax)
    based on how many axes you are going to plot"""

    if num_axes == 1:
        ax = fig.add_axes([0.16, 0.1, 0.6, 0.84])
        cax_height = 0.42
        cax_width = 0.035
        cax_left = 0.78
        cax_pos = fig.add_axes([cax_left, 0.52, cax_width, cax_height])
        cax_neg = fig.add_axes([cax_left, 0.1, cax_width, cax_height])

    elif num_axes == 2:
        ax_width = 0.32
        ax_height = 0.84
        cax_height = 0.42
        cax_width = 0.02

        if polarisation == 'xx':
            cax_xx_left = 0.42
            ax = fig.add_axes([0.09, 0.1, ax_width, ax_height])
            cax_pos = fig.add_axes([cax_xx_left, 0.52, cax_width, cax_height])
            cax_neg = fig.add_axes([cax_xx_left, 0.1, cax_width, cax_height])

        elif polarisation == 'yy':
            cax_yy_left = 0.88
            ax = fig.add_axes([0.55, 0.1, ax_width, ax_height])
            cax_pos = fig.add_axes([cax_yy_left, 0.52, cax_width, cax_height])
            cax_neg = fig.add_axes([cax_yy_left, 0.1, cax_width, cax_height])

    return ax, cax_pos, cax_neg


def do_2D_plot(chips_data):
    """Given the user supplied `args`, plot a 2D power spectrum"""

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 1e+12
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = 1e+3

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        # ##Read in data and convert to a 2D array for plotting
        # crosspower, weights = read_in_data(args, chips_data.parser_args.polarisation)
        # twoD_ps_array, extent = convert_to_2D_PS_array(crosspower, weights)

        twoD_ps_array, extent = chips_data.read_data_and_create_2Darray(chips_data.parser_args.polarisation)

        if chips_data.parser_args.colourscale == 'pos_and_negs':
            fig = plt.figure(figsize=(6,7))
            num_axes = 1
            ax, cax_pos, cax_neg = setup_ax_and_cax_for_double_colourbar(args,
                                                fig, num_axes, chips_data.parser_args.polarisation)

            plot_2D_on_ax_two_colour_bars(twoD_ps_array, extent, ax, fig,
                                cax_pos, cax_neg, chips_data.parser_args.polarisation, args,
                                args.title)
        else:

            fig, ax = plt.subplots(1,1,figsize=(6,7))
            plot_2D_on_ax(twoD_ps_array, extent, ax, fig, chips_data.parser_args.polarisation, args)
            plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_{chips_data.parser_args.chips_tag}_crosspower.png"

    elif chips_data.parser_args.polarisation == 'both':

        twoD_ps_array_xx, extent_xx = chips_data.read_data_and_create_2Darray('xx')
        twoD_ps_array_yy, extent_yy = chips_data.read_data_and_create_2Darray('yy')

        if chips_data.parser_args.colourscale == 'pos_and_negs':

            fig = plt.figure(figsize=(11,7))

            num_axes = 2

            ax_xx, cax_pos_xx, cax_neg_xx = setup_ax_and_cax_for_double_colourbar(args,
                                                                   fig, num_axes, 'xx')
            plot_2D_on_ax_two_colour_bars(twoD_ps_array_xx, extent_xx,
                                ax_xx, fig, cax_pos_xx, cax_neg_xx, 'xx', args,
                                args.title, hide_cbar_label=True)

            ax_yy, cax_pos_yy, cax_neg_yy = setup_ax_and_cax_for_double_colourbar(args,
                                                                   fig, num_axes, 'yy')
            plot_2D_on_ax_two_colour_bars(twoD_ps_array_yy, extent_yy,
                                ax_yy, fig, cax_pos_yy, cax_neg_yy, 'yy', args,
                                args.title, hide_k_perp_label=True)

        else:
            fig, axs = plt.subplots(1,2,figsize=(12,7))
            plot_2D_on_ax(twoD_ps_array_xx, extent_xx, axs[0], fig, 'xx', args,
                          hide_cbar_label=True)
            plot_2D_on_ax(twoD_ps_array_yy, extent_yy, axs[1], fig, 'yy', args,
                          hide_k_perp_label=True)

            plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_{chips_data.parser_args.chips_tag}_crosspower.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def make_2D_ratio(chips_data, polarisation):
    """Using input user arguments `args`, make a 2D Ratio array for the given
    `polarisation`"""
    twoD_ps_array_numer, extent_numer = chips_data.read_data_and_create_2Darray(polarisation,
                                                      chips_data.parser_args.chips_tag_one)

    twoD_ps_array_denom, extent_denom = chips_data.read_data_and_create_2Darray(polarisation,
                                                      chips_data.parser_args.chips_tag_two)
    ratio = twoD_ps_array_numer / twoD_ps_array_denom

    return ratio, extent_numer


def plot_2D_ratio_on_ax(twoD_ps_ratio_array, extent, ax, fig, polarisation,
                  args, hide_cbar_label=False, hide_k_perp_label=False):
    """Plot the 2D ratio PS data in `twoD_ps_ratio_array`, which covers the kper/k_par
    coords in `extent`, on the axes `ax` on figure `fig`"""

    norm = LogNorm(args.min_power,args.max_power)
#    cmap = mpl.cm.get_cmap('RdBu').copy()
    cmap = mpl.cm.get_cmap('RdBu')

    # if set_under_colour:
    # cmap.set_under('Grey')

    cmap.set_bad('Grey')

    im = ax.imshow(twoD_ps_ratio_array, cmap=cmap, origin='lower',
                   aspect='auto', extent=extent,
                   norm=norm, interpolation='none')
                   # vmin=args.min_power,vmax=args.max_power,

    ##Append a smaller axis to plot the colourbar on
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)

    # cb = fig.colorbar(im, cax=cax, format='%.0e',extend='both')
    cb = fig.colorbar(im, cax=cax, extend='both')

    if not hide_cbar_label:
        cax.set_ylabel('Ratio P(k) / P(k)',fontsize=14)
    cb.ax.tick_params(labelsize=11)

    if args.chips_tag_one_label and args.chips_tag_two_label:
        title = f'Ratio {chr(10)}{chips_data.parser_args.chips_tag_one_label} / {chips_data.parser_args.chips_tag_two_label}'
    else:
        title = 'Ratio'

    do_2D_axes_labels(ax, title, polarisation, hide_cbar_label, hide_k_perp_label)

def do_2D_ratio_plot(chips_data):
    """Given the user supplied `chips_data.parser_args`, plot a 2D power spectrum ratio"""

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 10
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = -10

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array, extent =  make_2D_ratio(chips_data, chips_data.parser_args.polarisation)


        fig, ax = plt.subplots(1,1,figsize=(6,7))
        plot_2D_ratio_on_ax(twoD_ps_ratio_array,
                            extent, ax, fig, chips_data.parser_args.polarisation, chips_data.parser_args)
        plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_ratio.png"

    elif chips_data.parser_args.polarisation == 'both':
        fig, axs = plt.subplots(1,2,figsize=(12,7))

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array_xx, extent =  make_2D_ratio(chips_data, 'xx')

        plot_2D_ratio_on_ax(twoD_ps_ratio_array_xx,
                            extent, axs[0], fig, 'xx', chips_data.parser_args)

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_ratio_array_yy, extent =  make_2D_ratio(chips_data, 'yy')

        plot_2D_ratio_on_ax(twoD_ps_ratio_array_yy,
                            extent, axs[1], fig, 'yy', chips_data.parser_args)

        plt.tight_layout()
        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_ratio.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def make_2D_diff(chips_data, polarisation):
    """Using input user arguments `args`, make a 2D Ratio array for the given
    `polarisation`"""

    twoD_ps_array_one, extent_one = chips_data.read_data_and_create_2Darray(polarisation,
                                                      args.chips_tag_one)

    twoD_ps_array_two, extent_two = chips_data.read_data_and_create_2Darray(polarisation,
                                                      args.chips_tag_two)

    difference = twoD_ps_array_one - twoD_ps_array_two

    return difference, extent_one

def do_2D_diff_plot(chips_data):
    """Given the user supplied `args`, plot a 2D power spectrum ratio"""

    if chips_data.parser_args.max_power == 0.0:
        chips_data.parser_args.max_power = 1e12
    if chips_data.parser_args.min_power == 0.0:
        chips_data.parser_args.min_power = 1e3

    if args.chips_tag_one_label and args.chips_tag_two_label:
        title = f'Difference{chr(10)}{chips_data.parser_args.chips_tag_one_label} - {chips_data.parser_args.chips_tag_two_label}'
    else:
        title = 'Difference'


    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array, extent =  make_2D_diff(chips_data, chips_data.parser_args.polarisation)

        fig = plt.figure(figsize=(6,7))
        num_axes = 1
        ax, cax_pos, cax_neg = setup_ax_and_cax_for_double_colourbar(chips_data.parser_args,
                                            fig, num_axes, chips_data.parser_args.polarisation)

        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array, extent, ax, fig,
                            cax_pos, cax_neg, chips_data.parser_args.polarisation, args,
                            title, cmap='BlueRed')

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_{chips_data.parser_args.polarisation}_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_diff.png"

    elif chips_data.parser_args.polarisation == 'both':

        fig = plt.figure(figsize=(11,7))

        num_axes = 2

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array_xx, extent_xx =  make_2D_diff(chips_data, 'xx')

        ##Read in data and convert to a 2D array for plotting
        twoD_ps_diff_array_yy, extent_yy =  make_2D_diff(chips_data, 'yy')

        ax_xx, cax_pos_xx, cax_neg_xx = setup_ax_and_cax_for_double_colourbar(args,
                                                               fig, num_axes, 'xx')
        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array_xx, extent_xx,
                            ax_xx, fig, cax_pos_xx, cax_neg_xx, 'xx', args,
                            title,cmap='BlueRed', hide_cbar_label=True)

        ax_yy, cax_pos_yy, cax_neg_yy = setup_ax_and_cax_for_double_colourbar(args,
                                                               fig, num_axes, 'yy')
        plot_2D_on_ax_two_colour_bars(twoD_ps_diff_array_yy, extent_yy,
                            ax_yy, fig, cax_pos_yy, cax_neg_yy, 'yy', args,
                            title,cmap='BlueRed', hide_k_perp_label=True)

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips2D_xx+yy_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_diff.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)

def plot_1D_on_ax(ax, oneD_k_modes, oneD_power_measured, oneD_delta_measured,
                  oneD_delta_2sig_noise,
                  min_power, max_power,
                  pol, colour_power='C0', marker_power='x',
                  label='', delta=False):
    """Plot the power and two sigma noise on the given axes"""

    notzero = np.where(oneD_k_modes != 0)
    plot_k_modes = oneD_k_modes[notzero]
    plot_delta = oneD_delta_measured[notzero]
    plot_power = oneD_power_measured[notzero]
    plot_noise = oneD_delta_2sig_noise[notzero]

    if pol == 'xx':
        pol_label = f'XX {label}'
    else:
        pol_label = f'YY {label}'

    if delta:
        ax.plot(plot_k_modes, plot_delta, drawstyle='steps-mid',
               color=colour_power, marker=marker_power, mfc='none', ms=4,
               label=f'{pol_label}Power')
        ax.set_ylabel(f'{chr(36)}{chr(92)}Delta{chr(36)} (mK{chr(36)}^2{chr(36)})', fontsize=16)

    else:
        ax.plot(plot_k_modes, plot_power, drawstyle='steps-mid',
               color=colour_power, marker=marker_power, mfc='none', ms=4,
               label=f'{pol_label}Power')
        ax.set_ylabel(f'P(k) mK{chr(36)}^2{chr(36)} {chr(36)}h^{-3}{chr(36)} Mpc{chr(36)}^3{chr(36)}',fontsize=16)

    np.savez_compressed(f"1D_power.npz", k_modes=plot_k_modes, power=plot_power,
                                        delta=plot_delta)



    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel(f'k ({chr(36)}h{chr(36)}Mpc{chr(36)}^{-1}{chr(36)})',fontsize=16)

    ax.set_xscale("log")
    ax.set_yscale("log")

    if max_power:
        ax.set_ylim(top=max_power)
    if min_power:
        ax.set_ylim(bottom=min_power)

    ax.legend(prop={"size":14}, loc='best')

def do_1D_plot(chips_data):
    """Plot a 1D power spectrum given the input arguments from the parser"""

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        ##Read in data and convert to a 1D array for plotting
        oneD_k_modes, oneD_delta_2sig_noise, oneD_power_measured, oneD_delta_measured =  chips_data.read_data_and_create_1Darray(chips_data.parser_args.polarisation)


        fig, ax = plt.subplots(1,1,figsize=(8,6))

        plot_1D_on_ax(ax, oneD_k_modes, oneD_power_measured,
                      oneD_delta_measured, oneD_delta_2sig_noise,
                      chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                      chips_data.parser_args.polarisation,
                      delta=chips_data.parser_args.plot_delta)

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_{chips_data.parser_args.polarisation}_{chips_data.parser_args.chips_tag}.png"

    elif chips_data.parser_args.polarisation == 'both':
        ##Read in data and convert to a 1D array for plotting
        oneD_k_modes_xx, oneD_delta_2sig_noise_xx, oneD_power_measured_xx, oneD_delta_measured_xx  =  chips_data.read_data_and_create_1Darray("xx")

        oneD_k_modes_yy, oneD_delta_2sig_noise_yy, oneD_power_measured_yy, oneD_delta_measured_yy =  chips_data.read_data_and_create_1Darray("yy")


        fig, axs = plt.subplots(2,1,figsize=(10,10))

        plot_1D_on_ax(axs[0], oneD_k_modes_xx, oneD_power_measured_xx,
                              oneD_delta_measured_xx, oneD_delta_2sig_noise_xx,
                              chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                              "xx", delta=chips_data.parser_args.plot_delta)

        plot_1D_on_ax(axs[1], oneD_k_modes_yy, oneD_power_measured_yy,
                              oneD_delta_measured_xx, oneD_delta_2sig_noise_yy,
                              chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                              "yy", delta=chips_data.parser_args.plot_delta)

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_xx+yy_{chips_data.parser_args.chips_tag}.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)


def save_or_plot(fig, output_plot_name, plot_mode):
    """Either save the figure to the given file name `output_plot_name`,
    or call plt.show(), based on whether plot_mode == 'png' or
    plot_mode == 'screen' or"""

    if plot_mode == 'png':
        print(f"Saving {output_plot_name}")
        plt.savefig(output_plot_name, bbox_inches='tight')
        plt.close()

    if plot_mode == 'screen':
        print("Printing to screen")
        plt.show()

def plot_1D_comparison(chips_data, pol, ax):
    """Given the axes `ax` and user supplied `chips_data.parser_args`, plot
    a comparison for the given polarisation `pol`"""

    ##Read in data and convert to a 1D array for plotting
    oneD_k_modes_1, oneD_delta_2sig_noise_1, oneD_power_measured_1, oneD_delta_measured_1 =  chips_data.read_data_and_create_1Darray(pol,
                                            chips_tag=chips_data.parser_args.chips_tag_one)


    if chips_data.parser_args.chips_tag_one_label:
        label_one = chips_data.parser_args.chips_tag_one_label + ' '
    else:
        label_one = chips_data.parser_args.chips_tag_one + ' '

    plot_1D_on_ax(ax, oneD_k_modes_1, oneD_power_measured_1,
                  oneD_delta_measured_1, oneD_delta_2sig_noise_1,
                  chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                  pol, label=label_one,
                  delta=chips_data.parser_args.plot_delta)


    if chips_data.parser_args.chips_tag_two_label:
        label_two = chips_data.parser_args.chips_tag_two_label + ' '
    else:
        label_two = chips_data.parser_args.chips_tag_two + ' '

    oneD_k_modes_2, oneD_delta_2sig_noise_2, oneD_power_measured_2, oneD_delta_measured_2 =  chips_data.read_data_and_create_1Darray(pol,
                                            chips_tag=chips_data.parser_args.chips_tag_two)

    plot_1D_on_ax(ax, oneD_k_modes_2, oneD_power_measured_2,
                  oneD_delta_measured_2, oneD_delta_2sig_noise_2,
                  chips_data.parser_args.min_power, chips_data.parser_args.max_power,
                  pol, label=label_two,
                  marker_power='o',colour_power='C1', delta=chips_data.parser_args.plot_delta)

def do_1D_comparison_plot(chips_data):
    """Given the user supplied `chips_data.parser_args`, plot 2 different
    1D power spectra on the same axes"""

    if chips_data.parser_args.polarisation == 'xx' or chips_data.parser_args.polarisation == 'yy':

        fig, ax = plt.subplots(1,1,figsize=(8,6))

        plot_1D_comparison(chips_data, chips_data.parser_args.polarisation, ax)

        plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_{chips_data.parser_args.polarisation}_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_comparison.png"

    elif chips_data.parser_args.polarisation == 'both':

        fig, axs = plt.subplots(2,1,figsize=(8,12))

        plot_1D_comparison(chips_data, 'xx', axs[0])
        plot_1D_comparison(chips_data, 'yy', axs[1])

        plt.tight_layout()

        output_plot_name = f"{chips_data.parser_args.outputdir}/chips1D_xx+yy_" \
                           f"{chips_data.parser_args.chips_tag_one}_" \
                           f"{chips_data.parser_args.chips_tag_two}_comparison.png"

    else:
        msg = f'--polarisation={chips_data.parser_args.polarisation} is not a valid argument{chr(10)}'
        'Must be one of either: xx, yy, both'
        sys.exit(msg)

    save_or_plot(fig, output_plot_name, chips_data.parser_args.plot_mode)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        args = get_args(sys.argv[1:])
    else:
        # is being called directly from nextflow
        args = get_args([
            "--basedir=${basedir}",
            "--polarisation=xx",
            "--chips_tag=${ext}",
            "--title=${title}",
            "--min_power=1e3",
            "--max_power=1e15"
        ])


    chips_data = ChipsDataProducts(args)

    if args.plot_type == '2D':
        do_2D_plot(chips_data)

    elif args.plot_type == '2D_Ratio' or args.plot_type == '2D_ratio':
        do_2D_ratio_plot(chips_data)

    elif args.plot_type == '2D_Diff' or args.plot_type == '2D_diff':
        do_2D_diff_plot(chips_data)

    elif args.plot_type == '1D':
        do_1D_plot(chips_data)

    elif args.plot_type == '1D_comp':
        do_1D_comparison_plot(chips_data)

    else:
        sys.exit(f"You entered --plot_type={args.plot_type}, which is not"
                 "a recognised type. You can enter either: 2D, 2D_ratio,"
                 " 2D_diffm or 1D")
