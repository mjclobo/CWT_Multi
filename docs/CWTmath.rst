=========================================
Basic CWT_Multi theory
=========================================

On this page we provide the basic theory required
to understand how CWT_Multi works.
This knowledge empowers the user to get the most out of
CWT_Multi for their specific application.
We present ideas and illustrations, rather than comprehensive maths.
Future pages in this documentation will cover the mathematic details.
For now, we refer the reader to Lobo et al. (2024) and references
therein.


Tidal species and constituents
~~~~~~~~~~~~~~~~~~~~~~~~~
Different *tidal constituents* describe unique frequencies that correspond
to different sets of Doodson numbers.
Doodson numbers describe the fundamental forcing temporal periods of the tides.
(The first Doodson number =s describes number of cycles per lunar day, the second Doodson
number describes number of cycles per sidereal year, etc.)
*Tidal species* describes groups of tidal constituents
that have the same number of cycles per lunar day, i.e, they
share a first Doodson number.


The Rayleigh criterion defined
~~~~~~~~~~~~~~~~~~~~~~~~~
The *Rayleigh criterion*, :math:`R_{C}`, describes the minimum analysis window
required to separate two signals with different frequencies.
That is, the minimum window length is

   .. math::
    L_{w} = R_{C} \, \left | f_{1} - f_{2} \right | ^{-1} \, ,

where we generally require :math:`R_{C} = 1`.
For example, in order to resolve the :math:`M_{2}`
and :math:`S_{2}` tidal constituents with :math:`R_{C}=1`,
we need to use an analysis window length of

   .. math::
     
     L_{w} 
     &= \left | f_{M_{2}} - f_{S_{2}} \right | ^{-1}  \\
     &= \left | 1/12.4206012 - 1/12 \right | ^{-1} \ \mathrm{hr}  \\
     &\approx 15 \ \mathrm{days} \, . 

Another way to word the Rayleigh criterion is:
if a signal is composed of sine waves at two different frequencies
we need the sine wave at the higher frequency to complete
at least one more cycle than the signal at the lower frequency,
in order to be able to tell the two frequencies apart.
This is shown graphically below.

.. image:: /images/RC_Lw.png
   :width: 600pt

The Rayleigh criterion in spectral space
~~~~~~~~~~~~~~~~~~~~~~~~~
A useful way to think of the Rayleigh criterion is in
terms of spectra.
A *spectrum* is a measure of a quantity as a function of
frequency, e.g., power spectrum, absorption spectrum.
Thus, we can draw a spectrum where the x-axis is
frequency, and the y-axis is some measure of tidal amplitude.

Below we show a the spectrum of a signal that is the sum of
equal-amplitude sine waves at the :math:`M_{2}` and :math:`S_{2}`
tidal frequencies.
We compute the energy spectrum for signals of three different
lengths: :math:`L_{w} / 2`, :math:`L_{w}`, and :math:`2 \, L_{w}`.
We find that with the shortest window we are not able to differentiate between
energy at the two frequencies (red line).
Once we analyze a signal that is at least the length :math:`L_{w}`,
we are able to resolve energy at the two frequencies (green line).

.. image:: /images/RC_spectra.png
   :width: 600pt

Note, however, that as the analyzed signal gets longer,
the peaks at the two frequencies become more distinct (yellow line).
If we had an infinitely long signal, the energy at the two frequencies would be represented by
vertical lines (hence the often-used term *line spectra*).
The apparent "spreading" of energy at frequencies around
:math:`M_{2}` and :math:`S_{2}` is an artifact of the finite-length
analysis window.

CWT_Multi application method for a full time series
~~~~~~~~~~~~~~~~~~~~~~~~~
The fundamental application of CWT_Multi is to *define
tidal amplitudes and phases that vary as functions of time*.
Here we provide a brief explanation of the framework used to accomplish this goal.

First, we note that CWT_Multi performs both a species and constituents analysis.
The *species analysis* defines time-varying amplitudes and phases for each tidal species,
i.e., diurnal (:math:`D_{1}`), semidiurnal (:math:`D_{2}`), etc.
This analysis can resolve time-changes in species amplitudes on the order of a couple/few days.

The *constituents analysis* defines time-varying amplitudes and phases for 7-9 individual tidal
constituents within the diurnal and semidiurnal tidal species bands.
Since constituents within the same species are fairly close together (below, we will detail how the
closeness of the :math:`M_{2}` and :math:`S_{2}` constituents affects our analysis, for example),
we resolve time-changes of constituent amplitudes on the order of one to two weeks.

The main steps that the CWT_Multi analysis is comprised of are:

1. Define the analysis window for a given time step, centered on time :math:`t_m`
2. Convolve each filter from the filter bank with data within the analysis window.
   (This step outputs a complex response.)
3. Solve the response coefficient matrix problem (detailed below).
4. Store complex solution for all frequencies that have corresponding filters at the time :math:`t_m`.
   (From this complex solution, one easily retrieves amplitude and phase.)
5. Move the analysis window forward to :math:`t_m \, + \, D_{f} \Delta t`, where :math:`D_{f}` is
   the decimation factor, i.e., the number of time steps between adjacent CWT_Multi analyses, and
   :math:`\Delta t` is the sampling period.
6. Repeat.


We now describe the maths behind the CWT_Multi process that occurs at each analysis time step,
centered on :math:`t_m`.


CWT_Multi filters
~~~~~~~~~~~~~~~~~~~~~~~~~
The spectra shown above were constructed using Fourier transforms.
The Fourier amplitude at a given frequency, :math:`f`, is essentially the magnitude of the convolution
of a complex sinusoid, of the form

   .. math::
    e^{i \, t \,2 \, \pi \, f}
    = \mathrm{cos}(2 \pi f t ) + i \, \mathrm{sin} (2 \pi f t )  \, ,

with the signal being analyzed, over the analysis window length.
The complex output then contains the information necessary to find
the amplitude and phase of the signal at the frequency :math:`f`.

CWT_Multi performs analogous convolutions using complex wavelet filters.
An example of such a filter is shown below.

.. image:: /images/M2_wavelet.png
   :width: 600pt

In short, the form of our wavelet maximizes the amount
of information one is able to extract from this convolution
given a finite analysis window length.
However, the optimal form of wavelets are a topic of active
research, and always require some trade-off (see Lilly and Ohelde 2012).


CWT_Multi defines wavelets at frequencies where tidal energy is
expected, and then constructs a matrix problem for the complex
convolution output.
This matrix problem allows for resolution of frequencies for
analysis windows of lengths that violate the Rayleigh criterion.
We will soon present the assumptions and methods of the response coefficient
matrix.
First, we must understand what a frequency response is, and how this
concept manifests in CWT_Multi.

Frequency response: A definition
~~~~~~~~~~~~~~~~~~~~~~~~~
From the spectrum plot above, we see that finite-length
complex sinusoids (and wavelet filters) within a given frequency
band, which we define as :math:`f \pm \Delta f`, will "respond" to
energy at the central frequency, :math:`f`.
Importantly, this *frequency response* is a function
of the analysis window length.
Shorter filters (equivalently, shorter analysis windows) will
increase the frequency range, :math:`\Delta f`, at which the filter
will respond to energy at adjacent frequencies.

**CWT_Multi leverages the frequency response of filters
centered on tidal frequencies to energy at adjacent tidal frequencies**
to construct a matrix problem.
We now present this matrix problem.


Response coefficient matrix: The problem
~~~~~~~~~~~~~~~~~~~~~~~~~
The response coefficient matrix problem is

   .. math::
    \vec{f} (t_m) = \boldsymbol{R} \, \vec{a}(t_m) \, ,

where:

- :math:`t_m` is the time at the center of the analysis window
- :math:`\vec{f}` is an :math:`N \times 1` column vector of the complex output from
  the :math:`N` complex wavelet filters (at frequency :math:`f_n`) with signal, centered on time :math:`t_m`
- :math:`\boldsymbol{R}` is the *response coefficient matrix* (RCM), which we describe in detail below
- :math:`\vec{a}(t_m)` is the :math:`N \times 1` column vector of the true amplitudes
  of the signal at the frequencies :math:`f_n`

The easiest way to understand the RCM is in terms of a simplified problem.
Consider a set of wavelet filters at the :math:`M_{2}` and :math:`S_{2}` frequencies,
where we would like to define the :math:`M_{2}` and :math:`S_{2}`
amplitudes as a function of time.
We thus define the RCM as

   .. math::
    \boldsymbol{R} =
    \begin{pmatrix}
    r_{M_{2}, \, M_{2}} & r_{M_{2}, \, S_{2}} \\
    r_{S_{2}, \, M_{2}} & r_{S_{2}, \, S_{2}}
    \end{pmatrix} \, ,

where :math:`r_{f_{1}, \, f_{2}}` describes the frequency of the :math:`f_{1}` filter
to energy at the :math:`f_{2}` frequency, with a maximum value of unity.
For example, :math:`r_{M_{2}, \, M_{2}} = 1`, since the :math:`M_{2}` filter will
respond to all of the energy at the :math:`M_{2}` frequency.

As noted above, the filter width in time (equivalently, the length of the analysis window),
will determine the width in frequency-space, :math:`\Delta f`, at which
the filter will respond to energy at adjacent frequencies.
We can now plot the frequency response for our simplified problem.
In particular, we show the filter responses for the two filters for two different
choices of wavelet filter length.

.. image:: /images/RCM_filter_response.png
   :width: 700pt

We show the frequency response for the :math:`M_{2}` (red)
and :math:`S_{2}` (blue) filters above, as a function of frequency.
For the narrower filters (panel (a)), the surrounding band of frequencies, for which the
respective filters respond to energy, is relatively wide.
In particular, :math:`r_{S_{2}, \, M_{2}} \approx 0.45` means that the :math:`S_{2}` filter
will include 45% of the energy that exists at the :math:`M_{2}` frequency in its estimate
of the amplitude of the :math:`S_{2}` component of the signal during the analysis window.
Though this may seem like a problem, we will explain how the RCM accounts for such overlap in the following section.
First, we review some salient aspects of the frequency response plot, and their connections to the RCM.

Here are some things to note for the frequency response figure above:

- We have :math:`r_{M_{2}, \, M_{2}} = 1` and :math:`r_{S_{2}, \, S_{2}} = 1`,
  as expected
- If the :math:`M_{2}` and :math:`S_{2}` filters are the same length, as above,
  then we have :math:`r_{S_{2}, \, M_{2}} = r_{M_{2}, \, S_{2}}`, and the RCM is a
  symmetric matrix
- The wider the filter in time, i.e., the longer the analysis window, the more narrow
  the frequency response is

The last point should be thought upon, as it is this feature of the RCM that guides
one's choice of filter lengths when using CWT_Multi.
**The user must choose a trade-off between having time-resolution (i.e., being able
to define a tidal amplitude that varies as a function of time) and frequency-resolution
(i.e., being able to distinguish energy between two frequencies.**

.. note::
    The reader might be wondering why the 15-day-long wavelet filters respond to nearby frequencies,
    whereas the Rayleigh criterion suggests that 15 days is long enough to resolve the :math:`M_{2}`
    and :math:`S_{2}` signals.
    This is because the wavelet filters are tapered, and carry about 80% of their energy in the middle
    half of the filter (see the plot of complex wavelet filter above).
    So the effective length of a wavelet filter, in terms of a Rayleigh criterion, is close to about half
    of the user-specified wavelet filter length.



Response coefficient matrix: The solution
~~~~~~~~~~~~~~~~~~~~~~~~~
We have defined the response coefficient matrix (RCM), and have hopefully
provided some insight into its meaning and its connection to CWT_Multi analysis.
As a final stop in our exposition of the theory that supports CWT_Multi analysis,
we consider the solution to the RCM problem.


The RCM problem (also defined above) is

   .. math::
    \vec{f} (t_m) = \boldsymbol{R} \, \vec{a}(t_m) \, ,

In the example currently under consideration, we consider filters
only at the :math:`M_{2}` and :math:`S_{2}` tidal frequencies.
Now, suppose that signal only has energy at the :math:`M_{2}` and :math:`S_{2}`
frequencies, each with unity amplitude.

For filters that are 15 days long (panel (a)) above, our RCM problem
becomes

    .. math::
     \begin{pmatrix}
     1.45 \\
     1.45 \\
     \end{pmatrix}
     =
     \begin{pmatrix}
     1.0 & 0.45 \\
     0.45 & 1.0 
     \end{pmatrix}
     \ \begin{pmatrix}
     a_{M_{2}} \\
     a_{S_{2}}
     \end{pmatrix} \, .

By multiplying both sides by :math:`\boldsymbol{R}^{-1}` we find

    .. math::
     \vec{a} =
     \begin{pmatrix}
     1.0 \\
     1.0
     \end{pmatrix} \, .

Thus we are able to recover our true amplitudes, :math:`\vec{a}`, from
(i) the response of our wavelet filters to the signal, and
(ii) the known response coefficient matrix.

Note that the RCM problem becomes trivial for
:math:`r_{S_{2}, \, M_{2}} = r_{M_{2}, \, S_{2}} \approx 0.0`,
where the filters do not respond to energy at the adjacent tidal frequency.



Additional reading
~~~~~~~~~~~~~~~~~~~~~~~~~
- See `Lobo et al., (2024) <https://journals.ametsoc.org/view/journals/atot/41/10/JTECH-D-23-0144.1.xml>`_
  for details on the information presented on this page.


