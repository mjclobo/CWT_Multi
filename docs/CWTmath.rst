=========================================
Basic CWT_Multi theory
=========================================

On this page we provide the basic theory required
to understand how CWT_Multi, and wavelet transforms more generally,
work.
We present ideas and illustrations, rather than comprehensive maths.
Future pages in this documentation will cover the mathematic details.
For now, we refer the reader to Lobo et al. (2024), and include
additional references at the bottom of this page.


Tidal species and constituents
~~~~~~~~~~~~~~~~~~~~~~~~~
Different *tidal constituents* describe unique frequencies that correspond
to different sets of Doodson numbers.
Doodson numbers describe the fundamental forcing temporal periods of the tides.
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
and :math":`S_{2}` tidal constituents with :math:`R_{C}=1`,
we need to use an analysis window length of

    .. math::
    L_{w} 
    &= \left | f_{M_{2}} - f_{S_{2}} | ^{-1} \\
    &= \left | 1/12.4206012 - 1/12 \right | ^{-1} \ \mathrm{hr} \\
    & \approx 15 \ \mathrm{days} \, . 

Another way to word the Rayleigh criterion is:
if a signal is composed of sine waves at two different frequencies
we need the sine wave at the higher frequency to complete
at least one more cycle than the signal at the lower frequency,
in order to be able to tell the two frequencies apart.
This is shown graphically below.

.. image:: /images/RC_Lw.png
   :width: 300pt

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
lengths: :math:`L_{w} / 2`, :math:`L_{w}`, and :math"`2 \, L_{w}`.
We find that with the shortest window we are not able to differentiate between
energy at the two frequencies (red line).
Once we analyze a signal that is at least the length :math:`L_{w}`,
we are able to resolve energy at the two frequencies (green line).

.. image:: /images/RC_spectra.png
   :width: 300pt

Note, however, that as the analyzed signal gets longer,
the peaks at the two frequencies become more distinct (yellow line).
If we had an infinitely long signal, the energy at the two frequencies would be represented by
vertical lines (hence the often-used term *line spectra*).
The apparent "spreading" of energy at frequencies around
:math:`M_{2}` and :math:`S_{2}` is an artifact of the finite-length
analysis window.


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
   :width: 300pt

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

Frequency response
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


Response coefficient matrix
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
Consider a signal that only has energy at the :math:`M_{2}` and :math:`S_{2}`
frequencies, where we would like to define the :math:`M_{2}` and :math:`S_{2}`
amplitudes as a function of time.
We thus define the RCM as

   .. math::
    \boldsymbol{R} =
    \begin{pmatrix}
    r_{M_{2}, \, M_{2}} & r_{M_{2}, \, S_{2}} \\
    r_{S_{2}, \, M_{2}} & r_{S_{2}, \, S_{2}}
    \end{pmatrix} \, ,

where :math:`r_{f_{1}, \, f_{2}}` describes the fre





From these definitions we can describe the RCM problem in words:
the true amplitudes are 



Additional reading
~~~~~~~~~~~~~~~~~~~~~~~~~



