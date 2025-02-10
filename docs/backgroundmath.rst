=========================================
Background math
=========================================

On this page we provide the math background required
to understand how CWT_Multi, and wavelet transforms more generally,
work.
We present ideas and illustrations, rather than comprehensive maths.
Future pages in this documentation will cover the mathematical details.
For now, we refer the reader to Lobo et al. (2024), and include
additional references at the bottom of this page.


The following sections aim to be self-contained and therefore contain
information of increasing level of complexity.
The reader is encouraged to scan the material until they
come to a section that contains information unfamiliar
to them.


Sinusoids and complex exponentials
~~~~~~~~~~~~~~~~~~~~~~~~~
A *sine wave* is a signal that takes the form

   .. math::
    s(t) = A \mathrm{sin} ( \omega t + \phi ) \, ,

where

- :math:`A` is the *amplitude* of the wave with the units of your signal (e.g., meters for water level data)
- :math:`\omega = 2 \pi f` is the *angular frequency* with units rad/s
- :math:`f` is the *frequency* with units of cycles/s, i.e., Hz
- :math:`t` is time with units of seconds
- :math:`\phi` is a phase offset

The total argument to the sine function is known as the *phase*,
a unitless quantity that determines the output from the sine function
(a values in the range :math:`-1` to :math:`1`).

The Fourier transform, in pictures
~~~~~~~~~~~~~~~~~~~~~~~~~
The *Fourier transform* provides the means to decompose a
signal into a linear superposition (i.e., a sum) of sine
waves.
In particular, a Fourier transform routine will provide the user
with the amplitude, frequency and phase offset for each sine wave
that contributes to the original signal.
Though there are many technical aspects to both the mathematics
and the practical application of Fourier transforms, the information
presented so far is sufficient for our purposes.


.. image:: /images/FT_drawing.png
   :width: 300pt



Convolution
~~~~~~~~~~~~~~~~~~~~~~~~~
One of the main math tools needed to compute a Fourier transform
is the *convolution*.
The ingredients of a convolution are two signals that are evenly-spaced
in time, which we call :math:`S_{1}(t)` and
:math:`S_{2}(t)` here for convenience.
Generally one signal is your signal of interest, e.g., a water level time series,
while the other signal is a piece of some analytic method, e.g., a sine wave.
A convolution is a function of time.

The steps to a convolution are:

- Line up :math:`S_{1}(t)` and a flipped version of :math:`S_{2}(t-\tau)` so that they overlap
  by one point in time; specifically, they overlap for the last point in time
  where you have a value for :math:`S_{1}(t)` (:math:`\tau`, loosely, tracks by how many points the two signals overlap) 
- Multiply every point in time where there is overlap and sum (note that the only non-zero parts
  of the sum are where the two signals overlap)
- Shift :math:`S_{2}(t - \tau)` one point to the left, so that the two signals
  overlap at two points in time and sum
- Repeat this process until :math:`S_{2}(t - \tau)` is on the left side
  of :math:`S_{1}(t)`, i.e., :math:`\tau > t_\mathrm{max}`
- You should now have a value for the convolution at every point in time where the two
  signals overlapped.

Complex numbers
~~~~~~~~~~~~~~~~~~~~~~~~~
We define a complex number as

    .. math::
     C \equiv a + i \, b \,  ,

where :math:`a` is the *real* component,
:math:`b` is the *imaginary* component,
and :math:`i = \sqrt{-1}` is the
*imaginary number*.
The power of complex numbers, arguably, lies in their ability to
represent the amplitude and phase of a signal.
In particular, we will associate a complex number with a given frequency,
such that we know the amplitude and phase of the signal at a given frequency,
similar to the Fourier transform, but as a function of time.

The equations for the amplitude and phase of a complex number are
found easily online.
We encourage the reader to view `this video <https://www.youtube.com/watch?v=YVLEsxq2kEA>`,
which illustrates exactly how a complex number represents a signal's amplitude and
phase.

A wavelet
~~~~~~~~~~~~~~~~~~~~~~~~~
In brief, wavelet analysis provides information on the amplitude
and phase of a signal as a function of both frequency and time.
This is different than Fourier analysis, where the latter assumes
stationarity of the amplitude and phase, i.e., there is only one amplitude
and phase value at a given frequency for the signal being analyzed.


Additional reading
~~~~~~~~~~~~~~~~~~~~~~~~~
- We recommend `this 3Blue1Brown video <https://www.youtube.com/watch?v=spUNpyF58BY>`
  for an intuitive introduction to the Fourier Transform.
- Jonathan Lilly has great `course material <http://jmlilly.net/course/index.html>`
  for more details on wavelet analysis and
  signal processing, more generally.


