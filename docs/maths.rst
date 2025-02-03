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
The steps to a convolution are:

- Line up :math:`S_{1}(t)` and :math:`S_{2}(t)` so that they overlap
  by one point in time
- Multiply every point in time and sum (note that the only non-zero parts
  of the sum are where the two signals overlap)
- Shift :math:`S_{2}(t)` one point to the left, so that the two signals
  overlap at two points in time and sum
- Repeat this process until :math:`S_{2}(t)` is on the left side
  of :math:`S_{1}(t)`


Importantly, the convolution of a sine wave at given frequency, :math:`f`,
with a signal, :math:`S(t)`, will return the amplitude of the sine wave
with frequency :math:`f` that contributes to the linear superposition of
sine waves that can be summed to reconstruct :math:`S(t)` from said sine waves.
For example, in the image above, a convolution of a sine wave at the frequency of the
red line with the black line would return the amplitude of the red line.


However, we are still left without information on the phase of the sine wave that
contributes to the signal.
For this, we turn to complex numbers.



Complex numbers
~~~~~~~~~~~~~~~~~~~~~~~~~
We define a complex number as

   .. math::
   C \equiv a + i \, b \, ,

where :math:`a` is the *real* component,
:math:`b` is the *imaginary* component,
and :math:`i = \sqrt{-1}` is the
*imaginary number*.
The power of complex numbers, arguably, lies in their ability to
represent sine waves.





A wavelet
~~~~~~~~~~~~~~~~~~~~~~~~~









A spectrogram
~~~~~~~~~~~~~~~~~~~~~~~~~









Response coefficient matrix
~~~~~~~~~~~~~~~~~~~~~~~~~









CWT_Multi filter bank
~~~~~~~~~~~~~~~~~~~~~~~~~









Additional reading
~~~~~~~~~~~~~~~~~~~~~~~~~



