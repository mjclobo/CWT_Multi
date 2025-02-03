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

.. figure:: /images/FT_drawing.pdf
   :width: 300pt

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




Convolution
~~~~~~~~~~~~~~~~~~~~~~~~~











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



