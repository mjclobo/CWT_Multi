Welcome to CWT_Multi's documentation!
===================================

**CWT_Multi** is a set of MATLAB routines that
perform a version of continuous wavelet transform (CWT) analysis
that is specialized for application to nonstationary tidal data.

See below to navigate to pages that detail program installation, background maths,
the theory of CWT_Multi, and how to use CWT_Multi in practice.

If the reader is already knowledgeable of basic signal processing
subjects such as complex sinusoids, convolution, and basic wavelet analysis,
we suggest that they skip to the section titled **Basic CWT_Multi theory**.

.. note::

   This project is under active development.

Introduction
--------
Signal processing tools fall into two broad categories: analysis and synthesis.
CWT_Multi provides analysis output in the form of tidal amplitudes and phases that
vary in time.
In addition, CWT_Multi reconstructs a time series of water level from this amplitude
and phase data that is assumed to exclusively include tide-related signal.
For some users, this tide-related signal will be of greatest interest.
For others, they can then subtract this tidal time series from their original data set
in order to focus on their signal of interest, e.g., storm surges.

Though the default parameters in CWT_Multi have proven sufficient for a variety of use-cases,
we encourage any reader that us seriously interested in using CWT_Multi to read through the documentation provided
here, taking specific note of the information in the *Basic CWT_Multi theory* section.

Any questions, whether related to the program or related to this documentation, may be
sent to *mattlobo@princeton.edu*.

Contents
--------

.. toctree::

   install
   backgroundmath
   CWTmath
   CWTex
   
