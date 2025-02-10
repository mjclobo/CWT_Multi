=========================================
A basic use case of CWT_Multi
=========================================

Having provided the user with a brief introduction
to the theory behind CWT_Multi, we now turn to the
fun part: how to use the program!

This code is almost exactly the same a that found in the file
CWT_MWE.m, which can be found in the CWT_Multi Github repository.

Firstly, the user adds the path to the CWT_Multi source code

.. code-block:: language

   addpath('./src')


A time series
~~~~~~~~~~~~~~~~~~~~~~~~~
We then load the data, which is also included in the main repository.

.. code-block:: matlab

   load('./data/Vancouver/vancouver_wl.mat')

For this example, we use a time series at NOAA Vancouver station number 9440083.
We now choose the time windoe over which we would like to perform the analysis.
With the default decimation factor of twenty time steps, and a sampling period of
one hour, two years of data will output a time series of tidal amplitudes and phases
that is about :math:`365.25 \times 2 \times (20/24) / 1 = 608` values long.
We define this window as

.. code-block:: matlab

    startDate = datetime(2010,06,03);
    endDate = datetime(2012,06,03);

    Van.start = find(Van.dates==startDate);
    Van.end = find(Van.dates==endDate);

    Van.dates = Van.dates(Van.start:Van.end);
    Van.wl = Van.wl(Van.start:Van.end);

    % set time zone
    Van.dates.TimeZone = 'America/Los_Angeles';


The main routine
~~~~~~~~~~~~~~~~~~~~~~~~~
The minimum amount of arguments to CWT_Multi is two: the timestamps in DateTime format
and the corresponding water level data.

.. code-block:: matlab

    [constits,species,~,~] = cwtMulti(Van.dates,Van.wl) 

There are many optional arguments to the routine, which are detailed in the header of the
cwtMulti.m file header.

The main output
~~~~~~~~~~~~~~~~~~~~~~~~~
The main output from CWT_Multi are two struct type variables:

- :code:`constits` returns constituent amplitudes for the respective :math:`D_{1}`, :math:`D_{2}`, and
  :math:`D_{3}+` analyses.
- :code:`species` returns the species amplitudes for the :math:`D_{1}` species, :math:`D_{2}` species, etc.

To access and plot the main :math:`D_{2}` constituents, for example, use

.. code-block:: matlab
    
    figure(); p=tiledlayout(1,1);

    ax0=nexttile;
    plot(constits.decTimesAll,constits.M2.amps,'linewidth',2,'DisplayName','M_2')
    hold(ax0,'on')
    plot(constits.decTimesAll,constits.S2.amps,'linewidth',2,'DisplayName','S_2')
    plot(constits.decTimesAll,constits.N2.amps,'linewidth',2,'DisplayName','N_2')
    hold(ax0,'off')
    grid on
    legend('location','northeast')
    title('Vancouver D_2', 'Units', 'normalized', 'Position', [0.5, 0.85, 0])
    set(gca,'xtick',[])
    ylabel(gca,'Amp. (m)')
    xlim([constits.decTimesAll(1) constits.decTimesAll(end)])
    ax = gca;
    ax.FontSize = 14;

.. image:: /images/Van_Jul2010_Jul2012.png
   :width: 600pt

  
A one-liner to plot the reconstruction for the constituents analysis is

.. code-block::

    plot(constits.alltimes,constits.reconstruction.reconAll,'linewidth',2,'DisplayName','constituent-based reconstruction')




