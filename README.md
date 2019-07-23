# VelFieSmo

Smoothing of Galactic Velocity Fields, a Summer 2019 Gradmap Study


## Installation

The following packages will be needed

     * python: should have astropy, specutils, pyspeckit. Usually via pip. There are also conda channels.
     * NEMO
     * ds9


## Examples

The various versions of cubespectrum*.py are a mess, as they are experiments using different fits access methods and more.
One of the goals is to make a nice clean version of this script.

    ./cubespectrum12.py ../data/ngc6503.fits
    ./cubespectrum12.py ../data/ngc6503.fits  244 182
    ./cubespectrum12.py ../data/ngc6503.fits  161 128
    ./cubespectrum12.py ../data/ngc6503.fits  77 85

    ./cubespectrum13.py ../data/ngc6503.fits



The main modeling script is *mkgalcube* and the models/Makefile has example how to create a series of models.
