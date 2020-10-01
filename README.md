# process_EC
An Eddy Covariance processing module written entirely in Fortran(90). 

This module was the main computational tool used in my master thesis (Aalstad, 2015) that I succesfully defended in June 2015 at the Section for Meteorology and Oceanography in the Department of Geosciences at the University of Oslo. The aim of this thesis was to look into the challenges of applying the eddy covariance technique to measure turbulent fluxes at the northernmost FLUXNET site at Bayelva (79○N) near Ny-Ålesund (Svalbard, Norway). This site was operated by Dr. Julia Boike and colleagues at AWI for more than a decade and has led to many publications (e.g. Westermann et al. 2009; Lüers et al. 2014). I only analyzed a few months of the raw 20 Hz time series, but the outputs of my module are arguably more in depth than what is usual from many popular plug&play eddy covariance modules. For example, the outputs include power spectra, cospectra, ogives, and structure functions for most of the measured turbulent variables (u,v,w,T_v,rho_v) and their products. The outputs also include detailed qualitfy flags and uncertainty estimates. 

It is worth noting that I wrote the module several years ago (between September 2014 and May 2015), and it has not received any edits or attention since then. If time/funding permits, I may decide to pick up this module to update and/or generalize it. After several years, I have decided to archive this code here in the hope that it could be useful for someone (including myself) at some point. I make no claim to as to the optimality of the code and there are probably some bugs in the >7000 lines of code. Nonetheless, the code is very explicit, so it should be (relatively) easy to follow the different processing steps by comparing the module (structure.f90 in src) to the steps outlined in Aalstad (2015). Being written in Fortran, it is also quite fast; around 100 times faster than a comparable implementation in an interpreted language (e.g. MATLAB). 

# Run requirements

How to run: The 

Programming language: Fortran90.

Operationg systems: The module has only been tested using Red Hat Enterprise Linux (RHEL) version 6; I'm quite confident it will also work with newer versions. It should also work on other Linux distributions such as Ubuntu or Fedora. I hope that after some edits it would be possible to also get it to run on a Mac or Windows based operating system. 

Compiler: The module has only been compiled with gfortran. It is possible that it could work with another compiler.

Libaries: The only external libraries required are netcdf (see: https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90_The-NetCDF-Fortran-90-Interface-Guide.html) and DFFTPACK. 


# Rerences

Aalstad, K. (2015): Applying the Eddy Covariance Method Under Difficult Conditions. Master thesis submitted for the MSc in Geoscience degree (Specialization: Meteorology and Oceanography) at the University of Oslo. Supervisor: Terje K. Berntsen. The best formatted version of this thesis is available through doi:10.13140/RG.2.1.3083.7844 (original version: https://www.duo.uio.no/handle/10852/45561). 

Lüers et al. (2014): Annual CO2 budget and seasonal CO2 exchange signals at a high Arctic permafrost site on Spitsbergen, Svalbard archipelago. Biogeosciences. doi:10.5194/bg-11-6307-2014

Westermann et al. (2009): The annual surface energy budget of a high-arctic permafrost siteon Svalbard, Norway. The Cryosphere. doi:10.5194/tc-3-245-2009
