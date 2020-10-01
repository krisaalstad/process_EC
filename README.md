# process_EC
An explicit (white box) Eddy Covariance processing module written entirely in Fortran(90). 

This module was the main computational tool used in my master thesis (Aalstad, 2015) that I succesfully defended in June 2015 at the Section for Meteorology and Oceanography in the Department of Geosciences at the University of Oslo. The aim of this thesis was to look into the challenges of applying the eddy covariance technique to measure turbulent fluxes at the northernmost FLUXNET site at Bayelva (79○N) near Ny-Ålesund (Svalbard, Norway). This site was operated by Dr. Julia Boike and colleagues at AWI for more than a decade and has led to many publications (e.g. Westermann et al. 2009; Lüers et al. 2014). I only analyzed a few months of the raw 20 Hz time series, but the outputs of my module are arguably more in depth than what is usual from many popular plug&play eddy covariance modules. For example, the outputs include power spectra, cospectra, ogives, and structure functions for most of the measured turbulent variables (u,v,w,T_v,rho_v) and their products. The outputs also include detailed qualitfy flags and uncertainty estimates. The turbulent flux estimates from my (nameless) module were in close agreement to those obtained from TK2 (Mauder and Foken, 2004) which is a "gold standard", widely used, and standardized eddy covariance processing package that is only available as an executable (i.e. source code is not provided).

It is worth noting that I wrote my module several years ago (between September 2014 and May 2015), and it has not received any edits or attention by me since then. If time/funding permits, I may decide to pick up this module to update and/or generalize it. After several years, I have decided to archive this code here in the hope that it could be useful for someone (including myself) at some point. I make no claim to as to the optimality of this module and there are probably some bugs in the >7000 lines of code. Nonetheless, the code is very explicit, so it should be (relatively) easy to follow the different processing steps by comparing the module (structure.f90 in src) to the steps outlined in Aalstad (2015). Being written in Fortran, it is also quite fast; around 100 times faster than a comparable implementation in an interpreted language. 

# Run requirements

Running: First compile the module using the Makefile. To do this cd into the obj directory and type "make" (without the quotes :D) in the terminal. This will generate an executable, called run which you can execute in the terminal by entering the command "./run". If everything works as intended, this will generate ouput netcdf and logfiles after the code is finished. 

Analysis: The processed output in the netcdf files can be easily analyzed and visualized in an interpreted language such as MATLAB or Python.

Navigation: The entire module is coded up in the structure.f90 script in the src directory. The "steering" script main.f90, also in the src directory, is used to run the routines in the module and specify run options (such as input and output files). Example input files are found in the input directory. Note that the actual raw input eddy covariance data needed for actually running this module for a real case is not provided as it is too large to store on github (a truncated example of the expected input is. Example output files are found in the output directory. 

Programming language: Fortran90.

Operating system: The module has only been tested using Red Hat Enterprise Linux (RHEL) version 6; I'm quite confident it will also work with newer versions. It should also work on other Linux distributions such as Ubuntu or Fedora. I hope that after some edits it would be possible to also get it to run on a Mac or Windows based operating system. 

Compiler: The module has only been compiled with gfortran. It is possible that it could work with another compiler.

Libaries: The only external libraries required are netcdf (see: https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90_The-NetCDF-Fortran-90-Interface-Guide.html) and DFFTPACK (a double precision clone of FFTPACK: see dp.tgz in http://www.netlib.org/fftpack/). 


# Rerences

Aalstad, K. (2015): Applying the Eddy Covariance Method Under Difficult Conditions. Master thesis submitted for the MSc in Geoscience degree (Specialization: Meteorology and Oceanography) at the University of Oslo. Supervisor: Terje K. Berntsen. The best formatted version of this thesis is available through doi:10.13140/RG.2.1.3083.7844 (original version: https://www.duo.uio.no/handle/10852/45561). 

Lüers et al. (2014): Annual CO2 budget and seasonal CO2 exchange signals at a high Arctic permafrost site on Spitsbergen, Svalbard archipelago. Biogeosciences. doi:10.5194/bg-11-6307-2014

Mauder and Foken (2004): Documentation and Instruction Manual of the Eddy Covariance Software Package TK2. Available at: https://core.ac.uk/download/pdf/33806389.pdf. Note: A newer version, TK3, is also available now, but we do not expect major differences between versions.

Westermann et al. (2009): The annual surface energy budget of a high-arctic permafrost siteon Svalbard, Norway. The Cryosphere. doi:10.5194/tc-3-245-2009
