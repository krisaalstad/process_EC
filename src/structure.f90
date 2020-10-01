! By Kristoffer Aalstad 2015 (Last revision 01.05.2015).

! -------------------------------------------------------------------------------------------------------------------------------

! An eddy covariance processing module which subsequently:

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reads & Allocates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 1) Reads in raw eddy covariance data (t,u,v,w,T,h) from a TOA5-xxxx.dat CR3000 type ASCII file.
! 2) Allocates the raw data to arrays defined in the derived data type 'eddycov'.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Screens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 3) Performs three despking routines: first setting physical limits then a median absolute deviation test (Mauder et al. 2013),
!	and finally calculates skewness and kurtosis for detrended block sonic temperature and vertical velocity checking values
!	against thresholds and flagging accordingly.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rotations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 4) Orients the sonic anemometer, computing wind directions while identifying instances where these fall within the closed sector.
! 5) Applies the planar fit coordinate rotation technique of Wilczak et al. (2001) if the series is longer than one week.
! 6) Rotates the plane winds into the natural ensemble streamline frame such that block averaged v = 0 and block averaged u = plane mean wind.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Supplement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 7) Reads in measurement height from SR50/BCS input files.
! 8) Reads in the relevant diurnal pressure field from an external netcdf file for the NYÅ met station; |p'|/p_0 << |d'|/d_0 where
! d is density is assumed. This allows use of slow meteorology for the allocation of pressure, the value of which is further assumed
! to be the same at NYÅ and the measurement site. Required for the WPL and SND corrections. Note that the input filename here
! needs to be changed if the code is to be used on another data set/period. The current pressure field runs from 01.01.2007 to
! 31.12.2010 with 4 entries per day (00:00 UTC, 06:00 UTC, 12:00 UTC and 18:00 UTC).
! 9) Calculates the air density using pressure and the absolute humidity, converts absolute to specific humidity and
! converts virtual temperature to temperature.


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flux Corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 10) Corrects for sensor separation by maximizing the square cross correlation between the LDT vertical velocity and absolute humidity.
! 11) Corrects for cospectral attenuation at high and low frequencies using the analytical method of Massman (2000).
! 12) SND correction, converting buoyancy heat flux to sensible heat flux.
! 13) WPL correction, diagnoses the block average vertical velocity based on the assumption of zero dry air mass flux.
! 14) Iterates corrections 11,12, and 13 until the flux varies by less than 1% from one iteration to the next.
! + FFSE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 11) Applies quality flags following Foken&Wichura (1996), Vickers&Mahrt (1997) and Mauder et al. (2013) for 30 minute block windows:
! Flags: 0=high quality, 1=medium quality, 2=bad quality (discarded).
! Routines:
!			-> Percentage of data consisting of spikes. If >10% Flag 2.
!			-> Checks for flow distortion based on wind direction and constancy ratio and hard flags accordingly.
!			-> Verifies that the block average vertical velocity satisfies |w|<0.15, otherwise hard flag (|w|>0.1 soft flag).
!			-> Stationarity test using a coefficient 's'. If |s|<=0.3 Flag 0, if 0.3<|s|<=1.0 Flag 1, if |s|>=1.0 Flag 2. 
!			-> Tests various vertical velocity integral turbulence characteristics for unstable and near neutral stabilities; If
!				|dITC|<=0.3 Flag 0 elseif |dITC|<=0.75 Flag 1 else Flag 2.
!			-> Combined all the individual flag categories (including skewness-kurtosis) into an overall block flag which takes
!			the maximum value of all the individual flags for a given block. If hard flagged (f=2) the block is discarded from
!			further analysis.
!	
! 12) Calculates the relative flux sampling uncertainty by combining the methods of Finkelstein and Sims (2001) and Bilesbach (2011).

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 12) Computes the 30 minute block averaged statistics including fluxes of momentum, buoyancy as well as sensible and latent heat.
! 13) Calculates autostatistics, the autocorrelations and the second order structure functions, for a user defined lag increment
! and segment length.
! 14) Produces (co)spectra of turbulence as well as Ogives by conditioning, fast fourier transform, folding, dealiasing and smoothing.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! 15) Writes the output arrays to a netcdf file with the same name as the input file but a .cdf extension.
! 16) Writes the terminal dialogue (warnings and messages) to a logfile with the same file name as the input file but 
! a log.txt extension.

! N.B. Processing about 6 months of data takes ~12 hours.

! -------------------------------------------------------------------------------------------------------------------------------

MODULE structure
	USE netcdf
	IMPLICIT NONE
	TYPE, PUBLIC				:: eddycov		! Derived data type.

		! Parameters to be specified by the user of the module. These include options.
		CHARACTER(LEN=80)		:: filename		! EC input ASCII-file name.

		CHARACTER(LEN=80)		:: outfilename		! Output netcdf-file name.

		CHARACTER(LEN=80)		:: logfilename		! Output logfile.

		CHARACTER(LEN=80)		:: pfile		! Input pressure file.

		CHARACTER(LEN=80)		:: tfile		! Temperature limits input file.

		CHARACTER(LEN=80)		:: zfile		! Input snoic height file name (OPTIONAL).
		
		INTEGER				:: latitude		! Latitude of site (nearest degree north).
		
		INTEGER				:: fs			! Sampling frequency (in Hz).

		INTEGER				:: avper		! Block averaging period (in minutes).

		INTEGER				:: spwin		! Spectral window (in minutes).

		INTEGER				:: drwin		! D & R, structure function  & autocorrelation, window (in minutes).

		REAL*8				:: lag1			! Initial time lag (sec) desired in D&R calculations, ensuing lags will be integer
									! products of this lag up to half the length of the window.

		LOGICAL				:: hfstat		! Status flag specifying inclusion (T) or exclusion (F) of high frequency data.

		LOGICAL				:: smooth		! Status flag specifying smoothing (T) or no smoothing of turbulent spectra.
									! Smoothing is highly recommended to reduce the size of output NetCDF files.

		LOGICAL				:: detrend		! Logical specifying if windows should be detrended prior to spectral analysis and
									! structure function computation.

		LOGICAL				:: taper		! Logical specifying if a bell taper should be applied as a filter to the
									! series prior to computing the fourier transform. 
		CHARACTER(LEN=2)		:: tapertype		! Type of taper to use ('B': Bell, 'H': Hamming).

		REAL*8				:: CSATbearing		! Bearing that the sonic is pointing towards, for wind direction calculations.

		REAL*8				:: LIbearing		! Bearing that the IRGA is pointing towards, for the cross-correlation technique.

		REAL*8				:: dsep			! Horizontal separation distance between the IRGA and the CSAT (in meters).

		INTEGER,ALLOCATABLE		:: disdir(:,:)		! Wind direction bins to remove due to flow distortion by upwind structures.
									! New row for each bin, with the columns representing the bounsd of the bin.
									! For example if disdir(1,:)=[100 105] the bin will discount any block where
									! the mean wind direction is in the range 100^deg <= MWD <= 105^deg.
		
		REAL*8				:: zbase		! Height of the sonic AGL (m).

		LOGICAL				:: hmes			! Logical specifying if changes in the height of the sonic due to snowfall or 										! sensor adjustment have occured. Should only be set to .T. if height measurements
									! are available.


		! Parameters for subsequent writing to logfiles within each subroutine.
		INTEGER				:: loglun=21	 	! logical unit number for the ASCII-logfile.
		INTEGER				:: logres	 	! result of file operation and end of file.
		

		! Calculated parameters from the above in the dimensions routine.
		REAL*8				:: dt			! Sample period (in seconds).
		INTEGER				:: blockl		! Averaging period length (number of points in time).
		INTEGER				:: spwinl		! Spectral window length (number of points in time).
		INTEGER				:: drwinl		! D & R window length (number of points in time).
		


		! Internal variables passed to and from the subroutines.
		INTEGER				:: nn						! Number of points in time.
		INTEGER				:: nvars=6					! Number of variables of interest in the input file.
		REAL*8				:: nanreal=-1e3					! NaN placeholder value.
		REAL*8				:: alpha					! Angle between due east and sonic x-axis.
		REAL*8, ALLOCATABLE		:: zmes(:)					! Measuring height as a function of time only allocated
												! IF hmes==True. (i.e. if input data is available).
		REAL*8, ALLOCATABLE		:: t(:,:)					! Time array (yyyy,mm,dd,hh,mm,ss,ms).
		REAL*8,	ALLOCATABLE		:: u(:),v(:),w(:),Ts(:),h(:),p(:) 		! Time series of measured variables considered (will be
												! despiked and planar fitted and flagged).
		REAL*8, ALLOCATABLE		:: unode(:),vnode(:),wnode(:),Tsnode(:),hnode(:)! Time series of measured variables considered (only 													despiked).
		REAL*8,	ALLOCATABLE		:: rho(:),q(:),Tv(:),TT(:)			! Auxillary variables calculated from the above.
		REAL*8, ALLOCATABLE		:: bt(:,:,:)					! Block timestamps.
		LOGICAL				:: spactive,dractive				! Logicals declaring if at least one window is active.
		REAL*8				:: cfH=SQRT(2.52),cfB=SQRT(1.16)		! Compensation factors for the tapers in spectral analysis.
												! restoring lost variance.
		
		! Unit vectors used in rotation routines.
		REAL*8,DIMENSION(3)		:: iunit,junit,kunit				! Unit vectors defining the planar frame: kunit
												! is normal to the streamline plane, iunit
												! and junit are in the streamline plane and have
												! components pointing east and north respectively.
		REAL*8,DIMENSION(3)		:: sepunit					! Separation unit vector, pointing from the sonic anemometer
												! to the open path IRGA.

		! Output variables, writen to netcdf file.
		REAL*8, ALLOCATABLE		:: su(:,:),sv(:,:),sw(:,:),sT(:,:),sh(:,:)	! Power spectra of the variables.
		REAL*8, ALLOCATABLE		:: sobukhovl(:),sz(:),sustar(:),sTstar(:)	! L, zmes, ustar and Tstar for the respective spectra.
		REAL*8, ALLOCATABLE	        :: meanuh(:)					! Mean 'horizontal' wind speed in the window.
		REAL*8, ALLOCATABLE		:: couw(:,:),covw(:,:),coTw(:,:),cohw(:,:)	! Cospectra of the variables.
		REAL*8, ALLOCATABLE		:: oguw(:,:),ogvw(:,:),ogTw(:,:),oghw(:,:)	! Ogives of the variables considered.
		REAL*8,	ALLOCATABLE		:: oguu(:,:),ogvv(:,:),ogww(:,:)	
		REAL*8,	ALLOCATABLE		:: ogTT(:,:),oghh(:,:)
		REAL*8, ALLOCATABLE		:: quuw(:,:),quvw(:,:),quTw(:,:),quqw(:,:)	! Quadrature spectra of the variables considered.
		REAL*8, ALLOCATABLE		:: d2uu(:,:),d2vv(:,:),d2ww(:,:),d2TT(:,:),d2qq(:,:)	! 2nd order structure functions
		REAL*8, ALLOCATABLE		:: d3uu(:,:),d3vv(:,:),d3ww(:,:),d3TT(:,:),d3qq(:,:) 	! 3rd order structure functions.
		REAL*8, ALLOCATABLE		:: d4uu(:,:),d4vv(:,:),d4ww(:,:),d4TT(:,:),d4qq(:,:) 	! 4th order structure functions.
		REAL*8, ALLOCATABLE		:: ruu(:,:),rvv(:,:),rww(:,:),rTT(:,:),rqq(:,:)	! Autocorrelation.
		REAL*8, ALLOCATABLE		:: drobukhovl(:),drz(:),druh(:)			! L,zmes and mean horz wind for autostatistics.
		REAL*8, ALLOCATABLE		:: tlag(:)
		REAL*8, ALLOCATABLE		:: fsp(:),fog(:)				! Frequency arrays.
		REAL*8, ALLOCATABLE		:: tspwin(:,:,:), tstwin(:,:,:)			! Time stamps for active windows.

		! Output block variables also writen to netcdf file.
		INTEGER,ALLOCATABLE		:: flag(:)							! Output simple combined block flag (0-2).
		INTEGER,ALLOCATABLE		:: flagnan(:),flagdir(:),flagstat(:),flagitcw(:),flagvert(:)	! Category flags.
		INTEGER,ALLOCATABLE		:: flagstatuw(:),flagstatTsw(:),flagstathw(:)	! Individual block stationarity flags.
		INTEGER,ALLOCATABLE		:: flagskewT(:),flagskeww(:),flagkurtT(:),flagkurtw(:),flagsk(:)! Kurtosis & Skewness flags.
		REAL*8,ALLOCATABLE		:: znoqc(:),obukhovlnoqc(:)			! Obukhovlength and mes. height prior to QC.
		REAL*8,ALLOCATABLE		:: zbar(:)							! Block average height.
		REAL*8,ALLOCATABLE		:: mwd(:),mws(:)						! Block average wind speed and direction.
		REAL*8,ALLOCATABLE		:: ubar(:),vbar(:),Tbar(:),wbar(:),qbar(:),pbar(:)			! Block averages.
		REAL*8,ALLOCATABLE		:: hbar(:),rhobar(:),Tvbar(:),Tsbar(:)					! -"-
		REAL*8,ALLOCATABLE		:: tblock(:,:,:)							! Block time stamps.
		REAL*8,ALLOCATABLE		:: uwcov(:),vwcov(:),Twcov(:),Tswcov(:),hwcov(:)			! Block covariances.
		REAL*8,ALLOCATABLE		:: Tvar(:),Tsvar(:),uvar(:),vvar(:),wvar(:),hvar(:),qvar(:)		! Block variances.
		REAL*8,ALLOCATABLE		:: obukhovl(:),tke(:),QE(:),QH(:),tau(:),ustar(:)			! Block fluxes.
		REAL*8,ALLOCATABLE		:: Fa_mom(:),Fa_lat(:),Fa_sen(:)					! Block flux attenuation factors.
		REAL*8,ALLOCATABLE		:: CFsep(:),CFSND(:),CFSND0(:),CFWPL(:)					! Flux correction factors.
		INTEGER,ALLOCATABLE		:: lagIRGA(:)								! Block lag index for the IRGA.
		REAL*8,ALLOCATABLE		:: niter(:)								! Number of iterations in flux
		REAL*8,ALLOCATABLE		:: csuw(:),csTw(:),cshw(:)						! Flux stationarity coefficients.
		REAL*8,ALLOCATABLE		:: cr(:)								! Constancy ratio.
		REAL*8,ALLOCATABLE		:: ffsetsw(:),ffsehw(:)							! Fractional flux sampling errors. 							
		
		
	CONTAINS

		! Combines each subroutine below, in the right order.
		PROCEDURE			:: doit => doitall
		
		! File operations.
		PROCEDURE			:: dimit => setdims	 ! Specifies dimensions based on user input.
		PROCEDURE			:: readf => readfromfile ! Read from raw data ASCII-file (TOA5-xxxx.dat)
		PROCEDURE			:: writef => writetofile ! Write processed data to netcdf file (TOA5-xxxx.nc)

	
		
		! Cleaning.
		PROCEDURE			:: limit => confinelims			! Physical limits algorithm (Vickers & Mahrt 1997)
		PROCEDURE			:: madit => madspike			! MAD despiking algorithm (Mauder et al. 2013).
		PROCEDURE			:: skit => skewkurt			! Set limits for the skewness and kurtosis of Ts and w
											! in each block and flag data accordingly (Vickers & Mahrt 1997).
		

		! Corrections & calculations.
		PROCEDURE			:: orient => orientnorth 		! Calculates block averaged wind direction and speed.
		PROCEDURE			:: fitit => planarfit			! Planar fit algorithm (Wilsczak et al. 2001).
		PROCEDURE			:: rotit => rotatestream		! Rotates into the plane mean wind (possible AFTER planar fit).
		! § The slow pressure from syncit is used in the remaining subroutines in this section; these perform conversions and corrections
		! an error is made in using the 'slow' pressure, but the error is considerably smaller than that without any corrections.
		PROCEDURE			:: syncp => syncpressure		! Synchronize a pressure array based on nearby eKlima data.
		PROCEDURE			:: syncz => syncsnowdepth		! Reads in the snow depth from the directory 'snow' IF hmes=True.
		PROCEDURE			:: diagit => diaganc			! Diagnoses ancillary variables (q,Tv,TT and rho).
		PROCEDURE			:: sepit => separation			! Corrects water vapor flux for sensor separation.
		PROCEDURE			:: attenit => attenuation		! Corrects all fluxes for flux attenuation.
		PROCEDURE			:: wplit => wplcorrection		! Corrects fluxes for density fluctuations (WPL).
		PROCEDURE			:: sndit => sndcorrection		! Corrects for the effects of moisture in the heat flux (SND).
		PROCEDURE			:: fluxit => fluxcorrections		! Applies all four flux corrections; the last three (attenuation,SND
											! & WPL) are applied iteratively.
		PROCEDURE			:: ustarit => frictionvelocity		! Calculates the friction velocity (ustar) given arrays
											! of the three velocity components. 
		PROCEDURE			:: olit=> obukhovlength			! Calculates the Obukhov length given arrays of
											! virtual potential temperature,vertical velocity and
											! friction velocity.		

		! Operations.
		PROCEDURE			:: barit => temporalmean		! Mean.
		PROCEDURE			:: devit => standarddeviation		! Standard deviation.
		PROCEDURE			:: covit => covariance			! Covariance.
		PROCEDURE			:: corrit => correlationcoefficient	! Correlation coefficient.
		PROCEDURE			:: detit => detrend			! Linear detrend.
		PROCEDURE			:: fillit => gapfill			! Fills gaps using linear interpolation.
		PROCEDURE			:: skewit => skewness			! Skewness.
		PROCEDURE			:: kurtit => kurtosis			! Kurtosis.
		PROCEDURE			:: statit => stationarity		! Stationarity coefficient (Foken & Wichura 1996).
		PROCEDURE			:: wITCit => wITC			! Vertical velocity integral turbulence characteristics (-"-).
		PROCEDURE			:: intermit => intermittency		! Intermittency coefficient (Mahrt et al. 1998).
		PROCEDURE			:: lagit => autoandstructure		! Autocorrelation and structure function (Stull 1988).	
		PROCEDURE			:: medit => median			! Median using a quicksort algorithm by Rew, J. (2003).
		

		
		! Spectral analysis.
		PROCEDURE			:: specit => spectrum          		! Fourier transform using the dfftpack (Swarztrauber 1985).
		PROCEDURE			:: crossit => crossspectrum		! Cross-spectral analysis.
		PROCEDURE			:: deali => dealias			! Low-pass filter applied in the frequency domain to
											! minimize aliasing (but not to Ogives).
		PROCEDURE			:: kosit  => kosmoothing		! Computes the smoothing matrix used
											! in the Konno-Ohmachi method for smoothing spectra on
											! a logarithmic scale.
		

		! Block calculations.
		PROCEDURE			:: blockerr => blockffse
		PROCEDURE			:: propit => propogate 
		PROCEDURE			:: blockf => blockflags
		PROCEDURE			:: blockit => blockvalues

		! Deallocation; used if multiple files are read through the module.
		PROCEDURE			:: deall => deallocateall
		
		
		


	END TYPE eddycov
	CONTAINS


		!---------------------------------------------------------------------------------------------------------------------
		! Run module.
		!---------------------------------------------------------------------------------------------------------------------
		
		! Runs all the ensuing subroutines in the intended order performing all the processing and QA/QC of the eddy covariance
		! data.
		SUBROUTINE doitall(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv

			PRINT*,'Details written to logfile=',ecv%logfilename
			PRINT*,'As well as the terminal.'
			! Open the logfile for subsequent write statements in the ensuing subroutines.
			OPEN(UNIT=ecv%loglun,FILE=ecv%logfilename,FORM="formatted",IOSTAT=ecv%logres)
 				IF(ecv%logres /= 0) THEN
 			PRINT *, "Error in opening output file, status: ", ecv%logres
 				STOP	! Stop the program.
		 	END IF



			CALL setdims(ecv)				! Setup the dimensions of the problem.

			CALL readfromfile(ecv)				! Read in the raw data from the desired TOA-______.dat file.

			CALL confinelims(ecv)				! Set physical limits and test for violation.

			CALL madspike(ecv)				! Perform the MAD despiking algorithm.

			CALL skewkurt(ecv)				! Confine higher-order moment limits.

			CALL orientnorth(ecv)				! Orient the sonic anemometer relative to north for the block wind directions.

			CALL planarfit(ecv)				! Apply the planar fit rotating into the (long term) mean streamlines.

			CALL syncpressure(ecv)				! Synchronize the slow pressure given in the pfile.

			CALL syncsnowdepth(ecv)				! Synchronize the measurement heights.

			CALL diaganc(ecv)				! Diagnose ancillary variables.

			CALL blockflags(ecv)				! Set the flags for the blocks and ignore the block if hard flag (flag=2) occurs.

			CALL blockffse(ecv)				! Calculate the block fractional flux sampling errors.

			CALL crossspectrum(ecv)				! Compute the cross-spectra after the QC.			

			CALL autoandstructure(ecv)			! Compute the autostatistics after the QC.		

			CALL blockvalues(ecv)				! Compute the block averages variables of interest after QC. All the flux
									! corrections are applied. Causing problems for thingy.
		
			CALL writetofile(ecv)				! Write desired processed output data to a netcdf file.

			CALL deallocateall(ecv)				! Deallocate the data type for the next pass.

			! Close the logfile.
			CLOSE(UNIT=ecv%loglun)
			
			RETURN
		END SUBROUTINE

		!----------------------------------------------------------------------------------------------------------------------
		! File operations
		!----------------------------------------------------------------------------------------------------------------------


		! Calculates the dimensions based on user specifications: number of points in time within each averaging block, spectral window
		! and D & R window.
		SUBROUTINE setdims(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv	
			INTEGER	      :: i
	              	CHARACTER(8)  :: date
        	      	CHARACTER(10) :: time
        	     	CHARACTER(5)  :: zone
        	      	INTEGER,DIMENSION(8) :: values
			! using keyword arguments
             		CALL DATE_AND_TIME(date,time,zone,values)
             		CALL DATE_AND_TIME(DATE=date,ZONE=zone)
             		CALL DATE_AND_TIME(TIME=time)
             		CALL DATE_AND_TIME(VALUES=values)

			ecv%dt=1.0/(1.0*ecv%fs)			! Sample period (in seconds).
			ecv%blockl=60*ecv%fs*ecv%avper		! Averaging period length (number of points in time).
			ecv%spwinl=60*ecv%fs*ecv%spwin		! Spectral window length (number of points in time).
			ecv%drwinl=60*ecv%fs*ecv%drwin		! D & R window length (number of points in time).

			! User specifications to the logfile.
			WRITE(ecv%loglun,*)'Eddy covariance processing module structure.f90.'
			WRITE(ecv%loglun,*)'Run on the:',date
			WRITE(ecv%loglun,*)'@:',time,zone
			WRITE(ecv%loglun,*)''
			WRITE(ecv%loglun,*)''
			WRITE(ecv%loglun,*)'////////////////////////////////////////////////////////////////////////'
			WRITE(ecv%loglun,*)'//////////////// Parameters ////////////////////////////////////////////'
			WRITE(ecv%loglun,*)"	Sampling frequency:",ecv%fs,' [Hz]'
			WRITE(ecv%loglun,*)"	Block averaging period:",ecv%avper,' [mins]'
			WRITE(ecv%loglun,*)" 	Sonic bearing:",ecv%CSATbearing,' [degrees from North]'
			WRITE(ecv%loglun,*)"	IRGA bearing:",ecv%LIbearing,' [degrees from North]'
			WRITE(ecv%loglun,*)"	Sensor separation:",ecv%dsep,' [m]'
			WRITE(ecv%loglun,*)"	Base sonic height:",ecv%zbase,' [m] AGL'
			DO i=1,SIZE(ecv%disdir,1)
				WRITE(ecv%loglun,*)"	Excluding winds from the bearing bin:",&
					ecv%disdir(i,1)," to ",ecv%disdir(i,2)," [degrees from North], due to upwind flow distortion" 
			END DO
			WRITE(ecv%loglun,*)"	Spectral window length:",ecv%spwin," [mins]"
			WRITE(ecv%loglun,*)"	Smallest resolved spectral period (1/f_ny):",2.0/(ecv%fs)," [s]"
			WRITE(ecv%loglun,*)"	Autocorrelation (R) and structure function (D) window length:",ecv%drwin," [mins]"
			WRITE(ecv%loglun,*)"	Base lag used in R&D:",ecv%lag1," [seconds]"
			WRITE(ecv%loglun,*)''
			WRITE(ecv%loglun,*)''
			WRITE(ecv%loglun,*)'////////////////////////////////////////////////////////////////////////'
			WRITE(ecv%loglun,*)'//////////////// Options    ////////////////////////////////////////////'
			WRITE(ecv%loglun,*)"	Monitored changes in sensor height AGL:",ecv%hmes
			WRITE(ecv%loglun,*)"	Write high frequency data to output file:",ecv%hfstat
			WRITE(ecv%loglun,*)"	Detrend spectra:",ecv%detrend
			WRITE(ecv%loglun,*)"	Taper	spectra:",ecv%taper
			WRITE(ecv%loglun,*)"	Taper	type:",ecv%tapertype
			WRITE(ecv%loglun,*)''
			WRITE(ecv%loglun,*)''


			RETURN
		END SUBROUTINE


		!// Reads in the EC data from a .dat (ASCII) file. A given file format, conformant with the CSAT3000 
		! data is assumed.
		! Entries before the first occurence of an exact hour/half hour entry are truncated to facilitate 30 minute
		! block-averaging. 
		SUBROUTINE readfromfile(ecv) 
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			! Internal declarations
			REAL*8, ALLOCATABLE,DIMENSION(:,:)	:: ttemp	! Temporary time array [yyyy mm	dd mi ss].
			REAL*8,	ALLOCATABLE,DIMENSION(:)	:: utemp, vtemp, wtemp, Tstemp, htemp ! Temporary series.
			INTEGER, DIMENSION(ecv%nvars)		:: posvar 	! Column position of a variable.
			REAL*8, DIMENSION(ecv%nvars-1)		:: temparray 	! Temporary array.
			CHARACTER(LEN=80)			:: temptime	! Temporary time string.
			CHARACTER(LEN=160), DIMENSION(ecv%nvars):: namevar 	! Variable name characters as appearing in the file.
			CHARACTER(LEN=160)			:: line,oldline ! Line entry.
			INTEGER					:: dlines 	! N lines of data.
			INTEGER					:: blockstart	! Occurence of the first hour/half hour entry.
			INTEGER, PARAMETER			:: hlines=4 	! N lines of headers.
			INTEGER					:: nanentries	! Number of nan entries.
			INTEGER, PARAMETER			:: nargs=8 	! N arguments in a column.
			CHARACTER(LEN=1), PARAMETER		:: deli=','	! Entry delimiter.
			CHARACTER(LEN=1)			:: sdeli	! String delimiter.
			INTEGER, PARAMETER			:: ilun=11 	! logical unit number for the ec-data file.
			INTEGER					:: res,eoff 	! result of file operation and end of file.
			INTEGER					:: i,j,k,l,m,n,o,ul,ll,npd 	! Counters.
			REAL*8					:: start,finish,nanreal 	! Timers and a NaN placeholder value.	


			PRINT*,'---------------------------------------------------------------------'
			PRINT*,'Reading from file:',ecv%filename
			WRITE(ecv%loglun,*) '---------------------------------------------------------------------'
			WRITE(ecv%loglun,*) 'Reading from file:',ecv%filename

                        
			
			! Allocation of variable names as these appear in the file.
			namevar(1)='"TIMESTAMP"'
			namevar(2)='"Ux"'
			namevar(3)='"Uy"'
			namevar(4)='"Uz"'
			namevar(5)='"Ts"'
			namevar(6)='"h2o"'

			!---------EC data file-----------------------------------------------------------------------------------
			! Open the file.
			OPEN(UNIT=ilun,FILE=ecv%filename,FORM='FORMATTED',&
				IOSTAT=res)
 			IF(res /= 0) THEN
				PRINT *, 'Error in opening b file, status:', res
  				STOP
 			END IF
			

			! Find the total number of lines in the file, as well as the
			! positions of any variables of interest by scanning the first
			! four header lines (i.e. m<4)
			eoff=0
			dlines=0
			i=0
			CALL cpu_time(start)
			DO WHILE(eoff.EQ.0)
				i=i+1
				READ(ilun,FMT='(A)',IOSTAT=res) line
				IF(i.LE.hlines) THEN
					m=INDEX(line,deli) ! Next occurence of ','
					IF(line(1:m-1).EQ.namevar(1))	THEN ! Check for timestamp occurence.
						k=1 ! Occurence counter.
						posvar(k)=k						
						line=line(m+1:) ! Truncate the line
						DO j=2,nargs
							IF(j.NE.nargs)	THEN ! Not last entry, ',' exists.
								m=INDEX(line,deli) ! Next occurence of ','
								IF(line(1:m-1).EQ.namevar(k+1))	THEN
									k=k+1
									posvar(k)=j
								END IF
								line=line(m+1:) ! Truncate the line.
							ELSE
								IF(line.EQ.namevar(k+1)) THEN
									k=k+1
									posvar(k)=j
								END IF
							END IF
						END DO			

					END IF
				END IF
				IF(line.EQ.oldline) THEN ! If a line is repeating itself we have reached the end of the file.
					CALL cpu_time(finish)
					eoff=1
					i=i-2
					PRINT '("-->Time elapsed scanning input file structure =",f5.2,"(mins)")',(finish-start)/60
					WRITE(ecv%loglun,*) '-->Time elapsed scanning input file structure =',(finish-start)/60,', minutes.'
				END IF
				oldline=line
			END DO
			dlines=i-hlines+1 ! Calculate the number of lines of data by excluding the header lines.

			!STOP

			! Rewind to the start of the file.
			REWIND(ilun)


			! Pre-allocate the temporary array dimensions.
			ALLOCATE(utemp(dlines))
			ALLOCATE(vtemp(dlines))
			ALLOCATE(wtemp(dlines))
			ALLOCATE(Tstemp(dlines))
			ALLOCATE(htemp(dlines))
			ALLOCATE(ttemp(dlines,6)) ! Y M D H M S


			! Read the header lines again so as to skip ahead to the data.
			DO i=1,hlines
				READ(ilun,FMT='(A)') line
			END DO
			
			! Assign the nan placeholder value to be an unrealistic real number.
			nanreal=ecv%nanreal
			PRINT'("-->NaN placeholder is nanreal=",I5)',NINT(ecv%nanreal)
			WRITE(ecv%loglun,*)"-->NaN placeholder is nanreal=",NINT(ecv%nanreal)
			! Count nanentries
			nanentries=0
			! Read the data lines (remaining lines) and allocate to arrays.
			CALL cpu_time(start)
			DO i=1,dlines ! Line counter.
				
				k=1 ! Variable counter.
				READ(ilun,FMT='(A)') line ! Read in the line as a character.
				sdeli='"'	! String delimiter.
				l=INDEX(line,sdeli) ! Position of the first occurence of '"'
				m=INDEX(line(l+1:),sdeli) ! -"- next '"'
				temptime=line(l+1:m) ! Allocate the temporary timestamp character.
				m=INDEX(line,deli) ! Position of next ','.
				line=line(m+1:) ! Truncate the line.
		
				! Time array treatment, form is [yyyy mm dd hh mn ss]:
				DO j=1,6	! Loop from year entry to seconds entry
					IF(j.NE.6) THEN ! If not seconds entry.
						IF(j.EQ.1.OR.j.EQ.2) THEN 		! Year or month position.
							sdeli='-'    			! Delimiter
						ELSEIF(j.EQ.3) THEN			! Day position.
							sdeli=' '    			! Delimiter
						ELSEIF(j.EQ.4.OR.j.EQ.5) THEN 		! Hour or minute position.
							sdeli=':'
						END IF
						m=INDEX(temptime,sdeli)
						READ(temptime(1:m-1),*) ttemp(i,j)
						temptime=temptime(m+1:)			! Truncation
					ELSE 						! Seconds position.
						READ(temptime,*) ttemp(i,j)
					END IF
				END DO

				! Reset the string delimiter.
				sdeli='"'
				! Variable array treatment.
				DO j=2,nargs
					IF(j.NE.nargs) THEN
						m=INDEX(line,deli) ! -"- next ','.
						IF(j.EQ.posvar(k+1)) THEN ! IF the position corresponds to that of variable k+1 							! "NaN" string handling.
							n=INDEX(line,sdeli) ! -"- next '"'.
							IF(n.GT.0) THEN ! "NaN" string occurence.
								nanentries=nanentries+1
								o=INDEX(line(n+1:),sdeli)
								temparray(k)=nanreal
								line=line(o+1:) ! Truncate the line.
							ELSE ! No "NaN" string occurence.
								! Ensure that there is no "" occurence.
								o=INDEX(line(1:m-1),sdeli)
								!n=INDEX(line(o:m-1),sdeli)
								IF(o.GT.0) THEN
									temparray(k)=nanreal	! Set to nan.
									line=line(m+1:)		! Truncate the line.
								ELSE 
									READ(line(1:m-1),*) temparray(k) ! Read the real to temp array.
									line=line(m+1:) ! Truncate the line.
								END IF
							END IF
							k=k+1 ! Update the temparray index.
						ELSE
							! Position does not correspond to that of variable number k+1.
							line=line(m+1:) ! Truncate the line.
						END IF
					ELSE ! j.EQ.nargs, end of the line.
						! "NaN" string handling.
						n=INDEX(line,sdeli) ! -"- next '"'.
						IF(n.GT.0) THEN ! "NaN" string occurence.
							o=INDEX(line(n+1:),sdeli)
							temparray(k)=nanreal
						ELSE ! No "NaN" string occurence.
							READ(line,*) temparray(k) ! Read the real to temp array.
						END IF
					END IF			
				END DO
				! Allocate the i-th element in the variable arrays to corresponding entries in the
				! temporary array.
				utemp(i)=temparray(1)
				vtemp(i)=temparray(2)
				wtemp(i)=temparray(3)
				Tstemp(i)=temparray(4)+273.15 ! From Celsius to Kelvin.
				htemp(i)=temparray(5)/1e3     ! From g m^-3 to kg m^-3.
				temparray(:)=0 ! Reset the temporary array.
			END DO

			! Truncate the series.
			n=0 ! Status counter.
			i=0 ! Level counter.
			DO WHILE(n.EQ.0)
				i=i+1
				IF(ttemp(i,5).EQ.0.AND.ttemp(i,6).EQ.0) THEN		! Top of the hour entry.
					blockstart=i+1					
					n=1
				ELSEIF(ttemp(i,5).EQ.30.AND.ttemp(i,6).EQ.0) THEN 	! Half hour entry.
					blockstart=i+1
					n=1
				END IF
			END DO
			
			ecv%nn=dlines-blockstart+1	! Effective series length (from the first hour/half hour occurence).
			
			! Allocate the dimensions of the time series and the measured variables.
			ALLOCATE(ecv%t(ecv%nn,6))
			ALLOCATE(ecv%u(ecv%nn))
			ALLOCATE(ecv%v(ecv%nn))
			ALLOCATE(ecv%w(ecv%nn))
			ALLOCATE(ecv%Ts(ecv%nn))
			ALLOCATE(ecv%h(ecv%nn))
			ALLOCATE(ecv%unode(ecv%nn))
			ALLOCATE(ecv%vnode(ecv%nn))
			ALLOCATE(ecv%wnode(ecv%nn))
			ALLOCATE(ecv%Tsnode(ecv%nn))
			ALLOCATE(ecv%hnode(ecv%nn))

			! Allocate the values of the respective series from the truncated temporary series:
			ecv%t=ttemp(blockstart:dlines,:)
			ecv%u=utemp(blockstart:dlines)
			ecv%v=vtemp(blockstart:dlines)
			ecv%w=wtemp(blockstart:dlines)
			ecv%Ts=Tstemp(blockstart:dlines)
			ecv%h=htemp(blockstart:dlines)
			ecv%unode=utemp(blockstart:dlines)
			ecv%vnode=vtemp(blockstart:dlines)
			ecv%wnode=wtemp(blockstart:dlines)
			ecv%Tsnode=Tstemp(blockstart:dlines)
			ecv%hnode=htemp(blockstart:dlines)

			CALL cpu_time(finish)
			PRINT '("-->Time elapsed in reading data from file =",f5.2," (mins)")',(finish-start)/60
			WRITE(ecv%loglun,*) "-->Time elapsed in reading data from file",(finish-start)/60,"mins"
			PRINT '("-->Data runs from ",I4," . ",I2," . ",I2," , ",I2," : ",I2," : ",f5.2)',&
						NINT(ecv%t(1,1)),NINT(ecv%t(1,2)),NINT(ecv%t(1,3)),&
						NINT(ecv%t(1,4)),NINT(ecv%t(1,5)),ecv%t(1,6)
			WRITE(ecv%loglun,*) "-->Data runs from",ecv%t(1,:)
			PRINT '("   and to ",I4," . ",I2," . ",I2," , ",I2," : ",I2," : ",f5.2)',&
						NINT(ecv%t(ecv%nn,1)),NINT(ecv%t(ecv%nn,2)),NINT(ecv%t(ecv%nn,3)),&
						NINT(ecv%t(ecv%nn,4)),NINT(ecv%t(ecv%nn,5)),ecv%t(ecv%nn,6)	
			WRITE(ecv%loglun,*) "	and to ",ecv%t(ecv%nn,:)
			PRINT '("   @ 20Hz corresponds to",I10," lines of data")',ecv%nn	
			WRITE(ecv%loglun,*)"   @ 20Hz corresponds to",ecv%nn," lines of data"	
			PRINT '("   with N NaNs=",I7)',nanentries
			WRITE(ecv%loglun,*)"   with N NaNs=",nanentries


			PRINT*,'**Done reading from file**'
			WRITE(ecv%loglun,*)"**Done reading from file**"
			PRINT*,'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'
					

			! Deallocation.
			
			! De-allocate the temporary arrays.
			DEALLOCATE(utemp)
			DEALLOCATE(vtemp)
			DEALLOCATE(wtemp)
			DEALLOCATE(Tstemp)
			DEALLOCATE(htemp)
			DEALLOCATE(ttemp)
					
			! Close the file.
			CLOSE(UNIT=ilun)
			RETURN	
		END SUBROUTINE	


		!// Writes the processed and QA/QCed data of interest to a netcdf file, with the same filename as the
		! input file but a '.nc' filetype extension
		SUBROUTINE writetofile(ecv)
			USE netcdf			
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER				:: i,j,k,l,m 			! Counter.
			INTEGER				:: ncid				! File id.

			INTEGER				:: specavail,dravail		! Finally include status integers specifying whether or not
											! D&R and spectra are available for the series, allowing
											! for a control before these are read (or not read if not
											! available) into a visualization program.


			! Combined dimension ids.
			INTEGER				:: tdimids(2)			! Integer array to hold time dimension identifiers.
			INTEGER				:: statdimids(2)		! Integer array to hold structure&autocorr array dim identifiers.
			INTEGER				:: specdimids(2)		! Integer array to hold spectral array dim identifiers.
			INTEGER				:: ogdimids(2)			! Integer array to hold ogive array dim identifiers.
			INTEGER				:: tspwindimids(3)		! Integer array to hold the spectral window time stamps ids.
			INTEGER				:: tstwindimids(3)		! Integer array to hold the statistical window time stamps ids.
			INTEGER				:: tbldimids(3)			! Integer array to hold the block time stamps ids.
			INTEGER				:: snbsid(2),ssnbsid(2)		! Integer array to hold the sub block ids.

			! Dimension ids.
			INTEGER				:: tid,utid			! Size of time dimension identifier, and unit of time identifier.
			INTEGER				:: statsid,statnsid		! Size of dimensions in structure&autocorr arrays identifier.
			INTEGER				:: specsid,specnsid		! Size of dimensions in spectral arrays identifier.
			INTEGER				:: ogsid			! Size of dimensions in ogive arrays identifier.
			INTEGER				:: tbsid			! Size of bounding time dimensions identifier.
			INTEGER				:: nbsid			! Block dimension identifier.
			!INTEGER				:: nsbsid			! Sub block dimension identifier.
			!INTEGER				:: nssbsid			! Sub-sub block dimension identifier.
			INTEGER				:: availid			! Availability dimension identifier.
			
			! Size of the respective dimensions.
			INTEGER				:: ts,uts,stats,statns,specs,specns,ogs,tbs,nbs,avails,nsbpb,nssbpb

			

			! Integers identifying the variables to 'put' in the netcdf file.
			INTEGER				:: vt				! Variable id for time array.
			INTEGER				:: vruu,vrvv,vrww,vrTT,vrqq	! Remaining variable ids.
			INTEGER				:: vduu,vdvv,vdww,vdTT,vdqq,vdrobukhovl,vdrz,vdruh,vtlag
			INTEGER				:: vd3uu,vd3vv,vd3ww,vd3qq,vd3TT
			INTEGER				:: vd4uu,vd4vv,vd4ww,vd4qq,vd4TT
			INTEGER				:: vsu,vsv,vsw,vsT,vsh,vsobukhovl,vsz,vsuh,vsTstar,vsustar
			INTEGER				:: vcouw,vcovw,vcoTw,vcohw
			INTEGER				:: vquuw,vquvw,vquTw,vquqw	
			INTEGER				:: voguw,vogvw,vogTw,voghw
			INTEGER				:: voguu,vogvv,vogww,vogTT,voghh
			INTEGER				:: vfog,vfsp,vtspwin,vtstwin
			INTEGER				:: vtblock
			INTEGER				:: vmwd,vmws,vcr
			INTEGER				:: vubar,vvbar,vwbar,vhbar,vqbar,vTbar
			INTEGER				:: vpbar,vrhobar,vTvbar,vTsbar
			INTEGER				:: vuwcov,vvwcov,vTwcov,vTswcov,vhwcov	
			INTEGER				:: vtvar,vTsvar,vuvar,vvvar,vwvar,vhvar,vqvar
			INTEGER				:: vobukhovl,vobukhovls,vobukhovlss
			INTEGER				:: vtke,vtau,vustar,vQE,vQH,vz	
			INTEGER				:: vfalat,vfamom,vfasen
			INTEGER				:: vCFsep,vCFSND,vCFSND0,vCFWPL
			INTEGER				:: vlagirga
			INTEGER				:: vcsuw,vcsvw,vcshw,vcsTw
			INTEGER				:: vflag
			INTEGER				:: vspecavail,vdravail	
			INTEGER				:: vfstsw,vfsuw,vfshw,vfstat,vfn,vfd,vfv,vfitc
			INTEGER				:: vfskew,vfskeT,vfkuw,vfkuT,vfsk
			INTEGER				:: volnoqc, vznoqc
			INTEGER				:: vffsetsw,vffsehw	
			

			! Integers identifying the 20Hz series for the processed variables, optional whether or not to include
			! these in the ncfile. Good for a visual inspection of the series, but should ultimately be excluded as 
			! these take up unecessary space in the nc-file. 
			INTEGER				:: vu,vv,vw,vTT,vq,vp,vuno,vvno,vwno,vTTno,vqno

			! Prompts and status to queery inclusion of the processed variables. Might want to do this via command line arguments
			! for block length, stat length and spec length as well in a separate dialogue subroutine.
						
			CHARACTER(LEN=80)		:: cdfname			! Filename.
			
			
			
			! Netcdf filename.
			cdfname=ecv%outfilename

			! Set status integers specifying availability of spectra and statistics.
			IF(ecv%spactive.EQV..TRUE.) THEN
				specavail=1
			ELSE
				specavail=0
			END IF
			IF(ecv%dractive.EQV..TRUE.) THEN
				dravail=1
			ELSE
				dravail=0
			END IF
			avails=1	! Size of availability status integers.
			


			PRINT*,'---------------------------------------------------------------------'
			PRINT*,'Writing to file:',cdfname
			PRINT*,'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'Writing to file:',cdfname
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'


			! Define the size of the dimensions.
			ts=SIZE(ecv%t,1)	! Number of points in high frequency series.
			uts=SIZE(ecv%t,2)	! Number of units of time.
			IF(ecv%dractive.EQV..TRUE.) THEN
				stats=SIZE(ecv%d2uu,1)	! Number of time lags in D & R windows.
				statns=SIZE(ecv%d2uu,2)	! Number of active D & R windows.
			END IF
			IF(ecv%spactive.EQV..TRUE.) THEN
				specs=SIZE(ecv%su,1)	! Number of frequencies in spectra.
				specns=SIZE(ecv%su,2)	! Number of active spectral windows.
				ogs=SIZE(ecv%oguw,1)	! Number of frequencies in the Ogives.
				tbs=SIZE(ecv%tspwin,1)	! Number of timestamps per block/window.
			END IF
			tbs=SIZE(ecv%tblock,1)	! Number of timestamps per block/window.
			nbs=SIZE(ecv%ubar,1)	! Number of blocks.
			!nsbpb=SIZE(ecv%uvars,2)	! Number of short sub-blocks per block.
			!nssbpb=SIZE(ecv%uvarss,2) ! Number of shortest sub-blocks per block.
			

			! Create the netCDF file. CLOBBER=overwrite, NF90_64BIT_OFFSET == large file support.
 			CALL CHECK( nf90_create(cdfname, OR(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid) )
			PRINT*,'Created ncfile'

	

  			! Define the dimensions.
  			CALL CHECK( nf90_def_dim(ncid,"t",ts,tid) )
			CALL CHECK( nf90_def_dim(ncid,"ut",uts,utid) )
			IF(ecv%dractive.EQV..TRUE.) THEN
				CALL CHECK( nf90_def_dim(ncid,"lag",stats,statsid) )
				CALL CHECK( nf90_def_dim(ncid,"statwindow",statns,statnsid) )
			END IF
			IF(ecv%spactive.EQV..TRUE.) THEN
				CALL CHECK( nf90_def_dim(ncid,"f",specs,specsid) )
				CALL CHECK( nf90_def_dim(ncid,"specwindow",specns,specnsid) )
				CALL CHECK( nf90_def_dim(ncid,"f_og",ogs,ogsid)	)
			END IF
			CALL CHECK( nf90_def_dim(ncid,"tbounds",tbs,tbsid)	)
			CALL CHECK( nf90_def_dim(ncid,"blocks",nbs,nbsid)	)
			CALL CHECK( nf90_def_dim(ncid,"available",avails,availid)	)
			!CALL CHECK( nf90_def_dim(ncid,"n_subblocks",nsbpb,nsbsid)	)
			!CALL CHECK( nf90_def_dim(ncid,"n_subsubblocks",nssbpb,nssbsid)	)


			PRINT*,'Defined dimensions'
			
			
 			
			! Multiple dimension ids.
			!snbsid=(/ nbsid,nsbsid /)
			!ssnbsid=(/ nbsid,nssbsid /)

			tdimids=(/ tid,utid /)
			IF(ecv%dractive.EQV..TRUE.) THEN
				statdimids=(/ statsid,statnsid /)
				tstwindimids=(/ tbsid,statnsid,utid /)
			END IF
			IF(ecv%spactive.EQV..TRUE.) THEN
				specdimids=(/ specsid,specnsid /)
				ogdimids=(/ ogsid,specnsid /)
				tspwindimids=(/ tbsid,specnsid,utid /)
			END IF
			tbldimids=(/ tbsid,nbsid,utid /)



			! Define the variables. 
			IF(ecv%spactive.EQV..TRUE.) THEN
				CALL CHECK( 	nf90_def_var(ncid,"couw",NF90_DOUBLE,specdimids,vcouw)	)
				CALL CHECK(	nf90_def_var(ncid,"covw",NF90_DOUBLE,specdimids,vcovw)	)
				CALL CHECK(	nf90_def_var(ncid,"coTw",NF90_DOUBLE,specdimids,vcoTw)	)
				CALL CHECK(	nf90_def_var(ncid,"cohw",NF90_DOUBLE,specdimids,vcohw)	)
				CALL CHECK(	nf90_def_var(ncid,"su",NF90_DOUBLE,specdimids,vsu)	)
				CALL CHECK(	nf90_def_var(ncid,"sv",NF90_DOUBLE,specdimids,vsv)	)
				CALL CHECK(	nf90_def_var(ncid,"sw",NF90_DOUBLE,specdimids,vsw)	)
				CALL CHECK(	nf90_def_var(ncid,"sT",NF90_DOUBLE,specdimids,vsT)	)
				CALL CHECK(	nf90_def_var(ncid,"sh",NF90_DOUBLE,specdimids,vsh)	)
				CALL CHECK(	nf90_def_var(ncid,"quuw",NF90_DOUBLE,specdimids,vquuw)	)
				CALL CHECK(	nf90_def_var(ncid,"quvw",NF90_DOUBLE,specdimids,vquvw)	)
				CALL CHECK(	nf90_def_var(ncid,"quTw",NF90_DOUBLE,specdimids,vquTw)	)
				CALL CHECK(	nf90_def_var(ncid,"quqw",NF90_DOUBLE,specdimids,vquqw)	)
				CALL CHECK(	nf90_def_var(ncid,"oguw",NF90_DOUBLE,ogdimids,voguw)	)
				CALL CHECK(	nf90_def_var(ncid,"ogvw",NF90_DOUBLE,ogdimids,vogvw)	)
				CALL CHECK(	nf90_def_var(ncid,"ogTw",NF90_DOUBLE,ogdimids,vogTw)	)
				CALL CHECK(	nf90_def_var(ncid,"oghw",NF90_DOUBLE,ogdimids,voghw)	)
				CALL CHECK(	nf90_def_var(ncid,"oguu",NF90_DOUBLE,ogdimids,voguu)	)
				CALL CHECK(	nf90_def_var(ncid,"ogvv",NF90_DOUBLE,ogdimids,vogvv)	)
				CALL CHECK(	nf90_def_var(ncid,"ogww",NF90_DOUBLE,ogdimids,vogww)	)
				CALL CHECK(	nf90_def_var(ncid,"ogTT",NF90_DOUBLE,ogdimids,vogTT)	)
				CALL CHECK(	nf90_def_var(ncid,"oghh",NF90_DOUBLE,ogdimids,voghh)	)
				CALL CHECK(	nf90_def_var(ncid,"obukhovl_spec",NF90_DOUBLE,specnsid,vsobukhovl)	)
				CALL CHECK(	nf90_def_var(ncid,"fog",NF90_DOUBLE,ogsid,vfog)	)
				CALL CHECK(	nf90_def_var(ncid,"fsp",NF90_DOUBLE,specsid,vfsp)	)
				CALL CHECK(	nf90_def_var(ncid,"twin_spec",NF90_DOUBLE,tspwindimids,vtspwin)	)
				CALL CHECK(	nf90_def_var(ncid,"U_spec",NF90_DOUBLE,specnsid,vsuh)	)
				CALL CHECK(	nf90_def_var(ncid,"ustar_spec",NF90_DOUBLE,specnsid,vsustar) )
				CALL CHECK(	nf90_def_var(ncid,"Tstar_spec",NF90_DOUBLE,specnsid,vsTstar) )
				CALL CHECK(	nf90_def_var(ncid,"z_spec",NF90_DOUBLE,specnsid,vsz)	)
			END IF
			IF(ecv%dractive.EQV..TRUE.) THEN
				CALL CHECK(	nf90_def_var(ncid,"ruu",NF90_DOUBLE,statdimids,vruu)	)
				CALL CHECK(	nf90_def_var(ncid,"rvv",NF90_DOUBLE,statdimids,vrvv)	)
				CALL CHECK(	nf90_def_var(ncid,"rww",NF90_DOUBLE,statdimids,vrww)	)
				CALL CHECK(	nf90_def_var(ncid,"rTT",NF90_DOUBLE,statdimids,vrTT)	)
				CALL CHECK(	nf90_def_var(ncid,"rhh",NF90_DOUBLE,statdimids,vrqq)	)
				CALL CHECK(	nf90_def_var(ncid,"d2uu",NF90_DOUBLE,statdimids,vduu)	)
				CALL CHECK(	nf90_def_var(ncid,"d2vv",NF90_DOUBLE,statdimids,vdvv)	)
				CALL CHECK(	nf90_def_var(ncid,"d2ww",NF90_DOUBLE,statdimids,vdww)	)
				CALL CHECK(	nf90_def_var(ncid,"d2TT",NF90_DOUBLE,statdimids,vdTT)	)
				CALL CHECK(	nf90_def_var(ncid,"d2hh",NF90_DOUBLE,statdimids,vdqq)	)
				CALL CHECK(	nf90_def_var(ncid,"d3uu",NF90_DOUBLE,statdimids,vd3uu)	)
				CALL CHECK(	nf90_def_var(ncid,"d3vv",NF90_DOUBLE,statdimids,vd3vv)	)
				CALL CHECK(	nf90_def_var(ncid,"d3ww",NF90_DOUBLE,statdimids,vd3ww)	)
				CALL CHECK(	nf90_def_var(ncid,"d3TT",NF90_DOUBLE,statdimids,vd3TT)	)
				CALL CHECK(	nf90_def_var(ncid,"d3hh",NF90_DOUBLE,statdimids,vd3qq)	)
				CALL CHECK(	nf90_def_var(ncid,"d4uu",NF90_DOUBLE,statdimids,vd4uu)	)
				CALL CHECK(	nf90_def_var(ncid,"d4vv",NF90_DOUBLE,statdimids,vd4vv)	)
				CALL CHECK(	nf90_def_var(ncid,"d4ww",NF90_DOUBLE,statdimids,vd4ww)	)
				CALL CHECK(	nf90_def_var(ncid,"d4TT",NF90_DOUBLE,statdimids,vd4TT)	)
				CALL CHECK(	nf90_def_var(ncid,"d4hh",NF90_DOUBLE,statdimids,vd4qq)	)
				CALL CHECK(	nf90_def_var(ncid,"obukhovl_dr",NF90_DOUBLE,statnsid,vdrobukhovl)	)
				CALL CHECK( 	nf90_def_var(ncid,"z_dr",NF90_DOUBLE,statnsid,vdrz)	)
				CALL CHECK(	nf90_def_var(ncid,"U_dr",NF90_DOUBLE,statnsid,vdruh)	)
				CALL CHECK(	nf90_def_var(ncid,"twin_dr",NF90_DOUBLE,tstwindimids,vtstwin)	)
				CALL CHECK(	nf90_def_var(ncid,"tlag",NF90_DOUBLE,statsid,vtlag)	)
			END IF
			! Define the processed high frequency sampled series if they are to be included.			
			IF(ecv%hfstat.EQV..TRUE.) THEN
				CALL CHECK(	nf90_def_var(ncid,"t",NF90_DOUBLE,tdimids, vt)	)
				CALL CHECK(	nf90_def_var(ncid,"uhf",NF90_DOUBLE,tid,vu)	)
				CALL CHECK(	nf90_def_var(ncid,"vhf",NF90_DOUBLE,tid,vv)	)
				CALL CHECK(	nf90_def_var(ncid,"whf",NF90_DOUBLE,tid,vw)	)
				CALL CHECK(	nf90_def_var(ncid,"Tshf",NF90_DOUBLE,tid,vTT)	)
				CALL CHECK(	nf90_def_var(ncid,"hhf",NF90_DOUBLE,tid,vq)	)
				CALL CHECK(	nf90_def_var(ncid,"uhfnode",NF90_DOUBLE,tid,vuno)	)
				CALL CHECK(	nf90_def_var(ncid,"vhfnode",NF90_DOUBLE,tid,vvno)	)
				CALL CHECK(	nf90_def_var(ncid,"whfnode",NF90_DOUBLE,tid,vwno)	)
				CALL CHECK(	nf90_def_var(ncid,"Tshfnode",NF90_DOUBLE,tid,vTTno)	)
				CALL CHECK(	nf90_def_var(ncid,"hhfnode",NF90_DOUBLE,tid,vqno)	)
			
			END IF
			
			! Define the block values: averages, fluxes and parameters.
			CALL CHECK(	nf90_def_var(ncid,"wd",NF90_DOUBLE,nbsid,vmwd)		)
			CALL CHECK(	nf90_def_var(ncid,"ws",NF90_DOUBLE,nbsid,vmws)		)
			CALL CHECK(	nf90_def_var(ncid,"CR",NF90_DOUBLE,nbsid,vcr)		)
			CALL CHECK(	nf90_def_var(ncid,"ubar",NF90_DOUBLE,nbsid,vubar)	)
			CALL CHECK(	nf90_def_var(ncid,"vbar",NF90_DOUBLE,nbsid,vvbar)	)
			CALL CHECK(	nf90_def_var(ncid,"wbar",NF90_DOUBLE,nbsid,vwbar)	)
			CALL CHECK(	nf90_def_var(ncid,"hbar",NF90_DOUBLE,nbsid,vhbar)	)
			CALL CHECK(	nf90_def_var(ncid,"qbar",NF90_DOUBLE,nbsid,vqbar)	)
			CALL CHECK(	nf90_def_var(ncid,"Tbar",NF90_DOUBLE,nbsid,vTbar)	)
			CALL CHECK(	nf90_def_var(ncid,"pbar",NF90_DOUBLE,nbsid,vpbar)	)
			CALL CHECK(	nf90_def_var(ncid,"rhobar",NF90_DOUBLE,nbsid,vrhobar)	)
			CALL CHECK(	nf90_def_var(ncid,"Tvbar",NF90_DOUBLE,nbsid,vTvbar)	)
			CALL CHECK(	nf90_def_var(ncid,"Tsbar",NF90_DOUBLE,nbsid,vTsbar)	)
			CALL CHECK(	nf90_def_var(ncid,"tblock",NF90_DOUBLE,tbldimids,vtblock)	)
			CALL CHECK(	nf90_def_var(ncid,"uwcov",NF90_DOUBLE,nbsid,vuwcov)	)
			CALL CHECK(	nf90_def_var(ncid,"vwcov",NF90_DOUBLE,nbsid,vvwcov)	)
			CALL CHECK(	nf90_def_var(ncid,"Twcov",NF90_DOUBLE,nbsid,vTwcov)	)
			CALL CHECK(	nf90_def_var(ncid,"Tswcov",NF90_DOUBLE,nbsid,vTswcov)	)
			CALL CHECK(	nf90_def_var(ncid,"hwcov",NF90_DOUBLE,nbsid,vhwcov)	)
			CALL CHECK( 	nf90_def_var(ncid,"Tvar",NF90_DOUBLE,nbsid,vTvar)	)
			CALL CHECK( 	nf90_def_var(ncid,"Tsvar",NF90_DOUBLE,nbsid,vTsvar)	)
			CALL CHECK(	nf90_def_var(ncid,"uvar",NF90_DOUBLE,nbsid,vuvar)	)
			CALL CHECK(	nf90_def_var(ncid,"vvar",NF90_DOUBLE,nbsid,vvvar)	)
			CALL CHECK(	nf90_def_var(ncid,"wvar",NF90_DOUBLE,nbsid,vwvar)	)
			CALL CHECK(	nf90_def_var(ncid,"hvar",NF90_DOUBLE,nbsid,vhvar)	)
			CALL CHECK(	nf90_def_var(ncid,"qvar",NF90_DOUBLE,nbsid,vqvar)	)	
			CALL CHECK(	nf90_def_var(ncid,"obukhovl",NF90_DOUBLE,nbsid,vobukhovl)	)
			CALL CHECK(	nf90_def_var(ncid,"zblock",NF90_DOUBLE,nbsid,vz)	)
			CALL CHECK(	nf90_def_var(ncid,"tke",NF90_DOUBLE,nbsid,vtke)		)
			CALL CHECK(	nf90_def_var(ncid,"QE",NF90_DOUBLE,nbsid,vqe)		)
			CALL CHECK(	nf90_def_var(ncid,"QH",NF90_DOUBLE,nbsid,vqh)		)
			CALL CHECK(	nf90_def_var(ncid,"tau",NF90_DOUBLE,nbsid,vtau)		)
			CALL CHECK(	nf90_def_var(ncid,"ustar",NF90_DOUBLE,nbsid,vustar)	)
			CALL CHECK(	nf90_def_var(ncid,"blockflag",NF90_INT,nbsid,vflag)	)
			CALL CHECK(	nf90_def_var(ncid,"rns_uw",NF90_DOUBLE,nbsid,vcsuw)	)
			CALL CHECK(	nf90_def_var(ncid,"rns_vw",NF90_DOUBLE,nbsid,vcsvw)	)
			CALL CHECK(	nf90_def_var(ncid,"rns_Tw",NF90_DOUBLE,nbsid,vcsTw)	)
			CALL CHECK(	nf90_def_var(ncid,"rns_qw",NF90_DOUBLE,nbsid,vcshw)	)
			CALL CHECK(	nf90_def_var(ncid,"Fa_mom",NF90_DOUBLE,nbsid,vfamom)	)
			CALL CHECK(	nf90_def_var(ncid,"Fa_lat",NF90_DOUBLE,nbsid,vfalat)	)
			CALL CHECK(	nf90_def_var(ncid,"Fa_sen",NF90_DOUBLE,nbsid,vfasen)	)
			CALL CHECK(	nf90_def_var(ncid,"CF_SND",NF90_DOUBLE,nbsid,vCFSND)	)
			CALL CHECK(	nf90_def_var(ncid,"CF_SND0",NF90_DOUBLE,nbsid,vCFSND0)	)
			CALL CHECK(	nf90_def_var(ncid,"CF_WPL",NF90_DOUBLE,nbsid,vCFWPL)	)
			CALL CHECK(	nf90_def_var(ncid,"CF_sep",NF90_DOUBLE,nbsid,vCFsep)	)
			CALL CHECK(	nf90_def_var(ncid,"IRGA_lagcount",NF90_INT,nbsid,vlagirga)	)
			CALL CHECK(	nf90_def_var(ncid,"spectra_status",NF90_INT,availid,vspecavail)	)
			CALL CHECK(	nf90_def_var(ncid,"autostruct_status",NF90_INT,availid,vdravail	)	)
			CALL CHECK(	nf90_def_var(ncid,"flagstatuw",NF90_INT,nbsid,vfsuw)	)
			CALL CHECK(	nf90_def_var(ncid,"flagstatTw",NF90_INT,nbsid,vfstsw)	)
			CALL CHECK(	nf90_def_var(ncid,"flagstathw",NF90_INT,nbsid,vfshw)	)
			CALL CHECK(	nf90_def_var(ncid,"flagstat",NF90_INT,nbsid,vfstat)	)
			CALL CHECK(	nf90_def_var(ncid,"flagdistor",NF90_INT,nbsid,vfd)	)
			CALL CHECK(	nf90_def_var(ncid,"flagfaulty",NF90_INT,nbsid,vfn)	)
			CALL CHECK(	nf90_def_var(ncid,"flagvertical",NF90_INT,nbsid,vfv)	)
			CALL CHECK(	nf90_def_var(ncid,"flagitcw",NF90_INT,nbsid,vfitc)	)
			CALL CHECK(	nf90_def_var(ncid,"flagskeww",NF90_INT,nbsid,vfskew)	)
			CALL CHECK(	nf90_def_var(ncid,"flagskewT",NF90_INT,nbsid,vfskeT)	)
			CALL CHECK(	nf90_def_var(ncid,"flagkurtw",NF90_INT,nbsid,vfkuw)	)
			CALL CHECK(	nf90_def_var(ncid,"flagkurtT",NF90_INT,nbsid,vfkuT)	)
			CALL CHECK(	nf90_def_var(ncid,"flagskewkurt",NF90_INT,nbsid,vfsk)	)
			!CALL CHECK(	nf90_def_var(ncid,"z_noqc",NF90_DOUBLE,nbsid,vznoqc)	)
			!CALL CHECK(	nf90_def_var(ncid,"obukhovl_noqc",NF90_DOUBLE,nbsid,volnoqc)	)
			CALL CHECK(	nf90_def_var(ncid,"ffseQH",NF90_DOUBLE,nbsid,vffsetsw)	)
			CALL CHECK(	nf90_def_var(ncid,"ffseQE",NF90_DOUBLE,nbsid,vffsehw)	)

			PRINT*,'Defined variables'
			


			! Define attributes.
			IF(ecv%spactive.EQV..TRUE.) THEN
				CALL CHECK(	nf90_put_att(ncid,vcouw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vcovw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vcoTw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vcohw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsh,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vquuw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vquvw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vquTw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vquqw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,voguw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vogvw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vogTw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,voghw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,voguu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vogvv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vogww,'missing_value',ecv%nanreal)	)	
				CALL CHECK(	nf90_put_att(ncid,vogTT,'missing_value',ecv%nanreal)	)	
				CALL CHECK(	nf90_put_att(ncid,voghh,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsobukhovl,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsuh,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsustar,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsTstar,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vsz,'missing_value',ecv%nanreal)	)
			END IF
			IF(ecv%dractive.EQV..TRUE.) THEN
				CALL CHECK(	nf90_put_att(ncid,vruu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vrvv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vrww,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vrTT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vrqq,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vduu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdvv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdww,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdTT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd3uu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd3vv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd3ww,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd3TT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd3qq,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd4uu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd4vv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd4ww,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd4TT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vd4qq,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdqq,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdrobukhovl,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdrz,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vdruh,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vtlag,'missing_value',ecv%nanreal)	)
			END IF
			! Define the attributes of the processed high frequency sampled series if they are to be included.			
			IF(ecv%hfstat.EQV..TRUE.) THEN
				CALL CHECK(	nf90_put_att(ncid,vt,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vu,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vv,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vw,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vTT,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vq,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vuno,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vvno,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vwno,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vTTno,'missing_value',ecv%nanreal)	)
				CALL CHECK(	nf90_put_att(ncid,vqno,'missing_value',ecv%nanreal)	)
			END IF
			! Define the arributes of the block values.
			CALL CHECK(	nf90_put_att(ncid,vmws,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vmwd,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vcr,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vubar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vvbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vwbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vhbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vqbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vpbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vrhobar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTvbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTsbar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vtblock,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vuwcov,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vvwcov,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTwcov,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTswcov,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vhwcov,'missing_value',ecv%nanreal)	)
			CALL CHECK( 	nf90_put_att(ncid,vTvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vTsvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vuvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vvvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vwvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vhvar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vqvar,'missing_value',ecv%nanreal)	)	
			CALL CHECK(	nf90_put_att(ncid,vobukhovl,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vtke,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vqe,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vqh,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vtau,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vustar,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vcsuw,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vcsvw,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vcsTw,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vcshw,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vflag,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vz,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vfamom,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vfalat,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vfasen,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vCFSND,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vCFSND0,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vCFWPL,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vCFsep,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vlagirga,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK( 	nf90_put_att(ncid,vfsk,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfkuT,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfkuw,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfskeT,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfskew,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfitc,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfv,'missing_value',NINT(ecv%nanreal)	) 	)
			CALL CHECK(	nf90_put_att(ncid,vfn,'missing_value',NINT(ecv%nanreal)	)	)
			CALL CHECK(	nf90_put_att(ncid,vfd,'missing_value',NINT(ecv%nanreal)	)	)
			CALL CHECK(	nf90_put_att(ncid,vfstat,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfshw,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfstsw,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vfsuw,'missing_value',NINT(ecv%nanreal) )	)
			CALL CHECK(	nf90_put_att(ncid,vffsetsw,'missing_value',ecv%nanreal)	)
			CALL CHECK(	nf90_put_att(ncid,vffsehw,'missing_value',ecv%nanreal)	)

			PRINT*,'Defined attributes'

			! End define mode.
			CALL CHECK( nf90_enddef(ncid) )  			
			
			! Write ('put') the variables to the file.
			IF(ecv%spactive.EQV..TRUE.) THEN
				CALL CHECK(	nf90_put_var(ncid,vcouw,ecv%couw)	)
				CALL CHECK(	nf90_put_var(ncid,vcovw,ecv%covw)	)
				CALL CHECK(	nf90_put_var(ncid,vcoTw,ecv%coTw)	)
				CALL CHECK(	nf90_put_var(ncid,vcohw,ecv%cohw)	)
				CALL CHECK(	nf90_put_var(ncid,vsu,ecv%su)	)
				CALL CHECK(	nf90_put_var(ncid,vsv,ecv%sv)	)
				CALL CHECK(	nf90_put_var(ncid,vsw,ecv%sw)	)
				CALL CHECK(	nf90_put_var(ncid,vsT,ecv%sT)	)
				CALL CHECK(	nf90_put_var(ncid,vsh,ecv%sh)	)
				CALL CHECK(	nf90_put_var(ncid,vsobukhovl,ecv%sobukhovl)	)
				CALL CHECK(	nf90_put_var(ncid,vsuh,ecv%meanuh)	)
				CALL CHECK(	nf90_put_var(ncid,vquuw,ecv%quuw)	)
				CALL CHECK(	nf90_put_var(ncid,vquvw,ecv%quvw)	)
				CALL CHECK(	nf90_put_var(ncid,vquTw,ecv%quTw)	)
				CALL CHECK(	nf90_put_var(ncid,vquqw,ecv%quqw)	)
				CALL CHECK(	nf90_put_var(ncid,voguw,ecv%oguw)	)
				CALL CHECK(	nf90_put_var(ncid,vogvw,ecv%ogvw)	)
				CALL CHECK(	nf90_put_var(ncid,vogTw,ecv%ogTw)	)			
				CALL CHECK(	nf90_put_var(ncid,voghw,ecv%oghw)	)
				CALL CHECK(	nf90_put_var(ncid,voguu,ecv%oguu)	)
				CALL CHECK(	nf90_put_var(ncid,vogvv,ecv%ogvv)	)
				CALL CHECK(	nf90_put_var(ncid,vogww,ecv%ogww)	)	
				CALL CHECK(	nf90_put_var(ncid,vogTT,ecv%ogTT)	)
				CALL CHECK(	nf90_put_var(ncid,voghh,ecv%oghh)	)
				CALL CHECK(	nf90_put_var(ncid,vfog,ecv%fog)	)
				CALL CHECK(	nf90_put_var(ncid,vfsp,ecv%fsp)	)
				CALL CHECK(	nf90_put_var(ncid,vtspwin,ecv%tspwin)	)
				CALL CHECK(	nf90_put_var(ncid,vsz,ecv%sz)	)
				CALL CHECK(	nf90_put_var(ncid,vsustar,ecv%sustar)	)
				CALL CHECK(	nf90_put_var(ncid,vsTstar,ecv%sTstar)	)
				PRINT*,'ok put sp'	
			END IF	
			IF(ecv%dractive.EQV..TRUE.) THEN		
				CALL CHECK(	nf90_put_var(ncid,vruu,ecv%ruu)	)
				CALL CHECK(	nf90_put_var(ncid,vrvv,ecv%rvv)	)
				CALL CHECK(	nf90_put_var(ncid,vrww,ecv%rww)	)
				CALL CHECK(	nf90_put_var(ncid,vrTT,ecv%rTT)	)
				CALL CHECK(	nf90_put_var(ncid,vrqq,ecv%rqq)	)
				CALL CHECK(	nf90_put_var(ncid,vduu,ecv%d2uu)	)
				CALL CHECK(	nf90_put_var(ncid,vdvv,ecv%d2vv)	)
				CALL CHECK(	nf90_put_var(ncid,vdww,ecv%d2ww)	)
				CALL CHECK(	nf90_put_var(ncid,vdTT,ecv%d2TT)	)
				CALL CHECK(	nf90_put_var(ncid,vdqq,ecv%d2qq)	)
				CALL CHECK(	nf90_put_var(ncid,vd3uu,ecv%d3uu)	)
				CALL CHECK(	nf90_put_var(ncid,vd3vv,ecv%d3vv)	)
				CALL CHECK(	nf90_put_var(ncid,vd3ww,ecv%d3ww)	)
				CALL CHECK(	nf90_put_var(ncid,vd3TT,ecv%d3TT)	)
				CALL CHECK(	nf90_put_var(ncid,vd3qq,ecv%d3qq)	)
				CALL CHECK(	nf90_put_var(ncid,vd4uu,ecv%d4uu)	)
				CALL CHECK(	nf90_put_var(ncid,vd4vv,ecv%d4vv)	)
				CALL CHECK(	nf90_put_var(ncid,vd4ww,ecv%d4ww)	)
				CALL CHECK(	nf90_put_var(ncid,vd4TT,ecv%d4TT)	)
				CALL CHECK(	nf90_put_var(ncid,vd4qq,ecv%d4qq)	)
				CALL CHECK(	nf90_put_var(ncid,vdrobukhovl,ecv%drobukhovl)	)
				CALL CHECK(	nf90_put_var(ncid,vtstwin,ecv%tstwin)	)
				CALL CHECK(	nf90_put_var(ncid,vdrz,ecv%drz)	)
				CALL CHECK(	nf90_put_var(ncid,vdruh,ecv%druh)	)
				CALL CHECK(	nf90_put_var(ncid,vtlag,ecv%tlag)	)
				PRINT*,'ok put dr'
			END IF
			! Put the processed high frequency sampled series if they are to be included.			
			IF(ecv%hfstat.EQV..TRUE.) THEN
				CALL CHECK(	nf90_put_var(ncid,vt,ecv%t)	)
				CALL CHECK(	nf90_put_var(ncid,vu,ecv%u)	)
				CALL CHECK(	nf90_put_var(ncid,vv,ecv%v)	)
				CALL CHECK(	nf90_put_var(ncid,vw,ecv%w)	)
				CALL CHECK(	nf90_put_var(ncid,vTT,ecv%Ts)	)	
				CALL CHECK(	nf90_put_var(ncid,vq,ecv%h)	)
				CALL CHECK(	nf90_put_var(ncid,vuno,ecv%unode)	)
				CALL CHECK(	nf90_put_var(ncid,vvno,ecv%vnode)	)
				CALL CHECK(	nf90_put_var(ncid,vwno,ecv%wnode)	)
				CALL CHECK(	nf90_put_var(ncid,vTTno,ecv%Tsnode)	)	
				CALL CHECK(	nf90_put_var(ncid,vqno,ecv%hnode)	)	
			END IF
			! Put the block timestamp.
			CALL CHECK(	nf90_put_var(ncid,vtblock,ecv%tblock)	)
			! Put the block values.
			CALL CHECK(	nf90_put_var(ncid,vmws,ecv%mws)		)
			CALL CHECK(	nf90_put_var(ncid,vmwd,ecv%mwd)		)
			CALL CHECK(	nf90_put_var(ncid,vcr,ecv%cr)		)
			CALL CHECK(	nf90_put_var(ncid,vubar,ecv%ubar)	)
			CALL CHECK(	nf90_put_var(ncid,vvbar,ecv%vbar)	)
			CALL CHECK(	nf90_put_var(ncid,vwbar,ecv%wbar)	)
			CALL CHECK(	nf90_put_var(ncid,vhbar,ecv%hbar)	)
			CALL CHECK(	nf90_put_var(ncid,vqbar,ecv%qbar)	)
			CALL CHECK(	nf90_put_var(ncid,vTbar,ecv%Tbar)	)
			CALL CHECK(	nf90_put_var(ncid,vpbar,ecv%pbar)	)
			CALL CHECK(	nf90_put_var(ncid,vrhobar,ecv%rhobar)	)
			CALL CHECK(	nf90_put_var(ncid,vTvbar,ecv%Tvbar)	)
			CALL CHECK(	nf90_put_var(ncid,vTsbar,ecv%Tsbar)	)
			CALL CHECK(	nf90_put_var(ncid,vuwcov,ecv%uwcov)	)
			CALL CHECK(	nf90_put_var(ncid,vvwcov,ecv%vwcov)	)
			CALL CHECK(	nf90_put_var(ncid,vTwcov,ecv%Twcov)	)
			CALL CHECK(	nf90_put_var(ncid,vTswcov,ecv%Tswcov)	)
			CALL CHECK(	nf90_put_var(ncid,vhwcov,ecv%hwcov)	)
			CALL CHECK( 	nf90_put_var(ncid,vTvar,ecv%Tvar)	)
			CALL CHECK( 	nf90_put_var(ncid,vTsvar,ecv%Tsvar)	)
			CALL CHECK(	nf90_put_var(ncid,vuvar,ecv%uvar)	)
			CALL CHECK(	nf90_put_var(ncid,vvvar,ecv%vvar)	)
			CALL CHECK(	nf90_put_var(ncid,vwvar,ecv%wvar)	)
			CALL CHECK(	nf90_put_var(ncid,vhvar,ecv%hvar)	)
			CALL CHECK(	nf90_put_var(ncid,vqvar,ecv%qvar)	)
			CALL CHECK(	nf90_put_var(ncid,vobukhovl,ecv%obukhovl)	)
			CALL CHECK(	nf90_put_vaR(ncid,vz,ecv%zbar)	)
			CALL CHECK(	nf90_put_var(ncid,vtke,ecv%tke)	)
			CALL CHECK(	nf90_put_var(ncid,vqe,ecv%qe)	)
			CALL CHECK(	nf90_put_var(ncid,vqh,ecv%qh)	)
			CALL CHECK(	nf90_put_var(ncid,vtau,ecv%tau)	)
			CALL CHECK(	nf90_put_var(ncid,vustar,ecv%ustar)	)
			CALL CHECK(	nf90_put_var(ncid,vcsuw,ecv%csuw)	)
			CALL CHECK(	nf90_put_var(ncid,vcsTw,ecv%csTw)	)
			CALL CHECK(	nf90_put_var(ncid,vcshw,ecv%cshw)	)
			CALL CHECK(	nf90_put_var(ncid,vflag,ecv%flag)	)
			CALL CHECK(	nf90_put_var(ncid,vspecavail,specavail)	)
			CALL CHECK(	nf90_put_var(ncid,vdravail,dravail)	)
			CALL CHECK(	nf90_put_var(ncid,vfamom,ecv%Fa_mom)	)
			CALL CHECK(	nf90_put_var(ncid,vfalat,ecv%Fa_lat)	)
			CALL CHECK(	nf90_put_var(ncid,vfasen,ecv%Fa_sen)	)
			CALL CHECK(	nf90_put_var(ncid,vCFSND,ecv%CFSND)	)
			CALL CHECK(	nf90_put_var(ncid,vCFSND0,ecv%CFSND0)	)
			CALL CHECK(	nf90_put_var(ncid,vCFWPL,ecv%CFWPL)	)
			CALL CHECK(	nf90_put_var(ncid,vCFsep,ecv%CFsep)	)
			CALL CHECK(	nf90_put_var(ncid,vlagirga,ecv%lagIRGA)	)
			CALL CHECK(	nf90_put_var(ncid,vfsuw,ecv%flagstatuw)	)
			CALL CHECK(	nf90_put_var(ncid,vfstsw,ecv%flagstatTsw))
			CALL CHECK(	nf90_put_var(ncid,vfshw,ecv%flagstathw)	)
			CALL CHECK(	nf90_put_var(ncid,vfstat,ecv%flagstat)	)
			CALL CHECK(	nf90_put_var(ncid,vfd,ecv%flagdir)	)
			CALL CHECK(	nf90_put_var(ncid,vfn,ecv%flagnan)	)
			CALL CHECK(	nf90_put_var(ncid,vfv,ecv%flagvert)	)
			CALL CHECK(	nf90_put_var(ncid,vfitc,ecv%flagitcw)	)
			CALL CHECK(	nf90_put_var(ncid,vfsk,ecv%flagsk)	)
			CALL CHECK(	nf90_put_var(ncid,vfskew,ecv%flagskeww) )
			CALL CHECK(	nf90_put_var(ncid,vfskeT,ecv%flagskewT)	)
			CALL CHECK(	nf90_put_var(ncid,vfkuT,ecv%flagkurtT)	)
			CALL CHECK(	nf90_put_var(ncid,vfkuw,ecv%flagkurtw)	)
			!CALL CHECK(	nf90_put_var(ncid,volnoqc,ecv%obukhovlnoqc)	)
			!CALL CHECK(	nf90_put_var(ncid,vznoqc,ecv%znoqc)	)
			CALL CHECK(	nf90_put_var(ncid,vffsetsw,ecv%ffsetsw)	)
			CALL CHECK(	nf90_put_var(ncid,vffsehw,ecv%ffsehw)	)

			PRINT*,'Put the variables'



			! Close the file, freeing all resources.
 			CALL check( nf90_close(ncid) )


		END SUBROUTINE




		SUBROUTINE check(IOSTAT) ! Check netcdf file status.
			USE netcdf
			IMPLICIT NONE
			INTEGER, INTENT(IN)	:: IOSTAT

			IF(IOSTAT.NE.nf90_noerr) THEN
				STOP "Write error, stopped"
			END IF
		END SUBROUTINE



		!----------------------------------------------------------------------------------------------------------------------
		! Cleaning, conversion and rotation.
		!----------------------------------------------------------------------------------------------------------------------


		! Automated routine for removing unphysical values from the data determined by a set of limits on each of the
		! measured variables. These limits are somewhat arbitrary, but are set so as to be relevant to the measurement
		! location and period. Such a limiting routine may be considered as a first step in a despiking algorithm,
		! and should thus either be called internally in any such algorithm or called first externally.
		SUBROUTINE confinelims(ecv)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv
			REAL*8					:: ul,vl,wl		! Wind limits.
			REAL*8,DIMENSION(12)			:: Tmll,Tmul		! Monthly temperature limits.
			REAL*8,DIMENSION(12)			:: hmul			! Monthly upper absolute humidity limit.
			REAL*8					:: hll			! Lower absolute humidity limit.
			REAL*8					:: setnan		! Flag.
			INTEGER					:: i,j			! Counter.
			INTEGER					:: uf,vf,wf,Tf,hf	! Flag counters.
			INTEGER					:: huf,hlf		! Flag counters upper lower limit.
			INTEGER					:: us,vs,ws,Ts,hs	! Flag status.

			! Temperature limits file declarations.
			CHARACTER*80	:: line
			INTEGER, PARAMETER 	:: tlun=11
			INTEGER			:: res

			

			!---------T limits file---------------------------
			! Open the b vector file.
 			OPEN(UNIT=tlun,FILE=ecv%tfile,FORM='FORMATTED',&
				IOSTAT=res)
 			IF(res /= 0) THEN
				PRINT *, 'Error in opening T file, status:', res
  				STOP
 			END IF
			! Read the dimension of b from the first line. 
		
			! Read the entries line by line.
			DO i=1,12
				READ(tlun,FMT='(A)') line
				j=INDEX(line,' ')	! Delimiter.
				READ(line(j+1:),*) Tmul(i)
				READ(line(1:j-1),*) Tmll(i)
			END DO

			! Convert to Kelvin from Celsius.
			Tmul=Tmul+273.16
			Tmll=Tmll+273.16

			! Set the absolute humidity limit  on saturation
			! at upper limit of temperature using Teten's equation
			! and a=0.21167*e/T from Foken (2008). Actual slightly larger than
			! the corresponding saturation absolute humidity since T_s ~= T_v >= T.
			DO i=1,12
				hmul(i)=(0.21167/Tmul(i))*(6.11*EXP(17.6294*(Tmul(i)-273.16)/(Tmul(i)-35.86)))
			END DO
			!hul=1e2!(0.21167/Tul)*(6.11*EXP(17.6294*(Tul-273.16)/(Tul-35.86)))	! About 17g/m^

			hll=0

			CLOSE(UNIT=tlun)




			! Print start of call.
			PRINT*,'***************************************'
			PRINT*,'Commencing physical plausability test confinelims()'
			WRITE(ecv%loglun,*)'***************************************'
			WRITE(ecv%loglun,*)'Commencing physical plausability test confinelims()'

			! Set the outer limit flag.
			setnan=ecv%nanreal

			! Set the absolute wind limits (ms^-1).
			ul=32.768!40
			vl=32.768!40
			wl=5
		
			

	
			! Ensure all the flag status and flag counter integers are initially zero.
			Ts=0
			us=0
			vs=0
			ws=0
			hs=0
			Tf=0
			uf=0
			vf=0
			wf=0
			hf=0
			hlf=0
			huf=0
			! Scan the variables and flag values exceeding limits.
			DO i=1,ecv%nn
				IF(ecv%u(i).EQ.setnan) THEN
					us=1
				ELSE IF(ABS(ecv%u(i)).GE.ul) THEN
					ecv%u(i)=setnan
					uf=uf+1
					us=1
				END IF
				IF(ecv%v(i).EQ.setnan) THEN
					vs=1
				ELSE IF(ABS(ecv%v(i)).GE.vl) THEN
						ecv%v(i)=setnan
						vf=vf+1
						vs=1
				END IF
				IF(ecv%w(i).EQ.setnan) THEN
					ws=1
				ELSE IF(ABS(ecv%w(i)).GE.wl) THEN
					ecv%w(i)=setnan
					wf=wf+1
					ws=1
				END IF
				IF(ecv%Ts(i).EQ.setnan) THEN
					Ts=1
				ELSE IF( ecv%Ts(i).GE.Tmul(NINT(ecv%t(i,2)))&
					.OR.ecv%Ts(i).LE.Tmll(NINT(ecv%t(i,2))) ) THEN
					ecv%Ts(i)=setnan
					Tf=Tf+1
					Ts=1
				END IF
				! Based on a visual inspection, the LI-7500 seems to report
				! erroneous values, both negative and above saturation at around
				! the same time. So we exclude these periods from the analysis.
				IF(ecv%h(i).EQ.setnan) THEN
					hs=1
				ELSEIF(ecv%h(i).LT.hll) THEN!.AND.ecv%confined.EQV..TRUE.) THEN	! <-- Only set h limits in second pass.
					ecv%h(i)=ecv%nanreal
					hlf=hlf+1
					hf=hf+1
					hs=1			
				ELSEIF(ecv%h(i).GT.hmul(NINT(ecv%t(i,2)))) THEN!.AND.ecv%confined.EQV..TRUE.) THEN	! <-- -"-.
					ecv%h(i)=ecv%nanreal
					huf=huf+1
					hf=hf+1
					hs=1
				END IF
				! Finally flag ALL variables if one or more variables is flagged.
				IF(us.EQ.1.OR.vs.EQ.1.OR.ws.EQ.1.OR.Ts.EQ.1.OR.hs.EQ.1) THEN
					ecv%u(i)=setnan
					ecv%v(i)=setnan
					ecv%w(i)=setnan
					ecv%Ts(i)=setnan
					ecv%h(i)=setnan
				END IF
				! Reset status flags.
				us=0
				vs=0
				ws=0
				Ts=0
				hs=0
			END DO

			! Set confined status to true if not arleady set.
			!IF(ecv%confined.EQV..FALSE.) THEN
			!	ecv%confined=.TRUE.
			!END IF


			PRINT '("-->ulim flags=",I10)',uf
			PRINT '("-->vlim flags=",I10)',vf
			PRINT '("-->wlim flags=",I10)',wf
			PRINT '("-->Tlim flags=",I10)',Tf
			PRINT '("-->hlim flags=",I10," lower limit %:",I10)',hf,NINT(100*(1.0*hlf)/(1.0*hlf+huf))
			
			WRITE(ecv%loglun,*)"-->ulim flags=",uf
			WRITE(ecv%loglun,*)"-->vlim flags=",vf
			WRITE(ecv%loglun,*)"-->wlim flags=",wf
			WRITE(ecv%loglun,*)"-->Tlim flags=",Tf
			WRITE(ecv%loglun,*)"-->hlim flags=",hf," lower limit %:",NINT(100*(1.0*hlf)/(1.0*hlf+huf))

			PRINT*,'**end of confinelims()**'
			PRINT*,'***************************************'

			WRITE(ecv%loglun,*)'**end of confinelims()**'
			WRITE(ecv%loglun,*)'***************************************'

			RETURN
		END SUBROUTINE
		
		! The median absolute deviation despiking algorithm of Mauder et al. (2013), is less strict and
		! faster (one pass and block windows) than the running sigma deviation despiking algorithm of Højstrup (1993). 
		! Furthermore the algorithm is less sensitive to the choice of window size (verified). For any window j
		! with measured variables X=[x_1,...,x_i,...,x_w] where w is the window length the mean absolute deviation
		! is computed based on:
		! MAD=< | x_i - <X> |_i >
		! where <> denotes the median, <X> is the median of the window and repeated indices imply summation,
		! such that |x_i-<X>|_i is the absolute deviation from the median for the entire window, and <> of this
		! is consequently the MAD. A data point in a given window, x_i, is considered not to be a spike if
		! <X>-q*MAD/c <= x_i <= <X>+q*MAD/c
		! Where q is a threshold value which through trial and error restricts the spiking to unphysical points
		! this is a key point as the sigma despiking algorithm is less discriminating and thus risks eliminating
		! unusual yet physical behaviour (which is often of particular interest). 
		! The factor c=0.6745 relates a MAD of unity to a single standard deviation in a Gaussian distribution. 
		! Since our distribution is unlikely to be Gaussian and the MAD is not unity the high q factor means 
		! that only values within extreme wings (where these exist) are eliminated. 
		! Moreover points which are already flagged by the limits or the instrument intself are NOT weighted
		! in the MAD.
		
		SUBROUTINE MADspike(ecv)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv	
			INTEGER					:: w			! Window length=block length.
			INTEGER					:: wl			! Effective window length.
			REAL*8					:: medw			! Median of a given window.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: adw			! Absolute difference of window value to window median.
			REAL*8					:: madw			! MAD of a given window.
			REAL*8					:: ltol,utol		! MAD tolerance bounds in a given window.
			REAL*8,	PARAMETER			:: g=7.0		! The MAD threshold value. If 7 only 5/121 4h windows (even if 0.05%flag).
			REAL*8, PARAMETER			:: c=0.6745		! Gaussian weight.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: var,wvar		! Variable array & windowed variable array.
			INTEGER					:: nvars		! Number of variables to control.
			INTEGER					:: npass		! Number of passes.
			INTEGER					:: i,j,k,n		! Counters.
			INTEGER					:: ll,ul,nn,nw,nwp	! Bounds on indices, series length and window length.
			INTEGER					:: intsave		! Saves an integer value.	
			INTEGER					:: spikecount		! Counts the number of points flagged as spikes.
			REAL*8					:: zero			! Zero reference for h in the file.
		
			! Set the dimensions.
			w=ecv%blockl
			ALLOCATE(adw(w))

			PRINT*,'***************************************'
			PRINT*,'Starting despiking algorithm MADspike()'
			WRITE(ecv%loglun,*)'***************************************'
			WRITE(ecv%loglun,*)'Starting despiking algorithm MADspike()'

			! Set the physical limits using the confining routine and set nans where necessary.
			!PRINT*,'First limit confinement:'
			!CALL	confinelims(ecv) 	


			! Initialize the bounds and medians to be zero.
			ltol=0
			utol=0
			medw=0
			adw(:)=0
			madw=0
		
			! Initialize array dimensions for allocation.
			nn=ecv%nn			! Series length.
			nvars=ecv%nvars-1 		! Number of variables, exclude time.
			wl=w				! Effective window length (discounting nans).
				
			! Allocate variable array dimensions.
			ALLOCATE(var(ecv%nn,nvars))	! Allocate space for a 2d variable array.
			ALLOCATE(wvar(wl,nvars))	! Allocate space a 2d windowed variable array.

			
			! Allocate the variables in the variable array. The third dimension temporarily flags
			! 'nan' by setting the value equal to nanreal, then checks neighbouring points.
			var(:,1)=ecv%u
			var(:,2)=ecv%v
			var(:,3)=ecv%w
			var(:,4)=ecv%Ts
			var(:,5)=ecv%h
			
			! Initialize the windowed variable array to be an array of zeros.
			wvar(:,:)=0
			
			
			! Integer number of windows in the time series.
			nw=FLOOR(1.0*nn/w)

			! Check if nw windows fit the entire series exactly.
			IF(MOD(nn,nw).EQ.0) THEN ! Yes, only loop over nw.
				nwp=nw
			ELSE ! No, add an extra window.
				nwp=nw+1 ! Add an extra window which includes the end of the series.
			END IF
			
			DO n=1,nwp ! Loop over windows.
				! Define the window.
				! If nw doesnt fit the series exactly AND we are at the last window then
				! redefine the last window, which will overlap with the second last
				! window. Note that due to the slack tolerance requirements the despiking 
				! is fairly insensitive to the arbitrary allocation of subsequent windows.
				IF(nw.NE.nwp.AND.n.EQ.nwp) THEN	
					ll=(nn-wl)+1
					ul=nn
				! Otherwise use block increments starting at the begining of the series.
				ELSE			
					ll=(n-1)*wl+1	! Start of window (ll=1 if window=1).
					ul=n*wl		! End of window   (ul=wl if window=1).
				END IF	
				! Allocate the windowed variables.
				wvar=var(ll:ul,:)
				DO j=1,nvars	! Loop over variables.
					! Check that the entire block isn't flagged as 'NaN'.
					IF((ANY(wvar.NE.ecv%nanreal)).EQV..TRUE.) THEN
						! Compute the median value for the given variable window (medw).
						CALL median(ecv,wl,wvar(:,j),medw)
					
						! Compute the absolute deviation for the given variable window.
						DO i=1,wl
							IF(wvar(i,j).NE.ecv%nanreal) THEN ! If not nan compute absdiff.
								adw(i)=ABS(wvar(i,j)-medw)
							ELSE
								adw(i)=ecv%nanreal
							END IF
						END DO
					
						! Compute the MAD for the given variable window.
						CALL median(ecv,wl,adw,madw)
			
						! Reset ad to an array of zeros.
						adw(:)=0
	
						! Compute the corresponding tolerance limits.
						ltol=medw-(g*madw/c)
						utol=medw+(g*madw/c)
					
						!-----------------------------------------------------------------------------------------
						! Despiking.
						!-----------------------------------------------------------------------------------------
						! Loop over the values of the windowed variable.
						! If tolerance is exceeded set to nan (flag as spike),
						! otherwise continue.
						DO i=1,wl 
							! First verify that the variable is not already flagged as a spike
							! by the instrument, the limits subroutine or in an earlier window.
							IF(wvar(i,j).EQ.ecv%nanreal) THEN
								! (Do nothing).
							
							! Check if the tolerance is exceeded if the windowed variable
							! is not already flagged as 'nan'.
							ELSEIF(wvar(i,j).GT.utol.OR.wvar(i,j).LT.ltol) THEN	
								wvar(i,j)=ecv%nanreal
							END IF
							! If none of the above occur, the variable at the given point
							! in time has passed the test, i.e. remains unchanged,
							! and is considered in the remaining analysis.
						END DO
						! Reset values for next variable pass.
						medw=0
						madw=0
					END IF
				END DO	
				! Finally assign the despiked window to the corresponding location in the
				! variable array.
				var(ll:ul,:)=wvar(1:w,:)


				! Reset window values to zero.
				wvar(:,:)=0
				ll=0
				ul=0
			END DO


			! Count the total number of spikes, and set all variables to nan (spike)
			! at a given level if any variables are flagged as a spike. May seem overly
			! strict but due to the 'kind' despiking algorithm typically less than
			! one % of the data is flagged, and if any one variable is a spike it will
			! contaminate neighbouring variables as it is indicative of instrument error
			! due to sensor drift, snow on the sensors etc...
			nn=ecv%nn
			spikecount=0
			DO i=1,nn
				IF(ANY(var(i,:).EQ.ecv%nanreal).EQV..true.) THEN
					spikecount=spikecount+1
					var(i,:)=ecv%nanreal
				END IF
			END DO

			! Verify the TOTAL ratio of spike occurence; if >10% print a warning to terminal
			! and logfile (which should always be checked!). Still allocate the despiked variables.
			IF(((1.0*spikecount)/nn).GT.0.1) THEN
				PRINT*,' --->% of spikes=',(100.0*spikecount)/nn
				PRINT*,' --->NB NB Bad quality data; exclude from analysis'
				WRITE(ecv%loglun,*)' --->% of spikes=',(100.0*spikecount)/nn
				WRITE(ecv%loglun,*)' --->Bad quality data; exclude from analysis'
				ecv%u=var(:,1)
				ecv%v=var(:,2)
				ecv%w=var(:,3)
				ecv%Ts=var(:,4)
				ecv%h=var(:,5)
			ELSE ! If <10% no warning and allocate despiked variables to the processed data.
				ecv%u=var(:,1)
				ecv%v=var(:,2)
				ecv%w=var(:,3)
				ecv%Ts=var(:,4)
				ecv%h=var(:,5)
				PRINT '(" -->Despiked=",f6.2," % of the data in the file")',(100.0*spikecount)/nn
				WRITE(ecv%loglun,*)" -->Despiked=",(100.0*spikecount)/nn," % of the data in the file"
			END IF	

			PRINT*,'Done despiking with MADspike()	       '
			PRINT*,'***************************************'
			WRITE(ecv%loglun,*)'Done despiking with MADspike()	       '
			WRITE(ecv%loglun,*)'***************************************'

			! Deallocate the temporary arrays.
			DEALLOCATE(var)
			DEALLOCATE(wvar)

			RETURN
		END SUBROUTINE

		! Flag blocks of raw data based on the values of the kurtosis and skewness (after a linear detrend) based on the
		! thresholds outlined in Vickers and Mahrt (1997).
		SUBROUTINE skewkurt(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			!  Upper and lower bounds for block skewness and kurtosis hard flags:
			REAL*8				:: skhfbl,skhfbu,kuhfbl,kuhfbu	
			!  Upper and lower bounds for block skewness and kurtosis soft flags:
			REAL*8				:: sksfbl,sksfbu,kusfbl,kusfbu	

			INTEGER				:: i,j,ii,jj,ul,ll	! Counters & bounds.
			REAL*8				:: skw,skT,kuw,kuT	! Temporary reals to hold the kurtosis and skeness values.
			REAL*8,DIMENSION(4)		:: allskflags
			! Track number of flags:
			INTEGER				:: nhfskw,nhfskT,nhfkuw,nhfkuT,nhfsk


			PRINT*,'***************************************'
			PRINT*,'Starting higher order moment analysis skewkurt()'
			WRITE(ecv%loglun,*)'***************************************'
			WRITE(ecv%loglun,*)'Starting higher order moment analysis skewkurt()'

			ii=ecv%blockl			! Number of points in a block.
			jj=FLOOR(1.0*ecv%nn/ii)		! Number of blocks of data.	

			! Allocate the flags:
			ALLOCATE(ecv%flagskewT(jj))
			ALLOCATE(ecv%flagskeww(jj))
			ALLOCATE(ecv%flagkurtT(jj))
			ALLOCATE(ecv%flagkurtw(jj))	
			ALLOCATE(ecv%flagsk(jj))

			! Set the bounds for the flags for kurtosis and sknewness values after VM97:
			skhfbl=-2
			skhfbu=2
			kuhfbl=1
			kuhfbu=8  
			sksfbl=-1
			sksfbu=1
			kusfbl=2
			kusfbu=5

			! Initialize hard flag counters as zeros.
			nhfskw=0
			nhfskT=0
			nhfkuw=0
			nhfkuT=0
			nhfsk=0

			! Loop through the file block by block.
			DO j=1,jj
				! Array bounds for the current block.
				ll=(j-1)*ii+1
				ul=j*ii

				! Make sure temporary values are zero initially.
				kuT=0
				kuw=0
				skT=0
				skw=0

				! Calculate the skewness and kurtosis for the raw (despiked) sonic temperature and
				! vertical velocity in the given block.
				CALL kurtosis(ecv,ii,ecv%Ts(ll:ul),kuT)
				CALL kurtosis(ecv,ii,ecv%w(ll:ul),kuw)
				CALL skewness(ecv,ii,ecv%Ts(ll:ul),skT)
				CALL skewness(ecv,ii,ecv%w(ll:ul),skw)
		
				! Flag accordingly:
				IF(kuT.GT.kusfbu.OR.kuT.LT.kusfbl) THEN
					IF(kuT.GT.kuhfbu.OR.kuT.LT.kuhfbl) THEN
						ecv%flagkurtT(j)=2
						nhfkuT=nhfkuT+1
					ELSE
						ecv%flagkurtT(j)=1
					END IF
				ELSE
					ecv%flagkurtT(j)=0
				END IF
				allskflags(1)=1.0*ecv%flagkurtT(j)

				IF(kuw.GT.kusfbu.OR.kuw.LT.kusfbl) THEN
					IF(kuw.GT.kuhfbu.OR.kuw.LT.kuhfbl) THEN
						ecv%flagkurtw(j)=2
						nhfkuw=nhfkuw+1
					ELSE
						ecv%flagkurtw(j)=1
					END IF
				ELSE
					ecv%flagkurtw(j)=0
				END IF
				allskflags(2)=1.0*ecv%flagkurtw(j)

				IF(skw.GT.sksfbu.OR.skw.LT.sksfbl) THEN
					IF(skw.GT.skhfbu.OR.skw.LT.skhfbl) THEN
						ecv%flagskeww(j)=2
						nhfskw=nhfskw+1
					ELSE
						ecv%flagskeww(j)=1
					END IF
				ELSE
					ecv%flagskeww(j)=0
				END IF
				allskflags(3)=1.0*ecv%flagskeww(j)

				IF(skT.GT.sksfbu.OR.skT.LT.sksfbl) THEN
					IF(skT.GT.skhfbu.OR.skT.LT.skhfbl) THEN
						ecv%flagskewT(j)=2
						nhfskT=nhfskT+1
					ELSE
						ecv%flagskewT(j)=1
					END IF
				ELSE
					ecv%flagskewT(j)=0
				END IF
				allskflags(4)=1.0*ecv%flagskewT(j)

				! Set the combined flag as the maximum of all the kurtosis and skewness flags:
				ecv%flagsk(j)=NINT(MAXVAL(allskflags))
				IF(ecv%flagsk(j).EQ.2) THEN	! Update flag count and set block to just spikes.
					nhfsk=nhfsk+1
					ecv%u(ll:ul)=ecv%nanreal
					ecv%v(ll:ul)=ecv%nanreal
					ecv%w(ll:ul)=ecv%nanreal
					ecv%Ts(ll:ul)=ecv%nanreal
					ecv%h(ll:ul)=ecv%nanreal
				END IF
				
			END DO
		
			! Print results to terminal and logfile.
			PRINT*,'-->% Blocks hard flagged for skewness(T)=',(100.0*nhfskT)/(1.0*jj)
 			PRINT*,'-->% Blocks hard flagged for skewness(w)=',(100.0*nhfskw)/(1.0*jj)
			PRINT*,'-->% Blocks hard flagged for kurtosis(T)=',(100.0*nhfkuT)/(1.0*jj)
			PRINT*,'-->% Blocks hard flagged for kurtosis(w)=',(100.0*nhfkuw)/(1.0*jj)
			PRINT*,'-->% Blocks hard flagged for kurtosis & skewness=',(100.0*nhfsk)/(1.0*jj)
			PRINT*,'eof skewkurt()'
			PRINT*,'***************************************'
			WRITE(ecv%loglun,*)'-->% Blocks hard flagged for skewness(T)=',(100.0*nhfskT)/(1.0*jj)
			WRITE(ecv%loglun,*)'-->% Blocks hard flagged for skewness(w)=',(100.0*nhfskw)/(1.0*jj)
			WRITE(ecv%loglun,*)'-->% Blocks hard flagged for kurtosis(T)=',(100.0*nhfkuT)/(1.0*jj)
			WRITE(ecv%loglun,*)'-->% Blocks hard flagged for kurtosis(w)=',(100.0*nhfkuw)/(1.0*jj)
			WRITE(ecv%loglun,*)'-->% Blocks hard flagged for kurtosis & skewness=',(100.0*nhfsk)/(1.0*jj)
			WRITE(ecv%loglun,*)'eof skewkurt()'
			WRITE(ecv%loglun,*)'***************************************'

			RETURN
		END SUBROUTINE

			




		! Orient: Calculates the block averaged wind direction and wind speed, performed before any block rotation. Assumes that the
		! sonic is propery leveled. The bearing (from due north) that the sonic is pointing INTO must be given by the user as a parameter.
		! The unit vector separating the LI and the CSAT (pointing from the CSAT to the LI) is also calculated in the NESW coordinate system.
		SUBROUTINE	orientnorth(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			REAL*8				:: alpha,phi,theta,dir	! Angles.
			REAL*8				:: pi
			INTEGER				:: ii,jj,ll,ul,i,j,k	! Counters & bounds.
			REAL*8				:: uc,vc		! Mean (horizontal) wind components in NESW coordinate system.
			REAL*8				:: us,vs		! Mean (horizontal) wind components in sonic coordinate system.
			REAL*8				:: calmlim=0.2		! Lowest wind speed for which a wind direction is defined.
			! Temporary reals used to compute the constancy ratio for each block.
			REAL*8			:: usum,vsum,miws,cr		! In order: block sum of u, block sum of v, 
										! sum of inst. wind speed and the constancy ratio. 
			LOGICAL				:: inc			! Whether or not to include block mws and mwd in statistics.
			REAL*8				:: nanlim		! Limit for fraction of block that can be a spike.
			INTEGER				:: nnan			! Number of spikes in the block.
			! Finally the angles needed to compute the unit separation vector in the NESW frame.
			REAL*8				:: wdtoLI	! Wind direction from which the wind is blowing directly from the CSAT
									! to the LI.
			REAL*8				:: wdtoCSAT	! Wind direction from which the wind is blowing directly from the LI
									! to the CSAT.
			REAL*8				:: normbearing	! Bearing of the line that transects the midpoint of the separation
									! vector between the two instruments.

			! Set the nanlimit:
			nanlim=0.1	! (10 %).


			PRINT*,'***************************************'
			PRINT*,'Orienting the sonic to compute horizontal wind directions'
			WRITE(ecv%loglun,*)'***************************************'
			WRITE(ecv%loglun,*)'Orienting the sonic to compute horizontal wind directions'

			ii=ecv%blockl			! Number of points in a block.
			jj=FLOOR(1.0*ecv%nn/ii)		! Number of blocks of data.
			
			! Compute pi.
			pi=4.0*ATAN(1.0)

			! Convert the bearing to radians.
			theta=2*pi*(1.0*ecv%CSATbearing)/360
					
			! Allocate space for the block averaged values.
			ALLOCATE(ecv%mwd(jj))
			ALLOCATE(ecv%mws(jj))

			! Calculate the angle between the due north axis and the sonic x-axis (positive u).
			IF(theta.LT.pi)	THEN
				phi=theta+pi
			ELSE
				phi=theta-pi
			END IF

			! Calculate the angle between due east and the sonic x-axis (positive u).
			IF(phi.LE.pi/2.0) THEN
				alpha=pi/2.0-phi
			ELSE
				alpha=2*pi-(phi-pi/2.0)
			END IF


			! Set the values to zero initially.
			ll=0
			ul=0
			uc=0
			vc=0
			us=0
			vs=0
			DO j=1,jj
				! Array bounds for the current block.
				ll=(j-1)*ii+1
				ul=j*ii

				! Set block inclusion to true initially.
				inc=.TRUE.

				! Calculate the average horizontal wind speeds in the sonic coordinate system
				! for the given block. NaN values are ignored as usual.
				CALL temporalmean(ecv,ii,ecv%u(ll:ul),us)
				CALL temporalmean(ecv,ii,ecv%v(ll:ul),vs)
				
				! Compute the constancy ratio to verify that the wind direction is well defined and test that
				! not more than 10% of the block is a nan value and that the mean wind speed is >=0.2 ms^-1.
				nnan=0
				usum=0
				vsum=0
				miws=0
				! Test if more than 10% of the data in the block is 'nan'.
				DO k=ll,ul
					IF(ecv%u(k).EQ.ecv%nanreal.OR.ecv%v(k).EQ.ecv%nanreal) THEN
						nnan=nnan+1
					ELSE
						usum=usum+ecv%u(k)				! Block sum of u components.
						vsum=vsum+ecv%v(k)				! Block sum of v components.
						miws=miws+SQRT(ecv%u(k)**2+ecv%v(k)**2)		! Sum of instantaneous wind speeds.
					END IF
				END DO
				! If nnan comprises more than 10% of the block then exclude the block from the orientation.
				IF(((1.0*nnan)/ii).GE.nanlim) THEN
					inc=.FALSE.
				END IF
				! Also test the constancy ration; if it is less than 0.7 the wind direction
				! is not well defined.
				
				CR=SQRT(usum**2+vsum**2)/miws
			
				
				IF(inc.EQV..TRUE.) THEN	! <-- If passed nanlim test check that CR>=0.7 AND mws>=0.2 ms^-1.
					IF(CR.LT.0.75.OR.SQRT(us**2+vs**2).LT.calmlim) THEN
						inc=.FALSE.
					END IF
				END IF



				! Calculate the corresponding wind speeds in the NESW coordinate system.
				! Ensure that the entire block is not missing/flagged as nan and that the direction is
				! well defined.
				IF(us.NE.ecv%nanreal.AND.vs.NE.ecv%nanreal.AND.inc.EQV..TRUE.) THEN
					uc=us*COS(alpha)-vs*SIN(alpha)
					vc=us*SIN(alpha)+vs*COS(alpha)
					uc=-1*uc	! Reverse the wind vector to point upwind instead of downwind.
					vc=-1*vc
				ELSE
					uc=ecv%nanreal
					vc=ecv%nanreal	
				END IF
				
				! Use the components in the NESW coordinate system to calculate the wind direction,
				! i.e. where the mean wind is blowing from.

				! First use arc tangent to calculate where the wind is coming from.
				IF(uc.NE.ecv%nanreal) THEN
					dir=ATAN2(uc,vc)
				ELSE
					dir=ecv%nanreal
				END IF

				! Take care of negative values so that the wind direction is output as a bearing.
				IF(dir.GE.0.AND.dir.LE.pi.AND.dir.NE.ecv%nanreal) THEN
					dir=dir
				ELSEIF(dir.NE.ecv%nanreal) THEN
					dir=2*pi+dir
				END IF
			
				! Allocate the mean wind direction in degrees.
				IF(dir.NE.ecv%nanreal) THEN
					dir=(dir)*360.0/(2*pi)
				END IF
				IF(dir.GT.360.AND.dir.NE.ecv%nanreal) THEN	! <-- In case of rounding error.
					dir=0
				END IF
				ecv%mwd(j)=dir

				! Finally compute the mean wind speed.
				IF(dir.NE.ecv%nanreal) THEN
					ecv%mws(j)=SQRT(uc**2+vc**2)
				ELSE
					ecv%mws(j)=ecv%nanreal
				END IF

				! Set values to zero for the next pass.
				dir=0
				uc=0
				vc=0
				us=0
				vs=0
			END DO

			! Add the alpha angle to the ecv structure to be used in the planar fit definition of 
			! the streamline plane unit vectors (ivec and jvec), such that the (system frame) horizontal component
			! of ivec is pointing due east and the (system frame) horizontal component of jvec is pointing due north.
			ecv%alpha=alpha 


			! Finally compute the bearing of the line that transects, i.e. passes through the midpoint and is normal to,
			! the separation vector between the instruments. And use this to compute the wind directions where
			! the LI/CSAT upwind of the CSAT/LI.
			! The module assumes that the separation angle between the instruments is at most 90deg since this is the
			! typical configuration.
			IF(ecv%CSATbearing.GT.ecv%LIbearing) THEN
					IF((ecv%CSATbearing-ecv%LIbearing).GT.90) THEN	! LI in first quadrant and CSAT in last.
						normbearing=ecv%CSATbearing+(ecv%LIbearing+(ecv%CSATbearing-360))/(2.0)
						wdtoLI=normbearing-90
						wdtoCSAT=(normbearing+90)-360
						IF(normbearing.GT.360) THEN	! If normbearing in first quadrant correct.
							wdtoLI=normbearing-90
							normbearing=normbearing-360
							wdtoCSAT=normbearing+90
							
						END IF
					ELSE	! Remaining configurations:
						normbearing=ecv%LIbearing+(ecv%CSATbearing-ecv%LIbearing)/(2.0)
						wdtoLI=normbearing+90
						IF(wdtoLI.GT.360) THEN
							wdtoLI=wdtoLI-360
						END IF
						wdtoCSAT=normbearing-90
						IF(wdtoCSAT.LT.0) THEN
							 wdtoCSAT=360+wdtoCSAT
						END IF
					END IF
			ELSE
					IF((ecv%LIbearing-ecv%CSATbearing).GT.90) THEN	! CSAT in first quadrant and LI in last.
						normbearing=ecv%LIbearing+(ecv%CSATbearing+(ecv%LIbearing-360))/(2.0)
						wdtoLI=(normbearing+90)-360
						wdtoCSAT=normbearing-90
						IF(normbearing.GT.360) THEN	! If normbearing in first quadrant correct.
							wdtoCSAT=normbearing-90
							normbearing=normbearing-360
							wdtoLI=normbearing+90
						END IF
					ELSE	! Remaining configurations:
						normbearing=ecv%CSATbearing+(ecv%LIbearing-ecv%CSATbearing)/(2.0)	
						wdtoCSAT=normbearing+90
						IF(wdtoCSAT.GT.360) THEN
							wdtoCSAT=wdtoCSAT-360
						END IF
						wdtoLI=normbearing-90
						IF(wdtoLI.LT.0) THEN
							wdtoLI=360+wdtoLI
						END IF
					END IF
			END IF

			! Use the wdtoCSAT bearing to define the separation unit vector pointing from the CSAT3 sonic anemometer
			! to the open path LI7500 IRGA.
			! First convert this angle into the angle from due east in the right handed (counter clockwise) sense.
			IF(wdtoCSAT.LE.90) THEN
				wdtoCSAT=90-wdtoCSAT
			ELSEIF(wdtoCSAT.LE.360.AND.wdtoCSAT.GT.90) THEN
				wdtoCSAT=450-wdtoCSAT
			END IF
			ecv%sepunit(3)=0
			ecv%sepunit(1)=COS(wdtoCSAT*(1.0/360.0)*2*pi)
			ecv%sepunit(2)=SIN(wdtoCSAT*(1.0/360.0)*2*pi)
		
			PRINT*,alpha
			PRINT'("--> Angle between due East and the sonic x-axis is",f12.4," deg")',(alpha/(2*pi))*360
			PRINT'("--> Wind direction where the wind is blowing directly from the IRGA to the sonic is",f12.4," deg")',wdtoCSAT
			PRINT*,'--> Separation unit vector (pointing to IRGA) is',ecv%sepunit(:)
			PRINT*,'eof orient()'
			PRINT*,'***************************************'
			WRITE(ecv%loglun,*)"--> Angle between due East and the sonic x-axis is",(alpha/(2*pi))*360," deg"
			WRITE(ecv%loglun,*)"--> Wind direction where the wind is blowing directly from the IRGA to the sonic is",wdtoCSAT," deg "
			WRITE(ecv%loglun,*)"--> Separation unit vector (pointing to IRGA) is",ecv%sepunit(:)
			WRITE(ecv%loglun,*)'eof orient()'
			WRITE(ecv%loglun,*)'***************************************'





			RETURN
		END SUBROUTINE

				
			
				

		! Performs the planar fit routine of (Wilczak et. al. 2001) rotating the coordinate system into
		! the mean stream lines such that in the new planar coordinate w_p; the mean vertical velocity
		! averaged over the entire (non flagged) series is zero.
		SUBROUTINE	planarfit(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			REAL*8				:: b0,b1,b2		! Coefficients (b0 is the vertical instrument offset).
			REAL*8				:: alpha,beta		! Rotation angles: pitch (rotation about y) and roll (rot about x).
			REAL*8				:: pi			! Pi.		
			INTEGER				:: nper			! Points in an averaging period.
			INTEGER				:: apernn		! Integer number number of averaging periods in the series. 
			INTEGER				:: nw			! Number of averaging periods that are weighted in the planar fit.
			REAL*8,DIMENSION(3,3)		:: Am,Aminv		! A matrix and its inverse used in the least squares regression.
			REAL*8,DIMENSION(3,3)		:: eye3,cont		! I the 3x3 identity matrix and a control matrix (should equal I). 
			REAL*8,DIMENSION(3)		:: b,c			! b array containing the b coefficients and c array defining b as the
										! solution of the inhomogeneous matrix equation Ab=c.
			REAL*8,DIMENSION(3)		:: ivec,jvec,kvec	! Planar orthogonal unit vectors expressed IN THE SYSTEM FRAME.
			REAL*8,DIMENSION(3)		:: eastvec,northvec	! Orthogonal unit vectors pointing due east and north IN THE SYSTEM FRAME.
			REAL*8				:: magnitude		! Temporary real to hold vector mangitudes. 
			LOGICAL				:: inc			! Logical declaring if block average is to be considered in the fit.	
			! Mean wind arrays in the sonic coordinates and the planar coordinates.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: mmwds0,mmwds
			INTEGER			:: i,j,k,nnan			! Counters.
			INTEGER			:: lb,ub			! Lower and upper limit of indices in a given averaging period.
			REAL*8			:: nanlim=0.1			! Limit on ratio of nans in a period for it to be weighted.
			! Reals needed for the least squares regression in computing the coefficients b0,b1,b2.
			REAL*8			:: sumum,sumvm,sumumum,sumvmvm,sumumvm,sumwm,sumumwm,sumvmwm,sumwmwm
			! Temporary reals.
			REAL*8			:: lhs,rhs,uu,vv,uv,vw,uw,utemp,vtemp,wtemp	
			! Arrays and reals to Control the final result.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: wpbar	! Block averaged vertical velocities in the planar frame.
			REAL*8					:: wplt,wmlt	! Long term mean vertical velocities in the two frames.	
			! Temporary reals used to compute the constancy ratio for each block.
			REAL*8			:: usum,vsum,miws,cr		! In order: block sum of u, block sum of v, 
										! sum of inst. wind speed and the constancy ratio. 


			PRINT*,'***************************************'
			PRINT*,'Applying the planar fit'
			WRITE(ecv%loglun,*)'***************************************'
			WRITE(ecv%loglun,*)'Applying the planar fit'
			
			! Set the number of points to an averaging period (already computed in setdims).
			nper=ecv%blockl
			
			! Compute pi.
			pi=4.0*ATAN(1.0)

			! Nearest integer number of averaging periods in the entire series.
			apernn=FLOOR(1.0*ecv%nn/nper)

			! Initially assign this to be the first dimension of the measured mean wind array.
			ALLOCATE(mmwds0(apernn,3))

			! Compute the mean winds in each averaging period.
			! Determine if there are enough non-flagged points to weight the given mean in the
			! planar fit.
			
			! Initialize interal counters as zeros.
			j=0
			nnan=0
			DO i=1,apernn	! Dummy counter.
				inc=.TRUE.		! <-- Set the block inclusion to true initially.
				nnan=0
				usum=0
				vsum=0	
				
				! Block bounds.
				lb=(i-1)*nper+1
				ub=i*nper
			
				! First compute the constancy ratio to verify that the wind direction is well defined and test that
				! not more than 10% of the block is a nan value.

				! Test if more than 10% of the data in the block is 'nan'.
				DO k=lb,ub
					IF(ecv%u(k).EQ.ecv%nanreal.OR.ecv%v(k).EQ.ecv%nanreal) THEN
						nnan=nnan+1
					ELSE
						usum=usum+ecv%u(k)				! Block sum of u components.
						vsum=vsum+ecv%v(k)				! Block sum of v components.
						miws=miws+SQRT(ecv%u(k)**2+ecv%v(k)**2)		! Sum of instantaneous wind speeds.
					END IF
				END DO
				! If nnan comprises more than 10% of the block then exclude the block from the planar fit.
				IF(((1.0*nnan)/nper).GE.nanlim) THEN
					inc=.FALSE.
				END IF
				! Also test for the wind direction; if it is within an
				! exclusion bin AND the constancy ratio is hight then set inc to false.
				CR=SQRT(usum**2+vsum**2)/miws

				IF(inc.EQV..TRUE.) THEN ! If the block is still included then:
					! Also test for the wind direction; if it is within an
					! exclusion bin AND the constancy ratio is high then set inc to false.
					
					CR=SQRT(usum**2+vsum**2)/miws	! <-- Compute the constancy ratio.

					DO k=1,SIZE(ecv%disdir,1)
						IF(ecv%mwd(i).GE.ecv%disdir(k,1).AND.ecv%mwd(i).LE.ecv%disdir(k,2)&
								.AND.CR.GT.0.75) THEN
							inc=.FALSE.
						END IF
					END DO
				END IF

				IF(inc.EQV..TRUE.) THEN
					j=j+1	! Weighted block.
					CALL temporalmean(ecv,nper,ecv%u(lb:ub),mmwds0(j,1))	! Block averaged u in sonic coordinates.
					CALL temporalmean(ecv,nper,ecv%v(lb:ub),mmwds0(j,2))	! Block averaged v in sonic coordinates.
					CALL temporalmean(ecv,nper,ecv%w(lb:ub),mmwds0(j,3))	! Block averaged w in sonic coordinates.
				END IF
				nnan=0	! Reset values for each pass.
				usum=0
				vsum=0	
				miws=0
				CR=0
			END DO

			! Set the weighted number of blocks.
			nw=j

			! Only include the weighted means in the measured mean wind array (j entries),
			! so first allocate space in the measured and planar wind arrays:
			ALLOCATE(mmwds(nw,3))
			! Then set the values of the measured wind array.
			mmwds=mmwds0(1:nw,:)


			!*******************************************************************************************
			! Least squares routine.

			! Calculate the entries of the matrix A explicitly.
			sumum=SUM(mmwds(:,1))
			sumvm=SUM(mmwds(:,2))
			sumwm=SUM(mmwds(:,3))
			sumumum=DOT_PRODUCT(mmwds(:,1),mmwds(:,1))
			sumvmvm=DOT_PRODUCT(mmwds(:,2),mmwds(:,2))
			sumwmwm=DOT_PRODUCT(mmwds(:,3),mmwds(:,3))
			sumumvm=DOT_PRODUCT(mmwds(:,1),mmwds(:,2))
			sumumwm=DOT_PRODUCT(mmwds(:,1),mmwds(:,3))
			sumvmwm=DOT_PRODUCT(mmwds(:,2),mmwds(:,3))

			! Populate the matrix symmetric matrix A.
			Am(1,1)=nw
			Am(1,2)=sumum
			Am(1,3)=sumvm
			Am(2,1)=Am(1,2)
			Am(3,1)=Am(1,3)
			Am(2,2)=sumumum
			Am(3,3)=sumvmvm
			Am(2,3)=sumumvm
			Am(3,2)=Am(2,3)

			! Verify that the 'error with matrix inverse' works correctly by setting A
			! to a matrix with no inverse. e.g. the zero matrix. Comment out once checked.
			! -> Am(:,:)=0
			! Checked.

			! Populate the array c.
			c(1)=sumwm
			c(2)=sumumwm
			c(3)=sumvmwm

			! Solve the matrix equation, first compute the inverse of A.
			CALL inverse(Am,Aminv,3)

			! Control that the inverse is correctly calculated using A^-1*A=I.
			! if the inverse is correct use it to calculate the array b
			! using b=A^-1*c
			eye3(:,:)=0
			DO j=1,3
				eye3(j,j)=1
			END DO
			
			cont=MATMUL(Aminv,Am)
			PRINT*,'cont=',cont
			PRINT*,'ISNAN(cont)',ANY(ISNAN(cont))
			IF(ANY(ISNAN(cont)).EQV..true.) THEN
				PRINT*,'ok working'
			END IF
			IF((ANY((ABS(cont-eye3)).GE.1e-3)).EQV..true.) THEN	! If any entry is far off 0.001 raise exception and stop.
				PRINT*,'Error with matrix inverse: >tolerance'
				STOP	
			ELSEIF((ANY(ISNAN(cont))).EQV..true.) THEN	! If any is nan raise exception and stop.
				PRINT*,'Error with matrix inverse: singular'
				STOP
			ELSE
				b(:)=0
				DO i=1,3
					DO j=1,3
						b(i)=b(i)+Aminv(i,j)*c(j)
					END DO
				END DO
			END IF
			! Extract the coefficients from the b array.
			b0=b(1)
			b1=b(2)
			b2=b(3)

			! eof least squares routine.
			!----------------------------------------------------------------------------------------

			!****************************************************************************************
			! Calculate the new ONB defining the planar fit coordinate frame.

			! Calculate the components of the new vertical unit vector expressed in the system frame.
			kvec(1)=-1*b1/(SQRT(b1**2+b2**2+1))
			kvec(2)=-1*b2/(SQRT(b1**2+b2**2+1))
			kvec(3)=1.0/(SQRT(b1**2+b2**2+1))

			! Verify that kvec is a unit vector.
			PRINT*,'|k|=',SQRT(kvec(1)**2+kvec(2)**2+kvec(3)**2)

			! Calculate the north and east unit vectors as defined in the system frame.
			! Note that alpha is the angle between due east and the x unit vector
			! in the system frame.
			eastvec(1)=COS(ecv%alpha)
			eastvec(2)=-SIN(ecv%alpha)
			eastvec(3)=0
			northvec(1)=SIN(ecv%alpha)
			northvec(2)=COS(ecv%alpha)
			northvec(3)=0

			
			! Calculate the entries of unit vectors which define the span of the streamline plane,
			! i.e. ivec and jvec, such that the horizontal (system frame) component of ivec points
			! due east and jvec due north.

			! That is jvec=(kvec X eastvec)/(|kvec X eastvec|).
			! First compute the components:
			jvec(1)=-1*kvec(3)*eastvec(2)	! eastvec(3)=0
			jvec(2)=kvec(3)*eastvec(1)	! eastvec(3)=0
			jvec(3)=kvec(1)*eastvec(2)-kvec(2)*eastvec(1)
			! Then normalize by the magnitude --> unit vector.
			magnitude=SQRT(DOT_PRODUCT(jvec,jvec))
			jvec(1)=jvec(1)/magnitude
			jvec(2)=jvec(2)/magnitude
			jvec(3)=jvec(3)/magnitude
			! Verify that jvec is a unit vector.
			PRINT*,'|j|=',SQRT(jvec(1)**2+jvec(2)**2+jvec(3)**2)

			! Use this to define the third and final unit vector that forms 
			! the new planar ONB, ivec=(jvec X kvec)/(|jvec X kvec|).
			! Again first compute the magnitude.
			ivec(1)=jvec(2)*kvec(3)-jvec(3)*kvec(2)
			ivec(2)=jvec(3)*kvec(1)-jvec(1)*kvec(3)
			ivec(3)=jvec(1)*kvec(2)-kvec(1)*jvec(2)
			! Then normalize by the magnitude --> unit vector.
			magnitude=SQRT(DOT_PRODUCT(ivec,ivec))
			ivec(1)=ivec(1)/magnitude
			ivec(2)=ivec(2)/magnitude
			ivec(3)=ivec(3)/magnitude
			! Verify that ivec is a unit vector.
			PRINT*,'|i|=',SQRT(ivec(1)**2+ivec(2)**2+ivec(3)**2)
			

			! NB. From the above routine ivec should be orthogonal to the 'northvec'
			! while having a component parallel to the 'eastvec'. (i.e. lies in the
			! eastvec, kvec plane).
			! while jvec should be orthogonal to the 'eastvec' while having
			! a component parallel to the 'northvec'. (i.e. lies in the northvec ,kvec plane).
			
			! Verify that the new ONB [i,j,k] defining the planar fit coordinate frame is truly orthogonal,
			! and also ensure that i lies in the [z,e] plane and that j lies in the [z,n] plane
			IF(ABS(DOT_PRODUCT(jvec,kvec)).GE.1e-4.OR.&
					ABS(DOT_PRODUCT(ivec,kvec)).GE.1e-4.OR.&
					 	ABS(DOT_PRODUCT(ivec,jvec)).GE.1e-4) THEN
					PRINT*, '~~~~~Obs! Planar frame is not an ONB'
					PRINT*, '~~~~~j dot k=',DOT_PRODUCT(jvec,kvec)
					PRINT*, '~~~~~i dot k=',DOT_PRODUCT(ivec,kvec)
					PRINT*, '~~~~~i dot j=',DOT_PRODUCT(ivec,jvec)
					PRINT*, '~~~~~i dot e=',DOT_PRODUCT(ivec,eastvec)
					PRINT*, '~~~~~j dot e=',DOT_PRODUCT(jvec,eastvec)
					WRITE(ecv%loglun,*)'~~~~~Obs! Planar frame is not an ONB'
					WRITE(ecv%loglun,*)'~~~~~j dot k=',DOT_PRODUCT(jvec,kvec)
					WRITE(ecv%loglun,*)'~~~~~i dot k=',DOT_PRODUCT(ivec,kvec)
					WRITE(ecv%loglun,*)'~~~~~i dot j=',DOT_PRODUCT(ivec,jvec)
					WRITE(ecv%loglun,*)'~~~~~i dot e=',DOT_PRODUCT(ivec,eastvec)
					WRITE(ecv%loglun,*)'~~~~~j dot e=',DOT_PRODUCT(jvec,eastvec)
			ELSE 
					PRINT*, 'Planar frame is an ONB'
					PRINT*, '-->j dot k=',DOT_PRODUCT(jvec,kvec)
					PRINT*, '-->i dot k=',DOT_PRODUCT(ivec,kvec)
					PRINT*, '-->i dot j=',DOT_PRODUCT(ivec,jvec)
					PRINT*, '-->i dot e=',DOT_PRODUCT(ivec,eastvec)
					PRINT*, '-->j dot e=',DOT_PRODUCT(jvec,eastvec)
					WRITE(ecv%loglun,*)'Planar frame is an ONB'
					WRITE(ecv%loglun,*)'~~~~~j dot k=',DOT_PRODUCT(jvec,kvec)
					WRITE(ecv%loglun,*)'~~~~~i dot k=',DOT_PRODUCT(ivec,kvec)
					WRITE(ecv%loglun,*)'~~~~~i dot j=',DOT_PRODUCT(ivec,jvec)
					WRITE(ecv%loglun,*)'~~~~~i dot e=',DOT_PRODUCT(ivec,eastvec)
					WRITE(ecv%loglun,*)'~~~~~j dot e=',DOT_PRODUCT(jvec,eastvec)
				
			END IF

			!***********************************************************************************************

			!-----------------------------------------------------------------------------------------------
			! Finally rotate the coordinate system into the ensemble mean streamline plane; i.e. the planar fit.
			
			DO i=1,ecv%nn
				IF(ecv%u(i).EQ.ecv%nanreal&	! As usual exclude spikes from the computation.
					.OR.ecv%v(i).EQ.ecv%nanreal.OR.ecv%w(i).EQ.ecv%nanreal) THEN
					ecv%u(i)=ecv%nanreal
					ecv%v(i)=ecv%nanreal
					ecv%w(i)=ecv%nanreal
				ELSE
					! Compute u in the planar frame, based on simple dot products:
					! up=i1*u + i2*v + i3*(w-b0)
					! Where we have taken offset into account and assumed that
					! for the horizontal components the offset is negligible.
					
					utemp=ivec(1)*ecv%u(i)+ivec(2)*ecv%v(i)+ivec(3)*(ecv%w(i)-b0)
					
					! Similarly for v:
					! vp=j1*u+j2*v+j3*(w-b0)

					vtemp=jvec(1)*ecv%u(i)+jvec(2)*ecv%v(i)+jvec(3)*(ecv%w(i)-b0)
				
					! And w:
					! wp=k1*u+k2*v+k3*(w-b0)
					wtemp=kvec(1)*ecv%u(i)+kvec(2)*ecv%v(i)+kvec(3)*(ecv%w(i)-b0)

					! Now that all the planar fit values are computed for the given point in time,
					! update to velocities to these values.
					ecv%u(i)=utemp
					ecv%v(i)=vtemp
					ecv%w(i)=wtemp


				END IF
			END DO

			! Last we control the result. That is verify how small the long term mean
			! of the block averaged PLANAR vertical velocity is.
			ALLOCATE(wpbar(nw))
			DO j=1,nw
				! Calculate the block averaged vertical velocity in the planar frame.
				wpbar(j)=kvec(1)*mmwds(j,1)+kvec(2)*mmwds(j,2)+kvec(3)*(mmwds(j,3)-b0)	
			END DO
			CALL temporalmean(ecv,nw,wpbar,wplt)		! Calculate the LT mean of planar block averaged vertical velocities. 
			CALL temporalmean(ecv,nw,mmwds(:,3),wmlt)	! Calculate the LT mean of measured block averaged vertical velocities.


			! Compute the pitch and the roll and print these to the terminal.
			alpha=(ATAN(-b1)/(2.0*pi))*360
			beta=(ATAN(b2)/(2.0*pi))*360
		
			PRINT*,'kvec=',kvec(:)
			PRINT"('-->Pitch is:',f12.4,' degrees (rotated about sonic y axis)')",alpha
			PRINT"('-->Roll is:',f12.4,' degrees (rotated about sonic x axis)')",beta
			PRINT"('-->Mean of block averaged system vertical velocities=',f12.4,' [m/s]')",wmlt
			PRINT"('-->Mean of block averaged planar vertical velocities=',f12.4,' [m/s]')",wplt
			PRINT"('-->Offset is:',f12.4,' [m/s]')",b0
			WRITE(ecv%loglun,*)'kvec=',kvec(:)
			WRITE(ecv%loglun,*)'-->Pitch is:',alpha,' degrees (rotated about sonic y axis)'
			WRITE(ecv%loglun,*)'-->Roll is:',beta,' degrees (rotated about sonic x axis)'
			WRITE(ecv%loglun,*)'-->Mean of block averaged system vertical velocities=',wmlt,' [m/s]'
			WRITE(ecv%loglun,*)'-->Mean of block averaged planar vertical velocities=',wplt,' [m/s]'
			WRITE(ecv%loglun,*)'-->Offset is:',b0,' [m/s]'


			PRINT*,'eof planarfit()'
			PRINT*,'***************************************'
			WRITE(ecv%loglun,*)'eof planarfit()'
			WRITE(ecv%loglun,*)'***************************************'
		
			RETURN
		END SUBROUTINE

		! Rotates the coordinate system in the streamwise plane (i.e. kvec
		! does not change) into the block or window averaged mean wind; i.e.
		! the mean of urot and vrot input vectors. urot and vrot are returned
		! as the alongstream (urot) and crosstream (vrot) winds respectively.
		! This is not included in the planar fit algorithm to enable
		! separate calls for this rotation whenver required; due
		! to different block lengths in spectra and statistical
		! calculations.
		SUBROUTINE rotatestream(ecv,ss,urot,vrot)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER,INTENT(IN)			:: ss		! Length of the input plane wind arrays (uin,vin).
			REAL*8,DIMENSION(ss),INTENT(INOUT)	:: urot,vrot
			REAL*8,DIMENSION(2)			:: svec,nvec	! Unit vectors to be calculated describing
										! the streamwise and crossstream directions
										! in the sonic reference frame.
			REAL*8					:: ubar,vbar	! Temporary reals to hold the mean winds in the
										! two directions defining the streamwise plane.
			REAL*8					:: m		! Magnitude of the mean wind.
			REAL*8					:: ute,vte	! Temporary reals to hold the old values of u and v.
			INTEGER					:: i		! Counter.
			
			! Compute the mean wind components.
			CALL temporalmean(ecv,ss,urot,ubar)
			CALL temporalmean(ecv,ss,vrot,vbar)
			
			! Compute the magnitude of the mean wind.
			m=SQRT(ubar**2+vbar**2)

			! Use the planar coordinate basis
			! to first express the unit vector pointing in the mean wind direction
			! in terms of the unit vectors i,j then use the cross product and
			! to compute the crossstream unit vector.
	
			! Alongstream unit vector:
			svec(1)=ubar/M
			svec(2)=vbar/M
			
			! Cross stream unit vector, using k x s = (-s2)i+(s1)j 
			nvec(1)=-1*svec(2)
			nvec(2)=svec(1)

			! Finally define the (high frequency) velocities according to these unit 
			! vectors.
			DO i=1,ss
				IF(urot(i).NE.ecv%nanreal.AND.vrot(i).NE.ecv%nanreal) THEN
					ute=urot(i)*svec(1)+vrot(i)*svec(2)	! Streamwise velocity: u_s = s dot (ui +vj)
					vte=urot(i)*nvec(1)+vrot(i)*nvec(2)	! Cross stream velocity: u_n= n dot (ui + vj)
				ELSE
					ute=ecv%nanreal
					vte=ecv%nanreal
				END IF

				! Now that these are computed, update the velocities to the rotated values.
				urot(i)=ute
				vrot(i)=vte

				ute=0
				vte=0
				
			END DO
			CALL temporalmean(ecv,ss,vrot,vbar)
			CALL temporalmean(ecv,ss,urot,ubar)

			ubar=0
			vbar=0
			svec(:)=0
			nvec(:)=0
	

			RETURN
		END SUBROUTINE
				

			
			
				

		! Slightly modified version of a subroutine by Alex G. made in 2009, which
		! computes the inverse of a square matrix using Doolittle LU factorization.
		! input: square matrix 'm', size (where size1=size2=) 'n'.
		! output: matrix inverse 'minv'.
  		SUBROUTINE inverse(m,minv,n)
			IMPLICIT NONE 
			INTEGER,INTENT(IN)	:: n		! Size of input square matrix; n x n.
			REAL*8,INTENT(IN)	:: m(n,n)	! Input square matrix n x n.
			REAL*8,INTENT(OUT)	:: minv(n,n)	! Inverse of the square n x n matrix.			
			REAL*8			:: a(n,n)	! Copy of m (changes during computation).
			REAL*8			:: c(n,n)	! Copy of minv.
			REAL*8			:: L(n,n), U(n,n), b(n), d(n), x(n)
			REAL*8			:: coeff	! Dummy coefficients.
			INTEGER			:: i, j, k	! Counters.

			! Allocate the work matrix used to find the inverse.
			a=m

			! step 0: initialization for matrices L and U and b
			! Fortran 90/95 allows such operations on matrices
			L=0.0
			U=0.0
			b=0.0

			! step 1: forward elimination
			DO k=1, n-1
   				DO i=k+1,n
      					coeff=a(i,k)/a(k,k)
      					L(i,k) = coeff
      					DO j=k+1,n
        			 		a(i,j) = a(i,j)-coeff*a(k,j)
      					END DO
   				END DO
			END DO

			! Step 2: prepare L and U matrices 
			! L matrix is a matrix of the elimination coefficient
			! + the diagonal elements are 1.0
			DO i=1,n
  				L(i,i) = 1.0
			END DO
			! U matrix is the upper triangular part of A
			DO j=1,n
  				DO i=1,j
    					U(i,j) = a(i,j)
  				END DO
			END DO
	
			! Step 3: compute columns of the inverse matrix C
			DO k=1,n
  				b(k)=1.0
  				d(1) = b(1)
				! Step 3a: Solve Ld=b using the forward substitution
  				DO i=2,n
    					d(i)=b(i)
    					DO j=1,i-1
      						d(i) = d(i) - L(i,j)*d(j)
    					END DO
  				END DO
				! Step 3b: Solve Ux=d using the back substitution
  				x(n)=d(n)/U(n,n)
  				DO i = n-1,1,-1
    					x(i) = d(i)
    					DO j=n,i+1,-1
      						x(i)=x(i)-U(i,j)*x(j)
    					END DO
    					x(i) = x(i)/u(i,i)
 				END DO
				! Step 3c: fill the solutions x(n) into column k of C
  				DO i=1,n
    					c(i,k) = x(i)
  				END DO
 					 b(k)=0.0
			END DO
			! Set the output inverse.
			minv=c
			RETURN
	
		END SUBROUTINE inverse




		! Reads the (slow-meteorology) observed pressure from the NYÅ weather station (eKlima .txt file)
		! that has been converted to a netcdf file in the folder ../input, this input file has
		! four entries per day @ UTC: 00:00,06:00,12:00,18:00 for the period 01.01.2007-31.12.2010.
		! Tracks the time in the ecv%t array and assigns synchronized pressure to the ecv%p array.
		! The observed pressure measured at hour XX:XX in the eKlima file is allocated to all ecv%p with
		! times in XX:XX-02:59.55 to XX:XX+03:00.00. That is the pressures are centered about the
		! observation times.
		SUBROUTINE syncpressure(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER,ALLOCATABLE,DIMENSION(:,:)	:: tp		! Time array to be read from .nc file.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: pOp		! Pressure field to be read from .nc file.
			REAL*8,PARAMETER			:: pst=1013.25	! Standard atmospheric pressure, placeholder value.
			CHARACTER(LEN=80)			:: pfile	! .nc filename.
			CHARACTER(LEN=*),PARAMETER		:: tname='t'	! Time array name in .nc file.
			CHARACTER(LEN=*),PARAMETER		:: pname='pO'	! Pressure array name in .nc file.
			CHARACTER(LEN=*),PARAMETER		:: rown='length'! Variable dimension name.
			CHARACTER(LEN=*),PARAMETER		:: coln='unit'	! Variable dimension name.
			INTEGER					:: ncols,nrows	! Size of arrays in the file (tbd).
			INTEGER					:: colid,rowid	! Dimension ids.
			INTEGER					:: ncid,v1,v2	! Variable and netCDF file ids.
			INTEGER					:: lyis,uyis	! Year of first and last ecv entry.
			INTEGER					:: lmis,umis	! Month of first and last ecv entry.
			INTEGER					:: ldis,udis	! Day of first and last ecv entry.
			INTEGER					:: lhis,uhis	! Hour of first and last ecv entry.
			INTEGER					:: lhc		! lhis in the (right) closed interval (lhc-3,lhc+3].
			INTEGER					:: uhc		! uhis in the (right) closed interval (uhc-3,uhc+3].
			INTEGER					:: i,j,k,m,l,n	! Counters.
			INTEGER					:: ll,ul	! Bounds.
			INTEGER					:: fl,fu	! Status for the bounds.
			INTEGER					:: np6		! Number of entries in 6 hours.
			INTEGER,DIMENSION(4)			:: hoursO	! UTC hours with observations.

		
			!------------------------------------------------------------------------------------------------
			! Reading from the observed pressure file.
			!------------------------------------------------------------------------------------------------

			! Set the filename and path.
			pfile=ecv%pfile


			PRINT*,'---------------------------------------------------------------------'
			PRINT*,'Reading from file:',pfile
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'Reading from file:',pfile
			

			! Open the .nc file.
			CALL CHECK( nf90_open(pfile,NF90_NOWRITE,ncid) )

			! Get the variable ids.
			CALL CHECK( nf90_inq_varid(ncid,tname,v1) )
			CALL CHECK( nf90_inq_varid(ncid,pname,v2) )

			! Get variable dimension id. Needs the name of the dimension, d1n (from e.g. ncview).
			CALL CHECK( nf90_inq_dimid(ncid,rown,rowid) )
			CALL CHECK( nf90_inq_dimid(ncid,coln,colid) )

			! Get variable dimensions.
			CALL CHECK( nf90_Inquire_Dimension(ncid,rowid,len=nrows) )
			CALL CHECK( nf90_Inquire_Dimension(ncid,colid,len=ncols) )

			! Allocate temp arrays.
			ALLOCATE(tp(nrows,ncols))
			ALLOCATE(pOp(nrows))

			! GET the variables.
			CALL CHECK( nf90_get_var(ncid,v1,tp) )
			CALL CHECK( nf90_get_var(ncid,v2,pOp) )

			! Close the file, freeing all resources.
 			CALL check( nf90_close(ncid) )

			!------------------------------------------------------------------------------------------------
			! Eddy covariance slow pressure allocation.
			!------------------------------------------------------------------------------------------------


			! Allocate the dimensions of the pressure array.
			ALLOCATE(ecv%p(ecv%nn))
		
			! Allocate the bounds of the ecv%p array, based on the first and last
			! occurence of years, months and days.
			lyis=NINT(ecv%t(1,1)) 		! ecv%t is type real*8, hence NINT is used.
			uyis=NINT(ecv%t(ecv%nn,1))
			lmis=NINT(ecv%t(1,2))
			umis=NINT(ecv%t(ecv%nn,2))
			ldis=NINT(ecv%t(1,3))
			udis=NINT(ecv%t(ecv%nn,3))
			lhis=NINT(ecv%t(1,4))
			uhis=NINT(ecv%t(ecv%nn,4))

			

			! Find the bounds for year,month day in the observed pressure array so as to synch
			! with the ecv pressure array.
			
			! Set found bounds status to not found (yet).
			fl=0			
			fu=0
			i=1			
			DO WHILE(fl.EQ.0.OR.fu.EQ.0)	! Loop through the observed pressure array.
				! Lower bounds:
				! If first entry in the correct year then:
				IF(i.NE.1.AND.fl.EQ.0.AND.tp(i,1).EQ.lyis) THEN
					k=0	! <-- Found month status.
					j=0	! <-- Entry counter.
					DO WHILE(k.EQ.0)
						! If first entry in the correct month: 
						IF(tp(i+j,2).EQ.lmis) THEN
							l=0	! <-- Found day status.
							m=0	! <-- Entry counter.
							DO WHILE(l.EQ.0)
								! If first entry in the correct day then:
								IF(tp(i+j+m,3).EQ.ldis) THEN
									ll=i+j+m	! Found the lower bound.
									l=1		! Set to found day.
									fl=1		! Set to found lower bound.
								END IF
								! Otherwise update the counter.
								m=m+1
							END DO
							! Set to found.
							k=1
						END IF
						! Otherwise update the counter.
						j=j+1
					END DO
				! If first index:
				ELSEIF(i.EQ.1.AND.tp(i,3).EQ.ldis.AND.tp(i,2).EQ.lmis.AND.tp(i,1).EQ.lyis) THEN
					ll=i
				END IF


				! Upper bounds:
				! If first entry in the correct year then:
				IF(i.NE.SIZE(pOp).AND.fu.EQ.0.AND.tp(i,1).EQ.uyis) THEN
					k=0	! <-- Found month status.
					j=0	! <-- Entry counter.
					DO WHILE(k.EQ.0)
						! If first entry in the correct month then:
						IF(tp(i+j,2).EQ.umis) THEN
							l=0	! <-- Found day status.
							m=0	! <-- Entry counter.
							DO WHILE(l.EQ.0)
								! If LAST entry in the correct day then:
								IF(tp(i+j+m,3).EQ.udis.AND.tp(i+j+m+1,3).NE.udis) THEN
									ul=i+j+m	! Found the upper bound.
									l=1		! Set to found day.
									fu=1		! Set to found upper bound.
								END IF
								! Otherwise update the counter.
								m=m+1
							END DO
							! Set to found.
							k=1
						END IF
						! Otherwise update the counter.
						j=j+1
					END DO
				! If last index:
				ELSEIF(i.EQ.SIZE(pOp).AND.tp(i,3).EQ.udis.AND.tp(i,2).EQ.umis.AND.tp(i,1).EQ.uyis) THEN
					ul=i
				END IF	
				! Update the counter.
				i=i+1
			END DO

			! Now synchronize the series.
			! Compute the number of entries in 6 hours, with f_s=20 Hz.
			np6=ecv%fs*6*60**2	! <--- is even.

			! Allocate the observation hours.
			hoursO(1)=0
			hoursO(2)=6
			hoursO(3)=12
			hoursO(4)=18

			! Update the lower bounds of the observed pressure array to account for the first hour entry
			! of the ecv time series. 
			DO j=1,SIZE(hoursO)
				IF(hoursO(j)-3.LT.lhis.AND.lhis.LE.hoursO(j)+3) THEN
					lhc=hoursO(j)
				END IF
				IF(hoursO(j)-3.LT.uhis.AND.uhis.LE.hoursO(j)+3) THEN
					uhc=hoursO(j)
				END IF
			END DO


			! Set found status to zero initially.
			fl=0
			fu=0
			! Update the bounds of the observed pressure array.
			i=0	! Counter.
			DO WHILE(fl.EQ.0.OR.fu.EQ.0)
				! Lower bounds.
				IF(tp(ll+i,4).EQ.lhc.AND.fl.EQ.0) THEN	! Found new lower bounds.
					ll=ll+i
					fl=1
				END IF
				! Upper bounds (scan from bottom).
				IF(tp(ul-i,4).EQ.uhc.AND.fu.EQ.0) THEN	! Found new upper bounds.
					ul=ul-i
					fu=1
				END IF
				! Update the counter for the next pass.
				i=i+1
			END DO

			! Search for the first occurence of one of the observation hours in the ecv%t series.
			i=0
			k=0
			n=0
			DO WHILE(k.EQ.0)
				i=i+1
				DO j=1,SIZE(hoursO)
					IF(ecv%t(i,5).EQ.0.AND.ecv%t(i,6).EQ.0) THEN 	! Top of the hour.
						IF(hoursO(j).EQ.NINT(ecv%t(i,4))) THEN	! Observation hour has occured.
							k=i
							IF(hoursO(j).NE.lhc) THEN	! <--- Unlikely to be within the first bin.
								n=1			! Account for this with a status counter.
							END IF
								
						END IF
					END IF
				END DO
			END DO

			! Allocate the observed pressure to the eddy covariance slow pressure array (synchronized).
			DO i=1,SIZE(pOp(ll:ul))
				! Set the bounds in the ecv pressure array, centered about entry k.
				l=k-(np6/2-1)
				m=k+np6/2
				! Allocate the values of the ecv pressure array.
				IF(i.EQ.1.AND.n.EQ.1.AND.l.GT.0) THEN ! First bin.
					ecv%p(1:l)=pOp(ll+(i-1))	! Not part of l,m bounds above.
				ELSEIF(i.EQ.1.AND.l.LE.0) THEN
					ecv%p(1:m)=pOp(ll+(i-1))	! Part of 1:m bounds above.
				ELSEIF(m.GT.ecv%nn.AND.l.GT.0)	THEN		
					ecv%p(l:ecv%nn)=pOp(ll+(i-1))	! m exceeds array bounds.
				ELSEIF(l.GT.0)	THEN
					ecv%p(l:m)=pOp(ll+(i-1))	! Part of l:m bounds above (standard case).	
				END IF
				! Update k for the next pass of the loop, i.e. the next point to center the pressure
				! array about.
				IF(i.NE.1.OR.l.LE.0) THEN
					k=k+np6
				END IF
			END DO

			! Due to the short periods where the NYÅ met station was malfunctioning implying missing pressure readings
			! which are set to -1000 (NaN) in the file we change these to standard atmospheric pressure. The reason
			! for treating pressure differently, i.e. if it were another variable all other variables would also be
			! flagged as faulty, is that we would risk discarding far too much data. One pressure reading comprises
			! 6 hours of data, i.e. 4.32x10^5 points at 20 Hz sampling frequency, and we do not wish to discard
			! that much just because the pressure is slightly off. Note this has no effect on most of the files as
			! the downtime of the NYÅ station is minimal.
			DO j=1,SIZE(ecv%p)
				IF(ecv%p(j).EQ.ecv%nanreal) THEN
					ecv%p(j)=pst
				END IF
			END DO


			! As with the other variables set to NaN if any of the other variables are NaN.
			DO j=1,SIZE(ecv%p)
				IF(ecv%u(j).EQ.ecv%nanreal) THEN
					ecv%p(j)=ecv%nanreal
				END IF
				IF(ecv%p(j).EQ.ecv%nanreal) THEN	! Double check that NaNs are properly assigned.
					ecv%u(j)=ecv%nanreal
					ecv%v(j)=ecv%nanreal
					ecv%w(j)=ecv%nanreal
					ecv%Ts(j)=ecv%nanreal
					ecv%h(j)=ecv%nanreal
				ELSE
					! Convert from hPa to Pa.
					ecv%p(j)=1e2*ecv%p(j)
				END IF
			END DO

			PRINT*,'->Done synchronizing the slow pressure'
			PRINT*,'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'->Done synchronizing the slow pressure'
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'

			! Deallocate
			DEALLOCATE(tp)
			DEALLOCATE(pOp)


			RETURN
		END SUBROUTINE


		!// Reads in the diurnaly interpolated measurement heights based on hourly measurement height measurements by the Campbell SR50 sonic 			! ranger or the Bayelva Climate station if such a file is available. These heights are then allocated to the array zmes with the same 			! length as the time array.
		SUBROUTINE syncsnowdepth(ecv) 
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			! Internal declarations
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: ztemp	! Measurement heights read in from the file.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: ttemp	! Time stamps read in from the file (one for each day). Format y m d.
			INTEGER					:: eoff		! End of file status (0,1).
			CHARACTER(LEN=4)			:: seoff
			INTEGER					:: i,j,k	! Counters.
			INTEGER					:: nlines	! Number of lines in file.
			CHARACTER(LEN=160)			:: line,subline ! Line entries.
			CHARACTER(LEN=80)			:: zfileis
			INTEGER					:: nn		! Number of lines of EC data.
			INTEGER,PARAMETER			:: zlun=11	! Logical unit number for the measurement file.
			INTEGER					:: res		! Result of file operation.
			CHARACTER(LEN=1)			:: deli=' '	! Delimiter in the zfile.
		
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'***Searching for measurement heights***'
			PRINT*,'---------------------------------------------------------------------'
			PRINT*,'***Searching for measurement heights***'

			!%%%%%%%%%%%%%% Case 1: Available measurement heights %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			IF(ecv%hmes.EQV..TRUE.) THEN	! < IF AVAILABLE
			zfileis=ecv%zfile
			WRITE(ecv%loglun,*)'-->Reading and synchronizing measurement heights from file:',zfileis
			PRINT*,'-->Reading and synchronizing measurement heights from file:',zfileis

			
			!---------z mesfile---------------------------
			! Open the zfile.
 			OPEN(UNIT=zlun,FILE=ecv%zfile,FORM='FORMATTED',&
				IOSTAT=res)
 			IF(res /= 0) THEN
				PRINT *, 'Error in opening z file, status:', res
  				STOP
 			END IF
		
			! Find the total number of lines in the file by reaching the end of file string 'eoff'.
			eoff=0
			seoff='eoff'
			i=0
			DO WHILE(eoff.EQ.0)
				i=i+1		! Line counter.
				READ(zlun,FMT='(A)',IOSTAT=res) line
				IF(TRIM(line).EQ.seoff) THEN
					eoff=1
					i=i-1	! Dont count the last line.
				END IF
			END DO
			nlines=i
			! Allocate the temporary arrays:
			ALLOCATE(ttemp(nlines,3))
			ALLOCATE(ztemp(nlines))
			
			REWIND(zlun)	! Rewind to the start of the file.
			
			! Now read in the time stamps and the heights in the file.
			DO i=1,nlines	! Loop over all days in the file.
				READ(zlun,FMT='(A)',IOSTAT=res) line
				DO k=1,3	! Read in the time stamps.
					j=INDEX(line,deli)
					subline=line(1:j)
					READ(subline,*) ttemp(i,k)
					line=line(j+1:)	! Truncate the line for the next pass.
				END DO
				READ(line,*)	ztemp(i)	! Read in the measurement height for the given time stamp.
			END DO

			CLOSE(UNIT=zlun) ! Close the zfile.


			! Now we synchronize the ztemp array with a zmes array with the same dimensions as the
			! other high frequency arrays (u,v,w,Ts,h and so on).
			nn=SIZE(ecv%u,1)
			ALLOCATE(ecv%zmes(nn))
			DO i=1,nn	! Loop over the EC data.
				IF(i.EQ.1) THEN	! Search for the same day as the first EC entry in the height measurements.
					k=0	! 'Found' integer.
					j=0	! Index in height measurement file.
					DO WHILE(k.EQ.0)
						j=j+1
						! Check if we are in the same day.
						IF(ttemp(j,1).EQ.ecv%t(i,1).AND.ttemp(j,2).EQ.ecv%t(i,2).AND.&
							ttemp(j,3).EQ.ecv%t(i,3)) THEN
							ecv%zmes(i)=ztemp(j)	! Set the measurement height.
							k=j	! Exit while loop and set position in height measurement file.
						END IF
					END DO
				ELSE
					IF(ecv%t(i,3).EQ.ecv%t(i-1,3)) THEN ! IF we are still in the same day.
						ecv%zmes(i)=ztemp(k)
					ELSE				    ! Update day in height measurement file.
						k=k+1
						ecv%zmes(i)=ztemp(k)
					END IF
				END IF
				IF(ecv%zmes(i).EQ.0.OR.ISNAN(ecv%zmes(i)).EQV..TRUE.) THEN
					ecv%zmes(i)=ecv%zbase
				END IF
			END DO
			WRITE(ecv%loglun,*)'-->Max in file z=',MAXVAL(ecv%zmes),' [m]'
			PRINT*,'-->Max zmes in file:',MAXVAL(ecv%zmes),' [m]'
			WRITE(ecv%loglun,*)'-->Min in file z=',MINVAL(ecv%zmes),' [m]'
			PRINT*,'-->Min zmes in file:',MINVAL(ecv%zmes),' [m]'

			DEALLOCATE(ttemp)
			DEALLOCATE(ztemp)

			ELSE
			!%%%%%%%%%%%%%% Case 2: Not available %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			WRITE(ecv%loglun,*)'-->No measurement height file found.'
			PRINT*,'-->No measurement height file found.'

			! If not available set to snow-free measurement height.
			ALLOCATE(ecv%zmes(SIZE(ecv%u,1)))
			ecv%zmes(:)=ecv%zbase

			PRINT*,'-->Set all heights to base (snow-free) height of z=',ecv%zbase
			WRITE(ecv%loglun,*)'-->Set all heights to base (snow-free) height of z=',ecv%zbase

			END IF

			PRINT*,'->Done allocating measurement heights'
			PRINT*,'---------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'->Done allocating measurement heights'
			WRITE(ecv%loglun,*)'---------------------------------------------------------------------'
			

				
			RETURN	
		END SUBROUTINE	





		SUBROUTINE	diaganc(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER					:: j,jj			! Counters.
			REAL*8,PARAMETER			:: Rd=287.04		! Gas constant for dry air.

			PRINT*,'----------------------------------------------------------------------------------'
			PRINT*,'diaganc(): Diagnosing ancillary variables (check for warnings).'	
			WRITE(ecv%loglun,*)'----------------------------------------------------------------------------------'
			WRITE(ecv%loglun,*)'diaganc(): Diagnosing ancillary variables (check for warnings).'		

			! Allocate the size of the instantaneous arrays to be diagnosed.
			jj=SIZE(ecv%Ts)
			ALLOCATE(ecv%q(jj))	! Specific humidty [kg/kg].
			ALLOCATE(ecv%Tv(jj))	! Virtual temperature [K].
			ALLOCATE(ecv%TT(jj))	! Absolute temperature [K].
			ALLOCATE(ecv%rho(jj))	! Air density [kgm^-3].


			! Loop through the data:			
			DO j=1,jj
				IF(ecv%Ts(j).EQ.ecv%nanreal.OR.ecv%p(j).EQ.ecv%nanreal&
					.OR.ecv%h(j).EQ.ecv%nanreal) THEN	! Set to spike if either Ts, p or h are flagged.
					ecv%q(j)=ecv%nanreal
					ecv%Tv(j)=ecv%nanreal
					ecv%TT(j)=ecv%nanreal
					ecv%rho(j)=ecv%nanreal
				ELSE	
					ecv%q(j)=(ecv%p(j)/(Rd*ecv%h(j)*ecv%Ts(j))-0.1)**(-1)
					ecv%TT(j)=ecv%Ts(j)/(1+ecv%q(j)*0.51)
					ecv%Tv(j)=ecv%TT(j)*(1+0.61*ecv%q(j))
					ecv%rho(j)=ecv%h(j)/ecv%q(j)
				END IF
			END DO

			PRINT*,'Done diaganc().'
			PRINT*,'----------------------------------------------------------------------------------'	
			WRITE(ecv%loglun,*)'Done diaganc().'
			WRITE(ecv%loglun,*)'----------------------------------------------------------------------------------'

			RETURN
		END SUBROUTINE		

		! Calculates the time lag between the LI7500 and the CSAT3 in a given block based on wind direction, with a limit on calm wind and
		! constancy, using the maximized cross-correlation technique (see e.g. Vickers&Mahrt 1997). Which of the series to lag according to
		! the other is based on the block averaged wind direction. Once the (block constant) time lag is determined the corresponding
		! covariance is output by the routine (hw) along with the number of indices that the h series is lagged according to the w series (lagli).
		SUBROUTINE separation(ecv,ss,u,wd,h,w,hw,lagli,CFsep)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: ss		! Length of input arrays.
			REAL*8,INTENT(IN)		:: u		! Input block average horizontal wind speed.
			REAL*8,INTENT(IN)		:: wd		! Input block average wind direction.
			REAL*8,DIMENSION(ss),INTENT(IN)	:: h,w		! Input arrays: high frequency absolute humidity and vertical velocity.
			REAL*8,INTENT(OUT)		:: hw		! Output lagged covariance.
			INTEGER,INTENT(OUT)		:: lagli	! Output lag index; if positive LI is downwind of CSAT and vise versa if negative.
			REAL*8,INTENT(OUT)		:: CFsep	! Correction factor.
			REAL*8				:: mR,Rte	! maximum (absolute) cross-correlation and temporary cross-correlation.
			INTEGER				:: i,j,k	! Counters.
			INTEGER				:: tlag		! Maximum number of indices to perturb the series in the cross-correlation
									! analysis.
			REAL*8,DIMENSION(3)		:: bunit	! Block unit wind vector; pointing in the direction that the block average
									! wind is blowing in the SLEC frame. 
			REAL*8				:: dp		! Dot product between block unit vector and the separation unit vector.
			REAL*8				:: madp		! Minimum absolute value of the dot product between the two vectors.			
			REAL*8				:: pi,wdxy	! Pi and wd in a SLEC coordinates but right handed sense.
			REAL*8				:: old		! Uncorrected flux value.
			REAL*8				:: ute		! Temporary real to hold u.
			REAL*8				:: nnans	! Temporary real to hold the number of nans.
			REAL*8				:: wbar,hbar	! Temporary reals to hold the block averages.
			REAL*8				:: sdw,sdh	! Temporary reals to hold the block standard deviations.
			REAL*8				:: covhw	! Temporary real to hold the covariance at zero lag.
			REAL*8				:: ccovhw	! Temporary real to hold the LDT cross covariance at a given lag.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: tempa ! Temporary array for gap filling.
			REAL*8,ALLOCATABLE,DIMENSION(:)	:: hs,ws	! Temporary arrays to hold h and w allowing for conditioning (LDT,gapfill).

			
			! Compute pi.
			pi=4.0*ATAN(1.0)
			
			! Assign 'right handed' wind direction and u as zero initially.
			wdxy=0
			ute=0
		
			! Convert block wind direction bearing to the angle between the east axis and
			! the block average wind direction (right handed sense: anticlockwise from east).
			IF(wd.LE.90.AND.wd.NE.ecv%nanreal) THEN
				wdxy=90-wd
			ELSEIF(wd.LE.360.AND.wd.GT.90.AND.wd.NE.ecv%nanreal) THEN
				wdxy=450-wd
			END IF
			! Calculate the block unit wind vector; sign reversal to convert from wind direction (bearing) to the
			! direction that the wind is blowing to.
			IF(wd.NE.ecv%nanreal) THEN	! Make sure that the wind direction is well defined (i.e. steady and not too calm a wind).
				bunit(1)=-1*COS(wdxy*(1.0/360.0)*2*pi)
				bunit(2)=-1*SIN(wdxy*(1.0/360.0)*2*pi)
				bunit(3)=0
			ELSE
				bunit(:)=ecv%nanreal
			END IF

			dp=DOT_PRODUCT(bunit,ecv%sepunit)	! Dot product between the block unit wind vector and the separation vector.
								! IF >0 the wind has a component blowing from the CSAT to the LI and if
								! <0 the wind has a component blowing from the LI to the CSAT.
			
			madp=0.5				! Threshold for the minimum value of the absolute dot product between the two
								! unit vectors; corresponds to an angle of <=60 deg between the separation vector
								! and the block wind vector. If exceeded no cross-correlation analysis is needed
								! as the wind is blowing almost directly into both instruments; i.e. the
								! component of the wind blowing between the two instruments is neglibible.

			! Set u to the calm wind threshold if it is below this threshold to ensure that tlag is bounded.
			IF(u.LT.0.2) THEN
				ute=0.2
			ELSE 
				ute=u
			END IF

			! Gap fill the arrays for convenience. [<= 10% spike entries so minimal effect].
			ALLOCATE(tempa(ss,2))
			ALLOCATE(hs(ss))
			ALLOCATE(ws(ss))
			hs=h
			ws=w
			tempa(:,1)=hs(:)
			tempa(:,2)=ws(:)
			CALL gapfill(ecv,ss,2,tempa)
			hs(:)=tempa(:,1)
			ws(:)=tempa(:,2)

			! Apply the linear detrend both the arrays.
			CALL detrend(ecv,ss,hs)
			CALL detrend(ecv,ss,ws)

			! Calculate block standard deviations and averages without calling the respective functions
			! to reduce the computational burden (i.e. in a single block loop).
			wbar=0
			hbar=0
			sdw=0
			sdh=0
			covhw=0
			nnans=0
			DO j=1,ss
				IF(w(j).NE.ecv%nanreal) THEN
					wbar=wbar+ws(j)
					hbar=hbar+hs(j)
					sdw=sdw+ws(j)**2
					sdh=sdh+hs(j)**2
					covhw=covhw+hs(j)*ws(j)
				ELSE
					nnans=nnans+1
				END IF
			END DO
			! Normalize:
			wbar=wbar/(1.0*(ss-nnans))
			hbar=hbar/(1.0*(ss-nnans))
			sdw=SQRT(sdw/(1.0*(ss-nnans))-wbar*wbar)
			sdh=SQRT(sdh/(1.0*(ss-nnans))-hbar*hbar)
			covhw=covhw/(1.0*(ss-nnans))-hbar*wbar

			! Set lag to zero initially; will remain zero if no new cross-correlation maximum is found or if no
			! cross-correlation analysis is applied in the case of undefined wind direction.
			lagli=0
			IF(ABS(dp).LT.madp.OR.wd.EQ.ecv%nanreal) THEN
				lagli=0				! No cross-correlation analysis if below threshold.

			ELSE					! Otherwise perform the analysis.
			
				tlag=CEILING((ecv%dsep/(ute*ABS(dp))))	! Threshold: Maximum number of seconds to perturb the series. Defined as twice the
								! time it takes for the component of the mean wind in the separation direction to travel
								! between the two instruments. Maximum is ~4 seconds.
				tlag=2*tlag*ecv%fs		! --> Convert to number of points in time. Maximum is ~80 points at 20 Hz sampling 									! frequency.

				! First calculate the correlation coefficient without lag.
				mR=covhw/(sdw*sdh)
				mR=ABS(mR) ! <-- Convert to absolute correlation coefficient.
				lagli=0	! Set to zero lag initially.

				! Now we use the SIGN of the dot product to determine which of the series to lag with respect to the other. If
				! dp>0 then the wind is blowing from the CSAT to the LI, so the h series must be lagged wrt. the w series. Vice
				! versa if dp<0.
				IF(dp.GT.0) THEN
					DO k=1,tlag	
						! Reset cross-covariance for each pass.
						ccovhw=0
						DO j=1,ss-k
							IF(hs(j+k).NE.ecv%nanreal.AND.ws(j).NE.ecv%nanreal) THEN
								ccovhw=ccovhw+(hs(j+k)-hbar)*(ws(j)-wbar)
							END IF
						END DO
						ccovhw=ccovhw/(1.0*(ss-nnans))
						Rte=ABS(ccovhw/(sdw*sdh))
						IF(Rte.GT.mR) THEN	! <-- If this exceeds the maximum correlation (so far) then update
									! the maximum and the 'lagli' index.
							lagli=k		! Set the lag to positive in that in this case the h series is considered later 
									! in time than the w series.
							mR=Rte
						END IF
					END DO
				ELSE
					DO k=1,tlag	! Same as above but this time the w series is lagged (h series is advanced) by 1 each pass.
						! Reset cross-covariance for each pass.
						ccovhw=0
						DO j=1,ss-k
							IF(hs(j).NE.ecv%nanreal.AND.ws(j+k).NE.ecv%nanreal) THEN
								ccovhw=ccovhw+(hs(j)-hbar)*(ws(j+k)-wbar)
							END IF
						END DO
						ccovhw=ccovhw/(1.0*(ss-nnans))
						Rte=ABS(ccovhw/(sdw*sdh))	
						IF(Rte.GT.mR) THEN
							lagli=-1*k
							mR=Rte
						END IF
					END DO
				END IF
			END IF

			! Finally compute the water vapor flux based on the lagtime determined by the cross-correlation analysis.
			! If the maximum occurs at lag index > tlag/1.25 the maximum is not well defined so keep the zero lag
			! covariance estimate.
			IF(lagli.NE.0) THEN
				IF(lagli.LT.0.AND.ABS(lagli).LT.NINT((1.0*tlag)/1.25)) THEN
					k=ABS(lagli)
					DO j=1,ss-k
						IF(h(j).NE.ecv%nanreal.AND.w(j+k).NE.ecv%nanreal) THEN
							hw=hw+(hs(j)-hbar)*(ws(j+k)-wbar)
						END IF
					END DO
					hw=hw/(1.0*(ss-nnans))
					!CALL covariance(ecv,ss+lagli,h(1:ss+lagli),w(1-lagli:ss),hw)
					!hw=hw*((1.0*ss)/(1.0*ss+lagli)) 	! Again correct the normalization.
				ELSEIF(lagli.GT.0.AND.ABS(lagli).LT.NINT((1.0*tlag)/1.25)) THEN
					k=ABS(lagli)
					DO j=1,ss-k
						IF(h(j+k).NE.ecv%nanreal.AND.w(j).NE.ecv%nanreal) THEN
							hw=hw+(hs(j+k)-hbar)*(ws(j)-wbar)
						END IF
					END DO
					hw=hw/(1.0*(ss-nnans))
					!CALL covariance(ecv,ss-lagli,h(lagli+1:ss),w(1:ss-lagli),hw)
					!hw=hw*((1.0*ss)/(1.0*ss-lagli)) 	! Again correct the normalization.
				ELSE
					lagli=0
					hw=covhw					
				END IF
			ELSE
				lagli=0
				hw=covhw
			END IF
			IF(u.EQ.0.OR.ISNAN(u).EQV..TRUE..OR.ISNAN(hw).EQV..TRUE.) THEN
				PRINT*,'ubar is zero or nan or hw is nan in separation'
				WRITE(ecv%loglun,*)'ubar is zero or nan or hw is nan in separation'	
			END IF		

			! Calculate the uncorrected value of the flux:
			old=covhw	
			! Calculate the flux correction factor:
			CFsep=hw/old

			! Check that the absolute value of the flux has increased, should always be the case and
			! the increase should be at most on the order of 20% (otherwise something has gone wrong in the calculation
			! e.g. small number division and problems with floating point arithmetic).:
			IF(CFsep.LT.1.OR.CFsep.GT.1.5) THEN
				hw=old
				CFsep=1
			END IF

			! Finally apply the correction.
			CALL covariance(ecv,ss,h,w,hw)
			hw=CFsep*hw

			! Deallocate temp arrays.
			DEALLOCATE(tempa)
			DEALLOCATE(ws)
			DEALLOCATE(hs)
				
			RETURN
		END SUBROUTINE

			
		


		! Simplified Flux attenuation corrections following Massman (2000,2001) henceforth M00 and M01 based on analytical solutions
		! to the integrated normalized Kaimal cospectra (see e.g. Kaimal & Finnigan 1994) multiplied by  
		! transfer functions for high and low frequency attenuation. INPUT: block averaged wind, mes height, MO-parameter and
		! block vertical eddy covariances for: longitudinal momentum, lateral momentum, sonic temperature and absolute humidity.
		! NB! Constants determined by instruments must be changed depending on setup.
		! Further all lengths are in meters, times in seconds and natural frequency in Hz.
		SUBROUTINE attenuation(ecv,ss,z,ubar,zL,&
					uw,vw,tsw,hw,famom,fasen,falat)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: ss			! Length of blocks on which input is based.
			REAL*8,INTENT(IN)		:: z,ubar		! Input mes. height and block average horizontal wind.
			REAL*8,INTENT(IN)		:: zL			! Obukhov length based on sonic temperature.
										! & corrected ustar.
			REAL*8,INTENT(INOUT)		:: uw,vw,tsw,hw 	! Input fluxes -> output corrected fluxes.
			REAL*8,INTENT(OUT)		:: famom,fasen,falat	! Flux attenuation factors as output (momentum, sensible heat and 
										! latent heat flux).
			REAL*8				:: alph			! Cospectral broadness factor.
			REAL*8				:: fxxs,fxxm		! Natural frequencies of maxima in cospectra.
			REAL*8				:: nxxs,nxxm		! Non-dimensional frequencies of maxima in cospectra.
			REAL*8				:: b,p			! Coefficients.
			REAL*8				:: tb			! Block averaging period.
			REAL*8				:: taumf,taubf,tauhf		! Combined equivalent first order time constants for the entire set
											! of LP-filters for momentum flux, buoyancy flux (sonic) and
											! vapor flux (IRGA & Sonic).
			REAL*8				:: tli,tsmh,tsmv,tss,tba 	! Equivalent first order time constants of first order filters
											! for LI-line averaging and sonic line averaging  
											! of u,w and Ts, as well as block averaging. [Table 1 M00].
			REAL*8				:: lli=0.125,lsv=0.1,lsh=0.058	! Path lengths for line averaging of the LI7500 and CSAT3 (vert
											! and horizontal).
			REAL*8				:: pi
			REAL*8				:: ute,zte,zLte			! Temporary reals to hold ubar, z and z/L.

			pi=4.0*ATAN(1.0)
			
			

			! Block averaging period in seconds:
			tb=(1.0*ss)/ecv%fs		

			! Set temporary reals to zero initially
			ute=0
			zte=0
			zLte=0

			! Set the equivalent first order time constants.

			IF(ubar.EQ.0.OR.ISNAN(ubar).EQV..TRUE.) THEN
				PRINT*,'ubar is zero or nan in atten'
				WRITE(ecv%loglun,*)'ubar is zero or nan in atten'
				ute=0.2
				tba=tb/2.8
				tli=lli/(4.0*ute)
				tsmh=lsh/(2.8*ute)
				tsmv=lsv/(5.7*ute)
				tss=lsv/(6.9*ute)
			ELSE
				ute=ubar
				tba=tb/2.8
				tli=lli/(4.0*ute)
				tsmh=lsh/(2.8*ute)
				tsmv=lsv/(5.7*ute)
				tss=lsv/(6.9*ute)
			END IF

			! Based on these calculate the combined equivalent first order time constants for the entire set of low pass transfer funcs
			! for each flux based on EQ.9 in M00 and the applicable LOW PASS transfer functions in each case.
			taumf=SQRT(tsmh**2+tsmv**2)
			taubf=tss
			tauhf=SQRT(tsmv**2+tli**2)

			! Check that z is not <1 or nan.
			IF(z.EQ.ecv%nanreal.OR.z.LT.1) THEN
				zte=1
			ELSE
				zte=z
			END IF


			! Impose bounds on zL in the stable case; not allowing a value over 5 to keep the attenuation
				
		
			! Calculate the maximum nd frequency based on the parametrization in M00 (in turn based on Kaimal & Finnigan 1994):
			IF(zL.LE.0) THEN ! Unstable case.
				nxxm=0.079
				nxxs=0.085
				alph=0.925	! Broadness factor.
			ELSE	! Stable case.
				IF(zL.GT.5) THEN
					zLte=5
				ELSE
					zLte=zL
				END IF
				nxxm=0.079*(1+7.9*zLte)**(3.0/4.0)
				nxxs=2.0-1.915/(1.0+0.5*zLte)
				alph=1.0		! Broadness factor.
			END IF
			! Convert to natural (dimensional) frequency in Hz:
			fxxs=nxxs*ute/z
			fxxm=nxxm*ute/z


			! Now correct the fluxes; calculate the flux attenuation factors first for each flux.
					
			! Momentum:
			! Constant associated with the high pass (boxcar) transfer function:
			b=2*pi*fxxm*tba	
			! Coefficient p changes in accordance with the equivalent low pass time constant:
			p=2*pi*fxxm*taumf

			! Apply equations (6) in M01 to compute the flux attenuation [NB. No recursive filtering is done so a^\alpha->infty].
			famom=((b**alph)/(b**alph+1))*((b**alph)/(b**alph+p**alph))*(1.0/(p**alph+1))

			! Sonic heat flux:
			b=2*pi*fxxs*tba	
			p=2*pi*fxxs*taubf

			! Apply equations (6) in M01 to compute the flux attenuation [NB. No recursive filtering is done so a^\alpha->infty].
			fasen=((b**alph)/(b**alph+1))*((b**alph)/(b**alph+p**alph))*(1.0/(p**alph+1))

		
			! Vapor flux:
			b=2*pi*fxxs*tba	
			p=2*pi*fxxs*tauhf

			! Apply equations (6) in M01 to compute the flux attenuation [NB. No recursive filtering is done so a^\alpha->infty].
			falat=((b**alph)/(b**alph+1))*((b**alph)/(b**alph+p**alph))*(1.0/(p**alph+1))

			! Finally correct the kinematic fluxes using the flux attenuation factors.
			IF(fasen.GT.1.OR.fasen.LT.0.5) THEN	! <-- Should be bounded as a correction between 1 and 1.5.
				fasen=1
				tsw=tsw
			ELSE
				tsw=tsw/fasen
			END IF
			IF(falat.GT.1.OR.falat.LT.0.5) THEN	! -"-.
				falat=1
				hw=hw
			ELSE
				hw=hw/falat
			END IF
			IF(famom.GT.1.OR.famom.LT.0.5) THEN	!-"-
				famom=1
				uw=uw
				vw=vw
			ELSE
				uw=uw/famom
				vw=vw/famom
			END IF
			

			RETURN
		END SUBROUTINE				
			
			
			
			
		! The SND correction after Schotanus et al. (1983), converting buoyancy flux to sensible heat flux.
		! NOTE: No cross-wind correction is required as this is done internally for the CSAT3.			
		!SUBROUTINE sndcorrection(ecv,Tsw,Tbar,rhobar,hw,Twout)
		SUBROUTINE sndcorrection(ecv,Tsw,hw,qw,Tsbar,Tbar,Pbar,hbar,Twout,CFSND,CFSND0)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			REAL*8,INTENT(IN)		:: Tsw,hw,qw,Tsbar,Tbar, Pbar,hbar ! Input block sonic heat flux, water vapor flux, spec h flux,
										! and block average: sonic temp, temp, pressure and absolute humidity. 
	
			REAL*8,INTENT(OUT)		:: Twout		! Output (kinematic) sensible heat flux corrected for sonic
										! temperature.
			REAL*8,INTENT(OUT)	:: CFSND,CFSND0	! Correction factor (SND) for our version and the usual version.

			REAL*8				:: a=0.51		! Constant.
			REAL*8				:: Rd=287.04 ! Gas constant for dry air [JK^-1kg^-1].
			REAL*8				:: old				! Temporary real used in the calculation of the correction factor.

			! Set the initial value of the sensible heat flux as that of the sonic heat flux.
			old=Tsw

			! Apply the SND correction.
			IF(pbar.EQ.0.OR.ISNAN(pbar).EQV..TRUE..OR.ISNAN(Tbar).EQV..TRUE..OR.ISNAN(hbar).EQV..TRUE.) THEN
				Twout=Tsw
				PRINT*,'NaN or zero dev issue in SND'
				WRITE(ecv%loglun,*)'NaN or zero dev issue in SND'
			ELSE
				Twout=Tsw*( 1.0- (a*Rd*hbar*Tbar/Pbar) ) - a*Rd*Tsbar*Tbar*hw/Pbar
			END IF

			! Calculate the correction factor for our version of the SND correction.
			CFSND=Twout/old
			
			! And for the usual version:
			CFSND0=(old-a*Tbar*qw)/old

			! Check that it is bounded; can be a problem when kinematic flux values are low.
			IF(CFSND.LT.0.5.OR.CFSND.GT.2) THEN
				CFSND=1
				CFSND0=1
				Twout=CFSND*old
			END IF

			
			
			RETURN
		END SUBROUTINE

		! The WPL correction after Webb et al. (1980); correcting the latent heat flux for the
		! effect of density fluctuations.
		SUBROUTINE wplcorrection(ecv,hw,rhobar,hbar,Tbar,Tw,Eout,CFWPL)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			REAL*8,INTENT(IN)		:: rhobar,hbar,Tbar,Tw	! Input block average density,
											! average absolute humidity, average absolute temperature and
											! kinematic sensible heat flux (post SND).
			REAL*8,INTENT(IN)		:: hw				! Input water vapor flux corrected.
			REAL*8,INTENT(OUT)		:: Eout				! Corrected kinematic water vapor flux.
			REAL*8,INTENT(OUT)		:: CFWPL			! Correction factor (WPL).
			REAL*8				:: mv,md			! Molar mass of water vapor and dry air (g/mol).
			REAL*8				:: mu,sigma			! Ratio of mv/md and block av. vapor and dry air density.
			REAL*8				:: old					! Uncorrected value.	
			 
			! Allocate the molar masses of dry air and water vapor as well as their ratio.
			mv=18.01528	! g/mol
			md=28.97	! g/mol
			mu=md/mv	! -

			! Allocate the ratio of vapor to dry air density (no liquid water assumed).
			sigma=hbar/(rhobar-hbar)

			! Allocate the uncorrected value.
			old=hw

			! Apply the WPL correction.
			IF(Tbar.NE.0.AND.ISNAN(Tbar).EQV..FALSE.) THEN
				Eout=(1+mu*sigma)*(hw+(hbar/Tbar)*Tw)
			ELSEIF(Tbar.EQ.0) THEN
				PRINT*,'Tbar iz zero in WPL'
				WRITE(ecv%loglun,*)'Tbar iz zero in WPL'
				Eout=hw
			ELSEIF(ISNAN(Tbar).EQV..TRUE.) THEN
				PRINT*,'Tbar is NaN in WPL'
				WRITE(ecv%loglun,*)'Tbar is NaN zero in WPl'
				Eout=hw
			END IF		
			
			! Calculate the correction factor.
			CFWPL=Eout/old

			! Check that it is bounded; can be a problem when kinematic flux values are low.
			IF(CFWPL.LT.0.5.OR.CFWPL.GT.4) THEN
				CFWPL=1
				Eout=CFWPL*old
			END IF
			
			RETURN
		END SUBROUTINE

		! Deals with all the flux corrections in one subroutine. Applying the last three corrections, i.e.
		! flux attenuation, SND and WPL iteratively until the fluxes differ by less than 0.1% from one iteration
		! to the next. Alot of arguments but this subroutine only needs to be called once.
		SUBROUTINE fluxcorrections(ecv,ss,u,v,w,Ts,h,q,rhobar,Tbar,pbar,zbar,wd,uw,vw,Tsw,Tw,hw,Eout,&
						lagli,famom,fasen,falat,ustar,ol,tau,niter,CFsep,CFSND,CFSND0,CFWPL)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: ss			! Length of a block.
			REAL*8,DIMENSION(ss),INTENT(IN)	:: u,v,w,Ts,h,q		! Input block series of high frequency variables.
			REAL*8,INTENT(IN)		:: rhobar,Tbar,pbar	! Input block average air density and absolute temperature.
			REAL*8,INTENT(IN)		:: wd,zbar		! Input block average (horizontal) wind direction and measurement height
			REAL*8,INTENT(OUT)		:: uw,vw,hw,Tsw,Tw	! Output corrected kinematic fluxes.
			REAL*8,INTENT(OUT)		:: Eout			! Output WPL corrected vapor flux.
			REAL*8,INTENT(OUT)		:: ol,ustar,tau		! Output corrected obukhov length, friction velocity and
										! Reynolds stress. 
			REAL*8,INTENT(OUT)		:: famom,fasen,falat	! Output flux attenuation factors.
			INTEGER,INTENT(OUT)		:: lagli		! Output block sample lag for the IRGA (>0: LI is lagged wrt. 
										! CSAT and <0: LI is advanced wrt to CSAT).
			REAL*8,INTENT(OUT)		:: niter		! Output number of iterations before flux corrections converge.
			REAL*8,INTENT(OUT)		:: CFsep,CFSND,CFSND0,CFWPL	! Correction factors.
			
			! Internal declarations:
			INTEGER				:: maxit=1e2		! Maximum number of iterations.
			INTEGER				:: i,j,k		! Counters.
			REAL*8				:: tol=1e-3		! Tolerance for convergence of iterations from one pass to the next.
			REAL*8				:: zLt			! MO-Stability parameter based on sonic temperature for the flux 
										! attenuation using the Kaimal flat-terrain cospectra.
			REAL*8				:: Tsbar,vptbar,hbar ! Block averaged Ts,T,VPT,AH and pressure.
			REAL*8				:: ubar,vbar		! Block averaged hor. wind speeds. 
			REAL*8		 		:: kappa=0.4,g=9.81	! Von-Karman constant and acceleration of gravity for zL.
			REAL*8				:: uwold,vwold,Tswold,Twold,hwold,qwold,Eold	! Kinematic fluxes from last iteration.
			REAL*8				:: uwnew,vwnew,Tswnew,Twnew,hwnew,qwnew,Enew	! Kinematic fluxes from newest iteration.
			REAL*8				:: uw0,vw0,Tsw0,Tw0,hw0
			REAL*8				:: duw,dvw,dTsw,dTw,dhw	,dE		! Normalized absolute difference of fluxes from one 													! iteration to the next.
			INTEGER				:: ok					! Status integer to break out of iterations.
			REAL*8				:: tmp					! Temporary real.

			
			! Zeroing.
			uw0=0
			vw0=0
			Tsw0=0
			Tw0=0
			hw0=0
			qwold=0

			! Calculate the kinematic buoyancy fluxes.
			CALL covariance(ecv,ss,Ts,w,Tsw0)

			! And specific humidity flux
			CALL covariance(ecv,ss,q,w,qwold)

			! And kinematic momentum fluxes.
			CALL covariance(ecv,ss,u,w,uw0)
			CALL covariance(ecv,ss,v,w,vw0)
			
			! And block averages.
			CALL temporalmean(ecv,ss,Ts,Tsbar)
			CALL temporalmean(ecv,ss,h,hbar)
			CALL temporalmean(ecv,ss,v,vbar)	! <--- In case a streamwise rotation has not been applied prior to the call on the flux corrections.
			CALL temporalmean(ecv,ss,u,ubar)
			ubar=SQRT(ubar**2+vbar**2)

			! Calculate the MO-stability parameter based on the above.
			ol=-1.0*((SQRT(SQRT(uw**2+vw**2)))**3)/((kappa*g/Tsbar)*Tsw)
			zLt=zbar/ol
		
			! Calculate the kinematic water vapor flux taking sensor separation into account.
			CALL separation(ecv,ss,ubar,wd,h,w,hw0,lagli,CFsep)

			! Start the iterative flux corrections.
			ok=0
			uwold=uw0
			vwold=vw0
			Tswold=Tsw0
			Twold=Tsw0
			hwold=hw0
			k=0	! Iteration counter.
			DO WHILE(k.LT.maxit.AND.ok.EQ.0)	! Start iterations.
				uwnew=uw0
				vwnew=vw0
				Tswnew=Tsw0
				Twnew=Tsw0
				hwnew=hw0

				! 1) Flux attenuation correction with MO-Stability parameter from last pass.
				CALL attenuation(ecv,ss,zbar,ubar,zLt,uwnew,vwnew,Tswnew,hwnew,famom,fasen,falat)
				
				! Also correct the specific humidity flux [Only used to compare our SND method with the usual one].
				qwnew=qwold*CFsep
				qwnew=qwnew/falat
		
				! 2) SND correction with flux attenuation corrections applied.	
				CALL sndcorrection(ecv,Tswnew,hwnew,qwnew,Tsbar,Tbar,Pbar,hbar,Twnew,CFSND,CFSND0)
	
				! 3) WPL correction iwth flux attenuation & SND correction applied.
				CALL wplcorrection(ecv,hwnew,rhobar,hbar,Tbar,Twnew,Enew,CFWPL)

				! Update the Obukhov length and MO-stability parameter based on the above for next pass (or output).
				ol=-1.0*((SQRT(SQRT(uwnew**2+vwnew**2)))**3)/((kappa*g/Tsbar)*Tswnew)
				zLt=zbar/ol	
				IF(zbar.EQ.0) THEN
					PRINT*,'zbar is zero in flux corrections!'
					WRITE(ecv%loglun,*)'zbar is zero in flux corrections!'
				ELSEIF(ISNAN(zLt).EQV..TRUE.) THEN
					PRINT*,'zL is nan in flux corrections!'
					WRITE(ecv%loglun,*)'zL is nan in flux corrections'			
				END IF				
	
				! Compute the normalized absolute deviations in the corrections compared to the last pass.
				IF(uwold.NE.0) THEN
					duw=ABS(1e3*(uwnew-uwold)/(1e3*uwold))
				ELSE
					duw=0
				END IF
				IF(vwold.NE.0) THEN
					dvw=ABS(1e3*(vwnew-vwold)/(1e3*vwold))
				ELSE
					dvw=0
				END IF
				IF(Tswold.NE.0) THEN
					dTsw=ABS(1e3*(Tswnew-Tswold)/(1e3*Tswold))
				ELSE
					dTsw=0
				END IF
				IF(Twold.NE.0) THEN
					dTw=ABS(1e3*(Twnew-Twold)/(1e3*Twold))
				ELSE
					dTw=0
				END IF
				IF(hwold.NE.0) THEN
					dhw=ABS(1e3*(hwnew-hwold)/(1e3*hwold))
				ELSE
					dhw=0	
				END IF
				IF(Eold.NE.0) THEN
					dE=ABS(1e3*(Enew-Eold)/(1e3*Eold))
				ELSE
					dE=0
				END IF

				IF(duw.LE.tol.AND.dvw.LE.tol.AND.dTsw.LE.tol&
						.AND.dTw.LE.tol.AND.dhw.LE.tol.AND.dE.LE.tol) THEN ! Break out of loop if all fluxes converge.  					
					k=k+1
					ok=1
				ELSE	! Otherwise swap new with old for next pass.
					uwold=uwnew
					vwold=vwnew
					Tswold=Tswnew
					Twold=Twnew
					hwold=hwnew
					Eold=Enew
					k=k+1	! Update iteration counter.
				END IF


			END DO
			! Set the number of iterations for convergence.
			niter=1.0*k

			! Update fluxes to the corrected values after convergence.
			uw=uwnew
			vw=vwnew
			Tw=Twnew
			Tsw=Tswnew
			hw=hwnew
			Eout=Enew
			ustar=SQRT(SQRT(uw**2+vw**2))
			tau=rhobar*(ustar**2)

			ol=zbar/zLt
			

			RETURN
		END SUBROUTINE
			
			
		

		! Calculates the friction velocity (ustar), given arrays of u,v,w of EQUAL length. The length (ss)
		! corresponds to the averaging period desired. If ss is set to 1 then uv and vw must be specified;
		! this is the case for the 'flux corrections' iterations.
		SUBROUTINE frictionvelocity(ecv,ss,u,v,w,ustar)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: ss	! Length of input arrays.
			REAL*8,DIMENSION(ss),INTENT(IN)	:: u,v,w! Input velocity arrays.
			REAL*8,INTENT(OUT)		:: ustar! Friction velocity: output.
			REAL*8				:: uw,vw! Covariances.
			
			ustar=0
			uw=0
			vw=0			
			
			CALL covariance(ecv,ss,u,w,uw)
			CALL covariance(ecv,ss,v,w,vw)
				
			! Compute friction velocity:
			ustar=SQRT(SQRT(uw**2+vw**2))

			RETURN
		END SUBROUTINE

		! Calculates the Obukhov length (L) given arrays of vpt,u,v and w of EQUAL length. The length (ss)
		! corresponds to the averaging period desired. 
		SUBROUTINE obukhovlength(ecv,ss,u,v,w,ts,ol)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: ss		! Length of input arrays.
			REAL*8,DIMENSION(ss),INTENT(IN)	:: u,v,w,ts	! Input arrays.
			REAL*8,INTENT(OUT)		:: ol		! Obukhov length: output. 
			REAL*8				:: ustar	! Friction velocity.
			REAL*8				:: tsw,tsbar	! vpt and w covariance and mean vpt.
			REAL*8,PARAMETER		:: kappa=0.4	! Von Karman constant.
			REAL*8,PARAMETER		:: g=9.81	! Acceleration of gravity.
			
			ustar=0
			ol=0
			
	
			CALL frictionvelocity(ecv,ss,u,v,w,ustar)
			CALL temporalmean(ecv,ss,ts,tsbar)
			CALL covariance(ecv,ss,ts,w,tsw)		
			
			! Compute the Obukhov length:
			ol=-1*(ustar**3)/(kappa*g*tsw/tsbar)

			RETURN
		END SUBROUTINE

		!----------------------------------------------------------------------------------------------------------------------
		! Statistics.
		!----------------------------------------------------------------------------------------------------------------------



		! Temporal mean, ignoring nan entries.
		SUBROUTINE	temporalmean(ecv,mm,a,mean) 
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: mm		! Segment length and weight.
			INTEGER				:: mmm		! Effective segment length (discounting nans).
			REAL*8,DIMENSION(mm),INTENT(IN)	:: a		! Array segment.
			REAL*8,INTENT(OUT)		:: mean		! Mean.
			INTEGER				:: i 		! Counter.
			
			mean=0 ! Initialize the mean to be zero (for multiple calls).
			mmm=mm
			DO i=1,mm
				IF(a(i).EQ.ecv%nanreal) THEN	! Ignore nan entries.
					mmm=mmm-1
				ELSE
					mean=mean+a(i)
				END IF
			END DO
			IF(mmm.NE.0) THEN
				mean=mean/mmm
			ELSE
				mean=ecv%nanreal
			END IF

			RETURN
		END SUBROUTINE	
	
		! Biased standard deviation, also ignores nan entries.
		SUBROUTINE	standarddeviation(ecv,mm,a,sd)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER,INTENT(IN)		:: mm		! Segment length		
			INTEGER				:: mmm		! Effective segment length (discounting nans).
			REAL*8,DIMENSION(mm),INTENT(IN)	:: a		! Array segment.
			REAL*8,INTENT(OUT)		:: sd		! Standard deviation.
			REAL*8				:: amean	! Mean.
			INTEGER				:: i 		! Counter.
			
			CALL temporalmean(ecv,mm,a,amean)
			CALL covariance(ecv,mm,a,a,sd)
			IF(sd.NE.ecv%nanreal) THEN
				sd=SQRT(sd)
			END IF
			
			RETURN
		END SUBROUTINE

		! Biased covariance between two series.
		SUBROUTINE covariance(ecv,mm,a,b,cov)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER, INTENT(IN)		:: mm		! Segment length.
			INTEGER				:: mmm		! Effective segment length.
			REAL*8, DIMENSION(mm),INTENT(IN):: a,b		! Array segments.
			REAL*8,INTENT(OUT)		:: cov		! Covariance between the array segments.
			REAL*8				:: amean,bmean	! Means.
			INTEGER				:: i 		! Counter.
			
			
			CALL temporalmean(ecv,mm,a,amean) ! Compute the mean of the a series/window.
			CALL temporalmean(ecv,mm,b,bmean) ! Compute the mean of the b series/window.
			
			cov=0	! Intialize the covariance to be zero (for multiple calls).
			mmm=mm
			IF(amean.NE.ecv%nanreal.AND.bmean.NE.ecv%nanreal) THEN
				DO i=1,mm
					IF(a(i).EQ.ecv%nanreal.OR.b(i).EQ.ecv%nanreal) THEN
						mmm=mmm-1
					ELSE
						cov=cov+(a(i)-amean)*(b(i)-bmean)
					END IF
				END DO
			ELSE
				cov=ecv%nanreal
			END IF
			IF(mmm.NE.0.AND.amean.NE.ecv%nanreal.AND.bmean.NE.ecv%nanreal) THEN
				cov=cov/mmm
			ELSE
				cov=ecv%nanreal
			END IF

			RETURN
		END SUBROUTINE

		! Correlation coefficient between two series.
		SUBROUTINE correlationcoefficient(ecv,mm,a,b,corrc)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv		
			INTEGER, INTENT(IN)		:: mm 		! Segment lengths.
			REAL*8, DIMENSION(mm),INTENT(IN):: a,b		! Array segments.
			REAL*8, INTENT(OUT)		:: corrc	! Correlation between the array segments.
			REAL*8				:: abcov,asd,bsd! Covariances and standard deviations.
			INTEGER				:: i 		! Counter.
			

			CALL standarddeviation(ecv,mm,a,asd)
			CALL standarddeviation(ecv,mm,b,bsd)
			CALL covariance(ecv,mm,a,b,abcov)
			
			! Calculate the linear correlation coefficient.
			IF(abcov.EQ.ecv%nanreal.OR.asd.EQ.ecv%nanreal.OR.bsd.EQ.ecv%nanreal) THEN
				corrc=ecv%nanreal
			ELSE
				corrc=abcov/(asd*bsd)
			END IF

			RETURN
		END SUBROUTINE

		! A routine for applying a simple linear detrend to segments of data. Following Gash & Culf 1996.
		! Note that this is NOT applied to the measured variables themselves, but only in the calculation
		! of higher moments (m>2) such as kurtosis and skewness which assume a random distribution where
		! we must at the very least remove any linear trend. Also included in spectra and structure functions
		! (optional).
		SUBROUTINE detrend(ecv,mm,a)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv
			INTEGER, INTENT(IN)			:: mm		! Segment length.
			INTEGER					:: mmm		! Effective segment length.
			REAL*8, DIMENSION(mm),INTENT(INOUT)	:: a		! Array segment.
			REAL*8, DIMENSION(mm)			:: ade 		! Skewness of array segment.
			REAL*8					:: beta,c,d	! Linear regression coefficients.
			REAL*8					:: e,f,g,h	! Placeholder reals.
			INTEGER					:: i,j		! Counters.
			
					


			! Assign the effective segment length.			
			mmm=mm
			! Set placeholder values to zero initially.
			e=0
			f=0
			g=0
			h=0
			ade(:)=0
			j=0
			! Evaluate the components determining the linear regression coefficients.
			DO i=1,mm ! Position in series.
				j=j+1 ! Update the actual counter.
				IF(a(i).EQ.ecv%nanreal) THEN	! Exclude nan entries.
					mmm=mmm-1
					j=j-1 ! Nan entries are not counted.
				ELSE ! We first use the c,d coefficients as placeholders.
					e=e+j**2
					f=f+j
					g=g+a(i)
					h=h+j*a(i)
				END IF
			END DO
			
			! Compute the linear regression coefficients.
			beta=mmm*e-f**2 
			d=(mmm*h-f*g)/beta
			c=(g-d*f)/mmm
			
			j=0	! <-- Counter (excluding nan entries).
			! Calculate the linear best fit line for the array.
			DO i=1,mm ! Position in series.
				j=j+1 ! Update the actual counter.
				IF(a(i).EQ.ecv%nanreal) THEN ! Exclude nan entries.
					j=j-1 ! Nan entries are not counted.
					ade(i)=0
				ELSE 
					ade(i)=c+d*j
				END IF
			END DO

			! Finally remove the linear trend from the series.
			ade=a-ade

			! Change the input/output array to the detrended version.
			a=ade
			
			! Reset values to make sure no errors occur on multiple passes.
			ade(:)=0
			e=0
			f=0
			g=0
			h=0
			beta=0
			d=0
			c=0
			j=0
			

			RETURN
		END SUBROUTINE

		! Fills the gaps in the array a with dimensions mm x nn, where mm is the segment length
		! and nn is the variable dimension. Simple gap filling using linear interpolation wherever
		! nans occur. NB! Assumes that all the series contained in 'ay' have synchronized nan entries. 
		SUBROUTINE gapfill(ecv,mm,nn,ay)
			TYPE(eddycov), INTENT(INOUT)		:: ecv
			INTEGER, INTENT(IN)			:: mm		! Segment length.
			INTEGER, INTENT(IN)			:: nn		! Number of variables in the input array.
			REAL*8,DIMENSION(mm,nn),INTENT(INOUT)	:: ay		! Input/output array.
			INTEGER					:: i,j,k,a,b	! Counters.
			REAL*8					:: nanis	! Nan (shorthand).
			
			nanis=ecv%nanreal
			a=0
			b=0
			k=0
			! Check end points of the segment for NaN-occurence.
			IF(ay(1,1).EQ.nanis) THEN	! If first point is nan then:
				k=1
				DO WHILE(ay(k+1,1).EQ.nanis)	! Scan for next non-nan occurence.
					k=k+1
				END DO
				a=k+1 				! Set next non-nan occurence index.
				IF(ay(a+1,1).EQ.nanis) THEN	! Check if next point is nan.
					k=a+1			
					DO WHILE(ay(k+1,1).EQ.nanis)	! If so scan for next non-nan occurence.
						k=k+1			
					END DO
					b=k+1				! Set next non-nan occurence index.
				ELSE				
					b=a+1			! If not set as next non-nan occurence index.
				END IF
				
				DO j=1,a-1			! Perform linear interpolation for all j<a.
					ay(j,:)=ay(b,:)-(ay(b,:)-ay(a,:))*(1.0*b-1.0*j)/(1.0*b-1.0*a)
				END DO
			END IF
			IF(ay(mm,1).EQ.nanis) THEN		! If last point is nan then.
				k=mm
				DO WHILE(ay(k-1,1).EQ.nanis)	! Scan for previous non-nan occurence.
					k=k-1
				END DO
				b=k-1 				! Set previous non-nan occurence index.
				IF(ay(b-1,1).EQ.nanis) THEN	! Check if previous point is nan.
					k=b-1			
					DO WHILE(ay(k-1,1).EQ.nanis)	! If so scan for previous non-nan occurence.
						k=k-1			
					END DO
					a=k-1				! Set previous non-nan occurence index.
				ELSE				
					a=b-1			! If not set as previous non-nan occurence index.
				END IF
				
				DO j=b+1,mm			! Perform linear interpolation for all j>b.
					ay(j,:)=ay(a,:)+(ay(b,:)-ay(a,:))*(1.0*j-1.0*a)/(1.0*b-1.0*a)
				END DO
			END IF
		
			DO i=2,mm-1
				a=0	! Position of previous non-nan occurence.
				b=0	! Position of next non-nan occurence.
				k=i	
				IF(ay(i,1).EQ.nanis) THEN	! Check for nan-occurence in the first variable at point i.
					a=i-1			! <-- Set last non-nan occurence.
					DO WHILE(ay(k+1,1).EQ.nanis)
						k=k+1
					END DO
					b=k+1			! <-- Set next non-nan occurence.
					DO j=a+1,b-1		! Loop through NaNs and interpolate. Vectorized.
						ay(j,:)=ay(a,:)+(ay(b,:)-ay(a,:))*(1.0*j-1.0*a)/(1.0*b-1.0*a)
					END DO
				END IF
			END DO
			a=0
			b=0
			k=0
			RETURN
		END SUBROUTINE
						
					
				
			
					
		! Flag from Vickers and Mahrt (1997).
		! Sample skewness, third order moment normalized by cubed standard deviation.
		! NB! Operates on a linearly detrended series to attempt to assure a random distribution.
		! Positive: PDF peak skewed to the left (\).
		! Negative: PDF peak skewed to the right (/).
		SUBROUTINE skewness(ecv,mm,a,skew)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv	
			INTEGER, INTENT(IN)			:: mm		! Segment length.
			INTEGER					:: mmm		! Effective segment length.
			REAL*8, DIMENSION(mm),INTENT(IN)	:: a		! Array segment to be detrended.
			REAL*8, INTENT(OUT)			:: skew 	! Skewness of array segment.
			REAL*8,DIMENSION(mm)			:: atemp	! Detrended version of a.
			REAL*8					:: asd,amean 	! Standard deviation and mean of array segment.
			INTEGER					:: i		! Counter.

			atemp=a
			CALL detrend(ecv,mm,atemp)
			CALL temporalmean(ecv,mm,atemp,amean)
			CALL standarddeviation(ecv,mm,atemp,asd)

			skew=0 ! Intialize as zero (for multiple calls).
			! First calculate the third order moment. (Mean of the cubed fluctuation).
			mmm=mm
			DO i=1,mm
				IF(a(i).EQ.ecv%nanreal) THEN
					mmm=mmm-1
				ELSE
					skew=skew+(atemp(i)-amean)**3
				END IF			
			END DO
			skew=skew/mmm
			
			! Then normalize by the cubed standard deviation to get the skewness.
			skew=skew/(asd**3)
		
			! Note that the third moment is recovered by multiplying the skewness by the cubed sd.
			RETURN
		END SUBROUTINE

		! Flag from Vickers and Mahrt (1997).
		! Sample kurtosis, fourth order moment normalized by square variance.
		! NB! Operates on a linearly detrended series to attempt to assure a random distribution.
		! Higher value => sharp peak & long wings (LaCasce).
		SUBROUTINE kurtosis(ecv,mm,a,kurt)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv	
			INTEGER, INTENT(IN)			:: mm		! Segment length.
			INTEGER					:: mmm		! Effective segment length.
			REAL*8, DIMENSION(mm),INTENT(IN)	:: a		! Array segment to be detrended.
			REAL*8, INTENT(OUT)			:: kurt		! Kurtosis of array segment.
			REAL*8,DIMENSION(mm)			:: atemp	! Detrended version of a.
			REAL*8					:: asd,amean 	! Standard deviation and mean of array segment.
			INTEGER					:: i		! Counter.		

			atemp=a
			CALL detrend(ecv,mm,atemp)
			CALL temporalmean(ecv,mm,atemp,amean)
			CALL standarddeviation(ecv,mm,atemp,asd)
			
			kurt=0	! Intialize as zero (for multiple calls).
			! First calculate the third order moment. (Mean of the cubed fluctuation).
			mmm=mm
			DO i=1,mm
				IF(a(i).EQ.ecv%nanreal) THEN
					mmm=mmm-1
				ELSE
					kurt=kurt+(atemp(i)-amean)**4
				END IF
			END DO
			kurt=kurt/mmm
			
			! Then normalize by the cubed standard deviation to get the skewness.
			kurt=kurt/(asd**4)
		
			! Note that the fourth moment is recovered by multiplying the skewness by the square variance.
			RETURN
		END SUBROUTINE



		! Calculates the stationarity coefficient as in Foken and Wichura (1996).
		! The stationarity coefficient is the difference between the mean of the covariance of
		! subsegments and the covariance of the segment normalized by the covariance of
		! the segment. As it approaches zero the segment approaches stationarity.
		! If a=b the stationarity coefficient is a measure of variance stationarity. Conversely
		! if a=/b and either a or b is the vertical velocity then the stationarity coefficient 
		! is a measure of flux stationarity.
		SUBROUTINE stationarity(ecv,mm,ss,a,b,rstat)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv	
			INTEGER, INTENT(IN)			:: mm			! Segment length.
			INTEGER, INTENT(IN)			:: ss			! Subsegment length, integer factor of the segment length.
			INTEGER					:: nss		 	! Number of subsegments.
			REAL*8, DIMENSION(mm),INTENT(IN)	:: a,b			! Array segment.
			REAL*8, ALLOCATABLE, DIMENSION(:,:)	:: as,bs		! Array subsegments.
			REAL*8, INTENT(OUT)			:: rstat		! Stationarity coefficient of array segment.
			REAL*8					:: abcov		! Covariance of array segment.
			REAL*8, ALLOCATABLE, DIMENSION(:)	:: abscov	 	! Covariance of array subsegments.
			REAL*8					:: mabscov		! Mean of array subsegment covariances.
			INTEGER					:: i,j,k		! Counters.

			! Control the segment length and the subsegment length.
			IF(1.0*mm/ss.NE.NINT(1.0*mm/ss)) THEN
				PRINT*,	'Error: Subsegment length must be an integer factor of the segment length'
				STOP
			ELSE ! Is an integer factor.
				nss=mm/ss ! Allocate the number of subsegments.
			END IF
			!PRINT*,'number of subsegments considered:',nss

			! Allocate the subsegment arrays.
			ALLOCATE(as(ss,nss))	
			ALLOCATE(bs(ss,nss))
			! Allocate the subsegment covariance array.
			ALLOCATE(abscov(nss))
			
			DO k=1,nss ! Loop over columns.
				! In general A_ij=a([j-1]*I+i)=a(k)=a_k where ss=I and nss=J.	
				i=((k-1)*ss+1) 	! Index in segment array of first row entry in the column of the subsegment array.
				j=i+ss-1	! Index in segment array of last row entry in the column of the subsegment array.
				as(1:ss,k)=a(i:j)	! Corresponding entry in subsegment array. 
				bs(1:ss,k)=b(i:j)	! -"-.
			END DO


			CALL covariance(ecv,mm,a,b,abcov)			! Covariance of the segment.
			DO j=1,nss
				CALL covariance(ecv,ss,as(:,j),bs(:,j),abscov(j))! Subsegment covariances.
			END DO
			CALL temporalmean(ecv,nss,abscov,mabscov)		! Mean of subsegment covariances.
			
			IF(mabscov.NE.ecv%nanreal.AND.abcov.NE.ecv%nanreal.AND.abcov.NE.0) THEN
				rstat=ABS((mabscov-abcov)/abcov)		! Stationarity coefficient.
			ELSE
				rstat=ecv%nanreal
			END IF

			i=0
			j=0
			mabscov=0
			abcov=0
			

			RETURN
		END SUBROUTINE

		! Calculates parametrized non-dimensional w variance in given stability ranges. Compares this to
		! the actual non-dimensional w variance and returns a flag based on the deviation of the latter from the former.
		! in the 'integral turbulence characteristics' assesment.
		SUBROUTINE wITC(ecv,ss,z,u,v,w,ts,fw)
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)	:: ecv
			INTEGER,INTENT(OUT)		:: fw			! Output ITC flag (0-2).
			INTEGER				:: ss			! Size of input arrays (should be block size).
			REAL*8,DIMENSION(ss),INTENT(IN)	:: z,u,v,w,ts		! Input arrays.
			REAL*8				:: zbar			! Mean measurement height.
			REAL*8				:: wvar			! Vertical velocity variance.
			REAL*8				:: ol			! Obukhov length.
			REAL*8				:: zoverL		! MO-Stability parameter.
			REAL*8				:: lat			! Latitude (in degrees N).
			REAL*8				:: ustar		! Friction velocity.
			REAL*8				:: nd
			REAL*8				:: nwvard,nwvarp	! Non-dimensional w-variance from the data and parametrization.
			REAL*8				:: pi,fcor,lnz		! For the scaling factor in near-neutral stability.
			REAL*8				:: reldev		! Relative deviation from the parametrization.
			REAL*8,PARAMETER		:: zp=1
			REAL*8,PARAMETER		:: a=0.24,b=2.9		! Coeffs for near-neutral parametrization.
			REAL*8,PARAMETER		:: c=1.9		! Coeff for unstable parametrization.
			REAL*8,PARAMETER		:: tol1=0.3,tol2=0.75	! Flag tolerances.


			lat=1.0*ecv%latitude				! Set the latitude.
			CALL temporalmean(ecv,ss,z,zbar)		! Calculate the mean measurement height.
			CALL frictionvelocity(ecv,ss,u,v,w,ustar)	! Friction velocity.
			CALL obukhovlength(ecv,ss,u,v,w,ts,ol)		! Obukhov length.
			zoverL=zbar/ol					! MO-Stability parameter.
			CALL covariance(ecv,ss,w,w,wvar)		! Vertical velocity variance.
	
			
			fw=0
			IF(zoverL.LT.0.4.AND.zoverL.GT.-3) THEN
				IF(zoverL.GE.-0.2) THEN	! Effectively -0.2<z/L<0.4 the near neutral range.
					! Compute pi.
					pi=4.0*ATAN(1.0)
					! Compute the Coriolis parameter for this latitude.
					fcor=2*(2*pi/(86400.0))*SIN(lat*((pi/2)/90))
					! Calculate the scaling factor for this parametrization.
					lnz=LOG(zp*fcor/ustar)
					! Calculate the parametrized non-dimensional w-variance for this stability
					! range according to Thomas&Foken (2002) with revised coefficients.
					nwvarp=a*lnz+b
				ELSEIF(zoverL.LT.-0.2) THEN ! Effectively -3<z/L<-0.2 the unstable range.
					nwvarp=1.3*((1+c*ABS(zoverL))**(1/3))
				END IF
				! And the actual non-dim wvariance from the data.
				nwvard=SQRT(wvar)/ustar
				! Compute the (absolute) relative deviation from the parametrization.
				reldev=ABS(nwvarp-nwvard)/nwvarp
				! Apply the flag:
				IF(reldev.GT.tol1) THEN
					IF(reldev.GT.tol2) THEN
						fw=2
					ELSE
						fw=1
					END IF
				ELSE
					fw=0
				END IF
			END IF
			! NB. Flag remains zero in the case of stable conditions or very unstable conditions.
			RETURN
		END SUBROUTINE
				
				

		! Calculates the flux intermittency as in Mahrt et al. 1998. Very similar to
		! the stationarity coefficient, defined as the ratio of standard deviation of
		! subsegment mean covariances to the segment mean covariance. As with the
		! stationarity coefficient if a=b we compute the intermittency of variance. Moreover
		! if a=/=b and either a or b is the vertical velocity we compute the flux intermittency.
		SUBROUTINE intermittency(ecv,mm,ss,a,b,rinterm)	
			IMPLICIT NONE
			TYPE(eddycov), INTENT(INOUT)		:: ecv	
			INTEGER, INTENT(IN)			:: mm			! Segment length.
			INTEGER, INTENT(IN)			:: ss			! Subsegment length, integer factor of the segment length.
			INTEGER					:: nss		 	! Number of subsegments.
			REAL*8, DIMENSION(mm),INTENT(IN)	:: a,b			! Array segment.
			REAL*8, ALLOCATABLE, DIMENSION(:,:)	:: as,bs		! Array subsegments.
			REAL*8, INTENT(OUT)			:: rinterm		! Stationarity coefficient of array segment.
			REAL*8					:: abcov		! Covariance of array segment.
			REAL*8, ALLOCATABLE, DIMENSION(:)	:: abscov	 	! Covariance of array subsegments.
			REAL*8					:: sdabscov		! Mean of array subsegment standard deviations.
			INTEGER					:: i,j,k		! Counters.

			! Control the segment length and the subsegment length.
			IF(1.0*mm/ss.NE.NINT(1.0*mm/ss)) THEN
				PRINT*,	'Error: Subsegment length must be an integer factor of the segment length'
				STOP
			ELSE ! Is an integer factor.
				nss=mm/ss ! Allocate the number of subsegments.
				PRINT*,nss
			END IF

			! Allocate the subsegment arrays.
			ALLOCATE(as(ss,nss))	
			ALLOCATE(bs(ss,nss))
			! Allocate the subsegment covariance array.
			ALLOCATE(abscov(nss))
			
			DO k=1,nss ! Loop over columns.
				! In general A_ij=a([j-1]*I+i)=a(k)=a_k where ss=I and nss=J.	
				i=((k-1)*ss+1) 	! Index in segment array of first row entry in the column of the subsegment array.
				j=i+ss-1	! Index in segment array of last row entry in the column of the subsegment array.
				as(1:ss,k)=a(i:j)	! Corresponding entry in subsegment array. 
				bs(1:ss,k)=b(i:j)	! -"-.
			END DO


			CALL covariance(ecv,mm,a,b,abcov)				! Covariance of the segment.
			DO j=1,nss
				CALL covariance(ecv,ss,as(:,j),bs(:,j),abscov(j))	! Subsegment covariances.
				IF(ISNAN(abscov(j))) THEN
					PRINT '("interm isnan at j=",I3)',nss-j
				END IF
			END DO
			CALL standarddeviation(ecv,nss,abscov,sdabscov)		! Standard deviation of subsegment covariances.
			rinterm=(sdabscov)/abcov				! Intermittency coefficient.

			! Deallocate temporary arrays.
			DEALLOCATE(as)
			DEALLOCATE(bs)
			DEALLOCATE(abscov)


			RETURN
		END SUBROUTINE


		SUBROUTINE	autoandstructure(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER					:: mm			! Size of input arrays.
			INTEGER					:: ss			! Size of segments (windows) on which to compute autocorrelation.
			INTEGER					:: ssj			! Number of lags to consider.
			INTEGER					:: nc			! Number of central lag bins (post smoothing).
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: tlag,tlagc		! Lag vectors.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: KOSM			! K-O smoothing matrix.
			REAL*8					:: lag1			! Initial (j=1) lag in time.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: us,vs,ws,Ts,qs	! Arrays split into segments of length ss.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: sobukhovl,ustar	! Obukhov length for the statistical windows.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: szmes		! Measurement height for the statistical windows.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: suh			! Mean horizontal wind for the statistical windows.
			REAL*8					:: tempreal		! Temporary real.
			REAL*8					:: teuw,tevw,tetw	! Temporary reals needed for computing the Obukhov length.
			INTEGER					:: i,j,k,l,ll,ul,di	! Counters.
			INTEGER					:: isnan,ok		! Gap filling counters.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: tempar		! Temporary array for gap filling.
			INTEGER					:: ii,jj		! Output dimensions.
			INTEGER					:: sspmm		! (NINT) segments in the series with length mm.
			REAL*8					:: stime,timeis,lag	! Total and current segment time and current time lag.
			REAL*8,ALLOCATABLE,DIMENSION(:,:,:)	:: twindow		! Timestamp to keep track of windows.
			! Temporary structure functions and autocorrelations containing nan entries:
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: d2uu,d2vv,d2ww,d2qq,d2TT	! <-- 2nd order Structure functions.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: d3uu,d3vv,d3ww,d3qq,d3TT	! <-- 3rd order structure functions.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: d4uu,d4vv,d4ww,d4qq,d4TT	! <-- 4th order structure functions.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: ruu,rvv,rww,rqq,rTT	! <-- Autocorrelations.
			! Placeholder reals for the laged and advanced means and standrad deviations:
			REAL*8					:: muk,mukj,suk,sukj
			REAL*8					:: mvk,mvkj,svk,svkj
			REAL*8					:: mwk,mwkj,swk,swkj
			REAL*8					:: mTk,mTkj,sTk,sTkj
			REAL*8					:: mqk,mqkj,shk,shkj


			! Set the dimensions as specified by the user and the input file.
			mm=ecv%nn
			ss=ecv%drwinl
			lag1=ecv%lag1
			

			PRINT*,'*******************************************************************************'
			PRINT*,'Starting structure function & autocorrelation computation'
			WRITE(ecv%loglun,*)'*******************************************************************************'
			WRITE(ecv%loglun,*)'Starting structure function & autocorrelation computation'

			sspmm=FLOOR(1.0*mm/ss)	! Number of segments in the series.
	
			! Allocate the dimensions of the segmented arrays.		
			ALLOCATE(us(ss,sspmm))
			ALLOCATE(vs(ss,sspmm))
			ALLOCATE(ws(ss,sspmm))
			ALLOCATE(Ts(ss,sspmm))
			ALLOCATE(qs(ss,sspmm))
			ALLOCATE(twindow(2,sspmm,SIZE(ecv%t,2)))
			ALLOCATE(sobukhovl(sspmm))
			sobukhovl(:)=0
			ALLOCATE(ustar(sspmm))
			ALLOCATE(szmes(sspmm))
			ALLOCATE(suh(sspmm))
			ALLOCATE(tempar(ss,5))	! <-- Allocate the gap filling array.



			DO j=1,sspmm
				ll=(j-1)*ss+1 		! Index lower limit.
				ul=j*ss			! Index upper limit.
				us(:,j)=ecv%u(ll:ul)
				vs(:,j)=ecv%v(ll:ul)
				ws(:,j)=ecv%w(ll:ul)
				Ts(:,j)=ecv%Ts(ll:ul)		! <-- Choice of Ts,Tv or TT.
				qs(:,j)=ecv%h(ll:ul)		! <-- Choice of specific or absolute humidity.

				twindow(1,j,:)=ecv%t(ll,:)	! <-- Allocate the temporal bounds for the windows.
				twindow(2,j,:)=ecv%t(ul,:)	

			
				! Verify that IF >10% of the segment lenght is NaN; if so exclude from analysis.
				ok=0				
				isnan=0
				k=1
				DO WHILE(ok.EQ.0)
					IF(us(k,j).EQ.ecv%nanreal) THEN
						isnan=isnan+1
						IF( ((1.0*isnan)/(1.0*ss)) .GE. 0.1 ) THEN	! IF > 10% nans flag window.
							us(:,j)=ecv%nanreal
							vs(:,j)=ecv%nanreal
							ws(:,j)=ecv%nanreal
							Ts(:,j)=ecv%nanreal
							qs(:,j)=ecv%nanreal
							sobukhovl(j)=ecv%nanreal	! <-- Effective inclusion flag (if zero not nanreal).
							szmes(j)=ecv%nanreal
							ok=1
						END IF
					END IF
					IF(k.GE.ss) THEN
						ok=1
					END IF
					k=k+1
				END DO
				
				
				! If included then:
				IF(sobukhovl(j).NE.ecv%nanreal) THEN
					
					IF(isnan.NE.0) THEN	! If some nans gap fill:
						! Interpolate flaged entries (when <10% of the segment) using the gapfill function. 
						tempar(:,1)=us(:,j)
						tempar(:,2)=vs(:,j)
						tempar(:,3)=ws(:,j)
						tempar(:,4)=Ts(:,j)
						tempar(:,5)=qs(:,j)
						CALL gapfill(ecv,ss,5,tempar)
						us(:,j)=tempar(:,1)
						vs(:,j)=tempar(:,2)
						ws(:,j)=tempar(:,3)
						Ts(:,j)=tempar(:,4)
						qs(:,j)=tempar(:,5)
					END IF
			

					! Rotate u into the (streamline plane) mean velocity for the given block,
					! and v into the cross-wind direction.
					CALL rotatestream(ecv,ss,us(:,j),vs(:,j))
	

					IF(ecv%detrend.EQV..TRUE.) THEN	! <-- Check if detrend is specified by the user.
					! If so detrend the subsegments before calculating the structure functions
					! and the autocorrelations.
						CALL detrend(ecv,ss,us(:,j))
						CALL detrend(ecv,ss,vs(:,j))
						CALL detrend(ecv,ss,ws(:,j))
						CALL detrend(ecv,ss,Ts(:,j))
						CALL detrend(ecv,ss,qs(:,j))
					END IF
					
					! Compute Obukhov length 'manually' to avoid dividing by detrended
					! block averaged sonic temperature which is near zero.
					tempreal=0
					tetw=0
					CALL temporalmean(ecv,ss,ecv%Ts(ll:ul),tempreal)
					CALL covariance(ecv,ss,Ts(:,j),ws(:,j),tetw)
					CALL frictionvelocity(ecv,ss,us(:,j),vs(:,j),ws(:,j),ustar(j))
					sobukhovl(j)=-1.0*(ustar(j)**3)/(0.4*9.81*tetw/tempreal)	
					
					! Measurement height for the autostatistical window.
					CALL temporalmean(ecv,ss,ecv%zmes(ll:ul),szmes(j))

					! Mean 'horizontal' velocity for the autostatistical window.
					tempreal=0
					CALL temporalmean(ecv,ss,ecv%u(ll:ul),tempreal)
					CALL temporalmean(ecv,ss,ecv%v(ll:ul),suh(j))
					suh(j)=SQRT(tempreal**2+suh(j)**2)

					tevw=0
					teuw=0
					tetw=0
					tempreal=0					
				END IF
			END DO
			
			ssj=NINT(5.0*60.0/lag)		! Number of lags to consider.
			! The i'th time lag is (i-1)*di
			di=NINT(lag1*ecv%fs)		! Number of points in the time lag for i=1.	
			stime=(ss-1)*ecv%dt		! Total time in the segment to be analyzed.
			ul=ssj				! Number of time lags to consider; must ensure lag <= segment length/2.
			ALLOCATE(tlag(ul))		! Time lag vector.
			DO j=1,ul
				tlag(j)=j*lag
			END DO

			! Allocate space for the structure function and autocorrelation arrays.
			ALLOCATE(d2uu(ul,sspmm))
			ALLOCATE(d3uu(ul,sspmm))
			ALLOCATE(d4uu(ul,sspmm))
			ALLOCATE(ruu(ul,sspmm))

			ALLOCATE(d2vv(ul,sspmm))
			ALLOCATE(d3vv(ul,sspmm))
			ALLOCATE(d4vv(ul,sspmm))
			ALLOCATE(rvv(ul,sspmm))
			
			ALLOCATE(d2ww(ul,sspmm))
			ALLOCATE(d3ww(ul,sspmm))
			ALLOCATE(d4ww(ul,sspmm))
			ALLOCATE(rww(ul,sspmm))
		
			ALLOCATE(d2TT(ul,sspmm))
			ALLOCATE(d3TT(ul,sspmm))
			ALLOCATE(d4TT(ul,sspmm))
			ALLOCATE(rTT(ul,sspmm))
			
			ALLOCATE(d2qq(ul,sspmm))
			ALLOCATE(d3qq(ul,sspmm))
			ALLOCATE(d4qq(ul,sspmm))
			ALLOCATE(rqq(ul,sspmm))


			DO j=1,sspmm
				IF(sobukhovl(j).NE.ecv%nanreal) THEN ! Only compute for complete (non flagged) segments.
					DO i=1,ul
						l=(i)*di	! i*di; alternately (i-1)*di to include zero lag.

						! Means of the laged subsegments.						
						CALL temporalmean(ecv,ss-l,us(1:ss-l,j),muk)
						CALL temporalmean(ecv,ss-l,vs(1:ss-l,j),mvk)
						CALL temporalmean(ecv,ss-l,ws(1:ss-l,j),mwk)
						CALL temporalmean(ecv,ss-l,Ts(1:ss-l,j),mTk)
						CALL temporalmean(ecv,ss-l,qs(1:ss-l,j),mqk)

						! Means of the advanced subsegments.
						CALL temporalmean(ecv,ss-l,us(l+1:ss,j),mukj)
						CALL temporalmean(ecv,ss-l,vs(l+1:ss,j),mvkj)
						CALL temporalmean(ecv,ss-l,ws(l+1:ss,j),mwkj)
						CALL temporalmean(ecv,ss-l,Ts(l+1:ss,j),mTkj)
						CALL temporalmean(ecv,ss-l,qs(l+1:ss,j),mqkj)

						DO k=1,ss-l
							! Components of the autuocorrelations.
							suk=suk+(us(k,j)-muk)**2
							sukj=sukj+(us(k+l,j)-mukj)**2
							ruu(i,j)=ruu(i,j)+(us(k,j)-muk)*(us(k+l,j)-mukj)
		
							svk=svk+(vs(k,j)-mvk)**2
							svkj=svkj+(vs(k+l,j)-mvkj)**2
							rvv(i,j)=rvv(i,j)+(vs(k,j)-mvk)*(vs(k+l,j)-mvkj)
							
							
							swk=swk+(ws(k,j)-mwk)**2
							swkj=swkj+(ws(k+l,j)-mwkj)**2
							rww(i,j)=rww(i,j)+(ws(k,j)-mwk)*(ws(k+l,j)-mwkj)

							sTk=sTk+(Ts(k,j)-mTk)**2
							sTkj=sTkj+(Ts(k+l,j)-mTkj)**2
							rTT(i,j)=rTT(i,j)+(Ts(k,j)-mTk)*(Ts(k+l,j)-mTkj)	

							shk=shk+(qs(k,j)-mqk)**2
							shkj=shkj+(qs(k+l,j)-mqkj)**2
							rqq(i,j)=rqq(i,j)+(qs(k,j)-mqk)*(qs(k+l,j)-mqkj)	
													


							! Components of the 2nd order structure functions.
							d2uu(i,j)=d2uu(i,j)+(us(k,j)-us(k+l,j))**2
							d2vv(i,j)=d2vv(i,j)+(vs(k,j)-vs(k+l,j))**2
							d2ww(i,j)=d2ww(i,j)+(ws(k,j)-ws(k+l,j))**2
							d2TT(i,j)=d2TT(i,j)+(Ts(k,j)-Ts(k+l,j))**2
							d2qq(i,j)=d2qq(i,j)+(qs(k,j)-qs(k+l,j))**2

							! Components of the 3rd order structure functions.
							d3uu(i,j)=d2uu(i,j)+(us(k,j)-us(k+l,j))**3
							d3vv(i,j)=d2vv(i,j)+(vs(k,j)-vs(k+l,j))**3
							d3ww(i,j)=d2ww(i,j)+(ws(k,j)-ws(k+l,j))**3
							d3TT(i,j)=d2TT(i,j)+(Ts(k,j)-Ts(k+l,j))**3
							d3qq(i,j)=d2qq(i,j)+(qs(k,j)-qs(k+l,j))**3
							

							! Components of the 4th order structure functions.
							d4uu(i,j)=d4uu(i,j)+(us(k,j)-us(k+l,j))**4
							d4vv(i,j)=d4vv(i,j)+(vs(k,j)-vs(k+l,j))**4
							d4ww(i,j)=d4ww(i,j)+(ws(k,j)-ws(k+l,j))**4
							d4TT(i,j)=d4TT(i,j)+(Ts(k,j)-Ts(k+l,j))**4
							d4qq(i,j)=d4qq(i,j)+(qs(k,j)-qs(k+l,j))**4


						END DO

						! Variance to standard deviation.
						suk=SQRT(suk)
						sukj=SQRT(sukj)
						svk=SQRT(svk)
						svkj=SQRT(svkj)
						swk=SQRT(swk)
						swkj=SQRT(swkj)
						sTk=SQRT(sTk)
						sTkj=SQRT(sTkj)
						shk=SQRT(shk)
						shkj=SQRT(shkj)

						! Value of the 2nd order structure function for this segment at lag (i-1)*lag1
						d2uu(i,j)=d2uu(i,j)/ss
						d2vv(i,j)=d2vv(i,j)/ss
						d2ww(i,j)=d2ww(i,j)/ss
						d2TT(i,j)=d2TT(i,j)/ss
						d2qq(i,j)=d2qq(i,j)/ss

						! Value of the 3rd order structure function at this time lag.
						d3uu(i,j)=d3uu(i,j)/ss
						d3vv(i,j)=d3vv(i,j)/ss
						d3ww(i,j)=d3ww(i,j)/ss
						d3TT(i,j)=d3TT(i,j)/ss
						d3qq(i,j)=d3qq(i,j)/ss


						! Value of the 4th order structure function at this time lag.
						d4uu(i,j)=d4uu(i,j)/ss
						d4vv(i,j)=d4vv(i,j)/ss
						d4ww(i,j)=d4ww(i,j)/ss
						d4TT(i,j)=d4TT(i,j)/ss
						d4qq(i,j)=d4qq(i,j)/ss

						! Value of autocorrelation for this segment at lag (i-1)*lag1
						ruu(i,j)=ruu(i,j)/(suk*sukj)
						rvv(i,j)=rvv(i,j)/(svk*svkj)
						rww(i,j)=rww(i,j)/(swk*swkj)
						rTT(i,j)=rTT(i,j)/(sTk*sTkj)
						rqq(i,j)=rqq(i,j)/(shk*shkj)

						! Set means and standard deviations to zero for the next pass.
						suk=0
						svk=0
						swk=0
						sTk=0
						shk=0
						sukj=0
						svkj=0
						swkj=0
						sTkj=0
						shkj=0
						muk=0
						mukj=0
						mvk=0
						mvkj=0
						mwk=0
						mwkj=0
						mTk=0
						mTkj=0
						mqk=0
						mqkj=0

					END DO
				ELSE	! If nan flag -> set to nan.
					d2uu(:,j)=ecv%nanreal
					d2vv(:,j)=ecv%nanreal
					d2ww(:,j)=ecv%nanreal
					d2TT(:,j)=ecv%nanreal
					d2qq(:,j)=ecv%nanreal

					d4uu(:,j)=ecv%nanreal
					d4vv(:,j)=ecv%nanreal
					d4ww(:,j)=ecv%nanreal
					d4TT(:,j)=ecv%nanreal
					d4qq(:,j)=ecv%nanreal

					ruu(:,j)=ecv%nanreal
					rvv(:,j)=ecv%nanreal
					rww(:,j)=ecv%nanreal
					rTT(:,j)=ecv%nanreal
					rqq(:,j)=ecv%nanreal
				END IF
			END DO

			! Scan through the arrays to determine the number of active windows.
			jj=0
			DO j=1,SIZE(d2uu,2)
				IF(d2uu(2,j).NE.ecv%nanreal) THEN
					jj=jj+1		! <-- Update number of active windows.
				END IF
			END DO
			ii=SIZE(d2uu,1)
			PRINT '("-->",I5," active statistical windows")',jj
			WRITE(ecv%loglun,*)"-->",jj," active statistical windows"
			IF(jj.EQ.0) THEN	! <-- Break out of the subroutine if no windows are active.
				ecv%dractive=.FALSE.
				RETURN
			ELSE
				ecv%dractive=.TRUE.
			END IF


			! Smoothing matrix:
			nc=90
			! Allocate the central time lag array
			ALLOCATE(tlagc(nc))
			! Allocate the Konno_Ohmachi smoothing matrix.
			ALLOCATE(KOSM(nc,ii))
			! Compute the Konno-Ohmachi smoothing matrix, KOSM, with nfc smoothing windows.
			CALL	kosmoothing(ecv,tlag,ii,tlagc,nc,10.0D0,KOSM)
			ALLOCATE(ecv%tlag(nc))
			ecv%tlag(:)=tlagc(:)	! <-- Pass the time lag array to the module.


			
			! Allocate the dimensions of the output structure functions and autocorrelations.
			ALLOCATE(ecv%ruu(nc,jj))
			ALLOCATE(ecv%rvv(nc,jj))
			ALLOCATE(ecv%rww(nc,jj))
			ALLOCATE(ecv%rTT(nc,jj))
			ALLOCATE(ecv%rqq(nc,jj))

			ALLOCATE(ecv%d2uu(nc,jj))
			ALLOCATE(ecv%d2vv(nc,jj))
			ALLOCATE(ecv%d2ww(nc,jj))
			ALLOCATE(ecv%d2TT(nc,jj))
			ALLOCATE(ecv%d2qq(nc,jj))

			ALLOCATE(ecv%d3uu(nc,jj))
			ALLOCATE(ecv%d3vv(nc,jj))
			ALLOCATE(ecv%d3ww(nc,jj))
			ALLOCATE(ecv%d3TT(nc,jj))
			ALLOCATE(ecv%d3qq(nc,jj))

			ALLOCATE(ecv%d4uu(nc,jj))
			ALLOCATE(ecv%d4vv(nc,jj))
			ALLOCATE(ecv%d4ww(nc,jj))
			ALLOCATE(ecv%d4TT(nc,jj))
			ALLOCATE(ecv%d4qq(nc,jj))

			! Allocate the statistical obukhov length and measurement height.
			ALLOCATE(ecv%drobukhovl(jj))
			ALLOCATE(ecv%drz(jj))
			ALLOCATE(ecv%druh(jj))

			! Allocate the statistical window timestamp.
			ALLOCATE(ecv%tstwin(2,jj,6))

			! Allocate the entries to the arrays, only including active (non nan) windows.
			j=0
			k=0
			i=0
			l=0
			DO WHILE(j.LT.jj)
				k=k+1	! <-- Update window counter.
				IF(d2uu(2,k).NE.ecv%nanreal) THEN
					j=j+1	! <-- Update active window counter.
					ecv%d2uu(:,j)=MATMUL(KOSM,d2uu(:,k))
					ecv%d2vv(:,j)=MATMUL(KOSM,d2vv(:,k))
					ecv%d2ww(:,j)=MATMUL(KOSM,d2ww(:,k))
					ecv%d2TT(:,j)=MATMUL(KOSM,d2TT(:,k))
					ecv%d2qq(:,j)=MATMUL(KOSM,d2qq(:,k))

					ecv%d3uu(:,j)=MATMUL(KOSM,d3uu(:,k))
					ecv%d3vv(:,j)=MATMUL(KOSM,d3vv(:,k))
					ecv%d3ww(:,j)=MATMUL(KOSM,d3ww(:,k))
					ecv%d3TT(:,j)=MATMUL(KOSM,d3TT(:,k))
					ecv%d3qq(:,j)=MATMUL(KOSM,d3qq(:,k))

					ecv%d4uu(:,j)=MATMUL(KOSM,d4uu(:,k))
					ecv%d4vv(:,j)=MATMUL(KOSM,d4vv(:,k))
					ecv%d4ww(:,j)=MATMUL(KOSM,d4ww(:,k))
					ecv%d4TT(:,j)=MATMUL(KOSM,d4TT(:,k))
					ecv%d4qq(:,j)=MATMUL(KOSM,d4qq(:,k))
		
					ecv%ruu(:,j)=MATMUL(KOSM,ruu(:,k))
					ecv%rvv(:,j)=MATMUL(KOSM,rvv(:,k))
					ecv%rww(:,j)=MATMUL(KOSM,rww(:,k))
					ecv%rTT(:,j)=MATMUL(KOSM,rTT(:,k))
					ecv%rqq(:,j)=MATMUL(KOSM,rqq(:,k))

					ecv%drobukhovl(j)=sobukhovl(k)
					ecv%drz(j)=szmes(k)
					ecv%druh(j)=suh(k)

					ecv%tstwin(:,j,:)=twindow(:,k,:)

				END IF
			END DO
					
			PRINT*,'End of structure function & autocorrelation computation'	
			PRINT*,'*******************************************************************************'
			WRITE(ecv%loglun,*)'End of structure function & autocorrelation computation'	
			WRITE(ecv%loglun,*)'*******************************************************************************'


			! Deallocate temporary arrays.
			DEALLOCATE(twindow)
			DEALLOCATE(tlag)
			DEALLOCATE(tlagc)
			DEALLOCATE(KOSM)
			DEALLOCATE(us)
			DEALLOCATE(vs)
			DEALLOCATE(ws)
			DEALLOCATE(Ts)
			DEALLOCATE(qs)
			DEALLOCATE(ruu)
			DEALLOCATE(rvv)
			DEALLOCATE(rww)
			DEALLOCATE(rqq)
			DEALLOCATE(rTT)
			DEALLOCATE(d2uu)
			DEALLOCATE(d2vv)
			DEALLOCATE(d2ww)
			DEALLOCATE(d2qq)
			DEALLOCATE(d2TT)
			DEALLOCATE(d3uu)
			DEALLOCATE(d3vv)
			DEALLOCATE(d3ww)
			DEALLOCATE(d3qq)
			DEALLOCATE(d3TT)
			DEALLOCATE(d4uu)
			DEALLOCATE(d4vv)
			DEALLOCATE(d4ww)
			DEALLOCATE(d4TT)
			DEALLOCATE(d4qq)
			DEALLOCATE(sobukhovl)
			DEALLOCATE(szmes)
			DEALLOCATE(tempar)

			RETURN
		END SUBROUTINE



		! Linked subroutines computing the median (intrinsic does not exist in FORTRAN), first the
		! quicksort algorithm is used to sort an input array from lowest to highest values recursively in
		! the quicksort subroutine in which the array is partitioned using the partition subroutine. Finally
		! once the array 'a' is sorted computing the median is straight forwards; for an even number of entries
		! it is the midpoint interpolation between the two midpoint values in the sorted array. Likewise if
		! the array is of odd length the median is the value of the (now defined) central point of the sorted array.

		! Note that the quicksort algorithm and the related parition algorithm were originally written in C
		! by Rew, J. and made fortran compatible by Brainerd, W (url: http://www.fortran.com/qsort_c.f95).
		! The author of this module simply adapted the algorithm, to the needs of this module, in particular the 
		! need for computing medians in the MAD despiking algorithm of Mauder et al. (2013). 
 
		
		! Can be called upon independently of the MADspike subroutine.
		SUBROUTINE	median(ecv,mm,a,medval)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: 		ecv
			INTEGER,INTENT(IN)		:: mm		! Segment length and weight.
			INTEGER				:: mmm		! Effective segment length (discounting nans).
			REAL*8,DIMENSION(mm),INTENT(IN) :: a		! Array segment.
			REAL*8,DIMENSION(mm)		:: asorted	! Sorted array.
			REAL*8,ALLOCATABLE,DIMENSION(:) :: temp(:)	! Temporary array (for discounting nans in the median).
			REAL*8,INTENT(OUT)		:: medval	! Output median value.
			INTEGER				:: i,j 		! Counter.
			INTEGER				:: eofnans	! Status integer.
			 
			! Initialize values to be zero.
			medval=0
			eofnans=0
			asorted(:)=0
			
			asorted=a ! Assign the unsorted values to the as array, as we dont want to sort the actual array in the structure.
			
			! Sort the array a in ascending order. Remember ecv%nanreal=-1e3 is the absolute minimum to expect.
			CALL quicksort(asorted)
			
			 
			! Check for nan entries in the now sorted array a.
			j=0
			DO WHILE(eofnans.EQ.0)
				IF(asorted(j+1).EQ.ecv%nanreal) THEN
					j=j+1
				ELSE
					eofnans=1	! No more nan entries.
				END IF
			END DO
			
			! Creat a (possibly) truncated version of array a, j is the entry of the last nan,
			! so the temporary truncated sorted array should have dimensions mm-j. (j can be zero).
			mmm=mm-j
			ALLOCATE(temp(mmm))
			temp(:)=0
			temp(1:mmm)=asorted(j+1:mm)	! Assign the sorted (and no nan) values to the temporary array.
			
			!PRINT*,temp
			
			! Finally compute the median value of the temporary array.
			IF (MOD(mmm,2).EQ.0) THEN ! Even array length
				medval=(temp(mmm/2)+temp(1+mmm/2))/2.0
			ELSE ! Odd array length.
				medval=temp((1+mmm)/2)
			END IF
			
			! Deallocate the temporary array.
			DEALLOCATE(temp)

			RETURN
		END SUBROUTINE median
		
		! The quicksort algorithm of Rew, J. (2003), used internally in the median subroutine.
		RECURSIVE SUBROUTINE quicksort(a)
			IMPLICIT NONE
			REAL*8, INTENT(INOUT), DIMENSION(:) :: a
			INTEGER :: iq

			IF(SIZE(A).GT.1) THEN
				CALL partition(a, iq)
				CALL quicksort(a(:iq-1))
				CALL quicksort(a(iq:))
			END IF
		END SUBROUTINE quicksort

		! The partition algorithm of Rew, J. (2003), used internally in the quicksort and hence
		! median subroutines.
		SUBROUTINE partition(a, marker)
			IMPLICIT NONE
			REAL*8, INTENT(INOUT), DIMENSION(:) :: a ! Array to be partitioned.
			INTEGER, INTENT(OUT) :: marker	! Marks the element in the array.
			INTEGER :: i, j	 ! Counters.
			REAL*8 :: temp	 ! Temporary real
			REAL*8 :: x      ! Pivot point
			x = a(1)	! First entry in the array as x.
			i= 0		! Left bounds of array.
			j= SIZE(a) + 1 	! Right bounds of array.

			DO 
				j = j-1	! Move left from right.
				DO
					IF (a(j).LE.x) EXIT ! Move further left from right.
					j = j-1 
				END DO
				i = i+1	! Move right from left.
				DO
					IF (a(i).GE.x) EXIT
					i = i+1
				END DO
				IF (i.LT.j) THEN
					! exchange a(i) and a(j)
					temp = a(i)
					a(i) = a(j)
					a(j) = temp
				ELSEIF (i.EQ.j) THEN
					marker = i+1
					RETURN
				ELSE
					marker = i
					RETURN
				END IF
			END DO

		END SUBROUTINE partition

		!----------------------------------------------------------------------------------------------------------------------
		! Spectral analysis.
		!----------------------------------------------------------------------------------------------------------------------

		! Calculates the Fourier spectrum (forward Fourier transform) of the series a in segments. We make use of the
		! fast fourier transform algorithm implemented in the dfftpack package (by H.C. Pumphrey), which is a double
		! precision clone of the fftpack by P.N.Swarztrauber (1985).
		! Computes the forward fourier transform of all of the measured/corrected variables at once to save
		! computational time involved in the initialization of the work array w required for the fourier transforms.
		SUBROUTINE spectrum(ecv,mm,ss,isodd,fu,fv,fw,fT,fq)
			IMPLICIT NONE
			TYPE(eddycov)			:: ecv

			EXTERNAL dffti
			EXTERNAL dfftf
			INTEGER,INTENT(IN)			:: mm,ss			! Series length, segment length.
			INTEGER					:: sspmm			! Nearest integer number of segments in the series.
			REAL*8,DIMENSION(mm)			:: u,v,w,T,q			! Arrays to be Fourier transformed.
			REAL*8,ALLOCATABLE,DIMENSION(:)		:: us,vs,ws,Ts,qs 		! A single segment of the above arrays.
			! Output forward Fourier transforms.			
			COMPLEX*16,ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)	:: fu,fv,fw,fT,fq
			REAL*8,DIMENSION(2*ss+15)		:: work				! Work array required for dfftf.
			INTEGER					:: i,j,k,l,ul,ll,isnan,ok	! Counters/Indices.
			INTEGER, INTENT(OUT)			:: isodd			! Status integer, for window length.
			INTEGER					:: iib				! Number of points in a block averaging period.
			REAL*8					:: pi
			REAL*8					:: taper			! Temporary real to hold the value of the taper
												! function at a given point in a series.
			REAL*8					:: compf			! Compensation factor for the taper.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)	:: tempb			! Temporary array for gap filling.
			REAL*8					:: tempmean

			! Define pi.
			pi=4.0*ATAN(1.0)
			

			! From dfft source code:
			! A call for dfftf(n,a,w) returns the forward fourier transform of 'a' as a new array 'a'
			! where:
			! a(1) is the sum from i=1 to n of a.

			! if n is even l=n/2, if odd l=(n+1)/2 
			! for k=2,...,l

			! a(2*k-2)=the sum from i=1 to n of r(i)*cos((k-1)*(i-1)*2pi/n)		---> real part of ft.
			! a(2*k-1)=the sum from i=1 to n of -r(i)*sin((k-1)*(i-1)*2pi/n)	---> complex part of ft.

			! if n is even then:
			!              r(n) = the sum from i = 1 to i = n of  (-1)**(i-1)*r(i)

			! NB: Remember to normalize by n, as the ft is unormalized here.
			
			! Compute the Nyquist index for the segment.
			IF(MOD(ss,2).NE.0) THEN	! The segment length is odd.
				l=(ss+1)/2
				isodd=1
			ELSE	! Segment length is even.
				l=ss/2
				isodd=0
			END IF
			
			iib=ecv%blockl
			
			!sspmm=FLOOR(1.0*mm/ss) ! Number of averaging periods in the entire series.
			sspmm=FLOOR((1.0*mm)/iib)-NINT((1.0*ss)/iib)+1 
			
			! Allocate the arrays to be transformed.
			u=ecv%u(1:mm)
			v=ecv%v(1:mm)
			w=ecv%w(1:mm)
			T=ecv%Ts(1:mm)	! <--- Choice of Ts,Tv or TT.
			q=ecv%h(1:mm)	! <--- Choice of absolute or specific humidity.

			! Allocation of dimensions.
			! A single segment of the arrays to be transformed.			
			ALLOCATE(us(ss))
			ALLOCATE(vs(ss))
			ALLOCATE(ws(ss))
			ALLOCATE(Ts(ss))
			ALLOCATE(qs(ss))

			ALLOCATE(tempb(ss,5))	! <-- Allocate the dimensions of the array used for gap filling.

			! Allocate the output spectrum.
			IF(isodd.EQ.1) THEN ! If odd series length the nyquist frequency is not resolved.
				ALLOCATE(fu(l,sspmm))
				ALLOCATE(fv(l,sspmm))
				ALLOCATE(fw(l,sspmm))
				ALLOCATE(fT(l,sspmm))
				ALLOCATE(fq(l,sspmm))
			ELSE	! If even series length the nyquist frequency is resolved, and we add an extra entry.
				ALLOCATE(fu(l+1,sspmm))
				ALLOCATE(fv(l+1,sspmm))
				ALLOCATE(fw(l+1,sspmm))
				ALLOCATE(fT(l+1,sspmm))
				ALLOCATE(fq(l+1,sspmm))
			END IF

			CALL dffti(ss,work) ! Initialize the fft (only need to do this once).
			
			! Initialize the spectrum as zeros, for multiple calls.
			fu(:,:)=DCMPLX(0,0)
			fv(:,:)=DCMPLX(0,0)
			fw(:,:)=DCMPLX(0,0)
			fT(:,:)=DCMPLX(0,0)
			fq(:,:)=DCMPLX(0,0)
			
			isnan=0 ! Intialize that nan flag as zero (not nan).
			DO j=1,sspmm
				ll=(j-1)*iib+1		! Index lower limit.
				ul=ll+ss-1		! Index upper limit.
				
				! Allocate the array segments.
				us=u(ll:ul)
				vs=v(ll:ul)
				ws=w(ll:ul)
				Ts=T(ll:ul)
				qs=q(ll:ul)

				ok=0
				isnan=0
				k=1
				DO WHILE(ok.EQ.0)
					IF(us(k).EQ.ecv%nanreal) THEN
						isnan=isnan+1
						IF( ((1.0*isnan)/(1.0*ss)) .GT. 0.1 ) THEN	! IF > 10% nans flag window.
							fu(:,j)=DCMPLX(ecv%nanreal,0)
							fv(:,j)=DCMPLX(ecv%nanreal,0)
							fw(:,j)=DCMPLX(ecv%nanreal,0)
							fT(:,j)=DCMPLX(ecv%nanreal,0)
							fq(:,j)=DCMPLX(ecv%nanreal,0)
							ok=1
						END IF
					END IF
					IF(k.GE.ss) THEN
						ok=1
					END IF
					k=k+1
				END DO

				
				! Only compute spectra for segments with nnans < 10%.
				IF(DBLE(fu(1,j)).NE.ecv%nanreal) THEN
					IF(isnan.NE.0) THEN
						! Interpolate flaged entries (when <1% but not 0% of the segment) using the gapfill function. 
						tempb(:,1)=us(:)
						tempb(:,2)=vs(:)
						tempb(:,3)=ws(:)
						tempb(:,4)=Ts(:)
						tempb(:,5)=qs(:)
						CALL gapfill(ecv,ss,5,tempb)
						us(:)=tempb(:,1)
						vs(:)=tempb(:,2)
						ws(:)=tempb(:,3)
						Ts(:)=tempb(:,4)
						qs(:)=tempb(:,5)
					END IF	
					tempb(:,:)=0






					! Rotate u into the (streamline plane) mean velocity for the given block,
					! and v into the cross-wind direction. Of course this does NOT change
					! meanuh(j).
					CALL rotatestream(ecv,ss,us,vs)
	
					IF(ecv%detrend.EQV..TRUE.) THEN	! <-- Check if detrend is specified by the user.
						CALL detrend(ecv,ss,us)
						CALL detrend(ecv,ss,vs)
						CALL detrend(ecv,ss,ws)
						CALL detrend(ecv,ss,Ts)
						CALL detrend(ecv,ss,qs)
					END IF
					! Add a taper to the series if specified by the user.
					IF(ecv%taper.EQV..TRUE.) THEN
						IF(ecv%tapertype.EQ.'B') THEN		! Use the Bell taper.
							DO i=1,NINT(0.1*ss)		! Apply Bell taper @ start of window.
								taper=SIN(5*pi*(i-1)/((1.0)*(ss-1)))**2
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper
							END DO
							DO i=NINT(0.9*ss),ss		! Apply Bell taper @ end of window.
								taper=SIN(5*pi*(i-1)/((1.0)*(ss-1)))**2
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper		
							END DO
						ELSEIF(ecv%tapertype.EQ.'H') THEN	! Use the Hamming window.
							DO i=1,ss
								taper=0.54-0.46*COS(2*pi*(i-1)/(ss-1))
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper
							END DO
						ELSE					! Make sure a taper is specified if taperedspectra=TRUE.
							PRINT*,'Please specify taper B or H'
							STOP
						END IF
						! Furthermore make sure to remove the residual mean introduced by the tapering operation.
						CALL temporalmean(ecv,ss,us,tempmean)
						us=us-tempmean
						CALL temporalmean(ecv,ss,vs,tempmean)
						vs=vs-tempmean
						CALL temporalmean(ecv,ss,ws,tempmean)
						ws=ws-tempmean
						CALL temporalmean(ecv,ss,Ts,tempmean)
						Ts=Ts-tempmean
						CALL temporalmean(ecv,ss,qs,tempmean)
						qs=qs-tempmean

						! Also scale the series to account for the variance lost by the tapering using
						! the compensation factors (after Kaimal et al. 1989).
						IF(ecv%tapertype.EQ.'B') THEN
							compf=ecv%cfB
						ELSEIF(ecv%tapertype.EQ.'H') THEN
							compf=ecv%cfH
						END IF
						us=compf*us
						vs=compf*vs
						ws=compf*ws
						Ts=compf*Ts
						qs=compf*qs
					END IF			
				


	
					CALL dfftf(ss,us,work)	! Compute the forward Fourier transform.
					CALL dfftf(ss,vs,work)
					CALL dfftf(ss,ws,work)
					CALL dfftf(ss,Ts,work)
					CALL dfftf(ss,qs,work)
					
					! Normalize the transforms by the segment length.
					us=us/ss
					vs=vs/ss
					ws=ws/ss
					Ts=Ts/ss
					qs=qs/ss

	
					! Allocate values to the output spectrum.
					fu(1,j)=DCMPLX(us(1),0)	! The mean component; zeroth mode.
					fv(1,j)=DCMPLX(vs(1),0)
					fw(1,j)=DCMPLX(ws(1),0)
					fT(1,j)=DCMPLX(Ts(1),0)
					fq(1,j)=DCMPLX(qs(1),0)


					! Remaining components:
					DO k=2,l 
						fu(k,j)=DCMPLX(us(2*k-2),us(2*k-1))
						fv(k,j)=DCMPLX(vs(2*k-2),vs(2*k-1))
						fw(k,j)=DCMPLX(ws(2*k-2),ws(2*k-1))
						fT(k,j)=DCMPLX(Ts(2*k-2),Ts(2*k-1))
						fq(k,j)=DCMPLX(qs(2*k-2),qs(2*k-1))
					END DO
					IF(isodd.EQ.0) THEN ! If even length segment allocate the missing last entry->f_ny resolved.
						fu(l+1,j)=DCMPLX(us(ss),0)
						fv(l+1,j)=DCMPLX(vs(ss),0)
						fw(l+1,j)=DCMPLX(ws(ss),0)
						fT(l+1,j)=DCMPLX(Ts(ss),0)
						fq(l+1,j)=DCMPLX(qs(ss),0)
					END IF
				END IF
				isnan=0 ! Reset nan status.
				ll=0
				ul=0
				us(:)=0
				vs(:)=0
				ws(:)=0
				Ts(:)=0
				qs(:)=0
			
			END DO


			RETURN
		END SUBROUTINE		






		! Calculates a variety of cross spectra for the two series a and b segment wise.
		! See Stull (1988) for a discussion of cross spectra; these are the cospectra, quadratures,
		! and finally the ogive spectrum.
		! The power spectrum (spectral density) of the respective series is also calculated and compared to the
		! biased variance as a control; this is also of interest in identifying the presence of an inertial subrange. 
		! These spectra are all passed to the ecv structure from the function call.
		! The time series that are analyzed can be detrended and tapered IF desired by the user. Moreover
		! a logarithmic smoothing is applied to the final spectra for presentation purposes and to reduce
		! memory requirements.
		SUBROUTINE 	crossspectrum(ecv)
			IMPLICIT NONE
			TYPE(eddycov)			:: ecv
			INTEGER						:: mm,ss			! Series length, segment length.
			INTEGER						:: sspmm			! Nearest integer number of segments in the series.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: sobukhovl			! Obukhov length, L, for the spectral windows.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: Tstar			! Dynamic temperature for spec windows. 
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: szmes			! Measurement height,z, for the spectral windows.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: u,v,w,T,q			! Arrays to be transformed.
			REAL*8,ALLOCATABLE,DIMENSION(:,:,:)		:: windowtime			! Keep track of times of active (clean) windows.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: us,vs,ws,Ts,qs,ustar	! A single segment of the input arrays.
			COMPLEX*16,ALLOCATABLE,DIMENSION(:,:)		:: fa,fb   			! Spectra; computed using the spectrum subroutine.
			COMPLEX*16,ALLOCATABLE,DIMENSION(:,:)		:: fu,fv,fw,fT,fq		! Spectra; computed using the spectrum subroutine.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: su,sv,sw,sT,sh		! Output power spectra (spectral density).
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: siu,siv,siw,siT,siq		! Intensity spectra (used for variance Ogives).
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: couw,covw,coTw,cohw		! Output cospectral density.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: ciuw,civw,ciTw,ciqw		! Cospectral intensity.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: quuw,quvw,quTw,quqw		! Output quadratures (quadrature density).
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: oguw,ogvw,ogTw,oghw		! Ogives.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: oguu,ogvv,ogww,ogTT,oghh	
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: freq				! Frequency array [Hz].
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: fogs				! Frequency for the ogives.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: varu,varv,varw,varT,varq	! Window statistical variances for control.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: svaru,svarv,svarw,svarT,svarq! Window spectral variances for control.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: covuw,covvw,covTw,covqw	! Window statistical covariances for control.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: meanuh			! Mean horizontal wind speed for each window.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: tempa,tempb			! Temporary arrays for gap filling if worthwhile.
			REAL*8						:: tsbar,tempvb,tempub,tempmean ! Temporary reals.
			INTEGER						:: i,j,k,l,ul,ll,isnan,a,b,ok	! Counters/Indices.
			INTEGER						:: ii,jj,iib			! Sizes of the output spectra.
			INTEGER						:: isodd			! Status: odd or even segment length.
			REAL*8						:: etol				! Error tolerance in the spectral analysis.
			REAL*8						:: dvar,dcov			! Differece between spectral&temporal covs&vars.
			INTEGER						:: awins			! Number of active (nanfree) spectral windows.
			! Counters for controlling the spectral variances and covariances.
			INTEGER						:: okvaru,okvarv,okvarw,okvarT,okvarq
			INTEGER						:: okcovuw,okcovvw,okcovTw,okcovqw
			REAL*8						:: pi
			REAL*8						:: taper			! Temporary real to hold the value of the taper
													! function at a given point in a series.
			REAL*8						:: compf			! Compensation factor for the taper.
			! Smoothing; see kosmoothing() for more details. The 'og' is for the Ogives which already have a different resolution than
			! the remaining spectra.
			REAL*8,ALLOCATABLE,DIMENSION(:,:)		:: KOSM,KOSMog			! Konno-Ohmachi smoothing matrix.
			REAL*8,ALLOCATABLE,DIMENSION(:)			:: fcen				! Central frequencies for smoothing; also
													! the output frequencies for spectra (not Ogives).
			REAL*8						:: deltaog			! Logarithmic frequency increments for Ogives.
													! spectra in the smoothing process.
			INTEGER						:: nfc,nfcog			! Output frequency array length.
			INTEGER						:: nf				! Number of frequencies to smooth.
			REAL*8						:: bwidth			! Bandwidth for smoothing.
			INTEGER,ALLOCATABLE,DIMENSION(:)		:: findi			! Index of high resolution frequencies corresponding
													! to Ogive frequencies.
													

			! Define pi.
			pi=4.0*ATAN(1.0)
			
			! Set dimensions as specified by the user and the input file.
			mm=ecv%nn
			ss=ecv%spwinl	! Spectral segment length [number of points in time: averaging period].
			iib=ecv%blockl	! Number of points in a block averaging period.
			ALLOCATE(u(mm))
			ALLOCATE(v(mm))
			ALLOCATE(w(mm))
			ALLOCATE(T(mm))
			ALLOCATE(q(mm))
						
			PRINT*,'*******************************************************************************'
			PRINT*,'Starting spectral analysis'
			WRITE(ecv%loglun,*)'*******************************************************************************'
			WRITE(ecv%loglun,*)'Starting spectral analysis'

			! Compute the spectra of u,v,w,T,q.	
			CALL spectrum(ecv,mm,ss,isodd,fu,fv,fw,fT,fq)
			
			! Get the length of the spectra, and verify that they are consistent.
			l=SIZE(fu,1)

			! Allocate the dimensions of the frequency array.
			ALLOCATE(freq(l))

			! Allocate the frequencies; entry i corresponds to frequency (i-1)/(ss*dt)
			DO i=1,l
				freq(i)=(1.0*i-1)/(ss*ecv%dt)
			END DO
			PRINT*,'Max frequency is:',MAXVAL(freq)
			PRINT*,'Min frequency is:',MINVAL(freq)	
			WRITE(ecv%loglun,*)'Max frequency is:',MAXVAL(freq)		
			WRITE(ecv%loglun,*)'Min frequency is:',MINVAL(freq)	

			!sspmm=FLOOR(1.0*mm/ss) ! Number of averaging periods in the entire series.
			! Number of spectral segment averaging periods in the entire series [advance one block averaging period at a time].
			sspmm=FLOOR((1.0*mm)/iib)-NINT((1.0*ss)/iib)+1 

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Power spectrum and variance control !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Allocate the arrays that are Fourier transformed for the variance/covariance control.
			u=ecv%u(1:mm)
			v=ecv%v(1:mm)
			w=ecv%w(1:mm)
			T=ecv%Ts(1:mm)		! Pick T or Tv or Ts.
			q=ecv%h(1:mm)		! Pick specific or absolute humidity.	

			! Allocate the segment/window temporal arrays.
			ALLOCATE(us(ss))
			ALLOCATE(vs(ss))
			ALLOCATE(ws(ss))
			ALLOCATE(Ts(ss))
			ALLOCATE(qs(ss))

			! Allocate a temporary array for gap filling.
			ALLOCATE(tempa(ss,5))			! 5 variables to gap fill.

			! Allocate space for the temporal variance and covariance arrays.
			ALLOCATE(varu(sspmm))
			ALLOCATE(varv(sspmm))
			ALLOCATE(varw(sspmm))
			ALLOCATE(varT(sspmm))
			ALLOCATE(varq(sspmm))

			ALLOCATE(covuw(sspmm))
			ALLOCATE(covvw(sspmm))
			ALLOCATE(covTw(sspmm))
			ALLOCATE(covqw(sspmm))
			ALLOCATE(ustar(sspmm))
			ALLOCATE(sobukhovl(sspmm))
			ALLOCATE(Tstar(sspmm))	
			ALLOCATE(szmes(sspmm))
	
			! And for the spectral window time stamp array.
			ALLOCATE(windowtime(2,sspmm,6))

			! And finally the mean horizontal wind speed.
			ALLOCATE(meanuh(sspmm))
			
			! Calculate window variances and covariances for controls. -> As with spectra IF >1% of the segment
			! is nans exclude. Otherwise interpolate.
			DO j=1,sspmm
				!ll=(j-1)*ss+1 		! Index lower limit.
				!ul=j*ss		! upper limit.
				ll=(j-1)*iib+1		! Index lower limit.
				ul=ll+ss-1		! Index upper limit
				
				! Array segments.
				us=u(ll:ul)
				vs=v(ll:ul)
				ws=w(ll:ul)
				Ts=T(ll:ul)		
				qs=q(ll:ul)

				! Allocate the start and end of the spectral windows:
				windowtime(1,j,:)=ecv%t(ll,:)
				windowtime(2,j,:)=ecv%t(ul,:)
			
				! Control IF > 1% nan flags in window; if so exclude from spectral analysis.
				ok=0
				isnan=0
				k=1
				DO WHILE(ok.EQ.0)
					IF(us(k).EQ.ecv%nanreal) THEN
						isnan=isnan+1
						IF( ((1.0*isnan)/(1.0*ss)) .GT. 0.1 ) THEN	! IF > 10% nans flag window.
							varu(j)=ecv%nanreal
							varv(j)=ecv%nanreal
							varw(j)=ecv%nanreal
							varT(j)=ecv%nanreal
							varq(j)=ecv%nanreal
							covuw(j)=ecv%nanreal
							covvw(j)=ecv%nanreal
							covTw(j)=ecv%nanreal
							covqw(j)=ecv%nanreal
							sobukhovl(j)=ecv%nanreal 
							ustar(j)=ecv%nanreal
							Tstar(j)=ecv%nanreal
							szmes(j)=ecv%nanreal
							ok=1
						END IF
					END IF
					IF(k.GE.ss) THEN
						ok=1
					END IF
					k=k+1
				END DO
				
				IF(varu(j).NE.ecv%nanreal) THEN
					IF(isnan.NE.0) THEN
						! Interpolate flaged entries (when <1% of the segment) using the gapfill function. 
						tempa(:,1)=us(:)
						tempa(:,2)=vs(:)
						tempa(:,3)=ws(:)
						tempa(:,4)=Ts(:)
						tempa(:,5)=qs(:)
						CALL gapfill(ecv,ss,5,tempa)
						us(:)=tempa(:,1)
						vs(:)=tempa(:,2)
						ws(:)=tempa(:,3)
						Ts(:)=tempa(:,4)
						qs(:)=tempa(:,5)

					END IF	
					tempa(:,:)=0
					

	

					! Compute the window variances and covariances and the mean 'horizontal'
					! wind speed.
					CALL temporalmean(ecv,ss,us,tempub)
					CALL temporalmean(ecv,ss,vs,tempvb)
					meanuh(j)=SQRT(tempub**2+tempvb**2)
					
					! Rotate u into the (streamline plane) mean velocity for the given block,
					! and v into the cross-wind direction. Of course this does NOT change
					! meanuh(j).
					CALL rotatestream(ecv,ss,us,vs)

					IF(ecv%detrend.EQV..TRUE.) THEN	! <-- Check if detrend is specified by the user.
						CALL detrend(ecv,ss,us)
						CALL detrend(ecv,ss,vs)
						CALL detrend(ecv,ss,ws)
						CALL detrend(ecv,ss,Ts)
						CALL detrend(ecv,ss,qs)
					END IF

					! Add a taper to the series if specified by the user.
					IF(ecv%taper.EQV..TRUE.) THEN
						IF(ecv%tapertype.EQ.'B') THEN		! Use the Bell taper.
							DO i=1,NINT(0.1*ss)		! Apply Bell taper @ start of window.
								taper=SIN(5*pi*(i-1)/((1.0)*(ss-1)))**2
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper
							END DO
							DO i=NINT(0.9*ss),ss		! Apply Bell taper @ end of window.
								taper=SIN(5*pi*(i-1)/((1.0)*(ss-1)))**2
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper		
							END DO
						ELSEIF(ecv%tapertype.EQ.'H') THEN	! Use the Hamming window.
							DO i=1,ss
								taper=0.54-0.46*COS(2*pi*(i-1)/(ss-1))
								us(i)=us(i)*taper
								vs(i)=vs(i)*taper
								ws(i)=ws(i)*taper
								Ts(i)=Ts(i)*taper
								qs(i)=qs(i)*taper
							END DO
						ELSE					! Make sure a taper is specified if taperedspectra=TRUE.
							PRINT*,'Please specify taper B or H'
							STOP
						END IF
						! Furthermore make sure to remove the residual mean introduced by the tapering operation.
						CALL temporalmean(ecv,ss,us,tempmean)
						us=us-tempmean
						CALL temporalmean(ecv,ss,vs,tempmean)
						vs=vs-tempmean
						CALL temporalmean(ecv,ss,ws,tempmean)
						ws=ws-tempmean
						CALL temporalmean(ecv,ss,Ts,tempmean)
						Ts=Ts-tempmean
						CALL temporalmean(ecv,ss,qs,tempmean)
						qs=qs-tempmean

						! Also scale the series to account for the variance lost by the tapering using
						! the compensation factors (after Kaimal et al. 1989).
						IF(ecv%tapertype.EQ.'B') THEN
							compf=ecv%cfB
						ELSEIF(ecv%tapertype.EQ.'H') THEN
							compf=ecv%cfH
						END IF
						us=compf*us
						vs=compf*vs
						ws=compf*ws
						Ts=compf*Ts
						qs=compf*qs
						
					END IF
					! NB. The above is important if the fourier analysis is performed on detrended and/or tapered 
					! series; as the statistical and spectral comparisson is only possible for identical series.


					CALL covariance(ecv,ss,us,us,varu(j))
					CALL covariance(ecv,ss,vs,vs,varv(j))
					CALL covariance(ecv,ss,ws,ws,varw(j))
					CALL covariance(ecv,ss,Ts,Ts,varT(j))
					CALL covariance(ecv,ss,qs,qs,varq(j))
					
					CALL covariance(ecv,ss,us,ws,covuw(j))
					CALL covariance(ecv,ss,vs,ws,covvw(j))
					CALL covariance(ecv,ss,Ts,ws,covTw(j))
					CALL covariance(ecv,ss,qs,ws,covqw(j))
				
					! Compute Obukhov length 'manually' to avoid dividing by detrended
					! block averaged sonic temperature which is near zero.
					tempmean=0
					CALL temporalmean(ecv,ss,T(ll:ul),tempmean)
					CALL frictionvelocity(ecv,ss,us,vs,ws,ustar(j))
					sobukhovl(j)=-1.0*(ustar(j)**3)/(0.4*9.81*covTw(j)/tempmean)				
									


					! Compute the obukhov length and ustar for the window.
					!  The friction velocity (ustar).
					!ustar(j)=SQRT(SQRT(covuw(j)**2+covvw(j)**2))
					!  The kinematic buoyancy flux bf.
					!CALL covariance(ecv,ss,vpts,ws,sobukhovl(j))
					! Correct the term above for L=-ustar^3/((k theta_v)*bf/g)
					! where u_star={covuw**2+covvw**2}^{1/4}.
					!CALL temporalmean(ecv,ss,vpts,vptbar)
					!sobukhovl(j)=(-1*ustar(j)**3)/(0.4*9.81*sobukhovl(j)/vptbar)	

					! Calculate the dynamic temperature for the spectral window.
					! Defined: T_*=-cov(Ts,w)/u_*
					CALL covariance(ecv,ss,ecv%Ts(ll:ul),ecv%w(ll:ul),Tstar(j))
					Tstar(j)=-1.0*(Tstar(j))/(ustar(j))
					

					! Calculate the measurement height for the spectral window.
					CALL temporalmean(ecv,ss,ecv%zmes(ll:ul),szmes(j))			

				END IF
				ll=0
				ul=0
				us(:)=0
				vs(:)=0
				ws(:)=0
				Ts(:)=0
				qs(:)=0

			END DO

			! Error tolerance for the spectra.
			etol=0.01!01

			! Allocate the dimensions of the power spectrum arrays; de and re-allocate if already
			! allocated.

			ALLOCATE(su(l,sspmm))
			ALLOCATE(sv(l,sspmm))
			ALLOCATE(sw(l,sspmm))
			ALLOCATE(sT(l,sspmm))
			ALLOCATE(sh(l,sspmm))
			ALLOCATE(siu(l,sspmm))
			ALLOCATE(siv(l,sspmm))
			ALLOCATE(siw(l,sspmm))
			ALLOCATE(siT(l,sspmm))
			ALLOCATE(siq(l,sspmm))


			! Initialize as an arrays of zeros.
			su(:,:)=0
			sv(:,:)=0
			sw(:,:)=0
			sT(:,:)=0
			sh(:,:)=0
			


			! Allocate the dimensions of the spectral variance arrays.
			ALLOCATE(svaru(sspmm))
			ALLOCATE(svarv(sspmm))
			ALLOCATE(svarw(sspmm))
			ALLOCATE(svarT(sspmm))
			ALLOCATE(svarq(sspmm))
			
			! Initialize as an arrays of zeros.
			svaru(:)=0
			svarv(:)=0
			svarw(:)=0
			svarT(:)=0
			svarq(:)=0
		
			! Initialize controls as zeros.
			okvaru=0
			okvarv=0
			okvarw=0
			okvarT=0
			okvarq=0
			isnan=0
			dvar=0
			awins=0

			PRINT*,'-->Starting power spectra'
			! Calculate the power spectrum; the energy spectrum normalized by frequency increment.
			DO j=1,sspmm
				IF(DBLE(fu(1,j)).EQ.ecv%nanreal) THEN	
					su(:,j)=ecv%nanreal
					sw(:,j)=ecv%nanreal
					sT(:,j)=ecv%nanreal
					sh(:,j)=ecv%nanreal

					siu(:,j)=su(:,j)
					siv(:,j)=sv(:,j)
					siw(:,j)=sw(:,j)
					siT(:,j)=sT(:,j)
					siq(:,j)=sh(:,j)

					isnan=isnan+1
				ELSE	
					! Update the number of active windows.
					awins=awins+1
					! The spectral energy at the fundamental frequency (i=1->f=0) does not contribute to the variance;
					! as the norm of this is equal to the mean of the window. Set value to nan
					su(1,j)=ecv%nanreal
					sv(1,j)=ecv%nanreal
					sw(1,j)=ecv%nanreal
					sT(1,j)=ecv%nanreal
					sh(1,j)=ecv%nanreal

					DO i=2,l! Fourier modes at remaining frequenceis contribute to the power spectrum, these are
						! deviations from the mean.

						! Compute the absolute value of the Fourier spectrum.
						su(i,j)=DBLE(fu(i,j))**2+DIMAG(fu(i,j))**2
						sv(i,j)=DBLE(fv(i,j))**2+DIMAG(fv(i,j))**2
						sw(i,j)=DBLE(fw(i,j))**2+DIMAG(fw(i,j))**2
						sT(i,j)=DBLE(fT(i,j))**2+DIMAG(fT(i,j))**2
						sh(i,j)=DBLE(fq(i,j))**2+DIMAG(fq(i,j))**2
	

						! Take into account folding about the Nyquist frequency:
						! i.e. include the equal energy from the corresponding 
						! (symmetric about f_ny) higher frequency.
						su(i,j)=2*su(i,j)
						sv(i,j)=2*sv(i,j)
						sw(i,j)=2*sw(i,j)
						sT(i,j)=2*sT(i,j)
						sh(i,j)=2*sh(i,j)
	
						! If the window length is even then the nyquist frequency is 
						! not folded.
						IF(i.EQ.l) THEN
							su(i,j)=su(i,j)-(1-isodd)*0.5*su(i,j)
							sv(i,j)=sv(i,j)-(1-isodd)*0.5*sv(i,j)
							sw(i,j)=sw(i,j)-(1-isodd)*0.5*sw(i,j)
							sT(i,j)=sT(i,j)-(1-isodd)*0.5*sT(i,j)
							sh(i,j)=sh(i,j)-(1-isodd)*0.5*sh(i,j)
			
						END IF
						! Add the spectral intensity to the total spectral window variance.
						svaru(j)=svaru(j)+su(i,j)
						svarv(j)=svarv(j)+sv(i,j)
						svarw(j)=svarw(j)+sw(i,j)
						svarT(j)=svarT(j)+sT(i,j)
						svarq(j)=svarq(j)+sh(i,j)

						! Save the spectral intensities for the variance Ogives.
						siu(i,j)=su(i,j)
						siv(i,j)=sv(i,j)
						siw(i,j)=sw(i,j)
						siT(i,j)=sT(i,j)
						siq(i,j)=sh(i,j)		

						! Normalize the spectral intensity by the frequency increment to get the
						! power spectrum (aka. spectral density).
						su(i,j)=su(i,j)/(freq(i)-freq(i-1))
						sv(i,j)=sv(i,j)/(freq(i)-freq(i-1))
						sw(i,j)=sw(i,j)/(freq(i)-freq(i-1))
						sT(i,j)=sT(i,j)/(freq(i)-freq(i-1))
						sh(i,j)=sh(i,j)/(freq(i)-freq(i-1))


					END DO
				
					! Compute the normalized differences between the spectral and temporal variances.

					! Compare the u-variances computed in the temporal and spectral domain:
					dvar=(svaru(j)-varu(j))/varu(j)
					IF(ABS(dvar).LE.etol) THEN
						okvaru=okvaru+1
					END IF
					dvar=0

					! Compare the v-variances computed in the temporal and spectral domain:
					dvar=(svarv(j)-varv(j))/varv(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarv=okvarv+1
					END IF
					dvar=0


					! Compare the w-variances computed in the temporal and spectral domain:
					dvar=(svarw(j)-varw(j))/varw(j)
					! Compare the variances computed in the temporal and spectral domain:
					IF(ABS(dvar).LE.etol) THEN
						okvarw=okvarw+1
					END IF
					dvar=0

					! Compare the T-variances computed in the temporal and spectral domain:
					dvar=(svarT(j)-varT(j))/varT(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarT=okvarT+1
					END IF
					dvar=0


					! Compare the Q-variances computed in the temporal and spectral domain:
					dvar=(svarq(j)-varq(j))/varq(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarq=okvarq+1
					END IF
					dvar=0
				END IF
			END DO
			!STOP

			! Print the number of windows where spectral and temporal variances matched for both series.
			!5 to 3 comment and change back
			PRINT '("------Active spectral windows:",I5)',awins
			PRINT '("------Matched u-variances for", f8.2, " % of active windows")',(100.0*okvaru)/(1*(sspmm-isnan))
			PRINT '("------Matched v-variances for", f8.2, " % of active windows")',(100.0*okvarv)/(1*(sspmm-isnan))
			PRINT '("------Matched w-variances for", f8.2, " % of active windows")',(100.0*okvarw)/(1*(sspmm-isnan))
			PRINT '("------Matched T-variances for", f8.2, " % of active windows")',(100.0*okvarT)/(1*(sspmm-isnan))
			PRINT '("------Matched h-variances for", f8.2, " % of active windows")',(100.0*okvarq)/(1*(sspmm-isnan))

			WRITE(ecv%loglun,*)"------Active spectral windows:",awins
			WRITE(ecv%loglun,*)"------Matched u-variances for", (100.0*okvaru)/(1*(sspmm-isnan)), " % of active windows"
			WRITE(ecv%loglun,*)"------Matched v-variances for", (100.0*okvarv)/(1*(sspmm-isnan)), " % of active windows"
			WRITE(ecv%loglun,*)"------Matched w-variances for", (100.0*okvarw)/(1*(sspmm-isnan)), " % of active windows"
			WRITE(ecv%loglun,*)"------Matched T-variances for", (100.0*okvarT)/(1*(sspmm-isnan)), " % of active windows"
			WRITE(ecv%loglun,*)"------Matched h-variances for", (100.0*okvarq)/(1*(sspmm-isnan)), " % of active windows"

			!STOP

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Cross-spectral calculations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! Set the dimensions of the ogive arrays. [360 logarithmically evenspaced entries].
			nfcog=90
			! Allocate the dimensions of the Ogive frequency array.
			ALLOCATE(fogs(nfcog))
			! And of the array containing the corresponding indices in the freq array.
			ALLOCATE(findi(nfcog))
			! Calculate the constant logarithmic frequency increment
			deltaog=LOG10(freq(l)/freq(2))*((1.0)/(1.0*nfcog-1.0))
			! Allocate the frequencies:
			DO i=1,nfcog
				fogs(i)=freq(2)*10**((i-1)*deltaog)
			END DO
				


			! Allocate dimensions of the cospectrum, quadrature and ogive arrays.
			ALLOCATE(couw(l,sspmm))
			ALLOCATE(ciuw(l,sspmm))
			ALLOCATE(quuw(l,sspmm))
			ALLOCATE(oguw(nfcog,sspmm))!ALLOCATE(oguw(l,sspmm))
			
			ALLOCATE(covw(l,sspmm))
			ALLOCATE(civw(l,sspmm))
			ALLOCATE(quvw(l,sspmm))
			ALLOCATE(ogvw(nfcog,sspmm))!ALLOCATE(ogvw(l,sspmm))
			
			ALLOCATE(coTw(l,sspmm))
			ALLOCATE(ciTw(l,sspmm))
			ALLOCATE(quTw(l,sspmm))
			ALLOCATE(ogTw(nfcog,sspmm))!ALLOCATE(ogTw(l,sspmm))
			
			ALLOCATE(cohw(l,sspmm))
			ALLOCATE(ciqw(l,sspmm))
			ALLOCATE(quqw(l,sspmm))
			ALLOCATE(oghw(nfcog,sspmm))!ALLOCATE(oghw(l,sspmm))

			ALLOCATE(oguu(nfcog,sspmm))!ALLOCATE(oguu(l,sspmm))
			ALLOCATE(ogvv(nfcog,sspmm))!ALLOCATE(ogvv(l,sspmm))
			ALLOCATE(ogww(nfcog,sspmm))!ALLOCATE(ogww(l,sspmm))
			ALLOCATE(ogTT(nfcog,sspmm))!ALLOCATE(ogTT(l,sspmm))
			ALLOCATE(oghh(nfcog,sspmm))!ALLOCATE(oghh(l,sspmm))

			
			PRINT*,'-->Starting cospectra'
			! Loop over windows.
			DO j=1,sspmm
				IF(DBLE(fu(1,j)).EQ.ecv%nanreal) THEN 	! If 'nan' flag then set crosspectra to nan.
					couw(:,j)=ecv%nanreal
					covw(:,j)=ecv%nanreal
					coTw(:,j)=ecv%nanreal
					cohw(:,j)=ecv%nanreal

					ciuw(:,j)=ecv%nanreal
					civw(:,j)=ecv%nanreal
					ciTw(:,j)=ecv%nanreal
					ciqw(:,j)=ecv%nanreal

					quuw(:,j)=ecv%nanreal
					quvw(:,j)=ecv%nanreal
					quTw(:,j)=ecv%nanreal
					quqw(:,j)=ecv%nanreal
		
					oguw(:,j)=ecv%nanreal
					ogvw(:,j)=ecv%nanreal
					ogTw(:,j)=ecv%nanreal
					oghw(:,j)=ecv%nanreal

					oguu(:,j)=ecv%nanreal
					ogvv(:,j)=ecv%nanreal
					ogww(:,j)=ecv%nanreal
					ogTT(:,j)=ecv%nanreal
					oghh(:,j)=ecv%nanreal
			


				ELSE
					DO i=2,l
						! Intensities.
						ciuw(i,j)=DBLE(fu(i,j))*DBLE(fw(i,j))+DIMAG(fu(i,j))*DIMAG(fw(i,j))	! Cospectrum.
						quuw(i,j)=DIMAG(fu(i,j))*DBLE(fw(i,j))-DBLE(fu(i,j))*DIMAG(fw(i,j))	! Quadrature.
						oguw(:,j)=0								! Initialize Ogive as zeros.

						! Convert to densities.
						couw(i,j)=ciuw(i,j)/(freq(i)-freq(i-1))
						quuw(i,j)=quuw(i,j)/(freq(i)-freq(i-1))


						! Intensities.
						civw(i,j)=DBLE(fv(i,j))*DBLE(fw(i,j))+DIMAG(fv(i,j))*DIMAG(fw(i,j))	! Cospectrum.
						quvw(i,j)=DIMAG(fv(i,j))*DBLE(fw(i,j))-DBLE(fv(i,j))*DIMAG(fw(i,j))	! Quadrature.
								

						! Convert to densities.
						covw(i,j)=civw(i,j)/(freq(i)-freq(i-1))
						quvw(i,j)=quvw(i,j)/(freq(i)-freq(i-1))


						! Intensities.
						ciTw(i,j)=DBLE(fT(i,j))*DBLE(fw(i,j))+DIMAG(fT(i,j))*DIMAG(fw(i,j))	! Cospectrum.
						quTw(i,j)=DIMAG(fT(i,j))*DBLE(fw(i,j))-DBLE(fT(i,j))*DIMAG(fw(i,j))	! Quadrature.

						! Convert to densities.
						coTw(i,j)=ciTw(i,j)/(freq(i)-freq(i-1))
						quTw(i,j)=quTw(i,j)/(freq(i)-freq(i-1))


						! Intensities.
						ciqw(i,j)=DBLE(fq(i,j))*DBLE(fw(i,j))+DIMAG(fq(i,j))*DIMAG(fw(i,j))	! Cospectrum.
						quqw(i,j)=DIMAG(fq(i,j))*DBLE(fw(i,j))-DBLE(fq(i,j))*DIMAG(fw(i,j))	! Quadrature.

						! Convert to densities.
						cohw(i,j)=ciqw(i,j)/(freq(i)-freq(i-1))
						quqw(i,j)=quqw(i,j)/(freq(i)-freq(i-1))

						! Fold the cospectra and the quadratures as was doen for the variance spectra.
						couw(i,j)=2*couw(i,j)
						covw(i,j)=2*covw(i,j)
						coTw(i,j)=2*coTw(i,j)
						cohw(i,j)=2*cohw(i,j)

						quuw(i,j)=2*quuw(i,j)
						quvw(i,j)=2*quvw(i,j)
						quTw(i,j)=2*quTw(i,j)
						quqw(i,j)=2*quqw(i,j)

	
						! If the window length is even then the nyquist frequency is 
						! not folded.
						IF(i.EQ.l) THEN
							couw(i,j)=couw(i,j)-(1-isodd)*0.5*couw(i,j)
							covw(i,j)=covw(i,j)-(1-isodd)*0.5*covw(i,j)
							coTw(i,j)=coTw(i,j)-(1-isodd)*0.5*coTw(i,j)
							cohw(i,j)=cohw(i,j)-(1-isodd)*0.5*cohw(i,j)

							quuw(i,j)=quuw(i,j)-(1-isodd)*0.5*quuw(i,j)
							quvw(i,j)=quvw(i,j)-(1-isodd)*0.5*quvw(i,j)
							quTw(i,j)=quTw(i,j)-(1-isodd)*0.5*quTw(i,j)
							quqw(i,j)=quqw(i,j)-(1-isodd)*0.5*quqw(i,j)
						END IF

					END DO
		
					! Initialize Ogives as zeros.
					oguw(:,j)=0
					ogvw(:,j)=0
					ogTw(:,j)=0
					oghw(:,j)=0
					oguu(:,j)=0
					ogvv(:,j)=0
					ogww(:,j)=0
					ogTT(:,j)=0
					oghh(:,j)=0


					! Set the initial value to nan as its undefined.
					couw(1,j)=ecv%nanreal
					quuw(1,j)=ecv%nanreal
					covw(1,j)=ecv%nanreal
					quvw(1,j)=ecv%nanreal
					coTw(1,j)=ecv%nanreal
					quTw(1,j)=ecv%nanreal
					cohw(1,j)=ecv%nanreal
					quqw(1,j)=ecv%nanreal

				END IF
			END DO
			
			! Compute the Ogives.
			
			! First initialize the covariance controls to be zero.
			okcovuw=0
			okcovvw=0
			okcovTw=0
			okcovqw=0	
			okvaru=0
			okvarv=0
			okvarw=0
			okvarT=0
			okvarq=0
			! Loop over windows.

			! Find the indices of the high resolution frequncies corresponding to the ogive frequencies.
			findi(1)=2	! The first and last ones are known a priori.
			findi(nfcog)=l
			DO i=2,nfcog-1	! Loop through Ogive frequencies.
				DO j=findi(i-1),l	! Loop through indices in high res frequency array, starting at previous index.
					IF(freq(j-1).LT.fogs(i).AND.freq(j).GT.fogs(i)) THEN
						findi(i)=j
						EXIT	! Stop looping when the index is found.
					END IF
				END DO
			END DO	
	
			DO j=1,sspmm
				IF(DBLE(fu(1,j)).NE.ecv%nanreal) THEN
					DO i=1,nfcog-1	! Loop through evenly spaced logarithmic frequencies
						oguw(i,j)=2*SUM(ciuw(findi(i):l,j),1)-(1-isodd)*ciuw(l,j)
						ogvw(i,j)=2*SUM(civw(findi(i):l,j),1)-(1-isodd)*civw(l,j)
						ogTw(i,j)=2*SUM(ciTw(findi(i):l,j),1)-(1-isodd)*ciTw(l,j)
						oghw(i,j)=2*SUM(ciqw(findi(i):l,j),1)-(1-isodd)*ciqw(l,j)
						oguu(i,j)=SUM(siu(findi(i):l,j),1) ! <-- Already folded.
						ogvv(i,j)=SUM(siv(findi(i):l,j),1)
						ogww(i,j)=SUM(siw(findi(i):l,j),1)
						ogTT(i,j)=SUM(siT(findi(i):l,j),1)
						oghh(i,j)=SUM(siq(findi(i):l,j),1)
					END DO
					! Last entry is zero by definition:
					oguw(nfcog,j)=0
					ogvw(nfcog,j)=0
					ogTw(nfcog,j)=0
					oghw(nfcog,j)=0
					oguu(nfcog,j)=0
					ogvv(nfcog,j)=0
					ogww(nfcog,j)=0
					ogTT(nfcog,j)=0
					oghh(nfcog,j)=0

					! Control the uw-covariance.
					dcov=0
					dcov=(oguw(1,j)-covuw(j))/covuw(j)
					IF(ABS(dcov).LE.etol) THEN
						okcovuw=okcovuw+1
					END IF
					dcov=0
		
					! Control the vw-covariance.
					dcov=(ogvw(1,j)-covvw(j))/covvw(j)
					IF(ABS(dcov).LE.etol) THEN
						okcovvw=okcovvw+1
					END IF
					dcov=0

					! Control the Tw-covariance.
					dcov=(ogTw(1,j)-covTw(j))/covTw(j)
					IF(ABS(dcov).LE.etol) THEN
						okcovTw=okcovTw+1
					END IF
					dcov=0

					! Control the qw-covariance.
					dcov=(oghw(1,j)-covqw(j))/covqw(j)
					IF(ABS(dcov).LE.etol) THEN
						okcovqw=okcovqw+1
					END IF
					dcov=0
		
					! Control the u-variance using the Ogive.
					dvar=(oguu(1,j)-varu(j))/varu(j)
					IF(ABS(dvar).LE.etol) THEN
						okvaru=okvaru+1
					END IF
					dvar=0
	
					! Control the v-variance using the Ogive.
					dvar=(ogvv(1,j)-varv(j))/varv(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarv=okvarv+1
					END IF
					dvar=0
			
					! Control the w-variance using the Ogive.
					dvar=(ogww(1,j)-varw(j))/varw(j)
					IF(ABS(dvar).LE.etol) THEN	
						okvarw=okvarw+1
					END IF
					dvar=0
				
					! Control the T-variance using the Ogive.
					dvar=(ogTT(1,j)-varT(j))/varT(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarT=okvarT+1
					END IF
					dvar=0
			
					! Control the q-variance using the Ogive.
					dvar=(oghh(1,j)-varq(j))/varq(j)
					IF(ABS(dvar).LE.etol) THEN
						okvarq=okvarq+1
					END IF
					dvar=0
					

				END IF
			END DO
		

			! Print the number of windows where spectral and temporal variances matched for both series.
			PRINT '("------Matched uw-covariances for" f8.2" % of active windows")',(100.0*okcovuw)/(sspmm-isnan)
			PRINT '("------Matched vw-covariances for" f8.2" % of active windows")',(100.0*okcovvw)/(sspmm-isnan)
			PRINT '("------Matched Tw-covariances for" f8.2" % of active windows")',(100.0*okcovTw)/(sspmm-isnan)
			PRINT '("------Matched hw-covariances for" f8.2" % of active windows")',(100.0*okcovqw)/(sspmm-isnan)
			WRITE(ecv%loglun,*)"------Matched uw-covariances for",(100.0*okcovuw)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched vw-covariances for",(100.0*okcovvw)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched Tw-covariances for",(100.0*okcovTw)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched hw-covariances for",(100.0*okcovqw)/(sspmm-isnan)," % of active windows"
			PRINT '("------Matched u-variances (Ogive) for" f8.2" % of active windows")',(100.0*okvaru)/(sspmm-isnan)
			PRINT '("------Matched v-variances (Ogive) for" f8.2" % of active windows")',(100.0*okvarv)/(sspmm-isnan)
			PRINT '("------Matched w-variances (Ogive) for" f8.2" % of active windows")',(100.0*okvarw)/(sspmm-isnan)
			PRINT '("------Matched T-variances (Ogive) for" f8.2" % of active windows")',(100.0*okvarT)/(sspmm-isnan)
			PRINT '("------Matched h-variances (Ogive) for" f8.2" % of active windows")',(100.0*okvarq)/(sspmm-isnan)
			WRITE(ecv%loglun,*)"------Matched u-variances (Ogive) for",(100.0*okvaru)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched v-variances (Ogive) for",(100.0*okvarv)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched w-variances (Ogive) for",(100.0*okvarw)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched T-variances (Ogive) for",(100.0*okvarT)/(sspmm-isnan)," % of active windows"
			WRITE(ecv%loglun,*)"------Matched h-variances (Ogive) for",(100.0*okvarq)/(sspmm-isnan)," % of active windows"


			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!! Smoothing and discarding nans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


			! Scan through the ogives and find the non-nan entries in EACH non-nan window.
			ii=0
			jj=0	
			k=0	
			! Loop through windows.	
			DO j=1,SIZE(oguw,2)
				! Check if the window is active.
				IF(oguw(2,j).NE.ecv%nanreal) THEN
					jj=jj+1	!  Window is active; update counter.
				END IF
			END DO


			! Now jj will be the number of active windows for all the spectra, while ii is the number of
			! active frequencies for the ogives. 



			! Dealiasing power spectra and cospectra.




			! We may now allocate space for the spectra that are fed to the ecv structure.
			PRINT*,'--> Restructuring spectral arrays, removing nan entries and Konno Ohmachi smoothing to save disk space'
			WRITE(ecv%loglun,*)'--> Restructuring spectral arrays, removing nan entries and Konno Ohmachi smoothing to save disk space'



			
			! Set the bandwidth.
			bwidth=10!10
			! Standard spectra:
			! Number of frequencies to smooth:
			nf=SIZE(su,1)-1		! Excluding f=0.
			! Calculate the corresponding number of frequencies in the smoothed frequency array.
			nfc=90
			! Allocate the central frequency array
			ALLOCATE(fcen(nfc))
			! Allocate the Konno_Ohmachi smoothing matrix.
			ALLOCATE(KOSM(nfc,nf))
			! Compute the Konno-Ohmachi smoothing matrix, KOSM, with nfc smoothing windows.
			CALL	kosmoothing(ecv,freq(2:nf+1),nf,fcen,nfc,bwidth,KOSM)

	


			IF(awins.EQ.0)	THEN	! <-- If no active spectral windows write dummy spectra to the output file.
				ecv%spactive=.FALSE.
				RETURN		! Break out of the subroutine.
			ELSE			! <-- Windows are active; Allocate space for these.
				IF(ecv%smooth.EQV..FALSE.) THEN
					DEALLOCATE(fcen)
					nfc=nf
					ALLOCATE(fcen(nfc))
					fcen=freq(2:nf+1)
				END IF
					
			
				ecv%spactive=.TRUE.
				! Allocate the dimensions:
				ALLOCATE(ecv%tspwin(2,jj,SIZE(ecv%t,2)))

				! Allocate the mean wind speed array.
				ALLOCATE(ecv%meanuh(jj))

				! Exclude the mean entry (f=0) for the spectra; as this is undefined.
				ALLOCATE(ecv%su(nfc,jj))
				ALLOCATE(ecv%sv(nfc,jj))
				ALLOCATE(ecv%sw(nfc,jj))
				ALLOCATE(ecv%sT(nfc,jj))
				ALLOCATE(ecv%sh(nfc,jj))

				ALLOCATE(ecv%couw(nfc,jj))
				ALLOCATE(ecv%covw(nfc,jj))
				ALLOCATE(ecv%coTw(nfc,jj))
				ALLOCATE(ecv%cohw(nfc,jj))

				ALLOCATE(ecv%quuw(nfc,jj))
				ALLOCATE(ecv%quvw(nfc,jj))
				ALLOCATE(ecv%quTw(nfc,jj))
				ALLOCATE(ecv%quqw(nfc,jj))

				! Allocate the Ogives that are sent to the module.
				ALLOCATE(ecv%oguw(nfcog,jj))
				ALLOCATE(ecv%ogvw(nfcog,jj))
				ALLOCATE(ecv%ogTw(nfcog,jj))
				ALLOCATE(ecv%oghw(nfcog,jj))
				ALLOCATE(ecv%oguu(nfcog,jj))
				ALLOCATE(ecv%ogvv(nfcog,jj))
				ALLOCATE(ecv%ogww(nfcog,jj))
				ALLOCATE(ecv%ogTT(nfcog,jj))
				ALLOCATE(ecv%oghh(nfcog,jj))
				ALLOCATE(ecv%fog(nfcog))
				ecv%fog=fogs
	
				

				! Allocate the dimensions of the spectral frequency array.
				ALLOCATE(ecv%fsp(nfc))

				! Allocate the spectral frequency.
				ecv%fsp=fcen	

				! And finally the dimensions of the obukhov length for the spectra.
				ALLOCATE(ecv%sobukhovl(jj))
				! and ustar.
				ALLOCATE(ecv%sustar(jj))
				! and Tstar.
				ALLOCATE(ecv%sTstar(jj))
				! as well as the measurement height for the spectra.
				ALLOCATE(ecv%sz(jj))
			END IF
					

			! Allocate the entries if at least one window is active:
			j=0	! <-- Active window counter.
			k=0	! <-- Actual window counter.
			IF(ecv%spactive.EQV..TRUE.) THEN	! <-- If statement not needed, but just to be safe to avoid allocation errors.
				DO WHILE(j.LT.jj)
					k=k+1	! <-- Update window counter.
					IF(oguw(2,k).NE.ecv%nanreal) THEN	! <-- Window is active.
						j=j+1				! <-- Update active window counter.

						! Set the number of frequencies resolved:
						ul=SIZE(su,1)
	
						! Allocate the spectral time stamps.
						ecv%tspwin(:,j,:)=windowtime(:,k,:)
	
						! Allocate spectral obukhov length.
						ecv%sobukhovl(j)=sobukhovl(k)	

						! Allocate spectral friction velocity.
						ecv%sustar(j)=ustar(k)
			
						! Allocate spectral dynamic temperature.
						ecv%sTstar(j)=Tstar(k)

						! Allocate the spectral measuring height.
						ecv%sz(j)=szmes(k)

						! Allocate the mean 'horizontal' wind speed for the active window.
						ecv%meanuh(j)=meanuh(k)	

						! If smoothing is specified then smooth and dealias the spectra:
						IF(ecv%smooth.EQV..TRUE.) THEN		
							! Allocate spectra.
							! § Remembering to exclude the first 'mean' entry.
							! Dealiasing before smoothing.
							CALL dealias(ecv,nf,freq(2:nf+1),su(2:ul,k))					
							ecv%su(:,j)=MATMUL(KOSM,su(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),sv(2:ul,k))	
							ecv%sv(:,j)=MATMUL(KOSM,sv(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),sw(2:ul,k))	
							ecv%sw(:,j)=MATMUL(KOSM,sw(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),sT(2:ul,k))	
							ecv%sT(:,j)=MATMUL(KOSM,sT(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),sh(2:ul,k))	
							ecv%sh(:,j)=MATMUL(KOSM,sh(2:ul,k))

							CALL dealias(ecv,nf,freq(2:nf+1),couw(2:ul,k))	
							ecv%couw(:,j)=MATMUL(KOSM,couw(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),covw(2:ul,k))	
							ecv%covw(:,j)=MATMUL(KOSM,covw(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),coTw(2:ul,k))	
							ecv%coTw(:,j)=MATMUL(KOSM,coTw(2:ul,k))
							CALL dealias(ecv,nf,freq(2:nf+1),cohw(2:ul,k))	
							ecv%cohw(:,j)=MATMUL(KOSM,cohw(2:ul,k))
			
							ecv%quuw(:,j)=MATMUL(KOSM,quuw(2:ul,k))
							ecv%quvw(:,j)=MATMUL(KOSM,quvw(2:ul,k))
							ecv%quTw(:,j)=MATMUL(KOSM,quTw(2:ul,k))
							ecv%quqw(:,j)=MATMUL(KOSM,quqw(2:ul,k))
						ELSE
							ecv%su(:,j)=su(2:ul,k)
							ecv%sv(:,j)=sv(2:ul,k)
							ecv%sw(:,j)=sw(2:ul,k)
							ecv%sT(:,j)=sT(2:ul,k)
							ecv%sh(:,j)=sh(2:ul,k)

							ecv%couw(:,j)=couw(2:ul,k)
							ecv%covw(:,j)=covw(2:ul,k)
							ecv%coTw(:,j)=coTw(2:ul,k)
							ecv%cohw(:,j)=cohw(2:ul,k)

							ecv%quuw(:,j)=quuw(2:ul,k)
							ecv%quvw(:,j)=quvw(2:ul,k)
							ecv%quTw(:,j)=quTw(2:ul,k)
							ecv%quqw(:,j)=quqw(2:ul,k)
						END IF

						ecv%oguw(:,j)=oguw(:,k)
						ecv%ogvw(:,j)=ogvw(:,k)
						ecv%ogTw(:,j)=ogTw(:,k)
						ecv%oghw(:,j)=oghw(:,k)
						ecv%oguu(:,j)=oguu(:,k)
						ecv%ogvv(:,j)=ogvv(:,k)
						ecv%ogww(:,j)=ogww(:,k)
						ecv%ogTT(:,j)=ogTT(:,k)
						ecv%oghh(:,j)=oghh(:,k)
					
					END IF
				END DO
			END IF
				
				


	
			PRINT*,'End of spectral analysis'
			PRINT*,'*******************************************************************************'
			WRITE(ecv%loglun,*)'End of spectral analysis'
			WRITE(ecv%loglun,*)'*******************************************************************************'

			! Deallocate temporary arrays.
			IF(ecv%spactive.EQV..TRUE.) THEN	! <-- If statement not needed, but just to be safe to avoid allocation errors.
				DEALLOCATE(windowtime)
				DEALLOCATE(us)
				DEALLOCATE(vs)
				DEALLOCATE(ws)
				DEALLOCATE(Ts)
				DEALLOCATE(qs)
				DEALLOCATE(fu)
				DEALLOCATE(fv)
				DEALLOCATE(fw)
				DEALLOCATE(fT)
				DEALLOCATE(fq)
				DEALLOCATE(su)
				DEALLOCATE(sv)
				DEALLOCATE(sw)
				DEALLOCATE(sT)
				DEALLOCATE(sh)
				DEALLOCATE(couw)
				DEALLOCATE(covw)
				DEALLOCATE(coTw)
				DEALLOCATE(cohw)
				DEALLOCATE(ciuw)
				DEALLOCATE(civw)
				DEALLOCATE(ciTw)
				DEALLOCATE(ciqw)
				DEALLOCATE(quuw)
				DEALLOCATE(quvw)
				DEALLOCATE(quTw)
				DEALLOCATE(quqw)
				DEALLOCATE(oguw)
				DEALLOCATE(ogvw)
				DEALLOCATE(ogTw)
				DEALLOCATE(oghw)
				DEALLOCATE(oguu)
				DEALLOCATE(ogvv)
				DEALLOCATE(ogww)
				DEALLOCATE(ogTT)
				DEALLOCATE(oghh)
				DEALLOCATE(freq)
				DEALLOCATE(varu)
				DEALLOCATE(varv)
				DEALLOCATE(varw)
				DEALLOCATE(varT)
				DEALLOCATE(varq)
				DEALLOCATE(svaru)
				DEALLOCATE(svarv)
				DEALLOCATE(svarw)
				DEALLOCATE(svarT)
				DEALLOCATE(svarq)
				DEALLOCATE(covuw)
				DEALLOCATE(covvw)
				DEALLOCATE(covTw)
				DEALLOCATE(covqw)
				DEALLOCATE(meanuh)
				DEALLOCATE(tempa)
				DEALLOCATE(fogs)
				DEALLOCATE(fcen)
				DEALLOCATE(KOSM)
				DEALLOCATE(sobukhovl)
				DEALLOCATE(szmes)
				DEALLOCATE(ustar)
				DEALLOCATE(Tstar)
				DEALLOCATE(findi)
			
			END IF

			RETURN
		END SUBROUTINE

		! Calculates and applies the digital low pass filter to an input spectrm to reduce aliasing. 
		! The filter has the form proposed by Gobbi et al. (2006).
		SUBROUTINE dealias(ecv,nf,f,spec)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER,INTENT(IN)			:: nf
			REAL*8,DIMENSION(nf),INTENT(IN)		:: f
			REAL*8,DIMENSION(nf),INTENT(INOUT)	:: spec
			REAL*8,DIMENSION(nf)			:: Hlp
			REAL*8					:: pi,dt
			INTEGER					:: i
	
			pi=4.0*ATAN(1.0)
			dt=ecv%dt
			DO i=1,nf
				Hlp(i)=((1.0+COS(pi*f(i)*dt))**2)/4.0
				spec(i)=Hlp(i)*spec(i)
			END DO
			
			RETURN
		END SUBROUTINE
			

		! Computes a smoothing matrix, warray, given a discrete frequency resolution set by the
		! frequency array f (length nf) based on the central frequencies in the fc array. 
		! The fc array is output by the subroutine and is defined as a logarithmically spaced
		! array with length nfc. Constrained to satisfy fc(1)=f(1) and fc(end)=f(end).
		! The parameter b sets the bandwidth of the smoothing windows defined
		! by the smoothing matrix warray. A low b (order 1) results in an overly smoothed spectra while
		! a high b (order >100) hardly smoothes the spectra at all, a value of around 50 is recommended.
		! This method was first proposed by Konno and Ohmachi (1998)
		! to improve presentation of spectra in loglog type plots. Here we also use the method
		! to reduce the dimensions (from nf to nfc) and reduce noise with hardly any loss of information. Only needs
		! to be called once at the end of the cross-spectral calculations. The smoothing matrix
		! is the same for all the spectra.
		SUBROUTINE	kosmoothing(ecv,f,nf,fc,nfc,b,warray)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)		:: ecv
			INTEGER,INTENT(IN)			:: nf,nfc		! Dimensions (input).
			REAL*8,INTENT(IN)			:: b			! Bandwidth.
			REAL*8,DIMENSION(nf),INTENT(IN)		:: f			! Input frequencies.
			REAL*8,DIMENSION(nfc),INTENT(OUT)	:: fc			! Input central frequencies, corresponds
											! to the resolution of the smoothed spectra.
			REAL*8,DIMENSION(nfc,nf),INTENT(OUT)	:: warray		! Output smoothing matrix.
			REAL*8					:: temp			! Temporary real.
			REAL*8					:: sumi			! Sum of the components in the ith window.
			REAL*8					:: delta		! Constant logarithmic increment.
			INTEGER					:: i,j			! Counters.
		
			! Define the central frequency array onto which the spectra are to be smoothed.
			delta=LOG10(f(nf)/f(1))/(1.0*(nfc-1))				! Compute the logarithmic increment:  												! log10[fc(i+1)]-log10[fc(i)]=delta for all i<nfc.
			DO i=1,nfc	! Allocate the central frequencies. Given by f(1)*10^((i-1)*delta)
				fc(i)=f(1)*10**((i-1)*delta)
			END DO
				

			sumi=0
			DO i=1,nfc
				DO j=1,nf
					IF(f(j).EQ.fc(i)) THEN
						warray(i,j)=1	! Avoid division by zero. Should be order unity when frequencies are coincident.
					ELSE
						temp=(1.0*b)*LOG10(f(j)/fc(i))			
						warray(i,j)=(SIN(temp)/temp)**4		! Compute the K&O smoothing functions
											! value for the i-th smoothing window at the j-th frequency.
					END IF
					sumi=sumi+warray(i,j)				! Update the sum of the ith window.
				END DO
				! Normalize window i. 
				warray(i,:)=warray(i,:)/sumi
				sumi=0							! Reset the sum for next window.
			END DO
			
			RETURN
		END SUBROUTINE
				
						
		
			
			
			



	
		!----------------------------------------------------------------------------------------------------------------------
		! Block calculations.
		!----------------------------------------------------------------------------------------------------------------------	

		! Calculations of the block fractional flux sampling error for the turbulent fluxes of sensible heat, water vapor and momentum
		! based on the method of Finkelstein and Sims (2001).
		SUBROUTINE 	blockffse(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER				:: i,j,k,ll,ul,n,ii,jj		! Counters and bounds.
			INTEGER				:: lagmax			! Maximum lag to shift (+/-).
			REAL*8				:: ffsetsw,ffsehw		! Temporary reals to hold the fractional flux sampling error for each block.
			REAL*8				:: covtsw,covhw			! Temporary reals to hold the covariances for each block.
			REAL*8				:: varw,varts,varh		! Temporary reals to hold the variances for each block.
			REAL*8,ALLOCATABLE,DIMENSION(:)	:: hs,ws,Tss			! Temporary arrays to hold the samples in each block.
			REAL*8,ALLOCATABLE,DIMENSION(:,:):: tempa			! Temporary array for gap filling.
			REAL*8				:: atst,awt,aht			! Temporary reals to hold the auto covariances at a given lag.
			REAL*8				:: ctswt,cwtst,chwt,cwht	! Temporary reals to hold the cross covariances at a given lag.
			REAL*8				:: hbar,tsbar,wbar		! Temporary reals to hold the block average.
			INTEGER				:: nnans			! Number of nan entries in the block.
			REAL*8				:: covtsw0,covhw0		! Temporary reals to hold the non LDT covariances.
			
			PRINT*,'*******************************************************************************'
			PRINT*,'starting ffse routine'		
			WRITE(ecv%loglun,*)'*******************************************************************************'
			WRITE(ecv%loglun,*)'starting ffse routine.'

			! Set the maximum lag index: i.e. the number of indices to shift the series against each other (or themselves).
			! F&S recommend at least 200 at 10 Hz. At 20 Hz we use 400 which corresponds to a 20 second shift. Nonetheless
			! in exploiting the symmetry of their expression we cut the number of computations in half [from 800 to 400
			! cross-covariance and auto-covariance estimates]. Despite this subroutine being short it is quite computationally
			! expensive and may be turned off if no error estimates are desired.
			lagmax=400

			! Set dimensions.
			ii=ecv%blockl			! Number of points in a given block.		
			jj=FLOOR(1.0*ecv%nn/ii)		! Number of blocks of data.

			! Allocate temporary arrays and the FFSEs
			ALLOCATE(hs(ii))
			ALLOCATE(ws(ii))
			ALLOCATE(Tss(ii))
			ALLOCATE(ecv%ffseTsw(jj))
			ALLOCATE(ecv%ffsehw(jj))
			ALLOCATE(tempa(ii,3))

			! Commence flagging of the respective blocks.	
			IF(ALLOCATED(ecv%flag).EQV..TRUE.) THEN
				DO j=1,jj	! Block loop.
					IF(ecv%flag(j).LT.2) THEN
						! Reset statistics on each pass.
						ffseTsw=0
						ffsehw=0
						covTsw=0
						covhw=0
						covTsw0=0
						covhw0=0
						varTs=0
						varh=0
						varw=0
						hbar=0
						Tsbar=0
						wbar=0
						nnans=0
						tempa(:,:)=0
						
						
							
						! Index bounds for the block.
						ll=(j-1)*ii+1
						ul=j*ii	
			
						! Set the values of the block arrays.
						hs=ecv%h(ll:ul)
						Tss=ecv%Ts(ll:ul)
						ws=ecv%w(ll:ul)
						
						! Interpolate flaged entries using the gapfill function. 
						tempa(:,1)=hs(:)
						tempa(:,2)=Tss(:)
						tempa(:,3)=ws(:)
						CALL gapfill(ecv,ii,3,tempa)
						hs(:)=tempa(:,1)
						Tss(:)=tempa(:,2)
						ws(:)=tempa(:,3)

						! Compute the non-LDT covariances.
						CALL covariance(ecv,ii,hs,ws,covhw0)
						CALL covariance(ecv,ii,Tss,ws,covTsw0)

						! Important to detrend these; particularly when computing the autocovariance
						! as it is not well defined in the presence of even a small trend.
						CALL detrend(ecv,ii,hs)
						CALL detrend(ecv,ii,Tss)
						CALL detrend(ecv,ii,ws)

						! Calculate the block averages,variances and covariances without calling routines
						! for efficiency (only need to traverse the block time series once).
						DO i=1,ii
							IF(hs(i).NE.ecv%nanreal.AND.&
							ws(i).NE.ecv%nanreal.AND.Tss(i).NE.ecv%nanreal) THEN
								hbar=hbar+hs(i)
								Tsbar=Tsbar+Tss(i)
								wbar=wbar+ws(i)
								covhw=covhw+hs(i)*ws(i)
								covTsw=covTsw+Tss(i)*ws(i)
								varh=varh+hs(i)*hs(i)
								varTs=varTs+Tss(i)*Tss(i)
								varw=varw+ws(i)*ws(i)
							ELSE
								nnans=nnans+1
							END IF
						END DO
						hbar=hbar/(1.0*(ii-nnans))
						Tsbar=Tsbar/(1.0*(ii-nnans))
						wbar=wbar/(1.0*(ii-nnans))
						covhw=(covhw/(1.0*(ii-nnans)))-1.0*hbar*wbar
						covTsw=(covTsw/(1.0*(ii-nnans)))-1.0*Tsbar*wbar
						varh=(varh/(1.0*(ii-nnans)))-1.0*hbar*hbar
						varTs=(varTs/(1.0*(ii-nnans)))-1.0*Tsbar*Tsbar
						varw=(varw/(1.0*(ii-nnans)))-1.0*wbar*wbar
						
						! Add these to the variance of the covariances.
						ffsetsw=(varts**2)*(varw**2)+(covtsw)**2
						ffsehw=(varh**2)*(varw**2)+(covhw)**2	
			
						! Calculate the block auto and cross covariances up to the maxlag index.
						DO k=1,lagmax	! Lag loop.
							! Set temporary reals to zero on each pass.
							atst=0
							awt=0
							aht=0
							ctswt=0
							cwtst=0
							chwt=0
							cwht=0 
							
							! Calculate auto and cross covariances at the given lag.
							DO i=1,ii-k
								IF(hs(i).NE.ecv%nanreal.AND.hs(i+k).NE.ecv%nanreal) THEN
									aht=aht+(hs(i)-hbar)*(hs(i+k)-hbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(ws(i).NE.ecv%nanreal.AND.ws(i+k).NE.ecv%nanreal) THEN
									awt=awt+(ws(i)-wbar)*(ws(i+k)-wbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(Tss(i).NE.ecv%nanreal.AND.Tss(i+k).NE.ecv%nanreal) THEN
									aTst=aTst+(Tss(i)-Tsbar)*(Tss(i+k)-wbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(Tss(i+k).NE.ecv%nanreal.AND.ws(i).NE.ecv%nanreal) THEN
									cwtst=cwtst+(Tss(i+k)-Tsbar)*(ws(i)-wbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(ws(i+k).NE.ecv%nanreal.AND.Tss(i).NE.ecv%nanreal) THEN
									ctswt=ctswt+(ws(i+k)-wbar)*(Tss(i)-Tsbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(hs(i+k).NE.ecv%nanreal.AND.ws(i).NE.ecv%nanreal) THEN
									cwht=cwht+(hs(i+k)-hbar)*(ws(i)-wbar)*(1.0/(1.0*(ii-nnans)))
								END IF
								IF(ws(i+k).NE.ecv%nanreal.AND.hs(i).NE.ecv%nanreal) THEN
									chwt=chwt+(ws(i+k)-wbar)*(hs(i)-hbar)*(1.0/(1.0*(ii-nnans)))
								END IF
							END DO

							! Add these to the variances of the covariances.
							ffsetsw=ffsetsw+2*(atst*awt)+2*(cwtst*ctswt)
							ffsehw=ffsehw+2*(aht*awt)+2*(cwht*chwt)
	
						END DO
		
						! Normalize the ffses to arrive at the variance of the covariances.
						ffsetsw=ffsetsw/(1.0*(ii-nnans))
						ffsehw=ffsehw/(1.0*(ii-nnans))
		
						! Take the square root and then normalize by the covariance to arrive at the ffses.
						ffsetsw=SQRT(ffsetsw)
						ffsehw=SQRT(ffsehw)
						
						! For (low) absolute fluxes the FFSE can become very large (i.e. >1) as both the
						! variance of the covariance and the covariance itself are small. For all practical
						! purposes we can simply set these FFSEs to 1 in that the fluxes here are of such
						! a small magnitude that they have little effect on the SEB, we still get
						! an accurate impression of the uncertainty associated with these fluxes.
						!IF(ffsetsw.GT.1) THEN
						!	ffsetsw=1
						!END IF
						!IF(ffsehw.GT.1) THEN
						!	ffsehw=1
						!END IF

						! Finally assign the block ffses:	
						ecv%ffsetsw(j)=ffsetsw
						ecv%ffsehw(j)=ffsehw
	
					ELSE
						ecv%ffsetsw(j)=ecv%nanreal
						ecv%ffsehw(j)=ecv%nanreal
					END IF
				
				END DO
			ELSE
				DO j=1,jj	! Block loop.
				
					! Reset statistics on each pass.
					ffseTsw=0
					ffsehw=0
					covTsw=0
					covhw=0
					covTsw0=0
					covhw0=0
					varTs=0
					varh=0
					varw=0
					hbar=0
					Tsbar=0
					wbar=0
					nnans=0
					tempa(:,:)=0
						
						
							
					! Index bounds for the block.
					ll=(j-1)*ii+1
					ul=j*ii	
			
					! Set the values of the block arrays.
					hs=ecv%h(ll:ul)
					Tss=ecv%Ts(ll:ul)
					ws=ecv%w(ll:ul)
						
					! Interpolate flaged entries using the gapfill function. 
					tempa(:,1)=hs(:)
					tempa(:,2)=Tss(:)
					tempa(:,3)=ws(:)
					CALL gapfill(ecv,ii,3,tempa)
					hs(:)=tempa(:,1)
					Tss(:)=tempa(:,2)
					ws(:)=tempa(:,3)

					! Compute the non-LDT covariances.
					CALL covariance(ecv,ii,hs,ws,covhw0)
					CALL covariance(ecv,ii,Tss,ws,covTsw0)

					! Important to detrend these; particularly when computing the autocovariance
					! as it is not well defined in the presence of even a small trend.
					CALL detrend(ecv,ii,hs)
					CALL detrend(ecv,ii,Tss)
					CALL detrend(ecv,ii,ws)

					! Calculate the block averages,variances and covariances without calling routines
					! for efficiency (only need to traverse the block time series once).
					DO i=1,ii
						IF(hs(i).NE.ecv%nanreal.AND.&
						ws(i).NE.ecv%nanreal.AND.Tss(i).NE.ecv%nanreal) THEN
							hbar=hbar+hs(i)
							Tsbar=Tsbar+Tss(i)
							wbar=wbar+ws(i)
							covhw=covhw+hs(i)*ws(i)
							covTsw=covTsw+Tss(i)*ws(i)
							varh=varh+hs(i)*hs(i)
							varTs=varTs+Tss(i)*Tss(i)
							varw=varw+ws(i)*ws(i)
						ELSE
							nnans=nnans+1
						END IF
					END DO
					hbar=hbar/(1.0*(ii-nnans))
					Tsbar=Tsbar/(1.0*(ii-nnans))
					wbar=wbar/(1.0*(ii-nnans))
					covhw=(covhw/(1.0*(ii-nnans)))-1.0*hbar*wbar
					covTsw=(covTsw/(1.0*(ii-nnans)))-1.0*Tsbar*wbar
					varh=(varh/(1.0*(ii-nnans)))-1.0*hbar*hbar
					varTs=(varTs/(1.0*(ii-nnans)))-1.0*Tsbar*Tsbar
					varw=(varw/(1.0*(ii-nnans)))-1.0*wbar*wbar
						
					! Add these to the variance of the covariances.
					ffsetsw=(varts**2)*(varw**2)+(covtsw)**2
					ffsehw=(varh**2)*(varw**2)+(covhw)**2	
			
					! Calculate the block auto and cross covariances up to the maxlag index.
					DO k=1,lagmax	! Lag loop.
						! Set temporary reals to zero on each pass.
						atst=0
						awt=0
						aht=0
						ctswt=0
						cwtst=0
						chwt=0
						cwht=0 
							
						! Calculate auto and cross covariances at the given lag.
						DO i=1,ii-k
							IF(hs(i).NE.ecv%nanreal.AND.hs(i+k).NE.ecv%nanreal) THEN
								aht=aht+(hs(i)-hbar)*(hs(i+k)-hbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(ws(i).NE.ecv%nanreal.AND.ws(i+k).NE.ecv%nanreal) THEN
								awt=awt+(ws(i)-wbar)*(ws(i+k)-wbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(Tss(i).NE.ecv%nanreal.AND.Tss(i+k).NE.ecv%nanreal) THEN
								aTst=aTst+(Tss(i)-Tsbar)*(Tss(i+k)-wbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(Tss(i+k).NE.ecv%nanreal.AND.ws(i).NE.ecv%nanreal) THEN
								cwtst=cwtst+(Tss(i+k)-Tsbar)*(ws(i)-wbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(ws(i+k).NE.ecv%nanreal.AND.Tss(i).NE.ecv%nanreal) THEN
								ctswt=ctswt+(ws(i+k)-wbar)*(Tss(i)-Tsbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(hs(i+k).NE.ecv%nanreal.AND.ws(i).NE.ecv%nanreal) THEN
								cwht=cwht+(hs(i+k)-hbar)*(ws(i)-wbar)*(1.0/(1.0*(ii-nnans)))
							END IF
							IF(ws(i+k).NE.ecv%nanreal.AND.hs(i).NE.ecv%nanreal) THEN
								chwt=chwt+(ws(i+k)-wbar)*(hs(i)-hbar)*(1.0/(1.0*(ii-nnans)))
							END IF
						END DO

						! Add these to the variances of the covariances.
						ffsetsw=ffsetsw+2*(atst*awt)+2*(cwtst*ctswt)
						ffsehw=ffsehw+2*(aht*awt)+2*(cwht*chwt)
	
					END DO
		
					! Normalize the ffses to arrive at the variance of the covariances.
					ffsetsw=ffsetsw/(1.0*(ii-nnans))
					ffsehw=ffsehw/(1.0*(ii-nnans))
		
					! Take the square root and then normalize by the covariance to arrive at the ffses.
					ffsetsw=SQRT(ffsetsw)
					ffsehw=SQRT(ffsehw)

					! For (low) absolute fluxes the FFSE can become very large (i.e. >1) as both the
					! variance of the covariance and the covariance itself are small. For all practical
					! purposes we can simply set these FFSEs to 1 in that the fluxes here are of such
					! a small magnitude that they have little effect on the SEB, we still get
					! an accurate impression of the uncertainty associated with these fluxes.
					

					! Finally assign the block ffses:	
					ecv%ffsetsw(j)=ffsetsw
					ecv%ffsehw(j)=ffsehw

					
				END DO
				
			END IF

			! Deallocate temporary arrays.
			DEALLOCATE(hs)
			DEALLOCATE(ws)
			DEALLOCATE(Tss)	
			DEALLOCATE(tempa)

			
			PRINT*,'--> Max ffse(tsw)=',MAXVAL(ecv%ffsetsw)		
			WRITE(ecv%loglun,*)'--> Max ffse(tsw)=',MAXVAL(ecv%ffsetsw)
			PRINT*,'--> Max ffse(hw)=',MAXVAL(ecv%ffsehw)		
			WRITE(ecv%loglun,*)'--> Max ffse(hw)=',MAXVAL(ecv%ffsehw)

			PRINT*,'eof ffse routine'	
			PRINT*,'*******************************************************************************'	
			WRITE(ecv%loglun,*)'eof ffse routine.'	
			PRINT*,'*******************************************************************************'
				
			RETURN
		END SUBROUTINE


		! Account for the propogation of uncertainty through the flux corrections for a given block following Billesbach (2011). NB. must be called after
		! the flux corrections routine and the blockffse routine!
		SUBROUTINE 	propogate(ecv,jblock)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER, INTENT(IN)		:: jblock			! Block index.
			REAL*8				:: mu=28.97/18.01528		! Molar mass ratio from WPL correction.
			REAL*8				:: Rd=287.04			! Gas constant for dry air [JK^-1kg^-1]
			REAL*8				:: sigma			! Ratio of block average water vapor to dry air density.
			REAL*8				:: unct,unch			! To hold the block variances of the covariances.
			REAL*8				:: temp1,temp2			! Temporary reals.
			INTEGER				:: j
			
			j=jblock			! Shorthand block index.
			unct=ecv%ffsetsw(j)		! Assign the variance of cov(T_s,w).
			unch=ecv%ffsehw(j)		! Assign the variance of cov(rho_v,w).
			sigma=ecv%hbar(j)/ecv%rhobar(j)	! Compute the block averaged density ratio.

		
			
			! --> These equations are detailed in appendix A.5.

			! First for the sensible heat flux: 
			temp1=((1.0/ecv%Fa_sen(j))*(1-0.51*Rd*ecv%hbar(j)*ecv%tbar(j)/ecv%pbar(j)))**2
			temp2=((ecv%CFsep(j)/ecv%Fa_lat(j))*(0.51*Rd*ecv%Tsbar(j)*ecv%Tbar(j)/ecv%pbar(j)))**2

			ecv%ffsetsw(j)=SQRT(temp1*(unct**2)+temp2*(unch**2))/ABS(ecv%Twcov(j)) ! <-- Sensible heat flux fractional uncertainty. 

			! Next the latent heat flux.

			temp1=((ecv%CFsep(j)/ecv%Fa_lat(j))*(1-0.51*Rd*ecv%tsbar(j)*ecv%hbar(j)/ecv%pbar(j)))**2
			temp2=((ecv%hbar(j)/(ecv%Fa_sen(j)*ecv%Tbar(j)))*(1-0.51*Rd*ecv%hbar(j)*ecv%Tbar(j)/ecv%pbar(j)))**2
			
			ecv%ffsehw(j)=SQRT(temp1*(unch**2)+temp2*(unct**2))/ABS(ecv%hwcov(j))  ! <-- Latent heat flux fractional uncertainty.


			! Done but reset values for next pass.
			j=0
			unct=0
			unch=0
			temp1=0
			temp2=0
			sigma=0
			
			RETURN

		END SUBROUTINE
		

		SUBROUTINE	blockflags(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER				:: i,j,k,ll,ul,n	! Counters and bounds.
			INTEGER				:: nstat,nnan,ndir,nitc,nv,n2! Hard flag counters.
			INTEGER				:: ii,jj,ss		! Block variable array dimensions.
			INTEGER				:: stuw,stvw,stqw,stTw	! Trackers for the four stationarity tests.
			INTEGER	,DIMENSION(3)		:: stis			! Array containing the trackers for the stationarity tests.
			INTEGER,ALLOCATABLE		:: flagstat(:)		! Stationarity flag.
			INTEGER,ALLOCATABLE		:: flagnan(:)		! Missing data flag.
			INTEGER,ALLOCATABLE		:: flagdir(:)		! Unfavorable mean wind directions due to upwind structures.
			INTEGER,ALLOCATABLE		:: flagvert(:)		! Flag high mean vertical wind speeds (after planar fit).
			INTEGER,ALLOCATABLE		:: flagitc(:)		! ITC flag.
			INTEGER,ALLOCATABLE		:: fsuw(:),fsvw(:),fsqw(:),fsTw(:),fiw(:)	! Stationarity and ITC flags.
			REAL*8				:: tewbar		! Temporary real to hold the block average vertical velocity.
			INTEGER				:: favper		! Corresponding frequency of the averaging period in Hz.
			REAL*8,ALLOCATABLE		:: us(:),vs(:)		! Temporary arrays holding block plane velocities allowing
										! for streamwise rotation.
			REAL*8				:: tempva,tempia	! Temporary values holding the vector and instantaneous
										! average wind speeds for the constancy ratio CR.	
					

			PRINT*,'*******************************************************************************'
			PRINT*,'starting flag routine.'		
			WRITE(ecv%loglun,*)'*******************************************************************************'
			WRITE(ecv%loglun,*)'starting flag routine.'	


			ii=ecv%blockl			! Number of points in a given block.
			ss=FLOOR(1.0*ii/6)		! Subsegment length, 6 per block. Corresponds to 5 mins for a 30 minute averaging period.		
			jj=FLOOR(1.0*ecv%nn/ii)		! Number of blocks of data.

			! Allocate space for the block flag array.
			ALLOCATE(ecv%flag(jj))	
				
			! Allocate the respective coefficients and individual flags.
			ALLOCATE(ecv%csuw(jj))
			ALLOCATE(ecv%csTw(jj))
			ALLOCATE(ecv%cshw(jj))
			ALLOCATE(ecv%cr(jj))
			ALLOCATE(ecv%flagstat(jj))
			ALLOCATE(ecv%flagnan(jj))
			ALLOCATE(ecv%flagdir(jj))
			ALLOCATE(ecv%flagvert(jj))
			ALLOCATE(ecv%flagitcw(jj))

			ALLOCATE(ecv%flagstatuw(jj))
			ALLOCATE(ecv%flagstatTsw(jj))
			ALLOCATE(ecv%flagstathw(jj))

			ALLOCATE(flagstat(jj))
			ALLOCATE(flagnan(jj))
			ALLOCATE(flagdir(jj))
			ALLOCATE(flagvert(jj))
			ALLOCATE(flagitc(jj))
			ALLOCATE(fsuw(jj))
			ALLOCATE(fsvw(jj))
			ALLOCATE(fsqw(jj))
			ALLOCATE(fsTw(jj))
			ALLOCATE(fiw(jj))
			ALLOCATE(ecv%znoqc(jj))
			ALLOCATE(ecv%obukhovlnoqc(jj))

			! Set all flags to zero initially.
			flagstat(:)=0
			flagnan(:)=0
			flagdir(:)=0
			flagvert(:)=0
			flagitc(:)=0
			fsuw(:)=0
			fsvw(:)=0
			fsTw(:)=0
			fsqw(:)=0
			fiw(:)=0

			ALLOCATE(us(ii))
			ALLOCATE(vs(ii))
			
			nstat=0
			nnan=0
			ndir=0
			nitc=0
			nv=0
			n2=0

			! Commence flagging of the respective blocks.	
			DO j=1,jj
				ll=(j-1)*ii+1
				ul=j*ii	
	
				us=ecv%u(ll:ul)
				vs=ecv%v(ll:ul)


				! Rotate u into the (streamline plane) mean velocity for the given block,
				! and v into the cross-wind direction. Of course this does NOT change
				! meanuh(j).
				CALL rotatestream(ecv,ii,us,vs)

				! Compute the Obukhov length and mean measurement height.
				!CALL obukhovlength(ecv,jj,us,vs,ecv%w(ll:ul),ecv%ts(ll:ul),ecv%obukhovlnoqc(j))
				!CALL temporalmean(ecv,jj,ecv%zmes(ll:ul),ecv%znoqc(j))

				! Compute the constancy ratio after Mahrt (1998) to verify that the wind direction
				! does not change much within a given block average. Due to the
				! rotation the mean of us is now the block averaged plane wind in the planar frame. To compute
				! CR we also need the time averaged instantaneous plane wind speed.
				ecv%cr(j)=0
				tempva=0	! <-- Temporary vector average*number of entries.
				tempia=0	! <-- Temporary instaneous average*number of entries.
				DO i=1,ii
					IF(us(i).NE.ecv%nanreal.AND.vs(i).NE.ecv%nanreal) THEN	! <-- Exclude nan entries from the analysis.
						tempva=tempva+us(i)
						tempia=tempia+SQRT(us(i)**2+vs(i)**2)
					END IF
				END DO
				IF(tempia.NE.0) THEN	! <-- Check that an instaneous average exists (not just nans).
					ecv%cr(j)=ABS(tempva)/tempia
				ELSE
					ecv%cr(j)=ecv%nanreal
				END IF
				! Reset values for next block.
				tempva=0
				tempia=0
			
					
				! By definition the constancy ratio is the ratio of the block vector averaged plane wind to the
				! block averaged instantaneous plane wind. Which is NOT equivalent to the expression above.

				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%% Stationarity flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				! Compute the four stationarity coefficients.
			
				CALL stationarity(ecv,ii,ss,us,ecv%w(ll:ul),ecv%csuw(j))

				CALL stationarity(ecv,ii,ss,ecv%Ts(ll:ul),ecv%w(ll:ul),ecv%csTw(j))

				CALL stationarity(ecv,ii,ss,ecv%h(ll:ul),ecv%w(ll:ul),ecv%cshw(j))



				! Set the stationarity flag from 0 to 2. If cs <=0.3 flag 0, elseif 0.3<cs<=1 flag 1, else flag 2.


				IF(ABS(ecv%csuw(j)).GT.0.3) THEN
					IF(ABS(ecv%csuw(j)).GT.1) THEN
						stuw=2
						fsuw(j)=2
						!PRINT*,'CSuw > 1 at j=',j
					ELSE
						stuw=1
						fsuw(j)=2
					END IF				
				ELSE
					stuw=0
					fsuw(j)=0
				END IF
				IF(ABS(ecv%cshw(j)).GT.0.3) THEN
					IF(ABS(ecv%cshw(j)).GT.1) THEN
						stqw=2
						fsqw(j)=2
						!PRINT*,'Cshw > 1 at j=',j
					ELSE
						stqw=1
						fsqw(j)=1
					END IF
				ELSE
					stqw=0
					fsqw(j)=0
				END IF
				IF(ABS(ecv%csTw(j)).GT.0.3) THEN
					IF(ABS(ecv%csTw(j)).GT.1) THEN
						stTw=2
						fsTw(j)=2
						!PRINT*,'CSTw > 1 at j=',j
					ELSE
						stTw=1
						fsTw(j)=1
					END IF
				ELSE
					stTw=0
					fsTw(j)=0
				END IF
				! Combine the stationarity flags for the given block.
				stis=(/	stuw,stqw,stTw /)
				! The combined flag is the maximum flag.
				IF(stis(1).GT.stis(2)) THEN
					IF(stis(1).GT.stis(3)) THEN
						IF(stis(2).GT.stis(3)) THEN
							flagstat(j)=stis(2)
						ELSE
							flagstat(j)=stis(3)
						END IF
					ELSE
						flagstat(j)=stis(1)
					END IF
				ELSE
					IF(stis(1).GT.stis(3)) THEN
						flagstat(j)=stis(1)
					ELSE
						IF(stis(3).GT.stis(2)) THEN
							flagstat(j)=stis(2)
						ELSE
							flagstat(j)=stis(3)
						END IF
					END IF
				END IF
				IF(flagstat(j).GT.1) THEN
					nstat=nstat+1
				END IF
				
				! Reset values for next pass.
				stis(:)=0
				stuw=0
				stqw=0
				stTw=0	

				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%% Missing data flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				! If spikes/internal flags exceed 10% for a given block flag the entire block (flag 2).
				n=0
				DO i=ll,ul					! First add up nan flags for the given block.
					IF(ecv%u(i).EQ.ecv%nanreal) THEN
						n=n+1
					END IF
				END DO	
				IF(((1.0*n)/ii).GE.0.1) THEN		! Then control the threshold.
					flagnan(j)=2
					nnan=nnan+1
				ELSE
					flagnan(j)=0
				END IF
				


				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%% Wind direction flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				DO i=1,SIZE(ecv%disdir,1)
					IF(ecv%mwd(j).GE.ecv%disdir(i,1).AND.ecv%mwd(j).LE.ecv%disdir(i,2)&
						.AND.ecv%cr(j).GT.0.75)	THEN	! Only flag if the constancy ratio exceeds
										! 0.7 such that the wind direction is well defined.
						flagdir(j)=2
						ndir=ndir+1
					END IF
				END DO
				
				! If no flag, flag as zero.
				IF(flagdir(j).NE.2) THEN
						flagdir(j)=0
				END IF


				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%% Integral Turbulence Characteristics Flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				! Vertical velocity integral turbulence characteristics flag. 
				CALL wITC(ecv,ii,ecv%zmes(ll:ul),us,vs,ecv%w(ll:ul),ecv%ts(ll:ul),fiw(j))
				IF(fiw(j).EQ.2) THEN
					nitc=nitc+1
				END IF
				
				

				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%% Vertical wind flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				CALL temporalmean(ecv,ii,ecv%w(ll:ul),tewbar)
				IF(ABS(tewbar).GT.0.1) THEN
					IF(ABS(tewbar).GT.0.15) THEN
						flagvert(j)=2
						nv=nv+1
					ELSE
						flagvert(j)=1
					END IF
				ELSE
					flagvert(j)=0
				END IF

				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%% Combine flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				

				ecv%flag(j)=flagstat(j)
				IF(flagnan(j).GT.ecv%flag(j)) THEN
					ecv%flag(j)=flagnan(j)
				END IF
				IF(flagdir(j).GT.ecv%flag(j)) THEN
					ecv%flag(j)=flagdir(j)
				END IF
				IF(fiw(j).GT.ecv%flag(j)) THEN
					ecv%flag(j)=fiw(j)
				END IF
				IF(flagvert(j).GT.ecv%flag(j)) THEN
					ecv%flag(j)=flagvert(j)
				END IF
				IF(ecv%flagsk(j).GT.ecv%flag(j)) THEN
					ecv%flag(j)=ecv%flagsk(j)
				END IF
		
				! If flag is 2, the block is not worth considering in further analysis.
				IF(ecv%flag(j).EQ.2) THEN
					n2=n2+1
					ecv%u(ll:ul)=ecv%nanreal
					ecv%v(ll:ul)=ecv%nanreal
					ecv%w(ll:ul)=ecv%nanreal
					ecv%TT(ll:ul)=ecv%nanreal
					ecv%Tv(ll:ul)=ecv%nanreal
					ecv%Ts(ll:ul)=ecv%nanreal
					ecv%rho(ll:ul)=ecv%nanreal
					ecv%h(ll:ul)=ecv%nanreal
					ecv%p(ll:ul)=ecv%nanreal
					ecv%q(ll:ul)=ecv%nanreal
				END IF

				! Set stationarity factors to nanreal if flagdir==2 and/or flagnan==2 since the
				! factors are poorly defined in such instances.
				IF(flagdir(j).EQ.2.OR.flagnan(j).EQ.2.OR.flagvert(j).EQ.2) THEN
					ecv%csuw(j)=ecv%nanreal
					ecv%csTw(j)=ecv%nanreal
					ecv%cshw(j)=ecv%nanreal
				END IF 

			END DO

			ecv%flagstatuw(:)=fsuw(:)
			ecv%flagstatTsw(:)=fsTw(:)
			ecv%flagstathw(:)=fsqw(:)

			ecv%flagstat(:)=flagstat(:)
			ecv%flagnan(:)=flagnan(:)
			ecv%flagdir(:)=flagdir(:)
			ecv%flagvert(:)=flagvert(:)
			ecv%flagitcw(:)=fiw(:)

			! Need to allocate individual flags to ecv structure and write to results!

			! Print a summary to the terminal.
			PRINT "('-->Set flags for',I6,' blocks')",jj
			PRINT "('-->Hard flagged (flag 2)',f8.2,' (%) of blocks for exceeding w-itc limits')",(100.0*nitc)/jj
			PRINT "('-->Hard flagged (flag 2)',f8.2,' (%) of blocks for non-stationarity')",(100.0*nstat)/jj
			PRINT "('-->Hard flagged (flag 2)',f8.2,' (%) of blocks for spike limit')",(100.0*nnan)/jj
			PRINT "('-->Hard flagged (flag 2)',f8.2,' (%) of blocks for possible flow distortion')",(100.0*ndir)/jj
			PRINT "('-->Hard flagged (flag 2)',f8.2,' (%) of blocks for excess block vertical wind speed')",(100.0*nv)/jj
			PRINT "('-->Removed',f8.2,' (%) of blocks due to hard flagging')", (100.0*n2)/jj

			WRITE(ecv%loglun,*)'-->Set flags for',jj,' blocks'
			WRITE(ecv%loglun,*)'-->Hard flagged (flag 2)',(100.0*nitc)/jj,' (%) of blocks for exceeding w-itc limits'
			WRITE(ecv%loglun,*)'-->Hard flagged (flag 2)',(100.0*nstat)/jj,' (%) of blocks for non-stationarity'
			WRITE(ecv%loglun,*)'-->Hard flagged (flag 2)',(100.0*nnan)/jj,' (%) of blocks for spike limit'
			WRITE(ecv%loglun,*)'-->Hard flagged (flag 2)',(100.0*ndir)/jj,' (%) of blocks for possible flow distortion'
			WRITE(ecv%loglun,*)'-->Hard flagged (flag 2)',(100.0*nv)/jj,' (%) of blocks for excess block vertical wind speed'
			WRITE(ecv%loglun,*)'-->Removed',(100.0*n2)/jj,' (%) of blocks due to hard flagging'

			PRINT*,'eof flag routine.'		
			PRINT*,'*******************************************************************************'

			WRITE(ecv%loglun,*)'eof flag routine.'		
			WRITE(ecv%loglun,*)'*******************************************************************************'

			RETURN
		END SUBROUTINE

				

				
			


		SUBROUTINE 	blockvalues(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			INTEGER				:: i,j,k,ll,ul	! Counters and bounds.
			INTEGER				:: lls,uls	! Bounds for sub blocks.
			INTEGER				:: llss,ulss	! Bounds for subsub blocks.
			INTEGER				:: ii,jj	! Block variable array dimensions.
			INTEGER				:: iis,iiss	! Sub (and subsub) block dimensions.
			INTEGER				:: ns,nss	! Number of subblocks and subsubblocks.
			REAL*8,PARAMETER		:: kappa=0.4	! Von Karman constant [-].
			REAL*8,PARAMETER		:: g=9.81	! Acceleration of gravity [ms^-2].
			REAL*8,PARAMETER		:: cpd=1004.67	! Specific heat at constant pressure for dry air [JK^-1kg^-1].
			REAL*8,PARAMETER		:: ld0=2501000	! Specific heat of evaparotation @ T=0 Celsius [Jkg^-1] based on TK3.
			REAL*8				:: cp		! Specific heat at constant pressure for moist air.	
			REAL*8				:: ld		! Specific heat of evaporation.
			REAL*8,ALLOCATABLE		:: us(:),vs(:)	! Block velocities, allowing for rotation.
			REAL*8				:: te1,te2,te3,Etemp 	! Temporary reals.
			INTEGER				:: maxits,meanits! Counters to keep track of flux correction iterations.


			PRINT*,'*******************************************************************************'
			PRINT*,'Starting block calculations'
			WRITE(ecv%loglun,*)'*******************************************************************************'
			WRITE(ecv%loglun,*)'Starting block calculations'
			
			
			ns=3	! Set 3 subblocks, corresponds to 10 minute sub block averages with a 30 minute total block average.
			nss=6	! Set 6 subsubblock, corresponds to 5 minute subsub block averages with a 30 minute total block average.
		
			ii=ecv%blockl			! Number of points in a block.
			iis=NINT((1.0*ii)/ns)		! Nuber of points in a 10 minute sub block.
			iiss=NINT((1.0*ii)/nss)		! Number of points in a 5 minute subsub block.
			jj=FLOOR(1.0*ecv%nn/ii)		! Number of blocks of data.

			! Allocate space for the block streamwise and cross-stream velocities.
			ALLOCATE(us(ii))
			ALLOCATE(vs(ii))

					
			! Allocate space for the block averaged values.
			ALLOCATE(ecv%ubar(jj))
			ALLOCATE(ecv%vbar(jj))
			ALLOCATE(ecv%wbar(jj))
			ALLOCATE(ecv%hbar(jj))
			ALLOCATE(ecv%qbar(jj))
			ALLOCATE(ecv%Tbar(jj))
			ALLOCATE(ecv%Tsbar(jj))
			ALLOCATE(ecv%Tvbar(jj))
			ALLOCATE(ecv%pbar(jj))
			ALLOCATE(ecv%rhobar(jj))
			ALLOCATE(ecv%zbar(jj))

			! Allocate space for block covariances.
			ALLOCATE(ecv%uwcov(jj))
			ALLOCATE(ecv%vwcov(jj))
			ALLOCATE(ecv%Twcov(jj))
			ALLOCATE(ecv%Tswcov(jj))
			ALLOCATE(ecv%hwcov(jj))
			!ALLOCATE(ecv%qwcov(jj))


			! Allocate space for block variances.
			ALLOCATE(ecv%uvar(jj))
			ALLOCATE(ecv%vvar(jj))
			ALLOCATE(ecv%wvar(jj))
			ALLOCATE(ecv%Tvar(jj))
			ALLOCATE(ecv%Tsvar(jj))
			ALLOCATE(ecv%hvar(jj))
			ALLOCATE(ecv%qvar(jj))

			! Allocate space for block tke.
			ALLOCATE(ecv%tke(jj))

			! Allocate space for the block friction velocity 
			ALLOCATE(ecv%ustar(jj))
			
			! Allocate space for the block momentum flux.
			ALLOCATE(ecv%tau(jj))
		
			! Allocate space for the block Obukhov length.
			ALLOCATE(ecv%obukhovl(jj))

			! Allocate space for the corrected fluxes of latent and sensible heat.
			ALLOCATE(ecv%QE(jj))
			ALLOCATE(ecv%QH(jj))

			! Allocate space for the flux attenuation factors:
			ALLOCATE(ecv%Fa_mom(jj))
			ALLOCATE(ecv%Fa_lat(jj))
			ALLOCATE(ecv%Fa_sen(jj))

			! And the correction factors:
			ALLOCATE(ecv%CFsep(jj))
			ALLOCATE(ecv%CFSND(jj))
			ALLOCATE(ecv%CFSND0(jj))
			ALLOCATE(ecv%CFWPL(jj))

			! And the lag array:
			ALLOCATE(ecv%lagIRGA(jj))
			ALLOCATE(ecv%niter(jj))

			! Allocate the block time stamps.
			ALLOCATE(ecv%tblock(2,jj,6))

			! Set the array bounds to zero initially.
			ll=0
			ul=0
			DO j=1,jj
				! Array bounds for the current block.
				ll=(j-1)*ii+1
				ul=j*ii

				us=ecv%u(ll:ul)
				vs=ecv%v(ll:ul)
	
				
				! Check that the block isn't hard flagged.
				IF(ABS(ecv%flag(j)).LT.2) THEN
					! Set block time stamps.
					ecv%tblock(1,j,:)=ecv%t(ll,:)
					ecv%tblock(2,j,:)=ecv%t(ul,:)

					! Rotate u into the (streamline plane) mean velocity for the given block,
					! and v into the cross-wind direction. Of course this does NOT change
					! meanuh(j).
					CALL rotatestream(ecv,ii,us,vs)
		
					! Compute block averaged values.
					CALL temporalmean(ecv,ii,us,ecv%ubar(j))
					CALL temporalmean(ecv,ii,vs,ecv%vbar(j))
					CALL temporalmean(ecv,ii,ecv%w(ll:ul),ecv%wbar(j))
					CALL temporalmean(ecv,ii,ecv%Ts(ll:ul),ecv%Tsbar(j))
					CALL temporalmean(ecv,ii,ecv%Tv(ll:ul),ecv%Tvbar(j))
					CALL temporalmean(ecv,ii,ecv%TT(ll:ul),ecv%Tbar(j))
					CALL temporalmean(ecv,ii,ecv%h(ll:ul),ecv%hbar(j))
					CALL temporalmean(ecv,ii,ecv%q(ll:ul),ecv%qbar(j))
					CALL temporalmean(ecv,ii,ecv%p(ll:ul),ecv%pbar(j))
					CALL temporalmean(ecv,ii,ecv%rho(ll:ul),ecv%rhobar(j))
					CALL temporalmean(ecv,ii,ecv%zmes(ll:ul),ecv%zbar(j))

					! Calculate block variances
					CALL covariance(ecv,ii,us,us,ecv%uvar(j))
					CALL covariance(ecv,ii,vs,vs,ecv%vvar(j))
					CALL covariance(ecv,ii,ecv%w(ll:ul),ecv%w(ll:ul),ecv%wvar(j))
					CALL covariance(ecv,ii,ecv%TT(ll:ul),ecv%TT(ll:ul),ecv%Tvar(j))
					CALL covariance(ecv,ii,ecv%Ts(ll:ul),ecv%Ts(ll:ul),ecv%Tsvar(j))
					CALL covariance(ecv,ii,ecv%h(ll:ul),ecv%h(ll:ul),ecv%hvar(j))
					CALL covariance(ecv,ii,ecv%q(ll:ul),ecv%q(ll:ul),ecv%qvar(j))
				
					! Calculate the block tke.
					IF(ecv%uvar(j).NE.ecv%nanreal) THEN
						ecv%tke(j)=0.5*(ecv%uvar(j)+ecv%vvar(j)+ecv%wvar(j))
					ELSE
						ecv%tke(j)=ecv%nanreal
					END IF

					! Apply the flux corrections.
					CALL fluxcorrections(ecv,ii,us,vs,ecv%w(ll:ul),ecv%Ts(ll:ul),ecv%h(ll:ul),ecv%q(ll:ul),&
						ecv%rhobar(j),ecv%Tbar(j),ecv%pbar(j),ecv%zbar(j),ecv%mwd(j),&
						ecv%uwcov(j),ecv%vwcov(j),ecv%Tswcov(j),ecv%Twcov(j),&
						ecv%hwcov(j),Etemp,ecv%lagIRGA(j),ecv%Fa_mom(j),ecv%Fa_sen(j),ecv%Fa_lat(j),&
						ecv%ustar(j),ecv%obukhovl(j),ecv%tau(j),ecv%niter(j),ecv%CFsep(j),ecv%CFSND(j),&
						ecv%CFSND0(j),ecv%CFWPL(j))

					! Calculate the final flux uncertanties accounting for error propagation.
					CALL propogate(ecv,j) 
	
					! Convert kinematic heat fluxes to dynamic heat fluxes.
					ecv%QH(j)=ecv%Twcov(j)
					ecv%QE(j)=ecv%hwcov(j)
					IF(ecv%qbar(j).NE.ecv%nanreal.AND.ecv%rhobar(j).NE.ecv%nanreal&
						.AND.ecv%QH(j).NE.ecv%nanreal) THEN
						cp=cpd*(1+0.84*ecv%qbar(j))
						ecv%QH(j)=ecv%rhobar(j)*cp*ecv%QH(j)
					ELSE
						ecv%QH(j)=ecv%nanreal
					END IF
					IF(ecv%Tbar(j).NE.ecv%nanreal.AND.Etemp.NE.ecv%nanreal) THEN
						ld=ld0-2360*(ecv%Tbar(j)-273.15)
						ecv%QE(j)=ld*Etemp
					ELSE
						ecv%QE(j)=ecv%nanreal
					END IF


					! Reset specific heats and temp E for next pass.
					cp=0
					ld=0
					Etemp=0

					! Calculate the obukhov length.
					CALL obukhovlength(ecv,ii,us,vs,ecv%w(ll:ul),ecv%Ts(ll:ul),ecv%obukhovl(j))


					! Calculate the block tke.
					IF(ecv%uvar(j).NE.ecv%nanreal) THEN
						ecv%tke(j)=0.5*(ecv%uvar(j)+ecv%vvar(j)+ecv%wvar(j))
					ELSE
						ecv%tke(j)=ecv%nanreal
					END IF

				ELSE
					! Set block time stamps.
					ecv%tblock(1,j,:)=ecv%t(ll,:)
					ecv%tblock(2,j,:)=ecv%t(ul,:)
		
					! Set the remaining blocks values to nan.
					ecv%ubar(j)=ecv%nanreal
					ecv%vbar(j)=ecv%nanreal
					ecv%vbar(j)=ecv%nanreal
					ecv%wbar(j)=ecv%nanreal
					ecv%Tsbar(j)=ecv%nanreal
					ecv%Tvbar(j)=ecv%nanreal
					ecv%Tbar(j)=ecv%nanreal
					ecv%hbar(j)=ecv%nanreal
					ecv%qbar(j)=ecv%nanreal
					ecv%pbar(j)=ecv%nanreal
					ecv%rhobar(j)=ecv%nanreal
					ecv%uwcov(j)=ecv%nanreal
					ecv%vwcov(j)=ecv%nanreal
					ecv%Twcov(j)=ecv%nanreal
					ecv%Tswcov(j)=ecv%nanreal
					ecv%hwcov(j)=ecv%nanreal
					ecv%uvar(j)=ecv%nanreal
					ecv%vvar(j)=ecv%nanreal
					ecv%wvar(j)=ecv%nanreal
					ecv%Tvar(j)=ecv%nanreal
					ecv%Tsvar(j)=ecv%nanreal
					ecv%hvar(j)=ecv%nanreal
					ecv%qvar(j)=ecv%nanreal
					ecv%tke(j)=ecv%nanreal
					ecv%ustar(j)=ecv%nanreal
					ecv%tau(j)=ecv%nanreal
					ecv%obukhovl(j)=ecv%nanreal
					ecv%zbar(j)=ecv%nanreal
					ecv%QE(j)=ecv%nanreal
					ecv%QH(j)=ecv%nanreal
					ecv%Fa_mom(j)=ecv%nanreal
					ecv%Fa_lat(j)=ecv%nanreal
					ecv%Fa_sen(j)=ecv%nanreal
					ecv%lagIRGA(j)=NINT(ecv%nanreal)
					ecv%niter(j)=ecv%nanreal
					ecv%CFsep(j)=ecv%nanreal
					ecv%CFSND(j)=ecv%nanreal
					ecv%CFSND0(j)=ecv%nanreal
					ecv%CFWPL(j)=ecv%nanreal
					ecv%ffsetsw(j)=ecv%nanreal
					ecv%ffsehw(j)=ecv%nanreal
				END IF
			END DO
			
			CALL temporalmean(ecv,jj,ecv%niter,te1)
			meanits=NINT(te1)
			te2=MAXVAL(ecv%niter)
			maxits=NINT(te2)
			PRINT "('-->Calculated block values for',I6,' blocks')",jj
			PRINT "('-->Block lengths correspond to an averaging period of',I5,' minutes.')",ecv%avper
			PRINT "('-->Flux corrections converged in an average of',I3,' iterations.')",meanits
			PRINT "('-->And maximum',I3,' iterations.')",maxits
			PRINT*,'eof block calculations'
			PRINT*,'*******************************************************************************'
			
			WRITE(ecv%loglun,*)'-->Calculated block values for',jj,' blocks'
			WRITE(ecv%loglun,*)'-->Block lengths correspond to an averaging period of',ecv%avper,' minutes.'
			WRITE(ecv%loglun,*)'-->Flux corrections converged in an average of',meanits,' iterations.'
			WRITE(ecv%loglun,*)'-->with a maximum of',maxits,' iterations.'
			WRITE(ecv%loglun,*)'eof block calculations'
			WRITE(ecv%loglun,*)'*******************************************************************************'

			RETURN
		END SUBROUTINE
				

		! Subroutine which deallocates an array when called IF the array has previously been allocated
 ! 		SUBROUTINE deallocaterealarray(a)

!			IMPLICIT NONE
			

		! Final subroutine for deallocation in the case that multiple files are passed to the module.
		SUBROUTINE deallocateall(ecv)
			IMPLICIT NONE
			TYPE(eddycov),INTENT(INOUT)	:: ecv
			
			! Measured and auxillary variables.
			IF(ALLOCATED(ecv%u).EQV..TRUE.) THEN
				DEALLOCATE(ecv%u)
			END IF
			IF(ALLOCATED(ecv%v).EQV..TRUE.) THEN
				DEALLOCATE(ecv%v)
			END IF
			IF(ALLOCATED(ecv%w).EQV..TRUE.) THEN
				DEALLOCATE(ecv%w)
			END IF
			IF(ALLOCATED(ecv%TT).EQV..TRUE.) THEN
				DEALLOCATE(ecv%TT)
			END IF
			IF(ALLOCATED(ecv%Ts).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Ts)
			END IF
			IF(ALLOCATED(ecv%Tv).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tv)
			END IF
			IF(ALLOCATED(ecv%q).EQV..TRUE.) THEN
				DEALLOCATE(ecv%q)
			END IF
			IF(ALLOCATED(ecv%h).EQV..TRUE.) THEN
				DEALLOCATE(ecv%h)
			END IF
			IF(ALLOCATED(ecv%p).EQV..TRUE.) THEN
				DEALLOCATE(ecv%p)
			END IF
			IF(ALLOCATED(ecv%rho).EQV..TRUE.) THEN
				DEALLOCATE(ecv%rho)
			END IF
			IF(ALLOCATED(ecv%t).EQV..TRUE.) THEN
				DEALLOCATE(ecv%t)
			END IF
			IF(ALLOCATED(ecv%unode).EQV..TRUE.) THEN
				DEALLOCATE(ecv%unode)
			END IF
			IF(ALLOCATED(ecv%vnode).EQV..TRUE.) THEN
				DEALLOCATE(ecv%vnode)
			END IF
			IF(ALLOCATED(ecv%wnode).EQV..TRUE.) THEN
				DEALLOCATE(ecv%wnode)
			END IF
			IF(ALLOCATED(ecv%Tsnode).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tsnode)
			END IF
			IF(ALLOCATED(ecv%hnode).EQV..TRUE.) THEN
				DEALLOCATE(ecv%hnode)
			END IF

			! Block varaibles.
			IF(ALLOCATED(ecv%flag).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flag)
			END IF
			IF(ALLOCATED(ecv%flagstatuw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagstatuw)
			END IF
			IF(ALLOCATED(ecv%flagstatTsw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagstatTsw)
			END IF
			IF(ALLOCATED(ecv%flagstathw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagstathw)
			END IF
			IF(ALLOCATED(ecv%flagstat).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagstat)
			END IF
			IF(ALLOCATED(ecv%flagnan).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagnan)
			END IF
			IF(ALLOCATED(ecv%flagdir).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagdir)
			END IF
			IF(ALLOCATED(ecv%flagvert).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagvert)
			END IF
			IF(ALLOCATED(ecv%flagitcw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagitcw)
			END IF
			IF(ALLOCATED(ecv%mwd).EQV..TRUE.) THEN
				DEALLOCATE(ecv%mwd)
			END IF
			IF(ALLOCATED(ecv%mws).EQV..TRUE.) THEN
				DEALLOCATE(ecv%mws)
			END IF
			IF(ALLOCATED(ecv%ubar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%ubar)
			END IF
			IF(ALLOCATED(ecv%vbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%vbar)
			END IF
			IF(ALLOCATED(ecv%wbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%wbar)
			END IF
			IF(ALLOCATED(ecv%Tbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tbar)
			END IF
			IF(ALLOCATED(ecv%qbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%qbar)
			END IF
			IF(ALLOCATED(ecv%pbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%pbar)
			END IF
			IF(ALLOCATED(ecv%rhobar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%rhobar)
			END IF
			IF(ALLOCATED(ecv%hbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%hbar)
			END IF
			IF(ALLOCATED(ecv%Tvbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tvbar)
			END IF
			IF(ALLOCATED(ecv%Tsbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tsbar)
			END IF
			IF(ALLOCATED(ecv%tblock).EQV..TRUE.) THEN
				DEALLOCATE(ecv%tblock)
			END IF
			IF(ALLOCATED(ecv%uwcov).EQV..TRUE.) THEN
				DEALLOCATE(ecv%uwcov)
			END IF
			IF(ALLOCATED(ecv%vwcov).EQV..TRUE.) THEN
				DEALLOCATE(ecv%vwcov)
			END IF
			IF(ALLOCATED(ecv%Twcov).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Twcov)
			END IF
			IF(ALLOCATED(ecv%Tswcov).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tswcov)
			END IF
			IF(ALLOCATED(ecv%hwcov).EQV..TRUE.) THEN
				DEALLOCATE(ecv%hwcov)
			END IF
			IF(ALLOCATED(ecv%Tvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tvar)
			END IF
			IF(ALLOCATED(ecv%Tsvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Tsvar)
			END IF
			IF(ALLOCATED(ecv%uvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%uvar)
			END IF
			IF(ALLOCATED(ecv%vvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%vvar)
			END IF
			IF(ALLOCATED(ecv%wvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%wvar)
			END IF
			IF(ALLOCATED(ecv%hvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%hvar)
			END IF
			IF(ALLOCATED(ecv%qvar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%qvar)
			END IF
			IF(ALLOCATED(ecv%obukhovl).EQV..TRUE.) THEN
				DEALLOCATE(ecv%obukhovl)
			END IF
			IF(ALLOCATED(ecv%tke).EQV..TRUE.) THEN
				DEALLOCATE(ecv%tke)
			END IF
			IF(ALLOCATED(ecv%QE).EQV..TRUE.) THEN
				DEALLOCATE(ecv%QE)
			END IF
			IF(ALLOCATED(ecv%QH).EQV..TRUE.) THEN
				DEALLOCATE(ecv%QH)
			END IF
			IF(ALLOCATED(ecv%tau).EQV..TRUE.) THEN
				DEALLOCATE(ecv%tau)
			END IF
			IF(ALLOCATED(ecv%ustar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%ustar)
			END IF
			IF(ALLOCATED(ecv%csuw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%csuw)
			END IF
			IF(ALLOCATED(ecv%cshw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%cshw)
			END IF
			IF(ALLOCATED(ecv%csTw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%csTw)
			END IF
			IF(ALLOCATED(ecv%cr).EQV..TRUE.) THEN
				DEALLOCATE(ecv%cr)
			END IF
			IF(ALLOCATED(ecv%zbar).EQV..TRUE.) THEN
				DEALLOCATE(ecv%zbar)
			END IF
			IF(ALLOCATED(ecv%zmes).EQV..TRUE.) THEN
				DEALLOCATE(ecv%zmes)
			END IF
			IF(ALLOCATED(ecv%Fa_lat).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Fa_lat)
			END IF
			IF(ALLOCATED(ecv%Fa_sen).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Fa_sen)
			END IF
			IF(ALLOCATED(ecv%Fa_mom).EQV..TRUE.) THEN
				DEALLOCATE(ecv%Fa_mom)
			END IF
			IF(ALLOCATED(ecv%lagIRGA).EQV..TRUE.) THEN
				DEALLOCATE(ecv%lagIRGA)
			END IF
			IF(ALLOCATED(ecv%niter).EQV..TRUE.) THEN
				DEALLOCATE(ecv%niter)
			END IF
			IF(ALLOCATED(ecv%CFsep).EQV..TRUE.) THEN
				DEALLOCATE(ecv%CFsep)
			END IF
			IF(ALLOCATED(ecv%CFSND).EQV..TRUE.) THEN
				DEALLOCATE(ecv%CFSND)
			END IF
			IF(ALLOCATED(ecv%CFSND0).EQV..TRUE.) THEN
				DEALLOCATE(ecv%CFSND0)
			END IF
			IF(ALLOCATED(ecv%CFWPL).EQV..TRUE.) THEN
				DEALLOCATE(ecv%CFWPL)
			END IF
			IF(ALLOCATED(ecv%flagskewT).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagskewT)
			END IF
			IF(ALLOCATED(ecv%flagskeww).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagskeww)
			END IF
			IF(ALLOCATED(ecv%flagkurtT).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagkurtT)
			END IF
			IF(ALLOCATED(ecv%flagkurtw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagkurtw)
			END IF
			IF(ALLOCATED(ecv%flagsk).EQV..TRUE.) THEN
				DEALLOCATE(ecv%flagsk)
			END IF
			IF(ALLOCATED(ecv%ffsetsw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%ffsetsw)
			END IF
			IF(ALLOCATED(ecv%ffsehw).EQV..TRUE.) THEN
				DEALLOCATE(ecv%ffsehw)
			END IF
			IF(ALLOCATED(ecv%znoqc).EQV..TRUE.) THEN
				DEALLOCATE(ecv%znoqc)
			END IF
			IF(ALLOCATED(ecv%obukhovlnoqc).EQV..TRUE.) THEN
				DEALLOCATE(ecv%obukhovlnoqc)
			END IF

			! Spectra.
			IF(ecv%spactive.EQV..TRUE.) THEN
				DEALLOCATE(ecv%su)
				DEALLOCATE(ecv%sv)
				DEALLOCATE(ecv%sw)
				DEALLOCATE(ecv%sT)
				DEALLOCATE(ecv%sh)
				DEALLOCATE(ecv%sobukhovl)
				DEALLOCATE(ecv%sustar)			
				DEALLOCATE(ecv%sTstar)
				DEALLOCATE(ecv%sz)
				DEALLOCATE(ecv%meanuh)
				DEALLOCATE(ecv%couw)
				DEALLOCATE(ecv%covw)
				DEALLOCATE(ecv%coTw)
				DEALLOCATE(ecv%cohw)
				DEALLOCATE(ecv%oguw)
				DEALLOCATE(ecv%ogvw)
				DEALLOCATE(ecv%ogTw)
				DEALLOCATE(ecv%oghw)
				DEALLOCATE(ecv%oguu)
				DEALLOCATE(ecv%ogvv)
				DEALLOCATE(ecv%ogww)
				DEALLOCATE(ecv%ogTT)
				DEALLOCATE(ecv%oghh)
				DEALLOCATE(ecv%quuw)
				DEALLOCATE(ecv%quvw)
				DEALLOCATE(ecv%quTw)
				DEALLOCATE(ecv%quqw)
				DEALLOCATE(ecv%fsp)
				DEALLOCATE(ecv%fog)
				DEALLOCATE(ecv%tspwin)
			END IF		
			! Structure functions and autocorrelations.
			IF(ecv%dractive.EQV..TRUE.) THEN
				DEALLOCATE(ecv%d2uu)
				DEALLOCATE(ecv%d2vv)
				DEALLOCATE(ecv%d2ww)
				DEALLOCATE(ecv%d2TT)
				DEALLOCATE(ecv%d2qq)
				DEALLOCATE(ecv%d3uu)
				DEALLOCATE(ecv%d3vv)
				DEALLOCATE(ecv%d3ww)
				DEALLOCATE(ecv%d3TT)
				DEALLOCATE(ecv%d3qq)
				DEALLOCATE(ecv%d4uu)
				DEALLOCATE(ecv%d4vv)
				DEALLOCATE(ecv%d4ww)
				DEALLOCATE(ecv%d4TT)
				DEALLOCATE(ecv%d4qq)
				DEALLOCATE(ecv%ruu)
				DEALLOCATE(ecv%rvv)
				DEALLOCATE(ecv%rww)
				DEALLOCATE(ecv%rTT)
				DEALLOCATE(ecv%rqq)
				DEALLOCATE(ecv%tstwin)
				DEALLOCATE(ecv%drobukhovl)
				DEALLOCATE(ecv%drz)
				DEALLOCATE(ecv%druh)
				DEALLOCATE(ecv%tlag)
			END IF

			RETURN	
		END SUBROUTINE


		 
END MODULE structure



