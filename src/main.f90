
PROGRAM main
	USE netcdf			! Load the netcdf module.
	USE structure 			! Load the eddy covariance module.
	IMPLICIT NONE
	TYPE(eddycov)	:: ecv 		! Declare the derived data type.

	! Filename (with path).
	INTEGER, PARAMETER		:: nfiles=20
	CHARACTER*80, DIMENSION(nfiles)	:: filesare,outfilesare,logfilesare
	CHARACTER*80			:: zfilesb1,zfilesb2
	INTEGER 			:: n,ns
	REAL*8				:: start,finish,tstart,tfinish
	
	! Time the entire run.
	CALL cpu_time(start)
	tstart=start

	! Clear the warning signs in the terminal associated with 'polymorphic entities'
	CALL SYSTEM('clear')

	! Site specifications.	

	ecv%latitude=79				! Latitude in degrees north (to the nearest degree).

	ecv%CSATbearing=215			! Orientation of the sonic relative to true north (deg).

	ecv%LIbearing=245			! Orientation of the IRGA relative to true north (deg).

	ecv%dsep=0.22				! Horizontal separation distance between the IRGA and the sonic.

	ecv%pfile='../input/pNYÅ2007-2010.nc'	! Input file for slow pressure measurements.

	ecv%tfile='../input/TNYÅ1984-2013in.txt'! Input file for monthly temperature limits.

	zfilesb1='../input/diurnalzmB1.txt'	! Measurement heights, batch 1 (2007).
	zfilesb2='../input/diurnalzmB2.txt'	! Measurement heights, batch 2 (2008-2009).

	

	ecv%zbase=2.75				! Base height of the sonic AGL (m).

	ecv%fs=20				! Sampling frequency (Hz).

	
	ALLOCATE(ecv%disdir(1,2))		! <-- Allocate and set bins where considerable upwind flow distortion occurs.
	ecv%disdir(1,1)=15			! <-- Lower limit of wind bin.
	ecv%disdir(1,2)=55			! <-- Upper limit of wind bin.

	! Parameters.
	
	ecv%avper=30				! Block averaging period to be used (minutes).

	ecv%spwin=210				! Spectral window lengths to be used (minutes).
	
	ecv%drwin=30				! Autostatistical window lengths to be used (minutes).

	ecv%lag1=30				! Autostatistical time lag increment (seconds).


	! Options.
	
	ecv%smooth=.TRUE.			! Smooth and dealias spectra status (T/F). True recommended.

	ecv%detrend=.TRUE.			! Linear detrend spectra and autostatistics (T/F).

	ecv%taper=.TRUE.			! Apply a taper to the series, hence avoiding spectral
						! leakage as a result of the discontinuity between the
						! begining and the end of the series.
	
	ecv%tapertype='H'			! If the previous is true set the taper type, either a Bell taper ('B') or a Hamming window
						! ('H'), the latter is recommended.

	ecv%hfstat=.TRUE.			! Write processed high frequency measurements to file (T/F).

	ecv%hmes=.TRUE.				! Are changes in the height of the sonic measured (T/F). False -> no data available.

	! -------------------------------------------------------------------------------------------------------------
	! Define input files.

	filesare(1)='../input/TOA5_070402_070411.dat'
	outfilesare(1)='../output/TOA5_070402_070411.nc'
	logfilesare(1)='../output/TOA5_070402_070411log.txt'

	filesare(2)='../input/TOA5_070411_070501.dat'
	outfilesare(2)='../output/TOA5_070411_070501.nc'
	logfilesare(2)='../output/TOA5_070411_070501log.txt'

	filesare(3)='../input/TOA5_070519_070527_0001.dat'
	outfilesare(3)='../output/TOA5_070519_070527_0001.nc'
	logfilesare(3)='../output/TOA5_070519_070527_0001log.txt'

	filesare(4)='../input/TOA5_070519_070527_0002.dat'
	outfilesare(4)='../output/TOA5_070519_070527_0002.nc'
	logfilesare(4)='../output/TOA5_070519_070527_0002log.txt'

	filesare(5)='../input/TOA5_070519_070527_0003.dat'
	outfilesare(5)='../output/TOA5_070519_070527_0003.nc'
	logfilesare(5)='../output/TOA5_070519_070527_0003log.txt'

	filesare(6)='../input/TOA5_070601_070703_0001.dat'
	outfilesare(6)='../output/TOA5_070601_070703_0001.nc'
	logfilesare(6)='../output/TOA5_070601_070703_0001log.txt'

	filesare(7)='../input/TOA5_070601_070703_0002.dat'
	outfilesare(7)='../output/TOA5_070601_070703_0002.nc'
	logfilesare(7)='../output/TOA5_070601_070703_0002log.txt'


	filesare(8)='../input/TOA5_070601_070703_0003.dat'
	outfilesare(8)='../output/TOA5_070601_070703_0003.nc'
	logfilesare(8)='../output/TOA5_070601_070703_0003log.txt'

	filesare(9)='../input/TOA5_070601_070703_0004.dat'
	outfilesare(9)='../output/TOA5_070601_070703_0004.nc'
	logfilesare(9)='../output/TOA5_070601_070703_0004log.txt'

	filesare(10)='../input/TOA5_070601_070703_0005.dat'
	outfilesare(10)='../output/TOA5_070601_070703_0005.nc'
	logfilesare(10)='../output/TOA5_070601_070703_0005log.txt'

	filesare(11)='../input/TOA5_070601_070703_0006.dat'
	outfilesare(11)='../output/TOA5_070601_070703_0006.nc'
	logfilesare(11)='../output/TOA5_070601_070703_0006log.txt'

	filesare(12)='../input/TOA5_070601_070703_0007.dat'
	outfilesare(12)='../output/TOA5_070601_070703_0007.nc'
	logfilesare(12)='../output/TOA5_070601_070703_0007log.txt'

	filesare(13)='../input/TOA5_070601_070703_0008.dat'
	outfilesare(13)='../output/TOA5_070601_070703_0008.nc'
	logfilesare(13)='../output/TOA5_070601_070703_0008log.txt'

	filesare(14)='../input/TOA5_070601_070703_0009.dat'
	outfilesare(14)='../output/TOA5_070601_070703_0009.nc'
	logfilesare(14)='../output/TOA5_070601_070703_0009log.txt'

	filesare(15)='../input/TOA5_070703_070914_0001.dat' 

	outfilesare(15)='../output/TOA5_070703_070914_0001.nc'
	logfilesare(15)='../output/TOA5_070703_070914_0001log.txt'

	filesare(16)='../input/TOA5_070703_070914_0002.dat'
	outfilesare(16)='../output/TOA5_070703_070914_0002.nc'
	logfilesare(16)='../output/TOA5_070703_070914_0002log.txt'

	filesare(17)='../input/TOA5_070703_070914_0003.dat'
	outfilesare(17)='../output/TOA5_070703_070914_0003.nc'
	logfilesare(17)='../output/TOA5_070703_070914_0003log.txt'

	filesare(18)='../input/TOA5_070703_070914l_0001.dat'
	outfilesare(18)='../output/TOA5_070703_070914l_0001.nc'
	logfilesare(18)='../output/TOA5_070703_070914l_0001log.txt'

	filesare(19)='../input/TOA5_1549.ts_data.dat'
	outfilesare(19)='../output/TOA5_1549.ts_data.nc'
	logfilesare(19)='../output/TOA5_1549.ts_datalog.txt'

	filesare(20)='../input/TOA5_1549.ts2_data.dat'
	outfilesare(20)='../output/TOA5_1549.ts2_data.nc'
	logfilesare(20)='../output/TOA5_1549.ts2_datalog.txt'

	!---------------------------------------------------------------------------------------------------------------------------
	ns=0	! File counter.
	! Loop through input files and perform eddy covariance analysis.
	DO n=15,15 !nfiles
		ns=ns+1
		CALL cpu_time(start)
		ecv%filename=filesare(n)		! <-- Set the current input file.
		ecv%logfilename=logfilesare(n)		! <-- Set the associated log file.
		ecv%outfilename=outfilesare(n)		! <-- Set the associated output file.
		IF(n.LE.18) THEN
			ecv%zfile=zfilesb1		! <-- Set the measurement height file if batch 1.
		ELSE
			ecv%zfile=zfilesb2		! <-- Same for batch 2.
		END IF
		CALL ecv%doit()				! <-- Perform the routines.
		CALL cpu_time(finish)
		PRINT*,'************************************************************************************************************'
		PRINT*,'************************************************************************************************************'
		PRINT*,'************************************************************************************************************'
		PRINT*,'Done with data in:',filesare(n)
		PRINT*,'Processed data written to:',outfilesare(n)
		PRINT*,'See logfile:',logfilesare(n)
		PRINT '("Time elapsed=",f12.3,"(min)","for file number=",I2,"out of",I2,"total.")',(finish-start)/(60),n,nfiles
		PRINT*,'************************************************************************************************************'
		PRINT*,'************************************************************************************************************'
		PRINT*,'************************************************************************************************************'
		PRINT*,''
		PRINT*,''
		PRINT*,''
		
	END DO

	! Time the entire run.
	CALL cpu_time(tfinish)
	PRINT '("Processed",I2," files in",f12.3," (hours).")',ns,(tfinish-tstart)/(60**2)
	PRINT*,'Done.'

END PROGRAM

