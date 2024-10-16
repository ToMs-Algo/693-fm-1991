C--**--CH4072--693--C:U--7:8:2000
C--**--CH4070--693--C:U--7:8:2000
C--**--CH3657--693--C:Fix--8:6:2000
C--**--CH3656--693--C:SU--8:6:2000
C--**--CH1935--693--P:Gen--2:9:1999
C--**--CH1934--693--P:MC--2:9:1999
C--**--CH1930--693--U:D--2:9:1999
      SUBROUTINE FMSET(NPREC)
C
C  Initialize the values in common which must be set before calling
C  other FM routines.
C
C  Base and precision will be set to give at least NPREC+3 decimal
C  digits of precision (giving the user three base ten guard digits).
C
C  JBASE is set to the largest permissible power of ten.
C  JFORM1 and JFORM2 are set to 1PE format displaying NPREC significant
C  digits.
C
C  The trace option is set off, and the mode for angles in trig
C  functions is set to radians.  The rounding mode is set to symmetric
C  rounding.
C
C  KW, the unit number for all FM output, is set to 6.
C
C  The size of all arrays is controlled by defining two parameters:
C  MXNDIG is the maximum value the user can set NDIG,
C  NBITS  is the number of bits per integer word.
C
C
C             Define the array sizes:
C
C
C  Here are all the common blocks used in FM.
C
C  /FMUSER/, /FM/, and /FMSAVE/ should also be declared in the main
C  program, because some compilers allocate and free space used for
C  labelled common which is declared only in subprograms.  This causes
C  the saved information to be lost.
C
C             FMUSER contains values which may need to be changed by
C             the calling program.
C
C
C             FM contains the work array used by the low-level
C             arithmetic routines, definitions for overflow and
C             underflow thresholds, etc.
C
C
C
C             FMSAVE contains information about saved constants.
C
C
C             MJSUMS is an array which can contain several FM numbers
C             being used to accumulate the concurrent sums in FMEXP2
C             and FMSIN2.  When MXNDIG is 256 eight is about the maximum
C             number of sums needed (but this depends on JBASE).  For
C             larger MXNDIG dimensioning MJSUMS to hold more than eight
C             FM numbers could speed the functions up.
C
C
C             MBUFF is a character array used by FMPRNT for printing
C             output from FMOUT.  This array may also be used for calls
C             to FMOUT from outside the FM package.
C
C
C
C             FM1 contains scratch arrays for temporary storage of FM
C             numbers while computing various functions.
C
C
C             FMPCK contains scratch arrays used to hold input arguments
C             in unpacked format when the packed versions of functions
C             are used.
C
C
C
C             KW is the unit number for all output from the FM package.
C                This includes trace output and error messages.
C
C     .. Parameters ..
      INTEGER MXNDIG,NBITS,LUNPCK,LMWA,LJSUMS,LMBUFF
      PARAMETER (MXNDIG=256,NBITS=32,LUNPCK= (6*MXNDIG)/5+20,
     +          LMWA=2*LUNPCK,LJSUMS=8*LUNPCK,
     +          LMBUFF= ((LUNPCK+3)* (NBITS-1)*301)/2000+6)
C     ..
C     .. Scalar Arguments ..
      INTEGER NPREC
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MJSUMS(LJSUMS),MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),
     +        MLN3(LUNPCK),MLN4(LUNPCK),MPA(LUNPCK),MPB(LUNPCK),
     +        MPISAV(LUNPCK),MWA(LMWA)
      CHARACTER MBUFF(LMBUFF)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER K,NPSAVE
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH
      REAL R1MACH
      INTEGER I1MACH
      EXTERNAL D1MACH,R1MACH,I1MACH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,LOG10,MAX,MIN,REAL,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMBUFF/MBUFF
      COMMON /FMPCK/MPA,MPB
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMSUMS/MJSUMS
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KW = I1MACH(2)
C
C             MAXINT should be set to a very large integer, possibly
C                    the largest representable integer for the current
C                    machine.  For most 32-bit machines MAXINT is set to
C                    2**31 - 1 = 2 147 483 647.
C
C                    Setting MAXINT to a smaller number is ok, but this
C                    unnecessarily restricts the permissible range of
C                    JBASE and MXEXP.  Too small a value of MAXINT will
C                    also slow the elementary functions like SIN, EXP,
C                    etc., since MXBASE = SQRT(MAXINT) is used to
C                    determine how many terms can be combined when
C                    summing series.
C
C                    The following code should set MAXINT to the largest
C                    representable number of the form 2**J - 1.
C
C          WARNING:  This loop causes integer overflow to occur, so it
C                    is a likely place for the program to fail when
C                    run on a different machine.  The loop below has
C                    been used successfully on a Fortran-77 compiler
C                    on each of these machines:
C                    IBM 3090, CDC 176, CRAY XMP, MAGNUSON, MIPS 800,
C                    SUN 4/260, IBM PC, MACINTOSH, MACINTOSH II,
C                    COMPAQ 386/20, TRS-80/16.
C                    However, even different versions of the same
C                    compiler may react differently, so check the value
C                    of MAXINT if there are problems installing FM on
C                    a new machine.
C
      MAXINT = I1MACH(9)
C
C             DPMAX should be set to a value near the machine's double
C                   precision overflow threshold, so that DPMAX and
C                   1.0D0/DPMAX are both representable in double
C                   precision.
C
      DPMAX = D1MACH(2)
C
C             SPMAX should be set to a value near the machine's single
C                   precision overflow threshold, so that 1.01*SPMAX
C                   and 1.0/SPMAX are both representable in single
C                   precision.
C
      SPMAX = R1MACH(2)
C
C             MXNDG2 is the maximum value for NDIG which can be used
C                    internally.  Several of the FM routines may raise
C                    NDIG above MXNDIG temporarily, in order to
C                    compute correctly rounded results.
C                    In the definition of LUNPCK The '6/5' condition
C                    allows for converting from a large base to the
C                    (smaller) largest power of ten base for output
C                    conversion.
C                    The '+ 20' condition allows for the need to carry
C                    many guard digits when using a small base like 2.
C
      MXNDG2 = LUNPCK - 1
C
C             MXBASE is the maximum value for JBASE.
C
      TEMP = MAXINT
      MXBASE = SQRT(TEMP)
C
C             JBASE is the currently used base for arithmetic.
C
      K = LOG10(REAL(MXBASE))
      JBASE = 10**K
C
C             NDIG is the number of digits currently being carried.
C
      NPSAVE = NPREC
      NDIG = 2 + (NPREC+2)/K
      IF (NDIG.LT.2 .OR. NDIG.GT.MXNDIG) THEN
          NDIG = MAX(2,MIN(MXNDIG,NDIG))
          WRITE (KW,FMT=9000) NPREC,NDIG
          NPSAVE = 0
      END IF
C
C             KFLAG is the flag for error conditions.
C
      KFLAG = 0
C
C             NTRACE is the trace switch.  Default is no printing.
C
      NTRACE = 0
C
C             LVLTRC is the trace level.  Default is to trace only
C                    routines called directly by the user.
C
      LVLTRC = 1
C
C             NCALL is the call stack pointer.
C
      NCALL = 0
C
C             Some constants which are often needed are stored in the
C             maximum precision which they have been computed with the
C             currently used base.  This speeds up the trig, log, power,
C             and exponential functions.
C
C             NDIGPI is the number of digits available in the currently
C                    stored value of pi (MPISAV).
C
      NDIGPI = 0
C
C             NJBPI is the value of JBASE for the currently stored
C                   value of pi.
C
      NJBPI = 0
C
C             NDIGE is the number of digits available in the currently
C                   stored value of e (MESAV).
C
      NDIGE = 0
C
C             NJBE is the value of JBASE for the currently stored
C                  value of e.
C
      NJBE = 0
C
C             NDIGLB is the number of digits available in the currently
C                    stored value of LN(JBASE) (MLBSAV).
C
      NDIGLB = 0
C
C             NJBLB is the value of JBASE for the currently stored
C                   value of LN(JBASE).
C
      NJBLB = 0
C
C             NDIGLI is the number of digits available in the currently
C                    stored values of the four logarithms used by FMLNI
C                    MLN1 - MLN4.
C
      NDIGLI = 0
C
C             NJBLI is the value of JBASE for the currently stored
C                   values of MLN1 - MLN4.
C
      NJBLI = 0
C
C             MXEXP  is the current maximum exponent.
C             MXEXP2 is the internal maximum exponent. This is used to
C                    define the overflow and underflow thresholds.
C
C             These values are chosen so that FM routines can raise the
C             overflow/underflow limit temporarily while computing
C             intermediate results, and so that EXP(MAXINT) is greater
C             than MXBASE**(MXEXP2+1).
C
C             The overflow threshold is JBASE**(MXEXP+1), and the
C             underflow threshold is JBASE**(-MXEXP-1).
C             This means the valid exponents in the first word of an FM
C             number can range from -MXEXP to MXEXP+1 (inclusive).
C
      TEMP = MXBASE
      MXEXP = (TEMP*TEMP)/ (2.0*LOG(TEMP)) - 1.0
      MXEXP2 = 2*MXEXP + MXEXP/100
C
C             KARGSW is a switch used by some of the elementary function
C                    routines to disable argument checking while doing
C                    calculations where no exceptions can occur.
C                    See FMARGS for a description of argument checking.
C                    KARGSW = 0 is the normal setting,
C                           = 1 means argument checking is disabled.
C
      KARGSW = 0
C
C             KEXPUN is the exponent used as a special symbol for
C                    underflowed results.
C
      KEXPUN = -MXEXP2 - 5*MXNDIG
C
C             KEXPOV is the exponent used as a special symbol for
C                    overflowed results.
C
      KEXPOV = -KEXPUN
C
C             KUNKNO is the exponent used as a special symbol for
C                    unknown FM results (1/0, SQRT(-3.0), etc).
C
      KUNKNO = KEXPOV + 5*MXNDIG
C
C             RUNKNO is returned from FM to real or double conversion
C                    routines when no valid result can be expressed in
C                    real or double precision.  On systems which provide
C                    a value for undefined results (e.g., Not A Number)
C                    setting RUNKNO to that value is reasonable.  On
C                    other systems set it to a value which is likely to
C                    make any subsequent results obviously wrong which
C                    use it.  In either case a KFLAG = -4 condition is
C                    also returned.
C
      RUNKNO = -1.01*SPMAX
C
C             IUNKNO is returned from FM to integer conversion routines
C                    when no valid result can be expressed as a one word
C                    integer.  KFLAG = -4 is also set.
C
      IUNKNO = -MXBASE*MXBASE
C
C             JFORM1 indicates the format used by FMOUT.
C
      JFORM1 = 1
C
C             JFORM2 indicates the number of digits used in FMOUT.
C
      JFORM2 = NPSAVE
C
C             KRAD = 1 indicates that trig functions use radians,
C                  = 0 means use degrees.
C
      KRAD = 1
C
C             KWARN = 0 indicates that no warning message is printed
C                       and execution continues when UNKNOWN or another
C                       exception is produced.
C                   = 1 means print a warning message and continue.
C                   = 2 means print a warning message and stop.
C
      KWARN = 1
C
C             KROUND = 1   Causes all results to be rounded to the
C                          nearest FM number, or to the value with
C                          an even last digit if the result is halfway
C                          between two FM numbers.
C                    = 0   Causes all results to be chopped.
C
      KROUND = 1
C
      RETURN

 9000 FORMAT (//' PRECISION OUT OF RANGE WHEN CALLING FMSET.','  NPREC',
     +       ' =',I20,/' THE NEAREST VALID NDIG WILL BE USED',' INSTEA',
     +       'D:   NDIG =',I6,//)
      END
      SUBROUTINE FMABS(MA,MB)
C
C  MB = ABS(MA)
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KWSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMNTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 1
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
C
      KFLAG = 0
      KWSV = KWARN
      KWARN = 0
      CALL FMEQU(MA,MB,NDIG,NDIG)
      MB(2) = ABS(MB(2))
      KWARN = KWSV
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMACOS(MA,MB)
C
C  MB = ACOS(MA)
C
C
C
C
C
C
C             Scratch array usage during FMACOS:   M01 - M06
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,MA2,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMATAN,FMDIV,FMENTR,FMEQU,FMEXIT,FMI2M,FMMPY,FMPI,
     +         FMSQRT,FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(2,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      MA2 = MA(2)
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Use ACOS(X) = ATAN(SQRT(1-X*X)/X)
C
      MB(2) = ABS(MB(2))
      CALL FMI2M(1,M05)
      CALL FMSUB(M05,MB,M03)
      CALL FMADD(M05,MB,M04)
      CALL FMMPY(M03,M04,M04)
      CALL FMSQRT(M04,M04)
      CALL FMDIV(M04,MB,MB)
C
      CALL FMATAN(MB,MB)
C
      IF (MA2.LT.0) THEN
          IF (KRAD.EQ.1) THEN
              CALL FMPI(M05)

          ELSE
              CALL FMI2M(180,M05)
          END IF

          CALL FMSUB(M05,MB,MB)
      END IF
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMADD(MA,MB,MC)
C
C  MC = MA + MB
C
C  This routine performs the trace printing for addition.
C  FMADD2 is used to do the arithmetic.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD2,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 3
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      CALL FMADD2(MA,MB,MC)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMADD2(MA,MB,MC)
C
C  Internal addition routine.  MC = MA + MB
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER J,JCOMP,JSIGN,K,KBIGMA,KRESLT,MA2,MB2,N1,N2,NGUARD,NLAST,
     +        NMWA
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADDN,FMADDP,FMARGS,FMEQU,FMMOVE,FMRND,FMRSLT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (KARGSW.NE.1) THEN
          CALL FMARGS(3,2,MA,MB,KRESLT)
          IF (KRESLT.NE.0) THEN
              CALL FMRSLT(MA,MB,MC,KRESLT)
              RETURN

          END IF

      ELSE
          IF (MA(2).EQ.0) THEN
              CALL FMEQU(MB,MC,NDIG,NDIG)
              KFLAG = 1
              RETURN

          END IF

          IF (MB(2).EQ.0) THEN
              CALL FMEQU(MA,MC,NDIG,NDIG)
              KFLAG = 1
              RETURN

          END IF

      END IF
C
      KFLAG = 0
      N1 = NDIG + 1
C
C             NGUARD is the number of guard digits used.
C
      IF (NCALL.GT.1) THEN
          NGUARD = 1

      ELSE
          B = JBASE
          NGUARD = 5.0*LOG(10.0)/LOG(B) + 2.0
          IF (NGUARD.GT.NDIG) NGUARD = NDIG
      END IF

      NMWA = N1 + NGUARD
C
C             Save the signs of MA and MB and then work with
C             positive numbers.
C             JSIGN is the sign of the result of MA + MB.
C
      JSIGN = 1
      MA2 = MA(2)
      MB2 = MB(2)
      MA(2) = ABS(MA(2))
      MB(2) = ABS(MB(2))
C
C             See which one is larger in absolute value.
C
      IF (MA(1).GT.MB(1)) THEN
          JCOMP = 1
          GO TO 20

      END IF

      IF (MB(1).GT.MA(1)) THEN
          JCOMP = 3
          GO TO 20

      END IF

      NLAST = NDIG + 1
C
      DO 10 J = 2,NLAST
          IF (ABS(MA(J)).GT.ABS(MB(J))) THEN
              JCOMP = 1
              GO TO 20

          END IF

          IF (ABS(MB(J)).GT.ABS(MA(J))) THEN
              JCOMP = 3
              GO TO 20

          END IF

   10 CONTINUE
C
      JCOMP = 2
C
   20 IF (JCOMP.LT.3) THEN
          KBIGMA = 1
          IF (MA2.LT.0) JSIGN = -1
          IF (MA2*MB2.GT.0) THEN
              CALL FMADDP(MA,MB,NMWA)

          ELSE
              CALL FMADDN(MA,MB,NGUARD,NMWA)
          END IF

      ELSE
          KBIGMA = 0
          IF (MB2.LT.0) JSIGN = -1
          IF (MA2*MB2.GT.0) THEN
              CALL FMADDP(MB,MA,NMWA)

          ELSE
              CALL FMADDN(MB,MA,NGUARD,NMWA)
          END IF

      END IF

      MA(2) = MA2
      MB(2) = MB2
C
C
C             Round the result.
C
      CALL FMRND(NDIG,NGUARD,0)
C
C             See if the result is equal to one of the input arguments.
C
      K = ABS(MA(1)-MB(1))
      IF (K.LT.NDIG) GO TO 50
      IF (K.GT.NDIG+1) THEN
          KFLAG = 1
          GO TO 50

      END IF
C
      N2 = NDIG + 4
      IF (KBIGMA.EQ.1) THEN
          DO 30 J = 3,N1
              IF (MWA(N2-J).NE.MA(N2-J)) GO TO 50
   30     CONTINUE
          IF (MWA(1).NE.MA(1)) GO TO 50
          IF (MWA(2).NE.ABS(MA(2))) GO TO 50

      ELSE
          DO 40 J = 3,N1
              IF (MWA(N2-J).NE.MB(N2-J)) GO TO 50
   40     CONTINUE
          IF (MWA(1).NE.MB(1)) GO TO 50
          IF (MWA(2).NE.ABS(MB(2))) GO TO 50
      END IF

      KFLAG = 1
C
C             Transfer to MC and fix the sign of the result.
C
   50 CALL FMMOVE(MC)
      IF (JSIGN.LT.0) MC(2) = -MC(2)
C
      RETURN

      END
      SUBROUTINE FMADDN(MA,MB,NGUARD,NMWA)
C
C  Internal addition routine.  MWA = MA - MB
C  The arguments are such that MA.GE.MB.GE.0.
C
C  NGUARD is the number of guard digits being carried.
C  NMWA is the number of words in MWA which will be used.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER NGUARD,NMWA
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,K,K2,KB,KL,KP,KP1,KP2,KPM1,KPT,KSH,L,N1,N2,NK,NK1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      N1 = NDIG + 1
C
C             Check for an insignificant operand.
C
      K = MA(1) - MB(1)
      IF (K.GE.NDIG+2) THEN
          DO 10 J = 1,N1
              MWA(J) = MA(J)
   10     CONTINUE
          MWA(N1+1) = 0
          RETURN

      END IF

      IF (NGUARD.LE.1) NMWA = N1 + 2
C
C             Subtract MB from MA.
C
      KP1 = MIN(N1,K+1)
      MWA(K+1) = 0
      DO 20 J = 1,KP1
          MWA(J) = MA(J)
   20 CONTINUE
      KP2 = K + 2
      IF (KP2.LE.N1) THEN
          DO 30 J = KP2,N1
              MWA(J) = MA(J) - MB(J-K)
   30     CONTINUE
      END IF

      N2 = NDIG + 2
      IF (N2-K.LE.1) N2 = 2 + K
      NK = MIN(NMWA,N1+K)
      IF (N2.LE.NK) THEN
          DO 40 J = N2,NK
              MWA(J) = -MB(J-K)
   40     CONTINUE
      END IF

      NK1 = NK + 1
      IF (NK1.LE.NMWA) THEN
          DO 50 J = NK1,NMWA
              MWA(J) = 0
   50     CONTINUE
      END IF
C
C             Normalize.  Fix the sign of any negative digit.
C
      IF (K.GT.0) THEN
          KB = NMWA - KP2 + 1
          K2 = NMWA + 1
          DO 60 J = 1,KB
              IF (MWA(K2-J).LT.0) THEN
                  MWA(K2-J) = MWA(K2-J) + JBASE
                  MWA(NMWA-J) = MWA(NMWA-J) - 1
              END IF

   60     CONTINUE
          KPT = KP2
   70     KPT = KPT - 1
          IF (MWA(KPT).LT.0 .AND. KPT.GE.3) THEN
              MWA(KPT) = MWA(KPT) + JBASE
              MWA(KPT-1) = MWA(KPT-1) - 1
              GO TO 70

          END IF

          GO TO 90

      END IF
C
      IF (K.EQ.0) THEN
          KP = N1 + 3
          KPM1 = KP - 1
          DO 80 J = 3,N1
              IF (MWA(KP-J).LT.0) THEN
                  MWA(KP-J) = MWA(KP-J) + JBASE
                  MWA(KPM1-J) = MWA(KPM1-J) - 1
              END IF

   80     CONTINUE
      END IF
C
C             Shift left if there are any leading zeros in the mantissa.
C
   90 DO 100 J = 2,NMWA
          KSH = J - 2
          IF (MWA(J).GT.0) GO TO 110
  100 CONTINUE
      MWA(1) = 0
      RETURN
C
  110 KL = NMWA - KSH
      IF (KSH.GT.0) THEN
          DO 120 J = 2,KL
              L = J + KSH
              MWA(J) = MWA(L)
  120     CONTINUE
          KL = KL + 1
          DO 130 J = KL,NMWA
              MWA(J) = 0
  130     CONTINUE
          MWA(1) = MWA(1) - KSH
      END IF
C
      RETURN

      END
      SUBROUTINE FMADDP(MA,MB,NMWA)
C
C  Internal addition routine.  MWA = MA + MB
C  The arguments are such that MA.GE.MB.GE.0.
C
C  NMWA is the number of words in MWA which will be used.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER NMWA
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,K,KB,KP,KP1,KP2,KPM1,KPT,KT,N1,N2,NK,NK1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      N1 = NDIG + 1
C
C             Check for an insignificant operand.
C
      K = MA(1) - MB(1)
      IF (K.GE.NDIG+1) THEN
          DO 10 J = 1,N1
              MWA(J) = MA(J)
   10     CONTINUE
          MWA(N1+1) = 0
          RETURN

      END IF
C
C             Add MA and MB.
C
      KP1 = K + 1
      DO 20 J = 1,KP1
          MWA(J) = MA(J)
   20 CONTINUE
      KP2 = K + 2
      IF (KP2.LE.N1) THEN
          DO 30 J = KP2,N1
              MWA(J) = MA(J) + MB(J-K)
   30     CONTINUE
      END IF

      N2 = NDIG + 2
      NK = MIN(NMWA,N1+K)
      IF (N2.LE.NK) THEN
          DO 40 J = N2,NK
              MWA(J) = MB(J-K)
   40     CONTINUE
      END IF

      NK1 = NK + 1
      IF (NK1.LE.NMWA) THEN
          DO 50 J = NK1,NMWA
              MWA(J) = 0
   50     CONTINUE
      END IF
C
C             Normalize.  Fix any digit not less than JBASE.
C
      IF (K.EQ.NDIG) RETURN
C
      IF (K.GT.0) THEN
          KB = N1 - KP2 + 1
          N2 = N1 + 1
          DO 60 J = 1,KB
              IF (MWA(N2-J).GE.JBASE) THEN
                  MWA(N2-J) = MWA(N2-J) - JBASE
                  MWA(N1-J) = MWA(N1-J) + 1
              END IF

   60     CONTINUE
          KPT = KP2
   70     KPT = KPT - 1
          IF (MWA(KPT).GE.JBASE .AND. KPT.GE.3) THEN
              MWA(KPT) = MWA(KPT) - JBASE
              MWA(KPT-1) = MWA(KPT-1) + 1
              GO TO 70

          END IF

          GO TO 90

      END IF
C
      IF (K.EQ.0) THEN
          KP = N1 + 3
          KPM1 = KP - 1
          DO 80 J = 3,N1
              IF (MWA(KP-J).GE.JBASE) THEN
                  MWA(KP-J) = MWA(KP-J) - JBASE
                  MWA(KPM1-J) = MWA(KPM1-J) + 1
              END IF

   80     CONTINUE
      END IF
C
C             Shift right if the leading digit is not less than JBASE.
C
   90 IF (MWA(2).GE.JBASE) THEN
  100     KP = NMWA + 4
          DO 110 J = 4,NMWA
              MWA(KP-J) = MWA(KP-J-1)
  110     CONTINUE
          KT = MWA(2)/JBASE
          MWA(3) = MWA(2) - KT*JBASE
          MWA(2) = KT
          MWA(1) = MWA(1) + 1
          IF (MWA(2).GE.JBASE) GO TO 100
      END IF
C
      RETURN

      END
      SUBROUTINE FMARGS(KROUTN,NARGS,MA,MB,KRESLT)
C
C  Check the input arguments to a routine for special cases.
C
C  KROUTN - ID number of the subroutine which was called
C  NARGS  - The number of input arguments (1 or 2)
C  MA     - First input argument
C  MB     - Second input argument (if NARGS is 2)
C  KRESLT - Result code returned to the calling routine.
C
C  Result codes:
C
C   0 - Perform the normal operation
C   1 - The result is the first input argument
C   2 - The result is the second input argument
C   3 - The result is -OVERFLOW
C   4 - The result is +OVERFLOW
C   5 - The result is -UNDERFLOW
C   6 - The result is +UNDERFLOW
C   7 - The result is -1.0
C   8 - The result is +1.0
C   9 - The result is -pi/2
C  10 - The result is +pi/2
C  11 - The result is 0.0
C  12 - The result is UNKNOWN
C  13 - The result is +pi
C  14 - The result is -pi/4
C  15 - The result is +pi/4
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KRESLT,KROUTN,NARGS
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER JBS,NCATMA,NCATMB,NDS
C     ..
C     .. Local Arrays ..
      INTEGER KACOS(15),KADD(15,15),KASIN(15),KATAN(15),KCOS(15),
     +        KCOSH(15),KDIV(15,15),KEXP(15),KLG10(15),KLN(15),
     +        KMPY(15,15),KPWR(15,15),KSIN(15),KSINH(15),KSQRT(15),
     +        KTAN(15),KTANH(15)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMCAT,FMIM,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
C     .. Data statements ..
C
C             These tables define the result codes to be returned for
C             given values of the input argument(s).
C
C             For example, in row 7 column 2 of this DATA statement
C             KADD(2,7) = 2 means that if the first argument in a call
C             to FMADD is in category 7 ( -UNDERFLOW ) and the second
C             argument is in category 2 ( near -OVERFLOW but
C             representable ) then the result code is 2 ( the value
C             of the sum is equal to the second input argument).
C             See routine FMCAT for descriptions of the categories.
C
C
C
C
C
      DATA KADD/3,3,3,3,3,3,3,3,3,3,3,3,3,12,12,3,0,0,0,0,0,1,1,1,0,0,0,
     +     0,0,12,3,0,0,0,0,0,1,1,1,0,0,0,0,0,4,3,0,0,0,0,0,1,1,1,0,0,0,
     +     0,0,4,3,0,0,0,0,0,1,1,1,0,0,0,0,0,4,3,0,0,0,0,0,12,1,12,0,0,
     +     0,0,0,4,3,2,2,2,2,12,12,5,12,12,2,2,2,2,4,3,2,2,2,2,2,5,2,6,
     +     2,2,2,2,2,4,3,2,2,2,2,12,12,6,12,12,2,2,2,2,4,3,0,0,0,0,0,12,
     +     1,12,0,0,0,0,0,4,3,0,0,0,0,0,1,1,1,0,0,0,0,0,4,3,0,0,0,0,0,1,
     +     1,1,0,0,0,0,0,4,3,0,0,0,0,0,1,1,1,0,0,0,0,0,4,12,0,0,0,0,0,1,
     +     1,1,0,0,0,0,0,4,12,12,4,4,4,4,4,4,4,4,4,4,4,4,4/
      DATA KMPY/4,4,4,4,12,12,12,11,12,12,12,3,3,3,3,4,0,0,0,0,0,12,11,
     +     12,0,0,1,0,0,3,4,0,0,0,0,0,12,11,12,0,0,1,0,0,3,4,0,0,0,0,0,
     +     6,11,5,0,0,1,0,0,3,12,0,0,0,0,0,6,11,5,0,0,1,0,0,12,12,0,0,0,
     +     0,0,6,11,5,0,0,1,0,0,12,12,12,12,6,6,6,6,11,5,5,5,5,12,12,12,
     +     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,5,5,5,
     +     5,11,6,6,6,6,12,12,12,12,0,0,0,0,0,5,11,6,0,0,1,0,0,12,12,0,
     +     0,0,0,0,5,11,6,0,0,1,0,0,12,3,2,2,2,2,2,5,11,6,2,2,2,2,2,4,3,
     +     0,0,0,0,0,12,11,12,0,0,1,0,0,4,3,0,0,0,0,0,12,11,12,0,0,1,0,
     +     0,4,3,3,3,3,12,12,12,11,12,12,12,4,4,4,4/
      DATA KDIV/12,12,12,4,4,4,4,12,3,3,3,3,12,12,12,12,0,0,0,0,0,4,12,
     +     3,0,0,1,0,0,12,12,0,0,0,0,0,4,12,3,0,0,1,0,0,12,6,0,0,0,0,0,
     +     4,12,3,0,0,1,0,0,5,6,0,0,0,0,0,12,12,12,0,0,1,0,0,5,6,0,0,0,
     +     0,0,12,12,12,0,0,1,0,0,5,6,6,6,6,12,12,12,12,12,12,12,5,5,5,
     +     5,11,11,11,11,11,11,11,12,11,11,11,11,11,11,11,5,5,5,5,12,12,
     +     12,12,12,12,12,6,6,6,6,5,0,0,0,0,0,12,12,12,0,0,1,0,0,6,5,0,
     +     0,0,0,0,12,12,12,0,0,1,0,0,6,5,0,0,0,0,0,3,12,4,0,0,1,0,0,6,
     +     12,0,0,0,0,0,3,12,4,0,0,1,0,0,12,12,0,0,0,0,0,3,12,4,0,0,1,0,
     +     0,12,12,12,12,3,3,3,3,12,4,4,4,4,12,12,12/
      DATA KPWR/12,12,0,5,12,12,12,8,12,12,12,3,0,12,12,12,12,0,0,12,12,
     +     12,8,12,12,12,1,0,12,12,12,12,0,0,12,12,12,8,12,12,12,1,0,12,
     +     12,12,12,0,0,12,12,12,8,12,12,12,1,0,12,12,12,12,0,0,12,12,
     +     12,8,12,12,12,1,0,12,12,12,12,0,0,12,12,12,8,12,12,12,1,0,12,
     +     12,12,12,0,3,12,12,12,8,12,12,12,5,0,12,12,12,12,12,12,12,12,
     +     12,12,11,11,11,11,11,11,11,4,4,4,4,12,12,12,8,12,12,12,6,6,6,
     +     6,4,4,0,0,0,8,8,8,8,0,0,1,0,6,6,4,4,0,0,0,8,8,8,8,0,0,1,0,6,
     +     6,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,6,6,0,0,0,8,8,8,8,8,0,1,0,4,
     +     4,6,6,0,0,0,8,8,8,8,8,0,1,0,4,4,6,6,6,6,12,12,12,8,12,12,12,
     +     4,4,4,4/
      DATA KSQRT/12,12,12,12,12,12,12,11,12,0,0,8,0,0,12/
      DATA KEXP/6,6,0,0,0,8,8,8,8,8,0,0,0,4,4/
      DATA KLN/12,12,12,12,12,12,12,12,12,0,0,11,0,0,12/
      DATA KSIN/12,12,0,0,0,0,5,11,6,0,0,0,0,12,12/
      DATA KCOS/12,12,0,0,0,8,8,8,8,8,0,0,0,12,12/
      DATA KTAN/12,12,0,0,0,0,5,11,6,0,0,0,0,12,12/
      DATA KASIN/12,12,12,9,0,0,5,11,6,0,0,10,12,12,12/
      DATA KACOS/12,12,12,13,0,10,10,10,10,10,0,11,12,12,12/
      DATA KATAN/9,9,0,14,0,0,5,11,6,0,0,15,0,10,10/
      DATA KSINH/3,3,0,0,0,1,5,11,6,1,0,0,0,4,4/
      DATA KCOSH/4,4,0,0,0,8,8,8,8,8,0,0,0,4,4/
      DATA KTANH/7,7,0,0,0,1,5,11,6,1,0,0,0,8,8/
      DATA KLG10/12,12,12,12,12,12,12,12,12,0,0,11,0,0,12/
C     ..
C
      KRESLT = 12
      KFLAG = -4
      IF (MA(1).EQ.KUNKNO) RETURN
      IF (NARGS.EQ.2) THEN
          IF (MB(1).EQ.KUNKNO) RETURN
      END IF

      KFLAG = 0
C
C             Check the validity of parameters if this is a user call.
C
      IF (NCALL.GT.1) GO TO 10
C
C             Check NDIG.
C
      IF (NDIG.LT.2 .OR. NDIG.GT.MXNDIG) THEN
          KFLAG = -1
          CALL FMWARN
          NDS = NDIG
          IF (NDIG.LT.2) NDIG = 2
          IF (NDIG.GT.MXNDIG) NDIG = MXNDIG
          WRITE (KW,FMT=9000) NDS,NDIG
          RETURN

      END IF
C
C             Check JBASE.
C
      IF (JBASE.LT.2 .OR. JBASE.GT.MXBASE) THEN
          KFLAG = -2
          CALL FMWARN
          JBS = JBASE
          IF (JBASE.LT.2) JBASE = 2
          IF (JBASE.GT.MXBASE) JBASE = MXBASE
          WRITE (KW,FMT=9010) JBS,JBASE
          RETURN

      END IF
C
C             Check exponent range.
C
      IF (MA(1).GT.MXEXP+1 .OR. MA(1).LT.-MXEXP) THEN
          IF (ABS(MA(1)).NE.KEXPOV .OR. ABS(MA(2)).NE.1) THEN
              CALL FMIM(0,MA)
              KFLAG = -3
              CALL FMWARN
              MA(1) = KUNKNO
              MA(2) = 1
              RETURN

          END IF

      END IF

      IF (NARGS.EQ.2) THEN
          IF (MB(1).GT.MXEXP+1 .OR. MB(1).LT.-MXEXP) THEN
              IF (ABS(MB(1)).NE.KEXPOV .OR. ABS(MB(2)).NE.1) THEN
                  CALL FMIM(0,MB)
                  KFLAG = -3
                  CALL FMWARN
                  MB(1) = KUNKNO
                  MB(2) = 1
                  RETURN

              END IF

          END IF

      END IF
C
C             Check for special cases.
C
   10 CALL FMCAT(MA,NCATMA)
      NCATMB = 0
      IF (NARGS.EQ.2) CALL FMCAT(MB,NCATMB)
C
      IF (KROUTN.EQ.3) THEN
          KRESLT = KADD(NCATMB,NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.29) THEN
          KRESLT = KMPY(NCATMB,NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.12) THEN
          KRESLT = KDIV(NCATMB,NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.34) THEN
          KRESLT = KPWR(NCATMB,NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.39) THEN
          KRESLT = KSQRT(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.15) THEN
          KRESLT = KEXP(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.21) THEN
          KRESLT = KLN(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.36) THEN
          KRESLT = KSIN(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.9) THEN
          KRESLT = KCOS(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.41) THEN
          KRESLT = KTAN(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.4) THEN
          KRESLT = KASIN(NCATMA)
          IF ((NCATMA.EQ.7.OR.NCATMA.EQ.9) .AND. KRAD.EQ.0) KRESLT = 12
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.2) THEN
          KRESLT = KACOS(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.5) THEN
          KRESLT = KATAN(NCATMA)
          IF ((NCATMA.EQ.7.OR.NCATMA.EQ.9) .AND. KRAD.EQ.0) KRESLT = 12
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.37) THEN
          KRESLT = KSINH(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.10) THEN
          KRESLT = KCOSH(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.42) THEN
          KRESLT = KTANH(NCATMA)
          GO TO 20

      END IF
C
      IF (KROUTN.EQ.20) THEN
          KRESLT = KLG10(NCATMA)
          GO TO 20

      END IF
C
      KRESLT = 0
      RETURN
C
   20 IF (KRESLT.EQ.12) THEN
          KFLAG = -4
          CALL FMWARN
      END IF

      IF (KRESLT.EQ.3 .OR. KRESLT.EQ.4) THEN
          IF (NCATMA.EQ.1 .OR. NCATMA.EQ.7 .OR. NCATMA.EQ.9 .OR.
     +        NCATMA.EQ.15 .OR. NCATMB.EQ.1 .OR. NCATMB.EQ.7 .OR.
     +        NCATMB.EQ.9 .OR. NCATMB.EQ.15) THEN
              KFLAG = -5

          ELSE
              KFLAG = -5
              CALL FMWARN
          END IF

      END IF

      IF (KRESLT.EQ.5 .OR. KRESLT.EQ.6) THEN
          IF (NCATMA.EQ.1 .OR. NCATMA.EQ.7 .OR. NCATMA.EQ.9 .OR.
     +        NCATMA.EQ.15 .OR. NCATMB.EQ.1 .OR. NCATMB.EQ.7 .OR.
     +        NCATMB.EQ.9 .OR. NCATMB.EQ.15) THEN
              KFLAG = -6

          ELSE
              KFLAG = -6
              CALL FMWARN
          END IF

      END IF

      RETURN

 9000 FORMAT (' NDIG WAS',I10,'.  IT HAS BEEN CHANGED TO',I10,'.')
 9010 FORMAT (' JBASE WAS',I10,'.  IT HAS BEEN CHANGED TO',I10,'.')
      END
      SUBROUTINE FMASIN(MA,MB)
C
C  MB = ARCSIN(MA)
C
C
C
C
C
C
C             Scratch array usage during FMASIN:   M01 - M06
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMATAN,FMDIV,FMENTR,FMEQU,FMEXIT,FMI2M,FMMPY,
     +         FMSQRT,FMSUB
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(4,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Use ASIN(X) = ATAN(X/SQRT(1-X*X))
C
      CALL FMI2M(1,M05)
      CALL FMSUB(M05,MB,M03)
      CALL FMADD(M05,MB,M04)
      CALL FMMPY(M03,M04,M04)
      CALL FMSQRT(M04,M04)
      CALL FMDIV(MB,M04,MB)
C
      CALL FMATAN(MB,MB)
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMATAN(MA,MB)
C
C  MB = ARCTAN(MA)
C
C
C
C
C
C
C             Scratch array usage during FMATAN:   M01 - M06
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X,XB,XM
      INTEGER J,KASAVE,KOVUN,KRESLT,KRSAVE,KST,KWSV,MA1,MA2,MXSAVE,
     +        NDSAV1,NDSAVE,NDSV
C     ..
C     .. Local Arrays ..
      INTEGER NSTACK(19)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIG,FMDIV,FMDIVI,FMDP2M,FMENTR,FMEQU,FMEXIT,
     +         FMI2M,FMM2DP,FMMPY,FMMPYI,FMPI2,FMSIN,FMSQRT,FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,LOG,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(5,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,M05,NDSAVE,NDIG)
C
C             If MA.GE.1 work with 1/MA.
C
      MA1 = MA(1)
      MA2 = MA(2)
      M05(2) = ABS(M05(2))
      IF (MA1.GE.1) THEN
          CALL FMI2M(1,MB)
          CALL FMDIV(MB,M05,M05)
      END IF
C
      KRSAVE = KRAD
      KRAD = 1
      KWSV = KWARN
C
      X = M05(1)
      XB = JBASE
      XM = MXBASE
C
C             In case pi has not been computed at the current precision
C             and will be needed here, get it to full precision first
C             to avoid repeated calls at increasing precision during
C             Newton iteration.
C
      IF (MA1.GE.1 .OR. KRSAVE.EQ.0) THEN
          IF (NJBPI.NE.JBASE .OR. NDIGPI.LT.NDIG) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMPI2(MPISAV)
              NJBPI = JBASE
              NDIGPI = NDIG
              NDIG = NDSV
          END IF

      END IF
C
C             If the argument is small, use the Taylor series,
C             otherwise use Newton iteration.
C
      IF (X*LOG(XB).LT.-5.0*LOG(XM)) THEN
          KWARN = 0
          CALL FMEQU(M05,MB,NDIG,NDIG)
          CALL FMMPY(M05,M05,M06)
          J = 3
          NDSAV1 = NDIG
C
   10     CALL FMMPY(M05,M06,M05)
          M05(2) = -M05(2)
          CALL FMDIVI(M05,J,M03)
          NDIG = NDSAV1
          CALL FMADD(MB,M03,MB)
          IF (KFLAG.EQ.1) THEN
              KFLAG = 0
              GO TO 30

          END IF

          NDIG = NDSAV1 - (MB(1)-M03(1))
          IF (NDIG.LT.2) NDIG = 2
          J = J + 2
          GO TO 10

      ELSE
C
          CALL FMM2DP(M05,X)
          X = ATAN(X)
          CALL FMDP2M(X,MB)
          CALL FMDIG(NSTACK,KST)
C
C             Newton iteration.
C
          DO 20 J = 1,KST
              NDIG = NSTACK(J)
              CALL FMSIN(MB,M06)
              CALL FMMPY(M06,M06,M03)
              CALL FMI2M(1,M04)
              CALL FMSUB(M04,M03,M03)
              CALL FMSQRT(M03,M04)
              CALL FMDIV(M06,M04,M04)
              CALL FMSUB(M04,M05,M04)
              CALL FMMPY(M03,M04,M04)
              CALL FMSUB(MB,M04,MB)
   20     CONTINUE
      END IF
C
C             If MA.GE.1 use pi/2 - ATAN(1/MA)
C
   30 IF (MA1.GE.1) THEN
          CALL FMDIVI(MPISAV,2,M06)
          CALL FMSUB(M06,MB,MB)
      END IF
C
C             Convert to degrees if necessary, round and return.
C
      KRAD = KRSAVE
      IF (KRAD.EQ.0) THEN
          CALL FMMPYI(MB,180,MB)
          CALL FMDIV(MB,MPISAV,MB)
      END IF

      IF (MA2.LT.0) MB(2) = -MB(2)
C
      IF (KFLAG.EQ.1) KFLAG = 0
      KWARN = KWSV
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMATN2(MA,MB,MC)
C
C  MC = ATAN2(MA,MB)
C
C  MC is returned as the angle between -pi and pi (or -180 and 180 if
C  degree mode is selected) for which TAN(MC) = MA/MB.  MC is an angle
C  for the point (MB,MA) in polar coordinates.
C
C
C
C
C
C
C             Scratch array usage during FMATN2:   M01 - M06
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER JQUAD,KASAVE,KOVUN,KRESLT,KWSV,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMATAN,FMDIV,FMDIVI,FMENTR,FMEQU,FMEXIT,FMI2M,FMIM,FMPI,
     +         FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(6,MA,MB,2,MC,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      KARGSW = 0
      KWSV = KWARN
      KWARN = 0
      CALL FMEQU(MA,M01,NDSAVE,NDIG)
      CALL FMEQU(MB,M02,NDSAVE,NDIG)
C
C             Check for special cases.
C
      IF (MA(1).EQ.KUNKNO .OR. MB(1).EQ.KUNKNO .OR.
     +    (MA(2).EQ.0.AND.MB(2).EQ.0)) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = -4
          GO TO 10

      END IF
C
      IF (MB(2).EQ.0 .AND. MA(2).GT.0) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(90,MC)

          ELSE
              CALL FMPI(MC)
              CALL FMDIVI(MC,2,MC)
          END IF

          GO TO 10

      END IF
C
      IF (MB(2).EQ.0 .AND. MA(2).LT.0) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(-90,MC)

          ELSE
              CALL FMPI(MC)
              CALL FMDIVI(MC,-2,MC)
          END IF

          GO TO 10

      END IF
C
      IF (MA(1).EQ.KEXPOV .AND. MB(1).LT.MXSAVE-NDIG-2) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(90,MC)

          ELSE
              CALL FMPI(MC)
              CALL FMDIVI(MC,2,MC)
          END IF

          IF (M01(2).LT.0) MC(2) = -MC(2)
          GO TO 10

      END IF
C
      IF (MA(1).EQ.KEXPUN .AND. (-MB(1)).LT.MXSAVE-NDIG-2 .AND.
     +    MB(2).LT.0) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(180,MC)

          ELSE
              CALL FMPI(MC)
          END IF

          IF (M01(2).LT.0) MC(2) = -MC(2)
          GO TO 10

      END IF
C
      IF (MB(1).EQ.KEXPOV .AND. MA(1).LT.MXSAVE-NDIG-2 .AND.
     +    MB(2).LT.0) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(180,MC)

          ELSE
              CALL FMPI(MC)
          END IF

          IF (M01(2).LT.0) MC(2) = -MC(2)
          GO TO 10

      END IF
C
      IF (MB(1).EQ.KEXPUN .AND. MA(2).EQ.0) THEN
          IF (MB(2).LT.0) THEN
              IF (KRAD.EQ.0) THEN
                  CALL FMI2M(180,MC)

              ELSE
                  CALL FMPI(MC)
              END IF

          ELSE
              CALL FMI2M(0,MC)
          END IF

          GO TO 10

      END IF
C
      IF (MB(1).EQ.KEXPUN .AND. (-MA(1)).LT.MXSAVE-NDIG-2) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(90,MC)

          ELSE
              CALL FMPI(MC)
              CALL FMDIVI(MC,2,MC)
          END IF

          IF (M01(2).LT.0) MC(2) = -MC(2)
          GO TO 10

      END IF
C
C             Determine the quadrant for the result, then use FMATAN.
C
      IF (MA(2).GE.0 .AND. MB(2).GT.0) JQUAD = 1
      IF (MA(2).GE.0 .AND. MB(2).LT.0) JQUAD = 2
      IF (MA(2).LT.0 .AND. MB(2).LT.0) JQUAD = 3
      IF (MA(2).LT.0 .AND. MB(2).GT.0) JQUAD = 4
C
      CALL FMDIV(M01,M02,MC)
      MC(2) = ABS(MC(2))
      CALL FMATAN(MC,MC)
C
      IF (JQUAD.EQ.2 .OR. JQUAD.EQ.3) THEN
          IF (KRAD.EQ.0) THEN
              CALL FMI2M(180,M05)
              CALL FMSUB(M05,MC,MC)

          ELSE
              CALL FMPI(M05)
              CALL FMSUB(M05,MC,MC)
          END IF

      END IF
C
      IF ((JQUAD.EQ.3.OR.JQUAD.EQ.4) .AND.
     +    MC(1).NE.KUNKNO) MC(2) = -MC(2)
C
C             Round the result and return.
C
   10 IF (KFLAG.EQ.1) KFLAG = 0
      KWARN = KWSV
      CALL FMEXIT(MC,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMBIG(MA)
C
C     MA = The biggest representable FM number using the current base
C          and precision.
C          The smallest positive number is then 1.0/MA.
C          Because of rounding, 1.0/(1.0/MA) will then overflow.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 7
C
      KFLAG = 0
      N1 = NDIG + 1
      DO 10 J = 2,N1
          MA(J) = JBASE - 1
   10 CONTINUE
      MA(1) = MXEXP + 1
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMCAT(MA,NCAT)
C
C  NCAT is returned as the category of MA.  This is used by the various
C  arithmetic routines to handle special cases such as:
C  'number greater than 1' + 'underflowed result' is the first argument,
C  'overflowed result' / 'overflowed result' is 'unknown'.
C
C  NCAT       range
C
C   1.         -OV                OV stands for overflowed results.
C   2.   (-OV   , -OVTH)             ( MA(1) .GE. MAXEXP+2 )
C   3.   (-OVTH ,    -1)
C   4.         -1                 OVTH stands for a representable
C   5.   (-1    , -UNTH)               number near the overflow
C   6.   (-UNTH ,   -UN)               threshold.
C   7.         -UN                     ( MA(1) .GE. MAXEXP-NDIG+1 )
C   8.          0
C   9.         +UN                UN stands for underflowed results.
C  10.   (+UN   , +UNTH)             ( MA(1) .LE. -MAXEXP-1 )
C  11.   (+UNTH ,    +1)
C  12.         +1                 UNTH stands for a representable
C  13.   (+1    , +OVTH)               number near the underflow
C  14.   (+OVTH ,   +OV)               threshold.
C  15.         +OV                     ( MA(1) .LE. -MAXEXP+NDIG-1 )
C  16.       UNKNOWN
C
C
C
C
C
C
C             Check for special symbols.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER NCAT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,MA2,NLAST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCAT = 16
      IF (MA(1).EQ.KUNKNO) RETURN
C
      IF (MA(1).EQ.KEXPOV) THEN
          NCAT = 15
          IF (MA(2).LT.0) NCAT = 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPUN) THEN
          NCAT = 9
          IF (MA(2).LT.0) NCAT = 7
          RETURN

      END IF
C
      IF (MA(2).EQ.0) THEN
          NCAT = 8
          RETURN

      END IF
C
C             Check for +1 or -1.
C
      MA2 = ABS(MA(2))
      IF (MA(1).EQ.1 .AND. MA2.EQ.1) THEN
          NLAST = NDIG + 1
          IF (NLAST.GE.3) THEN
              DO 10 J = 3,NLAST
                  IF (MA(J).NE.0) GO TO 20
   10         CONTINUE
          END IF

          NCAT = 12
          IF (MA(2).LT.0) NCAT = 4
          RETURN

      END IF
C
   20 IF (MA(1).GE.MXEXP-NDIG+1) THEN
          NCAT = 14
          IF (MA(2).LT.0) NCAT = 2
          RETURN

      END IF
C
      IF (MA(1).GE.1) THEN
          NCAT = 13
          IF (MA(2).LT.0) NCAT = 3
          RETURN

      END IF
C
      IF (MA(1).GE.-MXEXP+NDIG) THEN
          NCAT = 11
          IF (MA(2).LT.0) NCAT = 5
          RETURN

      END IF
C
      IF (MA(1).GE.-MXEXP) THEN
          NCAT = 10
          IF (MA(2).LT.0) NCAT = 6
          RETURN

      END IF
C
      RETURN

      END
      LOGICAL FUNCTION FMCOMP(MA,LREL,MB)
C
C  Logical comparison of FM numbers MA and MB.
C
C  LREL is a CHARACTER *2 description of the comparison to be done:
C  LREL = 'EQ' returns FMCOMP = .TRUE. if MA.EQ.MB
C       = 'NE', 'GE', 'GT', 'LE', 'LT' also work like a logical IF.
C
C  For comparisons involving 'UNKNOWN' or two identical special symbols
C  such as +OVERFLOW,'EQ',+OVERFLOW FMCOMP is returned FALSE and a
C  KFLAG = -4 error condition is returned.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      CHARACTER*2 LREL
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,JCOMP,NLAST
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNTRJ,FMPRNT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 8
C
      IF (NCALL.LE.LVLTRC .AND. ABS(NTRACE).GE.2) THEN
          WRITE (KW,FMT=9000)
C
          IF (NTRACE.GT.0) THEN
              CALL FMPRNT(MA)
              WRITE (KW,FMT=9010) LREL
              CALL FMPRNT(MB)

          ELSE
              CALL FMNTRJ(MA,NDIG)
              WRITE (KW,FMT=9010) LREL
              CALL FMNTRJ(MB,NDIG)
          END IF

      END IF
C
C             JCOMP will be 1 if MA.GT.MB
C                           2 if MA.EQ.MB
C                           3 if MA.LT.MB
C
C             Check for special cases.
C
      IF (LREL.NE.'EQ' .AND. LREL.NE.'NE' .AND. LREL.NE.'LT' .AND.
     +    LREL.NE.'GT' .AND. LREL.NE.'LE' .AND. LREL.NE.'GE') THEN
          KFLAG = -4
          FMCOMP = .FALSE.
          IF (KWARN.LE.0) GO TO 30
          WRITE (KW,FMT=9020) LREL
          GO TO 30

      END IF
C
      IF (MA(1).EQ.KUNKNO .OR. MB(1).EQ.KUNKNO) THEN
          KFLAG = -4
          FMCOMP = .FALSE.
          GO TO 30

      END IF
C
      IF (ABS(MA(1)).EQ.KEXPOV .AND. MA(1).EQ.MB(1) .AND.
     +    MA(2).EQ.MB(2)) THEN
          KFLAG = -4
          FMCOMP = .FALSE.
          IF (KWARN.LE.0) GO TO 30
          WRITE (KW,FMT=9030)
          GO TO 30

      END IF
C
C             Check for zero.
C
      KFLAG = 0
      IF (MA(2).EQ.0) THEN
          JCOMP = 2
          IF (MB(2).LT.0) JCOMP = 1
          IF (MB(2).GT.0) JCOMP = 3
          GO TO 20

      END IF

      IF (MB(2).EQ.0) THEN
          JCOMP = 1
          IF (MA(2).LT.0) JCOMP = 3
          GO TO 20

      END IF
C             Check for opposite signs.
C
      IF (MA(2).GT.0 .AND. MB(2).LT.0) THEN
          JCOMP = 1
          GO TO 20

      END IF

      IF (MB(2).GT.0 .AND. MA(2).LT.0) THEN
          JCOMP = 3
          GO TO 20

      END IF
C
C             See which one is larger in absolute value.
C
      IF (MA(1).GT.MB(1)) THEN
          JCOMP = 1
          GO TO 20

      END IF

      IF (MB(1).GT.MA(1)) THEN
          JCOMP = 3
          GO TO 20

      END IF

      NLAST = NDIG + 1
C
      DO 10 J = 2,NLAST
          IF (ABS(MA(J)).GT.ABS(MB(J))) THEN
              JCOMP = 1
              GO TO 20

          END IF

          IF (ABS(MB(J)).GT.ABS(MA(J))) THEN
              JCOMP = 3
              GO TO 20

          END IF

   10 CONTINUE
C
      JCOMP = 2
C
C             Now match the JCOMP value to the requested comparison.
C
   20 IF (JCOMP.EQ.1 .AND. MA(2).LT.0) THEN
          JCOMP = 3

      ELSE IF (JCOMP.EQ.3 .AND. MB(2).LT.0) THEN
          JCOMP = 1
      END IF
C
      FMCOMP = .FALSE.
      IF (JCOMP.EQ.1 .AND. (LREL.EQ.'GT'.OR.LREL.EQ.'GE'.OR.
     +    LREL.EQ.'NE')) FMCOMP = .TRUE.
C
      IF (JCOMP.EQ.2 .AND. (LREL.EQ.'EQ'.OR.LREL.EQ.'GE'.OR.
     +    LREL.EQ.'LE')) FMCOMP = .TRUE.
C
      IF (JCOMP.EQ.3 .AND. (LREL.EQ.'NE'.OR.LREL.EQ.'LT'.OR.
     +    LREL.EQ.'LE')) FMCOMP = .TRUE.
C
   30 IF (NTRACE.NE.0) THEN
          IF (NCALL.LE.LVLTRC .AND. ABS(NTRACE).GE.1) THEN
              IF (KFLAG.EQ.0) THEN
                  WRITE (KW,FMT=9040) NCALL,JBASE,NDIG

              ELSE
                  WRITE (KW,FMT=9050) NCALL,JBASE,NDIG,KFLAG
              END IF

              IF (FMCOMP) THEN
                  WRITE (KW,FMT=9060)

              ELSE
                  WRITE (KW,FMT=9070)
              END IF

          END IF

      END IF

      NCALL = NCALL - 1
      RETURN

 9000 FORMAT (' INPUT TO FMCOMP')
 9010 FORMAT (7X,'.',A2,'.')
 9020 FORMAT (/' ERROR OF TYPE KFLAG = -4 IN FM PACKAGE IN ROUTINE',
     +       ' FMCOMP',//1X,A,' IS NOT ONE OF THE SIX RECOGNIZED',' CO',
     +       'MPARISONS.',//' .FALSE. HAS BEEN RETURNED.',/)
 9030 FORMAT (/' ERROR OF TYPE KFLAG = -4 IN FM PACKAGE IN ROUTINE',
     +       ' FMCOMP',//' TWO NUMBERS IN THE SAME OVERFLOW OR',' UNDE',
     +       'RLOW CATEGORY CANNOT BE COMPARED.',//' .FALSE. HAS BEEN ',
     +       'RETURNED.',/)
 9040 FORMAT (' FMCOMP',15X,'CALL LEVEL =',I2,5X,'JBASE =',I10,5X,'NDI',
     +       'G =',I6)
 9050 FORMAT (' FMCOMP',6X,'CALL LEVEL =',I2,4X,'JBASE =',I10,4X,'NDIG',
     +       ' =',I6,4X,'KFLAG =',I3)
 9060 FORMAT (7X,'.TRUE.')
 9070 FORMAT (7X,'.FALSE.')
      END
      SUBROUTINE FMCOS(MA,MB)
C
C  MB = COS(MA)
C
C
C
C
C
C
C
C             Scratch array usage during FMCOS:   M01 - M04
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER JCOS,JSIN,JSWAP,KASAVE,KOVUN,KRESLT,MXSAVE,NDSAVE,NDSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMCOS2,FMDIVI,FMENTR,FMEQU,FMEXIT,FMMPY,FMPI2,FMRDC,
     +         FMSIN2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(9,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Reduce the argument, convert to radians if the input is
C             in degrees, and evaluate the function.
C
      CALL FMRDC(MB,MB,JSIN,JCOS,JSWAP)
      IF (MB(1).EQ.KUNKNO) GO TO 10
      IF (KRAD.EQ.0) THEN
          IF (NJBPI.NE.JBASE .OR. NDIGPI.LT.NDIG) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMPI2(MPISAV)
              NJBPI = JBASE
              NDIGPI = NDIG
              NDIG = NDSV
          END IF

          CALL FMMPY(MB,MPISAV,MB)
          CALL FMDIVI(MB,180,MB)
      END IF

      IF (MB(1).NE.KUNKNO) THEN
          IF (JSWAP.EQ.0) THEN
              CALL FMCOS2(MB,MB)

          ELSE
              CALL FMSIN2(MB,MB)
          END IF

      END IF
C
C             Append the sign, round, and return.
C
      IF (JCOS.EQ.-1) MB(2) = -MB(2)
   10 CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMCOS2(MA,MB)
C
C  Internal routine for computing MB = COS(MA) where 0.LE.MA.LE.1.
C
C
C
C
C
C
C             Scratch array usage during FMCOS2:   M01 - M03
C
C
C             Use COS(X) = SQRT(1-SIN(X)**2).
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMI2M,FMMPY,FMSIN2,FMSQRT,FMSUB
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMSIN2(MA,MB)
      CALL FMMPY(MB,MB,MB)
      CALL FMI2M(1,M02)
      CALL FMSUB(M02,MB,M02)
      CALL FMSQRT(M02,MB)
      RETURN

      END
      SUBROUTINE FMCOSH(MA,MB)
C
C  MB = COSH(MA)
C
C
C
C
C
C
C             Scratch array usage during FMCOSH:   M01 - M03
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,KWRNSV,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIV,FMDIVI,FMENTR,FMEQU,FMEXIT,FMEXP,FMI2M
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(10,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      KWRNSV = KWARN
      KWARN = 0
      KARGSW = 0
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
      CALL FMEXP(MB,MB)
      IF (MB(1).EQ.KEXPOV) THEN
          GO TO 10

      ELSE IF (MB(1).EQ.KEXPUN) THEN
          MB(1) = KEXPOV
          GO TO 10

      END IF

      CALL FMI2M(1,M01)
      CALL FMDIV(M01,MB,M01)
      CALL FMADD(MB,M01,MB)
      CALL FMDIVI(MB,2,MB)
C
C             Round and return.
C
   10 KWARN = KWRNSV
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMDIG(NSTACK,KST)
C
C  Compute the number of intermediate digits to be used in Newton
C  iteration.  This assumes that a starting approximation which is
C  accurate to double precision is used, and the root is simple.
C
C  KST is the number of iterations needed for final accuracy NDIG.
C  NSTACK(J) holds the value of NDIG to be used for the
C            Jth iteration.
C
C
C
C
C
C
C             Compute the approximate machine precision.  This will be
C             the initial accuracy obtained by using double precision
C             to compute the first approximation before starting Newton
C             iteration.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KST
C     ..
C     .. Array Arguments ..
      INTEGER NSTACK(19)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ONE,XBASE,Y,YT
      INTEGER J,JT,L,ND,NDT,NE
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      XBASE = JBASE
      ONE = 1.0
      Y = ONE
      NE = 1
C
   10 Y = Y/XBASE
      YT = ONE + Y
      IF (YT.GT.ONE) THEN
          NE = NE + 1
          GO TO 10

      END IF
C
C             Fill the intermediate digit stack (backwards).
C
      KST = 1
      ND = NDIG
      NSTACK(1) = ND
      IF (ND.LT.2*NE .OR. ND.LE.2) RETURN
C
   20 Y = ND
C
C             The 1.9 accounts for the fact that the number of correct
C             digits approximately doubles at each iteration.
C
      NDT = Y/1.9
      IF (2*NDT.LE.ND) NDT = NDT + 1
      ND = NDT
      KST = KST + 1
      NSTACK(KST) = ND
      IF (ND.GE.2*NE .AND. ND.GT.2) GO TO 20
C
C             Reverse the stack.
C
      L = KST/2
      DO 30 J = 1,L
          JT = NSTACK(J)
          NSTACK(J) = NSTACK(KST+1-J)
          NSTACK(KST+1-J) = JT
   30 CONTINUE
C
      RETURN

      END
      SUBROUTINE FMDIM(MA,MB,MC)
C
C  MC = DIM(MA,MB)
C
C  Positive difference.  MC = MA - MB  if MA.GE.MB,
C                           = 0        otherwise.
C
C
C
C
C
C
C             Scratch array usage during FMDIM:   M01 - M02
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,KWSV,MXSAVE,NDSAVE
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMENTR,FMEQU,FMEXIT,FMI2M,FMSUB
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(11,MA,MB,2,MC,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
      KARGSW = 0
      KWSV = KWARN
      KWARN = 0
      MXEXP = MXSAVE
C
      CALL FMEQU(MA,M01,NDSAVE,NDIG)
      CALL FMEQU(MB,M02,NDSAVE,NDIG)
C
      IF (FMCOMP(M01,'LT',M02)) THEN
          CALL FMI2M(0,MC)

      ELSE
          CALL FMSUB(M01,M02,MC)
      END IF
C
      IF (KFLAG.EQ.1) KFLAG = 0
      KWARN = KWSV
      CALL FMEXIT(MC,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMDIV(MA,MB,MC)
C
C  MC = MA / MB
C
C  This routine performs the trace printing for division.
C  FMDIV2 is used to do the arithmetic.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV2,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 12
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      CALL FMDIV2(MA,MB,MC)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMDIV2(MA,MB,MC)
C
C  Internal division routine.  MC = MA / MB
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B,X,XB,XBASE,XBLOG,XMWA
      INTEGER J,JB,JL,K,K1,KA,KB,KC,KDIFF,KL,KPTMWA,KRESLT,KT,LMAX,LQ,
     +        LT,MA2,MAXMWA,MB1,MB2,N1,N3,NG,NG1,NGAP,NGUARD,NL,NMBWDS,
     +        NONZPT,NORMPT,NZDMB
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMIM,FMMOVE,FMRND,FMRSLT,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (KARGSW.NE.1) THEN
          CALL FMARGS(12,2,MA,MB,KRESLT)
          IF (KRESLT.NE.0) THEN
              CALL FMRSLT(MA,MB,MC,KRESLT)
              RETURN

          END IF

      ELSE
          IF (MB(2).EQ.0) THEN
              CALL FMIM(0,MC)
              MC(1) = KUNKNO
              MC(2) = 1
              KFLAG = -4
              CALL FMWARN
              RETURN

          END IF

          IF (MA(2).EQ.0) THEN
              CALL FMIM(0,MC)
              RETURN

          END IF

      END IF

      KFLAG = 0
C
C             NGUARD is the number of guard digits used.
C
      IF (NCALL.GT.1) THEN
          NGUARD = 1

      ELSE
          B = JBASE
          NGUARD = 5.0*LOG(10.0)/LOG(B) + 2.0
      END IF

      N1 = NDIG + 1
      NG = N1 + NGUARD
C
C             Copy MA into the working array.
C
      DO 10 J = 3,N1
          MWA(J+1) = MA(J)
   10 CONTINUE
      MWA(1) = MA(1) - MB(1) + 1
      MWA(2) = 0
      N3 = NDIG + 3
      NL = N1 + NGUARD + 2
      DO 20 J = N3,NL
          MWA(J) = 0
   20 CONTINUE
C
C             Save the sign of MA and MB and then work only with
C             positive numbers.
C
      MA2 = MA(2)
      MB1 = MB(1)
      MB2 = MB(2)
      MA(2) = ABS(MA(2))
      MWA(3) = MA(2)
      MB(1) = 0
      MB(2) = ABS(MB(2))
C
C             NMBWDS is the number of words of MB used to compute
C             the trial quotient digit LT.
C
      NMBWDS = 2
      IF (JBASE.LT.100) NMBWDS = 4
C
C             XB is an approximation of MB used in selecting
C             trial quotients.
C
      XBASE = JBASE
      XB = 0
      JL = NMBWDS + 1
      DO 30 J = 2,JL
          IF (J.LE.N1) THEN
              XB = XB*XBASE + MB(J)

          ELSE
              XB = XB*XBASE
          END IF

   30 CONTINUE
      IF (JL+1.LE.N1) XB = XB + MB(JL+1)/XBASE
      XBLOG = LOG(XBASE)
C
C             LMAX determines when to normalize all of MWA.
C             NGAP determines how many digits must be normalized for
C             each trial quotient.
C
      IF (MXBASE**2/ (JBASE-1)**2.GE.NL) THEN
          LMAX = NL* (JBASE-1)

      ELSE
          LMAX = MXBASE**2/ (JBASE-1)
      END IF

      X = LMAX
      NGAP = LOG(X)/XBLOG + 2.0
C
C             Count the trailing zero digits of MB.
C
      DO 40 J = 1,NDIG
          IF (MB(NDIG+2-J).NE.0) THEN
              NZDMB = J - 1
              GO TO 50

          END IF

   40 CONTINUE
C
C             MAXMWA is an upper bound on the size of values in MWA
C             divided by JBASE-1.  It is used to determine whether
C             normalization can be postponed.
C
   50 MAXMWA = 0
C
C             KPTMWA points to the next digit in the quotient.
C
      KPTMWA = 2
C
C             This is the start of the division loop.
C
C             LT is the trial quotient digit.
C             Divide the first (NMBWDS+1) digits of MWA by MB so that
C             LT has a high probability of being the correct value.
C
   60 KL = KPTMWA + NMBWDS
      XMWA = 0
      DO 70 J = KPTMWA,KL
          IF (J.LE.NL) THEN
              XMWA = XMWA*XBASE + MWA(J)

          ELSE
              XMWA = XMWA*XBASE
          END IF

   70 CONTINUE
C
      LT = XMWA/XB
      IF (LT.GE.JBASE) LT = JBASE - 1
C
C             Subtract LT*MB from MWA.
C
      KA = KPTMWA + 1
      KB = MIN(KA+NDIG-1-NZDMB,NL)
      JB = KA - 2
      IF (LT.GT.0) THEN
          DO 80 J = KA,KB
              MWA(J) = MWA(J) - LT*MB(J-JB)
   80     CONTINUE
      END IF
C
C             Normalize enough digits in MWA so that LQ, the correct
C             quotient digit, can be determined.
C
      LQ = LT
      MAXMWA = MAXMWA + LT
C
C             NORMPT points to the last digit currently normalized.
C
      NORMPT = MIN(KA+NGAP+1,KB)
C
C             If the correct LQ has been found, then the current section
C             of MWA will be nonnegative and less than MB.
C             If MWA.LT.0 then too large a trial quotient was used and
C             a '+ MB' correction step is needed.
C             If MWA.GE.MB then too small a trial quotient was used and
C             a '- MB' correction step is needed.
C
   90 J = NORMPT
  100 IF (MWA(J).LT.0) THEN
          KT = (-MWA(J)-1)/JBASE + 1
          MWA(J-1) = MWA(J-1) - KT
          MWA(J) = MWA(J) + KT*JBASE
      END IF

      J = J - 1
      IF (J.GE.KA) GO TO 100
      IF (NORMPT.EQ.KB) MAXMWA = 0
C
C             NONZPT points to the first nonzero digit in MWA after the
C             subtraction and normalization have been done.
C
      IF (MWA(KA-1).NE.0) THEN
          NONZPT = KA - 1

      ELSE IF (MWA(KA).NE.0) THEN
          NONZPT = KA

      ELSE
          NONZPT = 0
          K1 = KA + 1
          IF (K1.LE.NORMPT) THEN
              DO 110 J = K1,NORMPT
                  IF (MWA(J).NE.0) THEN
                      NONZPT = J
                      GO TO 120

                  END IF

  110         CONTINUE
          END IF

      END IF
C
C             KDIFF points to the first digit in MWA which differs
C             from the corresponding digit in MB.
C
  120 IF (MWA(KA-1).NE.0) THEN
          KDIFF = KA - 1

      ELSE IF (MWA(KA).NE.MB(2)) THEN
          KDIFF = KA

      ELSE
          KDIFF = 0
          K1 = KA + 1
          IF (K1.LE.NORMPT) THEN
              DO 130 J = K1,NORMPT
                  IF (MWA(J).NE.MB(J-JB)) THEN
                      KDIFF = J
                      GO TO 140

                  END IF

  130         CONTINUE
          END IF

      END IF
C
C             Protect against cases like MWA = 2 0 0 0 ...
C             and MB = 1 9 9 9 ....  Since MWA may not be completely
C             normalized, MWA(KDIFF) might still be equal to
C             MB(KDIFF-JB).
C
  140 IF (KDIFF.NE.0 .AND. KDIFF.LT.NL .AND. KDIFF+1-JB.LE.N1) THEN
          IF (MWA(KDIFF+1).EQ.0 .AND. NORMPT.LT.KB .AND.
     +        MB(KDIFF+1-JB).EQ.JBASE-1) KDIFF = 0
      END IF
C
C             If NONZPT is zero then not enough digits were normalized
C             to tell whether MWA is negative.
C
      IF (NONZPT.EQ.0 .OR. NORMPT-NONZPT.LT.NGAP) THEN
          IF (NORMPT.EQ.KB) GO TO 150
          NORMPT = KB
          GO TO 90

      END IF
C
C             See if the current section of MWA is negative.
C
  150 IF (NONZPT.LE.0) GO TO 170
      IF (MWA(NONZPT).LT.0) THEN
C
C             At the end of the division LT could be off by more than
C             one because part of MB has shifted off the end of the
C             active part of MWA.  Skip the correction step for the
C             last guard digit.
C
          IF (KA.EQ.KB .AND. KB.EQ.NL) GO TO 210
C
C             Correct LQ by adding MB.
C
          LQ = LQ - 1
          MAXMWA = MAXMWA - 1
C
C             This is the only place where values greater than JBASE
C             can be introduced into MWA.  Normalize them here.
C
          K = KB
          KC = 0
  160     MWA(K) = MWA(K) + MB(K-JB) + KC
          IF (MWA(K).GE.JBASE) THEN
              KC = 1
              MWA(K) = MWA(K) - JBASE

          ELSE
              KC = 0
          END IF

          K = K - 1
          IF (K.GE.KA) GO TO 160
          IF (KC.EQ.1) MWA(KA-1) = MWA(KA-1) + 1
          GO TO 90

      END IF
C
C             If KDIFF is zero then not enough digits were normalized
C             to tell if MWA.GE.MB.
C
  170 IF (KDIFF.EQ.0) THEN
          IF (NORMPT.EQ.KB) GO TO 180
          NORMPT = KB
          GO TO 90

      END IF
C
C             See if the current section of MWA is greater than MB.
C
      IF (MWA(KDIFF).LT.MB(KDIFF-JB)) GO TO 210
C
C             If NORMPT-KDIFF.GE.NGAP then not enough digits were
C             normalized to tell if MWA.GE.MB.
C
      IF (NORMPT-KDIFF.LT.NGAP .AND. NORMPT.LT.KB) THEN
          NORMPT = KB
          GO TO 90

      END IF
C
C             Correct LQ by subtracting MB.
C
  180 IF (KA.EQ.KB .AND. KB.EQ.NL) GO TO 210
      LQ = LQ + 1
C
C             Due to part of MB shifting off the end of the active
C             part of MWA, LQ could be equal to JBASE here.  This
C             means the correct value of LQ is JBASE-1 for all
C             subsequent digits.
C
      IF (LQ.GE.JBASE) THEN
          NG1 = NG + 1
          DO 190 K = KPTMWA,NG1
              MWA(K) = JBASE - 1
  190     CONTINUE
          GO TO 230

      END IF
C
      MAXMWA = MAXMWA + 1
      DO 200 K = KA,KB
          MWA(K) = MWA(K) - MB(K-JB)
  200 CONTINUE
      GO TO 90
C
C             Here the correct value for LQ has been found.
C
  210 MWA(KPTMWA) = LQ
C
C             See if MWA must be normalized.
C
      IF (MAXMWA+JBASE-1.GE.LMAX) THEN
          J = KB
  220     IF (MWA(J).LT.0) THEN
              KT = (-MWA(J)-1)/JBASE + 1
              MWA(J-1) = MWA(J-1) - KT
              MWA(J) = MWA(J) + KT*JBASE
          END IF

          J = J - 1
          IF (J.GE.KA) GO TO 220
          MAXMWA = 0
      END IF
C
      KPTMWA = KPTMWA + 1
      IF (KPTMWA.LE.NG) GO TO 60
      IF (MWA(2).EQ.0 .AND. KPTMWA.LE.NG+1) GO TO 60
C
C             Round, affix the sign, and return.
C
  230 MA(2) = MA2
      MB(1) = MB1
      MB(2) = MB2
      IF (MWA(2).EQ.0) THEN
          CALL FMRND(NDIG,NGUARD,1)

      ELSE
          CALL FMRND(NDIG,NGUARD,0)
      END IF

      CALL FMMOVE(MC)
      IF (MA2*MB2.LT.0) MC(2) = -MC(2)
      RETURN

      END
      SUBROUTINE FMDIVI(MA,INT,MB)
C
C  MB = MA / INT
C
C  Divide FM number MA by one word integer INT.
C
C  This routine is faster than FMDIV when the divisor is less than
C  MXBASE.  When INT is not less than MXBASE, FMDIV2 is used.  In
C  this case if INT is known to be a product of two integers less than
C  MXBASE it is usually faster to make two calls to FMDIVI with half-
C  word factors than one call with their product.
C
C
C
C
C
C
C             Scratch array usage during FMDIVI:   M01
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER INTP,J,KA,KB,KL,KPT,KPTWA,KT,MA1,MA2,MODINT,N1,NGUARD
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV2,FMIM,FMMOVE,FMNTR,FMNTRI,FMRND,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 13
      IF (NTRACE.NE.0) THEN
          CALL FMNTR(2,MA,MA,1)
          CALL FMNTRI(2,INT,0)
      END IF

      KFLAG = 0
      N1 = NDIG + 1
C
C             Check for special cases.
C
      IF (MA(1).EQ.KUNKNO .OR. INT.EQ.0) THEN
          MA1 = MA(1)
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          IF (MA1.NE.KUNKNO) CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (ABS(INT).EQ.1) THEN
          DO 10 J = 1,N1
              MB(J) = MA(J)
   10     CONTINUE
          MB(2) = MA(2)*INT
          IF (MA(1).EQ.KEXPOV) KFLAG = -5
          IF (MA(1).EQ.KEXPUN) KFLAG = -6
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(2).EQ.0) THEN
          DO 20 J = 1,N1
              MB(J) = 0
   20     CONTINUE
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPUN) THEN
          MA2 = MA(2)
          CALL FMIM(0,MB)
          MB(1) = KEXPUN
          MB(2) = 1
          IF ((MA2.LT.0.AND.INT.GT.0) .OR.
     +        (MA2.GT.0.AND.INT.LT.0)) MB(2) = -1
          KFLAG = -6
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPOV) THEN
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             NGUARD is the number of guard digits used.
C
      IF (NCALL.GT.1) THEN
          NGUARD = 1

      ELSE
          B = JBASE
          NGUARD = 5.0*LOG(10.0)/LOG(B) + 2.0
      END IF

      N1 = NDIG + 1
C
C             If ABS(INT).GE.MXBASE use FMDIV.
C
      IF (ABS(INT).GT.MXBASE) THEN
          CALL FMIM(INT,M01)
          CALL FMDIV2(MA,M01,MB)
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Work with positive numbers.
C
      MA2 = MA(2)
      MA(2) = ABS(MA(2))
      INTP = ABS(INT)
C
C             Find the first significant digit of the quotient.
C
      KT = 0
      DO 30 J = 2,N1
          KPT = J
          KT = KT*JBASE + MA(J)
          IF (KT.GE.INTP) GO TO 50
   30 CONTINUE
C
   40 KPT = KPT + 1
      KT = KT*JBASE
      IF (KT.LT.INTP) GO TO 40
C
C             Do the rest of the division.
C
   50 KA = KPT + 1
      MWA(1) = MA(1) + 2 - KPT
      MWA(2) = KT/INTP
      MODINT = MOD(KT,INTP)
      KPTWA = 2
      IF (KA.LE.N1) THEN
          KL = 3 - KA
          DO 60 J = KA,N1
              KT = MODINT*JBASE + MA(J)
              MWA(KL+J) = KT/INTP
              MODINT = MOD(KT,INTP)
   60     CONTINUE
          KPTWA = KL + N1
      END IF
C
      KA = KPTWA + 1
      KB = N1 + NGUARD
      IF (KA.LE.KB) THEN
          DO 70 J = KA,KB
              KT = MODINT*JBASE
              MWA(J) = KT/INTP
              MODINT = MOD(KT,INTP)
   70     CONTINUE
      END IF
C
C             Round the result, put the sign on MB and return.
C
      MA(2) = MA2
      CALL FMRND(NDIG,NGUARD,0)
      CALL FMMOVE(MB)
      IF ((MA2.LT.0.AND.INT.GT.0) .OR.
     +    (MA2.GT.0.AND.INT.LT.0)) MB(2) = -MB(2)
      IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMDM(X,MA)
C
C  Internal routine for converting double precision to multiple
C  precision.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ONE,XBASE,Y,YT
      INTEGER J,K,MN,N1,NE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMIM2,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      N1 = NDIG + 1
      ONE = 1.0
      XBASE = JBASE
      Y = ONE
      K = 0
      NE = 2
C
C             Find the approximate machine precision.
C             NE-1 is the number of words at the current precision and
C             base roughly equal to machine precision.
C
      DO 10 J = 2,NDIG
          Y = Y/XBASE
          NE = NE + 1
          YT = ONE + Y
          IF (YT.LE.ONE) GO TO 20
   10 CONTINUE
   20 Y = X
      IF (X.LT.0.0) Y = -X
C
      IF (X.EQ.0.0) THEN
          DO 30 J = 1,N1
              MA(J) = 0
   30     CONTINUE
          GO TO 140

      END IF
C
C             Get the exponent.
C
      IF (Y.GT.ONE) THEN
          IF (Y/XBASE.LT.Y) THEN
   40         K = K + 1
              Y = Y/XBASE
              IF (Y.GT.ONE) GO TO 40
              IF (Y.LT.ONE) THEN
                  MA(1) = K
                  GO TO 80

              END IF

              GO TO 60

          ELSE
              CALL FMIM2(0,MA)
              MA(1) = KUNKNO
              MA(2) = 1
              KFLAG = -4
              CALL FMWARN
              RETURN

          END IF

      END IF
C
      IF (Y.LT.ONE) THEN
          IF (Y*XBASE.GT.Y) THEN
   50         K = K - 1
              Y = Y*XBASE
              IF (Y.LT.ONE) GO TO 50
              IF (Y.GT.ONE) THEN
                  K = K + 1
                  Y = Y/XBASE
                  MA(1) = K
                  GO TO 80

              END IF

          ELSE
              CALL FMIM2(0,MA)
              MA(1) = KUNKNO
              MA(2) = 1
              KFLAG = -4
              CALL FMWARN
              RETURN

          END IF

      END IF
C
   60 MA(1) = K + 1
      MA(2) = 1
      DO 70 J = 3,N1
          MA(J) = 0
   70 CONTINUE
      GO TO 140
C
C             Build the rest of the number.
C
   80 DO 90 J = 2,NE
          Y = Y*XBASE
          K = Y
          YT = K
          Y = Y - YT
          MA(J) = K
          IF (J.GE.N1) GO TO 110
   90 CONTINUE
      K = NE + 1
      DO 100 J = K,N1
          MA(J) = 0
  100 CONTINUE
C
C             Normalize.
C
  110 IF (ABS(MA(2)).GE.JBASE) THEN
          K = N1 + 1
          DO 120 J = 3,N1
              K = K - 1
              MA(K) = MA(K-1)
  120     CONTINUE
          MN = MA(2)/JBASE
          MA(3) = MA(2) - MN*JBASE
          MA(2) = MN
          MA(1) = MA(1) + 1
          GO TO 140

      END IF
C
      IF (MA(2).EQ.0) THEN
          DO 130 J = 2,NDIG
              MA(J) = MA(J+1)
  130     CONTINUE
          MA(1) = MA(1) - 1
          MA(N1) = 0
      END IF
C
  140 IF (X.LT.0.0) MA(2) = -MA(2)
      RETURN

      END
      SUBROUTINE FMDP2M(X,MA)
C
C  MA = X
C
C  Convert a double precision floating point number to FM format.
C
C  In general the relative accuracy of the number returned is only
C  the relative accuracy of a machine precision number.  This may be
C  true even if X can be represented exactly in the machine floating
C  point number system.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDM,FMNTR,FMNTRR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 14
      IF (NTRACE.NE.0) CALL FMNTRR(2,X,1)
C
      CALL FMDM(X,MA)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMENTR(NROUTN,MA,MB,NARGS,MC,KRESLT,NDSAVE,MXSAVE,
     +                  KASAVE,KOVUN)
C
C  Do the argument checking and increasing of precision, overflow
C  threshold, etc., upon entry to an FM routine.
C
C  NROUTN - routine number of calling routine
C  MA     - first input argument
C  MB     - second input argument (optional)
C  NARGS  - number of input arguments
C  MC     - result argument
C  KRESLT - returned nonzero if the input arguments give the result
C           immediately (e.g., MA*0 or OVERFLOW*MB)
C  NDSAVE - saves the value of NDIG after NDIG is increased
C  MXSAVE - saves the value of MXEXP
C  KASAVE - saves the value of KARGSW
C  KOVUN  - returned nonzero if an input argument is (+ or -) overflow
C           or underflow.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KASAVE,KOVUN,KRESLT,MXSAVE,NARGS,NDSAVE,NROUTN
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER K
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMDIVI,FMEQU,FMI2M,FMNTR,FMPI,FMRSLT,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,MAX
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = NROUTN
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,NARGS)
      CALL FMARGS(NROUTN,NARGS,MA,MB,KRESLT)
C
      KOVUN = 0
      IF (MA(1).EQ.KEXPOV .OR. MA(1).EQ.KEXPUN) KOVUN = 1
      IF (NARGS.EQ.2) THEN
          IF (MB(1).EQ.KEXPOV .OR. MB(1).EQ.KEXPUN) KOVUN = 1
      END IF
C
C             Increase the working precision.
C
      NDSAVE = NDIG
      IF (NCALL.EQ.1) THEN
          B = JBASE
          K = 5.0*LOG(10.0)/LOG(B) + 2.0
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              KRESLT = 12
          END IF

      END IF
C
      IF (KRESLT.NE.0) THEN
          IF (KRESLT.EQ.9 .OR. KRESLT.EQ.10 .OR. KRESLT.GE.13) THEN
              IF (KRAD.EQ.1) THEN
                  CALL FMPI(MC)

              ELSE
                  CALL FMI2M(180,MC)
              END IF

              IF (KRESLT.LE.10) CALL FMDIVI(MC,2,MC)
              IF (KRESLT.GE.14) CALL FMDIVI(MC,4,MC)
              CALL FMEQU(MC,MC,NDIG,NDSAVE)
              NDIG = NDSAVE
              IF (KRESLT.EQ.9 .OR. KRESLT.EQ.14) MC(2) = -MC(2)
              IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
              NCALL = NCALL - 1
              RETURN

          END IF
C
          NDIG = NDSAVE
          CALL FMRSLT(MA,MB,MC,KRESLT)
          IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Turn off argument checking for low-level routines while
C             a high-level call is in progress.
C
      KASAVE = KARGSW
      KARGSW = 1
C
C             Extend the overflow/underflow threshold.
C
      MXSAVE = MXEXP
      MXEXP = MXEXP2
      RETURN

      END
      SUBROUTINE FMEQU(MA,MB,NDA,NDB)
C
C  Set MB (having NDB digits) equal to MA (having NDA digits).
C
C  If MB has less precision than MA the result is rounded to NDB digits.
C
C  If MB has more precision the result has zero digits padded on the
C  right.
C
C
C
C
C
C
C             See if MA and MB are the same array.  If so, the loop
C             which copies MA to MB can be skipped.
C
C             The MIN calls try to obscure the test so that an
C             optimizing compiler will not remove the K = ... and
C             MA(3) = ... lines and sometimes change an input argument.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDA,NDB
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,K,KB,KSAME,KWRNSV,L,MB2,N1,NDA1,NDA2,NDB2,NDG
C     ..
C     .. External Subroutines ..
      EXTERNAL FMTRAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KSAME = 0
      K = MIN(MA(3),MXBASE)
      MB(3) = K + 1
      IF (MA(3).EQ.MB(3)) KSAME = 1
      MA(3) = MIN(K,MXBASE+MXEXP)
C
C             Check for special symbols.
C
      KFLAG = 0
      IF (ABS(MA(1)).GE.KEXPOV) THEN
          DO 10 J = 2,NDB
              MB(J+1) = 0
   10     CONTINUE
          MB(1) = MA(1)
          MB(2) = MA(2)
          GO TO 140

      END IF
C
      IF (NDB.EQ.NDA) GO TO 90
C
      IF (NDB.GT.NDA) GO TO 110
C
C             Round to NDB digits.
C
      NDG = NDB
      N1 = NDB + 1
      IF (KSAME.EQ.0) THEN
          DO 20 J = 1,N1
              MB(J) = MA(J)
   20     CONTINUE
      END IF

      IF (NDG.LT.1 .OR. (KROUND.EQ.0.AND.NCALL.LE.1)) GO TO 140
C
      L = NDB + 2
      IF (2* (MA(L)+1).LT.JBASE) GO TO 140
      IF (MOD(JBASE,2).EQ.0) THEN
          IF (2*MA(L).LT.JBASE) GO TO 140
          IF (2*MA(L).EQ.JBASE) THEN
              IF (L.LE.NDA) THEN
                  DO 30 J = L,NDA
                      IF (MA(J+1).GT.0) GO TO 50
   30             CONTINUE
              END IF
C
C                       Round to even.
C
              IF (MOD(MB(N1),2).EQ.0) GO TO 140
          END IF

      ELSE
          IF (2*MA(L)+1.EQ.JBASE) THEN
              IF (L.LE.NDA) THEN
                  DO 40 J = L,NDA
                      IF (2* (MA(J+1)+1).LT.JBASE) GO TO 140
                      IF (2*MA(J+1).GT.JBASE) GO TO 50
   40             CONTINUE
                  GO TO 140

              END IF

          END IF

      END IF
C
   50 MB(NDG+1) = MB(NDG+1) + 1
      MB(NDG+2) = 0
C
C             Check whether there was a carry in the rounded digit.
C
      MB2 = MB(2)
      MB(2) = ABS(MB(2))
      KB = NDG + 1
      IF (KB.GE.3) THEN
          K = KB + 1
          DO 60 J = 3,KB
              K = K - 1
              IF (MB(K).LT.JBASE) GO TO 80
              MB(K-1) = MB(K-1) + MB(K)/JBASE
              MB(K) = MOD(MB(K),JBASE)
   60     CONTINUE
      END IF
C
C             If there is a carry in the first digit then the exponent
C             must be adjusted and the number shifted right.
C
      IF (MB(2).LT.JBASE) GO TO 80
      IF (KB.GE.4) THEN
          K = KB + 1
          DO 70 J = 4,KB
              K = K - 1
              MB(K) = MB(K-1)
   70     CONTINUE
      END IF
C
      IF (KB.GE.3) MB(3) = MOD(MB(2),JBASE)
      MB(2) = MB(2)/JBASE
      MB(1) = MB(1) + 1
C
   80 IF (MB2.LT.0) MB(2) = -MB(2)
      GO TO 140
C
C             MA and MB have the same precision.
C
   90 IF (KSAME.EQ.0) THEN
          NDA1 = NDA + 1
          DO 100 J = 1,NDA1
              MB(J) = MA(J)
  100     CONTINUE
      END IF

      GO TO 140
C
C             Extend to NDB digits by padding with zeros.
C
  110 IF (KSAME.EQ.0) THEN
          NDA1 = NDA + 1
          DO 120 J = 1,NDA1
              MB(J) = MA(J)
  120     CONTINUE
      END IF

      NDA2 = NDA + 2
      NDB2 = NDB + 1
      DO 130 J = NDA2,NDB2
          MB(J) = 0
  130 CONTINUE
C
C             Check for overflow or underflow.
C
  140 IF (ABS(MB(1)).GT.MXEXP) THEN
          IF (MB(1).NE.KUNKNO .OR. MB(2).NE.1) THEN
              KWRNSV = KWARN
              KWARN = 0
              NCALL = NCALL + 1
              CALL FMTRAP(MB)
              NCALL = NCALL - 1
              KWARN = KWRNSV
          END IF

          IF (MB(1).EQ.KUNKNO) KFLAG = -4
      END IF
C
      RETURN

      END
      SUBROUTINE FMEXIT(MT,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
C
C  Upon exit from an FM routine the result MT (having precision NDIG)
C  is rounded and returned in MC (having precision NDSAVE).
C  The values of NDIG, MXEXP, and KARGSW are restored.
C  KOVUN is nonzero if one of the routine's input arguments was overflow
C  or underflow.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KASAVE,KOVUN,MXSAVE,NDSAVE
C     ..
C     .. Array Arguments ..
      INTEGER MC(LUNPCK),MT(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KFSAVE,KWSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMNTR,FMWARN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KWSV = KWARN
      KWARN = 0
      MXEXP = MXSAVE
      KFSAVE = KFLAG
      CALL FMEQU(MT,MC,NDIG,NDSAVE)
      IF (KFLAG.NE.-5 .AND. KFLAG.NE.-6) KFLAG = KFSAVE
      NDIG = NDSAVE
      KWARN = KWSV
      IF (KFLAG.EQ.1) KFLAG = 0
      IF ((MC(1).EQ.KUNKNO.AND.KFLAG.NE.-9) .OR.
     +    (MC(1).EQ.KEXPUN.AND.KOVUN.EQ.0) .OR.
     +    (MC(1).EQ.KEXPOV.AND.KOVUN.EQ.0)) CALL FMWARN
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      KARGSW = KASAVE
      RETURN

      END
      SUBROUTINE FMEXP(MA,MB)
C
C  MB = EXP(MA)
C
C
C
C
C
C
C
C             Scratch array usage during FMEXP:   M01 - M03
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL XBASE,XMA,XOV,XT
      INTEGER IEXTRA,J,KASAVE,KOVUN,KRESLT,KT,KWSV,MA1,MA2,MXSAVE,
     +        NDSAVE,NDSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMENTR,FMEQU,FMEXIT,FMEXP2,FMI2M,FMIM,FMINT,FMIPWR,
     +         FMM2I,FMMPY,FMNTR,FMSUB,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,MIN,REAL
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(15,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      MA1 = MA(1)
      MA2 = MA(2)
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Check for obvious underflow or overflow.
C             XOV is LN(LN(slightly above overflow))
C             XMA is LN(LN(EXP(MA))) approximately.
C
      XOV = LOG(1.01*REAL(MXEXP)) + LOG(LOG(REAL(JBASE)))
      XMA = LOG(REAL(MAX(ABS(MA2),1))) - LOG(REAL(JBASE)) +
     +      MA1*LOG(REAL(JBASE))
C
   10 IF (XMA.GE.XOV) THEN
          CALL FMIM(0,MB)
          IF (MA2.GT.0) THEN
              KFLAG = -5
              MB(1) = KEXPOV
              MB(2) = 1

          ELSE
              KFLAG = -6
              MB(1) = KEXPUN
              MB(2) = 1
          END IF

          NDIG = NDSAVE
          MXEXP = MXSAVE
          KARGSW = KASAVE
          CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Split MA into integer and fraction parts.
C             Work with a positive argument.
C             M02 = integer part of ABS(MA)
C             MB  = fraction part of ABS(MA)
C
      MB(2) = ABS(MB(2))
      CALL FMINT(MB,M02)
      CALL FMSUB(MB,M02,MB)
C
C             If the integer part is not zero, use FMIPWR to compute
C             E**(M02).  If M02 is too large to represent as a one word
C             integer, the definition of MXEXP insures that E**(M02)
C             overflows or underflows.
C
      KWSV = KWARN
      KWARN = 0
      CALL FMM2I(M02,KT)
      KWARN = KWSV
      IF (KFLAG.NE.0) THEN
          XMA = XOV
          GO TO 10

      END IF

      IF (KT.GT.0) THEN
C
C             Compute IEXTRA, the number of extra digits required in
C             order to get EXP(KT) correct to the current precision.
C             IEXTRA = (LN(KT) - 1)/LN(JBASE)
C
          XT = KT
          XBASE = JBASE
          IEXTRA = LOG(XT)/LOG(XBASE) + 0.5
          IF (IEXTRA.GT.0 .AND. NDIG+IEXTRA.LE.MXNDG2) THEN
              DO 20 J = 1,IEXTRA
                  MB(NDIG+1+J) = 0
   20         CONTINUE
          END IF

          NDIG = NDIG + IEXTRA
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MB(1) = KUNKNO
              MB(2) = 1
              DO 30 J = 2,NDSAVE
                  MB(J+1) = 0
   30         CONTINUE
              CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
              RETURN

          END IF
C
C             Check whether the current precision of e is large
C             enough.
C
          IF (NJBE.NE.JBASE .OR. NDIG.GT.NDIGE) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMI2M(1,MESAV)
              CALL FMEXP2(MESAV,MESAV)
              NJBE = JBASE
              NDIGE = NDIG
              NDIG = NDSV
          END IF
C
      END IF
C
C             Now do the fraction part of MA and combine the results.
C
      KWSV = KWARN
      KWARN = 0
      IF (MB(2).NE.0 .AND. KT.GT.0) THEN
          CALL FMEXP2(MB,MB)
          CALL FMIPWR(MESAV,KT,M03)
          CALL FMMPY(MB,M03,MB)

      ELSE IF (MB(2).NE.0 .AND. KT.EQ.0) THEN
          CALL FMEXP2(MB,MB)

      ELSE IF (MB(2).EQ.0 .AND. KT.GT.0) THEN
          CALL FMIPWR(MESAV,KT,MB)

      ELSE
          CALL FMI2M(1,MB)
      END IF
C
C             Invert if MA was negative.
C
      IF (MA2.LT.0) THEN
          CALL FMI2M(1,M02)
          CALL FMDIV(M02,MB,MB)
      END IF

      KWARN = KWSV
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMEXP2(MA,MB)
C
C  MB = EXP(MA)
C
C  Internal exponential routine (called with 0.LT.MA.LE.1).
C
C  Upon return M03 (in COMMON /FM1/) contains an accurate value for
C  EXP(MA)-1 which can be used to avoid loss of significance in other
C  routines if MA is close to zero.
C
C
C
C
C
C
C             Scratch array usage during FMEXP2:   M01 - M03
C
C
C
C             LJSUMS = 8*LUNPCK allows for up to eight concurrent sums.
C             Increasing this value will begin to improve the speed of
C             EXP when the base is large and precision exceeds about
C             1,500 decimal digits.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA,LJSUMS
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK,
     +          LJSUMS=8*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MJSUMS(LJSUMS),
     +        MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL ALOG2,ALOGB,ALOGT,B,T,XN
      INTEGER J,J2,K,K2,KPT,KTWO,L,L2,MAXVAL,N2,NBIG,NBOT,NDSAV1,NDSAVE,
     +        NTERM,NTOP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIVI,FMEQU,FMI2M,FMIPWR,FMMPY,FMMPYI,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,MAX,MIN,REAL,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSUMS/MJSUMS
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NDSAVE = NDIG
      IF (MA(1).EQ.1) THEN
C
C             Here the special case EXP(1.0) is computed.
C             Use the direct series  e = 1/0! + 1/1! + 1/2! + ...
C             Do as much of the work as possible using small integers
C             to minimize the number of FM calls.
C             Reduce NDIG while computing each term in the
C             sum as the terms get smaller.
C
          B = JBASE
          T = NDIG
          XN = T*LOG(B)/LOG(T)
          K = LOG(XN)/LOG(B)
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MB(1) = KUNKNO
              MB(2) = 1
              DO 10 J = 2,NDSAVE
                  MB(J+1) = 0
   10         CONTINUE
              NDIG = NDSAVE
              RETURN

          END IF

          NDSAV1 = NDIG
C
          CALL FMI2M(2,MB)
          CALL FMI2M(1,M02)
          J = 2
          NBIG = MXBASE
C
   20     NTOP = 1
          NBOT = J
   30     IF (NBOT.GT.NBIG/ (J+1)) GO TO 40
          J = J + 1
          NTOP = J*NTOP + 1
          NBOT = J*NBOT
          GO TO 30
C
   40     CALL FMDIVI(M02,NBOT,M02)
          IF (NTOP.GT.1) THEN
              CALL FMMPYI(M02,NTOP,M03)
              NDIG = NDSAV1
              CALL FMADD(MB,M03,MB)
              NDIG = NDSAV1 - (MB(1)-M03(1))

          ELSE
              NDIG = NDSAV1
              CALL FMADD(MB,M02,MB)
              NDIG = NDSAV1 - (MB(1)-M02(1))
          END IF

          IF (NDIG.LT.2) NDIG = 2
          IF (KFLAG.NE.1) THEN
              J = J + 1
              GO TO 20

          END IF

          NDIG = NDSAVE
          CALL FMI2M(-1,M02)
          CALL FMADD(MB,M02,M03)
          KFLAG = 0
          RETURN

      END IF
C
C             Here is the general case.  Compute EXP(MA) where
C             0 .LT. MA .LT. 1.
C
C             Use the direct series
C                  EXP(X) = 1 + X + X**2/2! + X**3/3! + ...
C
C             The argument will be halved K2 times before the series
C             is summed.  The series will be added as J2 concurrent
C             series.  The approximately optimal values of K2 and J2
C             are now computed to try to minimize the time required.
C             N2 is the approximate number of terms of the series which
C             will be needed, and L2 guard digits will be carried.
C
      B = JBASE
      K = 5.0*LOG(10.0)/LOG(B) + 2.0
      T = MAX(NDIG-K,2)
      ALOGB = LOG(B)
      ALOG2 = LOG(2.0)
      ALOGT = LOG(T)
      J2 = 0.081*ALOGB*T**0.3333 + 1.85
      J2 = MAX(1,MIN(J2,LJSUMS/MXNDG2))
      K2 = 0.7*SQRT(T*ALOGB/REAL(J2)) - 0.57*ALOGT + 0.5
C
      L = - (REAL(MA(1))*ALOGB+LOG(REAL(MA(2))/B+REAL(MA(3))/ (B*B)))/
     +    ALOG2 - 0.3
      K2 = K2 - L
      IF (K2.LT.0) THEN
          K2 = 0
          J2 = .43*SQRT(T*ALOGB/ (ALOGT+REAL(L)*ALOG2)) + .33
      END IF

      IF (J2.LE.1) J2 = 1
C
      N2 = T*ALOGB/ (ALOGT+REAL(L)*ALOG2)
      L2 = LOG(REAL(N2)+2.0**K2)/ALOGB
      NDIG = NDIG + L2
      IF (NDIG.GT.MXNDG2) THEN
          KFLAG = -9
          CALL FMWARN
          MB(1) = KUNKNO
          MB(2) = 1
          DO 50 J = 2,NDSAVE
              MB(J+1) = 0
   50     CONTINUE
          NDIG = NDSAVE
          RETURN

      END IF

      NDSAV1 = NDIG
C
C             Halve the argument K2 times.
C
      CALL FMEQU(MA,M02,NDSAVE,NDIG)
      KTWO = 1
      MAXVAL = MXBASE/2
      IF (K2.GT.0) THEN
          DO 60 J = 1,K2
              KTWO = 2*KTWO
              IF (KTWO.GT.MAXVAL) THEN
                  CALL FMDIVI(M02,KTWO,M02)
                  KTWO = 1
              END IF

   60     CONTINUE
          IF (KTWO.GT.1) CALL FMDIVI(M02,KTWO,M02)
      END IF
C
C             Sum the series X + X**2/2! + X**3/3! + ....
C             Split into J2 concurrent sums and reduce NDIG while
C             computing each term in the sum as the terms get smaller.
C
      CALL FMEQU(M02,MB,NDIG,NDIG)
      NTERM = 1
      DO 70 J = 1,J2
          CALL FMDIVI(MB,NTERM,MB)
          NTERM = NTERM + 1
          KPT = (J-1)* (NDIG+1) + 1
          CALL FMEQU(MB,MJSUMS(KPT),NDIG,NDIG)
   70 CONTINUE
      CALL FMIPWR(M02,J2,M03)
C
   80 CALL FMMPY(MB,M03,MB)
      DO 90 J = 1,J2
          CALL FMDIVI(MB,NTERM,MB)
          KPT = (J-1)* (NDSAV1+1) + 1
          NDIG = NDSAV1
          CALL FMADD(MJSUMS(KPT),MB,MJSUMS(KPT))
          IF (KFLAG.EQ.1) GO TO 100
          NDIG = NDSAV1 - (MJSUMS(KPT)-MB(1))
          IF (NDIG.LT.2) NDIG = 2
          NTERM = NTERM + 1
   90 CONTINUE
      GO TO 80
C
C             Next put the J2 separate sums back together.
C
  100 KFLAG = 0
      KPT = (J2-1)* (NDIG+1) + 1
      CALL FMEQU(MJSUMS(KPT),M03,NDIG,NDIG)
      IF (J2.GE.2) THEN
          DO 110 J = 2,J2
              CALL FMMPY(M02,M03,M03)
              KPT = (J2-J)* (NDIG+1) + 1
              CALL FMADD(M03,MJSUMS(KPT),M03)
  110     CONTINUE
      END IF
C
C             Reverse the effect of halving the argument to
C             compute EXP(MA).
C
      NDIG = NDSAV1
      IF (K2.GT.0) THEN
          CALL FMI2M(2,MB)
          DO 120 J = 1,K2
              CALL FMADD(M03,MB,M02)
              CALL FMMPY(M03,M02,M03)
  120     CONTINUE
      END IF
C
      CALL FMI2M(1,MB)
      CALL FMADD(M03,MB,MB)
      CALL FMEQU(MB,MB,NDSAV1,NDSAVE)
      NDIG = NDSAVE
C
      RETURN

      END
      SUBROUTINE FMI2M(INTEG,MA)
C
C  MA = INTEG
C
C  Convert an integer to FM format.
C
C  The conversion is exact if INTEG is less than JBASE**NDIG,
C  otherwise the result is an approximation.
C
C  This routine performs the trace printing for the conversion.
C  FMIM is used to do the arithmetic.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INTEG
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMIM,FMNTR,FMNTRI
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 16
      IF (NTRACE.NE.0) CALL FMNTRI(2,INTEG,1)
C
      CALL FMIM(INTEG,MA)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMIM(INTEG,MA)
C
C  MA = INTEG.  Internal integer conversion routine.
C
C  The conversion is exact if INTEG is less than JBASE**NDIG.
C  Otherwise FMDM is used to get an approximation.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INTEG
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER INT,J,JM2,K,KB,KB1,L,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      N1 = NDIG + 1
C
C              Check for INTEG equal to zero.
C
      INT = ABS(INTEG)
      IF (INT.EQ.0) THEN
          DO 10 J = 1,N1
              MA(J) = 0
   10     CONTINUE
          GO TO 50

      END IF
C
C             Compute and store the digits, right to left.
C
      MA(1) = 0
      J = NDIG + 1
C
   20 K = INT/JBASE
      L = INT - K*JBASE
      MA(1) = MA(1) + 1
      MA(J) = L
      IF (K.GT.0) THEN
          INT = K
          J = J - 1
          IF (J.GE.2) GO TO 20
C
C             Here INTEG cannot be expressed exactly.
C
          X = INTEG
          CALL FMDM(X,MA)
          GO TO 50

      END IF
C
C             Normalize MA.
C
      KB = N1 - J + 2
      JM2 = J - 2
      DO 30 J = 2,KB
          MA(J) = MA(J+JM2)
   30 CONTINUE
      KB1 = KB + 1
      IF (KB1.LE.N1) THEN
          DO 40 J = KB1,N1
              MA(J) = 0
   40     CONTINUE
      END IF
C
      IF (INTEG.LT.0) MA(2) = -MA(2)
C
   50 RETURN

      END
      SUBROUTINE FMIM2(INTEG,MA)
C
C  special case of fmim called with first argument zero
C  this routine is used for the two calls in FMDM to prevent
C  an apparent circular recursion (FMDM-> FMIM -> FMDM ...)
C  MA = INTEG.  Internal integer conversion routine.
C
C  The conversion is exact if INTEG is less than JBASE**NDIG.
C  Otherwise FMDM is used to get an approximation.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INTEG
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,N1
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      N1 = NDIG + 1
C
C              Check for INTEG equal to zero.
C
C     INT = ABS(INTEG)
      DO 10 J = 1,N1
          MA(J) = 0
   10 CONTINUE
      RETURN

      END
      SUBROUTINE FMINP(LINE,MA,NPT,LB)
C
C  Convert an A1 character string to floating point multiple precision
C  format.
C
C  LINE is an A1 character array of length LB to be converted
C       to FM format and returned in MA.
C  NPT is a pointer telling the routine where in the array to begin
C      the conversion.  This allows more than one number to be stored
C      in an array and converted in place.
C  LB is a pointer to the last character of the field for that number.
C
C  The input number may be in integer or any real format.
C  The convention is made that if no digits appear before 'E' then 1.0
C  is assumed.  For example 'E6' is converted as '1.0E+6'.
C  In exponential format the 'E' may also be 'D', 'Q', or 'M'.
C
C  So that FMINP will convert any output from FMOUT, LINE is tested
C  to see if the input is one of the special symbols +OVERFLOW,
C  -OVERFLOW, +UNDERFLOW, -UNDERFLOW, or UNKNOWN.
C  For user input the abbreviations OVFL, UNFL, UNKN may be used.
C
C
C
C
C
C
C             Scratch array usage during FMINP:   M01 - M05
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER LB,NPT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
      CHARACTER LINE(LB)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER J,JSTATE,K,K10PWR,KASAVE,KDFLAG,KEXP,KF1,KF2,KMN,KOF,
     +        KPOWER,KRSAVE,KSIGN,KSIGNX,KSTART,KSTOP,KTENEX,KTENF1,
     +        KTENF2,KTYPE,KUF,KUK,KVAL,KWSV,L,LARGE,N1,N2,NDSAV1,
     +        NDSAVE,NTRSAV
      CHARACTER KBLANK,KMINUS
C     ..
C     .. Local Arrays ..
      INTEGER JTRANS(8,4)
      CHARACTER KOVFL(4),KUNFL(4),KUNKN(4)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIV,FMEQU,FMI2M,FMIM,FMINP2,FMINPT,FMM2I,FMMPY,
     +         FMNTR,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,LOG10,MAX,MIN,MOD,REAL
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
C     .. Data statements ..
C
C  Simulate a finite-state automaton to scan the input line
C  and build the number.  States of the machine:
C
C  1.  Initial entry to the subroutine
C  2.  Sign of the number
C  3.  Scanning digits before a decimal point
C  4.  Decimal point
C  5.  Scanning digits after a decimal point
C  6.  E, D, Q, or M - precision indicator before the exponent
C  7.  Sign of the exponent
C  8.  Scanning exponent
C  9.  Syntax error
C
C  Character types recognized by the machine:
C
C  1.  Sign (+,-)
C  2.  Numeral (0,1,...,9)
C  3.  Decimal point (.)
C  4.  Precision indicator (E,D,Q,M)
C  5.  Illegal character for number
C
C  All blanks are ignored.  The analysis of the number proceeds as
C  follows:  If the simulated machine is in state JSTATE and a character
C  of type JTYPE is encountered the new state of the machine is given by
C  JTRANS(JSTATE,JTYPE).
C
C  In this DATA statement note the array is loaded by columns.
C
C          State   1  2  3  4  5  6  7  8
C  Type
C
      DATA JTRANS/2,9,9,9,9,7,9,9,3,3,3,5,5,8,8,8,4,4,4,9,9,9,9,9,6,6,6,
     +     6,6,9,9,9/
      DATA KBLANK/' '/,KMINUS/'-'/,KOVFL/'O','V','F','L'/,KUNFL/'U','N',
     +     'F','L'/,KUNKN/'U','N','K','N'/
C     ..
C
      NCALL = NCALL + 1
      KSTACK(NCALL) = 17
      IF (ABS(NTRACE).GE.2 .AND. NCALL.LE.LVLTRC) THEN
          WRITE (KW,FMT=9000) (LINE(L),L=NPT,LB)
      END IF

      NDSAVE = NDIG
      KASAVE = KARGSW
      KARGSW = 1
      KRSAVE = KROUND
      KROUND = 1
      KFLAG = 0
C
C             Since arithmetic tracing is not usually desired during
C             I/O conversion, disable tracing during this routine.
C
      NTRSAV = NTRACE
      NTRACE = 0
C
C             Check for special symbols.
C
      KMN = 1
      KOF = 1
      KUF = 1
      KUK = 1
      DO 10 J = NPT,LB
          IF (LINE(J).EQ.KMINUS) KMN = -1
          IF (LINE(J).EQ.KOVFL(KOF)) THEN
              KOF = KOF + 1
              IF (KOF.EQ.5) THEN
                  CALL FMIM(0,MA)
                  MA(1) = KEXPOV
                  MA(2) = KMN
                  GO TO 140

              END IF

          END IF

          IF (LINE(J).EQ.KUNFL(KUF)) THEN
              KUF = KUF + 1
              IF (KUF.EQ.5) THEN
                  CALL FMIM(0,MA)
                  MA(1) = KEXPUN
                  MA(2) = KMN
                  GO TO 140

              END IF

          END IF

          IF (LINE(J).EQ.KUNKN(KUK)) THEN
              KUK = KUK + 1
              IF (KUK.EQ.5) THEN
                  CALL FMIM(0,MA)
                  MA(1) = KUNKNO
                  MA(2) = 1
                  GO TO 140

              END IF

          END IF

   10 CONTINUE
C
C             Increase the working precision.
C
      IF (NCALL.EQ.1) THEN
          B = JBASE
          K = 5.0*LOG(10.0)/LOG(B) + 2.0
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MA(1) = KUNKNO
              MA(2) = 1
              DO 20 J = 2,NDSAVE
                  MA(J+1) = 0
   20         CONTINUE
              GO TO 140

          END IF

      END IF

      NDSAV1 = NDIG
      KSTART = NPT
      KSTOP = LB
      JSTATE = 1
      KSIGN = 1
      N1 = NDIG + 1
      DO 30 J = 1,N1
          M02(J) = 0
          M03(J) = 0
          M04(J) = 0
          M05(J) = 0
   30 CONTINUE
C
C             If JBASE is a power of ten then call FMINP2 for a much
C             faster input conversion.
C
      KPOWER = LOG10(REAL(JBASE)) + 0.5
      IF (JBASE.EQ.10**KPOWER) THEN
          CALL FMINP2(MA,LINE,KSTART,KSTOP,JTRANS,KBLANK,KPOWER)
          GO TO 130

      END IF
C
      N2 = 0
      KSIGNX = 1
      KF1 = 0
      KF2 = 0
      KEXP = 0
      KTENF1 = 1
      KTENF2 = 1
      KTENEX = 1
      K10PWR = 0
C
C             Large is a threshold used in order to do as much of the
C             conversion as possible in one-word integer arithmetic.
C
      LARGE = (MXBASE*MXBASE-10)/10
C
C             KDFLAG will be 1 if any digits are found before 'E'.
C
      KDFLAG = 0
C
C             Scan the number.
C
      DO 100 J = KSTART,KSTOP
          IF (LINE(J).EQ.KBLANK) GO TO 100
          CALL FMINPT(LINE(J),KTYPE,KVAL)
          IF (KTYPE.GE.5) GO TO 150
C
          JSTATE = JTRANS(JSTATE,KTYPE)
C
          GO TO (150,40,50,100,60,70,80,90,
     +           150) JSTATE
C
C             State 2.  Sign of the number.
C
   40     KSIGN = KVAL
          GO TO 100
C
C             State 3.  Digits before a decimal point.
C
   50     KDFLAG = 1
          KF1 = 10*KF1 + KVAL
          KTENF1 = 10*KTENF1
          IF (KTENF1.GT.LARGE) THEN
              IF (KTENF1.NE.K10PWR .AND. M03(2).NE.0) THEN
                  CALL FMI2M(KTENF1,MA)
                  K10PWR = KTENF1
              END IF

              IF (M03(2).EQ.0) THEN
                  CALL FMI2M(KF1,M03)

              ELSE
                  NDIG = MAX(2,MIN(M03(1)+MA(1),NDSAV1))
                  CALL FMMPY(M03,MA,M03)
                  NDIG = NDSAV1
                  CALL FMI2M(KF1,M02)
                  NDIG = MAX(2,MIN(MAX(M03(1),M02(1))+1,NDSAV1))
                  IF (KF1.NE.0) CALL FMADD(M03,M02,M03)
                  NDIG = NDSAV1
              END IF

              KF1 = 0
              KTENF1 = 1
          END IF

          GO TO 100
C
C             State 5.  Digits after a decimal point.
C
   60     KDFLAG = 1
          N2 = N2 + 1
          KF2 = 10*KF2 + KVAL
          KTENF2 = 10*KTENF2
          IF (KTENF2.GT.LARGE) THEN
              IF (KTENF2.NE.K10PWR .AND. M04(2).NE.0) THEN
                  CALL FMI2M(KTENF2,MA)
                  K10PWR = KTENF2
              END IF

              IF (M04(2).EQ.0) THEN
                  CALL FMI2M(KF2,M04)

              ELSE
                  NDIG = MAX(2,MIN(M04(1)+MA(1),NDSAV1))
                  CALL FMMPY(M04,MA,M04)
                  NDIG = NDSAV1
                  CALL FMI2M(KF2,M02)
                  NDIG = MAX(2,MIN(MAX(M04(1),M02(1))+1,NDSAV1))
                  IF (KF2.NE.0) CALL FMADD(M04,M02,M04)
                  NDIG = NDSAV1
              END IF

              KF2 = 0
              KTENF2 = 1
          END IF

          GO TO 100
C
C             State 6.  Precision indicator.
C
   70     IF (KDFLAG.EQ.0) CALL FMI2M(1,M03)
          GO TO 100
C
C             State 7.  Sign of the exponent.
C
   80     KSIGNX = KVAL
          GO TO 100
C
C             State 8.  Digits of the exponent.
C
   90     KEXP = 10*KEXP + KVAL
          KTENEX = 10*KTENEX
          IF (KTENEX.GT.LARGE) THEN
              IF (KTENEX.NE.K10PWR .AND. M05(2).NE.0) THEN
                  CALL FMI2M(KTENEX,MA)
                  K10PWR = KTENEX
              END IF

              IF (M05(2).EQ.0) THEN
                  CALL FMI2M(KEXP,M05)

              ELSE
                  NDIG = MAX(2,MIN(M05(1)+MA(1),NDSAV1))
                  CALL FMMPY(M05,MA,M05)
                  NDIG = NDSAV1
                  CALL FMI2M(KEXP,M02)
                  NDIG = MAX(2,MIN(MAX(M05(1),M02(1))+1,NDSAV1))
                  IF (KEXP.NE.0) CALL FMADD(M05,M02,M05)
                  NDIG = NDSAV1
              END IF

              KEXP = 0
              KTENEX = 1
          END IF
C
  100 CONTINUE
C
C             Form the number and return.
C             MA = KSIGN*(M03 + M04/10.0**N2)*10.0**M05
C
      IF (KTENF1.GT.1) THEN
          IF (KTENF1.NE.K10PWR .AND. M03(2).NE.0) THEN
              CALL FMI2M(KTENF1,MA)
              K10PWR = KTENF1
          END IF

          IF (M03(2).EQ.0) THEN
              CALL FMI2M(KF1,M03)

          ELSE
              NDIG = MAX(2,MIN(M03(1)+MA(1),NDSAV1))
              CALL FMMPY(M03,MA,M03)
              NDIG = NDSAV1
              CALL FMI2M(KF1,M02)
              NDIG = MAX(2,MIN(MAX(M03(1),M02(1))+1,NDSAV1))
              IF (KF1.NE.0) CALL FMADD(M03,M02,M03)
              NDIG = NDSAV1
          END IF

      END IF

      IF (KTENF2.GT.1) THEN
          IF (KTENF2.NE.K10PWR .AND. M04(2).NE.0) THEN
              CALL FMI2M(KTENF2,MA)
              K10PWR = KTENF2
          END IF

          IF (M04(2).EQ.0) THEN
              CALL FMI2M(KF2,M04)

          ELSE
              NDIG = MAX(2,MIN(M04(1)+MA(1),NDSAV1))
              CALL FMMPY(M04,MA,M04)
              NDIG = NDSAV1
              CALL FMI2M(KF2,M02)
              NDIG = MAX(2,MIN(MAX(M04(1),M02(1))+1,NDSAV1))
              IF (KF2.NE.0) CALL FMADD(M04,M02,M04)
              NDIG = NDSAV1
          END IF

      END IF

      IF (KTENEX.GT.1) THEN
          IF (KTENEX.NE.K10PWR .AND. M05(2).NE.0) THEN
              CALL FMI2M(KTENEX,MA)
              K10PWR = KTENEX
          END IF

          IF (M05(2).EQ.0) THEN
              CALL FMI2M(KEXP,M05)

          ELSE
              NDIG = MAX(2,MIN(M05(1)+MA(1),NDSAV1))
              CALL FMMPY(M05,MA,M05)
              NDIG = NDSAV1
              CALL FMI2M(KEXP,M02)
              NDIG = MAX(2,MIN(MAX(M05(1),M02(1))+1,NDSAV1))
              IF (KEXP.NE.0) CALL FMADD(M05,M02,M05)
              NDIG = NDSAV1
          END IF

      END IF
C
      IF (KSIGNX.EQ.-1) M05(2) = -M05(2)
      IF (M04(2).NE.0) THEN
          CALL FMI2M(10,M02)
          K = N2
          IF (MOD(K,2).EQ.0) THEN
              CALL FMI2M(1,MA)

          ELSE
              CALL FMEQU(M02,MA,NDIG,NDIG)
          END IF
C
  110     K = K/2
          NDIG = MAX(2,MIN(2*M02(1),NDSAV1))
          CALL FMMPY(M02,M02,M02)
          IF (MOD(K,2).EQ.1) THEN
              NDIG = MAX(2,MIN(M02(1)+MA(1),NDSAV1))
              CALL FMMPY(M02,MA,MA)
          END IF

          IF (K.GT.1) GO TO 110
          NDIG = NDSAV1
          CALL FMDIV(M04,MA,M04)
      END IF

      IF (M05(2).NE.0) THEN
          CALL FMI2M(10,M02)
          KWSV = KWARN
          KWARN = 0
          CALL FMM2I(M05,KEXP)
          KWARN = KWSV
          IF (KFLAG.NE.0) GO TO 150
          K = ABS(KEXP)
          IF (MOD(K,2).EQ.0) THEN
              CALL FMI2M(1,M05)

          ELSE
              CALL FMEQU(M02,M05,NDIG,NDIG)
          END IF
C
  120     K = K/2
          NDIG = MAX(2,MIN(2*M02(1),NDSAV1))
          CALL FMMPY(M02,M02,M02)
          IF (MOD(K,2).EQ.1) THEN
              NDIG = MAX(2,MIN(M02(1)+M05(1),NDSAV1))
              CALL FMMPY(M02,M05,M05)
          END IF

          IF (K.GT.1) GO TO 120
          NDIG = NDSAV1
          IF (KEXP.LT.0) THEN
              CALL FMI2M(1,M02)
              CALL FMDIV(M02,M05,M05)
          END IF

      END IF

      CALL FMADD(M03,M04,MA)
      IF (M05(2).NE.0) CALL FMMPY(MA,M05,MA)
      IF (KSIGN.EQ.-1) MA(2) = -MA(2)
  130 CALL FMEQU(MA,MA,NDIG,NDSAVE)
      IF (MA(1).EQ.KUNKNO) GO TO 150
C
  140 NDIG = NDSAVE
      KARGSW = KASAVE
      NTRACE = NTRSAV
      KROUND = KRSAVE
      IF (KFLAG.EQ.1) KFLAG = 0
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN
C
C             Error in converting the number.
C
  150 CALL FMIM(0,MA)
      MA(1) = KUNKNO
      MA(2) = 1
      KFLAG = -7
      CALL FMWARN
      GO TO 140

 9000 FORMAT (' INPUT TO FMINP',/ (14X,65A1))
      END
      SUBROUTINE FMINP2(MA,LINE,KSTART,KSTOP,JTRANS,KBLANK,KPOWER)
C
C  Internal routine for input conversion for a power of ten JBASE.
C
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KPOWER,KSTART,KSTOP
      CHARACTER KBLANK
C     ..
C     .. Array Arguments ..
      INTEGER JTRANS(8,4),MA(LUNPCK)
      CHARACTER LINE(KSTOP)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,JSTATE,KDFLAG,KEXP,KF1,KF1DIG,KF2,KF2DIG,KF2PT,KNZDIG,
     +        KSHIFT,KSIGN,KSIGNX,KTYPE,KVAL,LARGE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIVI,FMI2M,FMIM,FMINPT,FMMPY,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      JSTATE = 1
      KDFLAG = 0
      KSIGN = 1
      KSIGNX = 1
      KF1 = 0
      KNZDIG = 0
      KF1DIG = 0
      KF2 = 0
      KF2DIG = 0
      KF2PT = 2
      KEXP = 0
      LARGE = MXBASE*MXBASE/10
C
C             Scan the number.
C
      DO 70 J = KSTART,KSTOP
          IF (LINE(J).EQ.KBLANK) GO TO 70
          CALL FMINPT(LINE(J),KTYPE,KVAL)
          IF (KTYPE.GE.5) GO TO 80
C
          JSTATE = JTRANS(JSTATE,KTYPE)
C
          GO TO (80,10,20,70,30,40,50,60,
     +           80) JSTATE
C
C             State 2.  Sign of the number.
C
   10     KSIGN = KVAL
          GO TO 70
C
C             State 3.  Digits before a decimal point.
C
   20     KDFLAG = 1
          KF1 = 10*KF1 + KVAL
          IF (KVAL.GT.0 .OR. KNZDIG.NE.0) THEN
              KNZDIG = 1
              KF1DIG = KF1DIG + 1
          END IF

          IF (KF1DIG.EQ.KPOWER) THEN
              M03(1) = M03(1) + 1
              IF (M03(1).LT.NDIG) M03(M03(1)+1) = KF1
              KF1 = 0
              KF1DIG = 0
          END IF

          GO TO 70
C
C             State 5.  Digits after a decimal point.
C
   30     KDFLAG = 1
          IF (KF2PT.GT.NDIG+1) GO TO 70
          KF2 = 10*KF2 + KVAL
          KF2DIG = KF2DIG + 1
          IF (KF2DIG.EQ.KPOWER) THEN
              M04(KF2PT) = KF2
              IF (KF2.EQ.0 .AND. KF2PT.EQ.2) THEN
                  M04(1) = M04(1) - 1

              ELSE
                  KF2PT = KF2PT + 1
              END IF

              KF2 = 0
              KF2DIG = 0
          END IF

          GO TO 70
C
C             State 6.  Precision indicator.
C
   40     IF (KDFLAG.EQ.0) CALL FMI2M(1,M03)
          GO TO 70
C
C             State 7.  Sign of the exponent.
C
   50     KSIGNX = KVAL
          GO TO 70
C
C             State 8.  Digits of the exponent.
C
   60     IF (KEXP.GE.LARGE) THEN
              IF (M03(2).EQ.0 .AND. M04(2).EQ.0) THEN
                  CALL FMI2M(0,MA)
                  RETURN

              END IF

              CALL FMIM(0,MA)
              IF (KSIGNX.EQ.1) THEN
                  MA(1) = KEXPOV
                  KFLAG = -4

              ELSE
                  MA(1) = KEXPUN
                  KFLAG = -4
              END IF

              MA(2) = KSIGN
              CALL FMWARN
              RETURN

          END IF
C
          KEXP = 10*KEXP + KVAL
C
   70 CONTINUE
C
C             Form the number and return.
C             MA = KSIGN*(M03 + M04)*10.0**(KSIGNX*KEXP)
C
      IF (KF1DIG.NE.0) THEN
          M03(1) = M03(1) + 1
          KSHIFT = 10** (KPOWER-KF1DIG)
          IF (M03(1).LT.NDIG) M03(M03(1)+1) = KF1*KSHIFT
          IF (KSHIFT.GT.1) CALL FMDIVI(M03,KSHIFT,M03)
      END IF
C
      IF (KF2DIG.NE.0) THEN
          KSHIFT = 10** (KPOWER-KF2DIG)
          M04(KF2PT) = KF2*KSHIFT
      END IF

      IF (M04(2).EQ.0) M04(1) = 0
C
      IF (KEXP.NE.0) THEN
          IF (KSIGNX.EQ.1) THEN
              M05(1) = KEXP/KPOWER + 1
              M05(2) = 10** (MOD(KEXP,KPOWER))

          ELSE
              M05(1) = - ((KEXP-1)/KPOWER)
              KSHIFT = 10** (MOD(KEXP,KPOWER))
              IF (KSHIFT.GT.1) THEN
                  M05(2) = JBASE/KSHIFT

              ELSE
                  M05(2) = 1
              END IF

          END IF

      END IF
C
      CALL FMADD(M03,M04,MA)
      IF (KEXP.GT.0) CALL FMMPY(MA,M05,MA)
      MA(2) = KSIGN*MA(2)
C
      RETURN
C
C             Error in converting the number.
C
   80 CALL FMIM(0,MA)
      MA(1) = KUNKNO
      MA(2) = 1
      RETURN

      END
      SUBROUTINE FMINPT(KCHAR,KTYPE,KVAL)
C
C  Character KCHAR is looked up in the table of valid characters.
C  KTYPE is returned as its character type.
C  KVAL is returned as the associated value.
C
C  Types
C    1.  Sign (+,-)
C    2.  Numeral (0,1,...,9)
C    3.  Decimal point (.)
C    4.  Precision indicator (E,D,Q,M)
C    5.  Illegal character for number
C
C     .. Scalar Arguments ..
      INTEGER KTYPE,KVAL
      CHARACTER KCHAR
C     ..
C     .. Local Scalars ..
      INTEGER J
      CHARACTER*21 LCHARS
C     ..
C     .. Local Arrays ..
      INTEGER LTYPES(21),LVALS(21)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INDEX
C     ..
C     .. Data statements ..
      DATA LCHARS/'+-0123456789.EDQMedqm'/
      DATA LTYPES/1,1,10*2,3,8*4/
      DATA LVALS/1,-1,0,1,2,3,4,5,6,7,8,9,9*0/
C     ..
C
      J = INDEX(LCHARS,KCHAR)
      IF (J.GT.0) THEN
          KTYPE = LTYPES(J)
          KVAL = LVALS(J)

      ELSE
          KTYPE = 5
          KVAL = 0
      END IF

      RETURN

      END
      SUBROUTINE FMINT(MA,MB)
C
C  MB = INT(MA)
C
C  The integer part of MA is computed and returned in MB as a multiple
C  precision floating point number.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,KA,KB,KRESLT,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMNTR,FMRSLT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 18
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
      IF (KARGSW.NE.1) THEN
          CALL FMARGS(18,1,MA,MB,KRESLT)
          IF (KRESLT.NE.0) THEN
              CALL FMRSLT(MA,MA,MB,KRESLT)
              IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
              NCALL = NCALL - 1
              RETURN

          END IF

      END IF
C
      N1 = NDIG + 1
C
C             If MA is less than one in magnitude, return zero.
C
      IF (MA(1).LE.0) THEN
          DO 10 J = 1,N1
              MB(J) = 0
   10     CONTINUE
          GO TO 50

      END IF
C
C             If the radix point is off the right end of MA then MA is
C             already an integer.  Return MA.
C
      IF (MA(1).GE.NDIG) THEN
          DO 20 J = 1,N1
              MB(J) = MA(J)
   20     CONTINUE
          GO TO 50

      END IF
C
C             Here MA has both integer and fraction parts.  Replace
C             the digits right of the radix point by zeros.
C
      KA = MA(1) + 2
      KB = KA - 1
      DO 30 J = 1,KB
          MB(J) = MA(J)
   30 CONTINUE
C
      DO 40 J = KA,N1
          MB(J) = 0
   40 CONTINUE
C
   50 IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMIPWR(MA,INT,MB)
C
C  MB = MA ** INT
C
C  Raise an FM number to an integer power.
C  The binary multiplication method used requires an average of
C  1.5 * LOG2(INT) multiplications.  MA may be negative.
C
C
C
C
C
C
C             Scratch array usage during FMIPWR:   M01
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B,XINT
      INTEGER J,JSIGN,K,KWRNSV,KWSV,MA2,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMEQU,FMI2M,FMIM,FMMPY,FMNTR,FMNTRI,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 19
      IF (NTRACE.NE.0) THEN
          CALL FMNTR(2,MA,MA,1)
          CALL FMNTRI(2,INT,0)
      END IF
C
C             Check for special cases.
C
      IF (MA(1).EQ.KUNKNO .OR. (INT.LE.0.AND.MA(2).EQ.0)) THEN
          MA2 = MA(2)
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          IF (INT.LE.0 .AND. MA2.EQ.0) CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (INT.EQ.0) THEN
          CALL FMIM(1,MB)
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (ABS(INT).EQ.1) THEN
          KWRNSV = KWARN
          KWARN = 0
          IF (INT.EQ.1) THEN
              CALL FMEQU(MA,MB,NDIG,NDIG)

          ELSE
              CALL FMIM(1,M01)
              CALL FMDIV(M01,MA,MB)
          END IF

          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          KWARN = KWRNSV
          RETURN

      END IF
C
      IF (MA(2).EQ.0) THEN
          CALL FMIM(0,MB)
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPOV) THEN
          JSIGN = 1
          IF (MA(2).LT.0) JSIGN = -1
          CALL FMIM(0,MB)
          IF (INT.GT.0) THEN
              MB(1) = KEXPOV
              MB(2) = JSIGN**MOD(INT,2)
              KFLAG = -5

          ELSE
              MB(1) = KEXPUN
              MB(2) = JSIGN**MOD(INT,2)
              KFLAG = -6
          END IF

          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPUN) THEN
          JSIGN = 1
          IF (MA(2).LT.0) JSIGN = -1
          CALL FMIM(0,MB)
          IF (INT.GT.0) THEN
              MB(1) = KEXPUN
              MB(2) = JSIGN**MOD(INT,2)
              KFLAG = -6

          ELSE
              MB(1) = KEXPOV
              MB(2) = JSIGN**MOD(INT,2)
              KFLAG = -5
          END IF

          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Increase the working precision.
C
      NDSAVE = NDIG
      IF (NCALL.EQ.1) THEN
          B = JBASE
          XINT = ABS(INT)
          K = (5.0*LOG(10.0)+LOG(XINT))/LOG(B) + 2.0
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MB(1) = KUNKNO
              MB(2) = 1
              DO 10 J = 2,NDSAVE
                  MB(J+1) = 0
   10         CONTINUE
              NDIG = NDSAVE
              IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
              NCALL = NCALL - 1
              RETURN

          END IF

      END IF
C
C             Initialize.
C
      KWSV = KWARN
      KWARN = 0
      K = ABS(INT)
      CALL FMEQU(MA,M01,NDSAVE,NDIG)
      IF (MOD(K,2).EQ.0) THEN
          CALL FMI2M(1,MB)

      ELSE
          CALL FMEQU(M01,MB,NDIG,NDIG)
      END IF
C
C             This is the multiplication loop.
C
   20 K = K/2
      CALL FMMPY(M01,M01,M01)
      IF (MOD(K,2).EQ.1) CALL FMMPY(M01,MB,MB)
      IF (K.GT.1) GO TO 20
C
C             Invert if the exponent is negative.
C
      IF (INT.LT.0) THEN
          CALL FMI2M(1,M01)
          CALL FMDIV(M01,MB,MB)
      END IF

      KWARN = KWSV
C
C             Round the result and return.
C
      CALL FMEQU(MB,MB,NDIG,NDSAVE)
      NDIG = NDSAVE
      IF (KFLAG.LT.0) CALL FMWARN
      IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMLG10(MA,MB)
C
C  MB = LOG10(MA)
C
C
C
C
C
C
C             Scratch array usage during FMLG10:   M01 - M05
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMENTR,FMEQU,FMEXIT,FMLN,FMLNI
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(20,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
      CALL FMLN(MB,MB)
      CALL FMLNI(10,M03)
      CALL FMDIV(MB,M03,MB)
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMLN(MA,MB)
C
C  MB = LOG(MA)     (Natural logarithm)
C
C
C
C
C
C
C
C             Scratch array usage during FMLN:   M01 - M05
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION Y
      INTEGER INT,J,K2,K2EXP,KASAVE,KBOT,KM1,KOVUN,KRESLT,KSCALE,KST,
     +        KWSV,LAST,MA1,MXSAVE,N1,N3,NDSAV1,NDSAVE,NDSV
C     ..
C     .. Local Arrays ..
      INTEGER NSTACK(19)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIG,FMDIV,FMDIVI,FMDP2M,FMENTR,FMEQU,FMEXIT,
     +         FMEXP,FMI2M,FMLNI,FMM2DP,FMM2I,FMMPY,FMMPYI,FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(21,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
C             If MA is close to 1, use the Taylor series:
C                   LN(1+X) = X - X**2/2 + X**3/3 - ...
C             This is faster for small X and avoids cancellation error.
C
C             This method is faster for moderate sized NDIG, but is
C             asymptotically slower by a factor of NDIG**(2/3) than
C             using Newton and FMEXP.  For JBASE=10,000 the Taylor
C             series is faster for NDIG less than about 150 (and is
C             used only when MA is between .9999 and 1.0001).
C
      IF ((MA(1).EQ.1.AND.MA(2).EQ.1.AND.MA(3).EQ.0) .OR.
     +    (MA(1).EQ.0.AND.MA(2).EQ.JBASE-1)) THEN
          CALL FMEQU(MA,M03,NDSAVE,NDIG)
          CALL FMI2M(-1,M01)
          CALL FMADD(M03,M01,M03)
C
C             The sum will be done as two concurrent series.
C
          NDSAV1 = NDIG
          CALL FMEQU(M03,M04,NDIG,NDIG)
          CALL FMDIVI(M03,2,M05)
          CALL FMMPY(M03,M03,MB)
          CALL FMEQU(M03,M02,NDIG,NDIG)
          KBOT = 2
C
   10     KBOT = KBOT + 1
          CALL FMMPY(M02,MB,M02)
          CALL FMDIVI(M02,KBOT,M01)
          NDIG = NDSAV1
          CALL FMADD(M04,M01,M04)
          NDIG = MAX(2,NDSAV1- (M04(1)-M01(1)))
          KBOT = KBOT + 1
          CALL FMDIVI(M02,KBOT,M01)
          NDIG = NDSAV1
          CALL FMADD(M05,M01,M05)
          NDIG = MAX(2,NDSAV1- (M04(1)-M01(1)))
          IF (KFLAG.NE.1) GO TO 10
C
          NDIG = NDSAV1
          CALL FMMPY(M05,M03,M05)
          CALL FMSUB(M04,M05,MB)
          GO TO 60

      END IF
C
C             Check to see if the argument is a small integer.
C             If so use FMLNI.
C
      KM1 = 0
      MA1 = MA(1)
      CALL FMEQU(MA,M05,NDSAVE,NDIG)
      KWSV = KWARN
      KWARN = 0
      CALL FMM2I(M05,INT)
      KWARN = KWSV
      IF (KFLAG.EQ.0 .AND. INT.LT.MXBASE) THEN
          CALL FMLNI(INT,MB)
          GO TO 60

      END IF
C
C             See if the argument can be scaled to a small integer.
C
      N3 = NDIG + 3
      N1 = NDIG + 1
      DO 20 J = 2,N1
          IF (M05(N3-J).NE.0) THEN
              LAST = N3 - J - 1
              GO TO 30

          END IF

   20 CONTINUE
C
   30 KSCALE = MA1 - LAST
      M05(1) = LAST
      KWSV = KWARN
      KWARN = 0
      CALL FMM2I(M05,INT)
      KWARN = KWSV
      IF (KFLAG.EQ.0 .AND. INT.LT.MXBASE) THEN
          CALL FMLNI(INT,M04)
          IF (INT.EQ.1) KM1 = 1
          K2EXP = 0
          GO TO 50

      END IF
C
C             For the non-integer case scale the argument to lie between
C             e/2 and e to speed up the calls to FMEXP.
C
      M05(1) = 1
      KSCALE = MA1 - 1
      CALL FMM2DP(M05,Y)
      K2EXP = LOG(2.0*Y/2.71828)/0.693147
      IF (Y.LT.1.359141) THEN
          K2EXP = -1
          CALL FMMPYI(M05,2,M05)
          Y = 2.0*Y

      ELSE
          K2 = 2**K2EXP
          CALL FMDIVI(M05,K2,M05)
          Y = Y/K2
      END IF
C
C             Generate the initial approximation.
C
      Y = LOG(Y)
      CALL FMDP2M(Y,M04)
      CALL FMI2M(1,MB)
      CALL FMDIG(NSTACK,KST)
C
C             Newton iteration.
C
      DO 40 J = 1,KST
          NDIG = NSTACK(J)
          CALL FMEXP(M04,MB)
          CALL FMSUB(M05,MB,M02)
          CALL FMDIV(M02,MB,MB)
          CALL FMADD(M04,MB,M04)
   40 CONTINUE
C
C             Compute LN(JBASE**KSCALE).
C
   50 IF ((NJBLB.NE.JBASE.OR.NDIGLB.LT.NDIG) .AND. KSCALE.NE.0) THEN
          NDSV = NDIG
          NDIG = MIN(NDIG+2,MXNDG2)
          CALL FMLNI(JBASE,MLBSAV)
          NJBLB = JBASE
          NDIGLB = NDIG
          NDIG = NDSV
      END IF
C
      IF (KSCALE.NE.0 .AND. KM1.EQ.0) THEN
          CALL FMMPYI(MLBSAV,KSCALE,MB)
          CALL FMADD(M04,MB,MB)

      ELSE IF (KSCALE.NE.0 .AND. KM1.EQ.1) THEN
          CALL FMMPYI(MLBSAV,KSCALE,MB)

      ELSE IF (KSCALE.EQ.0 .AND. KM1.EQ.0) THEN
          CALL FMEQU(M04,MB,NDIG,NDIG)

      ELSE IF (KSCALE.EQ.0 .AND. KM1.EQ.1) THEN
          CALL FMI2M(0,MB)
      END IF
C
      IF (K2EXP.NE.0) THEN
          CALL FMLNI(2,M04)
          CALL FMMPYI(M04,K2EXP,M04)
          CALL FMADD(MB,M04,MB)
      END IF
C
C             Round the result and return.
C
   60 CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMLNI(INT,MA)
C
C  MA = LOG(INT)
C
C  Compute the natural logarithm of an integer INT.
C
C  If INT has only powers of 2, 3, 5, and 7 in its factorization then
C  FMLNI is faster than FMLN.  Otherwise, if INT.GE.MXBASE (i.e., INT
C  does not fit in 1/2 word) then FMLN is usually faster.
C
C  Use FMLN instead of FMLNI if 10*INT would cause integer overflow.
C
C
C
C
C
C
C
C             Scratch array usage during FMLNI:   M01 - M02
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B,XINT
      INTEGER INT2,J,J2,J3,J5,J7,JTEMP2,JTEMP3,JTEMP5,JTEMP7,K,K2,K3,K5,
     +        K7,KA,KASAVE,KB,KC,KD,KDELTA,LAST,ND,NDSAVE,NDSV,NT
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMEQU,FMI2M,FMIM,FMLNI2,FMMPYI,FMNTR,FMNTRI,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 22
      IF (NTRACE.NE.0) CALL FMNTRI(2,INT,1)
C
C             Check for special cases.
C
      IF (INT.LE.0) THEN
          CALL FMIM(0,MA)
          MA(1) = KUNKNO
          MA(2) = 1
          KFLAG = -4
          CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (INT.EQ.1) THEN
          CALL FMI2M(0,MA)
          IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Increase the working precision.
C
      NDSAVE = NDIG
      IF (NCALL.EQ.1) THEN
          B = JBASE
          K = 5.0*LOG(10.0)/LOG(B) + 2.0
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MA(1) = KUNKNO
              MA(2) = 1
              DO 10 J = 2,NDSAVE
                  MA(J+1) = 0
   10         CONTINUE
              NDIG = NDSAVE
              IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
              NCALL = NCALL - 1
              RETURN

          END IF

      END IF

      KASAVE = KARGSW
      KARGSW = 1
C
C             Find integers K2, K3, K5, and K7 such that
C                NT = 2**K2 * 3**K3 * 5**K5 * 7**K7
C             is a good approximation of INT.
C             KDELTA = ABS(INT - NT).
C
      INT2 = INT
      IF (INT.GT.MAXINT/100) INT2 = INT/100
      KDELTA = INT2
      NT = 0
      K2 = 0
      K3 = 0
      K5 = 0
      K7 = 0
C
C             Start the search loop.
C
      XINT = INT2
      LAST = LOG(XINT)/LOG(2.0) + 2.0
C
      JTEMP7 = 1
      DO 80 J7 = 1,LAST
          IF (JTEMP7.GT.INT2 .AND. ABS(JTEMP7-INT2).GT.
     +        KDELTA) GO TO 90
C
          JTEMP5 = JTEMP7
          DO 60 J5 = 1,LAST
              IF (JTEMP5.GT.INT2 .AND. ABS(JTEMP5-INT2).GT.
     +            KDELTA) GO TO 70
C
              JTEMP3 = JTEMP5
              DO 40 J3 = 1,LAST
                  IF (JTEMP3.GT.INT2 .AND. ABS(JTEMP3-INT2).GT.
     +                KDELTA) GO TO 50
C
                  JTEMP2 = JTEMP3
                  DO 20 J2 = 1,LAST
                      IF (ABS(JTEMP2-INT2).LE.KDELTA) THEN
                          IF (ABS(JTEMP2-INT2).EQ.KDELTA .AND.
     +                        JTEMP2.LT.INT2) GO TO 30
                          KDELTA = ABS(JTEMP2-INT2)
                          NT = JTEMP2
                          K2 = J2 - 1
                          K3 = J3 - 1
                          K5 = J5 - 1
                          K7 = J7 - 1
                          IF (KDELTA.EQ.0) GO TO 90
                      END IF

                      IF (JTEMP2.GT.INT2) GO TO 30
C
                      JTEMP2 = 2*JTEMP2
   20             CONTINUE
C
   30             JTEMP3 = 3*JTEMP3
   40         CONTINUE
C
   50         JTEMP5 = 5*JTEMP5
   60     CONTINUE
C
   70     JTEMP7 = 7*JTEMP7
   80 CONTINUE
C
C             If INT was too close to the integer overflow limit,
C             restore NT to an approximation of INT.
C
   90 IF (INT2.NE.INT) THEN
          IF (NT.LE.INT2) THEN
              NT = NT*100
              K2 = K2 + 2
              K5 = K5 + 2

          ELSE IF (NT.LE.INT/98) THEN
              NT = NT*98
              K2 = K2 + 1
              K7 = K7 + 2

          ELSE
              NT = NT*70
              K2 = K2 + 1
              K5 = K5 + 1
              K7 = K7 + 1
          END IF

      END IF
C
C             End of the search.  Now compute LN(NT) as a linear
C             combination of LN(125/126), LN(224/225), LN(2400/2401),
C             and LN(4374/4375).  Each of these arguments is near 1.0
C             and each logarithm is computed using only addition and
C             integer division operations -- both O(NDIG) speed.
C
      KA = -72*K2 - 114*K3 - 167*K5 - 202*K7
      KB = -27*K2 - 43*K3 - 63*K5 - 76*K7
      KC = 19*K2 + 30*K3 + 44*K5 + 53*K7
      KD = -31*K2 - 49*K3 - 72*K5 - 87*K7
C
      IF (JBASE.NE.NJBLI .OR. NDIG.GT.NDIGLI) THEN
          NDSV = NDIG
          NDIG = MIN(NDIG+2,MXNDG2)
          NJBLI = JBASE
          NDIGLI = NDIG
          CALL FMLNI2(1,126,MLN1)
          CALL FMLNI2(1,225,MLN2)
          CALL FMLNI2(1,2401,MLN3)
          CALL FMLNI2(1,4375,MLN4)
          NDIG = NDSV
      END IF
C
C             If NT.NE.INT then the final step is to compute LN(INT/NT)
C             and then use LN(INT) = LN(INT/NT) + LN(NT).
C
      IF (NT.NE.INT) THEN
          ND = NT - INT
          CALL FMLNI2(ND,NT,MA)
      END IF
C
      CALL FMMPYI(MLN1,KA,M02)
      CALL FMMPYI(MLN2,KB,M01)
      CALL FMADD(M02,M01,M02)
      CALL FMMPYI(MLN3,KC,M01)
      CALL FMADD(M02,M01,M02)
      CALL FMMPYI(MLN4,KD,M01)
      IF (NT.NE.INT) CALL FMADD(M02,MA,M02)
      CALL FMADD(M02,M01,MA)
C
C             Round and move the result to MA.
C
      CALL FMEQU(MA,MA,NDIG,NDSAVE)
      NDIG = NDSAVE
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      KARGSW = KASAVE
      RETURN

      END
      SUBROUTINE FMLNI2(INT1,INT2,MA)
C
C  MA = LN(1 - INT1/INT2)
C
C  Taylor series for computing the logarithm of a rational number
C  near 1.
C
C
C
C
C
C
C
C             Scratch array usage during FMLNI2:   M01 - M02
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT1,INT2
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIVI,FMEQU,FMI2M,FMMPYI
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMI2M(INT1,M02)
      CALL FMDIVI(M02,INT2,M02)
      CALL FMEQU(M02,MA,NDIG,NDIG)
      NDSAVE = NDIG
      J = 1
C
   10 J = J + 1
      IF (INT1.NE.1) CALL FMMPYI(M02,INT1,M02)
      CALL FMDIVI(M02,INT2,M02)
      CALL FMDIVI(M02,J,M01)
      NDIG = NDSAVE
      CALL FMADD(MA,M01,MA)
      NDIG = NDSAVE - (MA(1)-M01(1))
      IF (NDIG.LT.2) NDIG = 2
      IF (KFLAG.NE.1) GO TO 10
C
      MA(2) = -MA(2)
      NDIG = NDSAVE
      RETURN

      END
      SUBROUTINE FMM2DP(MA,X)
C
C  X = MA
C
C  Convert an FM number to double precision.
C
C  If KFLAG = -4 is returned for a value of MA which is in the range
C  of the machine's double precision number system, change the
C  definition of DPMAX in routine FMSET to reflect the current machine's
C  range.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KRESLT
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMMD,FMNTR,FMNTRR,FMWARN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 23
      CALL FMARGS(23,1,MA,MA,KRESLT)
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
      IF (KRESLT.NE.0) THEN
C
C             Here no valid result can be returned.  Set X to some
C             value which the user is likely to recognize as wrong.
C
          X = RUNKNO
          KFLAG = -4
          IF (MA(1).NE.KUNKNO) CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTRR(1,X,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      CALL FMMD(MA,X)
C
      IF (NTRACE.NE.0) CALL FMNTRR(1,X,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMM2I(MA,INT)
C
C  INT = MA
C
C  Convert an FM number to integer.
C
C  KFLAG =  0 is returned if the conversion is exact.
C        = -4 is returned if MA is larger than MAXINT in magnitude.
C             INT = IUNKNO is returned as an indication that INT
C             could not be computed without integer overflow.
C        =  2 is returned if MA is smaller than MXBASE**2 in magnitude
C             but MA is not an integer.  The next integer toward zero
C             is returned in INT.
C  It is sometimes convenient to call FMM2I in order to see if an FM
C  number can be represented as a one-word integer by checking KFLAG
C  upon return.  To avoid an unwanted error message being printed in
C  the KFLAG=-4 case, set KWARN=0 before the call to FMM2I and reset
C  it after the call.
C
C  This routine performs the trace printing for the conversion.
C  FMMI is used to do the arithmetic.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMI,FMNTR,FMNTRI
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 24
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
C
      CALL FMMI(MA,INT)
C
      IF (NTRACE.NE.0) CALL FMNTRI(1,INT,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMM2SP(MA,X)
C
C  X = MA
C
C  Convert an FM number to single precision.
C
C  MA is converted and the result is returned in X.
C
C  If KFLAG = -4 is returned for a value of MA which is in the range
C  of the machine's single precision number system, change the
C  definition of SPMAX in routine FMSET to reflect the current machine's
C  range.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION Y
      INTEGER KRESLT
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMMD,FMNTR,FMNTRR,FMWARN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 25
      CALL FMARGS(25,1,MA,MA,KRESLT)
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
      IF (KRESLT.NE.0) THEN
C
C             Here no valid result can be returned.  Set X to some
C             value which the user is likely to recognize as wrong.
C
          X = RUNKNO
          KFLAG = -4
          IF (MA(1).NE.KUNKNO) CALL FMWARN
          Y = X
          IF (NTRACE.NE.0) CALL FMNTRR(1,Y,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      CALL FMMD(MA,Y)
      X = Y
C
      Y = X
      IF (NTRACE.NE.0) CALL FMNTRR(1,Y,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMMAX(MA,MB,MC)
C
C  MC = MAX(MA,MB)
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KWSV
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMIM,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 26
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      KWSV = KWARN
      KWARN = 0
      IF (MA(1).EQ.KUNKNO .OR. MB(1).EQ.KUNKNO) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = -4

      ELSE IF (FMCOMP(MA,'LT',MB)) THEN
          CALL FMEQU(MB,MC,NDIG,NDIG)

      ELSE
          CALL FMEQU(MA,MC,NDIG,NDIG)
      END IF
C
      KWARN = KWSV
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMMD(MA,X)
C
C  X = MA
C
C  Internal routine for conversion to double precision.
C
C
C
C
C
C
C             Check to see if MA is in range for single or double
C             precision.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ONE,PMAX,RZERO,XBASE,Y,YT
      REAL DLOGDP,DLOGJB
      INTEGER J,KWSV,MA1,MA2,N1,NCASE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMI,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,LOG
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      PMAX = DPMAX
      IF (NCALL.GT.0) THEN
          IF (KSTACK(NCALL).EQ.25) PMAX = SPMAX
      END IF

      DLOGJB = LOG(DBLE(JBASE))
      DLOGDP = LOG(PMAX)
      MA1 = MA(1)
      NCASE = 0
      IF (DBLE(MA(1)-1)*DLOGJB.GT.DLOGDP) THEN
          KFLAG = -4
          X = RUNKNO
          CALL FMWARN
          RETURN

      ELSE IF (DBLE(MA(1)+1)*DLOGJB.GT.DLOGDP) THEN
          MA(1) = MA(1) - 2
          NCASE = 1

      ELSE IF (DBLE(MA(1)+1)*DLOGJB.LT.-DLOGDP) THEN
          KFLAG = -10
          X = 0.0D0
          CALL FMWARN
          RETURN

      ELSE IF (DBLE(MA(1)-1)*DLOGJB.LT.-DLOGDP) THEN
          MA(1) = MA(1) + 2
          NCASE = 2
      END IF
C
C             Try FMMI first so that small integers will be
C             converted exactly.
C
      KWSV = KWARN
      KWARN = 0
      CALL FMMI(MA,J)
      KWARN = KWSV
      IF (KFLAG.EQ.0) THEN
          X = J
          RETURN

      END IF

      KFLAG = 0
C
      MA2 = MA(2)
      MA(2) = ABS(MA2)
      RZERO = 0.0
      ONE = 1.0
      N1 = NDIG + 1
      XBASE = JBASE
      X = RZERO
      Y = ONE
      DO 10 J = 2,N1
          Y = Y/XBASE
          YT = MA(J)
          X = X + Y*YT
          YT = ONE + Y*XBASE
          IF (YT.LE.ONE) GO TO 20
   10 CONTINUE
C
   20 X = X*XBASE**MA(1)
      IF (MA2.LT.0) X = -X
      MA(2) = MA2
C
C             Check the result if it is near overflow or underflow.
C
      IF (NCASE.EQ.1) THEN
          IF (X.LE.PMAX/ (XBASE*XBASE)) THEN
              X = X*XBASE*XBASE

          ELSE
              KFLAG = -4
              X = RUNKNO
              CALL FMWARN
          END IF

      ELSE IF (NCASE.EQ.2) THEN
          IF (X.GE. (1.0D0/PMAX)*XBASE*XBASE) THEN
              X = X/ (XBASE*XBASE)

          ELSE
              KFLAG = -10
              X = 0.0D0
              CALL FMWARN
          END IF

      END IF

      MA(1) = MA1
      RETURN

      END
      SUBROUTINE FMMI(MA,INT)
C
C  INT = MA.  Internal FM to integer conversion routine.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,KA,KB,LARGE,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      N1 = NDIG + 1
      LARGE = MAXINT/JBASE
      INT = 0
      IF (MA(1).LE.0) THEN
          IF (MA(2).NE.0) KFLAG = 2
          RETURN

      END IF
C
      KB = MA(1) + 1
      INT = ABS(MA(2))
      IF (KB.GE.3) THEN
          DO 10 J = 3,KB
              IF (INT.GT.LARGE) THEN
                  KFLAG = -4
                  IF (MA(1).NE.KUNKNO) CALL FMWARN
                  INT = IUNKNO
                  RETURN

              END IF

              IF (J.LE.N1) THEN
                  INT = INT*JBASE
                  IF (INT.GT.MAXINT-MA(J)) THEN
                      KFLAG = -4
                      IF (MA(1).NE.KUNKNO) CALL FMWARN
                      INT = IUNKNO
                      RETURN

                  ELSE
                      INT = INT + MA(J)
                  END IF

              ELSE
                  INT = INT*JBASE
              END IF

   10     CONTINUE
      END IF
C
      IF (MA(2).LT.0) INT = -INT
C
C             Check to see if MA is an integer.
C
      KA = KB + 1
      IF (KA.LE.N1) THEN
          DO 20 J = KA,N1
              IF (MA(J).NE.0) THEN
                  KFLAG = 2
                  RETURN

              END IF

   20     CONTINUE
      END IF
C
      RETURN

      END
      SUBROUTINE FMMIN(MA,MB,MC)
C
C  MC = MIN(MA,MB)
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KWSV
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMIM,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 27
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      KWSV = KWARN
      KWARN = 0
      IF (MA(1).EQ.KUNKNO .OR. MB(1).EQ.KUNKNO) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = -4

      ELSE IF (FMCOMP(MA,'GT',MB)) THEN
          CALL FMEQU(MB,MC,NDIG,NDIG)

      ELSE
          CALL FMEQU(MA,MC,NDIG,NDIG)
      END IF
C
      KWARN = KWSV
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMMOD(MA,MB,MC)
C
C  MC = MA(MOD MB).
C
C
C
C
C
C
C             Scratch array usage during FMMOD:   M01 - M03
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,KASAVE,KOVUN,KRESLT,KWSV,MXSAVE,NDSAVE,NTRSAV
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIV,FMENTR,FMEQU,FMEXIT,FMINT,FMMPY,FMNTR,FMSUB,
     +         FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(28,MA,MB,2,MC,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
      KARGSW = 0
      KWSV = KWARN
      KWARN = 0
C
      IF (MB(1).GT.MA(1) .AND. MB(2).NE.0) THEN
          CALL FMEQU(MA,M01,NDIG,NDSAVE)

      ELSE
          IF (NCALL.LE.1 .AND. MA(2).NE.0) NDIG = NDIG + MA(1) - MB(1)
          IF (NDIG.GT.MXNDG2 .OR. MB(2).EQ.0) THEN
              KFLAG = -9
              IF (MA(1).EQ.KEXPOV .OR. MB(1).EQ.KEXPUN .OR.
     +            MB(2).EQ.0) KFLAG = -4
              KWARN = KWSV
              KARGSW = KASAVE
              MXEXP = MXSAVE
              CALL FMWARN
              MC(1) = KUNKNO
              MC(2) = 1
              DO 10 J = 2,NDSAVE
                  MC(J+1) = 0
   10         CONTINUE
              NDIG = NDSAVE
              IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
              NCALL = NCALL - 1
              RETURN

          END IF
C
          CALL FMEQU(MA,M02,NDSAVE,NDIG)
          CALL FMEQU(MB,M03,NDSAVE,NDIG)
          M02(2) = ABS(M02(2))
          M03(2) = ABS(M03(2))
          CALL FMDIV(M02,M03,M01)
          CALL FMINT(M01,M01)
          CALL FMMPY(M01,M03,M01)
          CALL FMSUB(M02,M01,M01)
C
C             Due to rounding, M01 may not be between 0 and MB here.
C
          NTRSAV = NTRACE
          NTRACE = 0
          IF (FMCOMP(M01,'GE',M03)) THEN
              NTRACE = NTRSAV
              CALL FMSUB(M01,M03,M01)
          END IF

          NTRACE = NTRSAV
          IF (M01(2).LT.0) CALL FMADD(M01,M03,M01)
          IF (MA(2).LT.0 .AND. M01(1).NE.KUNKNO) M01(2) = -M01(2)
      END IF
C
      IF (KFLAG.EQ.1) KFLAG = 0
      KWARN = KWSV
      CALL FMEXIT(M01,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMMOVE(MA)
C
C  Move a result from the work area (MWA) to MA.
C
C  If the result has MWA(2)=0, then it is shifted and the exponent
C  adjusted when it is moved to MA.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,N1,N2
C     ..
C     .. External Subroutines ..
      EXTERNAL FMTRAP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (MWA(2).NE.0) THEN
          N1 = NDIG + 1
          DO 10 J = 1,N1
              MA(J) = MWA(J)
   10     CONTINUE

      ELSE
          N2 = NDIG + 2
          DO 20 J = 3,N2
              MA(J-1) = MWA(J)
   20     CONTINUE
          IF (MA(2).NE.0) THEN
              MA(1) = MWA(1) - 1

          ELSE
              MA(1) = 0
          END IF

      END IF
C
      IF (ABS(MA(1)).GT.MXEXP) CALL FMTRAP(MA)
C
      RETURN

      END
      SUBROUTINE FMMPY(MA,MB,MC)
C
C  MC = MA * MB
C
C  When one of the numbers MA, MB is known to have more zero digits
C  (base JBASE) than the other, it is faster if MB is the one with
C  more zero digits.
C
C  This routine performs the trace printing for multiplication.
C  FMMPY2 is used to do the arithmetic.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMPY2,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 29
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      CALL FMMPY2(MA,MB,MC)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMMPY2(MA,MB,MC)
C
C  Internal multiplication routine.  MC = MA * MB
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER J,KRESLT,KSHIFT,MA2,MB2,N1,NGUARD,NZMA,NZMB
C     ..
C     .. External Subroutines ..
      EXTERNAL FMARGS,FMIM,FMMOVE,FMMPY3,FMRND,FMRSLT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (KARGSW.NE.1) THEN
          CALL FMARGS(29,2,MA,MB,KRESLT)
          IF (KRESLT.NE.0) THEN
              CALL FMRSLT(MA,MB,MC,KRESLT)
              RETURN

          END IF

      ELSE
          IF (MA(2).EQ.0 .OR. MB(2).EQ.0) THEN
              CALL FMIM(0,MC)
              RETURN

          END IF

      END IF

      KFLAG = 0
C
C             NGUARD is the number of guard digits used.
C
      IF (NCALL.GT.1) THEN
          NGUARD = 1

      ELSE
          B = JBASE
          NGUARD = 5.0*LOG(10.0)/LOG(B) + 2.0
          IF (NGUARD.GT.NDIG) NGUARD = NDIG
      END IF
C
C             Save the sign of MA and MB and then work only with
C             positive numbers.
C
      MA2 = MA(2)
      MB2 = MB(2)
      MA(2) = ABS(MA(2))
      MB(2) = ABS(MB(2))
C
      N1 = NDIG + 1
C
C             If there is a good chance of finding several zero digits,
C             see which number has more zero digits.
C
      NZMA = 0
      NZMB = 0
      IF (NDIG/JBASE.GE.6) THEN
          DO 10 J = 2,N1
              IF (MA(J).EQ.0) NZMA = NZMA + 1
              IF (MB(J).EQ.0) NZMB = NZMB + 1
   10     CONTINUE
      END IF
C
C             It is faster if the second argument is the one with
C             more zero digits.
C
      IF (NZMA.GT.NZMB) THEN
          CALL FMMPY3(MB,MA,NGUARD,KSHIFT)

      ELSE
          CALL FMMPY3(MA,MB,NGUARD,KSHIFT)
      END IF
C
C             The multiplication is complete.  Round the result,
C             move it to MC, and append the correct sign.
C
      MA(2) = MA2
      MB(2) = MB2
      CALL FMRND(NDIG,NGUARD,KSHIFT)
      CALL FMMOVE(MC)
      IF (MA2*MB2.LT.0) MC(2) = -MC(2)
      RETURN

      END
      SUBROUTINE FMMPY3(MA,MB,NGUARD,KSHIFT)
C
C  Internal multiplication of MA*MB.  The result is returned in MWA.
C  Both MA and MB are positive.
C
C  NGUARD is the number of guard digits which will be used.
C  KSHIFT = 1 is returned if a left shift is pending (i.e., MWA(2)=0).
C             The shift will be done in FMMOVE.  KSHIFT = 0 is returned
C             if no shift is pending.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHIFT,NGUARD
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,JA,JB,JBM1,JM1,K,KB,KI,KJ,KK,KL,KNZ,KT,KWA,L,LMAX,
     +        MAXMWA,MBJ,MBKJ,N1,N2,N3,NBNORM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN,MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      N1 = NDIG + 1
      MWA(1) = MA(1) + MB(1)
      L = NDIG + 1 + NGUARD
      MWA(L+1) = 0
C
C             The multiplication loop begins here.
C
C             NBNORM is the minimum number of digits which can be
C                    multiplied before normalization is required.
C             MAXMWA is an upper bound on the size of values in MWA
C                    divided by (JBASE-1).  It is used to determine
C                    whether to normalize before the next digit is
C                    multiplied.
C
      NBNORM = MXBASE**2/ (JBASE-1)**2
      JBM1 = JBASE - 1
      LMAX = (MXBASE*MXBASE)/JBM1 - JBM1
      IF (NBNORM.GT.1) THEN
          MBJ = MB(2)
C
C             Count the trailing zeros in MA.
C
          KNZ = N1
          N3 = NDIG + 3
          DO 10 J = 2,N1
              IF (MA(N3-J).NE.0) THEN
                  KNZ = N3 - J
                  GO TO 20

              END IF

   10     CONTINUE
C
   20     MWA(2) = 0
          N2 = NDIG + 2
          DO 30 K = N2,L
              MWA(K) = 0
   30     CONTINUE
          DO 40 K = 2,N1
              MWA(K+1) = MA(K)*MBJ
   40     CONTINUE
          MAXMWA = MBJ
          DO 70 J = 3,N1
              MBJ = MB(J)
              IF (MBJ.NE.0) THEN
                  MAXMWA = MAXMWA + MBJ
                  JM1 = J - 1
                  KL = MIN(KNZ,L-JM1)
C
C                       This is the inner loop for multiplication.
C                       KL.GE.2 since NGUARD.GE.1.
C
                  DO 50 K = 2,KL
                      MWA(JM1+K) = MWA(JM1+K) + MA(K)*MBJ
   50             CONTINUE
              END IF
C
              IF (MAXMWA.GT.LMAX) THEN
                  MAXMWA = 0
C
C                       Here normalization is only required for the
C                       range of digits currently changing in MWA.
C
                  JA = JM1 + 2
                  JB = JM1 + KL
                  KB = JB + 1
                  DO 60 K = JA,JB
                      KB = KB - 1
                      MWA(KB-1) = MWA(KB-1) + MWA(KB)/JBASE
                      MWA(KB) = MOD(MWA(KB),JBASE)
   60             CONTINUE
              END IF

   70     CONTINUE
C
C             Perform the final normalization.
C
          KB = L + 1
          DO 80 K = 3,L
              KB = KB - 1
              MWA(KB-1) = MWA(KB-1) + MWA(KB)/JBASE
              MWA(KB) = MOD(MWA(KB),JBASE)
   80     CONTINUE
C
      ELSE
C
C             If normalization must be done for each digit, combine
C             the two loops and normalize as the digits are multiplied.
C
          DO 90 J = 2,L
              MWA(J) = 0
   90     CONTINUE
          KJ = NDIG + 2
          DO 110 J = 2,N1
              KJ = KJ - 1
              MBKJ = MB(KJ)
              IF (MBKJ.EQ.0) GO TO 110
              KL = L - KJ + 1
              IF (KL.GT.N1) KL = N1
              KI = KL + 2
              KWA = KL + KJ + 1
              KK = 0
              DO 100 K = 2,KL
                  KT = MA(KI-K)*MBKJ + MWA(KWA-K) + KK
                  KK = KT/JBASE
                  MWA(KWA-K) = KT - JBASE*KK
  100         CONTINUE
              MWA(KWA-KL-1) = KK
  110     CONTINUE
C
      END IF
C
C             Set KSHIFT = 1 if a shift left is necessary.
C
      KSHIFT = 0
      IF (MWA(2).EQ.0) KSHIFT = 1
C
      RETURN

      END
      SUBROUTINE FMMPYI(MA,INT,MB)
C
C  MB = MA * INT
C
C  Multiply FM number MA by one word integer INT.
C
C  This routine is faster than FMMPY when INT*JBASE is a
C  one word integer.
C
C
C
C
C
C
C             Scratch array usage during FMMPYI:   M01
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER INTP,J,KA,KB,KC,KCARRY,KL,KM,KSHIFT,KT,MA2,N1,NGUARD
C     ..
C     .. External Subroutines ..
      EXTERNAL FMIM,FMMOVE,FMMPY2,FMNTR,FMNTRI,FMRND,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MOD,REAL
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 30
      IF (NTRACE.NE.0) THEN
          CALL FMNTR(2,MA,MA,1)
          CALL FMNTRI(2,INT,0)
      END IF

      KFLAG = 0
      N1 = NDIG + 1
C
C             Check for special cases.
C
      IF (MA(1).EQ.KUNKNO) THEN
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (INT.EQ.0) THEN
          CALL FMIM(0,MB)
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (ABS(INT).EQ.1) THEN
          DO 10 J = 1,N1
              MB(J) = MA(J)
   10     CONTINUE
          IF (MA(1).EQ.KEXPOV) KFLAG = -5
          IF (MA(1).EQ.KEXPUN) KFLAG = -6
          MB(2) = MA(2)*INT
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPOV) THEN
          MA2 = MA(2)
          CALL FMIM(0,MB)
          KFLAG = -5
          MB(1) = KEXPOV
          MB(2) = 1
          IF ((MA2.LT.0.AND.INT.GT.0) .OR.
     +        (MA2.GT.0.AND.INT.LT.0)) MB(2) = -1
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      IF (MA(1).EQ.KEXPUN) THEN
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          CALL FMWARN
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
C             Work with positive numbers.
C
      MA2 = MA(2)
      MA(2) = ABS(MA(2))
      INTP = ABS(INT)
C
C             To leave room for the normalization, shift the product
C             to the right KSHIFT places in MWA.
C
      KSHIFT = (LOG(REAL(MA(2)+1))+LOG(REAL(INTP)))/LOG(REAL(JBASE))
C
C             If INTP is too big use FMMPY.
C
      IF (KSHIFT.GT.NDIG .OR. INTP.GT.MXBASE**2/JBASE) THEN
          CALL FMIM(INT,M01)
          MA(2) = MA2
          CALL FMMPY2(MA,M01,MB)
          IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
          NCALL = NCALL - 1
          RETURN

      END IF
C
      MWA(1) = MA(1) + KSHIFT
      KA = 2 + KSHIFT
      KB = NDIG + 1 + KSHIFT
      KC = NDIG + 5
      IF (KB.LE.KC) THEN
          DO 20 J = KB,KC
              MWA(J) = 0
   20     CONTINUE
      END IF
C
      KCARRY = 0
      KL = KB + KA
      KM = KL - KSHIFT
C
C             This is the main multiplication loop.
C
      DO 30 J = KA,KB
          KT = MA(KM-J)*INTP + KCARRY
          MWA(KL-J) = MOD(KT,JBASE)
          KCARRY = KT/JBASE
   30 CONTINUE
C
C             Resolve the final carry.
C
      KB = KA - 1
      KL = KB + 2
      IF (2.LE.KB) THEN
          DO 40 J = 2,KB
              MWA(KL-J) = MOD(KCARRY,JBASE)
              KCARRY = KCARRY/JBASE
   40     CONTINUE
      END IF
C
C             Now the first significant digit in the product is in
C             MWA(2) or MWA(3).  Round the result and move it to MB.
C
      MA(2) = MA2
      IF (MWA(2).EQ.0) THEN
          NGUARD = KSHIFT - 1
          CALL FMRND(NDIG,NGUARD,1)

      ELSE
          CALL FMRND(NDIG,KSHIFT,0)
      END IF

      CALL FMMOVE(MB)
C
C             Put the sign on the result.
C
      IF ((INT.GT.0.AND.MA2.LT.0) .OR.
     +    (INT.LT.0.AND.MA2.GT.0)) MB(2) = -MB(2)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMNAME(NAME)
C
C  Convert the FM routine number in KSTACK(NCALL) to an A4 character
C  string for output.
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      CHARACTER*4 NAME
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER K
C     ..
C     .. Local Arrays ..
      CHARACTER*4 NROUT(43)
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
C     .. Data statements ..
C
      DATA NROUT/'ABS ','ACOS','ADD ','ASIN','ATAN','ATN2','BIG ',
     +     'COMP','COS ','COSH','DIM ','DIV ','DIVI','DP2M','EXP ',
     +     'I2M ','INP ','INT ','IPWR','LG10','LN  ','LNI ','M2DP',
     +     'M2I ','M2SP','MAX ','MIN ','MOD ','MPY ','MPYI','NINT',
     +     'OUT ','PI  ','PWR ','SIGN','SIN ','SINH','SP2M','SQRT',
     +     'SUB ','TAN ','TANH','ULP '/
C     ..
C
      K = KSTACK(NCALL)
      NAME = NROUT(K)
      RETURN

      END
      SUBROUTINE FMNINT(MA,MB)
C
C  MB = NINT(MA)  --  MB is returned as the nearest integer to MA.
C
C
C
C
C
C
C             Scratch array usage during FMNINT:   M01
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,KWSV,MA2,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIVI,FMENTR,FMEQU,FMEXIT,FMI2M,FMINT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(31,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      KWSV = KWARN
      KWARN = 0
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
      MA2 = MA(2)
      MB(2) = ABS(MB(2))
      CALL FMI2M(1,M01)
      CALL FMDIVI(M01,2,M01)
      CALL FMADD(MB,M01,MB)
      CALL FMINT(MB,MB)
      IF (MA2.LT.0) MB(2) = -MB(2)
      KWARN = KWSV
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMNTR(NTR,MA,MB,NARG)
C
C  Print FM numbers in base 10 format using FMOUT for conversion.
C  This is used for trace output from the FM routines.
C
C  NTR =  1 if a result of an FM call is to be printed.
C      =  2 to print input argument(s) to an FM call.
C
C  MA  -  the FM number to be printed.
C
C  MB  -  an optional second FM number to be printed.
C
C  NARG - the number of arguments.  NARG = 1 if only MA is to be
C         printed, and NARG = 2 if both MA and MB are to be printed.
C
C
C  NTRACE and LVLTRC (in COMMON /FMUSER/) control trace printout.
C
C  NTRACE = 0        No printout except warnings and errors.
C
C  NTRACE = 1        The result of each call to one of the routines
C                    is printed in base 10, using FMOUT.
C
C  NTRACE = -1       The result of each call to one of the routines
C                    is printed in internal base JBASE format.
C
C  NTRACE = 2        The input arguments and result of each call to one
C                    of the routines is printed in base 10, using FMOUT.
C
C  NTRACE = -2       The input arguments and result of each call to one
C                    of the routines is printed in base JBASE format.
C
C  LVLTRC defines the call level to which the trace is done.  LVLTRC = 1
C         means only FM routines called directly by the user are traced,
C         LVLTRC = 2 also prints traces for FM routines called by other
C         FM routines called directly by the user, etc.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER NARG,NTR
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      CHARACTER*4 NAME
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNAME,FMNTRJ,FMPRNT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (NTRACE.EQ.0) RETURN
      IF (NCALL.GT.LVLTRC) RETURN
      IF (NTR.EQ.2 .AND. ABS(NTRACE).EQ.1) RETURN
C
      IF (NTR.EQ.2) THEN
          CALL FMNAME(NAME)
          WRITE (KW,FMT=9000) NAME

      ELSE
          CALL FMNAME(NAME)
          IF (KFLAG.EQ.0) THEN
              WRITE (KW,FMT=9010) NAME,NCALL,JBASE,NDIG

          ELSE
              WRITE (KW,FMT=9020) NAME,NCALL,JBASE,NDIG,KFLAG
          END IF

      END IF
C
C             Check for base JBASE internal format trace.
C
      IF (NTRACE.LT.0) THEN
          CALL FMNTRJ(MA,NDIG)
          IF (NARG.EQ.2) CALL FMNTRJ(MB,NDIG)
      END IF
C
C             Check for base 10 trace using FMOUT.
C
      IF (NTRACE.GT.0) THEN
          CALL FMPRNT(MA)
          IF (NARG.EQ.2) CALL FMPRNT(MB)
      END IF
C
      RETURN

 9000 FORMAT (' INPUT TO FM',A4)
 9010 FORMAT (' FM',A4,15X,'CALL LEVEL =',I2,5X,'JBASE =',I10,5X,'NDIG',
     +       ' =',I6)
 9020 FORMAT (' FM',A4,6X,'CALL LEVEL =',I2,4X,'JBASE =',I10,4X,
     +       'NDIG =',I6,4X,'KFLAG =',I3)
      END
      SUBROUTINE FMNTRI(NTR,N,KNAM)
C
C  Internal routine for trace output of integer variables.
C
C  NTR = 1 for output values
C        2 for input values
C
C  N     Integer to be printed.
C
C  KNAM  is positive if the routine name is to be printed.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KNAM,N,NTR
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      CHARACTER*4 NAME
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (NTRACE.EQ.0) RETURN
      IF (NCALL.GT.LVLTRC) RETURN
      IF (NTR.EQ.2 .AND. ABS(NTRACE).EQ.1) RETURN
C
      IF (NTR.EQ.2 .AND. KNAM.GT.0) THEN
          CALL FMNAME(NAME)
          WRITE (KW,FMT=9000) NAME
      END IF

      IF (NTR.EQ.1 .AND. KNAM.GT.0) THEN
          CALL FMNAME(NAME)
          IF (KFLAG.EQ.0) THEN
              WRITE (KW,FMT=9010) NAME,NCALL,JBASE,NDIG

          ELSE
              WRITE (KW,FMT=9020) NAME,NCALL,JBASE,NDIG,KFLAG
          END IF

      END IF
C
      WRITE (KW,FMT=9030) N
C
      RETURN

 9000 FORMAT (' INPUT TO FM',A4)
 9010 FORMAT (' FM',A4,15X,'CALL LEVEL =',I2,5X,'JBASE =',I10,5X,'NDIG',
     +       ' =',I6)
 9020 FORMAT (' FM',A4,6X,'CALL LEVEL =',I2,4X,'JBASE =',I10,4X,
     +       'NDIG =',I6,4X,'KFLAG =',I3)
 9030 FORMAT (1X,I18)
      END
      SUBROUTINE FMNTRJ(MA,ND)
C
C  Print trace output in internal base JBASE format.  The number to
C  be printed is in MA.
C
C  ND is the number of base JBASE digits to be printed.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER ND
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,N1
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      N1 = ND + 1
C
      IF (JBASE.LE.10) THEN
          IF (ND.LE.25) THEN
              WRITE (KW,FMT=9000) (MA(J),J=1,N1)

          ELSE
              WRITE (KW,FMT=9010) (MA(J),J=1,N1)
          END IF

          RETURN

      END IF
C
      IF (JBASE.LE.100) THEN
          IF (ND.LE.15) THEN
              WRITE (KW,FMT=9020) (MA(J),J=1,N1)

          ELSE
              WRITE (KW,FMT=9030) (MA(J),J=1,N1)
          END IF

          RETURN

      END IF
C
      IF (JBASE.LE.1000) THEN
          IF (ND.LE.10) THEN
              WRITE (KW,FMT=9040) (MA(J),J=1,N1)

          ELSE
              WRITE (KW,FMT=9050) (MA(J),J=1,N1)
          END IF

          RETURN

      END IF
C
      IF (JBASE.LE.10000) THEN
          IF (ND.LE.10) THEN
              WRITE (KW,FMT=9060) (MA(J),J=1,N1)

          ELSE
              WRITE (KW,FMT=9070) (MA(J),J=1,N1)
          END IF

          RETURN

      END IF
C
      IF (JBASE.LE.100000) THEN
          IF (ND.LE.9) THEN
              WRITE (KW,FMT=9080) (MA(J),J=1,N1)

          ELSE
              WRITE (KW,FMT=9090) (MA(J),J=1,N1)
          END IF

          RETURN

      END IF
C
      IF (ND.LE.4) THEN
          WRITE (KW,FMT=9100) (MA(J),J=1,N1)

      ELSE
          WRITE (KW,FMT=9110) (MA(J),J=1,N1)
      END IF
C
      RETURN

 9000 FORMAT (1X,I19,I4,24I2)
 9010 FORMAT (1X,I19,I4,24I2,/ (22X,25I2))
 9020 FORMAT (1X,I19,I5,14I3)
 9030 FORMAT (1X,I19,I5,14I3,/ (22X,15I3))
 9040 FORMAT (1X,I19,I6,9I4)
 9050 FORMAT (1X,I19,I6,9I4,/ (22X,10I4))
 9060 FORMAT (1X,I19,I7,9I5)
 9070 FORMAT (1X,I19,I7,9I5,/ (22X,10I5))
 9080 FORMAT (1X,I19,I8,8I6)
 9090 FORMAT (1X,I19,I8,8I6,/ (22X,9I6))
 9100 FORMAT (1X,I19,I13,3I11)
 9110 FORMAT (1X,I19,I13,3I11,/ (22X,4I11))
      END
      SUBROUTINE FMNTRR(NTR,X,KNAM)
C
C  Internal routine for trace output of real variables.
C
C  NTR - 1 for output values
C        2 for input values
C
C  X   - Double precision value to be printed if NX.EQ.1
C
C  KNAM - Positive if the routine name is to be printed.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER KNAM,NTR
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      CHARACTER*4 NAME
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNAME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (NTRACE.EQ.0) RETURN
      IF (NCALL.GT.LVLTRC) RETURN
      IF (NTR.EQ.2 .AND. ABS(NTRACE).EQ.1) RETURN
C
      IF (NTR.EQ.2 .AND. KNAM.GT.0) THEN
          CALL FMNAME(NAME)
          WRITE (KW,FMT=9000) NAME
      END IF

      IF (NTR.EQ.1 .AND. KNAM.GT.0) THEN
          CALL FMNAME(NAME)
          IF (KFLAG.EQ.0) THEN
              WRITE (KW,FMT=9010) NAME,NCALL,JBASE,NDIG

          ELSE
              WRITE (KW,FMT=9020) NAME,NCALL,JBASE,NDIG,KFLAG
          END IF

      END IF
C
      WRITE (KW,FMT=9030) X
C
      RETURN

 9000 FORMAT (' INPUT TO FM',A4)
 9010 FORMAT (' FM',A4,15X,'CALL LEVEL =',I2,5X,'JBASE =',I10,5X,'NDIG',
     +       ' =',I6)
 9020 FORMAT (' FM',A4,6X,'CALL LEVEL =',I2,4X,'JBASE =',I10,4X,
     +       'NDIG =',I6,4X,'KFLAG =',I3)
 9030 FORMAT (1X,D30.20)
      END
      SUBROUTINE FMOUT(MA,LINE,LB)
C
C  Convert a floating multiple precision number to a character array
C  for output.
C
C  MA   is an FM number to be converted to an A1 character
C       array in base 10 format
C  LINE is the CHARACTER*1 array in which the result is returned.
C  LB   is the length of LINE.
C
C JFORM1 and JFORM2 (in COMMON) determine the format of LINE.
C
C JFORM1 = 0  normal setting  ( .314159M+6 )
C        = 1  1PE format      ( 3.14159M+5 )
C        = 2  F   format      ( 314159.000 )
C
C JFORM2 = number of significant digits to display (if JFORM1 = 0, 1)
C        = number of digits after the decimal point (if JFORM1 = 2)
C
C          If JFORM2.EQ.0 and JFORM1.NE.2 then a default number of
C          digits is chosen.  The default is roughly the full precision
C          of MA.
C
C          If JFORM2.EQ.0 and JFORM1.EQ.2 then the number is returned in
C          integer format with no decimal point.  Rounding is done as
C          with other settings, so the value displayed is the nearest
C          integer to MA.
C
C  If JFORM1.EQ.2 and MA is too large or too small to display in the
C  requested format, it is converted using JFORM1=0, JFORM2=0.
C
C  LINE should be dimensioned at least LOG10(JBASE)*NDIG + 15 on a
C  32-bit machine to allow for up to 10 digit exponents.  Replace
C  15 by 20 if 48-bit integers are used, 25 for 64-bit integers, etc.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER LB
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
      CHARACTER LINE(LB)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL X
      INTEGER J,JBSAVE,JDPT,JF1SAV,JF2SAV,K,K1,K2,KA,KASAVE,KB,KC,
     +        KDIGIT,KEXP,KEXP10,KEXPSH,KMS2SD,KMT,KPT,KRSAVE,L,MS1,MS2,
     +        MXSAVE,ND,NDE,NDE2,NDIGMS,NDS2,NDSAVE,NPOWER,NSD1,NSD2,
     +        NTRSAV,NVAL,NWORD,NWORD1,NWORD2
      CHARACTER KCHAR,KSTAR,M,MINUS,NBL,NDPT,NPLUS
C     ..
C     .. Local Arrays ..
      INTEGER MD(LUNPCK),MS(LUNPCK),MT(LUNPCK)
      CHARACTER NEXPOV(12),NEXPUN(12),NUMB(10),NUNKNO(12)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD2,FMDIV2,FMEQU,FMIM,FMMPY2,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG10,MAX,MIN,MOD,REAL
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
C     .. Data statements ..
C
      DATA NUMB/'0','1','2','3','4','5','6','7','8','9'/
      DATA NBL/' '/,NDPT/'.'/,M/'M'/,NPLUS/'+'/,MINUS/'-'/,KSTAR/'*'/
      DATA NUNKNO/' ',' ',' ','U','N','K','N','O','W','N',' ',' '/
      DATA NEXPOV/' ',' ',' ','O','V','E','R','F','L','O','W',' '/
      DATA NEXPUN/' ',' ',' ','U','N','D','E','R','F','L','O','W'/
C     ..
C
C             To avoid recursion, FMOUT calls only internal arithmetic
C             routines (FMADD2, FMMPY2, etc.), so no trace printout is
C             done during a call to FMOUT.
C
      NTRSAV = NTRACE
      NTRACE = 0
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 32
      DO 10 J = 1,LB
          LINE(J) = NBL
   10 CONTINUE
C
C             Check for special cases.
C
      IF (MA(1).EQ.KUNKNO) THEN
          DO 20 J = 1,12
              LINE(J) = NUNKNO(J)
   20     CONTINUE
          NCALL = NCALL - 1
          NTRACE = NTRSAV
          RETURN

      END IF

      IF (MA(1).EQ.KEXPOV) THEN
          DO 30 J = 1,12
              LINE(J) = NEXPOV(J)
   30     CONTINUE
          LINE(2) = NPLUS
          IF (MA(2).LT.0) LINE(2) = MINUS
          NCALL = NCALL - 1
          NTRACE = NTRSAV
          RETURN

      END IF

      IF (MA(1).EQ.KEXPUN) THEN
          DO 40 J = 1,12
              LINE(J) = NEXPUN(J)
   40     CONTINUE
          LINE(2) = NPLUS
          IF (MA(2).LT.0) LINE(2) = MINUS
          NCALL = NCALL - 1
          NTRACE = NTRSAV
          RETURN

      END IF

      KASAVE = KARGSW
      KARGSW = 1
      KRSAVE = KROUND
      KROUND = 1
      JF1SAV = JFORM1
      JF2SAV = JFORM2
C
C             ND is the number of base 10 digits required.
C
   50 ND = JFORM2
      IF (JFORM1.EQ.2 .AND. MA(1).GT.0) ND = JFORM2 +
     +    MA(1)*LOG10(REAL(JBASE)) + 1
      IF (ND.LE.1) THEN
          K = NDIG*LOG10(REAL(JBASE))
          ND = MAX(K,JFORM2)
      END IF

      IF (JFORM2.LE.0 .AND. JFORM1.LE.1) ND = 1.1 +
     +    (NDIG-1)*LOG10(REAL(JBASE))
      IF (ND.LT.2) ND = 2
C
      IF (LB.LT.ND+6) THEN
          IF (JFORM1.EQ.2) THEN
              JFORM1 = 0
              JFORM2 = 0
              GO TO 50

          END IF

          GO TO 270

      END IF
C
C             Convert to the base which is the largest power of 10
C             less than MXBASE and build the output number.
C
      NPOWER = LOG10(REAL(MXBASE))
      JBSAVE = JBASE
      NDSAVE = NDIG
      MXSAVE = MXEXP
      MXEXP = MXEXP2
      JBASE = 10**NPOWER
      NDIG = ND/NPOWER + 2
      IF (NDIG.LT.2) NDIG = 2
      IF (NDIG.GT.MXNDG2) THEN
          KFLAG = -9
          CALL FMWARN
          GO TO 270

      END IF
C
      IF (MA(2).EQ.0) THEN
          CALL FMIM(0,MS)
          GO TO 110

      END IF
C
C             Check to see if MA is already in a base which is a
C             power of ten.  If so the conversion can be skipped.
C
      K = NPOWER
      DO 60 J = 1,K
          JBASE = 10**J
          IF (JBASE.EQ.JBSAVE) THEN
              NPOWER = J
              NDIG = ND/NPOWER + 2
              IF (NDIG.LT.2) NDIG = 2
              IF (NDIG.GT.MXNDG2) THEN
                  KFLAG = -9
                  CALL FMWARN
                  GO TO 270

              END IF

              CALL FMEQU(MA,MS,NDSAVE,NDIG)
              MS(2) = ABS(MS(2))
              GO TO 110

          END IF

   60 CONTINUE
C
      CALL FMIM(JBSAVE,MD)
      NDS2 = NDSAVE + 1
      CALL FMIM(1,MT)
      KMT = 1
C
C             Convert the fraction part of MA to the new base.
C
      KPT = NDS2 + 1
      DO 70 J = 3,NDS2
          KPT = KPT - 1
          IF (MA(KPT).NE.0) GO TO 80
   70 CONTINUE
C
   80 KEXPSH = KPT - 1
      KDIGIT = ABS(MA(2))
      CALL FMIM(KDIGIT,MS)
      NDIGMS = NDIG
C
      DO 90 J = 3,KPT
          KDIGIT = MA(J)
          IF (JBSAVE.EQ.2) THEN
              NDIG = MIN(NDIGMS,MAX(2,MS(1)+1))
              CALL FMADD2(MS,MS,MS)

          ELSE
              NDIG = MIN(NDIGMS,MAX(2,MS(1)+MD(1)))
              CALL FMMPY2(MS,MD,MS)
          END IF
C
          IF (KDIGIT.GT.0) THEN
              IF (KMT.NE.KDIGIT) THEN
                  NDIG = MIN(NDIGMS,MAX(2,MD(1)))
                  CALL FMIM(KDIGIT,MT)
                  KMT = KDIGIT
              END IF

              NDIG = MIN(NDIGMS,MAX(2,MAX(MS(1),MT(1))+1))
              CALL FMADD2(MS,MT,MS)
          END IF

   90 CONTINUE
C
C             Convert the exponent.
C
      NDIG = NDIGMS
      CALL FMIM(1,MT)
      K = ABS(MA(1)-KEXPSH)
      IF (MOD(K,2).EQ.1) THEN
          CALL FMEQU(MD,MT,NDIG,NDIG)

      ELSE
          CALL FMIM(1,MT)
      END IF
C
  100 K = K/2
      NDIG = MIN(NDIGMS,MAX(2,MD(1)*2))
      IF (K.GT.0) CALL FMMPY2(MD,MD,MD)
      IF (MOD(K,2).EQ.1) THEN
          NDIG = MIN(NDIGMS,MAX(2,MT(1)+MD(1)))
          CALL FMMPY2(MT,MD,MT)
      END IF

      IF (K.GT.1) GO TO 100
C
      NDIG = NDIGMS
      IF (MA(1)-KEXPSH.LT.0) THEN
          CALL FMDIV2(MS,MT,MS)

      ELSE
          CALL FMMPY2(MS,MT,MS)
      END IF
C
C             Now MS is the value of MA converted to a
C             power of ten base.
C
C             Next convert it to a character string base 10 for output.
C
C             KEXP10 is the base 10 exponent.
C             KMS2SD is the number of base 10 significant digits
C                    in MS(2).
C
  110 MS1 = MS(1)
  120 KEXP10 = NPOWER*MS(1)
      KMS2SD = NPOWER
      K = JBASE
      DO 130 J = 1,NPOWER
          K = K/10
          IF (MS(2).LT.K .AND. MS(2).NE.0) THEN
              KEXP10 = KEXP10 - 1
              KMS2SD = KMS2SD - 1
          END IF

  130 CONTINUE
C
C             For printing using JFORM1 = 1, reduce the exponent to
C             account for the fact that the decimal point and first
C             significant digit will later be swapped.
C
      IF (JFORM1.EQ.1 .AND. MS(2).NE.0) KEXP10 = KEXP10 - 1
C
C             Find the position in the unpacked number for rounding.
C             NWORD is the word in which rounding is done, or zero if
C                   no rounding is necessary.
C                   NWORD is set to -1 if JFORM1 is 2 (F format) but no
C                   significant digits would be printed.  This case
C                   defaults to JFORM1 = 0.
C             NVAL gives the position within that word where rounding
C                  occurs.
C             NSD1 is the maximum number of base 10 S.D.'s in NWORD
C                  digits of base 10**NPOWER.
C             NSD2 is the number of base 10 S.D.'s needed to get ND
C                  base 10 digits after the decimal.
C
      NSD2 = ND
      IF (JFORM1.EQ.2) THEN
          NSD2 = JFORM2 + KEXP10
          IF (NSD2.GT.ND) NSD2 = ND
          NWORD = (NSD2-KMS2SD-1+NPOWER)/NPOWER + 2
          IF (NWORD.LT.2) NWORD = -1
          IF (NWORD.GT.NDIG) NWORD = 0
          IF (NWORD.GE.2 .AND. NSD2.LE.0) NWORD = -1

      ELSE
          NWORD = (ND-KMS2SD-1+NPOWER)/NPOWER + 2
      END IF

      NSD1 = KMS2SD + NPOWER* (NWORD-2)
      IF (NWORD.LT.2) THEN
          NVAL = 0

      ELSE
          NVAL = 10** (NSD1-NSD2)
      END IF
C
C             Now do the base 10 rounding.
C
      IF (NWORD.GE.2) THEN
          X = 0.0
          IF (NVAL.GT.1) X = MOD(MS(NWORD),NVAL)
          IF (NWORD.LT.NDIG+1) X = X + MS(NWORD+1)/REAL(JBASE)
          X = X/NVAL
          IF (X.LT.0.5) GO TO 140
          MS2 = MS(2)
          MS(NWORD) = (MS(NWORD)/NVAL)*NVAL
          MS(NWORD+1) = 0
          MS(NWORD+2) = 0
          MS(NWORD) = MS(NWORD) + NVAL
          IF (MS(NWORD).GE.JBASE) THEN
              NWORD1 = NWORD - 1
              NWORD2 = NWORD - 2
              IF (NWORD.GT.2) THEN
                  CALL FMEQU(MS,MS,NWORD1,NWORD2)

              ELSE
                  MS(1) = MS(1) + 1
                  MS(2) = MS(2)/JBASE
                  MS(3) = 0
              END IF

          END IF

          IF (MS(1).NE.MS1 .OR. MS(2).NE.MS2) GO TO 120
      END IF
C
C             Build the base 10 character string.
C
  140 IF (MA(2).LT.0) LINE(1) = MINUS
      LINE(2) = NDPT
      K = 10**KMS2SD
      L = 2
      IF (NWORD.EQ.-1) NSD2 = ND
      DO 150 J = 1,NSD2
          K = K/10
          IF (K.EQ.0) THEN
              K = JBASE/10
              L = L + 1
          END IF

          KDIGIT = MS(L)/K
          MS(L) = MOD(MS(L),K)
          LINE(J+2) = NUMB(KDIGIT+1)
  150 CONTINUE
C
      KA = NSD2 + 3
      KB = ND + 2
      IF (KB.GE.KA) THEN
          DO 160 J = KA,KB
              LINE(J) = NUMB(1)
  160     CONTINUE
      END IF
C
      LINE(ND+3) = M
      LINE(ND+4) = NPLUS
      IF (KEXP10.LT.0) LINE(ND+4) = MINUS
      IF (MA(2).EQ.0) LINE(ND+4) = NBL
C
C             Build the digits of the base 10 exponent backwards,
C             then reverse them.
C
      NDE = 1
      KEXP = ABS(KEXP10)
      DO 180 J = 1,LB
          KDIGIT = MOD(KEXP,10)
          LINE(ND+4+J) = NUMB(KDIGIT+1)
          KEXP = KEXP/10
          IF (KEXP.EQ.0) GO TO 190
C
          IF (ND+5+J.GT.LB) THEN
              DO 170 K = 1,LB
                  LINE(K) = KSTAR
  170         CONTINUE
              GO TO 210

          END IF
C
          NDE = NDE + 1
  180 CONTINUE
C
  190 NDE2 = NDE/2
      IF (NDE2.LT.1) GO TO 210
      K1 = ND + 4
      K2 = ND + 5 + NDE
      DO 200 J = 1,NDE2
          K1 = K1 + 1
          K2 = K2 - 1
          KCHAR = LINE(K1)
          LINE(K1) = LINE(K2)
          LINE(K2) = KCHAR
  200 CONTINUE
C
C             If JFORM1 is 1 put the first digit left of the decimal.
C
  210 IF (JFORM1.EQ.1) THEN
          KCHAR = LINE(2)
          LINE(2) = LINE(3)
          LINE(3) = KCHAR
      END IF
C
C             If JFORM1 is 2 put the number into fixed format.
C
      IF (JFORM1.EQ.2 .AND. JFORM2.GE.0) THEN
          IF (KEXP10.LE.-JFORM2 .OR. KEXP10+2.GT.LB) THEN
              JFORM1 = 0
              JFORM2 = 0
              JBASE = JBSAVE
              NDIG = NDSAVE
              MXEXP = MXSAVE
              DO 220 J = 1,LB
                  LINE(J) = NBL
  220         CONTINUE
              GO TO 50

          END IF

          KA = ND + 3
          DO 230 J = KA,LB
              LINE(J) = NUMB(1)
  230     CONTINUE
C
          IF (KEXP10.GT.0) THEN
              KEXP = KEXP10
              DO 240 J = 1,KEXP
                  LINE(J+1) = LINE(J+2)
  240         CONTINUE
              LINE(KEXP+2) = NDPT
          END IF
C
          IF (KEXP10.LT.0) THEN
              KEXP = -KEXP10
              KA = 3 + KEXP
              KB = LB + 1
              KC = KB - KEXP
              DO 250 J = KA,LB
                  KB = KB - 1
                  KC = KC - 1
                  LINE(KB) = LINE(KC)
                  LINE(KC) = NUMB(1)
  250         CONTINUE
          END IF
C
          JDPT = 0
          DO 260 J = 1,LB
              IF (LINE(J).EQ.NDPT) JDPT = J
              IF (JDPT.GT.0 .AND. J.GT.JDPT+JFORM2) LINE(J) = NBL
  260     CONTINUE
          IF (JFORM2.EQ.0 .AND. JDPT.GT.0) LINE(KEXP+2) = NBL
C
      END IF
C
C             Restore values and return
C
      GO TO 290
C
C             LINE is not big enough to hold the number
C             of digits specified.
C
  270 KFLAG = -8
      DO 280 J = 1,LB
          LINE(J) = KSTAR
  280 CONTINUE
      CALL FMWARN
C
  290 JBASE = JBSAVE
      NDIG = NDSAVE
      MXEXP = MXSAVE
      NCALL = NCALL - 1
      KARGSW = KASAVE
      NTRACE = NTRSAV
      KROUND = KRSAVE
      JFORM1 = JF1SAV
      JFORM2 = JF2SAV
      RETURN

      END
      SUBROUTINE FMPACK(MA,MP)
C
C  MA is packed two base NDIG digits per word and returned in MP.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MP(LPACK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,KP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KP = 2
      MP(1) = MA(1)
      MP(2) = ABS(MA(2))*JBASE + MA(3)
      IF (MA(2).LT.0) MP(2) = -MP(2)
      IF (NDIG.GE.4) THEN
          DO 10 J = 4,NDIG,2
              KP = KP + 1
              MP(KP) = MA(J)*JBASE + MA(J+1)
   10     CONTINUE
      END IF

      IF (MOD(NDIG,2).EQ.1) MP(KP+1) = MA(NDIG+1)*JBASE
      RETURN

      END
      SUBROUTINE FMPI(MA)
C
C  MA = pi
C
C
C
C
C
C
C
C             Scratch array usage during FMPI:   M01 - M04
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MESAV(LUNPCK),MLBSAV(LUNPCK),MLN1(LUNPCK),
     +        MLN2(LUNPCK),MLN3(LUNPCK),MLN4(LUNPCK),MPISAV(LUNPCK),
     +        MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL B
      INTEGER J,K,KASAVE,NDSAVE,NDSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMNTR,FMPI2,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,LOG,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 33
      IF (ABS(NTRACE).GE.2 .AND. NCALL.LE.LVLTRC) THEN
          WRITE (KW,FMT=9000)
      END IF

      KASAVE = KARGSW
      KARGSW = 1
C
C             Increase the working precision.
C
      NDSAVE = NDIG
      IF (NCALL.EQ.1) THEN
          B = JBASE
          K = 5.0*LOG(10.0)/LOG(B) + 2.0
          NDIG = MAX(NDIG+K,2)
          IF (NDIG.GT.MXNDG2) THEN
              KFLAG = -9
              CALL FMWARN
              MA(1) = KUNKNO
              MA(2) = 1
              DO 10 J = 2,NDSAVE
                  MA(J+1) = 0
   10         CONTINUE
              GO TO 20

          END IF

      END IF
C
C             Check to see if pi has previously been computed
C             in base JBASE with sufficient precision.
C
      IF (NJBPI.EQ.JBASE .AND. NDIGPI.GE.NDIG) THEN
          CALL FMEQU(MPISAV,MA,NDIGPI,NDSAVE)

      ELSE
          NDSV = NDIG
          NDIG = MIN(NDIG+2,MXNDG2)
          CALL FMPI2(MPISAV)
          NJBPI = JBASE
          NDIGPI = NDIG
          CALL FMEQU(MPISAV,MA,NDIG,NDSAVE)
          NDIG = NDSV
      END IF
C
   20 NDIG = NDSAVE
      KARGSW = KASAVE
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN

 9000 FORMAT (' INPUT TO FMPI')
      END
      SUBROUTINE FMPI2(MPI)
C
C  Internal routine to compute pi.
C  The formula used is due to S. Ramanujan:
C                                                (4n)!(1103+26390n)
C  1/pi = (sqrt(8)/9801) * sum(n=0 to infinity) --------------------
C                                               ((n!)**4)(396**(4n))
C  The result is returned in MPI.
C
C
C
C
C
C
C             Scratch array usage during FMPI2:   M01 - M04
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MPI(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER J,K,KST,LARGE,MX,N,NDSAVE
C     ..
C     .. Local Arrays ..
      INTEGER NSTACK(19)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIG,FMDIV,FMDIVI,FMDP2M,FMI2M,FMMPY,FMMPYI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NDSAVE = NDIG
      N = -1
      CALL FMI2M(1103,MPI)
      CALL FMI2M(1,M02)
      CALL FMI2M(26390,M03)
      CALL FMI2M(1103,M04)
      MX = MXBASE**2/JBASE
C
   10 N = N + 1
      LARGE = MX/ (4*N+3)
      J = 4*N + 1
      IF (J.GT.LARGE) THEN
          CALL FMMPYI(M02,J,M02)
          J = J + 1
          CALL FMMPYI(M02,J,M02)
          J = J + 1
          CALL FMMPYI(M02,J,M02)

      ELSE IF (J* (J+1).GT.LARGE) THEN
          K = J* (J+1)
          CALL FMMPYI(M02,K,M02)
          J = J + 2
          CALL FMMPYI(M02,J,M02)

      ELSE
          K = J* (J+1)* (J+2)
          CALL FMMPYI(M02,K,M02)
      END IF
C
      J = N + 1
      LARGE = MXBASE/J
      IF (J.GT.LARGE) THEN
          CALL FMDIVI(M02,J,M02)
          CALL FMDIVI(M02,J,M02)
          CALL FMDIVI(M02,J,M02)

      ELSE IF (J*J.GT.LARGE) THEN
          K = J*J
          CALL FMDIVI(M02,K,M02)
          CALL FMDIVI(M02,J,M02)

      ELSE
          K = J*J*J
          CALL FMDIVI(M02,K,M02)
      END IF
C
C             Break 4/396**4 into 1/(2178*2178*1296).
C
      J = 2178
      LARGE = MXBASE/J
      IF (J.GT.LARGE) THEN
          CALL FMDIVI(M02,J,M02)
          CALL FMDIVI(M02,J,M02)
          CALL FMDIVI(M02,1296,M02)

      ELSE
          K = J*J
          CALL FMDIVI(M02,K,M02)
          CALL FMDIVI(M02,1296,M02)
      END IF
C
      CALL FMADD(M03,M04,M04)
      CALL FMMPY(M02,M04,M01)
C
      NDIG = NDSAVE
      CALL FMADD(MPI,M01,MPI)
      NDIG = MAX(2,NDSAVE- (MPI(1)-M01(1)))
      IF (KFLAG.NE.1) GO TO 10
      NDIG = NDSAVE
C
      CALL FMI2M(8,M02)
      X = 8.0
      X = SQRT(X)
      CALL FMDP2M(X,M04)
      CALL FMDIG(NSTACK,KST)
      DO 20 J = 1,KST
          NDIG = NSTACK(J)
          CALL FMDIV(M02,M04,M01)
          CALL FMADD(M04,M01,M04)
          CALL FMDIVI(M04,2,M04)
   20 CONTINUE
      CALL FMI2M(9801,M03)
      CALL FMMPY(MPI,M04,MPI)
      CALL FMDIV(M03,MPI,MPI)
C
      RETURN

      END
      SUBROUTINE FMPRNT(MA)
C
C  Print MA in base 10 format.
C
C  FMPRNT can be called directly by the user for easy output
C  in M format.  MA is converted using FMOUT and printed.
C
C
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,NBITS,LUNPCK,LMWA,LMBUFF
      PARAMETER (MXNDIG=256,NBITS=32,LUNPCK= (6*MXNDIG)/5+20,
     +          LMWA=2*LUNPCK,LMBUFF= ((LUNPCK+3)* (NBITS-1)*301)/2000+
     +          6)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
      CHARACTER MBUFF(LMBUFF)
C     ..
C     .. Local Scalars ..
      INTEGER J,K,KSAVE,L,LAST,LB,ND,NEXP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMOUT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG10,MAX,MIN,REAL
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMBUFF/MBUFF
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KSAVE = KFLAG
      ND = NDIG*LOG10(REAL(JBASE)) + 1
      IF (ND.LT.2) ND = 2
      NEXP = 2.0*LOG10(REAL(MXBASE)) + 6
      LB = MAX(JFORM2+NEXP,ND+NEXP)
      LB = MIN(LB,LMBUFF)
      CALL FMOUT(MA,MBUFF,LB)
      KFLAG = KSAVE
      LAST = LB + 1
      DO 10 J = 1,LB
          IF (MBUFF(LAST-J).NE.' ' .OR. J.EQ.LB) THEN
              L = LAST - J
              WRITE (KW,FMT=9000) (MBUFF(K),K=1,L)
              RETURN

          END IF

   10 CONTINUE
      RETURN

 9000 FORMAT (6X,73A1)
      END
      SUBROUTINE FMPWR(MA,MB,MC)
C
C  MC = MA ** MB
C
C  If MB can be expressed exactly as a one word integer, then FMIPWR is
C  used.  This is much faster when MB is small, and using FMIPWR allows
C  MA to be negative.
C
C
C
C
C
C
C             Scratch array usage during FMPWR:   M01 - M06
C
C
C             Convert MB to an integer before changing NDIG.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER IEXTRA,INTMB,J,KASAVE,KFL,KOVUN,KRESLT,KWSV,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMENTR,FMEQU,FMEXIT,FMEXP,FMIM,FMIPWR,FMLN,FMMI,FMMPY,
     +         FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KWSV = KWARN
      KWARN = 0
      CALL FMMI(MB,INTMB)
      KWARN = KWSV
      KFL = KFLAG
C
      CALL FMENTR(34,MA,MB,2,MC,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
C             If the exponent is large, raise the precision.
C
      IEXTRA = MAX(0,MB(1))
      NDIG = NDIG + IEXTRA
      IF (NDIG.GT.MXNDG2) THEN
          KFLAG = -9
          CALL FMWARN
          MC(1) = KUNKNO
          MC(2) = 1
          DO 10 J = 2,NDSAVE
              MC(J+1) = 0
   10     CONTINUE
          CALL FMEXIT(MC,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
          RETURN

      END IF
C
C             If the exponent is a small integer, call FMIPWR.
C
      KWSV = KWARN
      KWARN = 0
      CALL FMEQU(MA,M06,NDSAVE,NDIG)
      IF (KFL.EQ.0) THEN
          CALL FMIPWR(M06,INTMB,MC)

      ELSE IF (M06(2).LE.0) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = -4

      ELSE
          CALL FMLN(M06,M06)
          CALL FMEQU(MB,M02,NDSAVE,NDIG)
          CALL FMMPY(M06,M02,M06)
          CALL FMEXP(M06,MC)
      END IF

      KWARN = KWSV
C
C             Round the result and return.
C
      CALL FMEXIT(MC,MC,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMRDC(MA,MB,JSIN,JCOS,JSWAP)
C
C  Reduce MA using various trigonometric identities to an equivalent
C  angle MB between 0 and 45 degrees.  The reduction is done in radians
C  if KRAD (in common /FMUSER/) is 1, in degrees if KRAD is 0.
C  JSIN and JCOS are returned +1 or -1 and JSWAP is returned to indicate
C  that the sin and cos functions have been interchanged as follows:
C
C  JSWAP = 0 means   SIN(MA) = JSIN*SIN(MB)
C                    COS(MA) = JCOS*COS(MB)
C
C  JSWAP = 1 means   SIN(MA) = JSIN*COS(MB)
C                    COS(MA) = JCOS*SIN(MB)
C
C
C
C
C
C
C
C             Scratch array usage during FMRDC:   M01 - M04
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER JCOS,JSIN,JSWAP
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X
      INTEGER J,NDSAVE,NDSV
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMDIVI,FMEQU,FMI2M,FMINT,FMM2DP,FMMPY,FMMPYI,FMPI2,
     +         FMSUB,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      JSIN = 1
      JCOS = 1
      JSWAP = 0
      NDSAVE = NDIG
      NDIG = NDIG + MAX(0,MA(1))
C
C             If the argument is too big, return UNKNOWN.
C
      IF (NDIG.GT.MXNDG2) THEN
          KFLAG = -9
          CALL FMWARN
          MB(1) = KUNKNO
          MB(2) = 1
          DO 10 J = 2,NDSAVE
              MB(J+1) = 0
   10     CONTINUE
          NDIG = NDSAVE
          RETURN

      END IF
C
C             If MA is less than 1/JBASE no reduction is needed.
C
      IF (MA(1).LT.0) THEN
          NDIG = NDSAVE
          CALL FMEQU(MA,MB,NDIG,NDIG)
          IF (MB(2).LT.0) THEN
              MB(2) = -MB(2)
              JSIN = -1
          END IF

          RETURN

      END IF
C
      J = 1
      IF (KRAD.EQ.1) THEN
   20     IF (NJBPI.NE.JBASE .OR. NDIGPI.LT.NDIG) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMPI2(MPISAV)
              NJBPI = JBASE
              NDIGPI = NDIG
              NDIG = NDSV
          END IF

          CALL FMEQU(MA,M04,NDSAVE,NDIG)
          IF (MA(2).LT.0) JSIN = -1
          M04(2) = ABS(M04(2))
          IF (M04(1).EQ.0) THEN
              CALL FMM2DP(M04,X)
              IF (X.LE.0.75) THEN
                  NDIG = NDSAVE
                  CALL FMEQU(M04,MB,NDIG,NDIG)
                  RETURN

              END IF

          END IF

          CALL FMMPYI(MPISAV,2,M02)
          IF (FMCOMP(M04,'GE',M02)) THEN
              CALL FMDIV(M04,M02,M01)
              CALL FMINT(M01,M01)
              CALL FMMPY(M01,M02,M01)
              CALL FMSUB(M04,M01,M04)
          END IF

          CALL FMEQU(MPISAV,M03,NDIG,NDIG)
          IF (FMCOMP(M04,'GE',M03)) THEN
              JSIN = -JSIN
              CALL FMSUB(M02,M04,M04)
          END IF

          CALL FMDIVI(M02,4,M02)
          IF (FMCOMP(M04,'GE',M02)) THEN
              JCOS = -JCOS
              CALL FMSUB(M03,M04,M04)
          END IF

          CALL FMDIVI(M03,4,M03)
          IF (FMCOMP(M04,'GE',M03)) THEN
              JSWAP = 1
              CALL FMSUB(M02,M04,M04)
          END IF
C
C             If the reduced argument is close to zero, then
C             cancellation has produced an inaccurate value.
C             Raise NDIG and do the reduction again.
C
          IF (J.EQ.1 .AND. M04(1).LT.0) THEN
              J = 2
              NDIG = NDIG - M04(1)
              IF (NDIG.GT.MXNDG2) THEN
                  KFLAG = -9
                  CALL FMWARN
                  MB(1) = KUNKNO
                  MB(2) = 1
                  DO 30 J = 2,NDSAVE
                      MB(J+1) = 0
   30             CONTINUE
                  NDIG = NDSAVE
                  RETURN

              END IF

              JSIN = 1
              JCOS = 1
              JSWAP = 0
              GO TO 20

          END IF
C
      ELSE
C
          CALL FMEQU(MA,M04,NDSAVE,NDIG)
          IF (MA(2).LT.0) JSIN = -1
          M04(2) = ABS(M04(2))
          IF (M04(1).EQ.0) THEN
              CALL FMM2DP(M04,X)
              IF (X.LE.44.0) THEN
                  NDIG = NDSAVE
                  CALL FMEQU(M04,MB,NDIG,NDIG)
                  RETURN

              END IF

          END IF

          CALL FMI2M(360,M02)
          IF (FMCOMP(M04,'GE',M02)) THEN
              CALL FMDIV(M04,M02,M01)
              CALL FMINT(M01,M01)
              CALL FMMPY(M01,M02,M01)
              CALL FMSUB(M04,M01,M04)
          END IF

          CALL FMI2M(180,M03)
          IF (FMCOMP(M04,'GE',M03)) THEN
              JSIN = -JSIN
              CALL FMSUB(M02,M04,M04)
          END IF

          CALL FMI2M(90,M02)
          IF (FMCOMP(M04,'GE',M02)) THEN
              JCOS = -JCOS
              CALL FMSUB(M03,M04,M04)
          END IF

          CALL FMI2M(45,M03)
          IF (FMCOMP(M04,'GE',M03)) THEN
              JSWAP = 1
              CALL FMSUB(M02,M04,M04)
          END IF
C
      END IF
C
C             Round the result and return.
C
      CALL FMEQU(M04,MB,NDIG,NDSAVE)
      NDIG = NDSAVE
      RETURN

      END
      SUBROUTINE FMRND(ND,NGUARD,KSHIFT)
C
C  Round MWA to ND digits (base JBASE).
C
C  MWA is non-negative and has ND+NGUARD+KSHIFT digits.
C
C  NGUARD is the number of guard digits carried.
C  KSHIFT is 1 if a left shift is pending when MWA(2)=0.
C
C  Round to position MWA(ND+1+KSHIFT) using the guard digits
C  MWA(ND+2+KSHIFT), ..., MWA(ND+1+NGUARD+KSHIFT).
C
C  This routine is designed to be called only from within the FM
C  package.  The user should call FMEQU to round numbers.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHIFT,ND,NGUARD
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,K,KB,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (KROUND.EQ.0 .AND. NCALL.LE.1) RETURN
      L = ND + 2 + KSHIFT
      IF (2* (MWA(L)+1).LT.JBASE) RETURN
      IF (2*MWA(L).GT.JBASE) THEN
          MWA(L-1) = MWA(L-1) + 1
          MWA(L) = 0
          IF (MWA(L-1).LT.JBASE) RETURN
          GO TO 40

      END IF
C
C             If the first guard digit gives a value close to 1/2 then
C             further guard digits must be examined.
C
      IF (MOD(JBASE,2).EQ.0) THEN
          IF (2*MWA(L).LT.JBASE) RETURN
          IF (2*MWA(L).EQ.JBASE) THEN
              IF (NGUARD.GE.2) THEN
                  DO 10 J = 2,NGUARD
                      IF (MWA(L+J-1).GT.0) GO TO 30
   10             CONTINUE
              END IF
C
C                       Round to even.
C
              IF (MOD(MWA(L-1),2).EQ.0) RETURN
          END IF

      ELSE
          IF (2*MWA(L)+1.EQ.JBASE) THEN
              IF (NGUARD.GE.2) THEN
                  DO 20 J = 2,NGUARD
                      IF (2* (MWA(L+J-1)+1).LT.JBASE) RETURN
                      IF (2*MWA(L+J-1).GT.JBASE) GO TO 30
   20             CONTINUE
                  RETURN

              END IF

          END IF

      END IF
C
   30 MWA(L-1) = MWA(L-1) + 1
      MWA(L) = 0
C
C             Check whether there was a carry in the rounded digit.
C
   40 KB = L - 1
      IF (KB.GE.3) THEN
          K = KB + 1
          DO 50 J = 3,KB
              K = K - 1
              IF (MWA(K).LT.JBASE) RETURN
              MWA(K-1) = MWA(K-1) + MWA(K)/JBASE
              MWA(K) = MOD(MWA(K),JBASE)
   50     CONTINUE
      END IF
C
C             If there is a carry in the first digit then the exponent
C             must be adjusted and the number shifted right.
C
      IF (MWA(2).GE.JBASE) THEN
          IF (KB.GE.4) THEN
              K = KB + 1
              DO 60 J = 4,KB
                  K = K - 1
                  MWA(K) = MWA(K-1)
   60         CONTINUE
          END IF
C
          IF (KB.GE.3) MWA(3) = MOD(MWA(2),JBASE)
          MWA(2) = MWA(2)/JBASE
          MWA(1) = MWA(1) + 1
      END IF
C
      IF (KSHIFT.EQ.1 .AND. MWA(2).NE.0) KSHIFT = 0
C
      RETURN

      END
      SUBROUTINE FMRSLT(MA,MB,MC,KRESLT)
C
C  Handle results which are special cases, such as overflow,
C  underflow, unknown, etc.
C
C  MA and MB are the input arguments to an FM subroutine
C
C  MC is the result which is returned
C
C  KRESLT is the result code from FMARGS.  Result codes handled here:
C
C   0 - Perform the normal operation
C   1 - The result is the first input argument
C   2 - The result is the second input argument
C   3 - The result is -OVERFLOW
C   4 - The result is +OVERFLOW
C   5 - The result is -UNDERFLOW
C   6 - The result is +UNDERFLOW
C   7 - The result is -1.0
C   8 - The result is +1.0
C  11 - The result is 0.0
C  12 - The result is 'UNKNOWN'
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      INTEGER KRESLT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KFSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMIM
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFSAVE = KFLAG
      IF (KRESLT.EQ.1) THEN
          CALL FMEQU(MA,MC,NDIG,NDIG)
          IF (KSTACK(NCALL).EQ.3 .OR. KSTACK(NCALL).EQ.40) THEN
              KFLAG = 1

          ELSE
              KFLAG = KFSAVE
          END IF

          RETURN

      END IF
C
      IF (KRESLT.EQ.2) THEN
          CALL FMEQU(MB,MC,NDIG,NDIG)
          IF (KSTACK(NCALL).EQ.3 .OR. KSTACK(NCALL).EQ.40) THEN
              KFLAG = 1

          ELSE
              KFLAG = KFSAVE
          END IF

          RETURN

      END IF
C
      IF (KRESLT.EQ.3 .OR. KRESLT.EQ.4) THEN
          CALL FMIM(0,MC)
          MC(1) = KEXPOV
          MC(2) = 1
          IF (KRESLT.EQ.3) MC(2) = -1
          KFLAG = KFSAVE
          RETURN

      END IF
C
      IF (KRESLT.EQ.5 .OR. KRESLT.EQ.6) THEN
          CALL FMIM(0,MC)
          MC(1) = KEXPUN
          MC(2) = 1
          IF (KRESLT.EQ.5) MC(2) = -1
          KFLAG = KFSAVE
          RETURN

      END IF
C
      IF (KRESLT.EQ.7) THEN
          CALL FMIM(-1,MC)
          KFLAG = KFSAVE
          RETURN

      END IF
C
      IF (KRESLT.EQ.8) THEN
          CALL FMIM(1,MC)
          KFLAG = KFSAVE
          RETURN

      END IF
C
      IF (KRESLT.EQ.11) THEN
          CALL FMIM(0,MC)
          KFLAG = KFSAVE
          RETURN

      END IF
C
      IF (KRESLT.EQ.12 .OR. KRESLT.LT.0 .OR. KRESLT.GT.15) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = KFSAVE
          RETURN

      END IF
C
      RETURN

      END
      SUBROUTINE FMSIGN(MA,MB,MC)
C
C  MC = SIGN(MA,MB)
C
C  MC is set to ABS(MA) if MB is positive or zero,
C     or -ABS(MA) if MB is negative.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KWSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMIM,FMNTR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 35
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      KWSV = KWARN
      KWARN = 0
      IF (MA(1).EQ.KUNKNO .OR. MB(1).EQ.KUNKNO) THEN
          CALL FMIM(0,MC)
          MC(1) = KUNKNO
          MC(2) = 1
          KFLAG = -4

      ELSE IF (MB(2).GE.0) THEN
          CALL FMEQU(MA,MC,NDIG,NDIG)
          MC(2) = ABS(MC(2))

      ELSE
          CALL FMEQU(MA,MC,NDIG,NDIG)
          MC(2) = -ABS(MC(2))
      END IF
C
      KWARN = KWSV
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMSIN(MA,MB)
C
C  MB = SIN(MA)
C
C
C
C
C
C
C
C             Scratch array usage during FMSIN:   M01 - M04
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER JCOS,JSIN,JSWAP,KASAVE,KOVUN,KRESLT,MXSAVE,NDSAVE,NDSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMCOS2,FMDIVI,FMENTR,FMEQU,FMEXIT,FMMPY,FMPI2,FMRDC,
     +         FMSIN2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(36,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Reduce the argument, convert to radians if the input is
C             in degrees, and evaluate the function.
C
      CALL FMRDC(MB,MB,JSIN,JCOS,JSWAP)
      IF (MB(1).EQ.KUNKNO) GO TO 10
      IF (KRAD.EQ.0) THEN
          IF (NJBPI.NE.JBASE .OR. NDIGPI.LT.NDIG) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMPI2(MPISAV)
              NJBPI = JBASE
              NDIGPI = NDIG
              NDIG = NDSV
          END IF

          CALL FMMPY(MB,MPISAV,MB)
          CALL FMDIVI(MB,180,MB)
      END IF

      IF (MB(1).NE.KUNKNO) THEN
          IF (JSWAP.EQ.0) THEN
              CALL FMSIN2(MB,MB)

          ELSE
              CALL FMCOS2(MB,MB)
          END IF

      END IF
C
C             Append the sign, round, and return.
C
      IF (JSIN.EQ.-1 .AND. MB(1).NE.KUNKNO) MB(2) = -MB(2)
   10 CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMSIN2(MA,MB)
C
C  Internal subroutine for MB = SIN(MA) where 0.LE.MA.LE.1.
C
C
C
C
C
C
C             Scratch array usage during FMSIN2:   M01 - M03
C
C
C
C             LJSUMS = 8*LUNPCK allows for up to eight concurrent sums.
C             Increasing this value will begin to improve the speed of
C             SIN when the base is large and precision exceeds about
C             1,500 decimal digits.
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA,LJSUMS
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK,
     +          LJSUMS=8*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MJSUMS(LJSUMS),
     +        MWA(LMWA)
C     ..
C     .. Local Scalars ..
      REAL ALOG3,ALOGB,ALOGT,B,T
      INTEGER J,J2,K,K2,KPT,KTHREE,KWSV,L,L2,MAXVAL,N2,NBOT,NDSAV1,
     +        NDSAVE,NTERM
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIVI,FMEQU,FMI2M,FMIPWR,FMMPY,FMMPYI,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,MAX,MIN,REAL,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSUMS/MJSUMS
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (MA(2).EQ.0) THEN
          CALL FMI2M(0,MA)
          RETURN

      END IF

      NDSAVE = NDIG
      KWSV = KWARN
      KWARN = 0
C
C             Use the direct series
C                  SIN(X) = X - X**3/3! + X**5/5! - ...
C
C             The argument will be divided by 3**K2 before the series
C             is summed.  The series will be added as J2 concurrent
C             series.  The approximately optimal values of K2 and J2
C             are now computed to try to minimize the time required.
C             N2/2 is the approximate number of terms of the series
C             which will be needed, and L2 guard digits will be carried.
C
      B = JBASE
      K = 5.0*LOG(10.0)/LOG(B) + 2.0
      T = MAX(NDIG-K,2)
      ALOGB = LOG(B)
      ALOG3 = LOG(3.0)
      ALOGT = LOG(T)
      J2 = 0.081*ALOGB*T**0.3333 + 1.85
      J2 = MAX(1,MIN(J2,LJSUMS/MXNDG2))
      K2 = 0.1*SQRT(T*ALOGB/REAL(J2)) - 0.09*ALOGT + 0.5
C
      L = - (REAL(MA(1))*ALOGB+LOG(REAL(MA(2))/B+REAL(MA(3))/ (B*B)))/
     +    ALOG3 - 0.3
      K2 = K2 - L
      IF (K2.LT.0) THEN
          K2 = 0
          J2 = .43*SQRT(T*ALOGB/ (ALOGT+REAL(L)*ALOG3)) + .33
      END IF

      IF (J2.LE.1) J2 = 1
C
      N2 = T*ALOGB/ (ALOGT+REAL(L)*ALOG3)
      L2 = LOG(REAL(N2)+3.0**K2)/ALOGB
      NDIG = NDIG + L2
      IF (NDIG.GT.MXNDG2) THEN
          KFLAG = -9
          CALL FMWARN
          MB(1) = KUNKNO
          MB(2) = 1
          DO 10 J = 2,NDSAVE
              MB(J+1) = 0
   10     CONTINUE
          NDIG = NDSAVE
          KWARN = KWSV
          RETURN

      END IF

      NDSAV1 = NDIG
C
C             Divide the argument by 3**K2.
C
      CALL FMEQU(MA,M02,NDSAVE,NDIG)
      KTHREE = 1
      MAXVAL = MXBASE/3
      IF (K2.GT.0) THEN
          DO 20 J = 1,K2
              KTHREE = 3*KTHREE
              IF (KTHREE.GT.MAXVAL) THEN
                  CALL FMDIVI(M02,KTHREE,M02)
                  KTHREE = 1
              END IF

   20     CONTINUE
          IF (KTHREE.GT.1) CALL FMDIVI(M02,KTHREE,M02)
      END IF
C
C             Split into J2 concurrent sums and reduce NDIG while
C             computing each term in the sum as the terms get smaller.
C
      CALL FMEQU(M02,M03,NDIG,NDIG)
      NTERM = 1
      DO 30 J = 1,J2
          NBOT = NTERM* (NTERM-1)
          IF (NBOT.GT.1) CALL FMDIVI(M03,NBOT,M03)
          NTERM = NTERM + 2
          KPT = (J-1)* (NDIG+1) + 1
          CALL FMEQU(M03,MJSUMS(KPT),NDIG,NDIG)
          M03(2) = -M03(2)
   30 CONTINUE
      CALL FMMPY(M02,M02,M02)
      CALL FMIPWR(M02,J2,MB)
C
   40 CALL FMMPY(M03,MB,M03)
      DO 50 J = 1,J2
          NBOT = NTERM* (NTERM-1)
          IF (NTERM.GT.MXBASE .OR. NBOT.GT.MXBASE) THEN
              CALL FMDIVI(M03,NTERM,M03)
              NBOT = NTERM - 1
              CALL FMDIVI(M03,NBOT,M03)

          ELSE
              CALL FMDIVI(M03,NBOT,M03)
          END IF

          KPT = (J-1)* (NDSAV1+1) + 1
          NDIG = NDSAV1
          CALL FMADD(MJSUMS(KPT),M03,MJSUMS(KPT))
          IF (KFLAG.EQ.1) GO TO 60
          NDIG = NDSAV1 - (MJSUMS(KPT)-M03(1))
          IF (NDIG.LT.2) NDIG = 2
          M03(2) = -M03(2)
          NTERM = NTERM + 2
   50 CONTINUE
      GO TO 40
C
C             Next put the J2 separate sums back together.
C
   60 KFLAG = 0
      KPT = (J2-1)* (NDIG+1) + 1
      CALL FMEQU(MJSUMS(KPT),MB,NDIG,NDIG)
      IF (J2.GE.2) THEN
          DO 70 J = 2,J2
              CALL FMMPY(M02,MB,MB)
              KPT = (J2-J)* (NDIG+1) + 1
              CALL FMADD(MB,MJSUMS(KPT),MB)
   70     CONTINUE
      END IF
C
C             Reverse the effect of reducing the argument to
C             compute SIN(MA).
C
      NDIG = NDSAV1
      IF (K2.GT.0) THEN
          CALL FMI2M(3,M02)
          DO 80 J = 1,K2
              CALL FMMPY(MB,MB,M03)
              CALL FMMPYI(M03,-4,M03)
              CALL FMADD(M02,M03,M03)
              CALL FMMPY(M03,MB,MB)
   80     CONTINUE
      END IF
C
      CALL FMEQU(MB,MB,NDSAV1,NDSAVE)
      NDIG = NDSAVE
      KWARN = KWSV
C
      RETURN

      END
      SUBROUTINE FMSINH(MA,MB)
C
C  MB = SINH(MA)
C
C
C
C
C
C
C             Scratch array usage during FMSINH:   M01 - M03
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,KWRNSV,MA1,MA2,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIV,FMDIVI,FMENTR,FMEQU,FMEXIT,FMEXP,FMI2M,FMMPY,
     +         FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(37,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      MA1 = MA(1)
      MA2 = MA(2)
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Work with a positive argument.
C
      MB(2) = ABS(MB(2))
      KARGSW = 0
      KWRNSV = KWARN
      KWARN = 0
C
      CALL FMEXP(MB,MB)
      IF (MB(1).EQ.KEXPOV) THEN
          GO TO 10

      END IF

      IF (MB(1).EQ.KEXPUN) THEN
          MB(1) = KEXPOV
          MB(2) = -1
          GO TO 10

      END IF

      CALL FMI2M(1,M01)
C
C             If MA is small use the accurate value of EXP(MA) - 1
C             which is computed by FMEXP and is stored in M03 to
C             avoid cancellation error.
C
      IF (MA1.GT.0) THEN
          CALL FMDIV(M01,MB,M01)
          CALL FMSUB(MB,M01,MB)

      ELSE
          CALL FMADD(MB,M01,M01)
          CALL FMMPY(M03,M01,M01)
          CALL FMDIV(M01,MB,MB)
      END IF

      CALL FMDIVI(MB,2,MB)
C
C             Round and return.
C
   10 IF (MA2.LT.0) MB(2) = -MB(2)
      KWARN = KWRNSV
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMSP2M(X,MA)
C
C  MA = X
C
C  Convert a single precision number to FM format.
C
C  In general the relative accuracy of the number returned is only
C  the relative accuracy of a machine precision number.  This may be
C  true even if X can be represented exactly in the machine floating
C  point number system.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION Y
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDM,FMNTR,FMNTRR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 38
      Y = X
      IF (NTRACE.NE.0) CALL FMNTRR(2,Y,1)
C
      CALL FMDM(Y,MA)
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MA,MA,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMSQRT(MA,MB)
C
C  MB = SQRT(MA)
C
C
C
C
C
C
C             Scratch array usage during FMSQRT:   M01 - M02
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION X,XB
      INTEGER J,KASAVE,KE,KOVUN,KRESLT,KST,MA1,MXSAVE,NDSAVE
C     ..
C     .. Local Arrays ..
      INTEGER NSTACK(19)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIG,FMDIV,FMDIVI,FMDP2M,FMENTR,FMEQU,FMEXIT,
     +         FMM2DP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD,SQRT
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(39,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      MA1 = MA(1)
      CALL FMEQU(MA,M02,NDSAVE,NDIG)
C
C             Generate the first approximation.
C
      M02(1) = 0
      CALL FMM2DP(M02,X)
      X = SQRT(X)
      KE = MA1/2
      IF (MOD(ABS(MA1),2).EQ.1) THEN
          XB = JBASE
          X = X*SQRT(XB)
          KE = (MA1-1)/2
      END IF

      CALL FMDP2M(X,MB)
      MB(1) = MB(1) + KE
C
C             Initialize.
C
      M02(1) = MA1
      CALL FMDIG(NSTACK,KST)
C
C             Newton iteration.
C
      DO 10 J = 1,KST
          NDIG = NSTACK(J)
          CALL FMDIV(M02,MB,M01)
          CALL FMADD(MB,M01,MB)
          CALL FMDIVI(MB,2,MB)
   10 CONTINUE
C
C             Round the result and return.
C
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMSUB(MA,MB,MC)
C
C  MC = MA - MB
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK),MC(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER K1,K2,KFLG1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD2,FMNTR
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NCALL = NCALL + 1
      KSTACK(NCALL) = 40
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MB,2)
C
      KFLG1 = 0
      IF (MB(1).GT.MA(1) .OR. MA(2).EQ.0) KFLG1 = 1
      IF (MB(2).EQ.0) KFLG1 = 0
      K1 = MB(2)
      MB(2) = -MB(2)
C
      CALL FMADD2(MA,MB,MC)
C
C             If MB and MC are distinct arrays then the sign of MB must
C             be restored.  If they refer to identical memory locations
C             then MB now contains the result of the subtraction, so it
C             should not be restored.  Detecting which of these two
C             cases is true is tricky because some optimizing compilers
C             remove lines 1 and 3 of the following:
C                   K2 = MC(2)
C                   MB(2) = K1
C                   MC(2) = K2
C             The code used below tries to do the same thing in an
C             obscure way which defeats the optimizer.
C             If this still fails for some compiler, copy MB to a local
C             array and negate that as the argument to FMADD2.  This
C             adds one unpacked array to the program size and makes
C             subtraction run noticeably slower than addition.
C             See the Portability Notes before FMSET for a test to
C             detect this failure.
C
      IF (MB(2).NE.MC(2)) THEN
          MB(2) = K1

      ELSE
          K2 = (MB(2)+MC(2))/2
          MB(2) = K1
          MC(2) = K2
      END IF
C
C             If MA was smaller than MB then KFLAG = 1 returned from
C             FMADD means the result from FMSUB is the opposite of the
C             input argument of larger magnitude, so reset KFLAG.
C
      IF (KFLAG.EQ.1 .AND. KFLG1.EQ.1) KFLAG = 0
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MC,MC,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMTAN(MA,MB)
C
C  MB = TAN(MA)
C
C
C
C
C
C
C
C             Scratch array usage during FMTAN:   M01 - M04
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NDIGE,NDIGLB,NDIGLI,NDIGPI,NJBE,
     +        NJBLB,NJBLI,NJBPI,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MESAV(LUNPCK),
     +        MLBSAV(LUNPCK),MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     +        MLN4(LUNPCK),MPISAV(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER JCOS,JSIN,JSWAP,KASAVE,KOVUN,KRESLT,MXSAVE,NDSAVE,NDSV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMDIVI,FMENTR,FMEQU,FMEXIT,FMI2M,FMIM,FMMPY,FMPI2,
     +         FMRDC,FMSIN2,FMSQRT,FMSUB,FMWARN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMSAVE/NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     +       MPISAV,MESAV,MLBSAV,MLN1,MLN2,MLN3,MLN4
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(41,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Reduce the argument, convert to radians if the input is
C             in degrees, and evaluate the function.
C
      CALL FMRDC(MB,MB,JSIN,JCOS,JSWAP)
      IF (MB(1).EQ.KUNKNO) GO TO 10
      IF (MB(2).EQ.0) THEN
          IF (JSWAP.EQ.1) THEN
              CALL FMIM(0,MB)
              MB(1) = KUNKNO
              MB(2) = 1
              KFLAG = -4
              CALL FMWARN
          END IF

          GO TO 10

      END IF

      IF (KRAD.EQ.0) THEN
          IF (NJBPI.NE.JBASE .OR. NDIGPI.LT.NDIG) THEN
              NDSV = NDIG
              NDIG = MIN(NDIG+2,MXNDG2)
              CALL FMPI2(MPISAV)
              NJBPI = JBASE
              NDIGPI = NDIG
              NDIG = NDSV
          END IF

          CALL FMMPY(MB,MPISAV,MB)
          CALL FMDIVI(MB,180,MB)
      END IF

      IF (MB(1).NE.KUNKNO) THEN
          CALL FMSIN2(MB,MB)
          CALL FMMPY(MB,MB,M03)
          CALL FMI2M(1,M02)
          CALL FMSUB(M02,M03,M03)
          CALL FMSQRT(M03,M03)
          IF (JSWAP.EQ.0) THEN
              CALL FMDIV(MB,M03,MB)

          ELSE
              CALL FMDIV(M03,MB,MB)
          END IF

      END IF
C
C             Append the sign, round, and return.
C
      IF (JSIN*JCOS.EQ.-1) MB(2) = -MB(2)
   10 CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMTANH(MA,MB)
C
C  MB = TANH(MA)
C
C
C
C
C
C
C             Scratch array usage during FMTANH:   M01 - M03
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),M01(LUNPCK),M02(LUNPCK),M03(LUNPCK),
     +        M04(LUNPCK),M05(LUNPCK),M06(LUNPCK),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER KASAVE,KOVUN,KRESLT,KWRNSV,MA1,MA2,MXSAVE,NDSAVE
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMDIV,FMENTR,FMEQU,FMEXIT,FMEXP,FMI2M,FMMPY,FMSUB
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FM1/M01,M02,M03,M04,M05,M06
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      CALL FMENTR(42,MA,MA,1,MB,KRESLT,NDSAVE,MXSAVE,KASAVE,KOVUN)
      IF (KRESLT.NE.0) RETURN
C
      KARGSW = 0
      KWRNSV = KWARN
      KWARN = 0
      MA1 = MA(1)
      MA2 = MA(2)
      CALL FMEQU(MA,MB,NDSAVE,NDIG)
C
C             Work with a positive argument.
C
      MB(2) = ABS(MB(2))
C
      CALL FMEXP(MB,MB)
      IF (MB(1).EQ.KEXPOV) THEN
          CALL FMI2M(1,MB)
          GO TO 10

      END IF

      IF (MB(1).EQ.KEXPUN) THEN
          CALL FMI2M(-1,MB)
          GO TO 10

      END IF

      CALL FMI2M(1,M01)
C
C             If MA is small use the accurate value of EXP(MA) - 1
C             which is computed by FMEXP and is stored in M03 to
C             avoid cancellation error.
C
      IF (MA1.GT.0) THEN
          CALL FMDIV(M01,MB,M01)
          CALL FMSUB(MB,M01,M02)
          CALL FMADD(MB,M01,M01)
          CALL FMDIV(M02,M01,MB)

      ELSE
          CALL FMADD(MB,M01,M01)
          CALL FMMPY(M03,M01,M01)
          CALL FMI2M(2,M02)
          CALL FMADD(M01,M02,MB)
          CALL FMDIV(M01,MB,MB)
      END IF
C
C             Round and return.
C
   10 IF (MA2.LT.0) MB(2) = -MB(2)
      KWARN = KWRNSV
      CALL FMEXIT(MB,MB,NDSAVE,MXSAVE,KASAVE,KOVUN)
      RETURN

      END
      SUBROUTINE FMTRAP(MA)
C
C  If MA has overflowed or underflowed, replace it by the appropriate
C  symbol.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMWARN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (NCALL.LE.0) RETURN
      IF (MA(1).GT.MXEXP+1) THEN
          MA(1) = KEXPOV
          IF (MA(2).GT.0) THEN
              MA(2) = 1

          ELSE
              MA(2) = -1
          END IF

          KFLAG = -5
          CALL FMWARN
      END IF

      IF (MA(1).LT.-MXEXP) THEN
          MA(1) = KEXPUN
          IF (MA(2).GT.0) THEN
              MA(2) = 1

          ELSE
              MA(2) = -1
          END IF

          KFLAG = -6
          CALL FMWARN
      END IF
C
      RETURN

      END
      SUBROUTINE FMULP(MA,MB)
C
C  MB = The value of one Unit in the Last Place of MA at the current
C       base and precision.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MB(LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER J,KWRNSV,MA1,N1
C     ..
C     .. External Subroutines ..
      EXTERNAL FMIM,FMMOVE,FMNTR,FMWARN
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KFLAG = 0
      NCALL = NCALL + 1
      KSTACK(NCALL) = 43
      IF (NTRACE.NE.0) CALL FMNTR(2,MA,MA,1)
C
      MA1 = MA(1)
      N1 = NDIG + 1
      DO 10 J = 3,N1
          MWA(J) = 0
   10 CONTINUE
      MWA(2) = 1
      IF (MA(2).LT.0) MWA(2) = -1
      MWA(1) = MA(1) - NDIG + 1
      IF (MA(2).EQ.0 .OR. MA(1).GE.KEXPOV) THEN
          CALL FMIM(0,MB)
          MB(1) = KUNKNO
          MB(2) = 1
          KFLAG = -4
          IF (MA1.NE.KUNKNO) CALL FMWARN

      ELSE
          KWRNSV = KWARN
          IF (MA1.EQ.KEXPUN) KWARN = 0
          CALL FMMOVE(MB)
          KWARN = KWRNSV
      END IF
C
      IF (NTRACE.NE.0) CALL FMNTR(1,MB,MB,1)
      NCALL = NCALL - 1
      RETURN

      END
      SUBROUTINE FMUNPK(MP,MA)
C
C  MP is unpacked and the value returned in MA.
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LUNPCK),MP(LPACK)
C     ..
C     .. Scalars in Common ..
      INTEGER JBASE,JFORM1,JFORM2,KFLAG,KRAD,KROUND,KW,KWARN,LVLTRC,
     +        NDIG,NTRACE
C     ..
C     .. Local Scalars ..
      INTEGER J,KP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MOD
C     ..
C     .. Common blocks ..
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      KP = 2
      MA(1) = MP(1)
      MA(2) = ABS(MP(2))/JBASE
      MA(3) = MOD(ABS(MP(2)),JBASE)
      IF (MP(2).LT.0) MA(2) = -MA(2)
      IF (NDIG.GE.4) THEN
          DO 10 J = 4,NDIG,2
              KP = KP + 1
              MA(J) = MP(KP)/JBASE
              MA(J+1) = MOD(MP(KP),JBASE)
   10     CONTINUE
      END IF

      IF (MOD(NDIG,2).EQ.1) MA(NDIG+1) = MP(KP+1)/JBASE
      RETURN

      END
      SUBROUTINE FMWARN
C
C  Called by one of the FM routines to print a warning message
C  if any error condition arises in that routine.
C
C
C
C
C
C
C     .. Parameters ..
      INTEGER MXNDIG,LUNPCK,LMWA
      PARAMETER (MXNDIG=256,LUNPCK= (6*MXNDIG)/5+20,LMWA=2*LUNPCK)
C     ..
C     .. Scalars in Common ..
      DOUBLE PRECISION DPMAX
      REAL RUNKNO,SPMAX
      INTEGER IUNKNO,JBASE,JFORM1,JFORM2,KARGSW,KEXPOV,KEXPUN,KFLAG,
     +        KRAD,KROUND,KUNKNO,KW,KWARN,LVLTRC,MAXINT,MXBASE,MXEXP,
     +        MXEXP2,MXNDG2,NCALL,NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER KSTACK(19),MWA(LMWA)
C     ..
C     .. Local Scalars ..
      INTEGER NCS
      CHARACTER*4 NAME
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNAME
C     ..
C     .. Common blocks ..
      COMMON /FM/MWA,NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,KUNKNO,
     +       IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK,MAXINT,SPMAX,DPMAX
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      IF (KFLAG.GE.0 .OR. NCALL.NE.1 .OR. KWARN.LE.0) RETURN
      NCS = NCALL
      CALL FMNAME(NAME)
      WRITE (KW,FMT=9000) KFLAG,NAME
C
   10 NCALL = NCALL - 1
      IF (NCALL.GT.0) THEN
          CALL FMNAME(NAME)
          WRITE (KW,FMT=9010) NAME
          GO TO 10

      END IF
C
      IF (KFLAG.EQ.-1) THEN
          WRITE (KW,FMT=9020) MXNDIG

      ELSE IF (KFLAG.EQ.-2) THEN
          WRITE (KW,FMT=9030) MXBASE

      ELSE IF (KFLAG.EQ.-3) THEN
          WRITE (KW,FMT=9040)
          WRITE (KW,FMT=9050)

      ELSE IF (KFLAG.EQ.-4 .OR. KFLAG.EQ.-7) THEN
          WRITE (KW,FMT=9060)
          WRITE (KW,FMT=9050)

      ELSE IF (KFLAG.EQ.-5) THEN
          WRITE (KW,FMT=9070)

      ELSE IF (KFLAG.EQ.-6) THEN
          WRITE (KW,FMT=9080)

      ELSE IF (KFLAG.EQ.-8) THEN
          WRITE (KW,FMT=9090)

      ELSE IF (KFLAG.EQ.-9) THEN
          WRITE (KW,FMT=9100)
          WRITE (KW,FMT=9110) NDIG,MXNDG2
          WRITE (KW,FMT=9050)

      ELSE IF (KFLAG.EQ.-10) THEN
          IF (KSTACK(NCS).EQ.25) THEN
              WRITE (KW,FMT=9120)

          ELSE
              WRITE (KW,FMT=9130)
          END IF

          WRITE (KW,FMT=9140)
      END IF
C
      NCALL = NCS
      IF (KWARN.GE.2) STOP
      RETURN

 9000 FORMAT (/' ERROR OF TYPE KFLAG =',I3,' IN FM PACKAGE IN ROUTINE ',
     +       'FM',A4,/)
 9010 FORMAT (' CALLED FROM FM',A4)
 9020 FORMAT (' NDIG MUST BE BETWEEN 2 AND',I10,/)
 9030 FORMAT (' JBASE MUST BE BETWEEN 2 AND',I10,/)
 9040 FORMAT (' AN INPUT ARGUMENT IS NOT A VALID FM NUMBER.','  ITS EX',
     +       'PONENT IS OUT OF RANGE.',/)
 9050 FORMAT (' UNKNOWN HAS BEEN RETURNED.',/)
 9060 FORMAT (' INVALID INPUT ARGUMENT FOR THIS ROUTINE.',/)
 9070 FORMAT (' THE RESULT HAS OVERFLOWED.',/)
 9080 FORMAT (' THE RESULT HAS UNDERFLOWED.',/)
 9090 FORMAT (' THE RESULT ARRAY IS NOT BIG ENOUGH TO HOLD THE',' OUTP',
     +       'UT CHARACTER STRING',/' IN THE CURRENT FORMAT.',/' THE R',
     +       'ESULT ''***...***'' HAS BEEN RETURNED.',/)
 9100 FORMAT (' PRECISION COULD NOT BE RAISED ENOUGH TO PROVIDE ALL',
     +       ' REQUESTED GUARD DIGITS.',/)
 9110 FORMAT (' REQUESTED NDIG=',I7,'.  MAXIMUM AVAILABLE NDIG=',I7,/)
 9120 FORMAT (' AN FM NUMBER WAS TOO SMALL IN MAGNITUDE TO ','CONVERT ',
     +       'TO SINGLE PRECISION.',/)
 9130 FORMAT (' AN FM NUMBER WAS TOO SMALL IN MAGNITUDE TO ','CONVERT ',
     +       'TO DOUBLE PRECISION.',/)
 9140 FORMAT (' ZERO HAS BEEN RETURNED.',/)
      END
C
C  Here are the routines which work with packed FM numbers.  All names
C  are the same as unpacked versions with 'FM' replaced by 'FP'.
C
C  To convert a program using the FM package from unpacked calls to
C  packed calls make these changes to the program:
C  '(LUNPCK)' to '(LPACK)' in dimensions.
C  'CALL FM' to 'CALL FP'
C  'FMCOMP' to 'FPCOMP'.
C
      SUBROUTINE FPABS(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMABS,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMABS(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPACOS(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMACOS,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMACOS(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPADD(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMADD,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMADD(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPASIN(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMASIN,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMASIN(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPATAN(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMATAN,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMATAN(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPATN2(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMATN2,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMATN2(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPBIG(MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMBIG,FMPACK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMBIG(MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      LOGICAL FUNCTION FPCOMP(MA,LREL,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      CHARACTER*2 LREL
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Functions ..
      LOGICAL FMCOMP
      EXTERNAL FMCOMP
C     ..
C     .. External Subroutines ..
      EXTERNAL FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      FPCOMP = FMCOMP(MX,LREL,MY)
      RETURN

      END
      SUBROUTINE FPCOS(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMCOS,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMCOS(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPCOSH(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMCOSH,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMCOSH(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPDIG(NSTACK,KST)
C     .. Scalar Arguments ..
      INTEGER KST
C     ..
C     .. Array Arguments ..
      INTEGER NSTACK(19)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIG
C     ..
      CALL FMDIG(NSTACK,KST)
      RETURN

      END
      SUBROUTINE FPDIM(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIM,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMDIM(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPDIV(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIV,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMDIV(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPDIVI(MA,INT,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDIVI,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMDIVI(MX,INT,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPDP2M(X,MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMDP2M,FMPACK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMDP2M(X,MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPEQU(MA,MB,NDA,NDB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER NDA,NDB
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Scalars in Common ..
      INTEGER JBASE,JFORM1,JFORM2,KFLAG,KRAD,KROUND,KW,KWARN,LVLTRC,
     +        NDIG,NTRACE
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. Local Scalars ..
      INTEGER NDASAV,NDBSAV,NDGSAV
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEQU,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
      COMMON /FMUSER/NDIG,JBASE,JFORM1,JFORM2,KRAD,KW,NTRACE,LVLTRC,
     +       KFLAG,KWARN,KROUND
C     ..
      NDGSAV = NDIG
      NDASAV = NDA
      NDBSAV = NDB
      NDIG = NDASAV
      CALL FMUNPK(MA,MX)
      CALL FMEQU(MX,MY,NDASAV,NDBSAV)
      NDIG = NDBSAV
      CALL FMPACK(MY,MB)
      NDA = NDASAV
      NDB = NDBSAV
      NDIG = NDGSAV
      RETURN

      END
      SUBROUTINE FPEXP(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMEXP,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMEXP(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPI2M(INT,MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMI2M,FMPACK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMI2M(INT,MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPINP(LINE,MA,NPT,LB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER LB,NPT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
      CHARACTER LINE(LB)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMINP,FMPACK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMINP(LINE,MX,NPT,LB)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPINT(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMINT,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMINT(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPIPWR(MA,INT,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMIPWR,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMIPWR(MX,INT,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPLG10(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMLG10,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMLG10(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPLN(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMLN,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMLN(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPLNI(INT,MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMLNI,FMPACK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMLNI(INT,MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPM2DP(MA,X)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMM2DP,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMM2DP(MX,X)
      RETURN

      END
      SUBROUTINE FPM2I(MA,INT)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMM2I,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMM2I(MX,INT)
      RETURN

      END
      SUBROUTINE FPM2SP(MA,X)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMM2SP,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMM2SP(MX,X)
      RETURN

      END
      SUBROUTINE FPMAX(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMAX,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMMAX(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPMIN(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMIN,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMMIN(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPMOD(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMOD,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMMOD(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPMPY(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMPY,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMMPY(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPMPYI(MA,INT,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER INT
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMMPYI,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMMPYI(MX,INT,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPNINT(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMNINT,FMPACK,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMNINT(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPOUT(MA,LINE,LB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      INTEGER LB
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
      CHARACTER LINE(LB)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMOUT,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMOUT(MX,LINE,LB)
      RETURN

      END
      SUBROUTINE FPPI(MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMPI
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMPI(MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPPRNT(MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPRNT,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMPRNT(MX)
      RETURN

      END
      SUBROUTINE FPPWR(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMPWR,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMPWR(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPSET(NPREC)
C     .. Scalar Arguments ..
      INTEGER NPREC
C     ..
C     .. External Subroutines ..
      EXTERNAL FMSET
C     ..
      CALL FMSET(NPREC)
      RETURN

      END
      SUBROUTINE FPSIGN(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSIGN,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMSIGN(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPSIN(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSIN,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMSIN(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPSINH(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSINH,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMSINH(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPSP2M(X,MA)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSP2M
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMSP2M(X,MX)
      CALL FMPACK(MX,MA)
      RETURN

      END
      SUBROUTINE FPSQRT(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSQRT,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMSQRT(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPSUB(MA,MB,MC)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK),MC(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMSUB,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMUNPK(MB,MY)
      CALL FMSUB(MX,MY,MX)
      CALL FMPACK(MX,MC)
      RETURN

      END
      SUBROUTINE FPTAN(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMTAN,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMTAN(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPTANH(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMTANH,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMTANH(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN

      END
      SUBROUTINE FPULP(MA,MB)
C     .. Parameters ..
      INTEGER MXNDIG,LPACK,LUNPCK
      PARAMETER (MXNDIG=256,LPACK= (MXNDIG+1)/2+1,
     +          LUNPCK= (6*MXNDIG)/5+20)
C     ..
C     .. Array Arguments ..
      INTEGER MA(LPACK),MB(LPACK)
C     ..
C     .. Arrays in Common ..
      INTEGER MX(LUNPCK),MY(LUNPCK)
C     ..
C     .. External Subroutines ..
      EXTERNAL FMPACK,FMULP,FMUNPK
C     ..
C     .. Common blocks ..
      COMMON /FMPCK/MX,MY
C     ..
      CALL FMUNPK(MA,MX)
      CALL FMULP(MX,MX)
      CALL FMPACK(MX,MB)
      RETURN
C             This is the end of the FM package.  This is line 9352.
      END
