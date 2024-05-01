# Algorithm 693: a FORTRAN package for floating-point multiple-precision arithmetic

> [FREE ACCESS]
> ACM Transactions on Mathematical Software,
> Volume 17, Issue 2, pp 273–283, https://doi.org/10.1145/108556.108585


## FM package

```for
C     FM 1.0                  David M Smith                  1-06-90
C
C
C  The FM package performs floating point multiple precision arithmetic.
C
C  Before calling any routine in the package, several variables in the
C  common blocks /FMUSER/, /FM/, and /FMSAVE/ must be initialized.
C  These three common blocks contain information which must be saved
C  between calls, so they should be declared in the main program.
C
C  Subroutine FMSET initializes these variables to default values and
C  defines all machine-dependent values in the package.  After calling
C  FMSET once at the start of a program, the user may sometimes want
C  to reset some of the variables in common block /FMUSER/.
C
C
C  JBASE is the base in which the arithmetic is done.
C        JBASE must be bigger than one, and less than or equal
C        to the square root of the largest representable integer.
C        For best efficiency JBASE should be about 1/4 to 1/2 as big as
C        the square root of the largest representable integer.
C        Input and output conversion are much faster when JBASE is a
C        power of ten.
C
C  NDIG  is the number of base JBASE digits which are carried in the
C        multiple precision numbers.  NDIG must be at least two.
C        The upper limit for NDIG is defined in the PARAMETER statement
C        at the top of each routine and is restricted only by the amount
C        of memory available.
C
C
C  There are two representations for a floating multiple precision
C  number.  The unpacked representation used by the routines while
C  doing the computations is base JBASE and is stored in NDIG+1 words.
C  A packed representation is available to store the numbers in the
C  user's program in compressed form.  In this format, the NDIG
C  (base JBASE) digits of the mantissa are packed two per word to
C  conserve storage.  Thus the external, packed form of a number
C  requires (NDIG+1)/2+1 words.
C
C  The unpacked format of a floating multiple precision number is as
C  follows.  The number is kept in an integer array with the first
C  element containing the exponent and each of the succeeding NDIG
C  locations containing one digit of the mantissa, expressed in base
C  JBASE.  The exponent is a power of JBASE and the implied radix point
C  is immediately before the first digit of the mantissa.  Every nonzero
C  number is normalized so that the second array element (the first
C  digit of the mantissa) is nonzero.
C
C  In both representations the sign of the number is carried on the
C  second array element only.  Elements 3,4,... are always nonnegative.
C  The exponent is a signed integer and may be as large in magnitude as
C  MXEXP (defined in FMSET).
C
C  For JBASE = 10,000 and NDIG = 4 the number -pi would have these
C  representations:
C                   Word 1         2         3         4         5
C
C         Unpacked:      1        -3      1415      9265      3590
C         Packed:        1    -31415  92653590
C
C  Because of normalization, the equivalent number of base 10
C  significant digits for an FM number may be as small as
C  LOG10(JBASE)*(NDIG-1) + 1.
C
C
C  Subroutine FMOUT performs conversion of an FM number to base 10 and
C  formats it for output as a character array.
C  The user sets JFORM1 and JFORM2 to determine the output format.
C
C  JFORM1 = 0     E   format       ( .314159M+6 )
C         = 1     1PE format       ( 3.14159M+5 )
C         = 2     F   format       ( 314159.000 )
C
C  JFORM2 = Number of significant digits to display (if JFORM1 = 0, 1).
C           If JFORM2.EQ.0 then a default number of digits is chosen.
C           The default is roughly the full precision of the number.
C  JFORM2 = Number of digits after the decimal point (if JFORM1 = 2).
C           See the FMOUT documentation for more details.
C
C
C  KW is the unit number to be used for all output from the package,
C     including error and warning messages, and trace output.
C
C
C  NTRACE and LVLTRC control trace printout from the package.
C
C  NTRACE =  0   No printout except warnings and errors.
C         =  1   The result of each call to one of the routines
C                   is printed in base 10, using FMOUT.
C         = -1   The result of each call to one of the routines
C                   is printed in internal base JBASE format.
C         =  2   The input arguments and result of each call to one
C                   of the routines is printed in base 10, using FMOUT.
C         = -2   The input arguments and result of each call to one
C                   of the routines is printed in base JBASE format.
C
C
C  LVLTRC defines the call level to which the trace is done.  LVLTRC = 1
C         means only FM routines called directly by the user are traced,
C         LVLTRC = 2 also prints traces for FM routines called by other
C         FM routines called directly by the user, etc.
C
C  In the above description internal JBASE format means the number is
C  printed as it appears in the array as an exponent followed by NDIG
C  base JBASE digits.
C
C
C  KFLAG is a condition parameter returned by the package.  Negative
C        values indicate conditions for which a warning message will
C        be printed unless KWARN = 0.  Positive values indicate
C        conditions which may be of interest but are not errors.
C        No warning message is printed if KFLAG is nonnegative.
C
C    KFLAG =  0     Normal operation.
C
C          =  1     One of the operands in FMADD or FMSUB was
C                       insignificant with respect to the other, so
C                       that the result was equal to the argument of
C                       larger magnitude.
C          =  2     In converting an FM number to a one word integer
C                       in FMM2I the FM number was not exactly an
C                       integer.  The next integer toward zero was
C                       returned.
C
C          = -1     NDIG was less than 2 or more than MXNDIG.
C          = -2     JBASE was less than 2 or more than MXBASE.
C          = -3     An exponent was out of range.
C          = -4     Invalid input argument(s) to an FM routine.
C                        UNKNOWN was returned.
C          = -5     + or - OVERFLOW was generated as a result from an
C                        FM routine.
C          = -6     + or - UNDERFLOW was generated as a result from an
C                        FM routine.
C          = -7     The input string to FMINP was not legal.
C          = -8     The output array for FMOUT was not large enough.
C          = -9     Precision could not be raised enough to provide all
C                        requested guard digits.  UNKNOWN was returned.
C          = -10    An FM input argument was too small in magnitude to
C                        convert in FMM2SP or FMM2DP.  Zero was
C                        returned.
C
C  When a negative KFLAG condition is encountered the routine calls
C  FMWARN, which uses the value of KWARN as follows.
C
C  KWARN = 0     Execution continues and no message is printed.
C        = 1     A warning message is printed and execution continues.
C        = 2     A warning message is printed and execution stops.
C
C  When an overflow or underflow is generated for an operation in which
C  an input argument was already an overflow or underflow, no additional
C  message is printed.  When an unknown result is generated and an input
C  argument was already unknown, no additional message is printed.  In
C  these cases the negative KFLAG value is still returned.
C
C
C  KRAD = 0     Causes all angles in the trigonometric functions and
C                  inverse functions to be measured in degrees.
C       = 1     Causes all angles to be measured in radians.
C
C
C  KROUND = 0   Causes all final results to be chopped (rounded toward
C                  zero).  Intermediate results are rounded.
C         = 1   Causes all results to be rounded to the nearest FM
C                  number, or to the value with an even last digit if
C                  the result is halfway between two FM numbers.
C
C
C  Here is a list of the routines in FM which are designed to be called
C  by the user.  All are subroutines except logical function FMCOMP.
C  MA, MB, MC refer to FM format numbers.
C
C  FMABS(MA,MB)         MB = ABS(MA)
C  FMACOS(MA,MB)        MB = ACOS(MA)
C  FMADD(MA,MB,MC)      MC = MA + MB
C  FMASIN(MA,MB)        MB = ASIN(MA)
C  FMATAN(MA,MB)        MB = ATAN(MA)
C  FMATN2(MA,MB,MC)     MC = ATAN2(MA,MB)
C  FMBIG(MA)            MA = Biggest FM number less than overflow.
C  FMCOMP(MA,LREL,MB)   Logical comparison of MA and MB.  LREL is a
C                            CHARACTER*2 value identifying which
C                            comparison is made.
C                            Example:  IF (FMCOMP(MA,'GE',MB)) ...
C  FMCOS(MA,MB)         MB = COS(MA)
C  FMCOSH(MA,MB)        MB = COSH(MA)
C  FMDIG(NSTACK,KST)    Find a set of precisions to use during Newton
C                            iteration for finding a simple root
C                            starting with about double precision
C                            accuracy.
C  FMDIM(MA,MB,MC)      MC = DIM(MA,MB)
C  FMDIV(MA,MB,MC)      MC = MA/MB
C  FMDIVI(MA,INT,MB)    MB = MA/INT for one word integer INT.
C  FMDP2M(X,MA)         MA = X conversion from double precision to FM.
C  FMEQU(MA,MB,NDA,NDB) MB = MA where MA has NDA digits and MB has
C                            NDB digits.
C  FMEXP(MA,MB)         MB = EXP(MA)
C  FMI2M(INT,MA)        MA = INT conversion from one word integer to FM.
C  FMINP(LINE,MA,LA,LB) MA = LINE input conversion of LINE(LA) through
C                            LINE(LB) from characters to FM.
C  FMINT(MA,MB)         MB = INT(MA) integer part of MA.
C  FMIPWR(MA,INT,MB)    MB = MA**INT raise FM number to a one word
C                            integer power.
C  FMLG10(MA,MB)        MB = LOG10(MA)
C  FMLN(MA,MB)          MB = LOG(MA)
C  FMLNI(INT,MA)        MA = LOG(INT) natural log of one word integer.
C  FMM2DP(MA,X)         X  = MA conversion from FM to double precision.
C  FMM2I(MA,INT)        INT = MA conversion from FM to one word integer.
C  FMM2SP(MA,X)         X  = MA conversion from FM to single precision.
C  FMMAX(MA,MB,MC)      MC = MAX(MA,MB)
C  FMMIN(MA,MB,MC)      MC = MIN(MA,MB)
C  FMMOD(MA,MB,MC)      MC = MA mod MB
C  FMMPY(MA,MB,MC)      MC = MA*MB
C  FMMPYI(MA,INT,MB)    MB = MA*INT multiplication by one word integer.
C  FMNINT(MA,MB)        MB = NINT(MA) nearest integer.
C  FMOUT(MA,LINE,LB)    LINE = MA conversion from FM to character.  LB
C                              is the size of array LINE.
C  FMPI(MA)             MA = pi
C  FMPRNT(MA)           Print MA using current format.
C  FMPWR(MA,MB,MC)      MC = MA**MB
C  FMSET(NPREC)         Set default values and machine-dependent
C                            variables to give at least NPREC base 10
C                            digits plus three base 10 guard digits.
C  FMSIGN(MA,MB,MC)     MC = SIGN(MA,MB) sign transfer.
C  FMSIN(MA,MB)         MB = SIN(MA)
C  FMSINH(MA,MB)        MB = SINH(MA)
C  FMSP2M(X,MA)         MA = X conversion from single precision to FM.
C  FMSQRT(MA,MB)        MB = SQRT(MA)
C  FMSUB(MA,MB,MC)      MC = MA - MB
C  FMTAN(MA,MB)         MB = TAN(MA)
C  FMTANH(MA,MB)        MB = TANH(MA)
C  FMULP(MA,MB)         MB = 1 Unit in the Last Place of MA.
C
C  For each of these routines there is also a version available for
C  which the argument list is the same but all FM numbers are in packed
C  format.  The packed versions have the same names except 'FM' is
C  replaced by 'FP' at the start of each name.
C
C
C  NOTES ON ARRAY DIMENSIONS.
C
C  The dimensions of the arrays in the FM package are defined using
C  a PARAMETER statement at the top of each routine.  The size of
C  these arrays depends on the values of parameters MXNDIG and NBITS.
C  MXNDIG is the maximum value the user may set for NDIG.
C  NBITS is the number of bits used to represent integers.
C
C  FM numbers in packed format have size LPACK, and those in unpacked
C  format have size LUNPCK.
C
C
C  PORTABILITY NOTES.
C
C  In routines FMEQU and FMSUB there is code which attempts to
C  determine if two input arrays refer to identical memory locations.
C  Some optimizing compilers assume the arrays must be distinct and
C  may remove code which would then be redundant.  This code removal
C  could cause errors, so the tests are done in a way which should
C  keep the compiler from removing code. The current version works
C  correctly on all compilers tested.  Computing SIN(1.0) in radian
C  mode should reveal whether other compilers handle it correctly.
C  If there is a problem, SIN(1) gives 0.999... instead of 0.841....
C  To fix such a problem, MB can be copied to a local temporary array
C  and then negated in FMSUB before calling FMADD2.  For FMEQU
C  simply set KSAME = 0 after the section which tries to determine if
C  MA and MB are the same array. This makes both routines run slower.
C  A simpler fix which often works is to re-compile at a lower
C  optimization (e.g., OPT=0).
C
C  In FMSET there is machine-dependent code which attempts to
C  approximate the largest one word integer value.  The current code
C  works on all machines tested, but if an FM run fails, check the
C  MAXINT loop in FMSET in addition to the three routines mentioned
C  above.  Values for SPMAX and DPMAX are also defined which should
C  be set to values near overflow for single precision and double
C  precision.
C
C  Using the CFT77 compiler on a Cray X-MP computer there are
C  problems using a large value for JBASE when the default 46-bit
C  integer arithmetic mode is used.  In particular, FMSET chooses
C  too large a JBASE value since some of the arithmetic in the
C  MAXINT loop is done with 64-bit precision.  Setting JBASE = 10**6
C  or less may be ok, but the preferred solution is to select the
C  64-bit integer compiler option.  Then JBASE = 10**9 can be used.
```
