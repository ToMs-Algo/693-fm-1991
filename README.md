# Algorithm 693: a FORTRAN package for floating-point multiple-precision arithmetic

> David M. Smith. 1991.  
> Algorithm 693: a FORTRAN package for floating-point multiple-precision arithmetic.  
> ACM Trans. Math. Softw. 17, 2 (June 1991), 273–283.  
> https://doi.org/10.1145/108556.108585

- For latest version see: [Multiple Precision Computation](https://dmsmith.lmu.build/)


## FM Package

```for
C     FM 1.0                  David M Smith                  1-06-90
C
C  The FM package performs floating point multiple precision arithmetic.
```

### List of the Routines

```for
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
```


## Copyright and License Agreement

See: [copyright.htm](copyright.htm)

### ACM Copyright Notice

Permission to make digital or hard copies of part or all of this work for personal or classroom use is granted without fee provided that copies are not made or distributed for profit or commercial advantage and that copies bear this notice and the full citation on the first page.
Copyrights for components of this work owned by others than ACM must be honored. Abstracting with credit is permitted.
To copy otherwise, to republish, to post on servers, or to redistribute to lists, requires prior specific permission and/or a fee.  
Request permissions from Publications Dept, ACM Inc., fax +1 (212) 869-0481, or permissions@acm.org.

Collected Algorithms, Special Edition CD-ROM  
Copyright © 2000 by the Association for Computing Machinery, Inc. 1-58113-338-5/00/11 ... $5.00
