
      PROGRAM ROOTS
C                                                           3-21-90
C
C  Test program to compute the roots of the 10th degree Legendre
C  polynomial to 50 places.  Since the function is even, only the
C  positive roots are printed.
C
C  Some integer coefficients used in this program have magnitudes
C  about one million, so it should not be run on a machine with less
C  than 20-bit words unless the integers are broken into small factors.
C
C  Local arrays in this routine are kept in unpacked format.
C
C
C  The output from this program should be:
C
C
C ROOT NUMBER 1:
C       .14887433898163121088482600112971998461756485942069
C
C ROOT NUMBER 2:
C       .43339539412924719079926594316578416220007183765625
C
C ROOT NUMBER 3:
C       .67940956829902440623432736511487357576929471183481
C
C ROOT NUMBER 4:
C       .86506336668898451073209668842349304852754301496533
C
C ROOT NUMBER 5:
C       .97390652851717172007796401208445205342826994669238
C
C
      PARAMETER ( MXNDIG=256 , NBITS=32 ,
     *          LPACK  = (MXNDIG+1)/2 + 1 , LUNPCK = (6*MXNDIG)/5 + 20,
     *          LMWA   = 2*LUNPCK         , LJSUMS = 8*LUNPCK,
     *          LMBUFF = ((LUNPCK+3)*(NBITS-1)*301)/2000 + 6  )
C
      CHARACTER *80 RESULT(5)
      DIMENSION MX(LUNPCK),MXSQ(LUNPCK),MFX(LUNPCK),MFPX(LUNPCK),
     *          MC(LUNPCK),NSTACK(19)
      DOUBLE PRECISION X,XSQ,FX,FPX,GUESS(5),C(6)
C
      COMMON /FMUSER/ NDIG,JBASE,JFORM1,JFORM2,KRAD,
     *                KW,NTRACE,LVLTRC,KFLAG,KWARN,KROUND
C
      DOUBLE PRECISION DPMAX
C
      COMMON /FM/ MWA(LMWA),NCALL,MXEXP,MXEXP2,KARGSW,KEXPUN,KEXPOV,
     *            KUNKNO,IUNKNO,RUNKNO,MXBASE,MXNDG2,KSTACK(19),
     *            MAXINT,SPMAX,DPMAX
C
      COMMON /FMSAVE/ NDIGPI,NJBPI,NDIGE,NJBE,NDIGLB,NJBLB,NDIGLI,NJBLI,
     *                MPISAV(LUNPCK),MESAV(LUNPCK),MLBSAV(LUNPCK),
     *                MLN1(LUNPCK),MLN2(LUNPCK),MLN3(LUNPCK),
     *                MLN4(LUNPCK)
C
      DATA GUESS/ .14887D0, .43340D0, .67941D0, .86506D0, .97391D0/
      DATA C/ 46189.D0, -109395.D0, 90090.D0, -30030D0, 3465.D0, -63.D0/
C
      RESULT(1) = ' .14887433898163121088482600112971998461756485942069'
      RESULT(2) = ' .43339539412924719079926594316578416220007183765625'
      RESULT(3) = ' .67940956829902440623432736511487357576929471183481'
      RESULT(4) = ' .86506336668898451073209668842349304852754301496533'
      RESULT(5) = ' .97390652851717172007796401208445205342826994669238'
      NOGOOD = 0
C
C             Set precision to give at least 50 significant digits, and
C             initialize all FM parameters.  Set the print format to
C             give fixed point with 50 places.
C
      CALL FPSET(50)
      JFORM1 = 2
      JFORM2 = 50
C
      KW = 11
      OPEN(KW,FILE='ROOTS.OUT')
C
      DO 160 J = 1, 5
C
C             Generate the initial FM approximation using double
C             precision Newton iteration.
C
         X = GUESS(J)
         DO 110 K = 1, 4
            XSQ = X*X
            FX = ((((C(1)*XSQ+C(2))*XSQ+C(3))*XSQ+C(4))*XSQ+C(5))*XSQ
     *           +C(6)
            FPX = ((((10.0D0*C(1)*XSQ+8.0D0*C(2))*XSQ+6.0D0*C(3))*XSQ
     *            +4.0D0*C(4))*XSQ+2.0D0*C(5))*X
            IF (FPX.NE.0.0D0) X = X - FX/FPX
  110    CONTINUE
C
C             Convert to FM numbers and use increasing precision
C             until 50 digits are correct.
C
         NDSAVE = NDIG
         CALL FMDIG(NSTACK,KST)
         CALL FMDP2M(X,MX)
         DO 140 K = 1, KST
            NDIG = NSTACK(K)
            CALL FMMPY(MX,MX,MXSQ)
C
C             Convert the coefficients to integers and then to FM
C             numbers to avoid double precision rounding errors.
C
            KC = NINT(C(1))
            CALL FMI2M(KC,MFX)
            DO 120 L = 2, 6
               CALL FMMPY(MXSQ,MFX,MFX)
               KC = NINT(C(L))
               CALL FMI2M(KC,MC)
               CALL FMADD(MFX,MC,MFX)
  120       CONTINUE
C
            KC = 10*NINT(C(1))
            CALL FMI2M(KC,MFPX)
            DO 130 L = 2, 5
               CALL FMMPY(MXSQ,MFPX,MFPX)
               KC = (12-2*L)*NINT(C(L))
               CALL FMI2M(KC,MC)
               CALL FMADD(MFPX,MC,MFPX)
  130       CONTINUE
            CALL FMMPY(MFPX,MX,MFPX)
C
            CALL FMDIV(MFX,MFPX,MC)
            CALL FMSUB(MX,MC,MX)
            NDIG = NDSAVE
  140    CONTINUE
C
         WRITE (KW,150) J
  150    FORMAT(/' ROOT NUMBER',I2,':')
         CALL FMPRNT(MX)
         CALL CHECK(RESULT(J),MX,NOGOOD,KW)
  160 CONTINUE
C
      IF (NOGOOD.EQ.0) THEN
          WRITE (KW,170)
  170     FORMAT(//' ROOTS COMPLETED.  NO ERRORS.'//)
      ELSE IF(NOGOOD.EQ.1) THEN
          WRITE (KW,180)
  180     FORMAT(//' ROOTS COMPLETED.  1 ERROR.'//)
      ELSE
          WRITE (KW,190) NOGOOD
  190     FORMAT(//' ROOTS COMPLETED.',I4,' ERRORS.'//)
      ENDIF
C
      STOP
      END
      SUBROUTINE CHECK(RESULT,MX,NOGOOD,KW)
C
C  Check to see that the computed result MX agrees with the
C  correct output RESULT.
C
C  NOGOOD counts the number of cases where errors were found.
C
      PARAMETER ( MXNDIG=256 , NBITS=32 ,
     *          LPACK  = (MXNDIG+1)/2 + 1 , LUNPCK = (6*MXNDIG)/5 + 20,
     *          LMWA   = 2*LUNPCK         , LJSUMS = 8*LUNPCK,
     *          LMBUFF = ((LUNPCK+3)*(NBITS-1)*301)/2000 + 6  )
C
      CHARACTER *80 RESULT
      CHARACTER LINE(80)
      DIMENSION MX(LUNPCK)
C
      CALL FMOUT(MX,LINE,80)
      DO 120 J = 1, 80
         IF (LINE(J).NE.RESULT(J:J)) THEN
             NOGOOD = NOGOOD + 1
             WRITE (KW,110) RESULT
  110        FORMAT(/' ERROR IN ROOTS.  THE CORRECT RESULT SHOULD BE:'
     *              //6X,A/)
             RETURN
         ENDIF
  120 CONTINUE
      RETURN
      END
