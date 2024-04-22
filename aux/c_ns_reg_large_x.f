!     Large-x expansion of the regular part of the third order non-singlet
!     coefficient functions in DIS. They are coded to be used with the
!     Fortran files of Moch, Vermaseren, and Vogt given in
!
!     S. Moch, J. Vermaseren and A. Vogt, hep-ph/0209100
!     J. Vermaseren, A. Vogt and S. Moch, hep-ph/0504242
!     A. Vogt, J. Vermaseren and S. Moch, arXiv:0812.4168 [hep-ph]
!
!     and replace the routines X2NP3A, X3NM3A, and XLNP3A in the large-x
!     limit. In practice the agreement between this expansion and the
!     exact expressions are at a relative 10^-4 starting from 1-y =
!     10^-7.
!
!     Please cite V. Bertone, A. Karlberg, arXiv:2404:XXXXX if using the
!     code.

      FUNCTION X2NP3A_large_x (Y, DL, NF)
      IMPLICIT NONE
      INTEGER NF
      DOUBLE PRECISION Y, DL, DL1
      DOUBLE PRECISION X2NP3A_large_x
      DOUBLE PRECISION DL1VAL
      
      DL1 = DL1VAL(Y, DL)
      
      X2NP3A_large_x = 5894.634952596363d0*DL1 + 2319.655820717043d0
     $     *DL1**2 -1787.0418273217867d0*DL1**3 + 348.44444444444446d0
     $     *DL1**4 -18.962962962962962d0*DL1**5 +(1199.690656381538d0
     $     *DL1 -787.5420539087113d0*DL1**2 +146.69958847736626d0*DL1
     $     **3 -7.901234567901234d0*DL1**4) *NF +(-65.15652656374832d0
     $     *DL1 +14.617283950617283d0*DL1 **2 -0.7901234567901234d0*DL1
     $     **3)*NF**2
      
      RETURN
      END FUNCTION

      FUNCTION X3NM3A_large_x (Y, DL, NF)
      IMPLICIT NONE
      INTEGER NF
      DOUBLE PRECISION Y, DL, DL1
      DOUBLE PRECISION X3NM3A_large_x
      DOUBLE PRECISION DL1VAL
      
      DL1 = DL1VAL(Y, DL)
      
      X3NM3A_large_x = 6889.89378009207d0*DL1 + 1581.208087862925d0*DL1
     $     **2 -1609.6420585153533d0*DL1**3 + 329.48148148148147d0*DL1
     $     **4 -18.962962962962962d0*DL1**5 +(859.380948706157d0*DL1 -
     $     675.1899801124716d0*DL1**2 +134.0576131687243d0*DL1**3 -
     $     7.901234567901234d0*DL1**4)*NF +(-50.144180884735974d0*DL1 +
     $     12.246913580246913d0*DL1**2 -0.7901234567901234d0*DL1**3)*NF
     $     **2
      
      END FUNCTION

      FUNCTION XLNP3A_large_x (Y, DL, NF)
      IMPLICIT NONE
      INTEGER NF
      DOUBLE PRECISION Y, DL, DL1
      DOUBLE PRECISION XLNP3A_large_x
      DOUBLE PRECISION DL1VAL
      
      DL1 = DL1VAL(Y, DL)
      
      XLNP3A_large_x = -995.2588274957099d0*DL1 + 738.4477328541179d0
     $     *DL1**2- 177.39976880643366d0*DL1**3 + 18.962962962962962d0
     $     *DL1**4+(340.3097076753812d0*DL1 - 112.35207379623978d0*DL1
     $     **2 +12.641975308641975d0*DL1**3)*NF + (-15.012345679012345d0
     $     *DL1 +2.3703703703703702d0*DL1**2)*NF**2
      
      END FUNCTION

!     For Y values close to 1, use the series expansion close 
!     to 1-Y=0 instead of full value, for numerical convergence
!     
      FUNCTION DL1VAL (Y, DL)
      IMPLICIT NONE
      DOUBLE PRECISION D2, D24, D2880, Y, DL, DL1VAL
      PARAMETER (D2=0.5D0, D24=1.0D0/24.0D0, D2880=1.0D0/2880.0D0)

      IF (ABS(DL).LT.1D-4) THEN
         DL1VAL = LOG(-DL) +  DL*D2 + DL**2*D24 - DL**4*D2880
      ELSE
         DL1VAL = LOG(1.0D0 - Y)
      ENDIF

      RETURN
      END FUNCTION


