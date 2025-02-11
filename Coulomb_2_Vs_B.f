C     Este progrma calcula la Integral de Coulomb para un QD eje-simÇtrico
C     El c†lculo solo usa el estado fundamental para electr¢n y hueco (l=0, l: momentum angular)
C     Se usan las funciones ELLIPTICAS COMPLETAS DE PRIMERA CLASE
C     De ese modo se hace la integral angular para phi_e y phi_h entre 0 y 2 pi.

      implicit double precision (A-H,O-Z)
      parameter(np=19881)
      dimension A_e(np,3),A_h(np,3)
      dimension A_ee(np,3),A_hh(np,3)
      dimension AMatrix_e(np,43),AMatrix_h(np,43)

      OPEN(1,file='electron_Vs_B_F_mas_40.txt',status='unknown')
      OPEN(2,file='hueco_Vs_B_F_mas_40.txt',status='unknown')
      OPEN(3,file='tabla_mas_40.dat',status='unknown')

      OPEN(10,file='duque_ee_v3.txt',status='unknown')
      OPEN(11,file='duque_hh_v3.txt',status='unknown')

      pi =4.0d0*datan(1.0d0)
      dro_e = 0.5d-9
      dz_e  = 0.5d-9
      dro_h = 0.5d-9
      dz_h  = 0.5d-9


C     Leyendo la matrix de electrones y huecos para diferentes campos magnÇticos (desde 0 a 20 en pasos de 0.5, o sea 41 valores
      do 500 i=1,np
      read(1,*)(AMatrix_e(i,j),j=1,43)
      read(2,*)(AMatrix_h(i,j),j=1,43)
  500 continue


C     en el siguiente loop se corren los 41 valores de campo magnÇtico
      BE = 0.0d0
      do 1000 k=1,41
C     lectura de posiciones y funciones de onda del estado fundamental de electr¢n
C     tambiÇn se identifica la funci¢n de onda m†xima de la lista
      elec_max = 0.0d0
      do 1 i=1,np
      A_ee(i,1) = AMatrix_e(i,1)
      A_ee(i,2) = AMatrix_e(i,2)
      A_ee(i,3) = AMatrix_e(i,2+k)
      elec_max = max(elec_max,A_ee(i,3))
  1   continue
c      write(*,*)elec_max

C     selecci¢n de las funciones de onda que son mayores o iguales que el 1%
C     del valor m†ximo
      nt=1
      do 100 i=1,np
      if(A_ee(i,3).ge.0.01d0*elec_max)then
      A_e(nt,1)=A_ee(i,1)
      A_e(nt,2)=A_ee(i,2)
      A_e(nt,3)=A_ee(i,3)
      nt=nt+1
      else
      endif
  100 continue
      n1 = nt-1
c      write(*,*)n1

      do 30 i=1,n1
      write(10,*)A_e(i,1),A_e(i,2),A_e(i,3)
  30   continue

C     lectura de posiciones y funciones de onda del estado fundamental de hueco
C     tambiÇn se identifica la funci¢n de onda m†xima de la lista
      hole_max = 0.0d0
      do 10 i=1,np
      A_hh(i,1) = AMatrix_h(i,1)
      A_hh(i,2) = AMatrix_h(i,2)
      A_hh(i,3) = AMatrix_h(i,2+k)
      hole_max = max(hole_max,A_hh(i,3))
  10   continue
c      write(*,*)hole_max

C     selecci¢n de las funciones de onda que son mayores o iguales que el 1%
C     del valor m†ximo
      nq=1
      do 200 i=1,np
      if(A_hh(i,3).ge.0.01d0*hole_max)then
      A_h(nq,1)=A_hh(i,1)
      A_h(nq,2)=A_hh(i,2)
      A_h(nq,3)=A_hh(i,3)
      nq=nq+1
      else
      endif
  200 continue
      n2 = nq-1
c      write(*,*)n2
c      read(*,*)absbs

      do 40 i=1,n2
      write(11,*)A_h(i,1),A_h(i,2),A_h(i,3)
  40   continue


C     n1: n£mero de funciones de onda de electr¢n que son mayores
C     o iguales que el 1% del valor m†ximo
C     n2: n£mero de funciones de onda de hueco que son mayores
C     o iguales que el 1% del valor m†ximo


C     normalizar las funciones de onda    Int. Psi^2*2*pi*ro*dro*dz
      suma_e = 0.0d0
      do 2 i=1,n1
      ro_e = A_e(i,1)*1.0d-9
      dVe  = 2.0d0*pi*ro_e*dro_e*dz_e
      psi_e= A_e(i,3)
      suma_e = suma_e+psi_e**2*dVe
c      write(*,*)A_e(i,1),A_e(i,2),A_e(i,3),psi_e**2*dVe/2.0d0
  2   continue
      xnorm_e = suma_e


      suma_h = 0.0d0
      do 3 i=1,n2
      ro_h = A_h(i,1)*1.0d-9
      dVh  = 2.0d0*pi*ro_h*dro_h*dz_h
      psi_h= A_h(i,3)
      suma_h = suma_h+psi_h**2*dVh
  3   continue
      xnorm_h = suma_h

C     Calculemos la interacci¢n de Coulomb
C      go to 20
      suma_C = 0.0d0
      do 4 i=1,n1
      do 5 j=1,n2
      ro_e = A_e(i,1)*1.0d-9
      ro_h = A_h(j,1)*1.0d-9
      z_e = A_e(i,2)*1.0d-9
      z_h = A_h(j,2)*1.0d-9
      
      r = dsqrt((ro_e-ro_h)**2+(z_e-z_h)**2)
      if(r.eq.0.0d0)go to 5
      psi_e = A_e(i,3)/dsqrt(xnorm_e)
      psi_h = A_h(j,3)/dsqrt(xnorm_h)
      
c      dVe = ro_e*2.0d0*pi*dro_e*dz_e
c      dVh = ro_h*2.0d0*pi*dro_h*dz_h
c      suma_C = suma_C+psi_e**2*psi_h**2*dVe*dVh/r
C     las tres lineas anterior las vamos a cambiar por las siguientes lineas
C     donde la integral angular en phi_e y phi_h  de  (1/r) entre 0 y 2 pi queda expresada
C     en tÇrminos de funciones ELIPTICAS
      dVe = ro_e*dro_e*dz_e        !notese que he quitado el factor 2*pi
      dVh = ro_h*dro_h*dz_h        !notese que he quitado el factor 2*pi
C     SUBROUTINE COMELP(XK,CK,CE)  Calcula la funci¢n Elptica
      rp = 4.0d0*ro_e*ro_h/r
      XK = rp/(1.0d0+rp)
      call COMELP(XK,CK,CE)
      xEliptica = CK
      xIntegral_Angular = 8.0d0*pi/r*xEliptica/dsqrt(1.0d0+rp)
      suma_C = suma_C+psi_e**2*psi_h**2*dVe*dVh*xIntegral_Angular

  5   continue
  4   continue
      xcoulomb = suma_C
  20  continue



C     Calculo de la posici¢n z_e del electr¢n, z_h del hueco y la integral de overlap
      suma_ze = 0.0d0
      suma_zh = 0.0d0
      suma_overlap = 0.0d0
      do 41 j=1,np
      ro_e = A_ee(j,1)*1.0d-9
      ro_h = A_hh(j,1)*1.0d-9
      z_e = A_ee(j,2)*1.0d-9
      z_h = A_hh(j,2)*1.0d-9

      psi_e = A_ee(j,3)/dsqrt(xnorm_e)
      psi_h = A_hh(j,3)/dsqrt(xnorm_h)

      dVe = ro_e*2.0d0*pi*dro_e*dz_e
      dVh = ro_h*2.0d0*pi*dro_h*dz_h
      suma_ze = suma_ze+psi_e**2*z_e*dVe
      suma_zh = suma_zh+psi_h**2*z_h*dVh
      suma_overlap = suma_overlap+psi_e*psi_h*dVe
  41   continue
      xIntegral_overlap=suma_overlap**2

C     constantes y valores usados en los c†lculos
      h = 6.62606876d-34             ! Planck's constant [J*s]
      hbar = h/(2.0d0*pi)            ! Planck's constant reduced [J*s]
      epsilon0 = 8.854187817d-12     ! Vacuum permittivity [F/m] = [A*s/(V*m)]
      xme0 = 9.10938188d-31          ! Electron mass [kg]
      c = 299792458                  ! speed of light in m/s
      e = 1.60217662d-19             ! elementary charge in coulombs
      xJeV = 1.0d0/1.60217646d-19    ! Joule in eV

C     GaAs related constants,
      xme_GaAs = 0.067d0*xme0        ! GaAs effective electron mass [kg]
      xmhh_GaAs = 0.51d0*xme0        ! GaAs effective heavy hole mass [kg]
      epsilon_GaAs = 13.1d0*epsilon0 ! GaAs permittivity

      xm_Ex = xme_GaAs*xmhh_GaAs/(xme_GaAs+xmhh_GaAs)          ! Exciton effective mass
      r_Ex_GaAs = 4.0d0*pi*epsilon_GaAs*hbar**2/(xm_Ex*e**2)   ! GaAs exciton Bohr radius in [m]
      E_Ry_GaAs = xm_Ex*e**4/(2.0d0*(2.0d0*epsilon_GaAs*h)**2) ! Rydberg constant in GaAs bulk in [J]
C     *******************************************
      constante = e**2/(4.0d0*pi*epsilon_GaAs)*xJeV*1.0d3

      E_eh = constante*xcoulomb
      E_eh = E_eh/2.0d0         ! ESTE PASO NO LO ENTIENDO. NO SE POR QUê CHRISTIAN DIVIDE POR DOS

      write(3,*)BE,E_eh,suma_ze,suma_zh,xIntegral_overlap
      write(*,*)BE,E_eh,suma_ze,suma_zh,xIntegral_overlap
      BE = BE+0.5d0

 1000 continue
      WRITE(*,*)'DEME CUALQUIER COSA PARA TERMINAR EL PROGRAMA Coulomb2'
      read(*,*)SSSSS
      stop
      end


C       ==================================================
        SUBROUTINE COMELP(XK,CK,CE)
C       COMENTARIO: lo que en este programa se evalua poniendo la variable como x
C       en el mathematica hay que poner x^2. Por esa raz¢n se redefine la variabble
C       en la usando la raiz cuadrada. Ver la l°nea insertada despuÇs del
C       IMPLICIT DOUBLE PRECISION

C
C
C       ==================================================
C       Purpose: Compute complete elliptic integrals K(k)
C                and E(k)
C       Input  : K  --- Modulus k ( 0 Û k Û 1 )
C       Output : CK --- K(k)
C                CE --- E(k)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        HK = dsqrt(XK)

        PK=1.0D0-HK*HK
        IF (HK.EQ.1.0) THEN
           CK=1.0D+300
           CE=1.0D0
        ELSE
           AK=(((.01451196212D0*PK+.03742563713D0)*PK
     &        +.03590092383D0)*PK+.09666344259D0)*PK+
     &        1.38629436112D0
           BK=(((.00441787012D0*PK+.03328355346D0)*PK+
     &        .06880248576D0)*PK+.12498593597D0)*PK+.5D0
           CK=AK-BK*DLOG(PK)
           AE=(((.01736506451D0*PK+.04757383546D0)*PK+
     &        .0626060122D0)*PK+.44325141463D0)*PK+1.0D0
           BE=(((.00526449639D0*PK+.04069697526D0)*PK+
     &        .09200180037D0)*PK+.2499836831D0)*PK
           CE=AE-BE*DLOG(PK)
        ENDIF
        RETURN
        END


