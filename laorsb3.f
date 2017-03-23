c
c     laorsb3.f
c
c     Subroutine for calculating laor3 relativistically broadened
c     emission lines from accretion discs with a twice-broken power
c     law emissivity profile
c
c     v2.0, November 2014, D.R. Wilkins
c
      SUBROUTINE laorsb3(ear, ne, param, nenrgy, nrad, nteta, efil, 
     &                  ebin, rad, incl, trs, bin, start, end, fstart, 
     &                  fend, photar)


      INTEGER ne, nenrgy, nrad, nteta

      REAL param(9), ear(0:ne), photar(ne)
      REAL efil(0:*), ebin(*), rad(*), incl(*), trs(*), bin(*)
      REAL start(*), end(*), fstart(*), fend(*)

C Subroutine to do the actual calculation of the LAOR3 model
C Parameters :
C       ear      r           i: energies for output
C       ne       i           i: number of energies
C       param    r           i: model parameters
C       nenrgy   i           i: number of tabulated energies
C       nrad     i           i: number of tabulated radii
C       nteta    i           i: number of tabulated inclination angles
C       efil     r           i: energies on which model is calculated
C       ebin     r           i: tabulated energies (relative to line center)
C       rad      r           i: tabulated radii
C       incl     r           i: tabulated inclination angles
C       trs      r           i: tabulated transfer function
C       bin      r           w: workspace for accumulating model spectrum
C       start    r           w: workspace array used for rebinning
C       end      r           w: workspace array used for rebinning
C       fstart   r           w: workspace array used for rebinning
C       fend     r           w: workspace array used for rebinning
C       photar   r           r: final model spectrum

      REAL smo(3)
      REAL inclin, r, wlow, whgh, ems, spm, area, r0

      INTEGER i, j, ioff, inclow, iro, iri, ismo, ir0


c Initialize the flux array and set the observed energies.
c Note the use of the high energy bin and the low energy
c bin to suppress error messages from the interpolation
c program.

      efil(0) = 0.
      efil(1) = (3*ebin(1)-ebin(2)) * param(1) / 2.
      DO i = 2, nenrgy
         efil(i) = (ebin(i-1)+ebin(i)) * param(1) / 2.
      ENDDO
      efil(nenrgy+1) = (3*ebin(nenrgy)-ebin(nenrgy-1)) * param(1) / 2.
      efil(nenrgy+2) = 1.e6
      DO i = 1, nenrgy+2
         bin(i) = 0.
      ENDDO

c  selecting the integer values of the radial bins which are at the inner and
c  outer disk boundary (according to inner and outer radii)

      if(param(6).lt.param(3))param(6)=param(3)
      if(param(6).gt.param(4))param(6)=param(4)
      iro=1
      iri=nrad

      DO i=1,nrad-1
         IF(param(4).LE.rad(i).AND.param(4).GT.rad(i+1)) iro=i
         IF(param(3).LE.rad(i).AND.param(3).GT.rad(i+1)) iri=i
         IF(param(6).LE.rad(i).AND.param(6).GT.rad(i+1)) ir0=i
         IF(param(8).LE.rad(i).AND.param(8).GT.rad(i+1)) ir1=i
      ENDDO


c Set the cos(inclination) variable.

      inclin = cos(param(5)*3.14159/180.)

c  find the tabulated inclinations above and below that requested
c  and the associated weights

      i = 1
      DO WHILE( inclin.GT.incl(i) .AND. i.LT.nteta )
         i = i + 1
      ENDDO
      inclow = i - 1
      wlow = (incl(inclow+1)-inclin)/(incl(inclow+1)-incl(inclow))
      whgh = (inclin-incl(inclow))/(incl(inclow+1)-incl(inclow))

c Summing the profiles of individual rings according according to a given
c weight function

      area=0.
      r0=rad(ir0)
      r1=rad(ir1)
      DO j=iro,iri
         r=rad(j)
         IF(j.EQ.iro) THEN
            r=(rad(iro)+rad(iro+1))/2.
            area=param(4)**2-r**2
         ELSE
            IF(j.EQ.iri) THEN
               r=(rad(iri-1)+rad(iri))/2.
               area=r**2-param(3)**2
            ELSE
               area=((rad(j-1)+rad(j))/2.)**2-
     &           ((rad(j)+rad(j+1))/2.)**2
            ENDIF
         ENDIF

c  ems-the radial weight function
c      print *,'ir0,rad(ir0)',ir0,rad(ir0)
       if(r.lt.rad(ir0))  ems=area*r**(-param(2))
       if(r.eq.rad(ir0))  ems=(param(6)-rad(ir0+1))/(r-rad(ir0+1))
     &      *area*r**(-param(2))
     &      + (1.-(param(6)-rad(ir0+1))/(r-rad(ir0+1)))
     &      * area*param(6)**(-param(2)+param(7))*r**(-param(7))
       if(r.gt.rad(ir0).and.r.lt.rad(ir1))  
     &      ems=area*param(6)**(-param(2)+param(7))*r**(-param(7))
       if(r.eq.rad(ir1)) ems=(param(8)-rad(ir1+1))/(r-rad(ir1+1))
     &      * area*param(6)**(-param(2)+param(7))*r**(-param(7))
     &      + (1.-(param(8)-rad(ir1+1))/(r-rad(ir1+1)))
     &      * area* 
     &     param(6)**(-param(2)+param(7))*param(8)**(-param(7)+param(9))
     &     *r**(-param(9))
       if(r.gt.rad(ir1)) ems=area
     &     *param(6)**(-param(2)+param(7))*param(8)**(-param(7)+param(9)
     &     )*r**(-param(9))


c  sum in the required transfer functions

         ioff = ((j-1)*nteta + inclow - 1)*nenrgy

         DO i=1,nenrgy
            bin(i+1) = bin(i+1) + wlow*trs(ioff+i)*ems
         ENDDO

         ioff = ioff + nenrgy
         DO i=1,nenrgy
            bin(i+1) = bin(i+1) + whgh*trs(ioff+i)*ems
         ENDDO

      ENDDO

c  Perform 3-pt smoothing to reduce the small wiggles which result
c  from the finite number of radial bins

      DO i = 1, 2
         smo(i) = bin(i+1)
      ENDDO
      ismo = 3
      DO i = 2, nenrgy-1
         smo(ismo) = bin(i+2)
         ismo = ismo + 1
         IF(ismo.EQ.4) ismo = 1
         bin(i+1) = 0.
         DO j = 1, 3
            bin(i+1) = bin(i+1) + smo(j)
         ENDDO
         bin(i+1) = bin(i+1)/3.
      ENDDO

c  Rebin onto passed energies

      CALL inibin((nenrgy+2), efil, ne, ear, start, end,
     &            fstart, fend, 0)
      CALL erebin((nenrgy+2), bin, ne, start, end, fstart, fend, photar)

c normalise values to total

      spm = 0.
      DO i = 1, ne
         spm = spm + photar(i)
      ENDDO
      IF ( (spm .NE. 0.) .AND. (spm .NE. 1.) ) THEN
         DO i = 1, ne
            photar(i) = photar(i)/spm
         ENDDO
      ENDIF

      RETURN
      END




