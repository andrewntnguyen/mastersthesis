      program fft

      INCLUDE 'fft.dim'
c
      integer p,q
c
      dimension data(2*npts1*npts2*npts3),nn(3),nn2(2)
      dimension xk(3),xk2(3)
      dimension a(3,3),g1(3,3,3,3)
c
      dimension aux6(6),aux33(3,3)
c
      dimension minv1(3),minv2(3)
c
      dimension ddisgrad(3,3),ddisgradim(3,3)
cx      dimension ddisgradav(3,3)
c
      dimension work(3,3,NPTS1,NPTS2,NPTS3)
      dimension workim(3,3,NPTS1,NPTS2,NPTS3)
c
      dimension epav(3,3),edotpav(3,3)
c
      dimension disgradmacroactual(3,3),disgradmacrot(3,3)
      dimension velgradmacro(3,3)
c
      dimension cg(3,3,3,3),cg66aux(6,6)
caps
      dimension cinvg(3,3,3,3),eel(3,3)
caps
c
      dimension fbar(3,3)
cx      dimension velmax(3)
c
      pi=4.*atan(1.)
c
      open(55,file='vm.out',status='unknown')
      open(56,file='str_str.out',status='unknown')
      open(21,file='err.out',status='unknown')
      open(25,file='conv.out',status='unknown')
cc      open(35,file='activ.out',status='unknown')

      write(21,*) 'IT     ERRE       ERRS       SVM'
cc     #    ACT ...'
c
      nn(1)=npts1
      nn(2)=npts2
      nn(3)=npts3
c
      nn2(1)=npts1
      nn2(2)=npts2
c
      prodnn=float(nn(1)*nn(2)*nn(3))
      wgt=1./prodnn
c
      ur0=0
      open(ur0,file='fft_X.in',status='old')
      UR1=1      ! FILECRYS
      UR2=2      ! FILETEXT
c
      call vpsc_input
ch
      IF(IWFIELDS.EQ.1) then
       open(91,file='efield.out',status='unknown')
       open(92,file='sfield.out',status='unknown')
       open(93,file='dfield.out',status='unknown')
caps
       open(94,file='elfield.out',status='unknown')
c       open(95,file='sgtot.out',status='unknown')
caps
      ENDIF
ch
cc
cc     INITIALIZATION
cc
      do ii=1,3
      do jj=1,3
      disgradmacrot(ii,jj)=0.
      disgradmacro(ii,jj)=udot(ii,jj)*tdot
      enddo
      enddo

      do 777 i=1,npts1
      do 777 j=1,npts2
      do 777 k=1,npts3
cth
      jgr=jgrain(i,j,k)
cth
      do ii=1,3
      do jj=1,3
      edotp(ii,jj,i,j,k)=0.
cth      ept(ii,jj,i,j,k)=0.
      if(jgr.gt.0) then
       ept(ii,jj,i,j,k)=eth(ii,jj,jgr)
      else
       ept(ii,jj,i,j,k)=0.
      endif
cth
chardcx      disgrad(ii,jj,i,j,k)=0.
      disgrad(ii,jj,i,j,k)=0.
cth
cth      if(jgr.gt.0) then
cth       disgrad(ii,jj,i,j,k)=eini(ii,jj,jgr)
cth      else
cth       disgrad(ii,jj,i,j,k)=0.
cth      endif
cth
cx      sg(ii,jj,i,j,k)=0.
      enddo
      enddo
cx
cx    ELASTIC INITIAL GUESS FOR THE STRESS
cx
cg
      jph=jphase(i,j,k)
      if(igas(jph).eq.0) then
cg
       do ii=1,6
       do jj=1,6
        cg66aux(ii,jj)=cg66(ii,jj,i,j,k)
       enddo
       enddo
c
      call chg_basis(aux6,aux33,cg66aux,cg,3,6)
c
       do ii=1,3
       do jj=1,3
       sg(ii,jj,i,j,k)=0.
       do kk=1,3
       do ll=1,3
       sg(ii,jj,i,j,k)=sg(ii,jj,i,j,k)+
     #                cg(ii,jj,kk,ll)*disgradmacro(kk,ll)
cth       sg(ii,jj,i,j,k)=sg(ii,jj,i,j,k)+
cth     #                cg(ii,jj,kk,ll)*disgrad(kk,ll,i,j,k)
       enddo
       enddo
       enddo
       enddo
cg
       else     !   igas else
cg
       do ii=1,3
       do jj=1,3
       sg(ii,jj,i,j,k)=0.
       enddo
       enddo
cg
       endif    ! igas endif
cg
777   continue

cc      evm=0.

      evm=dvm*tdot

      write(55,*)
     #'EVM         EVMP        DVM         DVMP        SVM         SVM1'

cx      write(56,'(a234)')
cx     #    ' EVM           E11         E22         E33       EPVM
cx     #EP11         EP22         EP33         DVM         D11         D22
cx     #     D33        DPVM        DP11      DP22      DP33
cx     #SVM         S11         S22         S33'

c      write(95,format(5f7,3(e11.3,1x))) sg23avg,sg31avg,sg12avg
      write(56,'(t2,a,t14,a,t26,a,t38,a,t50,a,t62,a,t74,a,t86,a,t98,a,
     #           t110,a,t122,a,t134,a,t146,a,t158,a,t170,a,t182,a,
cg     #           t194,a,t206,a,t218,a,t230,a)') 'EVM','E11','E22',
     #    t194,a,t206,a,t218,a,t230,a,t242,a,t254,a,t266,a,t278,a)')
     #'EVM','E11','E22',
     #'E33','EPVM','EP11','EP22','EP33','DVM','D11','D22','D33','DPVM',
     #'DP11','DP22','DP33','SVM','S11','S22','S33','S23','S31','S12'
cg
      do 3000 imicro=1,nsteps

      write(*,*) '***************************************************'
      write(*,*) 'STEP = ',imicro
      if(nsteps.ne.1) write(21,*) 'STEP = ',imicro
      write(25,*) 'STEP = ',imicro
      sg23tot=0
      sg31tot=0
      sg12tot=0
ch
chh      IF(IWFIELDS.EQ.1) then
      IF(IWFIELDS.EQ.1.AND.imicro/iwstep*iwstep.eq.imicro) then
chh      IF(IWFIELDS.EQ.1.AND.
chh     #   (imicro-1)/iwstep*iwstep.eq.(imicro-1)) then
       write(91,*) 'STEP = ',imicro
       write(92,*) 'STEP = ',imicro
       write(93,*) 'STEP = ',imicro
caps
       write(94,*) 'STEP = ',imicro
caps
      ENDIF
cc
      do 7771 i=1,npts1
      do 7771 j=1,npts2
      do 7771 k=1,npts3
c
      do ii=1,3
      do jj=1,3
cg
cg    in velgrad we store disgrad at t
cg
      velgrad(ii,jj,i,j,k)=disgrad(ii,jj,i,j,k)
cg
      disgrad(ii,jj,i,j,k)=
     #   disgrad(ii,jj,i,j,k)+udot(ii,jj)*tdot
chard      disgrad(ii,jj,i,j,k)=disgradmacro(ii,jj)
      enddo
      enddo

7771  continue

      if(imicro.eq.1.or.iupdate.eq.1) call update_schmid
cc
      do ii=1,3
      do jj=1,3
cth
      ddisgradmacro(ii,jj)=0.
cth
      ddisgradmacroacum(ii,jj)=0.
      enddo
      enddo
cc
      iter=0

cc      err2mod=2*error
cc      do while(iter.lt.itmax.and.err2mod.gt.error)

      erre=2*error
      errs=2*error
      do while(iter.lt.itmax.and.(errs.gt.error.or.erre.gt.error))
cc      do while(iter.lt.itmax.and.erre.gt.error)

      iter=iter+1
c
      write(*,*)'ITER = ',iter
c
      write(25,*)'ITER = ',iter

cx       write(*,*) 'UPDATE STRESS FIELD'
cx       call evp(imicro,iter)
cx       call get_smacro

      write(*,*) 'DIRECT FFT OF STRESS FIELD'
c
      do 300 ii=1,3
      do 300 jj=1,3
c
      k1=0
      do 5 k=1,npts3
c
cw      write(*,'(1H+,a,i2,2(a,i4))')
cw     #   'STRESS - COMPONENT',ii,jj,'  -   Z = ',k,'   OUT OF ',npts
c
      do 5 j=1,npts2
      do 5 i=1,npts1
c
      k1=k1+1
      data(k1)=sg(ii,jj,i,j,k)
      k1=k1+1
      data(k1)=0.
5     continue
c
      if(npts3.gt.1) then
      call fourn(data,nn,3,1)
      else
      call fourn(data,nn2,2,1)
      endif
c
      k1=0
      do 6 kzz=1,npts3
      do 6 kyy=1,npts2
      do 6 kxx=1,npts1
c
      k1=k1+1
      work(ii,jj,kxx,kyy,kzz)=data(k1)
c
      k1=k1+1
      workim(ii,jj,kxx,kyy,kzz)=data(k1)
c
6     continue
c
300   continue

cx      write(*,*) 'Re(^SGij)(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (work(ii,jj,2,2,2),jj=1,3)
cx      enddo
cx      write(*,*) 'Im(^SGij)(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (workim(ii,jj,2,2,2),jj=1,3)
cx      enddo
cx      pause
cx

      write(*,*) 'CALCULATING G^pq,ij : SG^ij ...'
c
      do 1 kzz=1,npts3
      do 1 kyy=1,npts2
      do 1 kxx=1,npts1

      if(kxx.le.npts1/2) kx=kxx-1
      if(kxx.gt.npts1/2) kx=kxx-npts1-1
c
      if(kyy.le.npts2/2) ky=kyy-1
      if(kyy.gt.npts2/2) ky=kyy-npts2-1
c
      if(kzz.le.npts3/2) kz=kzz-1
      if(kzz.gt.npts3/2) kz=kzz-npts3-1
c
      xk(1)=kx/(delt(1)*nn(1))
      xk(2)=ky/(delt(2)*nn(2))
c
      if(npts3.gt.1) then
      xk(3)=kz/(delt(3)*nn(3))
      else
      xk(3)=0.
      endif
c
      xknorm=sqrt(xk(1)**2+xk(2)**2+xk(3)**2)
c
      if (xknorm.ne.0.) then
      do i=1,3
      xk2(i)=xk(i)/(xknorm*xknorm*2.*pi)
      xk(i)=xk(i)/xknorm
      enddo
      endif
ck
      if(kxx.eq.npts1/2+1.or.kyy.eq.npts2/2+1.or.
     #   (npts3.gt.1.and.kzz.eq.npts3/2+1)) then
ck      if(kxx.eq.npts1/2.or.kyy.eq.npts2/2.or.
ck     #   (npts3.gt.1.and.kzz.eq.npts3/2)) then
      g1=-s0
ck
      else
ck
      do 2 i=1,3
      do 2 k=1,3
      a(i,k)=0.
      do 2 j=1,3
      do 2 l=1,3
      a(i,k)=a(i,k)+c0(i,j,k,l)*xk(j)*xk(l)
2     continue

      call minv(a,3,det,minv1,minv2)
c
c      if(det.eq.0) then
c      write(*,*) kx,ky,kz,'  --> SINGULAR SYSTEM'
c      stop
c      pause
c      endif
c
      do 3 p=1,3
      do 3 q=1,3
      do 3 i=1,3
      do 3 j=1,3
      g1(p,q,i,j)=-a(p,i)*xk(q)*xk(j)
3     continue
ck
      endif
ck
c######################################################
ct      do 33 i=1,3
ct      do 33 j=1,3
ct      do 33 k=1,3
ct      g2(i,j,k)=-xk2(j)*a(k,i)
ct33     continue
c######################################################

cx      if(kxx.eq.2.and.kyy.eq.2.and.kzz.eq.2) then
cx      write(*,*) 'G1='
cx      write(*,*) g1
cx      write(*,*) 'WORK='
cx      do ii=1,3
cx      write(*,*) (work(ii,jj,kxx,kyy,kzz),jj=1,3)
cx      enddo
cx      write(*,*) 'WORKIM='
cx      do ii=1,3
cx      write(*,*) (workim(ii,jj,kxx,kyy,kzz),jj=1,3)
cx      enddo
cx      pause
cx      endif

      do 4 i=1,3
      do 4 j=1,3

cc      ddisgrad(i,j,kxx,kyy,kzz)=0.
cc      ddisgradim(i,j,kxx,kyy,kzz)=0.

      ddisgrad(i,j)=0.
      ddisgradim(i,j)=0.
c
      if(kx.eq.0.and.ky.eq.0.and.kz.eq.0.) goto 4
cx      if(kx.eq.0.and.ky.eq.0.and.kz.eq.0.) goto 44
c
      do k=1,3
      do l=1,3

cc      ddisgrad(i,j,kxx,kyy,kzz)=
cc     #    ddisgrad(i,j,kxx,kyy,kzz)+g1(i,j,k,l)*sg(k,l,kxx,kyy,kzz)
cc      ddisgradim(i,j,kxx,kyy,kzz)=
cc     #  ddisgradim(i,j,kxx,kyy,kzz)+g1(i,j,k,l)*sgim(k,l,kxx,kyy,kzz)

      ddisgrad(i,j)=
     #    ddisgrad(i,j)+g1(i,j,k,l)*work(k,l,kxx,kyy,kzz)
      ddisgradim(i,j)=
     #  ddisgradim(i,j)+g1(i,j,k,l)*workim(k,l,kxx,kyy,kzz)

cx      if(kxx.eq.2.and.kyy.eq.2.and.kzz.eq.2) then
cx       write(*,*) 'i,j,k,l=',i,j,k,l
cx       write(*,*) 'ddisgrad(i,j)'
cx       write(*,*) ddisgrad(i,j)
cx       write(*,*) 'g1(i,j,k,l)'
cx       write(*,*) g1(i,j,k,l)
cx       write(*,*) 'work(k,l,2,2,2)'
cx       write(*,*) work(k,l,kxx,kyy,kzz)
cx       write(*,*) 'workim(k,l,2,2,2)'
cx       write(*,*) workim(k,l,kxx,kyy,kzz)
cx       pause
cx      endif

      enddo
      enddo

cx44    work(i,j,kxx,kyy,kzz)=ddisgrad(i,j)
cx      workim(i,j,kxx,kyy,kzz)=ddisgradim(i,j)

4     continue

      do i=1,3
      do j=1,3
      work(i,j,kxx,kyy,kzz)=ddisgrad(i,j)
      workim(i,j,kxx,kyy,kzz)=ddisgradim(i,j)
      enddo
      enddo

1     continue

cx      write(*,*) 'Re(^DdisGRADij)(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (work(ii,jj,2,2,2),jj=1,3)
cx      enddo
cx      write(*,*) 'Im(^DdisGRADij)(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (workim(ii,jj,2,2,2),jj=1,3)
cx      enddo
cx      pause

cc      call equilibrium(snormfft,snormfftim,nn,err2mod)

      write(*,*) 'INVERSE FFT TO GET STRAIN FIELD'
c
      do 51 m=1,3
      do 51 n=1,3

cxx      ddisgradav(m,n)=0.

cc
cc      disgradmacro(m,n)=disgradmacro(m,n)+ddisgradmacro(m,n)
cc
      k1=0
c
      do 50 k=1,npts3
      do 50 j=1,npts2
      do 50 i=1,npts1
c
      k1=k1+1
      data(k1)=work(m,n,i,j,k)
c
      k1=k1+1
      data(k1)=workim(m,n,i,j,k)
c
50     continue
c
      if(npts3.gt.1) then
      call fourn(data,nn,3,-1)
      else
      call fourn(data,nn2,2,-1)
      endif
c
      do i=1,2*npts1*npts2*npts3
      data(i)=data(i)/prodnn
      enddo
c
      k1=0
      do 16 kzz=1,npts3
      do 16 kyy=1,npts2
      do 16 kxx=1,npts1
c
      k1=k1+1
c      write(*,*) 'REAL PART =',kxx,kyy,kzz,data(k1)
cc
cc      disgrad(m,n,kxx,kyy,kzz)=disgrad(m,n,kxx,kyy,kzz)+data(k1)
cc
      disgrad(m,n,kxx,kyy,kzz)=disgrad(m,n,kxx,kyy,kzz)+
     #  ddisgradmacro(m,n)+data(k1)

cxx      ddisgradav(m,n)=ddisgradav(m,n)+abs(data(k1))*wgt

      k1=k1+1
c      write(*,*) 'IMAGINARY PART =',kxx,kyy,kzz,data(k1)
c      pause
16    continue
c
51    continue
c
c     SYMMETRIZATION
c
cg      do 167 kzz=1,npts3
cg      do 167 kyy=1,npts2
cg      do 167 kxx=1,npts1
cgc
cg      do ii=1,3
cg      do jj=1,3
cg      aux33(ii,jj)=(disgrad(ii,jj,kxx,kyy,kzz)+
cg     #              disgrad(jj,ii,kxx,kyy,kzz))/2.
cg      enddo
cg      enddo
cgc
cg      do ii=1,3
cg      do jj=1,3
cg      disgrad(ii,jj,kxx,kyy,kzz)=aux33(ii,jj)
cg      enddo
cg      enddo
c
167   continue

cxx      erre2=tnorm(ddisgradav,3,3)

       write(*,*) 'UPDATE STRESS FIELD'
       call evpal(imicro)
       call get_smacro
c
      erre=erre/evm
cxx      erre2=erre2/evm
      errs=errs/svm

cc       write(*,*) 'STRESS FIELD ERROR =',errs/svm
cc       write(*,*) 'STRAIN FIELD ERROR =',erre/evm

       write(*,*) 'STRAIN FIELD ERROR =',erre
cxx       write(*,*) 'STRAIN FIELD ERROR (ALT) =',erre2
       write(*,*) 'STRESS FIELD ERROR =',errs

      write(21,101) iter,erre,errs,svm

c     STATISTIC ON ACTIVITIES
c
cc      call statactiv
c
cc      if(nph.eq.2) then
ccc
cc      if(igas(2).eq.1) then
cc      write(21,101) iter,errd/dvm,errs/svm,dvm,svm,
cc     # ((gavmod(k1,1)),k1=1,nmodes(1))
cc      else
cc      write(21,101) iter,errd/dvm,errs/svm,dvm,svm,
cc     # (((gavmod(k1,kph)),k1=1,nmodes(kph)),kph=1,2)
cc      endif
ccc
cc      else
ccc
cc      write(21,101) iter,errd/dvm,errs/svm,dvm,svm,
cc     # ((gavmod(k1,1)),k1=1,nmodes(1))
ccc
cc      endif
c
101   format(i3,4(1x,e10.4),10(1x,F7.4))
c
c     ENDDO ... WHILE
c
      enddo
cth
cth   kinematic hardening (ithermo.eq.1)
cth
      if(ithermo.eq.1.and.imicro.eq.1) call kinhard_param
cth
cg
cg    velgrad, which contained disgrad at t, is updated
cg
      velgrad=(disgrad-velgrad)/tdot
cg
cc      do ii=1,3
cc      do jj=1,3
cc      disgradmacro(ii,jj)=disgradmacro(ii,jj)+ddisgradmacro(ii,jj)
cc      enddo
cc      enddo

      do ii=1,3
      do jj=1,3
c
      disgradmacroactual(ii,jj)=
     #    disgradmacro(ii,jj)+ddisgradmacroacum(ii,jj)
c
cg      velgradmacro(ii,jj)=disgradmacroactual(ii,jj)-
cg     #                    disgradmacrot(ii,jj)
      velgradmacro(ii,jj)=(disgradmacroactual(ii,jj)-
     #                    disgradmacrot(ii,jj))/tdot
c
      enddo
      enddo

      write(*,*)
     #  'disgradmacro(1,1),disgradmacro(2,2),disgradmacro(3,3)'
      write(*,*) disgradmacroactual(1,1),disgradmacroactual(2,2),
     #            disgradmacroactual(3,3)
c
      write(*,*) 'disgradmacro(1,1)/disgradmacro(3,3)'
      write(*,*) disgradmacroactual(1,1)/
     #           disgradmacroactual(3,3)

cg      write(*,*) 'scauav(1,1),scauav(2,2),scauav(3,3)'
cg      write(*,*) scauav(1,1),scauav(2,2),scauav(3,3)
      write(*,*) 'scauav1(1,1),scauav1(2,2),scauav1(3,3)'
      write(*,*) scauav1(1,1),scauav1(2,2),scauav1(3,3)
c
c     TOTAL (EL+PL) VM
c
      evm=vm(disgradmacroactual)
cc      dvm=vm(disgradmacroactual-disgradmacrot)/tdot
cg      dvm=vm(velgradmacro)/tdot
      dvm=vm(velgradmacro)
      disgradmacrot=disgradmacroactual
c
c     INITIAL GUESS OF DISGRADMACRO AT t+dt ALWAYS ELASTIC
c
      do ii=1,3
      do jj=1,3
chard      disgradmacro(ii,jj)=disgradmacro(ii,jj)+udot(ii,jj)*tdot
      disgradmacro(ii,jj)=disgradmacrot(ii,jj)+udot(ii,jj)*tdot
      enddo
      enddo
c
      ept=ept+edotp*tdot
c
c     PLASTIC VM
c
      epav=0.
      edotpav=0.

      do 4777 i=1,npts1
      do 4777 j=1,npts2
      do 4777 k=1,npts3
       do ii=1,3
       do jj=1,3
        epav(ii,jj)=epav(ii,jj)+ept(ii,jj,i,j,k)*wgt
        edotpav(ii,jj)=edotpav(ii,jj)+edotp(ii,jj,i,j,k)*wgt
       enddo
       enddo
4777   continue

      evmp=0.
      dvmp=0.
      do ii=1,3
      do jj=1,3
       evmp=evmp+epav(ii,jj)**2
       dvmp=dvmp+edotpav(ii,jj)**2
      enddo
      enddo

      evmp=sqrt(2./3.*evmp)
      dvmp=sqrt(2./3.*dvmp)

cth      IF(IUPDATE.EQ.1) THEN
      if(iupdate.eq.1.and.(ithermo.ne.1.or.imicro.gt.1)) then
ccc
ccc
ccc     VELMAX
ccc
cc      velmax(1)=dsim(1,1)*delt(1)*(npts1-1)
cc      velmax(2)=dsim(2,2)*delt(2)*(npts2-1)
cc      velmax(3)=dsim(3,3)*delt(3)*(npts3-1)
ccc
ccc     UPDATE ORIENTATIONS
ccc
      call update_orient
ccc
ccc     UPDATE DELT
ccc
cc      delt(1)=(delt(1)*(npts1-1)+velmax(1)*tdot)/(npts1-1)
cc      delt(2)=(delt(2)*(npts2-1)+velmax(2)*tdot)/(npts2-1)
cc      if(npts3.gt.1) then
cc      delt(3)=(delt(3)*(npts3-1)+velmax(3)*tdot)/(npts3-1)
cc      endif
ccc
      ENDIF                !  IUPDATE ENDIF
c
cth      IF(IUPHARD.EQ.1) CALL HARDEN
      if(iuphard.eq.1.and.(ithermo.ne.1.or.imicro.gt.2))
cth      if(iuphard.eq.1.and.(ithermo.ne.1.or.imicro.gt.1))
     #  call harden
c
      write(55,315) evm,evmp,dvm,dvmp,svm,svm1
      do k=1,npts3
        do j=1,npts2
          do i=1,npts1
            sg23tot=sg23tot+sg(2,3,i,j,k)
            sg31tot=sg31tot+sg(3,1,i,j,k)
            sg12tot=sg12tot+sg(1,2,i,j,k)
          end do
        end do
      end do
c
      sg23avg=sg23tot/(npts1*npts2*npts3)
      write(*,*) sg23avg
      sg31avg=sg31tot/(npts1*npts2*npts3)
      write(*,*) sg31avg
      sg12avg=sg12tot/(npts1*npts2*npts3)
      write(*,*) sg12avg
      write(56,315) evm,disgradmacroactual(1,1),
     #   disgradmacroactual(2,2),disgradmacroactual(3,3),
     #   evmp,epav(1,1),epav(2,2),epav(3,3),
     #   dvm,velgradmacro(1,1),velgradmacro(2,2),velgradmacro(3,3),
     #   dvmp,edotpav(1,1),edotpav(2,2),edotpav(3,3),
     #   svm,scauav(1,1),scauav(2,2),scauav(3,3),
     #   sg23avg,sg31avg,sg12avg
cg

cg
315   format(23(e11.4,1x))
c
chh      IF(IWFIELDS.EQ.1) THEN
      IF(IWFIELDS.EQ.1.AND.imicro/iwstep*iwstep.eq.imicro) then
chh      IF(IWFIELDS.EQ.1.AND.
chh     #   (imicro-1)/iwstep*iwstep.eq.(imicro-1)) then
ch
ch      IF(IWFIELDS.EQ.1.AND.
ch     #  (IMICRO.EQ.1.OR.IMICRO.EQ.40.OR.IMICRO.EQ.100) ) THEN
ch
c
c      FIELDS.OUT
c
ch      IF(IMICRO.EQ.1) THEN
ch      open(91,file='fields1.out',status='unknown')
ch      ELSE IF(IMICRO.EQ.40) THEN
ch      open(91,file='fields40.out',status='unknown')
ch      ELSE IF(IMICRO.EQ.100) THEN
ch      open(91,file='fields100.out',status='unknown')
ch      ENDIF
ch
      write(91,*) delt
      write(92,*) delt
      write(93,*) delt
caps
      write(94,*) delt
caps
c
      write(91,*)
     # '   x    y    z  ngr   ph  ELOC (11,22,33,23,31,12)'
      write(92,*)
     # '   x    y    z  ngr   ph  SLOC (11,22,33,23,31,12)'
      write(93,*)
     # '   x    y    z  ngr   ph  DPLOC (11,22,33,23,31,12)'
caps
      write(94,*)
     # '   x    y    z  ngr   ph  ELLOC (11,22,33,23,31,12)'
caps

      do k=1,npts3
      do j=1,npts2
      do i=1,npts1
caps
      jph=jphase(i,j,k)
      if(igas(jph).eq.0) then
cg
       do ii=1,6
       do jj=1,6
        cg66aux(ii,jj)=cg66(ii,jj,i,j,k)
       enddo
       enddo
c
      call lu_inverse(cg66aux,6)
      call chg_basis(aux6,aux33,cg66aux,cinvg,3,6)
c
       do ii=1,3
       do jj=1,3
       eel(ii,jj)=0.
       do kk=1,3
       do ll=1,3
       eel(ii,jj)=eel(ii,jj)+
     #                cinvg(ii,jj,kk,ll)*sg(kk,ll,i,j,k)
       enddo
       enddo
       enddo
       enddo
cg
       else     !   igas else
cg
       do ii=1,3
       do jj=1,3
caps       eel(ii,jj)=0.
       eel(ii,jj)=-100.
       enddo
       enddo
cg
       endif    ! igas endif
caps
      write(91,305) i,j,k,jgrain(i,j,k),jphase(i,j,k),
     # disgrad(1,1,i,j,k),disgrad(2,2,i,j,k),disgrad(3,3,i,j,k),
     # disgrad(2,3,i,j,k),disgrad(3,1,i,j,k),disgrad(1,2,i,j,k)

      write(92,305) i,j,k,jgrain(i,j,k),jphase(i,j,k),
     # sg(1,1,i,j,k),sg(2,2,i,j,k),sg(3,3,i,j,k),
     # sg(2,3,i,j,k),sg(3,1,i,j,k),sg(1,2,i,j,k)

caps
      write(93,305) i,j,k,jgrain(i,j,k),jphase(i,j,k),
     # edotp(1,1,i,j,k),edotp(2,2,i,j,k),edotp(3,3,i,j,k),
     # edotp(2,3,i,j,k),edotp(3,1,i,j,k),edotp(1,2,i,j,k)
caps
      write(94,305) i,j,k,jgrain(i,j,k),jphase(i,j,k),
     # eel(1,1),eel(2,2),eel(3,3),eel(2,3),eel(3,1),eel(1,2)


caps
305   format(5i7,18(e11.3,1x))
c
      enddo
      enddo
      enddo


c
ch      close(91)

      ENDIF
c
3000  CONTINUE
cw
      IF(IWTEX.EQ.1) THEN
c
c     TEX.OUT
c
      open(24,file='tex.out',status='unknown')
c
      one=1.
c
      do i=1,3
      do j=1,3
      fbar(i,j)=0.
      enddo
      enddo
c
      do i=1,3
      fbar(i,i)=delt(i)
      enddo
c
      write(24,*)'TEX.OUT from FFT'
      write(24,*)'Formatted to plot with POLE'
      write(24,111)((fbar(i,j),j=1,3),i=1,3)
      write(24,'(a1,i10)') 'B',npts1*npts2*npts3
111   format(9f7.3)
c
      ig=0
      do k=1,npts3
      do j=1,npts2
      do i=1,npts1
      ig=ig+1
c
      do ii=1,3
      do jj=1,3
      aux33(ii,jj)=ag(jj,ii,i,j,k)
      enddo
      enddo
c
      call euler(1,ph,th,om,aux33)
c
       write(24,3305) ph,th,om,one,i,j,k,jgrain(i,j,k),jphase(i,j,k)
3305  format(4(f9.3,1x),5i5)
c
      enddo
      enddo
      enddo
c
      ENDIF
c
ccc###########################################################
cc
ccc      STRESS.OUT
ccc
cc      if(isave.eq.1) then
ccc
cc      open(40,file='stress.out',status='unknown',
cc     #     access='sequential',form='unformatted')
ccc
cc      do i=1,5
cc      write(40) (xlsec(i,j),j=1,5)
cc      enddo
ccc
cc      do kk=1,npts3
cc      do jj=1,npts2
cc      do ii=1,npts1
ccc
cc      write(40) (sg(i,ii,jj,kk),i=1,5)
ccc
cc      enddo
cc      enddo
cc      enddo
ccc
cc      endif
ccc
      end
c
      SUBROUTINE fourn(data,nn,ndim,isign)
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then

            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif

          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr

                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue
      return
      END
c
      SUBROUTINE VOIGT(C2,C4,IOPT)
C
C *** TRANSFORMS SECOND ORDER MATRIX C2 INTO FOURTH ORDER TENSOR C4 IF
C *** IOPT=1 AND VICEVERSA IF IOPT=2. IF IOPT=3,TRANSFORMS WITH INV.FACT.
C *** IOPT=4 FOR GO FROM 6x6 TO 3x3x3x3 WITH Aijkl ANTISYMMETRY
C
      DIMENSION C2(6,6),C4(3,3,3,3),IJV(6,2),F(6,6)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
C
      IF(IOPT.EQ.1) THEN
      DO 10 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 10 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
   10 C4(I2,I1,J2,J1)=C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.2) THEN
      DO 20 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 20 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
   20 C2(I,J)=C4(I1,I2,J1,J2)
      ENDIF
c
      IF(IOPT.EQ.3) THEN
      DO 9 I=1,6
      DO 9 J=1,6
      F(I,J)=1.
      IF(I.GT.3) F(I,J)=0.5
      IF(J.GT.3) F(I,J)=0.5*F(I,J)
9     CONTINUE
C
      DO 101 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 101 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      C4(I1,I2,J1,J2)=F(I,J)*C2(I,J)
      C4(I2,I1,J1,J2)=F(I,J)*C2(I,J)
      C4(I1,I2,J2,J1)=F(I,J)*C2(I,J)
101   C4(I2,I1,J2,J1)=F(I,J)*C2(I,J)
      ENDIF
C
      IF(IOPT.EQ.4) THEN
      DO 17 I=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      DO 17 J=1,6
      J1=IJV(J,1)
      J2=IJV(J,2)
      IF(I.LE.3) THEN
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
      C4(I2,I1,J2,J1)=C2(I,J)
      ELSE
      C4(I1,I2,J1,J2)=C2(I,J)
      C4(I2,I1,J1,J2)=-C2(I,J)
      C4(I1,I2,J2,J1)=C2(I,J)
      C4(I2,I1,J2,J1)=-C2(I,J)
      ENDIF
17    CONTINUE
      ENDIF
      RETURN
      END
c
      subroutine euler(iopt,ph,th,tm,a)
      dimension a(3,3)
      pi=4.*atan(1.d0)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES OF ca REFERRED TO sa.
c
      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph)
        cph=cos(ph)
        sth=sin(th)
        cth=cos(th)
        stm=sin(tm)
        ctm=cos(tm)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end
c
C ********************************************************************
C     SUBROUTINE VPSC_INPUT      --->      VERSION 31/jan/99
C
C     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
C     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
C     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
C     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
C     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
C     OPENS AND CLOSES INPUT FILES.   OPENS OUTPUT FILES.
C
C     MODIFIED 21/07/98 by CNT:
C     INITIALIZATION RELATED TO 'ELEMENTS' IS DONE INSIDE A SINGLE BLOCK.
C *****************************************************************************
C
      SUBROUTINE VPSC_INPUT

      INCLUDE 'fft.dim'

      dimension aux55(5,5),aux66(6,6),aux3333(3,3,3,3)
      dimension aux6(6),aux33(3,3)
      dimension dsimdev(3,3)

      DIMENSION IJV(6,2)
      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/

C *********   INITIALIZATION BLOCK   ***************************
C
      PI=4.*ATAN(1.)
C
cc      IREADSET=0      ! used as control in CUBCOMP to open unit
C
C     CALCULATES TENSORS OF THE SYMMETRIC BASIS 'B(3,3,6)'
cth      CALL CHG_BASIS(DUM1,DUM2,DUM3,DUM4,0,6)
      CALL CHG_BASIS(AUX6,AUX33,AUX66,AUX3333,0,6)
C
C     SEED FOR RANDOM NUMBER GENERATOR (RAN2) (USED FOR TWINNING AND RX)
      JRAN=-1
C
      READ(UR0,*) nph

C     THE FOLLOWING REQUIRED FOR SEVERAL ROUTINES WITH 'do iph=iphbot,iphtop'
        iphbot=1
        iphtop=nph
c
      if(nph.gt.nphmx) then
        write(*,'('' number of phases exceeds multiphase dimens !!'')')
        write(*,'('' --> increase parameter NPHMX to'',i5)') nph
        stop
      endif
c
      READ(UR0,*) ngr
c
      if(ngr.ne.npts1*npts2*npts3) then
      write(*,*) 'NUMBER OF FPs DECLARED IN FFT.IN = ',ngr
      write(*,*) 'DOES NOT MATCH WITH NPTS1*NPTS2*NPTS3 IN FFT.DIM =',
     #           npts1*npts2*npts3
      stop
      endif
c
c     RVE DIMENSIONS
c
      READ(UR0,*) DELT
      DELTVOL3=(DELT(1)*DELT(2)*DELT(3))**(1./3.)
c
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filetext

c ***************************************************************************
c     LOOP OVER PHASES
c ***************************************************************************

      DO IPH=1,NPH
c
      READ(UR0,'(a)') prosa
      READ(UR0,*) igas(iph)
      READ(UR0,'(a)') prosa
      READ(UR0,'(a)') filecryspl
cevp
      READ(UR0,'(a)') filecrysel
cevp
c
c     READS SLIP AND TWINNING MODES FOR THE PHASE
c
      if(igas(iph).eq.0) then
      OPEN (unit=UR1,file=filecryspl,status='old')
        call data_crystal(iph)
      CLOSE(unit=UR1)
cevp
      OPEN (unit=UR1,file=filecrysel,status='old')
        call data_crystal_elast(iph)
      CLOSE(unit=UR1)
cevp
      endif
c
      ENDDO     ! END OF DATA INPUT LOOP OVER ALL PHASES
c
c     READS INITIAL TEXTURE FROM FILETEXT
cevp  and calculates the local elastic stiffness
c
      OPEN(unit=UR2,file=filetext,status='old')
        call data_grain
      CLOSE(unit=UR2)
ccc
ccc     SEARCH FOR NRSMIN (NEEDED TO GET TAUMAX INSIDE NR SUBROUTINE)
ccc
cc      nrsmin=nrs(1,1)
cc      DO IPH=1,NPH
cc        do is=1,nsyst(iph)
cc          if(nrs(is,iph).lt.nrsmin) nrsmin=nrs(is,iph)
cc        enddo
cc      ENDDO
c
C ****************************************************************************
C     READ BOUNDARY CONDITIONS ON OVERALL STRESS AND STRAIN-RATE
C ****************************************************************************
C
      READ(UR0,'(A)') PROSA
      READ(UR0,'(A)') PROSA
C
      do i=1,3
        READ(UR0,*) (iudot(i,j),j=1,3)
      enddo
c
      if(iudot(1,1)+iudot(2,2)+iudot(3,3).eq.2) then
        write(*,*) 'CHECK DIAGONAL BOUNDARY CONDITIONS IUDOT'
        write(*,*) 'CANNOT ENFORCE ONLY TWO DEVIATORIC COMPONENTS'
        stop
      endif
c
      do i=1,3
      do j=1,3
        if(i.ne.j.and.iudot(i,j)+iudot(j,i).eq.0) then
          write(*,*) 'CHECK OFF-DIAGONAL BOUNDARY CONDITIONS IUDOT'
          stop
        endif
      enddo
      enddo
c
      READ(UR0,*)
      DO I=1,3
        READ(UR0,*) (UDOT(I,J),J=1,3)
      ENDDO
C
c     SYMMETRIC STRAIN-RATE, ANTISYMMETRIC ROTATION-RATE TENSORS
c     AND INDICES OF IMPOSED COMPONENTS
c
      do i=1,3
      do j=1,3
        dsim(i,j)=(udot(i,j)+udot(j,i))/2.
        tomtot(i,j)=(udot(i,j)-udot(j,i))/2.
      enddo
      enddo
c
      do i=1,3
        idsim(i)=iudot(i,i)
      enddo
c
      idsim(4)=0
      if(iudot(2,3).eq.1.and.iudot(3,2).eq.1) idsim(4)=1
      idsim(5)=0
      if(iudot(1,3).eq.1.and.iudot(3,1).eq.1) idsim(5)=1
      idsim(6)=0
      if(iudot(1,2).eq.1.and.iudot(2,1).eq.1) idsim(6)=1
c
c     WRITES STRAIN RATE DSIM(I,J) IN b-BASIS AS A 5-DIM VECTOR DBAR(K)
c
      call chg_basis(dbar6,dsim,aux66,aux3333,2,6)

      do i=1,5
      dbar5(i)=dbar6(i)
      enddo

      call chg_basis(dbar5,dsimdev,aux55,aux3333,1,5)
c
      dvm=0.
      do i=1,3
      do j=1,3
      dvm=dvm+dsimdev(i,j)**2
      enddo
      enddo
      dvm=sqrt(2./3.*dvm)
cw
cw      evm=0.
cw
      READ(UR0,*)
      READ(UR0,*) iscau(1),iscau(6),iscau(5)
      READ(UR0,*) iscau(2),iscau(4)
      READ(UR0,*) iscau(3)
c
      do i=1,6
        if(iscau(i)*idsim(i).ne.0.or.iscau(i)+idsim(i).ne.1) then
          WRITE(*,*) ' CHECK BOUNDARY CONDITS ON STRAIN-RATE AND STRESS'
          WRITE(*,'('' IDSIM = '',6I3)') IDSIM
          WRITE(*,'('' ISCAU = '',6I3)') ISCAU
          STOP
        endif
      enddo
c
      READ(UR0,*)
      READ(UR0,*) scauchy(1,1),scauchy(1,2),scauchy(1,3)
      READ(UR0,*) scauchy(2,2),scauchy(2,3)
      READ(UR0,*) scauchy(3,3)
c
      scauchy(3,2)=scauchy(2,3)
      scauchy(3,1)=scauchy(1,3)
      scauchy(2,1)=scauchy(1,2)
c
c     IF SOME OFF-DIAG COMPS OF SCAUCHY ARE IMPOSED
c            WRITE THE CORRESP COMP OF SBAR
c
cc      if(iscau(4).eq.1) sbar(3)=sqrt(2.d0)*scauchy(2,3)
cc      if(iscau(5).eq.1) sbar(4)=sqrt(2.d0)*scauchy(1,3)
cc      if(iscau(6).eq.1) sbar(5)=sqrt(2.d0)*scauchy(1,2)
c
      READ(UR0,*)
      READ(UR0,*) DUMMY
      READ(UR0,*) ICTRL
c
      if(ictrl.eq.-1) tdot=dummy
      if(ictrl.eq.0)  tdot=dummy/dvm
      if(ictrl.gt.0)  then
      write(*,*) 'ICTRL>0 NOT IMPLEMENTED'
      stop
      endif
c
c      if(ictrl.gt.0)  eijincr=dummy
c
c      if(ictrl.gt.0) then
c        if(.not.(idsim(ictrl).eq.1.and.
c     #      dsim(ijv(ictrl,1),ijv(ictrl,2)).ne.0.)) then
c         write(*,*) 'ICTRL        =',ictrl
c         write(*,*) 'IDSIM(ICTRL) =',idsim(ictrl)
c         write(*,*) 'DSIM(ICRTL)  =',dsim(ijv(ictrl,1),ijv(ictrl,2))
c         write(*,*) 'ICTRL MUST BE -1 TO IMPOSE DSIM*TDOT'
c         write(*,*) 'OR IT MUST BE 0 TO CONTROL VON MISES INCREMENT'
c         write(*,*) 'OR IT MUST IDENTIFY A NON-ZERO STRAIN COMPONENT'
c         stop
c        endif
c      endif
c
      READ(UR0,*)
      READ(UR0,*) NSTEPS
      READ(UR0,*) ERROR
      READ(UR0,*) ITMAX
      READ(UR0,*) IRECOVER

      if(irecover.eq.1) open(50,file='stress.in',status='old',
     #     access='sequential',form='unformatted')
c
      READ(UR0,*) ISAVE
      READ(UR0,*) IUPDATE
c
cw      if(iupdate.eq.1) open(26,file='update.out',status='unknown')
c
      READ(UR0,*) IUPHARD
      READ(UR0,*) IWTEX
chh      READ(UR0,*) IWFIELDS
      READ(UR0,*) IWFIELDS,IWSTEP
cw      READ(UR0,*) FACT2
cw      FACT2=1.
cw
c
c     INITIALIZE CRSS AND ACCUM SHEAR FOR GRAINS
c
c      do KKK=ngr(iph-1)+1,ngr(iph)
c
c        gtotgr(kkk)=0.
c
cw      DO IPH=1,NPH
c
        do ii=1,npts1
        do jj=1,npts2
        do kk=1,npts3
c
        gacumgr(ii,jj,kk)=0.
c
        iph=jphase(ii,jj,kk)
cth
        if(igas(iph).eq.0) then
cth
         do i=1,nsyst(iph)
c
cw        fact=1.
cw        if(jphase(ii,jj,kk).eq.2) fact=fact2
cw
cw          crss(i,1,ii,jj,kk)=fact*tau(I,1,iph)
cw          crss(i,2,ii,jj,kk)=fact*tau(I,2,iph)
cw
cw          crss(i,1)=fact*tau(I,1,iph)
cw          crss(i,2)=fact*tau(I,2,iph)

          crss(i,1,ii,jj,kk)=tau(I,1,iph)
          crss(i,2,ii,jj,kk)=tau(I,2,iph)
chard
          TRIALTAU(I,1,II,JJ,KK)=tau(I,1,iph)
          TRIALTAU(I,2,II,JJ,KK)=tau(I,2,iph)
chard
cth
          XKIN(I,II,JJ,KK)=0.
cth
c
         enddo
cth
        endif
cth
      enddo
      enddo
      enddo
c
cw      ENDDO
cth
      READ(UR0,*) ITHERMO
      if(ithermo.eq.1) then
c
      READ(UR0,'(a)') FILETHERMO
cx       open(49,file='thermo.in',status='old')
       open(49,file=filethermo,status='old')
c
       do ig=1,1118
        read(49,*) eth(1,1,ig),eth(2,2,ig),
     #             eth(3,3,ig),eth(2,3,ig),
     #             eth(3,1,ig),eth(1,2,ig)
cx     #             ,eini(1,1,ig),eini(2,2,ig),
cx     #             eini(3,3,ig),eini(2,3,ig),
cx     #             eini(3,1,ig),eini(1,2,ig)

        eth(3,2,ig)=eth(2,3,ig)
        eth(1,3,ig)=eth(3,1,ig)
        eth(2,1,ig)=eth(1,2,ig)
cx        eini(3,2,ig)=eini(2,3,ig)
cx        eini(1,3,ig)=eini(3,1,ig)
cx        eini(2,1,ig)=eini(1,2,ig)
       enddo
cthh
cx      eth=eth*0.1
cx      eini=eini*0.1
cthh
      else
       eth=0.
      endif
cth
      RETURN
      END
c
C *****************************************************************************
C     SUBROUTINE DATA_CRYSTAL        --->      VERSION 03/FEB/2000
C *****************************************************************************

      SUBROUTINE DATA_CRYSTAL(IPH)

      INCLUDE 'fft.dim'

      DIMENSION ISN(12,4),ISB(12,4),SN(3),SB(3),CDIM(3)
      DIMENSION aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)
      DIMENSION HSELFX(10),HLATEX(10,10),MODE(10)
c
      READ(UR1,'(a)') prosa
      READ(UR1,'(a)') icryst(iph)
      READ(UR1,*)     (cdim(i),i=1,3)
      covera=cdim(3)/cdim(1)
      READ(UR1,*)     nmodesx
      READ(UR1,*)     nmodes(iph)
      READ(UR1,*)     (mode(i),i=1,nmodes(iph))
C
      IF(NMODES(IPH).GT.NMODMX) THEN
        WRITE(*,'('' NMODES IN PHASE'',I3,'' IS'',I3)') IPH,NMODES(IPH)
        WRITE(*,'('' CHANGE PARAMETER NMODMX IN fft.dim'')')
        STOP
      ENDIF
c
      NTWMOD(iph)=0
      NSYST(iph) =0
      NTWSYS(iph)=0
      KOUNT=1
c
c     START READING DEFORMATION MODES FROM FILECRYS
c
      do 100 nm=1,nmodesx
c
        READ(UR1,'(a)') prosa
        READ(UR1,*)     modex,nsmx,nrsx,gamd0x,twshx,isectwx
        READ(UR1,*)     tau0xf,tau0xb,tau1x,thet0x,thet1x
        READ(UR1,*)     hselfx(nm),(hlatex(nm,jm),jm=1,nmodesx)
c
c     SKIPS nsmx LINES IF THE MODE IS NOT IN THE LIST.
        if(modex.ne.mode(kount)) then
          do iz=1,nsmx
            READ(UR1,*)
          enddo
          go to 100
        endif
C
        IF(THET0X.LT.THET1X) THEN
          WRITE(*,'('' INITIAL HARDENING LOWER THAN FINAL HARDENING FOR
     #      MODE'',I3,''  IN PHASE'',I3)') KOUNT,IPH
          STOP
        ENDIF
C
C     CASE TAU1=0 CORRESPONDS TO LINEAR HARDENING AND IS INDEPENDENT OF TAU0.
C     AVOID DIVISION BY ZERO
        IF(TAU1X.LE.1.E-6) THEN
          TAU1X=1.E-6
          THET0X=THET1X
        ENDIF
c
c     REORDER HARDENING COEFFICIENTS TO ACCOUNT ONLY FOR ACTIVE MODES
        hselfx(kount)=hselfx(nm)
        do i=1,nmodes(iph)
          hlatex(kount,i)=hlatex(nm,mode(i))
        enddo
c
c     SYSTEMS GIVEN IN FOUR INDEX NOTATION: HEXAGONALS AND TRIGONALS
c     SYSTEMS GIVEN IN THREE INDEX NOTATION: CUBIC AND ORTHORHOMBIC
c
      IF(icryst(iph).EQ.'HEX' .OR. icryst(iph).EQ.'TRI' .OR.
     #   icryst(iph).EQ.'hex' .OR. icryst(iph).EQ.'tri') THEN
        DO J=1,NSMX
          READ(UR1,*) (ISN(J,K),K=1,4),(ISB(J,K),K=1,4)
        ENDDO
      ELSE IF(icryst(iph).EQ.'CUB' .OR. icryst(iph).EQ.'ORT' .OR.
     #        icryst(iph).EQ.'cub' .OR. icryst(iph).EQ.'ort') THEN
        DO J=1,NSMX
          READ(UR1,*)(ISN(J,K),K=1,3),(ISB(J,K),K=1,3)
        ENDDO
      ELSE
        WRITE(*,'('' CANNOT IDENTIFY THE CRYSTAL SYMMETRY OF PHASE '',
     #            I3)') IPH
        STOP
      ENDIF
C
      NSM(kount,iph)=nsmx
      IF(TWSHX.NE.0) NTWMOD(iph)=NTWMOD(iph)+1
C
      IF(NTWMOD(IPH).GT.NTWMMX) THEN
        WRITE(*,'('' NTWMOD IN PHASE'',I3,'' IS'',I3)') IPH,NTWMOD(IPH)
        WRITE(*,'('' CHANGE PARAMETER NTWMMX IN VPSC.DIM'')')
        STOP
      ENDIF
C
      DO 205 JS=1,NSM(kount,iph)
C
      NSYST(iph)=NSYST(iph)+1
      NSYSX=NSYST(iph)
      IF(TWSHX.NE.0) NTWSYS(iph)=NTWSYS(iph)+1
C
      IF(NSYST(IPH).GT.NSYSMX) THEN
        WRITE(*,'('' NSYST IN PHASE'',I3,'' IS'',I3)') IPH,NSYST(IPH)
        WRITE(*,'('' CHANGE PARAMETER NSYSMX IN VPSC.DIM'')')
        STOP
      ENDIF
C
C   DEFINES RATE SENSITIVITY AND CRSS FOR EACH SYSTEM IN THE MODE
C
      GAMD0(NSYSX,iph) =GAMD0X
      NRS(NSYSX,iph)   =NRSX
      TWSH(NSYSX,iph)  =TWSHX
      TAU(NSYSX,1,iph) =TAU0XF
      TAU(NSYSX,2,iph) =TAU0XB
      TAU(NSYSX,3,iph) =TAU1X
      THET(NSYSX,0,iph)=THET0X
      THET(NSYSX,1,iph)=THET1X
c
      isectw(NSYSX,iph)=isectwx
C
      IF(icryst(iph).EQ.'HEX' .OR. icryst(iph).EQ.'TRI' .OR.
     #   icryst(iph).EQ.'hex' .OR. icryst(iph).EQ.'tri') THEN
        SN(1)= ISN(JS,1)
        SN(2)=(ISN(JS,1)+2.*ISN(JS,2))/SQRT(3.)
        SN(3)= ISN(JS,4)/COVERA
        SB(1)= 3./2.*ISB(JS,1)
        SB(2)=(ISB(JS,1)/2.+ISB(JS,2))*SQRT(3.)
        SB(3)= ISB(JS,4)*COVERA
      ELSE IF(icryst(iph).EQ.'CUB' .OR. icryst(iph).EQ.'ORT' .OR.
     #        icryst(iph).EQ.'cub' .OR. icryst(iph).EQ.'ort') THEN
        DO M=1,3
          SN(M)=ISN(JS,M)/CDIM(M)
          SB(M)=ISB(JS,M)*CDIM(M)
        ENDDO
      ENDIF
C
C *** NORMALIZES SYSTEM VECTORS AND CHECKS NORMALITY
C
      SNOR=SQRT(SN(1)*SN(1)+SN(2)*SN(2)+SN(3)*SN(3))
      QNOR=SQRT(SB(1)*SB(1)+SB(2)*SB(2)+SB(3)*SB(3))
      PROD=0.
      DO J=1,3
        DNCA(J,NSYSX,iph)=SN(J)/SNOR
        DBCA(J,NSYSX,iph)=SB(J)/QNOR
        IF(ABS(DNCA(J,NSYSX,iph)).LT.1.E-03) DNCA(J,NSYSX,iph)=0.
        IF(ABS(DBCA(J,NSYSX,iph)).LT.1.E-03) DBCA(J,NSYSX,iph)=0.
        PROD=PROD+DNCA(J,NSYSX,IPH)*DBCA(J,NSYSX,IPH)
      ENDDO
      IF(PROD.GE.1.E-3) THEN
        WRITE(*,'(''SYSTEM'',I4,''  IN MODE'',I4,'' IN PHASE'',I4,
     #       ''  IS NOT ORTHOGONAL !!'')') JS,NM,IPH
        STOP
      ENDIF
C
C   DEFINE SCHMID VECTOR IN CRYSTAL AXES FOR EACH SYSTEM
C
      DO I=1,3
      DO J=1,3
        AUX33(i,j)=(DNCA(i,NSYSX,IPH)*DBCA(j,NSYSX,IPH)+
     #              DNCA(j,NSYSX,IPH)*DBCA(i,NSYSX,IPH))/2.
      ENDDO
      ENDDO
c
      call chg_basis(aux5,aux33,aux55,aux3333,2,5)
c
      DO I=1,5
        SCHCA(I,NSYSX,IPH)=AUX5(I)
      ENDDO

  205 CONTINUE    ! END OF LOOP OVER A GIVEN DEFORMATION MODE

      kount=kount+1

  100 CONTINUE    ! END OF LOOP OVER ALL MODES IN A GIVEN PHASE

C     INITIALIZE SELF & LATENT HARDENING COEFS FOR EACH SYSTEM OF THE PHASE.
C     ABSOLUTE UNITS ARE ACCOUNTED FOR BY MODULATING FACTOR IN HARDENING LAW.

      I=0
      DO IM=1,NMODES(iph)
      DO IS=1,NSM(IM,iph)
        I=I+1
        J=0
        DO JM=1,NMODES(iph)
        DO JS=1,NSM(JM,iph)
          J=J+1
          HARD(I,J,IPH)=HLATEX(IM,JM)
        ENDDO
        ENDDO
        HARD(I,I,IPH)=HSELFX(IM)
      ENDDO
      ENDDO

cc      do i=1,2
cc      write(*,*) (hard(i,j,iph),j=1,2)
cc      enddo
cc      pause

C     LATENT HARDENING OF SLIP AND TWINNING BY TWINNING IS BASED ON THE
C     RELATIVE DIRECTIONS OF THE SHEAR DIRECTION AND THE TWIN PLANE
C     THIS APPROACH IS STILL BEING TESTED (30/4/99)

C     NSLSYS=NSYST(IPH)-NTWSYS(IPH)
C     DO IS=1,NSYST(IPH)
C       DO JT=NSLSYS+1,NSYST(IPH)
C         IF(IS.NE.JT) THEN
C           COSA=DBCA(1,IS,IPH)*DNCA(1,JT,IPH)+
C    #           DBCA(2,IS,IPH)*DNCA(2,JT,IPH)+
C    #           DBCA(3,IS,IPH)*DNCA(3,JT,IPH)
C           COSA=ABS(COSA)
C           HARD(IS,JT,IPH)=HARD(IS,JT,IPH)*(0.5+1.0*COSA)
C         ENDIF
C       ENDDO
C     ENDDO
C
C     WRITE(10,'(''  HARDENING MATRIX FOR PHASE'',I3)') IPH
C     DO I=1,NSYST(IPH)
C       WRITE(10,'(24F5.1)') (HARD(I,J,IPH),J=1,NSYST(IPH))
C     ENDDO

C     VERIFICATION OF TWINNING DATA TO BE SURE PROGRAM WILL RUN PROPERLY

      IF (NMODES(IPH). GT. 1) THEN
        DO I=2,NSYST(IPH)
          IF(TWSH(I,IPH).EQ.0. .AND. TWSH(I-1,IPH).NE.0.) THEN
            WRITE(*,*) ' WARNING! THE TWINNING MODES MUST FOLLOW THE'
            WRITE(*,*) ' SLIP MODES   -->   REORDER CRYSTAL FILE'
            STOP
          ENDIF
        ENDDO
      ENDIF
C
      RETURN
      END
C
C *****************************************************************************
C     SUBROUTINE DATA_GRAIN        --->      VERSION 31/mar/99
C *****************************************************************************

      SUBROUTINE DATA_GRAIN

      INCLUDE 'fft.dim'

      DIMENSION AA(3,3)

cx      dimension caux3333(3,3,3,3),caux66(6,6),c066(6,6),s066(6,6)
      dimension caux3333(3,3,3,3),caux66(6,6),s066(6,6)
cx
      dimension aux6(6),aux33(3,3)
cx
      do i=1,6
      do j=1,6
       c066(i,j)=0.
      enddo
      enddo
cg
      nph1=0
cg
cth      do kkk=1,ngr
      DO 188 KKK=1,NGR

      READ(UR2,*) ph,th,om,ii,jj,kk,jgr,jph
cg
      if(jph.eq.1) nph1=nph1+1
cg
      jgrain(ii,jj,kk)=jgr
      jphase(ii,jj,kk)=jph
cth
      if(igas(jph).eq.1) goto 188
cth
C     CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
C     SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.

cbug        CALL EULER(2,ph,th,om,aa)
        CALL EULER(2,ph*pi/180,th*pi/180,om*pi/180,aa)
c
        DO J=1,3
        DO K=1,3
          AG(J,K,ii,jj,kk)=AA(K,J)
        ENDDO
        ENDDO
cw
cw        DO J=1,3
cw        DO K=1,3
cw          AG(J,K,jgrain(ii,jj,kk))=AA(K,J)
cw        ENDDO
cw        ENDDO
cw

      do 1 i1=1,3
      do 1 j1=1,3
      do 1 k1=1,3
      do 1 l1=1,3
      dum=0.
      do 2 i2=1,3
      do 2 j2=1,3
      do 2 k2=1,3
      do 2 l2=1,3
      dum=dum+
cg     #    aa(i2,i1)*aa(j2,j1)*aa(k2,k1)*aa(l2,l1)*cc(i2,j2,k2,l2)
     #    aa(i2,i1)*aa(j2,j1)*aa(k2,k1)*aa(l2,l1)*cc(i2,j2,k2,l2,jph)
2     continue
      caux3333(i1,j1,k1,l1)=dum
1     continue
c
      call chg_basis(aux6,aux33,caux66,caux3333,4,6)
c
      do i=1,6
      do j=1,6
      cg66(i,j,ii,jj,kk)=caux66(i,j)
      c066(i,j)=c066(i,j)+caux66(i,j)*wgt
      enddo
      enddo

cth      ENDDO
188   CONTINUE
cg
      wph1=(1.*nph1)/(npts1*npts2*npts3)
cg
cx      write(*,*) 'c066='
cx      do i=1,6
cx      write(*,*)(c066(i,j),j=1,6)
cx      enddo
cx      pause
c
      s066=c066
      call lu_inverse(s066,6)
C
      call chg_basis(aux6,aux33,c066,c0,3,6)
      call chg_basis(aux6,aux33,s066,s0,3,6)
c
      RETURN
      END

C ************************************************************************
C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
C
C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
C     (modif. 10/MAY/01 - KDIM version - R.L.)
C
C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
C **************************************************************************

      SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)

c      PARAMETER (SQR2=1.41421356237309   )
      PARAMETER (RSQ2=0.70710678118654744)
      PARAMETER (RSQ3=0.57735026918962584)
      PARAMETER (RSQ6=0.40824829046386304)

      DIMENSION CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3)

C     DIMENSION B(3,3,6)
C     DATA B /RSQ6,0,   0,   0,   RSQ6,0,   0,   0,  -2*RSQ6,
C    #        RSQ2,0,   0,   0,  -RSQ2,0,   0,   0,   0,
C    #        0,   0,   0,   0,   0,   RSQ2,0,   RSQ2,0,
C    #        0,   0,   RSQ2,0,   0,   0,   RSQ2,0,   0,
C    #        0,   RSQ2,0,   RSQ2,0,   0,   0,   0,   0,
C    #        RSQ3,0,   0,   0,   RSQ3,0,   0,   0,   RSQ3/

      COMMON/BASIS/ B(3,3,6)

      IF(IOPT.EQ.0) THEN
C *** CALCULATES BASIS TENSORS B(N)

        DO I=1,3
          DO J=1,3
            DO N=1,6
              B(I,J,N)=0.0
            ENDDO
          ENDDO
        ENDDO

        B(1,1,2)=-RSQ6
        B(2,2,2)=-RSQ6
        B(3,3,2)= 2.D0*RSQ6

        B(1,1,1)=-RSQ2
        B(2,2,1)= RSQ2

        B(2,3,3)=RSQ2
        B(3,2,3)=RSQ2

        B(1,3,4)=RSQ2
        B(3,1,4)=RSQ2

        B(1,2,5)=RSQ2
        B(2,1,5)=RSQ2

        B(1,1,6)=RSQ3
        B(2,2,6)=RSQ3
        B(3,3,6)=RSQ3

      ENDIF

C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
      IF(IOPT.EQ.1) THEN
        DO 40 I=1,3
        DO 40 J=1,3
        C2(I,J)=0.0
        DO 40 N=1,KDIM
   40   C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
      ENDIF

C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
      IF(IOPT.EQ.2) THEN
        DO 50 N=1,KDIM
        CE2(N)=0.0
        DO 50 I=1,3
        DO 50 J=1,3
   50   CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
      ENDIF

C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
      IF(IOPT.EQ.3) THEN
        DO 20 I=1,3
        DO 20 J=1,3
        DO 20 K=1,3
        DO 20 L=1,3
        C4(I,J,K,L)=0.0
        DO 20 N=1,KDIM
        DO 20 M=1,KDIM
   20   C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
      ENDIF

C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
      IF(IOPT.EQ.4) THEN
        DO 30 N=1,KDIM
        DO 30 M=1,KDIM
        CE4(N,M)=0.0
        DO 30 I=1,3
        DO 30 J=1,3
        DO 30 K=1,3
        DO 30 L=1,3
   30   CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
      ENDIF

      RETURN
      END
c
C **************************************************************************
C     SUBROUTINE UPDATE_SCHMID
C
C     ROTATES SCHMID TENSORS OF EACH GRAIN FROM CRYSTAL TO SAMPLE AXES
C **************************************************************************
      SUBROUTINE UPDATE_SCHMID
c
      INCLUDE 'fft.dim'
c
      DIMENSION aux5(5),aux33(3,3),aux55(5,5),aux3333(3,3,3,3)
      dimension aux33r(3,3)
c
      do ii=1,npts1
      do jj=1,npts2
      do kk=1,npts3
c
      jph=jphase(ii,jj,kk)
cg
      if(igas(jph).eq.0) then
cg
      do 120 is=1,nsyst(jph)
c
      do j=1,5
        aux5(j)=schca(j,is,jph)
      enddo
c
      call chg_basis(aux5,aux33,aux55,aux3333,1,5)
c
      do 140 i=1,3
      do 140 j=1,3
        aux33r(i,j)=0.
      do 140 i1=1,3
      do 140 j1=1,3
        aux33r(i,j)=aux33r(i,j)+
     #      ag(i,i1,ii,jj,kk)*ag(j,j1,ii,jj,kk)*aux33(i1,j1)
cw
cw      aux33r(i,j)=aux33r(i,j)+
cw     #ag(i,i1,jgrain(ii,jj,kk))*ag(j,j1,jgrain(ii,jj,kk))*aux33(i1,j1)
cw
140   continue
c
      call chg_basis(aux5,aux33r,aux55,aux3333,2,5)
c
      do j=1,5
        sch(j,is,ii,jj,kk)=aux5(j)
      enddo
c
120   CONTINUE
cg
      endif    !  igas endif
cg
      ENDDO
      ENDDO
      ENDDO
c
      return
      end

C *************************************************************************
      SUBROUTINE LU_INVERSE(A,N)
c
c   INVERTS A MATRIX USING LU DECOMPOSITION
c
      DIMENSION A(N,N),Y(N,N),INDX(N)
c
      DO I=1,N
        DO J=1,N
          Y(I,J)=0.
        ENDDO
        Y(I,I)=1.
      ENDDO
c
      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
c
      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO
c
      DO I=1,N
      DO J=1,N
       A(I,J)=Y(I,J)
      ENDDO
      ENDDO
c
      RETURN
      END
c
c *************************************************************************
c
      SUBROUTINE LU_EQSYSTEM(A,B,N,ISINGULAR)
c
c     SOLVES A*X=B USING LU DECOMPOSITION
c
      DIMENSION A(N,N),B(N),INDX(N)
c
      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
c
      IF(ISINGULAR.EQ.1) RETURN
c
      CALL LUBKSB(A,N,N,INDX,B)
c
      RETURN
      END
c
C *****************************************************************************
c
      SUBROUTINE ludcmp(a,n,np,indx,d,isingular)
      INTEGER n,np,indx(n),NMAX
c      REAL d,a(np,np),TINY
c      PARAMETER (NMAX=500,TINY=1.0e-20)
      REAL d,a(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,imax,j,k,isingular
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
c
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
c
        if(aamax.eq.0.) then
        isingular=1
        return
        endif
c
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
c
c        if(a(j,j).eq.0.) a(j,j)=TINY
c
        if(a(j,j).eq.0.) then
        isingular=1
        return
        endif
c
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
c
      isingular=0
c
      return
      END
c
c *****************************************************************************
c
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NROWSxNCOLS MATRICES
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
C     OF BOTH DATA.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION tmismatch (v1,v2,nrows,ncols)
C
COLD      DIMENSION v1(nrows,ncols),v2(nrows,ncols)
COLD      DIMENSION v_dif(6,6),v_ave(6,6)
      DIMENSION v1(36),v2(36)
      DIMENSION v_dif(36),v_ave(36)
C
      do i=1,nrows*ncols
        v_dif(i)=v1(i)-v2(i)
        v_ave(i)=0.5d0*(v1(i)+v2(i))
      enddo
      tmismatch=tnorm(v_dif,nrows,ncols)/tnorm(v_ave,nrows,ncols)
C
      RETURN
      END

C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION tnorm(v,nrows,ncols)
C
COLD      DIMENSION v(nrows,ncols)
      DIMENSION v(36)
C
      tnorm=0.d0
      do i=1,nrows*ncols
        tnorm=tnorm+v(i)*v(i)
      enddo
      tnorm=sqrt(tnorm)
C
      RETURN
      END
c
      SUBROUTINE MINV (A,N,D,L,M)
c
      DIMENSION A(*),L(*),M(*)
C
C       SEARCH FOR LARGEST ELEMENT
C
      D=1.d0
      NK=-N
      DO 180 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF (ABS(BIGA)-ABS(A(IJ))) 10,20,20
   10 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C       INTERCHANGE ROWS
C
      J=L(K)
      IF (J-K) 50,50,30
   30 KI=K-N
      DO 40 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   40 A(JI)=HOLD
C
C       INTERCHANGE COLUMNS
C
   50 I=M(K)
      IF (I-K) 80,80,60
   60 JP=N*(I-1)
      DO 70 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   70 A(JI)=HOLD
C
C       DIVIDE COLUMN BY MINUS PIVOT (BIGA)
C
   80 IF (ABS(BIGA).LT.1.e-10) THEN
   90 D=0.0
      RETURN
      ENDIF
  100 DO 120 I=1,N
      IF (I-K) 110,120,110
  110 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
  120 CONTINUE
C
C       REDUCE MATRIX
C
      DO 150 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 150 J=1,N
      IJ=IJ+N
      IF (I-K) 130,150,130
  130 IF (J-K) 140,150,140
  140 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
  150 CONTINUE
C
C       DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 170 J=1,N
      KJ=KJ+N
      IF (J-K) 160,170,160
  160 A(KJ)=A(KJ)/BIGA
  170 CONTINUE
C
C       PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C       REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.d0/BIGA
  180 CONTINUE
C
C       FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  190 K=(K-1)
      IF (K) 260,260,200
  200 I=L(K)
      IF (I-K) 230,230,210
  210 JQ=N*(K-1)
      JR=N*(I-1)
      DO 220 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  220 A(JI)=HOLD
  230 J=M(K)
      IF (J-K) 190,190,240
  240 KI=K-N
      DO 250 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  250 A(JI)=HOLD
      GO TO 190
  260 RETURN
      END
c
cc      subroutine statactiv
cc
cc      include 'fft.dim'
ccc
cc      do 4000 ii=1,npts1
cc      do 4000 jj=1,npts2
cc      do 4000 kk=1,npts3
ccc
cc      jph=jphase(ii,jj,kk)
ccc
cc      if(igas(jph).eq.0) then
ccc
cc      is=0
cc      gavgr(ii,jj,kk)=0.
cc
cc      do 4001 k1=1,nmodes(jph)
cc        gtotmodgr(k1,ii,jj,kk)=0.
cc        do 4002 k2=1,nsm(k1,jph)
cc          is=is+1
cc          gabs=abs(gamdot(is,ii,jj,kk))
cc          gtotmodgr(k1,ii,jj,kk)=gtotmodgr(k1,ii,jj,kk)+gabs
cc4002    continue
cccw        gtot(k1,ii,jj,kk,jph)=gtot(k1,ii,jj,kk,jph)*tdot
cc        gavgr(ii,jj,kk)=gavgr(ii,jj,kk)+gtotmodgr(k1,ii,jj,kk)
ccch        gacumgr(ii,jj,kk)=gacumgr(ii,jj,kk)+
ccch     #                    gtotmodgr(k1,ii,jj,kk)*tdot
cc4001  continue
ccc
cc      endif   !  igas endif
ccc
cc4000  continue
cc
cc      gavpc=0.
ccc
cc      do kph=1,nph
cc      if(igas(kph).eq.0) then
cc       do k1=1,nmodes(kph)
cc         gavmod(k1,kph)=0.
cc       enddo
cc      endif
cc      enddo
ccc
cc      do 4005 ii=1,npts1
cc      do 4005 jj=1,npts2
cc      do 4005 kk=1,npts3
ccc
cc      jph=jphase(ii,jj,kk)
cc      if(igas(jph).eq.0) then
ccc
cc      do 4006 k1=1,nmodes(jph)
cc      gavmodgr(k1,ii,jj,kk)=gtotmodgr(k1,ii,jj,kk)/gavgr(ii,jj,kk)
cc4006  gavmod(k1,jph)=gavmod(k1,jph)+gtotmodgr(k1,ii,jj,kk)*wgt
cc      gavpc=gavpc+gavgr(ii,jj,kk)*wgt
ccc
cc      endif
ccc
cc4005  continue
cc
cc      do kph=1,nph
cc      if(igas(kph).eq.0) then
cc       do k1=1,nmodes(kph)
cc         gavmod(k1,kph)=gavmod(k1,kph)/gavpc
cc       enddo
cc      endif
cc      enddo
cc
cccc      write(35,*) (((gavmod(k1,kph)),k1=1,nmodes(iph)),kph=1,2)
cc
cc      return
cc      end
c
      subroutine get_smacro
c
      include 'fft.dim'
c
      dimension sav6(6),sav5(5)
      dimension aux55(5,5),aux66(6,6),aux3333(3,3,3,3)
      DIMENSION IJV(6,2)
cc      dimension scauav(3,3)

      DATA ((IJV(N,M),M=1,2),N=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
c
c     OVERALL STRESS
c
      do ii=1,3
      do jj=1,3
      scauav(ii,jj)=0.
      scauav1(ii,jj)=0.
      do k=1,npts3
      do j=1,npts2
      do i=1,npts1
      scauav(ii,jj)=scauav(ii,jj)+sg(ii,jj,i,j,k)*wgt
      if(jphase(i,j,k).eq.1) then
       scauav1(ii,jj)=scauav1(ii,jj)+sg(ii,jj,i,j,k)*wgt/wph1
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
cw
cw      enddo
cw      enddo
cc
cc    MIXED BC
cc
      do i=1,6
c
      ii=ijv(i,1)
      jj=ijv(i,2)
c
       ddisgradmacro(ii,jj)=0.
c
      if(idsim(i).eq.0) then
c
       do k=1,6
c
       kk=ijv(k,1)
       ll=ijv(k,2)
c
c       write(*,*) 'I,II,JJ,K,KK,LL',i,ii,jj,k,kk,ll
c       pause
c
       ddisgradmacro(ii,jj)=ddisgradmacro(ii,jj)+
     #    s0(ii,jj,kk,ll)*iscau(k)*(scauchy(kk,ll)-scauav(kk,ll))
c
       enddo
c
      endif
c
      enddo
c
      do ii=1,3
      do jj=1,3
      ddisgradmacroacum(ii,jj)=
     #   ddisgradmacroacum(ii,jj)+ddisgradmacro(ii,jj)
      enddo
      enddo
c
      write(*,*) 'DDISGRADMACRO(1,1),(2,2)=',
     #         ddisgradmacro(1,1),ddisgradmacro(2,2)
cc      write(*,*) 'DDISGRADMACROACUM(1,1),(2,2)=',
cc     #         ddisgradmacroacum(1,1),ddisgradmacroacum(2,2)
c
cc      call chg_basis(sav6,scauchy,aux66,aux3333,2,6)
      call chg_basis(sav6,scauav,aux66,aux3333,2,6)
c
      do ii=1,5
      sav5(ii)=sav6(ii)
      enddo
c
      call chg_basis(sav5,sdeviat,aux55,aux3333,1,5)
c
      svm=0.
      do i=1,3
      do j=1,3
      svm=svm+sdeviat(i,j)*sdeviat(i,j)
      enddo
      enddo
      svm=sqrt(3./2.*svm)
cg
      call chg_basis(sav6,scauav1,aux66,aux3333,2,6)
      do ii=1,5
      sav5(ii)=sav6(ii)
      enddo
      call chg_basis(sav5,sdeviat,aux55,aux3333,1,5)
      svm1=0.
      do i=1,3
      do j=1,3
      svm1=svm1+sdeviat(i,j)*sdeviat(i,j)
      enddo
      enddo
      svm1=sqrt(3./2.*svm1)
cg
cc      write(*,*) 'SVM=',svm
cc      pause
      return
      end
c
      subroutine update_orient
c
      include 'fft.dim'
c
      dimension aa(3,3),distor(3,3)
      dimension dnsa(3),dbsa(3)
      dimension rotslip(3,3),rotloc(3,3),rot(3,3)
c
      RSLBAR=0.
      RLCBAR=0.
c
C     MASTER DO
C
      do 2000 k=1,npts3
      do 2000 j=1,npts2
      do 2000 i=1,npts1
c
      iph=jphase(i,j,k)
c
      if(igas(iph).eq.0) then
c
c     LOCAL ROTATION RATE: ANTISYM(VELGRAD)
c
      do 100 ii=1,3
      do 100 jj=1,3
      rotloc(ii,jj)=(velgrad(ii,jj,i,j,k)-velgrad(jj,ii,i,j,k))/2.
cg      rotloc(ii,jj)=0.
  100 continue
c
c     SLIP ROTATION RATE
c
      do ii=1,3
      do jj=1,3
        aa(ii,jj)=ag(ii,jj,i,j,k)
        distor(ii,jj)=0.
      enddo
      enddo

cww        do is=1,nsyst(1)
        do is=1,nsyst(iph)
c
        do ii=1,3
          dnsa(ii)=0.
          dbsa(ii)=0.
          do jj=1,3
cw            dnsa(ii)=dnsa(ii)+aa(ii,jj)*dnca(jj,is,1)
cw            dbsa(ii)=dbsa(ii)+aa(ii,jj)*dbca(jj,is,1)
            dnsa(ii)=dnsa(ii)+aa(ii,jj)*dnca(jj,is,iph)
            dbsa(ii)=dbsa(ii)+aa(ii,jj)*dbca(jj,is,iph)
          enddo
        enddo
c
        do ii=1,3
        do jj=1,3
         distor(ii,jj)=distor(ii,jj)+dbsa(ii)*dnsa(jj)*gamdot(is,i,j,k)
        enddo
        enddo
c
      enddo
C
      do ii=1,3
      do jj=1,3
        rotslip(ii,jj)=(distor(ii,jj)-distor(jj,ii))/2.
      enddo
      enddo
C
C     AVERAGE ROTATION RATE
C
      rslbar=rslbar+sqrt(rotslip(3,2)**2+rotslip(1,3)**2+
     #       rotslip(2,1)**2)*wgt
      rlcbar=rlcbar+sqrt(rotloc(3,2)**2+rotloc(1,3)**2+
     #       rotloc(2,1)**2)*wgt
c
c    TOTAL ROTATION
c
      do ii=1,3
      do jj=1,3
        rot(ii,jj)=(tomtot(ii,jj)+rotloc(ii,jj)-rotslip(ii,jj))*tdot
      enddo
      enddo
c
c     REORIENTATION
c
      call orient(aa,rot)
c
c     UPDATE ORIENTATION MATRIX
c
      do ii=1,3
      do jj=1,3
        ag(ii,jj,i,j,k)=aa(ii,jj)
      enddo
      enddo
c
      endif   !  igas endif
c
2000  continue
c
      write(*,*)
      write(*,*) 'AVERAGE PLASTIC ROTATION =',rslbar
      write(*,*) 'AVERAGE LOCAL ROTATION =',rlcbar
      write(*,*)
c
      RETURN
      END
ccc
      subroutine orient(a,c)
      dimension a(3,3),c(3,3),th2(3,3),v(3),vbar(3)
      dimension th(3,3)
      dimension rot(3,3),anew(3,3)

c     BUILD ROTATION TENSOR BASED ON RODRIGUES FORMULA

      v(1)=c(3,2)
      v(2)=c(1,3)
      v(3)=c(2,1)
      snorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      snorm1=tan(snorm/2.)
      if(snorm.gt.1.e-06) go to 97
      snorm=1.
97    do 20 i=1,3
      vbar(i)=snorm1*v(i)/snorm
20    continue
      snorm=vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3)
      th(3,2)=vbar(1)
      th(1,3)=vbar(2)
      th(2,1)=vbar(3)
      th(2,3)=-vbar(1)
      th(3,1)=-vbar(2)
      th(1,2)=-vbar(3)
      do 40 i=1,3
40    th(i,i)=0.
      do 30 i=1,3
      do 30 j=1,3
      th2(i,j)=0.
      do 50 k=1,3
50    th2(i,j)=th2(i,j)+th(i,k)*th(k,j)
30    continue
      do 60 i=1,3
      do 60 j=1,3
60    rot(i,j)=(i/j)*(j/i)+2.*(th(i,j)+th2(i,j))/(1.+snorm)
      do 70 i=1,3
      do 70 j=1,3
      anew(i,j)=0.
      do 80 k=1,3
80    anew(i,j)=anew(i,j)+rot(i,k)*a(k,j)
70    continue
      do 90 i=1,3
      do 90 j=1,3
90    a(i,j)=anew(i,j)
      return
      end

C *****************************************************************************
C     SUBROUTINE UPDATE_CRSS_VOCE     --->      VERSION OF 23/APR/00
C
C     A VOCE LAW FUNCTION OF THE ACCUMULATED SHEAR IN EACH GRAIN IS ADDED
C     AS A MULTIPLICATIVE FACTOR THAT MODULATES THE ORIGINAL LINEAR HARDENING.
C     THE UNITS (AND THE STRENGTH) OF THE HARDENING ARE CARRIED BY THE
C     MULTIPLICATIVE FACTOR 'VOCE' (THE SAME FOR EVERY MODE).
C     THE SELF & LATENT COUPLING COEFFICIENTS 'HARD' ARE DIMENSIONLESS
C     CONSTANTS RELATIVE TO THE FACTOR 'VOCE'.
*******************************************************************************

      SUBROUTINE HARDEN
c
      INCLUDE 'fft.dim'
c
      do i=1,npts1
      do j=1,npts2
      do k=1,npts3
c
      iph=jphase(i,j,k)
c
       if(igas(iph).eq.0) then
c
          GAMTOTX=GACUMGR(I,J,K)
          DELTGAM=0.0
          DO IS=1,NSYST(IPH)
            DELTGAM=DELTGAM+ABS(GAMDOT(IS,I,J,K))*TDOT
          ENDDO

          DO IS=1,NSYST(IPH)
            DTAU=0.
            DO JS=1,NSYST(IPH)
              DTAU=DTAU+HARD(IS,JS,IPH)*ABS(GAMDOT(JS,I,J,K))*TDOT
            ENDDO
            TAU0 =TAU (IS,1,IPH)
            TAU1 =TAU (IS,3,IPH)
            THET0=THET(IS,0,IPH)
            THET1=THET(IS,1,IPH)
            TINY=1.E-4*TAU0

            VOCE=0.0
            IF(ABS(THET0).GT.TINY) THEN
              VOCE=THET1*DELTGAM
              IF(ABS(TAU1).GT.TINY) THEN
                FACT=ABS(THET0/TAU1)
                EXPINI=EXP(-GAMTOTX*FACT)
                EXPDEL=EXP(-DELTGAM*FACT)
                VOCE  =VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
     #            (EXPDEL-1.)-THET1/FACT*EXPINI*
     #            (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
              ENDIF
            ENDIF
            CRSS(IS,1,I,J,K)=CRSS(IS,1,I,J,K)+DTAU*VOCE/DELTGAM
            CRSS(IS,2,I,J,K)=CRSS(IS,2,I,J,K)+DTAU*VOCE/DELTGAM
chard
             TRIALTAU(IS,1,I,J,K)=CRSS(IS,1,I,J,K)
             TRIALTAU(IS,2,I,J,K)=CRSS(IS,2,I,J,K)
chard
          ENDDO
          GACUMGR(I,J,K)=GAMTOTX+DELTGAM
ch
ch        DO IS=1,NSYST(IPH)
chc
ch          dtau=0.
ch            do js=1,nsyst(iph)
ch              dtau=dtau+hard(is,js,iph)*abs(gamdot(js,i,j,k))*tdot
ch            enddo
ch
ch          thet0=thet(is,0,iph)
ch          thet1=thet(is,1,iph)
ch          fact =gacumgr(i,j,k)*thet0/tau(is,3,iph)
ch          voce=thet1+(thet0-thet1+thet1*fact)*exp(-fact)
chc
ch          crss(is,1,i,j,k)=crss(is,1,i,j,k)+dtau*voce
ch          crss(is,2,i,j,k)=crss(is,2,i,j,k)+dtau*voce
chc
ch        ENDDO
ch
      endif   !  igas endif
c
      enddo
      enddo
      enddo
c
      END
cevp
      SUBROUTINE DATA_CRYSTAL_ELAST(IPH)

      INCLUDE 'fft.dim'
c
	dimension dde(3,3),xid4(3,3,3,3)
        dimension cc66v(6,6),ccaux(3,3,3,3)
c
c       UNITARY TENSORS
c
	do i=1,3
	do j=1,3
	dde(i,j)=0.d0
	if(i.eq.j) dde(i,j)=1.d0
	enddo
	enddo
c
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	xid4(i,j,k,l)=(dde(i,k)*dde(j,l)+dde(i,l)*dde(j,k))/2.d0
	enddo
	enddo
	enddo
	enddo
c
	read(UR1,*) iso
c
	if(iso.eq.0) then
c
	do i=1,6
	read(UR1,*)(cc66v(i,j),j=1,6)
	enddo
cg	call voigt(cc66v,cc,1)
	call voigt(cc66v,ccaux,1)
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	cc(i,j,k,l,iph)=ccaux(i,j,k,l)
	enddo
	enddo
	enddo
	enddo
cg
c
	else
c
cw	read(UR1,*) tmu,tnu
	read(UR1,*) young,tnu
        tmu=young/(2.*(1.+tnu))
	tla=2.d0*tmu*tnu/(1.d0-2.d0*tnu)
c       bulk=2.d0*tmu*(1+tnu)/3.d0/(1.d0-2d0*tnu)
c
c	rho=3.d0*(1.d0-2.d0*tnu)/(1.d0+tnu)
c
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
cg	cc(i,j,k,l)=tla*dde(i,j)*dde(k,l)+2.d0*tmu*xid4(i,j,k,l)
	cc(i,j,k,l,iph)=tla*dde(i,j)*dde(k,l)+2.d0*tmu*xid4(i,j,k,l)
	enddo
	enddo
	enddo
	enddo
c
     	endif
c
      RETURN
      END
c
      subroutine evpal(imicro)

      INCLUDE 'fft.dim'

      dimension xlambda(3,3),xlambda6(6)
c
      dimension sgaux(3,3),sg6(6),sg6old(6)

cx      dimension dsg(3,3),dse(6)

      dimension eptaux(3,3),ept6(6)
      dimension edotpaux(3,3),edotp6(6)

      dimension dedotp66(6,6)

cg      dimension disgradaux(3,3),disgrad6(6)
      dimension strainaux(3,3),strain6(6)

      dimension aux66(6,6),aux3333(3,3,3,3)

      dimension sg66(6,6)

cg      dimension disgradceq(3,3),disgradceq6(6)
      dimension strainceq(3,3),strainceq6(6)
      dimension res(6)
      dimension xjacobinv(6,6)
c
      ERRE=0.
      ERRS=0.
cth
      write(*,*)
cth
      do 877 i=1,npts1
cth
caps      write(*,'(1H+,i5,a)') int(100.*i/npts1),'% COMPLETED'
cth
      do 877 j=1,npts2
      do 877 k=1,npts3
cg
       jph=jphase(i,j,k)
cg
      if(igas(jph).eq.0) then
cg
       do ii=1,6
       do jj=1,6
        sg66(ii,jj)=cg66(ii,jj,i,j,k)
       enddo
       enddo
c
      call lu_inverse(sg66,6)
c
cx      write(*,*) 'SG66='
cx      do ii=1,6
cx      write (*,*) (sg66(ii,jj),jj=1,6)
cx      enddo
cx      pause
c
       do ii=1,3
       do jj=1,3
        xlambda(ii,jj)=sg(ii,jj,i,j,k)
        sgaux(ii,jj)=sg(ii,jj,i,j,k)
        eptaux(ii,jj)=ept(ii,jj,i,j,k)
cg        disgradaux(ii,jj)=disgrad(ii,jj,i,j,k)
        strainaux(ii,jj)=(disgrad(ii,jj,i,j,k)+disgrad(jj,ii,i,j,k))/2.
       enddo
       enddo
c
      call chg_basis(xlambda6,xlambda,aux66,aux3333,2,6)

      call chg_basis(sg6,sgaux,aux66,aux3333,2,6)
      call chg_basis(ept6,eptaux,aux66,aux3333,2,6)
cg      call chg_basis(disgrad6,disgradaux,aux66,aux3333,2,6)
      call chg_basis(strain6,strainaux,aux66,aux3333,2,6)
c
      sgnorm=0.
      dgnorm=0.
c
      do ii=1,3
      do jj=1,3
       sgnorm=sgnorm+xlambda(ii,jj)**2
cg       dgnorm=dgnorm+disgradaux(ii,jj)**2
       dgnorm=dgnorm+strainaux(ii,jj)**2
      enddo
      enddo
c
      sgnorm=sqrt(sgnorm)
      dgnorm=sqrt(dgnorm)
c
cx      write(*,*) 'SGNORM,DGNORM=',sgnorm,dgnorm
cx      pause
c
      erroral=0.0000001
      itmaxal=100
      iter1=0
cxx      erral=2*error
      erral=2*erroral
c
cxx      do while(iter1.lt.itmax.and.erral.gt.error)
      do while(iter1.lt.itmaxal.and.erral.gt.erroral)
      iter1=iter1+1

      sg6old=sg6

cth      call edotp_sg_eval(sg6,edotp6,dedotp66,i,j,k,jph)
cthh      if(ithermo.ne.1)
      if(ithermo.ne.1.or.imicro.ne.1) then
        call edotp_sg_eval(sg6,edotp6,dedotp66,i,j,k,jph)
      else
        edotp6=0.
      endif
c
      do ii=1,6
cth
cthh      if(ithermo.eq.1) then
      if(ithermo.eq.1.and.imicro.eq.1) then
cx
cx    THERMOELASTIC
cx
cth       strainceq6(ii)=0.
       strainceq6(ii)=ept6(ii)
      else
       strainceq6(ii)=ept6(ii)+edotp6(ii)*tdot
      endif
cth
       do jj=1,6
        strainceq6(ii)=strainceq6(ii)+sg66(ii,jj)*sg6(jj)
       enddo
      enddo

cgcx      write(*,*) 'disgrad6=',disgrad6
cgcx      write(*,*) 'disgradceq6=',disgradceq6
cgcx      pause
cx      write(*,*) 'strain6=',strain6
cx      write(*,*) 'strainceq6=',strainceq6
cx      pause

cg      call chg_basis(disgradceq6,disgradceq,aux66,aux3333,1,6)
      call chg_basis(strainceq6,strainceq,aux66,aux3333,1,6)

cx      if(i.eq.2.and.j.eq.2.and.k.eq.2) then
cx      write(*,*) 'disgradceq='
cx      do ii=1,3
cx      write(*,*) (disgradceq(ii,jj),jj=1,3)
cx      enddo
cx      pause
cx      endif

      do ii=1,6
       res(ii)=sg6(ii)-xlambda6(ii)
      do jj=1,6
cg       res(ii)=res(ii)+c066(ii,jj)*(disgradceq6(jj)-disgrad6(jj))
       res(ii)=res(ii)+c066(ii,jj)*(strainceq6(jj)-strain6(jj))
      enddo
      enddo

cx      write(*,*) 'res=',res
cx      pause

      do ii=1,6
      do jj=1,6
cbug       xjacobinv(ii,jj)=(ii/jj)*(jj*ii)
       xjacobinv(ii,jj)=(ii/jj)*(jj/ii)
       do kk=1,6
cthh        if(ithermo.eq.1) then
        if(ithermo.eq.1.and.imicro.eq.1) then
cx
cx    ELASTIC
cx
         xjacobinv(ii,jj)=xjacobinv(ii,jj)+
     #                  c066(ii,kk)*sg66(kk,jj)
        else
         xjacobinv(ii,jj)=xjacobinv(ii,jj)+
     #                  c066(ii,kk)*(sg66(kk,jj)+dedotp66(kk,jj)*tdot)
        endif
       enddo
      enddo
      enddo
c
      call lu_inverse(xjacobinv,6)
c
      do ii=1,6
      do jj=1,6
      sg6(ii)=sg6(ii)-xjacobinv(ii,jj)*res(jj)
      enddo
      enddo
c
      dsgnorm1=0.
      dsgnorm2=0.
      ddgnorm=0.

      do ii=1,6
       dsgnorm1=dsgnorm1+(sg6(ii)-sg6old(ii))**2
       dsgnorm2=dsgnorm2+(sg6(ii)-xlambda6(ii))**2
cg       ddgnorm=ddgnorm+(disgradceq6(ii)-disgrad6(ii))**2
       ddgnorm=ddgnorm+(strainceq6(ii)-strain6(ii))**2
      enddo
c
      erral=dsgnorm1/sgnorm
      errald=ddgnorm/dgnorm

cx      if(i.eq.2.and.j.eq.2.and.k.eq.2) then
cx      write(*,*) 'FP',i,j,k,'ITER1,ERRAL,ERRALD=',iter1,erral,errald
cx      pause
cx      endif
chard
cth      call get_trialtau(i,j,k,jph)
      if(iuphard.eq.1.and.(ithermo.ne.1.or.imicro.gt.2))
cth      if(iuphard.eq.1.and.(ithermo.ne.1.or.imicro.gt.1))
     #   call get_trialtau(i,j,k,jph)
chard
      enddo  ! end do while (iter1)
c
      call chg_basis(sg6,sgaux,aux66,aux3333,1,6)
      call chg_basis(edotp6,edotpaux,aux66,aux3333,1,6)
c
       do ii=1,3
       do jj=1,3
        sg(ii,jj,i,j,k)=sgaux(ii,jj)
        edotp(ii,jj,i,j,k)=edotpaux(ii,jj)
       enddo
       enddo

       ERRS=ERRS+dsgnorm2*WGT
       ERRE=ERRE+ddgnorm*WGT

cx      if(i.eq.2.and.j.eq.2.and.k.eq.2) then
cx      write(*,*) 'SGij(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (sg(ii,jj,i,j,k),jj=1,3)
cx      enddo
cx      write(*,*) 'DISGRADij(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (disgradaux(ii,jj),jj=1,3)
cx      enddo
cx      write(*,*) 'DGCEQij(2,2,2)'
cx      do ii=1,3
cx      write(*,*) (disgradceq(ii,jj),jj=1,3)
cx      enddo
cx      pause
cx      endif
cg
      else    !   igas else
cg
       do ii=1,3
       do jj=1,3
        sg(ii,jj,i,j,k)=0.
        edotp(ii,jj,i,j,k)=0.    !  ????
       enddo
       enddo
cg
      endif    !   igas endif
cg
877   continue

      write(*,*) 'ERRE=',erre
      write(*,*) 'ERRS=',errs
cx      pause

      return
      end
c
      subroutine edotp_sg_eval(sg6,edotp6,dedotp66,i,j,k,jph)

      INCLUDE 'fft.dim'
c
      dimension sg6(6)
      dimension edotp6(6),dedotp66(6,6)
c
      dimension rss(NSYSMX)
      dimension rss1(NSYSMX),rss2(NSYSMX)
      dimension sc(5,NSYSMX),taux(NSYSMX,2),nrsx(NSYSMX)
cth
      dimension xkinaux(NSYSMX)
cth
c
      DO is=1,nsyst(jph)

        nrsx(is)=nrs(is,jph)
c
chard        taux(is,1)=crss(is,1,i,j,k)
chard        taux(is,2)=crss(is,2,i,j,k)
        taux(is,1)=trialtau(is,1,i,j,k)
        taux(is,2)=trialtau(is,2,i,j,k)
chard
cth
        xkinaux(is)=xkin(is,i,j,k)
cth
        DO jj=1,5
          sc(jj,is)=sch(jj,is,i,j,k)
        ENDDO
      ENDDO

C     GET RESOLVED SHEAR STRESSES 'rss' AND SHEAR RATES 'gamdot'.
C     SIGN(GAMDOT)=SIGN(RSS).
C     NRS CAN BE EVEN OR ODD.
C     RSS1 IS ALWAYS > 0 AND IS USED TO CALCULATE VISCOUS COMPLIANCE.

      do is=1,nsyst(jph)
        rss(is)=sc(1,is)*sg6(1)+sc(2,is)*sg6(2)+sc(3,is)*sg6(3)+
     #          sc(4,is)*sg6(4)+sc(5,is)*sg6(5)
        isign=1
cth        if(rss(is).lt.0.) isign=2
cth        rss(is)=rss(is)/taux(is,isign)
        if((rss(is)-xkinaux(is)).lt.0.) isign=2
        rss(is)=(rss(is)-xkinaux(is))/taux(is,isign)
        rss1(is)=gamd0(is,jph)*nrsx(is)*
     #           abs(rss(is)**(nrsx(is)-1))/taux(is,isign)
        rss2(is)=
     #  gamd0(is,jph)*abs(rss(is)**(nrsx(is)))*sign(1.,rss(is))
chard
      gamdot(is,i,j,k)=rss2(is)
chard
      enddo
c
      DO II=1,5
      edotp6(II)=0.
        DO K1=1,NSYST(jph)
          edotp6(II)=edotp6(II)+SC(II,K1)*RSS2(K1)
        ENDDO
      ENDDO
c
      edotp6(6)=0.
c
      DO II=1,5
      DO JJ=1,5
       dedotp66(II,JJ)=0.
       DO K1=1,NSYST(jph)
         dedotp66(II,JJ)=
     #   dedotp66(II,JJ)+SC(II,K1)*SC(JJ,K1)*RSS1(K1)
       ENDDO
      ENDDO
      ENDDO

      do ii=1,6
       dedotp66(II,6)=0.
       dedotp66(6,II)=0.
      enddo

      RETURN
      END
c
      function vm(dtensor)
c
c     CALCULATES THE VM EQUIVALENT OF A NON-SYMMETRIC, NON-TRACELESS
c     (displacement GRADIENT OR VELOCITY GRADIENT) TENSOR
c
      dimension dtensor(3,3),dt(3,3)
c
      trace=dtensor(1,1)+dtensor(2,2)+dtensor(3,3)
c
      do i=1,3
      do j=1,3
       dt(i,j)=(dtensor(i,j)+dtensor(j,i))/2.-(i/j)*(j/i)*trace/3.
      enddo
      enddo
c
      vm=0.
      do i=1,3
      do j=1,3
      vm=vm+dt(i,j)**2
      enddo
      enddo
      vm=sqrt(2./3.*vm)
c
      return
      end

      subroutine get_trialtau(i,j,k,jph)

      INCLUDE 'fft.dim'

      GAMTOTX=GACUMGR(I,J,K)
      DELTGAM=0.0
      DO IS=1,NSYST(JPH)
       DELTGAM=DELTGAM+ABS(GAMDOT(IS,I,J,K))*TDOT
      ENDDO

      DO IS=1,NSYST(JPH)

       DTAU=0.
       DO JS=1,NSYST(JPH)
         DTAU=DTAU+HARD(IS,JS,JPH)*ABS(GAMDOT(JS,I,J,K))*TDOT
       ENDDO
       TAU0 =TAU (IS,1,JPH)
       TAU1 =TAU (IS,3,JPH)
       THET0=THET(IS,0,JPH)
       THET1=THET(IS,1,JPH)
       TINY=1.E-4*TAU0

       VOCE=0.0
       IF(ABS(THET0).GT.TINY) THEN
        VOCE=THET1*DELTGAM
        IF(ABS(TAU1).GT.TINY) THEN
         FACT=ABS(THET0/TAU1)
         EXPINI=EXP(-GAMTOTX*FACT)
         EXPDEL=EXP(-DELTGAM*FACT)
         VOCE=VOCE-(FACT*TAU1-THET1)/FACT*EXPINI*
     #   (EXPDEL-1.)-THET1/FACT*EXPINI*
     #   (EXPDEL*((GAMTOTX+DELTGAM)*FACT+1.)-(GAMTOTX*FACT+1.))
        ENDIF
      ENDIF

      TRIALTAU(IS,1,I,J,K)=CRSS(IS,1,I,J,K)+DTAU*VOCE/DELTGAM
      TRIALTAU(IS,2,I,J,K)=CRSS(IS,2,I,J,K)+DTAU*VOCE/DELTGAM

      ENDDO

      RETURN
      END
c
      subroutine kinhard_param

      INCLUDE 'fft.dim'

      dimension sc5(5),sg5(5),sgx(3,3)
      dimension aux55(5,5),aux3333(3,3,3,3)

      do 7277 i=1,npts1
      do 7277 j=1,npts2
      do 7277 k=1,npts3

       jph=jphase(i,j,k)

       if(igas(jph).eq.0) then

        do ii=1,3
        do jj=1,3
         sgx(ii,jj)=sg(ii,jj,i,j,k)
        enddo
        enddo
c
        call chg_basis(sg5,sgx,aux55,aux3333,2,5)
c
        do is=1,nsyst(jph)

         DO jj=1,5
           sc5(jj)=sch(jj,is,i,j,k)
         ENDDO
c
         xkin(is,i,j,k)=sc5(1)*sg5(1)+sc5(2)*sg5(2)+
     #      sc5(3)*sg5(3)+sc5(4)*sg5(4)+sc5(5)*sg5(5)
        enddo

       endif

7277  continue
c
      return
      end
