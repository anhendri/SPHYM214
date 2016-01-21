!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMME PRINCIPAL
program ELYTRE

  implicit none
  integer,parameter              :: R =  4

  complex(kind = R), parameter   :: i = cmplx(0,1)
  complex(kind = R), allocatable :: EXPO(:)
  complex(kind = R)              :: Wu(2,2), Wd(2,2), Vu(2), Vd(2), CSTE, A
  real(kind = R), parameter      :: c = 299792458, mu0 = 12.5663706E-7, pi = 3.1415926535897932
  real(kind = R), allocatable    :: LGR0(:), LNGR(:), EPS0(:), kz(:), EPS(:)
  real(kind = R)                 :: EPSI, EPSS, LMBD, LMIN, LMAX, PASL, Jx, JUPS ,w, JUP0, SGNE(2,2), SGNK(2,2)
  integer                        :: NPER, NCCH, INDX, NMBR, IDXL, POSS

    open (10,file="INPUT.dat")
    open (20,file="OUTPUT.dat")

    read (10,*) LMIN, LMAX, PASL, EPSI, EPSS, NCCH, NPER

    allocate(LGR0(NPER),EPS0(NPER),EPS(0:NCCH+1),kz(0:NCCH+1),LNGR(NCCH))

    read (10,*) (LGR0(INDX), EPS0(INDX), INDX=1,NPER)

    EPS(0) = EPSI
    do INDX=1,NCCH
        EPS(INDX)  = EPS0(mod(INDX-1,NPER)+1)
        LNGR(INDX) = LGR0(mod(INDX-1,NPER)+1)*1.E-9
    enddo
    EPS(NCCH+1) = EPSS

    read (10,*) POSS, Jx

    SGNE = reshape([+1,-1,+1,-1],[2,2])
    SGNK = reshape([+1,-1,-1,+1],[2,2])

    do IDXL=0,int((LMAX-LMIN)/PASL)
        LMBD = LMIN+IDXL*PASL

        w    = 2.E9*pi*c/LMBD
        CSTE = -i*mu0*Jx*c**2/(EPS(POSS)*w)
        kz   = sqrt((w/c)**2*EPS)
        EXPO = exp(i*kz(1:NCCH)*LNGR)

        Wu = (.5+SGNK*kz(NCCH+1)/(2.*kz(POSS)))/EXPO(POSS)**SGNE
        Wd = (.5+SGNK*kz(0)/(2.*kz(POSS)))
        Vu = -CSTE/2.*EXPO(POSS)**[-1,+1]
        Vd = -CSTE/2.

        A    = (Wd(2,2)*(Vd(1)-Vu(1))+Wd(1,2)*(Vu(2)-Vd(2)))/(Wd(2,2)*Wu(1,1)-Wd(1,2)*Wu(2,1))
        JUP0 = kz(NCCH+1)*cabs(A)**2/(2.*mu0*w)

        Wu = (.5+SGNK*kz(POSS+1)/(2.*kz(POSS)))/EXPO(POSS)**SGNE
        Wd = (.5+SGNK*kz(POSS-1)/(2.*kz(POSS)))

        do INDX=POSS+1,NCCH
            Wu = Matmul(Wu,(.5+SGNK*kz(INDX+1)/(2.*kz(INDX)))/EXPO(INDX)**SGNE)
        enddo

        do INDX=POSS-1,1,-1
            Wd = Matmul(Wd,(.5+SGNK*kz(INDX-1)/(2.*kz(INDX)))*(EXPO(INDX)**SGNE))
        enddo

        A    = (Wd(2,2)*(Vd(1)-Vu(1))+Wd(1,2)*(Vu(2)-Vd(2)))/(Wd(2,2)*Wu(1,1)-Wd(1,2)*Wu(2,1))
        JUPS = kz(NCCH+1)*cabs(A)**2/(2.*mu0*w)

        write (20,'(2E13.5)') LMBD, JUPS/JUP0
    enddo

end program ELYTRE