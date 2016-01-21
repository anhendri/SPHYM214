!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMME PRINCIPAL
program TRANSFERT

  implicit none
  integer,parameter              :: R =  4

  complex(kind = R), parameter   :: i = cmplx(0,1)
  complex(kind = R), allocatable :: EXPO(:), EPS1(:), kz(:), EPS(:)
  complex(kind = R)              :: Wd(2,2), w, ky, f
  real(kind = R), parameter      :: c = 299792458, pi = 3.1415926535897932
  real(kind = R), allocatable    :: LNGR(:), LGR0(:), EPS0(:)
  real(kind = R)                 :: EPSI, EPSK, EPSS, LMBD, LMIN, LMAX, PASL, ANGL, SGNE(2,2), SGNK(2,2), RFLC,TRNS,BAND
  integer                        :: NCCH, INDX, NMBR, IDXL, TETM, NPER

    open (10,file="INPUT.dat")
    open (20,file="OUTPUT.dat")

    read (10,*) ANGL, LMIN, LMAX, PASL, TETM, NCCH, NPER, EPSI, EPSS

    allocate(eps(0:NCCH+1),eps1(0:NCCH+1), LNGR(NCCH), kz(0:NCCH+1), LGR0(NPER),EPS0(NPER))

    read (10,*) (LGR0(INDX), EPS0(INDX), INDX=1,NPER)

    EPS(0) = cmplx(EPSI,0)
    do INDX=1,NCCH
        EPS(INDX)  = EPS0(mod(INDX-1,NPER)+1)
        LNGR(INDX) = LGR0(mod(INDX-1,NPER)+1)*1.E-9
    enddo
    EPS(NCCH+1) = cmplx(EPSS,0)
    EPS1 = EPS**TETM

    SGNE = reshape([+1,-1,+1,-1],[2,2])
    SGNK = reshape([+1,-1,-1,+1],[2,2])

    do IDXL = 0,int((LMAX-LMIN)/PASL)
        LMBD = LMIN+IDXL*PASL
        w    = cmplx(2.E9*pi*c/LMBD,0.)
        ky   = (w/c)*sin(ANGL*360./(2.*pi))*sqrt(eps(0))
        kz   = sqrt((w/c)**2.*eps-ky**2)
        where((real(kz)<0.and.aimag(kz)==0).OR.aimag(kz)<0) kz = - kz

        EXPO = exp(i*kz(1:NCCH)*LNGR)

        Wd = (.5+SGNK*eps1(0)*kz(1)/(2.*eps1(1)*kz(0)))
        do INDX = 1,NCCH
            Wd = Matmul(Wd,(.5+SGNK*eps1(INDX)*kz(INDX+1)/(2.*eps1(INDX+1)*kz(INDX)))/(EXPO(INDX)**SGNE))
        enddo

        RFLC = cabs(Wd(2,1)/Wd(1,1))**2
        TRNS = (kz(NCCH+1)/kz(0))*(eps1(0)/eps1(NCCH+1))*cabs(1/Wd(1,1))**2

        f = real(((eps1(2)*kz(1))**2+(eps1(1)*kz(2))**2)/(2*eps1(1)*kz(1)*eps1(2)*kz(2)))
        BAND = 1./(LGR0(1)+LGR0(2))*acos(cos(kz(1)*lngr(1))*cos(kz(2)*lngr(2))-f*sin(kz(1)*lngr(1))*sin(kz(2)*lngr(2)))

        write (20,*) LMBD, RFLC*100, TRNS*100, 1-RFLC-TRNS, BAND
    enddo

end program TRANSFERT