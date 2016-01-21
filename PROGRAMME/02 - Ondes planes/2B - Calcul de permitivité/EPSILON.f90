!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMME PRINCIPAL
program STRUCTURE
  implicit none

  integer,parameter              :: Rk =  8

  real(kind=rk) :: a=1,f=0.74,Vc,eps1=1.,eps2=2.25,c=299792458,pi=3.1415926535, normeG, eps, nnum1,nnum2,nnum,epsm1,epsm2,&
                   epsm,nth, R
  real(kind=rk) :: e1(3),e2(3),e3(3),a1(3),a2(3),a3(3),b1(3),b2(3),b3(3), G1(3), vr(3), epsrg(50,50),epsr
  real(kind=rk),allocatable :: G2(:,:), G(:,:), k(:,:),h(:,:), Vp(:,:), Vp1(:), kn(:), vcp(:,:), w(:,:),fv1(:),fv2(:)
  integer :: Nmax = 0, N = 0, Gmax = 15, i=0, nk = 150, i1, i2, i3, ng, ng1,m, size1, size2, info, j
  integer, allocatable     :: WORK(:)
  character(len=1) :: VECT = "N", STOC = "U";


    Vc = a**3/4;

    R=(f*a**3.*3./(16.*pi))**(1./3.);

    e1=[1,0,0];
    e2=[0,1,0];
    e3=[0,0,1];

    a1=a/2*[0,1,1];
    a2=a/2*[1,0,1];
    a3=a/2*[1,1,0];

    b1=2*pi*cross(a2,a3)/(dot_product(a1,cross(a2,a3)));
    b2=2*pi*cross(a3,a1)/(dot_product(a1,cross(a2,a3)));
    b3=2*pi*cross(a1,a2)/(dot_product(a1,cross(a2,a3)));
    
    N = floor(Gmax/norm2(b1))
    i = count((/((/((/(norm2(i1*b1+i2*b2+i3*b3)<Gmax,i1=-N,N)/),i2=-N,N)/),i3=-N,N)/))

    allocate(G2(i,4),G(i,3))
    i = 0
    do i1=-N,N
        do i2=-N,N
            do i3=-N,N
                G1=i1*b1+i2*b2+i3*b3;
                if (norm2(G1)<Gmax) then
                    i=i+1;
                    G2(i,:)=[G1,norm2(G1)];
                endif
            enddo
        enddo
    enddo

    G2 = G2(QuickSort((/(j,j=1,i)/)),:);

    nG = size(G2,1);
    G  = G2(:,1:3);

    do i1=1,10
        vr(1) = a*i1/20.
        do i2=1,50
            vr(2) = a*i2/50.
            do i3=1,50
                vr(3) = a*i3/50.
                epsr = 0;
                do j=1,nG
                    normeG=norm2(G(j,:));
                    if (normeG==0) then
                        eps=1/eps1+(1/eps2-1/eps1)*f;
                    else
                        eps=(1/eps2-1/eps1)/Vc*(4*pi*(sin(normeG*R)-normeG*R*cos(normeG*R)))/(normeG**3);
                    endif
                    epsr = epsr+eps*exp(cmplx(0,dot_product(G(j,:),vr)));
                enddo
                epsrg(i2,i3)=real(1/epsr);
            enddo
        enddo
    enddo
    open(10,file='OUTPUT.dat')
    write(10,'(1E13.5)') real(size(epsrg)),epsrg

    contains

    function cross(VCT1,VCT2) result(CRSS)
      real(kind = rk) :: CRSS(3)
      real(kind = rk), intent(in)  :: VCT1(:),VCT2(:)

        CRSS = VCT1([2,3,1])*VCT2([3,1,2])-VCT1([3,1,2])*VCT2([2,3,1])
    end function cross

    recursive function QuickSort(InList) result(OutList)
      integer, dimension(:) :: InList
      integer, allocatable  :: Suplist(:), InfList(:), OutList(:)
      real                  :: pivot

        if(size(InList,1)<2) then
            OutList = Inlist
        else
            pivot = G2(InList(1),4)
            InfList = QuickSort(pack(Inlist(2:), G2(InList(2:),4) < Pivot))
            SupList = QuickSort(pack(Inlist(2:), G2(InList(2:),4) >= Pivot))
            OutList = [InfList, InList(1), SupList]
        end if
    end function QuickSort

end program STRUCTURE