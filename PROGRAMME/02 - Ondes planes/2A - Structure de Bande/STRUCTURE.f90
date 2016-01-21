!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROGRAMME PRINCIPAL
program STRUCTURE
  implicit none

  integer,parameter              :: Rk =  8

  real(kind=8) :: a,f,Vc,eps1,eps2,c=299792458,pi=3.1415926535, normeG, eps, nnum1,nnum2,nnum,epsm1,epsm2,&
                   epsm,nth, R
  real(kind=8) :: e1(3),e2(3),e3(3),a1(3),a2(3),a3(3),b1(3),b2(3),b3(3), G1(3),eta1(3),eta2(3),eta3(3),eta4(3)
  real(kind=8) :: kX(3),kU(3),kL(3),kG(3),kW(3),kK(3),k0(7,3)
  real(kind=8),allocatable :: G2(:,:), G(:,:), k(:,:),h(:,:), Vp(:,:), Vp1(:), kn(:), vcp(:,:), w(:,:),fv1(:),fv2(:)
  integer                  :: Nmax = 0, N = 0, Gmax, i=0, nk, i1, i2, i3, ng, ng1,m, size1, size2, info, j
  integer, allocatable     :: IWORK(:)
  real(kind=8), allocatable :: WORK2(:)
  complex(kind=8), allocatable :: WORK(:)
  character(len=1) :: VECT = "N", STOC = "U";

    open(10,file = 'INPUT.dat')
    read(10,*) a, f, eps1, eps2, Gmax, nk
    close(10)

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

    kX=pi/a*[0.,2.,0.];
    kU=pi/a*[1./2,2.,1./2];
    kL=pi/a*[1.,1.,1.];
    kG=pi/a*[0.,0.,0.];
    kW=pi/a*[1.,2.,0.];
    kK=pi/a*[3./2,3./2,0.];

    k0=transpose(reshape([kX,kU,kL,kG,kX,kW,kK],[3,7]))

    allocate(k(nk+1,3))

    do i1 = 0,5
        i=i1*nk/6+1;
        k(i,:)=k0(i1+1,:);
        do i2=1,nk/6-1
            i=i+1;
            k(i,:)=k0(i1+1,:)+i2*(k0(i1+2,:)-k0(i1+1,:))/(nk/6.);
        enddo
    enddo

    k(nk+1,:)=k0(7,:);

    allocate(H(2*nG,2*nG),Vp1(2*nG),Vp(nk+1,10),kn(nk+1),WORK(6*nG), vcp(2*ng,2*ng),fv1(2*ng),fv2(2*ng))
    allocate(WORK2(6*nG), IWORK(1))

    size1 = 2*nG
    size2 = 6*nG

    do j=1,nk+1
        m=1;
        do i1=1,2*nG-1,2
            if(all((k(j,:)+G(m,:))==0)) then
                eta1=e1;
                eta2=e2;
            elseif (abs(k(j,1)+G(m,1))<abs(k(j,2)+G(m,2))) then
                eta1=cross(e1,(k(j,:)+G(m,:)));
                eta2=cross(eta1,(k(j,:)+G(m,:)));
            else
                eta1=cross(e2,(k(j,:)+G(m,:)));
                eta2=cross(eta1,(k(j,:)+G(m,:)));
            endif
             
            eta1=eta1/(norm2(eta1));
            eta2=eta2/(norm2(eta2));

            n=1;
            do i2=1,2*nG-1,2
                if (all((k(j,:)+G(n,:))==0)) then
                    eta3=e1;
                    eta4=e2;
                elseif (abs(k(j,1)+G(n,1))<abs(k(j,2)+G(n,2))) then
                    eta3=cross(e1,(k(j,:)+G(n,:)));
                    eta4=cross(eta3,(k(j,:)+G(n,:)));
                else
                    eta3=cross(e2,(k(j,:)+G(n,:)));
                    eta4=cross(eta3,(k(j,:)+G(n,:)));
                endif

                eta3=eta3/(norm2(eta3));
                eta4=eta4/(norm2(eta4));

                normeG=norm2(G(m,:)-G(n,:));

                if (normeG==0) then
                    eps=1./eps1+(1./eps2-1./eps1)*f;
                else
                    eps=(1./eps2-1./eps1)/Vc*(4.*pi*(sin(normeG*R)-normeG*R*cos(normeG*R)))/(normeG**3.);
                endif

                H(i1,i2)=eps*dot_product(cross(k(j,:)+G(m,:),eta1),cross(k(j,:)+G(n,:),eta3));
                H(i1,i2+1)=eps*dot_product(cross(k(j,:)+G(m,:),eta1),cross(k(j,:)+G(n,:),eta4));
                H(i1+1,i2)=eps*dot_product(cross(k(j,:)+G(m,:),eta2),cross(k(j,:)+G(n,:),eta3));
                H(i1+1,i2+1)=eps*dot_product(cross(k(j,:)+G(m,:),eta2),cross(k(j,:)+G(n,:),eta4));

                n=n+1;
            enddo
            m=m+1;
        enddo
        call dsyev('V','U', size1, H, size1, Vp1, WORK, size2, INFO)
        Vp(j,:)=Vp1(1:10);
    enddo

    kn = norm2(k,2)*a/(2*pi)

    w=reshape(sqrt(Vp)/(2*pi),[size(Vp),1]);

    nnum1=(kn((3*nk)/5)-kn(nk/2))/(w((3*nk)/5,1)-w(nk/2,1))
    nnum2=(kn((2*nk)/5)-kn(nk/2))/(w((2*nk)/5,1)-w(nk/2,1))
    nnum=(nnum1+nnum2)/2

    epsm1=f*eps2+(1-f)*eps1;
    epsm2=1/(f/eps2+(1-f)/eps1);
    epsm=(epsm1+epsm2)/2;
    nth=sqrt(epsm)

    write(*,*) nth, nnum

    open(100,file='OUTPUT.dat')
    w = sqrt(Vp)/(2*pi);
    write(100,'(I3,10E13.5)') (i-1,w(i,1:10),i=1,nk+1)
    close(100)

    deallocate(k,H,Vp1,Vp,kn,WORK,G2,G)

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