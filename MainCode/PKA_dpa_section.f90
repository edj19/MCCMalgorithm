module constant
	implicit none
	real,parameter::Pi=3.14159,c=3.0d8
end module constant

program photon
	use constant
	implicit real*8(a-z)
	real,parameter ::me=9.11d-31,q=1.6d-16,&
	                 mc=1.99d-26,Z=22.,r0=2.8d-15,Ed=0.025
!  parameter,Ed
	real :: Ep,Ee0,Em0,Ee,Ee1,T1,T,value,M=47.87,Ec
	real :: n1,P,P1,M1,sex,sigmaD,sigmaD1,sigmaD2,Tav,sigE
        
    integer ::i,number
	real ::vector,sigmaF,eglos        !function需要在主函数前进行声明
	real:: a0,b0,v,egd
	open(10,file="photon.txt")
! 	call solve(fun1,s,0d0,2d0,1d-7,n)
 !
!   Ep=1250
! do while(Ep<600)
!  		Ep=Ep+10
 do Ee=20,6000,100       !电子的能量分布
    Ee0=me*(c**2)/q
!    	Ee=Ep/(1+Ee0/(2*Ep))  !发生正碰撞时光子传递给电子的最大能量
	a0=Z/137
	v=vector(Ee0,Ee,me)
!	print*,v
! 	pause
 	b0=v/c
	M1=mc*M/12
    Em0=mc*M*c**2/(12*q)
    T1=2*Ee*(Ee+2*Ee0)/Em0  !发生正碰撞电子传输给原子的能量
!     Ec=Sqrt(Ee0**2+Em0*Ed/-Ee0
    egd=eglos(r0,Ee0,b0,Z,Ee)
    !    print*,T1
 !   pause
!    inter: do j=0,180,5
! 			theta=j*Pi/180
! 			T=T1*sin(theta/2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!位移取代原子截面的计算
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

      call sigmaE(T1,Ed,a0,b0,r0,Z,Ee0,sigE,Tav)
!            sigmaF0=-sigmaF(T1,Ed,a0,b0,Z,Ee0)
!积分结果，详细结果参见毕业论文    
!     write(10,*)Ee,sigmaE(T1,Ed,a0,b0,r0,Z)
    if(Tav<=Ed)then
    	sigmaD=0
    	else if(Tav>Ed.and.Tav<=2*Ed/0.8)then
      		call sigd1(T1,Ed,a0,b0,Z,sigmaD1)
         		    sigmaD=sigE
          		    ! sigmaD=sigE
      	else
      		call sigd2(T1,Ed,a0,b0,Z,sigmaD2)
                	! sigmaD=sigmaD2
                    sigmaD=0.8*sigE*Tav/(2*Ed)
     	end if
      	write(10,*)Ee,sigE,sigmaD
!      	write(10,*)Ee,T1,Tav
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!计算PKA产生的离位原子数目
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if(T1<=Ed)then
!     	number=0
!     	else if(T1>Ed.and.T1<=2*Ed/0.8)then
!     		number=1	
!      	else
!       		number=0.8*T1/(2*Ed)
!      	end if
!     	write(10,*)T1,Ep,number
! !     end do inter 
    
end do 
	end program photon

!!!!!!!!!!!!!!!!!!!!!!!!
!分布均匀的电子能量的平均值

subroutine qsf(Ep,Ee0,Ee)
	implicit none
	real,intent(in)::Ep,Ee0
	real,intent(out):: Ee
	integer :: j,counter
	real :: theta,a
	real,parameter::Pi=3.14159
	counter=0
	Ee=0.0
	do j=0,180,1
			theta=j*Pi/180
			counter=counter+1
			a=Ep**2*(1-cos(theta))
			a=a/(Ee0+Ep*(1-cos(theta)))
			Ee=Ee+a
		end do
		Ee=Ee/counter
	end subroutine qsf
	!---------------------------------------------------
	!sigd1  在第一个区间积分所得位移取代截面的大小
	!sigd2  在第二区间积分所得位移取代截面的大小
	!---------------------------------------------------

	subroutine sigd1(T1,Ed,a0,b0,Z,sigmaD1)
		use constant
		implicit none
		real,intent(in)::T1,Ed
		real,intent(in)::a0,b0,Z
		real,intent(out)::sigmaD1
		real::T4,T5,C3
		T4=T1/Ed
		T5=2*(sqrt(T4)-1)-log(T4)
		C3=0.3136*(Z**2)*(1-b0**2)/(b0**4)
         sigmaD1=Pi*C3*(T4-(b0**2)*log(T4)-1+Pi*a0*b0*T5)/4
!         sigmaD1=Pi*C3*(T4-(b0**2)*log(T4)-1+Pi*(a0/b0)*T5)/4
    end subroutine sigd1
 

    !-----------------------------------------------

    subroutine sigd2(T1,Ed,a0,b0,Z,sigmaD2)
		use constant
		implicit none
		real,intent(in)::T1,Ed
		real,intent(in)::a0,b0,Z
		real,intent(out)::sigmaD2
		real::T4,T5,C3
		T4=T1/(2.5*Ed)
		T5=T4-2*Sqrt(T4)+1
		C3=0.3136*(Z**2)*(1-b0**2)/(2*b0**4)
 		sigmaD2=Pi*C3*(T4*log(T4)-b0**2*(T4-1)+Pi*a0*b0*T5)/4+Pi*C3*(1.5-b0**2*log(2.5)+Pi*a0*b0*(2*sqrt(2.5)-2-log(2.5)))/4
!         sigmaD2=Pi*C3*(T4*log(T4)-b0**2*(T4-1)+Pi*(a0/b0)*T5)/4+Pi*C3*(1.5-b0**2*log(2.5)+Pi*(a0/b0)*(2*sqrt(2.5)-2-log(2.5)))/4
    end subroutine sigd2



!---------------------------------------------------
!PKA cross section in Original paper
!--------------------------------------

 subroutine sigmaE(T1,Ed,a0,b0,r0,Z,Ee0,sigE,Tav)
    use constant
	! calculate the displacement cross section and the average atom recoil kinetic energy Tav 
	implicit none
	real,intent(in)::T1,Ed,Ee0
	real,intent(in)::a0,b0,r0,Z 
	real,intent(out)::sigE,Tav
	real :: T2,T3,C1,sigE1
	real,parameter::e=1.6d-16
	T2=T1/Ed
	T3=2*(sqrt(T2)-1)-log(T2)
! 	C1=4*Pi*Z**2*e**4*(1-b0**2)/(b0**4*Ee0**2)
	C1=0.3136*(Z**2)*(1-b0**2)/(b0**4)
! 	sigmaE=C1*(T2-1-b0**2*log(T2)+Pi*a0*b0*T3)
     sigE=Pi*C1*(T2-1-(b0**2)*log(T2)+Pi*a0*b0*T3)/4
!     sigE=Pi*C1*(T2-1-(b0**2)*log(T2)+Pi*(a0/b0)*T3)/4
    sigE1=4*sigE/(Pi*C1)
    Tav=T1*(log(T2)-b0**2*(1-1/T2)+Pi*a0*b0*(sqrt(1/T2)-1)**2)/sigE1

    end subroutine sigmaE
    
   
!相对论电子速度vector

    real function vector(Ee0,Ee,me)
    implicit none
    real,intent(in)::Ee0,Ee
    real,intent(in)::me
    real,parameter::e=1.6d-16
   ! print*,Ee0,Ee
   real::p2,v2
   real,parameter::c=3.0d8
!    p2=(e*Ee)**2/(c**2)+2*(e*Ee)*me
!    v2=p2*c**2/(p2+Ee0*e*me)
   p2=(Ee0*e/(Ee*e+Ee0*e))**2
   v2=(1-p2)*c**2
    vector=Sqrt(v2)
!     vector=Sqrt(2*Ee*e/me)
 end function vector

  !-------------------------------------------
  !ao-Boher radius
  ! Er=ε2/2a0=13.6eV -Rydberg energy
  !the cross section of PkA from the book--the fundamentals
  !-----------------------------------------------

  real function sigmaF(T1,Ed,a0,b0,Z,Ee0)
  implicit none
  real,intent(in)::T1,Ed,Ee0
  real,intent(in)::a0,b0,Z
  real,parameter::Er=0.0136,ao=5.29d-11,Pi=3.14159
  real::T2,C1,T3
  T2=T1/Ed
  C1=4*Pi*ao**2*Z**2*Er**2/(Ee0**2)
  T3=(1-b0**2)*(T2-1)/(b0**4)
  sigmaF=C1*T3-log(T2)*b0**2+2*a0*b0*Sqrt(T2)-1-log(T2)
end function sigmaF



  !----------------------------------------------
  ! the electron stopping power,calculated
  !following the expression given by Bethe and Ashkin 
  !-------------------------------------------------

Real function eglos(r0,Ee0,b0,Z,Ee)
implicit none
real,intent(in)::r0,Ee0,b0
real,intent(in)::Z,Ee
real::Zm,I,ga,ga1,E1,C2
real,parameter::Pi=3.14159,Na=4.23,q=1.6d-16
I=9.76*Z+58.8*(Z**(-0.19))
C2=2*Pi*Na*Z**2*r0**2*Ee0*q/(b0**2)
ga=sqrt(1/(1-b0**2))
ga1=(1-b0**2)*(1+(2*ga-1)*log(2.0)+(ga-1)**2/8)
E1=log(Ee0*q*b0**2*ga**2*Ee/(2*I**2))
eglos=C2*(E1-ga1)
! print*,I,C2,ga,ga1,E1
! pause

end function eglos

			
