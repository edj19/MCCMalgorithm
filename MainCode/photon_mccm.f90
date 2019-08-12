module constant
	implicit none
	real,parameter::Pi=3.1415926,c=3.0d8,Ed=0.019,me=9.11d-31,&
	                q=1.6d-16,mc=1.99d-26,Z=8,r0=2.8d-15,M=15.99,&
	                Nm=26.62
end module constant


module autoSimpson
	use constant
contains
	subroutine solve(func,s,a,b,tol,n)
		implicit real*8(a-z)
		external func
		integer::n,i
		n=40
		do i=1,20
			call simp(func,s,a,b,n)
			n=n*2
			call simp(func,s1,a,b,n)
			if(del<tol) exit
		end do
		s=s1
	end subroutine solve
	subroutine simp(func,s,a,b,n)
		implicit real*8(a-z)
		external func
		integer::n,k
		s=0d0
		h=(b-a)/n/2d0
		call func(f1,a)
		call func(f2,b)
		s=f1+f2
		call func(f1,a+h)
		s=s+4d0*f1
		do k=1,n-1
			t1=a+(2d0*k+1)*h
			t2=a+2d0*k*h
! 			print*,t1,t2
! 		     pause
			call func(f3,t1)
			call func(f4,t2)
			s=s+f3*4d0+f4*2d0
		end do
		s=s*h/3d0
	end subroutine simp
	subroutine fun1(f,E1)
		implicit real*8(a-z)
		real::Tm,b0,a0,T1,spka0,spka1,spka,g1,g2,I,Ee0,sigE,sigP,Tav,TavE,TavP
		real::dEdx1,dEdx
	    a0=Z/137
	    Ee0=me*c**2/q
	    M1=M*mc/12
	    Em0=M1*c**2/q
		a0=Z/137
		Tm=2*E1*(E1+2*Ee0)/Em0
		call vector(vec,Ee0,E1,me)
	    v=vec        
	    b0=v/c
        T1=Tm/Ed
!         call sigmaE(Tm,Ed,a0,b0,r0,Z,Ee0,sigE,Tav)
        call sigmaE_electron(Tm,Ed,a0,b0,r0,Z,Ee0,sigE,TavE)
        call sigmaE_positron(Tm,Ed,a0,b0,r0,Z,Ee0,sigP,TavP)
        spka0=T1-1-b0**2*Log(T1)+Pi*a0*b0*(2*(sqrt(T1)-1)-Log(T1))        
        Tav = TavE
        if (Tav<Ed) then
        	spka1=0
        	else if (Tav<2*Ed/0.8) then
        		spka1=spka0
        	else
!         		spka1=(T1-1-b0**2*Log(T1)+Pi*a0*b0*(2*(sqrt(T1)-1)-Log(T1)))*Tm/(2*Ed)
        		! T1=T1/2.5
        		! spka1=1.5-b0**2*Log(2.5)+Pi*a0*b0*(2*(sqrt(2.5)-1)-Log(2.5)) &
        		!        & +T1*Log(T1)-b0**2*(T1-1)+Pi*a0*b0*(T1-2*sqrt(T1)+1)
        		 spka1=spka0*0.8*Tav/Ed
        	end if  
!         spka1=T1-1-b0**2*Log(T1)+Pi*a0*b0*(2*(sqrt(T1)-1)-Log(T1))
        g2=1/(1-b0**2)
        g1=sqrt(g2)
        spka=Pi*Z**2*r0**2*spka1/(b0**4*g2)
        !I=0.00976*Z+0.0588*Z**(-0.19)
        I=0.013*Z
        dEdx1=Log(Ee0*b0*g2*E1/2*I**2)-(1+(2*g1-1)*Log(2.0)+(g1-1)**2/8)/g2
        dEdx=2*Pi*Nm*Z**2*r0**2*Ee0/(b0**2)*dEdx1
        f=Nm*spka/dEdx
        ! f=Nm*spka
!               print*,Ee0,Em0,v 
!               pause
	end subroutine fun1

	subroutine vector(vec,Ee0,E1,me)
		implicit real*8(a-z)
        real ::me,Ee0
        real,parameter::e=1.6d-16
   ! print*,Ee0,Ee
       real::p2,v2
       real,parameter::c=3.0d8
!    p2=(e*Ee)**2/(c**2)+2*(e*Ee)*me
!    v2=p2*c**2/(p2+Ee0*e*me)
     p2=(Ee0*e/(E1*e+Ee0*e))**2
     v2=(1-p2)*c**2
!       p2=((E1*e)**2-(Ee0*e)**2)/c**2
!       v2=p2/(p2+Ee0*me)
       vec=Sqrt(v2)
!      print*,Ee0,E1,me,p2,vec
!      pause

end subroutine vector
!   


! 	real function vector(Ee0,E1,me)
!     implicit real*8(a-z)
!     real,intent(in)::Ee0,E1
!     real,intent(in)::me
!     real,parameter::e=1.6d-16
!    ! print*,Ee0,Ee
!    real::p2,v2
!    real,parameter::c=3.0d8
! !    p2=(e*Ee)**2/(c**2)+2*(e*Ee)*me
! !    v2=p2*c**2/(p2+Ee0*e*me)
! !    p2=(Ee0*e/(Ee*e+Ee0*e))**2
! !    v2=(1-p2)*c**2
!    p2=((E1*e)**2-(Ee0*e)**2)/c**2
!    v2=p2/(p2+Ee0*me)
!     vector=Sqrt(v2)
! !     print*,Ee0,E1,me,p2,vector
! !     pause
! !     vector=Sqrt(2*Ee*e/me)
!  end function vector

end module autoSimpson

program photonvec
	use constant
	use autoSimpson
	implicit real*8(a-z)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!注意单位  能量：keV  电子静止质量：9.11d-31kg 元电荷：1.6d-19 碳原子质量：1.99d-26kg
	!          玻尔半径(a0)：5.29d-11m  电子经典半径：2.8d-15m
	!单位换算  1barn=10d-28 m2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real,parameter ::Ed=0.073,me=9.11d-31,q=1.6d-16,&
! 	                 mc=1.99d-26,Z=22.,r0=2.8d-15,c=3.0d8,&
! 	                 Pi=3.1415926
 	real :: Ep,Ee0,Em0,Ee1,T1,T,value
	real :: n1,P,P1,M1,sex,sigmaD,sigmaD1,sigmaD2,Tav,sigE,sigmaE_electron,sigmaE_positron
    real :: Tm  
    integer ::i,number,n
	real ::sigmaF,eglos        !function需要在主函数前进行声明
	real:: a0,b0,v,egd,g1,g2,spka,spka1
	open(15,file='photonvec.txt')
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!the displacement per atom cross section 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Do E1=1000,4000,10
	   a0=Z/137
	   Ee0=me*c**2/q
	   M1=M*mc/12
	   Em0=M1*c**2/q
! 	Tm=2*Ee*(Ee+2*Ee0)/Em0
	   Ec=sqrt(Ee0**2+Em0*Ed/2.0)-Ee0
	   if (E1<Ec) then
	   	     s=0
	   	else
	        call solve(fun1,s,Ec,E1,1d-7,n)
	     end if 
	   write(15,*)E1,s
! 	   102 format(T5,'total sequences:',I5,/,T5,'the results:',F18.10)
   end do 
end program photonvec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!相对论电子速度vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     real function vector(Ee0,E1,me)
!     implicit none
!     real,intent(in)::Ee0,E1
!     real,intent(in)::me
!     real,parameter::e=1.6d-16
!    ! print*,Ee0,Ee
!    real::p2,v2
!    real,parameter::c=3.0d8
! !    p2=(e*Ee)**2/(c**2)+2*(e*Ee)*me
! !    v2=p2*c**2/(p2+Ee0*e*me)
! !    p2=(Ee0*e/(Ee*e+Ee0*e))**2
! !    v2=(1-p2)*c**2
!    p2=((E1*e)**2-(Ee0*e)**2)/c**2
!    v2=p2/(p2+Ee0*me)
!     vector=Sqrt(v2)
!     print*,Ee0,E1,me,p2,vector
!     pause
! !     vector=Sqrt(2*Ee*e/me)
!  end function vector

subroutine sigmaE(T1,Ed,a0,b0,r0,Z,Ee0,sigE,Tav)
    ! use constant
	! calculate the displacement cross section and the average atom recoil kinetic energy Tav 
	implicit none
	real,intent(in)::T1,Ed,Ee0
	real,intent(in)::a0,b0,r0,Z 
	real,intent(out)::sigE,Tav
	real :: T2,T3,C1,sigE1
	real,parameter::e=1.6d-16,Pi=3.1415926
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


 subroutine sigmaE_electron(T1,Ed,a0,b0,r0,Z,Ee0,sigE,TavE)
!     use constant
	! calculate the displacement cross section and the average atom recoil kinetic energy Tav
	implicit none
	real,intent(in)::T1,Ed,Ee0
	real,intent(in)::a0,b0,r0,Z
	real,intent(out)::sigE,TavE
	real :: T2,T3,C1,sigE1
	real,parameter::e=1.6d-16,Pi=3.1415926
	T2=T1/Ed
	T3=2*(sqrt(T2)-1)-log(T2)
! 	C1=4*Pi*Z**2*e**4*(1-b0**2)/(b0**4*Ee0**2)
	C1=0.3136*(Z**2)*(1-b0**2)/(b0**4)
! 	sigmaE=C1*(T2-1-b0**2*log(T2)+Pi*a0*b0*T3)
     sigE=Pi*C1*(T2-1-(b0**2)*log(T2)+Pi*a0*b0*T3)/4
!     sigE=Pi*C1*(T2-1-(b0**2)*log(T2)+Pi*(a0/b0)*T3)/4
    sigE1=4*sigE/(Pi*C1)
    TavE=T1*(log(T2)-b0**2*(1-1/T2)+Pi*a0*b0*(sqrt(1/T2)-1)**2)/sigE1

    end subroutine sigmaE_electron

!   ! positron interaction cross section forms

 subroutine sigmaE_positron(T1,Ed,a0,b0,r0,Z,Ee0,sigP,TavP)
!     use constant
	! calculate the displacement cross section and the average atom recoil kinetic energy Tav
	implicit none
	real,intent(in)::T1,Ed,Ee0
	real,intent(in)::a0,b0,r0,Z
	real,intent(out)::sigP,TavP
	real :: T2,T3,C1,sigP1
	real,parameter::e=1.6d-16,Pi=3.1415926
	T2=T1/Ed
	T3=2*(sqrt(T2)-1)-log(T2)
! 	C1=4*Pi*Z**2*e**4*(1-b0**2)/(b0**4*Ee0**2)
	C1=0.3136*(Z**2)*(1-b0**2)/(b0**4)
! 	sigmaE=C1*(T2-1-b0**2*log(T2)+Pi*a0*b0*T3)
     sigP=Pi*C1*(T2-1-(b0**2)*log(T2)-Pi*a0*b0*T3)/4
!     sigE=Pi*C1*(T2-1-(b0**2)*log(T2)+Pi*(a0/b0)*T3)/4
    sigP1=4*sigP/(Pi*C1)
    TavP=T1*(log(T2)-b0**2*(1-1/T2)-Pi*a0*b0*(sqrt(1/T2)-1)**2)/sigP1

    end subroutine sigmaE_positron


 
