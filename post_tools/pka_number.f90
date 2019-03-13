subroutine pka_number(Z,M,number)
	implicit none
	real,parameter ::me=9.11d-31,q=1.6d-16,&
	                 mc=1.99d-26,r0=2.8d-15,Ed=0.073,T1=0.5
	real :: Z,M
	real :: number
	!f2py intent(in) Z,M
	!f2py intent(out) number
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculation the PKA numbers of defects
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(T1<=Ed)then
    	number=0
    else if(T1>Ed.and.T1<=2*Ed/0.8)then
    	number=1
    else
      	number=0.8*T1/(2*Ed)
    end if
end subroutine 
