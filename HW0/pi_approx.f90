program pi_approximation

    implicit none

    real :: threshold1,threshold2,threshold3,threshold4, diff, pi_true, pi_aprx,n
    real, dimension(4) :: thresholds
    integer :: i
    
    n=0
    pi_aprx = 0
    threshold1 = 1.e-4
    threshold2 = 1.e-8
    threshold3 = 1.e-12
    threshold4 = 1.e-16 
    thresholds = (/threshold1,threshold2,threshold3,threshold4/)
    pi_true = acos(-1.d0)
    diff = abs(pi_true - pi_aprx)
   
    do i = 1,4
        do while (diff > thresholds(i))
            pi_aprx = pi_aprx + (16**(-1*n))*(4/(8*n+1) - (2/(8*n+4)) -(1/(8*n+5)) - (1/(8*n+6)))
            n = n+1
            diff = abs(pi_true - pi_aprx)
        end do

        print*, "--------------------------------------------"
        print*, "For Threshold",i,":"
        print*, "Approximation of pi:",pi_aprx
        print*, " diff =",diff
        print*, "N =",N
    end do

end program pi_approximation