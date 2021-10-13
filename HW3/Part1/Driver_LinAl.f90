Program Driver_LinAl

!(i) This driver is compiled via the Makefile submitted alongside it. To compile and run, simply type "make -f Makefile".
!    After the make command has been executed, the output of the program is stored in "output.txt".
!    To remove the files created by the make command, type "make -f Makefile clean".

!(ii) All the coding homework problem outputs are output from this single driver file. The inputs
!     and outputs are defined as necessary for the subroutines that are called from  the module LinAl.f90
!     to be called properly. 

  use LinAl

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j,msize,nsize,p
  real :: traceA, normVal, pi
  real,allocatable,dimension(:,:) :: A,B,L,R, squareR,shortB,Q,X, coeffMat, vandMat, gMat, curveVals,err,Id
  real,allocatable, dimension(:) :: v
  integer,allocatable, dimension(:) :: s
  logical :: singular

!Problem 1 (Cholesky stuff)--------------------------------------------------------------------
  print*, '===================================================================================='
  print*, '//////////////////////////////////PROBLEM 1\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
  print*, '===================================================================================='
  !Read in matrix of Atkinson data.

  myFileName = 'atkinson.dat'
  call readMat(myFileName,A,msize,nsize)

  p = 3 !--------> Change degree of interpolant HERE <--------

  !Construct Vandermonde matrix
  call vander(A, vandMat, p)
  allocate(coeffMat(p+1,1)) ! Allocation of coefficient vector A^Tb in normal equation A^TAx = A^Tb, and
  allocate(gMat(p+1,p+1))   ! Gram matrix A^TA.
  allocate(B(msize,1))
  allocate(err(msize,1)) !Initialize error vector as an m x 1 matrix
  allocate(v(msize)) 
  allocate(curveVals(msize,1)) ! A m x 1 matrix to store values from (vandMat)*gMat
  B(:,1) = A(:,2)
  gMat = matmul(transpose(vandMat), vandMat)
  coeffMat = matmul(transpose(vandMat), B)
  
  !Perform the Cholesky decomposition on the Gram matrix, (vandMat)^T(vandMat), and solve the respective system.
  call CholeskyDecomp(gMat,singular)
  !Returns the solution to the system, stored in coeffMat.
  call CholeskyBacksub(gMat, coeffMat)
  print*, 'Solution vector with degree:',p
  call printMat(coeffMat, p+1,1)

  !The vector vandMat*coeffMat results in the fitted curve values.
  curveVals = matmul(vandMat, coeffMat)
  
  !Calculates the 2-norm error between the fitted curve values and given data.
  err = B - curveVals
  v = err(:,1)
  call stdNorm(v, msize, normVal)
  print*,'2-norm error between curve and data: ',normVal

  !Print fitted curve data to file 'curve_Chol.dat'
  call fileMat(curveVals, 'curve_Chol.dat')

  !Deallocate our variables for later use in Problem 2.
  deallocate(A)
  deallocate(coeffMat)
  deallocate(gMat)
  deallocate(B)
  deallocate(err)
  deallocate(v)
  deallocate(curveVals)
  deallocate(vandMat)




!Problem 2 (QR stuff) --------------------------------------------------------------------------
  print*, '===================================================================================='
  print*, '//////////////////////////////////PROBLEM 2\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
  print*, '===================================================================================='
  
  p = 3 !--------> Change degree of interpolant HERE <--------

  !Read in atkinson data and construct the Vandermonde matrix interpolant degree 3.
  myFileName = 'atkinson.dat'
  call readMat(myFileName,A,msize,nsize)
  call vander(A,vandMat,p)
  
  !Construct identity matrix for later use in computing Q^TQ - I.
  allocate(Id(msize,msize))
  Id = 0.0
  do i = 1,msize
    Id(i,i) = 1.0
  end do

  !Allocation of RHS vector B, the error matrix err = A - QR, a truncated R and B, named "squareR" and "shortB" 
  !respectively. We did this to solve (squareR)x = Q^T(shortB). Also allocate curveVals to store fitted curve values.
  allocate(B(msize,1))
  allocate(err(msize,p+1))
  allocate(X(p+1,1))
  allocate(squareR(p+1,p+1))
  allocate(shortB(p+1,1))
  allocate(curveVals(msize,1))
  allocate(v(msize))
  B(:,1) = A(:,2)
  
  !Perform QR decomposition on the Vandermonde matrix, resulting in Q and R
  call HHQR(vandMat,Q,R)
  err = vandMat - matmul(Q,R) !Error matrix vandMat - QR
  print*, 'Printing A - QR:'
  call printMat(err,msize,p+1)
  call FrobeniusNorm(err, normVal)
  print*, 'Frobenius norm of A - QR: ', normVal
  deallocate(err)

  allocate(err(msize,msize))
  err = 0.0
  err = matmul(transpose(Q), Q) - Id
  print*,'Printing Q^TQ - I:'
  call printMat(err, msize,msize)
  call FrobeniusNorm(err, normVal)
  print*, 'Frobenius norm of Q^TQ - I: ',normVal

  !RHS of Rx = Q^Tb
  B = matmul(transpose(Q),B)
  
  !Construct squareR from R, and shortB from B.
  do i = 1,p+1
    squareR(i,:) = R(i,:)
  end do

  do i = 1,p+1
    shortB(i,1) = B(i,1)
  end do

  !Set B back to being the second column of atkinson.dat.
  B(:,1) = A(:,2)

  print*, 'Printing R:'
  call printMat(squareR,p+1,p+1)
  call backsub(squareR, shortB, X, p+1, 1)
  print*, 'Solution vector for degree:',p
  call printMat(X,p+1,1)

  !Computes fitted curve values and prints them to file 'curve_QR.dat'.
  curveVals = matmul(vandMat, X)
  call fileMat(curveVals, 'curve_QR.dat')

  err = B - curveVals
  v = err(:,1)
  call stdNorm(v, msize, normVal)
  print*,'2-norm error between curve and data: ',normVal
  


End Program Driver_LinAl
 