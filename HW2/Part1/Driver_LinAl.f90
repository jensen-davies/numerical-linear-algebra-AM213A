Program Driver_LinAl

!(i) This driver is compiled via the Makefile submitted alongside it. To compile and run, simply type "make -f Makefile".
!    After the make command has been executed, the output of the program is stored in "output.txt".
!    To remove the files created by the make command, type "make -f Makefile clean".

!(ii) All the coding homework problem outputs are output from this single driver file. The inputs
!     and outputs are defined as necessary for the subroutines that are called from  the module LinAl.f90
!     to be called properly. 

  use LinAl, only: readMat, traceMat, stdNorm, printMat, GEPP, backsub, LUdecomp, LUbacksub

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j,msize,nsize
  real :: traceA, normVal, pi
  real,allocatable,dimension(:,:) :: A,B,X, As, Bs, E, Alu, Blu, U, L
  real,allocatable,dimension(:,:) :: pMat, pVec, pSol !matrix for system for Problem 5, Part 1
  real,allocatable, dimension(:) :: v
  integer,allocatable, dimension(:) :: s
  logical :: singular

!Read in matrix A-----------------------------------
  myFileName = 'Amat.dat'

  call readMat(myFileName, A, msize, nsize)

!Warm Up Exercises Start Here----------------------
  allocate(v(msize))
  print*, 'Original A:'
  call printMat(A,msize,nsize)
  call traceMat(A,msize, traceA)
  print*, 'Trace of A is ',traceA
  do i = 1,nsize
    v = A(:,i)
    call stdNorm(v, msize, normVal)
    write(*,*) 'For column ',i,' the 2-norm is:', normVal
  end do


  !Allocating matrices As (for error) and Alu (for LU).
  allocate(As(msize,nsize))
  allocate(Alu(msize,nsize))
  !Saving a copy of A for computing error and calling LUdecomp later.
  As = A 
  Alu = A
  print*, ''

!Read in matrix B, and allocate lots of vars-------------------
  myFileName = 'Bmat.dat'

  call readMat(myFileName, B, msize, nsize)
  print*,'Original B:'
  do i = 1, msize
    write(*,*) (B(i,j) , j = 1, nsize)
  end do
  !Allocating the numerous arrays we pass to our subroutines in LinAl.
  allocate(Bs(msize,nsize))
  allocate(E(msize,nsize))
  allocate(Blu(msize,nsize))
  allocate(U(msize,msize))
  allocate(L(msize,msize))

  !Saving a copy of B for computer error and calling LUdecomp later.
  Bs = B 
  Blu = B
  print*,''

!3. Gaussian Elimination with Partial Pivoting
  !Performs Gaussian elimination on AX=B
  call GEPP(A,B,msize,nsize,singular)
  print*,'Row reduced A:'
  !A is now row reduced, and B transformed
  call printMat(A,msize,nsize)
  print*,'Transformed B:'
  call printMat(B,msize,nsize)
  print*,'Solution Matrix X from GE:'
  call backsub(A,B,X,msize,nsize)
  call printMat(X,msize,nsize)

  print*, 'Error matrix:'
  !calculates error matrix E, defined as E = AsX - Bs
  E = matmul(As,X) - Bs
  call printMat(E,msize,nsize)

  !prints the 2-norm normVal of each column of the error matrix E
  print*, 'Norm values of columns of E:'
  do i = 1, nsize
    v = E(:,i)
    call stdNorm(v,msize,normVal)
    print*,normVal
  end do
   
print*, '============================================='
print*, '--------------LU DECOMPOSITION--------------- '
print*, '============================================='

  Print*, 'Matrix A prior to LU decomposition:'
  call printMat(As,msize,msize)
!4. LU decomposition with Partial Pivoting
  print*, 'LU decomposition of A (in one matrix):'
  call LUdecomp(As,msize,singular,s)
  call printMat(As,msize,msize)

  print*, 'Upper triangular matrix U:'
  
  do i = 1,msize
    do j = 1,msize
      if (j < i) then
        U(i,j) = 0
      else
        U(i,j) = As(i,j)
      end if
    end do
  end do

  call printMat(U,msize,msize)

  print*, 'Lower Triangular matrix L:'

  do i = 1,msize
    do j = 1,msize
      if (i == j) then
        L(i,i) = 1
      else if (j > i) then
        L(i,j) = 0
      else
        L(i,j) = As(i,j)
      end if
    end do
  end do

  call printMat(L,msize,msize)
  call LUbacksub(As,Bs,msize,s)
  print*, 'Solution matrix X from LU decomp:'
  call printMat(Bs, msize, msize)

  !calculates error matrix E, defined as E = AsX - Bs
  E = matmul(Alu,Bs) - Blu
  print*,'Error matrix:'
  call printMat(E, msize,nsize)

  print*, 'Norm values of error matrix E columns:'
  do i = 1, nsize
    v = E(:,i)
    call stdNorm(v,msize,normVal)
    print*,normVal
  end do

!5. This point marks the beginning of problem 5, "A very basic application."

  ! We're given three points, A = (1,2,3), B = (-3,2,5), C = (pi, e, -\sqrt(2)).
  ! In general, if we have 3 arbitrary points that aren't collinear, this problem
  ! boils down to finding the normal vector, (a,b,c) of the plane
  !                   ax + by + cz + d = 0.
  ! Let d = -1 (which we can do, assuming the plane doesn't pass through the origin.)
  ! Hence, we can solve the linear system
  !   [x_1   y_1   z_1] [a]  = [1]
  !   [x_2   y_2   z_2] [b]  = [1]
  !   [x_3   y_3   z_3] [c]  = [1]
  ! For the direction of the plane, (a,b,c).

  allocate(pMat(3,3))
  allocate(pVec(3,1))
  pVec = 0.0
  pVec(:,1) = 1.0
  pMat = 0.0
  pi = acos(-1.0)

  !definine matrix for the above system
  pMat = reshape((/1.0,2.0,3.0,-3.0,2.0,5.0,pi,exp(1.0),-sqrt(2.0)/), shape(pMat))
  pMat = transpose(pMat)

  call GEPP(pMat, pVec, 3,1, singular)
  call backsub(pMat, pVec, pSol, 3,1)
  print*, 'The solution vector is:'
  call printMat(pSol,3,1)


  End Program Driver_LinAl
 