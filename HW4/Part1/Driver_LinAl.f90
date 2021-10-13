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
  integer :: i,j,msize,nsize
  real,allocatable,dimension(:,:) :: A1,A2,A3, eigvec
  real,allocatable, dimension(:) :: v
  integer,allocatable, dimension(:) :: s
  logical :: singular

  

!Problem 1 ==========================================================
  print*, '===================================================================================='
  print*, '//////////////////////////////////PROBLEM 1\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
  print*, '===================================================================================='

    myFileName = 'A1_mat.dat'

    call readMat(myFileName, A1, msize, nsize)

    print*,'Printing matrix: '
    call printMat(A1,msize,msize)

    call SymmToTri(A1)
    call printMat(A1,msize,msize)
    
  
!Problem 2 ==========================================================
    print*, '===================================================================================='
    print*, '//////////////////////////////////PROBLEM 2\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
    print*, '===================================================================================='

  !(i) Approximates eigenvalues of A via QR algorithm without shift.
     myFileName = 'A2_mat.dat'

     call readMat(myFileName,A2,msize,nsize)

    print*,'Printing matrix: '
    call printMat(A2,msize,msize)
    
    call QR_eigenvalues(A2)
    call printMat(A2,msize,msize)

    print*, 'The QR algorithm without shift eigenvalues: '
    do i = 1,msize
      print*, A2(i,i)
    end do

  !(ii) Approximates eigenvalues of A via QR algorithm with shift and deflation.

    call readMat(myFileName,A2,msize,nsize)

    print*,'Printing matrix: '
    call printMat(A2,msize,msize)
    
    call QR_eigenvalues_shift_deflation(A2)
  

    print*, 'The QR algorithm without shift & deflation eigenvalues: '
    do i = 1,msize
      print*, A2(i,i)
    end do
    

!Problem 3 ==========================================================
    print*, '===================================================================================='
    print*, '//////////////////////////////////PROBLEM 3\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
    print*, '===================================================================================='

    myFileName = 'A3_mat.dat'

    call readMat(myFileName,A3,msize,nsize)

    allocate(v(msize))
    v = (/-8.0286 ,7.9329,5.6689,-1.5731/)

    print*,'Printing matrix:'
    call printMat(A3,msize,nsize)

    do i = 1,size(v)
      call Inverse_Iteration(A3,v(i),eigvec)
      print*, 'Eigenvector associated with approximate eigenvalue', v(i)
      call printMat(eigvec, size(v),1)
    end do


End Program Driver_LinAl
 