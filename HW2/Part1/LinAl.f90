module LinAl

contains

  !********************************************************

subroutine readMat(filename, mat, msize, nsize)

    implicit none
    character(len=*), intent(IN) :: filename
    real, allocatable, dimension(:,:), intent(OUT) :: mat
    integer, intent(OUT) :: msize,nsize
    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    ! Read Matrix Dimensions
    open(10,file=filename)
    read(10,*) msize,nsize
    allocate(mat(msize,nsize))
    mat = 0.0
    
    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)

end subroutine readMat



subroutine traceMat(mat,msize,c)

    implicit none
    
    real, dimension(:,:), allocatable, intent(IN) :: mat
    real, intent(OUT) :: c
    integer, intent(IN) :: msize
    integer :: i
  
    ! This routine calculates the trace, i.e. sum of diagonal elements, of an
    ! input msize by msize matrix, mat. The trace is returned in variable c.  


    c = 0 
  
    !sum diagonal elements into c
    do i = 1,msize
      c = c + mat(i,i)
    end do

end subroutine traceMat
  

subroutine stdNorm(v, vdim, vnorm)

    implicit none

    integer :: i 
    integer, intent(IN) :: vdim
    real, dimension(:), allocatable, intent(IN) :: v
    real, intent(OUT) :: vnorm

    ! This subroutine calculates the 2-norm of an input vector v, with dimension vdim.
    ! The output is the scalar vnorm.
    
    vnorm = 0
    do i=1,vdim
      vnorm = vnorm + v(i)**2
    end do

    vnorm = sqrt(vnorm)

end subroutine stdNorm


subroutine printMat(mat, msize, nsize)

    implicit none

    real, dimension(:,:), allocatable, intent(IN) :: mat
    integer, intent(IN) :: msize, nsize
    integer :: i,j
    
    ! This subroutine prints an input (msize) x (nsize) matrix, mat, in a readable manner

    print*,'dimensions:',msize, 'x',nsize
    do i= 1,msize
          write(*,*) mat(i,:)
    end do
    print*,''
end subroutine printMat

  
subroutine GEPP(A, B, m, n, singular)

    implicit none

    real, dimension(:,:), allocatable, intent(INOUT) :: A, B
    real, dimension(:), allocatable :: Atemp, Btemp
    integer, intent(IN) :: m, n
    logical, intent(OUT):: singular
    real :: factor
    integer :: i,j,k,l


    ! This subroutine performs Gaussian Elimination with partial pivoting on
    ! the matrix A augmented with matrix B = [b1 | b2 | ... | b_n].
    ! The bulk of this algorithm is followed step-by-step from the PDF lecture notes provided.
    ! Note: the user must specify B as a matrix, even if they want to solve Ax = b, where b is an (m x 1) vector.
    
    allocate(Atemp(m))
    allocate(Btemp(n))
    Atemp = 0.0

    do j = 1,m-1
      k = maxloc(abs(A(j:m, j)),1) + j-1 ! Index associated with largest abs value in column j, below the diagonal.
      if (k .NE. j) then !swap rows of A and B if needed
        !A swaps
        Atemp = A(k,:)
        A(k,:) = A(j,:)
        A(j,:) = Atemp

        !B swaps
        Btemp = B(k,:)
        B(k,:) = B(j,:)
        B(j,:) = Btemp

      end if

      if (A(j,j) == 0.0) then
        singular = .TRUE.
        stop 'Matrix is gonna explode, abort mission'
      end if

      do i = j+1,m
        factor = A(i,j)/A(j,j)
        A(i,:) = A(i,:) - (factor)*A(j,:)
        B(i,:) = B(i,:) - (factor)*B(j,:)
      end do
    end do
     
end subroutine GEPP
    

subroutine backsub(U,B,X, m, n)
      
      implicit none

      real,allocatable,dimension(:,:), intent(IN) :: U,B
      real,allocatable,dimension(:,:), intent(OUT):: X
      integer, intent(IN):: m,n
      integer :: i,j
      logical :: singular
      real,allocatable,dimension(:) :: sum

      !Performs backsubstitution on the system of the form UX = B, where U is an upper
      !triangular (m x m) matrix, X is the (m x n) solution matrix, and B is a matrix of n column m-vectors [b1 | b2 | ... | bn].
      !So, the first column of X is the solution to the system Ux1 = b1, etc.

      allocate(X(m,n))
      allocate(sum(n))
      X = 0.0
      sum = 0.0

      if (U(m,m) == 0) then
        stop 'Matrix is gonna explode, abort mission'
      end if

      X(m,:) = B(m,:)/U(m,m)
      do i = m-1,1,-1
        if (U(i,i) == 0) then
          singular = .TRUE.
          stop 'Matrix is gonna explode, abort mission'
        else
          singular = .FALSE.
        end if
        sum = 0.0
        do j = i+1,m
          sum = sum + U(i,j)*X(j,:)
        end do
        X(i,:) = (B(i,:) - sum)/U(i,i)
      end do

end subroutine backsub


subroutine LUdecomp(A, msize, singular, s)

      implicit none

      real,allocatable,dimension(:,:), intent(INOUT) :: A
      integer,intent(IN)::msize
      logical, intent(OUT) :: singular
      integer,allocatable,dimension(:), intent(OUT) :: s
      real,allocatable,dimension(:) :: Atemp
      integer :: Stemp
      integer :: i,j,k,z

      ! This subroutine performs an LU decomposition on input matrix A, returning
      ! LU where L is lower triangular and U is upper triangular.
      ! The algorithm is directly coded from pseudocode provided in lecture notes.

      allocate(Atemp(msize))
      allocate(s(msize))

      do j = 1,msize
        s(j) = j
      end do

      do j = 1,msize
        k = maxloc(abs(A(j:msize,j)), 1) + j-1 ! Index associated with largest abs value in column j, below the diagonal.
        if (k .NE. j) then
          !swap rows k and j of A
          Atemp = A(k,:)
          A(k,:) = A(j,:)
          A(j,:) = Atemp

          !swap elements s(k) and s(j) of s
          Stemp = s(k)
          s(k) = s(j)
          s(j) = Stemp
        end if
        if (A(j,j) == 0.0) then
          singular = .TRUE.
          stop 'Matrix is gonna explode, abort mission'
        else 
          singular = .FALSE.
        end if
        do i = j+1,msize
          A(i,j) = A(i,j)/A(j,j) !overwrites subdiagonal elts with values for L
          do z = j+1,msize
            A(i,z) = A(i,z) - (A(i,j)*A(j,z)) !updates A
          end do
        end do
      end do

end subroutine LUdecomp

subroutine LUbacksub(A, B, msize, s)

      implicit none

      real,allocatable,dimension(:,:),intent(IN):: A
      real,allocatable,dimension(:,:),intent(INOUT):: B
      integer,allocatable,dimension(:),intent(INOUT) :: s
      integer, intent(IN) :: msize
      real,allocatable,dimension(:,:) :: Y
      real,allocatable,dimension(:) :: sum
      integer:: i,j,k,nsize

      ! This subroutine solves the systems LUX = B. It does this first by solving
      ! LY = B, and then UX = Y.
      ! The algorithm follows from the in-class pseudo code.
      
      nsize = size(B,2)
      
      allocate(y(msize,nsize))
      allocate(sum(nsize))
      Y=0.0
      sum = 0.0

      do j=1,msize
        Y(j,:) = B(s(j),:)
      end do

      do j=1,msize-1
        do i = j+1,msize
          Y(i,:) = Y(i,:) - Y(j,:)*A(i,j)
        end do
      end do
      do i = msize,1,-1
        if (A(i,i)==0.0) then
          stop 'Matrix is gonna explode, abort mission'
        end if
        sum = 0.0
        do k = i+1,msize
          sum = sum + A(i,k)*B(k,:)
        end do
        B(i,:) = (Y(i,:)-sum)/A(i,i)
      end do

end subroutine LUbacksub

subroutine CholeskyDecomp(A,singular)
  implicit none

  real, allocatable, dimension(:,:), intent(INOUT) :: A
  logical, intent(INOUT) :: singular
  integer :: msize,i,j,k

  msize = size(A,1)

  do j = 1,msize
    do k = 1,j-1
      A(j,j) = A(j,j) - A(j,k)*A(j,k)
    end do
    A(j,j) = sqrt(A(j,j))
    do i = j+1,m
      do k=1,j-1
        A(i,j) = A(i,j) - A(i,k)*A(j,k)
      end do
      A(i,j) = A(i,j)/A(j,j)
    end do
  end do

end subroutine CholeskyDecomp



  end module LinAl
