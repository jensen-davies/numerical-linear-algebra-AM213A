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

subroutine fileMat(mat, filename)
  implicit none

  real, allocatable, dimension(:,:), intent(IN) :: mat
  character(len=*), intent(IN) :: filename
  integer:: msize,nsize,i,j

  !Writes the matrix mat to a file, with filename specified via user input.

  msize = size(mat,1)
  nsize = size(mat,2)

  !Write matrix to file 'filename'.
  open(11, file=filename)
  do i=1,msize
    write(11,*) ( mat(i,j), j=1,nsize)
  end do
  close(11)

end subroutine fileMat

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

subroutine FrobeniusNorm(mat, normVal)   
  implicit none

  real,allocatable,dimension(:,:),intent(IN) :: mat
  real, intent(OUT) :: normVal
  integer::i,j,msize,nsize

  !Calculates the Frobienius norm of the input matrix, mat. Outputs the value of the norm, normVal.

  msize = size(mat,1)
  nsize = size(mat,2)
  normVal = 0.0

  !Calculate the Frobenius norm of A.
  do i = 1,msize
    do j = 1,nsize
    normVal = normVal + mat(i,j)**2
    end do
  end do

  normVal = sqrt(normVal)

end subroutine FrobeniusNorm

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

      if (U(m,m) == 0.0) then
        stop 'Matrix is gonna explode, abort mission'
      end if

      X(m,:) = B(m,:)/U(m,m)
      do i = m-1,1,-1
        if (U(i,i) == 0.0) then
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
      
      allocate(Y(msize,nsize))
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

subroutine CholeskyDecomp(A,pd)
  implicit none

  !Performs the Cholesky decomposition on input matrix A. The second input, pd, indicates whether or not the matrix is
  !positive definite, or singular.

  real, allocatable, dimension(:,:), intent(INOUT) :: A
  logical, intent(INOUT) :: pd
  integer :: msize,i,j,k

  msize = size(A,1)
  pd = .TRUE.

  do j = 1,msize
    !Calculate the new diagonal elements
    do k = 1,j-1
      A(j,j) = A(j,j) - A(j,k)*A(j,k)
    end do
    !Check if it's positive definite/singular
    if (A(j,j) <= 1.0E-15) then
    pd = .FALSE.
    stop 'Input matrix is either singular or not positive definite.' 
    end if
    A(j,j) = sqrt(A(j,j)) 
    !Calculate the elements below the diagonal
    do i = j+1,msize
      do k=1,j-1
        A(i,j) = A(i,j) - A(i,k)*A(j,k)
      end do
      A(i,j) = A(i,j)/A(j,j)
    end do
  end do

end subroutine CholeskyDecomp

subroutine CholeskyBacksub(A,B)
  implicit none
  real, allocatable, dimension(:,:), intent(INOUT) :: A, B
  real,allocatable,dimension(:,:) :: Y
  real,allocatable,dimension(:) :: sum
  integer :: msize,nsize,i,j,k

  !This subroutine performs Cholesky backsubstitution on the system AX = B, where both A and B are inputs and, once modified, outputs.

  msize = size(A,1)
  nsize = size(B,2)
      
  allocate(Y(msize,nsize))
  allocate(sum(nsize))
  Y=B
  sum = 0.0

  !Forward substitution solving Ly = b
  do i = 1,msize
    sum = B(i,:)
    do j = 1,i-1
      sum = sum - Y(j,:)*A(i,j)
    end do
    Y(i,:) = sum/A(i,i)
  end do

  !Backward substitution solving L*x = y
  do i = msize,1,-1
    if (A(i,i) <= 1.0E-15) then
      stop 'Singular, abort'
    end if
    do k = i+1,msize
      Y(i,:) = Y(i,:) - A(k,i)*B(k,:)
    end do
    B(i,:) = Y(i,:)/A(i,i)
  end do

end subroutine CholeskyBacksub

subroutine HHQR(A,Q,R)
  implicit none

  real,allocatable,dimension(:,:), intent(IN):: A
  real,allocatable,dimension(:,:), intent(OUT) :: R,Q
  real,allocatable,dimension(:,:) :: V, outerprod,Id
  real,allocatable,dimension(:) :: s
  real :: normA,normV
  integer :: msize, nsize, i,j,k

  !Calculates the full QR decomposition of input matrix A, outputting an msize by msize orthogonal matrix Q, and
  !padded upper triangular matrix R, of dimension msize by nsize. The algorithm below follows directly from
  !pseudocode provided in the Chapter 3 lecture notes.

  msize = size(A,1)
  nsize = size(A,2)
  allocate(s(nsize))
  allocate(V(msize,nsize))
  allocate(outerprod(msize,msize))
  allocate(R(msize,nsize))
  allocate(Q(msize,msize))
  allocate(Id(msize,msize))
  
  !Construct identity matrix for later computing the recursive product that results in Q.
  Id = 0.0
  do i = 1,msize
    Id(i,i) = 1.0
  end do
  
  !Initialize our norms, Q, and V. We set R = A because A is purely an input value, and R is output.
  Q = Id
  normA = 0.0
  normV = 0.0
  V = 0.0
  R = A
  

  do j = 1,nsize
    !Calculates the 2-norm of A.
    do i = j,msize
      normA = normA + R(i,j)**2
    end do
    normA = sqrt(normA)
    
    s(j) = sign(normA, R(j,j))
    !Sets the appropriate value for V(j,j) prior to entering next do-loop.
    V(j,j) = R(j,j) + s(j)
    
    !Assign the remaining values for V below the diagonal entry.
    do k = j+1,msize
      V(k,j) = R(k,j)
    end do

    !Calculates the 2-norm of V.
    do i = j,msize
      normV = normV + V(i,j)**2
    end do
    normV = sqrt(normV)
    
    V(:,j) = V(:,j)/normV

    !Constructs the outer-product matrix vv^T
    do i = 1,msize
      do k = 1,msize
      outerprod(k,i) = V(k,j)*V(i,j)
      end do
    end do

    !Constructs the Householder reflector H = R - 2(vv^T)R/(v^Tv)
    R = R - 2.0*matmul(outerprod,R)

    !Reset norm values to calculate future norms correctly
    normA = 0.0
    normV = 0.0

    !Recursively multiply this quantity over and over, constructing Q with each iteration.
    Q = matmul(Q, Id - 2.0*outerprod)

  end do

  !for some reason Q and R are the opposite sign they should be, so as a tempory fix...
  Q = -Q
  R = -R

end subroutine HHQR

subroutine vander(A, vandMat, p)
  implicit none

  real, dimension(:,:), intent(IN) :: A
  real, allocatable, dimension(:,:), intent(OUT) :: vandMat
  integer, intent(INOUT) ::  p
  integer :: i,j,msize
  
  !This subroutine simply constructs the Vandermonde matrix. The input is a real matrix A, and a degree p. The output is the associated
  !Vandermonde matrix, and the associated interpolant degree p.

  msize = size(A,1)
  allocate(vandMat(msize, p+1)) !allocate vandermonde matrix to be msize tall, and the degree of the polynomial + 1 wide to interpolate data with.
  vandMat = 0.0

  !Construct Vandermonde matrix using matrix A.
  vandMat(:,1) = 1.0
  do j = 1,p
    do i = 1,msize
      vandMat(i,j+1) = A(i,1)**j    
    end do
  end do

end subroutine vander

end module LinAl
