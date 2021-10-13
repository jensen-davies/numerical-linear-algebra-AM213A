Program Driver_LinAl

  use LinAl, only: mat, msize, nsize, readMat, traceMat, stdNorm, printMat

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real, allocatable, dimension(:) :: v
  real :: traceA, vnorm

  
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0

  allocate(v(msize))
  v = 0.0
  
  call readMat(myFileName)

  do i = 1, msize
     write(*,*) (mat(i,j) , j = 1, nsize )
  end do
  
  call traceMat(mat, msize, traceA)
  print*, 'Trace of matrix is:',traceA
  print*, 'Printing matrix...'
  call printMat(mat, msize, nsize)

  do i =1,nsize
    v = mat(:,i)
    call stdNorm(v,msize,vnorm)
    write(*,*) 'For column ',i,', the 2-norm is:'
    write(*,*), vnorm
  end do

  deallocate(mat)


End Program Driver_LinAl
