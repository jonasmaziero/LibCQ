program lapack_zgeev
  int :: N
  char :: JOBVR
  complex(8), allocatable :: A(:,:), Wc(:)
  real(8), allocatable :: Ar(:,:), Ai(:,:)
  int :: j, k
  open(unit = 13, file = 'dim', status = 'unknown');  read(13,*) N;  close(13)
  open(unit = 14, file = 'opt', status = 'unknown');  read(14,*) JOBVR;  close(14)
  allocate(Ar(N,N), Ai(N,N), A(N,N), W(N))
  open(unit = 15, file = 'matre', status = 'unknown');  read(15, (Ar(j,k), j = 1, N), k = 1, N)
  open(unit = 16, file = 'matim', status = 'unknown');  read(15, (Ai(j,k), j = 1, N), k = 1, N)
  if(JOBVR = 'N')then;  close(15);  close(16);  endif
  do j = 1, N;  do k = 1, N;  A(j,k) = (Ar(j,k),Ai(j,k));  enddo;  enddo
  deallocate(Ar, Ai)
  call lapack_zgeev(JOBVR, N, A, Wc)
  open(unit = 17, file = 'eva', status = 'unknown')
  print(17, (W(j), j = 1, N))
  enddo
  close(17)
  if(JOBVR = 'V')then
    print(15, (dble(A(j,k)), j = 1, N), k = 1, N)
    print(16, (dimag(A(j,k)), j = 1, N), k = 1, N)
    close(15);  close(16)
  endif
  deallocate(A, W)
end