program Vortex
use fgsl
use, intrinsic :: iso_c_binding

implicit none
integer*8, allocatable :: m(:,:),A_G(:,:),A_X(:),m_new(:), k(:)
integer(c_size_t), pointer :: sample(:)
integer*8 :: i,n,nrep,nmax,seed,i1,j1,node,comp,ti,tlim,tini,nl
integer*8 :: ngh1,ngh2,nghx,neigh1,neigh2,neigh3,ngh0,nlines
real*8 :: magx, magy,m2,T,Jc,Jc1, E_new, E_old, HCar_new, HCar_old, HCross_new, HCross_old, Heis_x, Heis_y
real*8 :: HHeis_new, HHeis_old, Hw, Hd, m3
real*8 :: my, mvor, mx, mx2, mmed,var,varv,mvor2, mxav,myav,mmedav,mvor2av,mvorav,mmedav2
character :: dato*128,fmt1*8,salida*128,salida2*128,dato2*128,dato3*128,dato4*128, dato5*128
type(fgsl_rng) :: rng
type(fgsl_rng_type) :: ini_fgsl

!Read by terminal the arguments.

call get_command_argument(1, salida2) !Number of rhombus
read(salida2,*) nl

call get_command_argument(2, salida) !Temperature
read(salida,*) T

call get_command_argument(3, dato2) !Heisenberg interaction. The DM interaction is fixed at 1.0 by defect.
read(dato2,*) Jc

call get_command_argument(4, dato3) !Network adjacency matrix
read(dato3,*) 

call get_command_argument(5, dato4) !X and Y interactions to calculate the Jc1/2 strength
read(dato4,*) 

call get_command_argument(6, dato5) !Seed to perform the simulation
read(dato5,*) seed

!Start Gsl random numbers generator 
ini_fgsl = fgsl_rng_env_setup()
fgsl_rng_default_seed=seed
ini_fgsl = fgsl_rng_mt19937 !Mersenne Twister pseudo-random generator to ensure quasi-true randomness :))
rng = fgsl_rng_alloc(ini_fgsl)

N=4*nl*nl
!Allocate all the variables, considering the number of nodes
allocate(m(N,2),A_G(N,3),A_X(N),m_new(2),sample(N),k(N))

!Initialize variables
tlim=10000*N !Time limit (in MC steps)
tini=10*N !Time to start the measures
nmax=1!00 !Number of total averages
Jc1=1.0 !DM interaction strength
A_G=0
A_X=0


!Read the adjacency matrices
open (28, file=dato3, status='old') 
read(28,*) n,nlines

k=1
do i=1, nlines
read(28,*) i1,j1
A_G(i1,k(i1))=j1
A_G(j1,k(j1))=i1
k(i1)=k(i1)+1
k(j1)=k(j1)+1
enddo
close(28)

open (28, file=dato4, status='old')
read(28,*) n,nlines
do i=1, nlines
read(28,*) i1,j1
A_X(i1)=j1
A_X(j1)=i1
enddo
close(28)

!CHECK
!do i=1,N
!  write(*,*) i-1,A_G(i,:)-1,A_X(i)-1
!enddo

!Initialize variables to measure averages
mvor2av=0.0d0
mvorav=0.0d0
mxav=0.0d0
myav=0.0d0
mmedav=0.0d0
mmedav2=0.0d0

do nrep=1, nmax

  !Open the file to write the results on each average
  fmt1 = '(ES9.3)' ! an integer of width 5 with zeros at the left
  write (salida,fmt1) 1.0d0*N ! converting integer to string using a 'internal file'
  write (dato,fmt1)  ! converting integer to string using a 'internal file'
  write (dato2,fmt1) Jc ! converting integer to string using a 'internal file'
  write (dato3,fmt1) T ! converting integer to string using a 'internal file'
  salida2='OrderP_Size'//trim(salida)//'_str'//trim(dato2)//'_temp'//trim(dato3)//'_'
  open (unit = 25, file = salida2)

  write(*,*) 'Iteration',nrep

  !Initialize variables to measure the order of the system
  m=0.0d0
  mx=0.0d0
  my=0.0d0
  mmed=0.0d0
  mvor=0.0d0
  mvor2=0.0d0
  mx2=0.0d0
  !Initialize a state with random magnetization
  do i=1, N
    i1=fgsl_rng_uniform_int(rng, int(2,8))+1
    m(i,i1)=2d0*fgsl_rng_uniform_int(rng, int(2,8))-1
  enddo
  
  !Prepare the sample
  do i=1, N
    sample(i) = i
  enddo            

!Start the time (on MC steps)
do ti=1,tlim
  call fgsl_ran_shuffle(rng,sample,N) !Shuffle the vector to randomly update the nodes
  m3=0.0d0
  magx=0.0d0
  magy=0.0d0
  m2=0.0d0
  do i=1,N
    node=sample(i)
    
    !Flip the spin ensuring we will have a different arrow
    comp=fgsl_rng_uniform_int(rng, int(2,8))+1
    m_new=(0,0)
    if (m(node,comp)/=0) then
        m_new(comp)=-m(node,comp)
    else
    	m_new(comp)=2*fgsl_rng_uniform_int(rng, int(2,8))-1
    endif
    !Metropolis algorithm
    
    !Heisenberg interaction
    neigh1=A_G(node,1)
    neigh2=A_G(node,2)
    neigh3=A_G(node,3)

    Heis_x=m(neigh1,1)+m(neigh2,1)+m(neigh3,1)
    Heis_y=m(neigh1,2)+m(neigh2,2)+m(neigh3,2)
    HHeis_old=Jc*(-m(node,1)*Heis_x-m(node,2)*Heis_y)
    HHeis_new=Jc*(-m_new(1)*Heis_x-m_new(2)*Heis_y)
    
    !Dzyaloshinskii–Moriya interaction
    do i1=1,3
      if (A_G(node, i1) == (node + 3) .and. (A_G(node, i1) /= A_x(node))) ngh0=node+3 ! 
      if ((A_G(node,i1)==(node-1)).and.(A_G(node,i1)/=A_x(node))) ngh0=node-1
      if ((A_G(node,i1)==(node+1)).and.(A_G(node,i1)/=A_x(node))) ngh1=node+1
      if ((A_G(node,i1)==(node-3)).and.(A_G(node,i1)/=A_x(node))) ngh1=node-3
    enddo
    nghx=A_x(node)
    
    Hd=(m(node,1)*m(nghx,2)-m(node,2)*m(nghx,1))
    Hw=(m_new(1)*m(nghx,2)-m_new(2)*m(nghx,1))
    HCross_old=Jc1*(0.5*Hd+m(node,1)*m(ngh1,2)-m(node,2)*m(ngh1,1)-m(node,1)*m(ngh0,2)+m(node,2)*m(ngh0,1))
    HCross_new=Jc1*(0.5*Hw+m_new(1)*m(ngh1,2)-m_new(2)*m(ngh1,1)-m_new(1)*m(ngh0,2)+m_new(2)*m(ngh0,1))
    
    !Charge exclusion term on the rhombus
    do i1=1,3
      if ((A_G(ngh1,i1)/=node).and.(A_G(ngh1,i1)/=A_x(ngh1))) ngh2=A_G(ngh1,i1)
    enddo
    
     if (mod(node,2) == 0) then
       HCar_old=abs(m(node,2)-m(ngh1,1)+m(ngh0,1)-m(ngh2,2))+abs(m(node,1)-m(ngh1,2)+m(ngh0,2)-m(ngh2,1))
       HCar_new=abs(m_new(2)-m(ngh1,1)+m(ngh0,1)-m(ngh2,2))+abs(m_new(1)-m(ngh1,2)+m(ngh0,2)-m(ngh2,1))
     else
       HCar_old=abs(m(node,2)+m(ngh1,1)-m(ngh0,1)-m(ngh2,2))+abs(m(node,1)+m(ngh1,2)-m(ngh0,2)-m(ngh2,1))
       HCar_new=abs(m_new(2)+m(ngh1,1)-m(ngh0,1)-m(ngh2,2))+abs(m_new(1)+m(ngh1,2)-m(ngh0,2)-m(ngh2,1))
     endif       
       
     !Computation of the total energy
     E_old=-HCross_old+HHeis_old+HCar_old
     E_new=-HCross_new+HHeis_new+HCar_new
     m3=HCross_old
     if ((E_new<E_old).or.(fgsl_rng_uniform(rng)<(exp(-(E_new-E_old)/T)))) then
       m(node,:)=m_new(:)
       m3=m3-HCross_old
       m3=m3+HCross_new  
     endif 
     
     !Perform measures only when the system thermalize                 
     if (ti>=tini) then
       magx=magx+m(node,1)
       magy=magy+m(node,2)
       m2=m2+m3
     endif
   enddo

!Normalize the measures and accumulate. We measure total magnetization, the squared order parameter and the vortex order parameter.
magx=magx/(1.0d0*N)
magy=magy/(1.0d0*N)
m2=m2/(2.0d0*N)
mx=mx+magx
my=my+magy
mmed=mmed+sqrt(magx*magx+magy*magy)
mx2=mx2+sqrt(magx*magx+magy*magy)*sqrt(magx*magx+magy*magy)
mvor2=mvor2+m2*m2
mvor=mvor+m2

write(45,*) magx,magy,sqrt(magx*magx+magy*magy), m2, ti
enddo

!Average over runs
!rsq=rsq+dsqrt(rmed-rkmed*rkmed)
!ravrg=ravrg+dsqrt(rmed-kurR*kurR)
!sq=sq+(rmed-kurR*kurR)
mxav=mxav+mx/(1.0d0*tlim-1.0d0*tini)
myav=myav+mx/(1.0d0*tlim-1.0d0*tini)
mvorav=mvorav+mvor/(1.0d0*tlim-1.0d0*tini)
mmedav=mmedav+mmed/(1.0d0*tlim-1.0d0*tini)
mmedav2=mmedav2+mx2/(1.0d0*tlim-1.0d0*tini)
mvor2av=mvor2av+mvor2/(1.0d0*tlim-1.0d0*tini)

var=(mmedav2/(1.0d0*nrep))-((mmedav/(1.0d0*nrep))*(mmedav/(1.0d0*nrep)))
varv=(mvor2av/(1.0d0*nrep))-((mvorav/(1.0d0*nrep))*(mvorav/(1.0d0*nrep)))
write(25,*) mxav/(1.0d0*nrep),myav/(1.0d0*nrep),mmedav/(1.0d0*nrep),var,mvorav/(1.0d0*nrep),varv,nrep,T
!Write on file for each average all the measures
!write(45,*) ravrg/(1.0d0*nrep),sq/(1.0d0*nrep),rk/(1.0d0*nrep),rsq/(1.0d0*nrep),sigma,w0,s,nrep

close(25)
enddo !Sobre nreps

stop
end
