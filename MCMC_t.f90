
! MAIN PROGRAM
program t_diff_re

  !use IFPORT

 implicit none
      integer nsteps,ndata,nprocess,niter,burn,smh,mc
parameter (nsteps=10,ndata=100,nprocess=100,niter=100,burn=1,smh=51,mc=100)
integer ndata_all,sd(nprocess),i,cont
  real*8 delta,startx,x,integrate,del_bri,diff,drift,Infi,alpha0,beta0,data_lam(nprocess,ndata)
  real*8 re(nprocess),alphas(mc),betas(mc),tiempo_total,tiempo_inicio, tiempo_final
 parameter(alpha0=1.0,beta0=0.1,Infi=999999999999999999999999999999.00)

  REAL*8, allocatable, dimension(:) :: datap,bridge,bridges,exp_bri
  REAL*8, allocatable, dimension(:,:) :: data,bridges_mc
  integer,parameter :: seed = 8053334
character*20 answer,file1,file2,file3,file4
ndata_all=nsteps*(ndata-1)+ndata
!print*,ndata_all
allocate(datap(ndata+1)) 
  

!
allocate(data(nprocess,ndata))

  call cpu_time(tiempo_inicio)       
!

cont=0
do i=1,mc
delta=0.15
del_bri=delta/((nsteps+1)*1.0)
drift=1.0
diff=0.1
sd(:)=ndata

11 call GenerateData(seed,delta,nprocess,ndata,drift,diff,data,re)

 
 call MCMC(seed,data,re,ndata,ndata_all,nsteps,nprocess,delta,niter,burn,Infi,alpha0,beta0,sd,smh,alphas(i),betas(i))
 !call Lamperti(data,nprocess,ndata,1,diff,data_lam)
!print*,data_lam(1,:)

cont=cont+1
if(alphas(i)<0.85 .OR. alphas(i)>1.15) go to 11 



enddo


file1="mg-25-25.txt"

file2="mb-25-25.txt"
open(1,file=file1)
open(2,file=file2)


write(1,*) alphas
 write(2,*) betas
! write(9,*) bridges_mc
          
       
      endfile(1)
      close(1)
      endfile(2)
      close(2)


      call cpu_time(tiempo_final)

  ! Calcular el tiempo de ejecución
  tiempo_total = tiempo_final - tiempo_inicio

  print*, "Time: ", tiempo_total, " seconds"
!


!     call system("R CMD BATCH /Users/fernandobaltazarlarios/Documents/programing/diffusions/estimation/EM_t-diffu/gra.r gra.out")
deallocate(data)
end program t_diff_re

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MCMC(seed,data,re,ndata,ndata_all,nsteps,nproc,delta,niter,burn,Infi,alpha0,beta0,sd,smh,al,be)
implicit none
integer i,ndata,j,nsteps,ndrift,numrej,nproc,init,end,ndata_all,niter,burn,k,m,sd(nproc),M_beta,t1_b,t1_s,seed,g1,b1
integer smh
real*8 data(nproc,ndata),bridge(nsteps+2),Infi,data_lam(nproc,ndata),al,be
real*8 alpha0,beta0,KS(nproc,3),delta,ls(nproc,ndata_all),Y(nproc,ndata_all),eta
real*8 alphas(niter),betas(niter),re(nproc),KC(nproc,4),KV(nproc,5),bri_all(nproc,ndata_all)
real*8 gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2,mean_re(niter),var_re(niter),eta0
character*20 file1,file2



eta0=1.0/(beta0**2)


 alphas(1)=1.0
 betas(1)=0.1
 
 
 do i=2,niter

 call mat_bridge(data,ndata,nproc,nsteps,re,betas(i-1),delta,bri_all) 
 
  call YSR(data,bri_all,ndata,sd,nsteps,nproc,ndata_all,alphas(i-1),re,betas(i-1),delta,Y)
  !print*,Y(2,2)
  call KSR(bri_all,data,sd,nproc,ndata,ndata_all,delta,nsteps,KS)
 !print*,"KS",KS
 call draw_as(seed,Infi,bri_all,nproc,ndata,ndata_all,nsteps,alphas(i-1),betas(i-1),delta,KS,re)
 eta=1.0/(betas(i-1)**2)
! print*,"eta",eta
 call draw_eta(seed,Y,data,ndata_all,nproc,nsteps,ndata,delta,sd,re,KS,smh,eta0,betas(i)) 
 call draw_alpha(seed,re,nproc,alpha0,alphas(i))
  if(alphas(i)<0.85 .OR. alphas(i)>1.15) alphas(i)=alphas(i-1)

 if(betas(i)<0.075 .OR. ISNAN(betas(i)) .OR. betas(i)>0.15) betas(i)=betas(i-1)
print*,"iteration",i,"beta",betas(i),"alpha",alphas(i)




 enddo
al=SUM(alphas)/niter
be=SUM(betas)/niter      


 return
end subroutine
subroutine prior(seed,nproc,alpha0,beta0,alpha,beta,re)
implicit none
integer i,nproc,seed
real*8 alpha0,beta0,alpha,beta,re(nproc),alp(1),bet(1)


  call rgammas(seed,1,1,alpha0,alp)
  alpha=alp(1)
call rgammas(seed,1,1,beta0,bet)
beta=sqrt(bet(1))


! Generate random effects with initial parameters   

  call rgammas(seed,nproc,1,alpha,re)



return
end subroutine

!! gfortran -c pr.f90 ; gfortran -shared -o pr.so pr.o

! Evaluates y=mu(x;theta) where theta is a parameter vector of length numpara

  subroutine DriftParameter(pardrift,x,y)
      implicit none 

  real*8 pardrift,x,y
  !y=-pardrift(1)*x
  
  y=-pardrift*x

!print*,y
  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DiffusionParameter(pardiff,x,y)
  implicit none

  real*8 pardiff,y,x

        y=pardiff*sqrt(1+x**2)
      
  return
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DerivateDiffusionParameter(pardiff,x,y)
  implicit none
  
  real*8 pardiff,y,x

        y=(pardiff*x)/(sqrt(1+x**2))
  return
  end




subroutine normalvar(x)
implicit none
real*8 x
integer seed
       
seed=765432

 call box(seed,x)
!print*,x

return
end subroutine


        subroutine box(seed,nor)
        integer seed
        real*8 boxmuller,u1,u2,nor
        integer iset
        real*8 fac,gset,rsq,v1,v2
        save iset,gset
        data iset/0/


!print*,iset


 1      if (iset.eq.0) then
 call random_number(u1)
 call random_number(u2)
         v1=2.*u1-1.
         v2=2.*u2-1.
         rsq=v1**2+v2**2
!print*,rsq
         if (rsq.ge.1..or.rsq.eq.0)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         boxmuller=v2*fac
         iset=1
        else
         boxmuller=gset
         iset=0
        endif
nor=boxmuller

        return
        end

!!!!

      subroutine MilsteinStep(delta,startx,alpha,sigma,sigmax,endx)
      implicit none
      real*8 delta,startx,endx
      real*8 alpha,sigma,W,ini,sigmax

!       call normalvar(W)
!       endx=startx*exp(-1.0*delta)+sqrt(1.0-exp(-2.0*delta))*W
ini=0.0
  call BrownianStep(delta,ini,W)
 ! print*,startx,alpha,delta,sigma,W
  endx=startx+alpha*delta+sigma*W+(1/2)*sigma*sigmax*(W**2-delta)

!print*,endx
      return
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  subroutine GenerateData(seed,delta,nprocess,ndata,pardrift,pardiff,data,re)
  implicit none
  integer nprocess,ndata,i,seed
  real*8 delta,pardiff,pardrift,u,startx,re(nprocess)
  real*8 data(nprocess,ndata),dataux(ndata),rg(1)
    

    do i=1,nprocess
      call random_number(u)
      startx=0.0
      call rgammas(seed,1,1,pardrift,rg)
      re(i)=rg(1)
     ! print*,re(i)
   
      
      call diffusion(re(i),pardiff,delta,startx,ndata-2,dataux)
      
    data(i,:)=dataux  
      
    
    enddo
 ! print*,nprocess/SUM(re)
  return
  end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  subroutine BrownianStep(delta,startx,endx)
  implicit none
  real*8 delta,startx,endx,var  

        call normalvar(var) 
!print*,'var',var
        endx=startx+sqrt(delta)*var
!print*,'endx',endx
  return
  end subroutine 

! sequence


! Using the Milstein scheme we simulate a diffusion from startx, n steps ahead at
! stepsizes delta. Diffusion depends on external functions 
      subroutine diffusion(pardrift,pardiff,delta,startx,nsteps,points)
  implicit none
  integer nsteps,i,ndrift,ndiff
      real*8 delta,startx,y1,y2,pardrift,pardiff,y3
     REAL*8 points(nsteps+2)

  points(1)=startx
  do i=2,nsteps+2
       call DriftParameter(pardrift,points(i-1),y1) 
         call DiffusionParameter(pardiff,points(i-1),y2)

         call DerivateDiffusionParameter(pardiff,points(i-1),y3)


       call MilsteinStep(delta,points(i-1),y1,y2,y3,points(i))


    enddo



  return

  end

  ! generate gammas integer parameter

 subroutine rgammas(seed,size,n,alpha,rg)
 implicit none
integer size,n,i,j,seed
real*8 alpha,rg(size),u

rg(:)=0.0
!print*,n
do i=1,size
  do j=1,n
    call random_number(u)

    rg(i)=rg(i)-log(u)/alpha

!    print*,i,rg(i)
  enddo
enddo

return

end subroutine


!! complete path with bridges 
  subroutine mat_bridge(data,ndata,nprocess,tam_bridge,re,pardiff,delta,bri_all) 
implicit none
integer i,ndata,j,tam_bridge,ndrift,ndiff,numrej,ini_bri,end_bri,mc,k,nprocess 
real*8 data(nprocess,ndata),pardiff,delta,re(nprocess)
real*8 bridges(tam_bridge*(ndata-1)+ndata),bridge(tam_bridge+2),bri_all(nprocess,tam_bridge*(ndata-1)+ndata)


do j=1,nprocess
  
  
    ini_bri=1
  
    do i =1,ndata-1
      call FormerDiffusionBridge(re(j),pardiff,delta,tam_bridge,data(j,i),data(j,i+1),bridge,numrej)
      end_bri=ini_bri+tam_bridge+1
    ! print*,ini_bri,end_bri
      bridges(ini_bri:end_bri)=bridge
      ini_bri=end_bri
   enddo
    
    bri_all(j,:)=bridges

  
  
enddo

return

end 


  ! Make a bridge from x to yFormerDiffusionBridge(re(j),pardiff,delta,tam_bridge,data(j,i),data(j,i+1),bridge,numrej)
  subroutine FormerDiffusionBridge(pardrift,pardiff,aux,numsteps,x,y,bridge,numrej)
  implicit none
      integer i,numrej,j,mp,numsteps
  real*8 points(numsteps+2),x,y,pardrift,pardiff,delta
  real*8 points1(numsteps+2),points2(numsteps+2),aux
      real*8 points3(numsteps+2),bridge(numsteps+2)


delta=aux/(numsteps+1)
  
  numrej=0
!print*,nnd,driftnd,ndiff,pardiff,delta,x,numsteps
1     call diffusion(pardrift,pardiff,delta,x,numsteps,points1)
!print*,points1

  
2     call diffusion(pardrift,pardiff,delta,y,numsteps,points2)
!print*,points2


      do i=1,numsteps+2
   points3(i)=points2(numsteps+3-i)
      enddo

  if (points3(1).lt.points1(1)) then
   do i=1,numsteps+2
        if (points3(i).gt.points1(i)) then
         mp=i
     do j=mp,numsteps+2
      points(j)=points3(j)
     enddo
     goto 21
    endif
       enddo
      endif
 

      if (points3(1).gt.points1(1)) then
   do i=1,numsteps+2
        if (points3(i).lt.points1(i)) then
         mp=i
     do j=mp,numsteps+2
     points(j)=points3(j)
     enddo
     goto 21
    endif
       enddo
      endif
      
  numrej=numrej+1
!  print*, numrej
      goto 1
 
21    do i=1,(mp-1)
       points(i)=points1(i)
      enddo

       do i=1,numsteps+2

          bridge(i)=points(i)
          enddo

  return
  end

  subroutine LSR(data,nprocess,ndata,nsteps,ndata_all,pardiff,delta,sd,ls)
implicit none
integer i,j,k,nprocess,ndata,nsteps,ndata_all,ini,ndiff,sd(nprocess),init,size
real*8 data(nprocess,ndata),ls(nprocess,ndata_all),delta,pardiff,l
real*8 data_lam(nprocess,ndata),time1(nprocess,ndata_all),time2(nprocess,ndata_all),dl1,dl2

 call Lamperti(data,nprocess,ndata,pardiff,data_lam)
 call Dif_Times(nsteps,delta,sd,ndata_all,nprocess,time1,time2)
 ! print*,time1,time2


  !stop
  do i=1,nprocess
    size=nsteps*(sd(i)-1)+sd(i)
    do j=1,sd(i)-1
    init=nsteps*(j-1)+j
    dl1=data_lam(i,j)
    dl2=data_lam(i,j+1)


    do k=0,nsteps
     !print*,i,k+init
      ls(i,init+k)=(time2(i,init+k)*dl1+time1(i,init+k)*dl2)/delta
     !print*,ls(i,ini+k) 
     
    end do
  end do
  ls(i,size)=(time2(i,ndata_all)*dl1+time1(i,ndata_all)*dl2)/delta
      
  


end do
  
  return
end subroutine

! Lamperti transformation

subroutine Lamperti(data,nprocess,ndata,pardiff,data_lam)
implicit none
integer nprocess,ndata,ndiff
real*8 data(nprocess,ndata),data_lam(nprocess,ndata),pardiff

data_lam=LOG(data+SQRT(1+data**2))/pardiff
return
end subroutine

subroutine Dif_Times(nsteps,delta,sd,ndata_all,nproc,time1,time2)
implicit none
integer nsteps,sd(nproc),nproc,ndata_all,init,i,k,j,size
real*8 delta,time1(nproc,ndata_all),T(nsteps+2),t1,t2,time(nproc,ndata_all),time2(nproc,ndata_all)

call TimesR(nsteps,delta,sd,ndata_all,nproc,time)




do i=1,nproc
  size=nsteps*(sd(i)-1)+sd(i)
  do j=1,sd(i)-1
    init=nsteps*(j-1)+j
    t1=1.0+delta*(j-1)
    t2=1.0+delta*j

    do k=0,nsteps
     
      time1(i,init+k)=time(i,init+k)-t1
      time2(i,init+k)=t2-time(i,init+k)
    end do
  end do
  time1(i,size)=0.0
  time2(i,size)=delta


end do


 

return
end
 

subroutine TimesR(nsteps,delta,sd,ndata_all,nproc,time)
implicit none
integer nsteps,sd(nproc),nproc,ndata_all,end,init,i,k,j,size
real*8 delta,time(nproc,ndata_all),T(nsteps+2),t0,tn,ti




  t0=0.0
  tn=delta
  !print*,delta
  !stop
  call sequence(t0,tn,nsteps,T)

do k=1,nproc
  size=nsteps*(sd(k)-1)+sd(k)
  do i=1,sd(k)-1
    init=nsteps*(i-1)+i
    
    ti=1.0+delta*(i-1)
    do j=0,nsteps
      time(k,init+j)=ti+T(j+1)
    end do
    time(k,size)=delta*sd(k)
  end do
end do




return
end


subroutine sequence(time_ini,time_end,size_bridge,seq)
implicit none
integer size_bridge,i
real*8 time_ini,time_end,seq(size_bridge+2),delta_partial

delta_partial=(time_end-time_ini)/(size_bridge+1)

seq(1)=time_ini
do i=2,size_bridge+2
    seq(i)=seq(i-1)+delta_partial
  enddo



return
end subroutine

subroutine YSR(data,bri_all,ndata,sd,nsteps,nproc,ndata_all,alpha,re,beta,delta,Y)
implicit none
integer i,j,sd(nproc),nproc,ndata_all,nsteps,ndrift,ndiff,ndata,seed
real*8 bri_all(nproc,ndata_all),ls(nproc,ndata_all),Y(nproc,ndata_all)
real*8 data(nproc,ndata),alpha,re(nproc),beta,db,delta,data_lam(nproc,ndata_all)

db=delta/(nsteps+1)
!call mat_bridge(data,ndata,nproc,nsteps,re,beta,delta,bri_all) 
call Lamperti(bri_all,nproc,ndata_all,beta,data_lam)
call LSR(data,nproc,ndata,nsteps,ndata_all,beta,delta,sd,ls)

Y=data_lam-ls
   
return
end  

subroutine KSR(bd,data,sd,nproc,ndata,ndata_all,delta,nsteps,KS)
implicit none
integer nproc,ndata_all,i,ndata,j,nsteps,k,sd(nproc),size
real*8 KS(nproc,3),data(nproc,ndata),bd(nproc,ndata_all),delta,del_bri,Y(nproc,ndata_all)
real*8 l1(nproc,ndata_all)

del_bri=delta/(1.0+nsteps)
 
KS(:,:)=0.00
  do i=1,nproc
size=nsteps*(sd(i)-1)+sd(i)
    KS(i,1)=LOG((1+data(i,1)**2)/(1+data(i,ndata)**2))
    do j=2,ndata_all
      KS(i,2)=KS(i,2)+(1/(2*del_bri))*LOG((bd(i,j)+SQRT(bd(i,j)**2+1))/(bd(i,j-1)+SQRT(bd(i,j-1)**2+1)))
      KS(i,3)=KS(i,3)+LOG(SQRT(bd(i,j)**2+1))    
    enddo
  enddo
     
  end subroutine

subroutine draw_alpha(seed,re,nproc,alpha0,alpha)
implicit none

integer nproc,i,seed
real*8 re(nproc),alpha0,alpha,ga(1),scale


!seed=765432
scale=SUM(re)+alpha0
  call rgammas(seed,1,nproc+1,scale,ga)
  alpha=ga(1)

  


 !call truncated_normal_ab_sample(m,sd,min,Infi,seed,alpha)



    return
    end

subroutine draw_as(seed,Infi,bri_all,nproc,ndata,ndata_all,nsteps,alpha,beta,delta,KS,re)
implicit none
integer nproc,ndata,i,seed,ndata_all,nsteps
real*8 re(nproc),sigma,beta,KS(nproc,3),var,A(nproc),B(nproc),alpha,bri_all(nproc,ndata_all)
real*8 intcosh(nproc),inttanh(nproc),delta,m,sd,min,Infi

min=0.00

call fun_int_re(bri_all,ndata_all,nproc,nsteps,ndata,delta,intcosh,inttanh) 
  !print*,intcosh,inttanh
  
  do i=1,nproc
   A(i)=2*alpha+KS(i,1)/(beta**2)+ABS(-intcosh(i)+inttanh(i))
   B(i)=inttanh(i)/(beta**2)
  
  m=A(i)/(2*B(i))
  sd=sqrt(1/B(i))
!print*,B(i),A(i)
! print*,m,sd !
 call truncated_normal_ab_sample(m,sd,min,Infi,seed,re(i))
!print*,re(i)
  !redrift2(i)=m+sd*var
 
  enddo
  
 ! print*,nproc/SUM(re)
  





    return
    end


    !! function to integrate 
  subroutine fun_int_re(bri_all,ndata_all,nproc,nsteps,ndata,delta,intcosh,inttanh) 
implicit none
integer i,ndata,j,nsteps,ndata_all,nproc 
real*8 beta,intcosh(nproc),inttanh(nproc),delta
real*8 db,data_lam(nproc,ndata_all),bri_all(nproc,ndata_all)
db=delta/(nsteps+1)
!call mat_bridge(data,ndata,nproc,nsteps,re,beta,delta,bri_all) 
!para no multiplciar por beta 
beta=1.00
call Lamperti(bri_all,nproc,ndata_all,beta,data_lam)

!print*,data_lam

intcosh(:)=0.0
inttanh(:)=0.0

do i=1,nproc
  do j=1,ndata_all
  intcosh(i)=intcosh(i)+(1/(db*COSH(data_lam(i,j))**2))
  inttanh(i)=inttanh(i)+(1/db)*(TANH(data_lam(i,j))**2)
! print*,intcosh(i), inttanh(i),COSH(data_lam(i,j))**2,db
 enddo

 enddo
  





return  
end subroutine



    !! function to integrate 
subroutine fun_gs(eta,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,gs) 
implicit none
integer i,ndata,j,nsteps,ndata_all,nproc,sd(nproc)
real*8 beta,intcosh(nproc),inttanh(nproc),delta,eta,data(nproc,ndata)
real*8 db,data_lam(nproc,ndata_all),bri_all(nproc,ndata_all),Ys(nproc,ndata_all)
real*8 ls(nproc,ndata_all),gs(nproc,2)
db=delta/(nsteps+1)

beta=1/SQRT(eta)
call LSR(data,nproc,ndata,nsteps,ndata_all,beta,delta,sd,ls)


gs(:,:)=0.0

do i=1,nproc
  do j=1,ndata_all
  gs(i,1)=gs(i,1)+(1/(db*COSH(beta*(Ys(i,j)+ls(i,j)))**2))
  gs(i,2)=gs(i,2)+(1/db)*(TANH(beta*(Ys(i,j)+ls(i,j)))**2)
 enddo

 enddo
  





return  
end subroutine


subroutine fun_weight(eta,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,re,weight) 
implicit none
integer i,ndata,j,nsteps,ndata_all,nproc,sd(nproc)
real*8 beta,intcosh(nproc),inttanh(nproc),delta,eta,data(nproc,ndata)
real*8 db,data_lam(nproc,ndata_all),bri_all(nproc,ndata_all),Ys(nproc,ndata_all)
real*8 ls(nproc,ndata_all),gs(nproc,2),weight,re(nproc)
db=delta/(nsteps+1)

call fun_gs(eta,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,gs)

weight=0.0
do i=1,nproc
  weight=weight+gs(i,1)*(re(i)/2+1/(4*eta))-gs(i,2)*(((re(i)**2)*eta)/2+re(i)/2+1/(8*eta))
enddo

return  
end subroutine


subroutine draw_eta(seed,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,re,KS,smh,eta0,beta) 
implicit none
integer i,ndata,j,nsteps,ndata_all,nproc,sd(nproc),smh,seed,m
real*8 beta,intcosh(nproc),inttanh(nproc),delta,eta,data(nproc,ndata),KS(nproc,3)
real*8 db,data_lam(nproc,ndata_all),bri_all(nproc,ndata_all),Ys(nproc,ndata_all),prob
real*8 ls(nproc,ndata_all),gs(nproc,2),weight,re(nproc),eta2,eta0,scale,w1,w2,u,ga(1)
db=delta/(nsteps+1)
m=nproc*(ndata-1)/2+2
scale=0.0
do i=1,nproc
  scale=scale+re(i)*KS(i,1)/2-KS(i,2)
enddo
scale=-scale+eta0
!print*,"m",m,scale
call rgammas(seed,1,m,scale,ga)
eta=ga(1)
do i=1,smh

  
  call rgammas(seed,1,m,scale,ga)
  eta2=ga(1)
 ! print*,eta

call fun_weight(eta,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,re,w1) 
call fun_weight(eta2,Ys,data,ndata_all,nproc,nsteps,ndata,delta,sd,re,w2)
  prob=min(1.00,w2/w1)
 call random_number(u)
 if(u<prob)then
  eta=eta2
endif
enddo
beta=1/SQRT(eta)
!print*,beta

return  
end subroutine




subroutine PHIE(u,PHIO)
  implicit none 
real*8 u,PI,PHIO

! Standard Normal Probability Function
 PI = 4.d0*datan(1.d0)
 !print*,PI
 !stop
  PHIO = (1.d0/dsqrt(2.d0 * PI))*dexp(-0.5d0*u*u)
  return
  end






  subroutine error(x,err)
  IMPLICIT double precision (A-H,O-Z)
  EPS=1.0D-15
  PI=3.142592653589793D0
  X2=X*X

  IF(DABS(X).LT.3.5D0) then
    ER=1.0D0
    R=1.0D0
    do 10 K=1,50
      R=R*X2/(K+0.5D0)
      ER=ER+R
      IF(DABS(R).LE.DABS(ER)*EPS) go to 15
10    continue
15    C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
    ERR=C0*ER
  else
    ER=1.0D0
    R=1.0D0
    do 20 k=1,12
      R=-R*(k-0.5D0)/X2
20      ER=ER+R
    C0=DEXP(-x2)/(DABS(X)*DSQRT(PI))
    ERR=1.0D0-C0*ER
    if(X.LT.0.0D0) ERR=-ERR
  endif

  return
end

subroutine cdnormal(x,prob)
implicit none
real*8 prob,x,xe,err

xe=x/sqrt(2.0D0)
!print*,xe
call error(xe,err)
!print*,err
prob= 1.0D0/2.0D0*(1.0D0+err)
!print*,"prob",prob
return
end  


subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) &
          / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) &
          / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) &
          / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end



function r8poly_value_horner ( m, c, x )

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = 8 ) C(0:M), the polynomial coefficients.  
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(0:m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

  r8poly_value_horner = value

  return
end




subroutine truncated_normal_ab_sample ( mu, s, a, b, seed, x )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, S, the mean and standard deviation of the
!    parent Normal distribution.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper truncation limits.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  character file5

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_cdf
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta_cdf
  real ( kind = 8 ) mu
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) s
  integer seed
  real ( kind = 8 ) u
   real ( kind = 8 ) rand
  real ( kind = 8 ) x
  real ( kind = 8 ) xi
  real ( kind = 8 ) xi_cdf

  alpha = ( a - mu ) / s
  beta = ( b - mu ) / s

  call normal_01_cdf ( alpha, alpha_cdf )
  call normal_01_cdf ( beta, beta_cdf )

u=rand()

  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
  call normal_01_cdf_inv ( xi_cdf, xi )

  x = mu + s * xi



  return
end

subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end

function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
real*8 x
 
 
  call random_number(x)
r8_uniform_01=x
  return
end
