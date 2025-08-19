program EM_OU

 implicit none
      integer numpara,nsteps,ndata,nproc,ndrift,ndiff,num_iter,burn,mc,ndatat,nm
      real*8 beta_0,alpha_0,Infi,delta
  parameter (numpara=4,nproc=100,nsteps=100,ndatat=100,ndata=100,ndrift=2,ndiff=1,num_iter=50,mc=100,nm=500)
  integer ndata_all,i
  real*8 startx,x,integrate,lambda,a(2),b(2,2)
  real*8 pardrift(ndrift),del,pardiff(ndiff),redrift(nproc)
  parameter(Infi=999999999999999999999999999999.00,delta=0.5)
  REAL*8, allocatable, dimension(:) :: da 
  REAL*8, allocatable, dimension(:,:) :: datare,bridges_all
  real*8 tiempo_inicio, tiempo_final
  real*8  tiempo_total,alpha,beta,data(nproc,ndata),del_obs,gammas(nm),betas(nm)
integer,parameter :: seed = 805333
character*20 answer,file1,file2,file3,file4
allocate(datare(nproc,ndatat)) 
allocate(bridges_all(nproc,ndata_all))

call cpu_time(tiempo_inicio)




!print*,"delta_obs",del_obs
!print*,"delta_bri",del_bri
ndata_all=nsteps*(ndata-1)+ndata
do i=1,nm
  pardrift(1)=5.0
pardrift(2)=0.0
pardiff(1)=1.0
beta_0=3.0
alpha_0=1.0

del=delta/10.0
 11 call GenerateData(del,ndatat,ndrift,pardrift,ndiff,pardiff,nproc,datare)



 
call EM(Infi,seed,datare,ndata,nsteps,nproc,ndata_all,mc,ndrift,pardrift,ndiff,pardiff,del,num_iter,&
alpha_0,beta_0,gammas(i),betas(i))

if(ISNAN(gammas(i))) goto 11
print*,i,"gamma",gammas(i),sum(gammas(1:i))/i,"beta",betas(i),sum(betas(1:i))/i

enddo
file1="mg.txt"

file2="mb.txt"
open(1,file=file1)
open(2,file=file2)


write(1,*) gammas
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
!print*,data(1),data(2)
!call FormerDiffusionBridge(ndrift,para_diff,ndiff,pa_diff,delta,nsteps,data(1),data(2),bridge,numrej)
!print*,bridge
!call mat_bridge(data,ndata,nproc,mc,nsteps,ndrift,para_diff,ndiff,para_diff,del_bri,bridge_mc,bridges) 


!     call system("R CMD BATCH /Users/fernandobaltazarlarios/Documents/programing/diffusions/estimation/EM_t-diffu/gra.r gra.out")
     

end program EM_OU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine SEMV(nproc,ndata,data,delta,alpha,beta)
 implicit none
integer nproc,ndata,i,j
real*8 data(nproc,ndata),delta,alpha,beta,alphas(nproc),betas(nproc)


do i=1,nproc
 alphas(i)=-LOG(SUM(data(i,2:ndata)*data(i,1:(ndata-1)))/SUM(data(i,1:(ndata-1))**2))/delta
 betas(i)=(2.0*alphas(i))/(ndata*(1.0-EXP(-2.0*delta*alphas(i))))
 betas(i)=betas(i)*SUM((data(i,2:ndata)-data(i,1:(ndata-1))*EXP(-delta*alphas(i)))**2)
enddo

alpha=SUM(alphas)/nproc
beta=SUM(betas)/nproc


return

end subroutine


subroutine Sdata_reduction(nproc,ndatat,datat,ndata,data)
 implicit none
integer nproc,ndata,i,j,inc,ndatat,elem(ndata+1)
real*8 data(nproc,(ndata+1)),datat(nproc,ndatat)

inc=ndatat/ndata
elem(1)=1
do i=1,ndata
 elem(i+1)=i*inc
enddo
!print*,"ele",elem
data(:,:)=datat(:,elem)


return

end subroutine
       
  subroutine GenerateData(delta,ndata,ndrift,pardrift,ndiff,pardiff,nproc,datare)
  implicit none
  integer ndata,ndrift,ndiff,i,nproc
  real*8 delta,pardiff(ndiff),pardrift(ndrift),u1,u2,startx,lambda
  real*8 datare(nproc,ndata),data(ndata),rg(1),re(nproc)
    
 
    lambda=pardrift(1)
    
    
    do i=1,nproc
  call random_number(u1)

    startx=0.0
     call rgammas(1,1,lambda,re(i))
    

  !re(i)=rg(1)
      pardrift(1)=re(i)

     call diffusion(ndrift,pardrift,ndiff,pardiff,delta,startx,(ndata-2),data)
     
   !   call sim_ou_exac(startx,ndata,pardrift,ndrift,pardiff,ndiff,delta,data)
      datare(i,:)=data
      enddo
      
  !  print*,"real_alpha",re(1)  
      
  ! print*,"mean obs", SUM(re)/nproc
    
  
  return
  end subroutine 


  subroutine sim_ou_exac(startx,ndata,pardrift,ndrift,pardiff,ndiff,delta,data)
implicit none
integer i,ndata,ndrift,ndiff
real*8 pardrift(ndrift),pardiff(ndiff),delta,alpha,m,sd,data(ndata),var,beta,startx
  
  data(1)=startx
  alpha=pardrift(1)
  beta=pardiff(1)
  
  do i=2,ndata

    call normalvar(var)  
    sd=beta*sqrt((1-exp(-2*alpha*delta))/(2*alpha))
    m=exp(-alpha*delta)*data(i-1)
    data(i)=m+sd*var

  enddo
  
 
  return
end subroutine

subroutine STS(nproc,ndata,data,beta,delta,bb)
 implicit none
integer nproc,ndata,i,j
real*8 data(nproc,ndata),delta,bb(nproc),beta




bb(:)=(1.0/2.0)*(ndata-1.0)*delta-(1.0/(2.0*(beta**2)))*(data(:,ndata)**2-data(:,1)**2)!


return

end subroutine

subroutine SBS(nproc,ndata,ndata_all,matls,delbri,beta,Ys,BS)
 implicit none
integer nproc,ndata_all,i,j,ndata
real*8 matls(nproc,ndata_all),delbri,Ys(nproc,ndata_all),BS(nproc),beta


BS(:)=0.0


do i=1,nproc
  do j=2,ndata_all
    
    BS(i)=BS(i)+((Ys(i,j-1)+matls(i,j-1)/beta)**2)*delbri
  enddo
  

enddo

return

end subroutine

 subroutine SG1(nproc,ndata,data,delta,G1)
 implicit none
integer nproc,ndata,i,j
real*8 data(nproc,ndata),delta,G1

G1=0.0

do i=1,nproc
  do j=2,ndata
    
    G1=G1+(data(i,j)-data(i,j-1))**2/(2*delta)

    
  enddo
  

enddo

return

end subroutine

subroutine SG2(nproc,ndata,data,re,G2)
 implicit none
integer nproc,ndata,i
real*8 data(nproc,ndata),G2,re(nproc)

G2=0.0

do i=1,nproc
  
  G2=G2+re(i)*(data(i,ndata)**2-data(i,1)**2)/2.0

enddo

return

end subroutine


subroutine SES(nproc,ndata,ndata_all,matls,delbri,re,Ys,ES)
 implicit none
integer nproc,ndata_all,i,j,ndata
real*8 matls(nproc,ndata_all),delbri,ES(2),re(ndata),Ys(nproc,ndata_all)

ES(:)=0.0



do i=1,nproc
  do j=2,ndata_all
    
    ES(1)=ES(1)+(matls(i,j-1)**2)*delbri
    ES(2)=ES(2)+matls(i,j-1)*Ys(i,j-1)*delbri
!  print*,matls(i,j),Ys(i,j),matls(i,j)*Ys(i,j) 
  enddo
  ES(1)=ES(1)+ES(1)*(re(i)**2)/2.0
  ES(2)=ES(2)+ES(2)*(re(i)**2)
!  print*,ES
enddo
ES(2)=-ES(2)
!print*,ES
return

end subroutine



subroutine EM(Infi,seed,datare,ndata,nsteps,nproc,ndata_all,mc,ndrift,pardrift,ndiff,pardiff,delta,num_iter,&
  alpha_0,beta_0,alpha,beta)
implicit none
integer i,ndata,j,nsteps,ndrift,ndiff,numrej,nproc,k,mc,init,end,ndata_all,ns,burn,num_iter,seed
real*8 datare(nproc,ndata),pardrift(ndrift),pardiff(ndiff),del_bri,bridges_all(nproc,ndata_all),bridge(nsteps+2)
real*8 alphas(num_iter),betas(num_iter),delta,alpha_0,beta_0,E2,beta,bb(nproc),BS(nproc),min,as(nproc),abar,alpha
real*8 bridges_MC(mc*nproc,ndata_all),re(nproc),bri_lam(mc*nproc,ndata_all),IMC,Infi,fd(2),rem(nproc,mc),amc(mc)
real*8 G1,G2,E1,GS(2),seq(nsteps+2),ti,matls(nproc,ndata_all),delbri,data_lam(nproc,ndata),Ys(nproc,ndata_all)
character*20 answer,file1,file2,file3,file4
min=0.0
delbri=delta/((nsteps+1)*1.0)
ti=0.0

alphas(1)=5.0
betas(1)=1.0

    call SG1(nproc,ndata,datare,delta,G1)
    
    call SELES(nproc,ndata,datare,nsteps,ndata_all,delta,matls)


 do k=2,num_iter
 
    
 
    call MCMC(Infi,seed,mc,nproc,ndata,datare,nsteps,ndata_all,ndrift,alphas(1),ndiff,betas(1),delta,matls,abar,G2,E1,E2)
    
    
    GS(1)=G1
    GS(2)=G2
    alphas(k)=1.0/abar
    call Sbetas(nproc,ndata,E1,E2,GS,betas(k))
 !  print*,k,"beta",betas(k),"alpha",alphas(k)
  
  enddo

alpha=alphas(num_iter)
beta=betas(num_iter)


return 
   

end subroutine

!!! sequence 
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


 !!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
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
      subroutine diffusion(ndrift,pardrift,ndiff,pardiff,delta,startx,nsteps,points)
  implicit none
  integer nsteps,i,ndrift,ndiff
      real*8 delta,startx,y1,y2,pardrift(ndrift),pardiff(ndiff),y3
     REAL*8 points(nsteps+2)

!      print*,'startx',startx
!print*,'nsteps',nsteps

  points(1)=startx
  do i=2,nsteps+2
!print*,i
       call DriftParameter(ndrift,pardrift,points(i-1),y1) 

!print*,'i,drift',y1
         call DiffusionParameter(ndiff,pardiff,points(i-1),y2)

         call DerivateDiffusionParameter(ndiff,pardiff,points(i-1),y3)
!print*,'i,diff',y2


       call MilsteinStep(delta,points(i-1),y1,y2,y3,points(i))
!print*,'step2',points(i)


    enddo
!print*,points(1)



  return

  end
  ! Make a bridge from x to y
      subroutine FormerDiffusionBridge(ndrift,pardrift,ndiff,pardiff,aux,numsteps,x,y,bridge,numrej)
  implicit none
      integer i,numrej,j,mp,numsteps,ndrift,ndiff
  real*8 points(numsteps+2),x,y,pardrift(ndrift),pardiff(ndiff),delta
  real*8 points1(numsteps+2),points2(numsteps+2),aux
      real*8 points3(numsteps+2),bridge(numsteps+2)


!delta=aux/(numsteps+1)
delta=aux
  
  numrej=0
!print*,nnd,driftnd,ndiff,pardiff,delta,x,numsteps
1     call diffusion(ndrift,pardrift,ndiff,pardiff,delta,x,numsteps,points1)
!print*,points1

  
2     call diffusion(ndrift,pardrift,ndiff,pardiff,delta,y,numsteps,points2)
!print*,points2
!stop

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
     goto 20
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
     goto 20
    endif
       enddo
      endif
      
  numrej=numrej+1
 ! print*, numrej
      goto 1
 
20    do i=1,(mp-1)
       points(i)=points1(i)
      enddo

       do i=1,numsteps+2

          bridge(i)=points(i)
          enddo

  return
  end


!! complete path with bridges (data,ndata,nsteps,ndrift,pardrift,ndiff,pardiff,del_bri,bridges) 
subroutine mat_bridge(datare,ndata,nsteps,nproc,ndata_all,ndrift,redrift,ndiff,beta,delta,bridges_all) 
implicit none
integer i,ndata,j,nsteps,ndrift,ndiff,numrej,nproc,k,mc,init,end,ndata_all
real*8 datare(nproc,ndata),pardrift(ndrift),pardiff(ndiff),delta,bridges_all(nproc,ndata_all),bridge(nsteps+2)
real*8 beta,redrift(nproc),t0,tn,del_bri

del_bri=delta/(nsteps+1)

    pardiff(1)=beta
    do i=1,nproc
    pardrift(1)=redrift(i)

    do j=1,ndata-1
     init=nsteps*(j-1)+j
     end=init+nsteps+1

t0=0.0
tn=delta
    call bridge_ou_exact(t0,tn,datare(i,j),datare(i,j+1),nsteps,pardrift(1),beta,bridge)
    
 !    call FormerDiffusionBridge(ndrift,pardrift,ndiff,pardiff,del_bri,nsteps,datare(i,j),datare(i,j+1),bridge,numrej)
      bridges_all(i,init:end)=bridge
    ! print*,j,init,end 
   !  print*,"bridge",bridges_all(i,init:end)

    enddo
    enddo
   
 


return  
end subroutine

     subroutine DriftParameter(ndrift,pardrift,x,y)
      implicit none 
  integer ndrift
  real*8 pardrift(ndrift),x,y
  !y=-pardrift(1)*x
  
  y=pardrift(2)-pardrift(1)*x

!print*,y
  return
  end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DiffusionParameter(ndiff,pardiff,x,y)
  implicit none
  integer ndiff
  real*8 pardiff(ndiff),y,x

       ! y=pardiff(1)*sqrt(1+x**2)
       y=pardiff(1)
  return
  end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine DerivateDiffusionParameter(ndiff,pardiff,x,y)
  implicit none
  integer ndiff
  real*8 pardiff(ndiff),y,x

        y=0.00
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

  !!! sequence  SELES(nproc,ndata,datare,nsteps,ndata_all,delta,matls)
subroutine SELES(nproc,ndata,datare,nsteps,ndata_all,delta,matls)
implicit none
integer nproc,ndata,ndata_all,i,j,nsteps,ini,fin
real*8 time_ini,delta,matls(nproc,ndata_all),seq(nsteps+2),datare(nproc,ndata)

time_ini=0.0

call  sequence(time_ini,delta,nsteps,seq)
!print*,delta,seq
!stop
 do i=1,nproc

  do j=2,ndata
    ini=(nsteps+1)*(j-2)+1
    fin=(nsteps+1)*(j-1)+1
    matls(i,ini:fin)=((delta-seq)*datare(i,j-1)+seq*datare(i,j))/delta

   !stop
   enddo
 enddo
!print*,matls

!stop
return
end subroutine

  !!! sequence  SELES(nproc,ndata,datare,nsteps,ndata_all,delta,matls)
subroutine SELES_lam(nproc,ndata,datare,nsteps,ndata_all,beta,delta,matls)
implicit none
integer nproc,ndata,ndata_all,i,j,nsteps,ini,fin
real*8 time_ini,delta,matls(nproc,ndata_all),seq(nsteps+2),datare(nproc,ndata),beta,data_lam(nproc,ndata)

time_ini=0.0

call  sequence(time_ini,delta,nsteps,seq)
call SHlamp(nproc,ndata,datare,beta,data_lam)
 do i=1,nproc
  do j=2,ndata
    ini=(nsteps+1)*(j-2)+1
    fin=(nsteps+1)*(j-1)+1
    matls(i,ini:fin)=((delta-seq)*data_lam(i,j-1)+seq*data_lam(i,j))/delta
 !   print*,i,ini,fin
   !stop
   enddo
 enddo
!print*,matls


return
end subroutine

 !generate gammas integer parameter 

 subroutine rgammas(size,n,lambda,rg)
 implicit none
integer size,n,i,j
real*8 lambda,rg(size),u

rg(:)=0.0

do i=1,size
  do j=1,n
    call random_number(u)
    
    rg(i)=rg(i)-log(u)/lambda
  enddo
enddo

return

end subroutine


subroutine bridge_ou_exact(t0,tn,x0,xn,n,alpha,sigma,Z) 
implicit none
integer i,n,j
real*8 t0,tn,x0,xn,T(n+2),alpha,sigma,X(n+2),Z(n+2),norval,sd

  call sequence(t0,tn,n,T)
 
  X(1)=x0
  Z(1)=x0
 ! print*,"parameters",t0,tn,x0,xn,n,alpha,sigma
  
  i=2
  do while(i<=n+2)
    sd=sigma*sqrt((1-exp(-2*alpha*(T(i)-T(i-1))))/(2*alpha))
    call normalvar(norval)
    X(i)=exp(-alpha*(T(i)-T(i-1)))*X(i-1)+norval*sd
    i=i+1  
 !   print*,X(i)
  enddo
  
  j=2
 do while(j<=n+2)

    Z(j)=X(j)+(xn-X(n+2))*(exp(alpha*T(j))-exp(-alpha*T(j)))/(exp(alpha*T(n+2))-exp(-alpha*T(n+2)))
    j=j+1
  enddo
  
 return 
end subroutine


subroutine SHlamp(nproc,ndata,datare,beta,data_lam)
implicit none
integer nproc,ndata
real*8 datare(nproc,ndata),data_lam(nproc,ndata),beta

data_lam=datare/beta

 return 
end subroutine


 !!! sequence  SELES(nproc,ndata,datare,nsteps,ndata_all,delta,matls)
subroutine SYS(re,nproc,ndata,datare,nsteps,ndata_all,ndrift,ndiff,beta,delta,Ys)
implicit none
integer nproc,ndata,ndata_all,i,j,nsteps,ndrift,ndiff
real*8 pardrift(2),beta,delta,matls(nproc,ndata_all),datare(nproc,ndata),bridges_all(nproc,ndata_all)
real*8 YS(nproc,ndata_all),re(nproc),datalam(nproc,ndata),nbeta


datalam=datare/beta

call mat_bridge(datalam,ndata,nsteps,nproc,ndata_all,ndrift,re,ndiff,beta,delta,bridges_all) 
!call SHlamp(nproc,ndata,datare,beta,data_lam)
call  SELES(nproc,ndata,datalam,nsteps,ndata_all,delta,matls)


 YS=bridges_all-matls




return
end subroutine



!!! sequence  
subroutine MCMC(Infi,seed,mc,nproc,ndata,datare,nsteps,ndata_all,ndrift,alpha,ndiff,beta,delta,matls,ahat,G2,E1,E2)
implicit none
integer nproc,ndata,ndata_all,i,j,nsteps,ndrift,ndiff,mc,seed
real*8 pardrift(2),beta,delta,matls(nproc,ndata_all),datare(nproc,ndata),bridges_all(nproc,ndata_all)
real*8 YS(nproc,ndata_all),re(nproc),datalam(nproc,ndata),nbeta,GS2(mc),ESMC(2,mc),alpha
real*8 delbri,Infi,min,bb(nproc),BS(nproc),ES(2),as(mc,nproc),G2,ahat,E1,E2
delbri=delta/((nsteps+1)*1.0)


call rgammas(nproc,1,alpha,re)
!print*,nproc/SUM(re)

do i=1,mc

      call SG2(nproc,ndata,datare,re,GS2(i))
     ! print*,GS2(i)
      
      call SYS(re,nproc,ndata,datare,nsteps,ndata_all,ndrift,ndiff,beta,delta,Ys)
      call SES(nproc,ndata,ndata_all,matls,delbri,re,Ys,ES)
      ESMC(1,i)=ES(1)
      ESMC(2,i)=ES(2)
     ! print*, ESMC(:,i)
      
      call STS(nproc,ndata,datare,beta,delta,bb)
     ! print*,bb

      call SBS(nproc,ndata,ndata_all,matls,delbri,beta,Ys,BS)
      !print*,BS
      do j=1,nproc
!print*,bb(j),(bb(j)-alpha)/BS(j)
        call truncated_normal_ab_sample(((bb(j)-alpha)/BS(j)),SQRT(1.0/BS(j)),min,Infi,seed,re(j))
       
       as(i,j)=re(j)
   
  !     print*,bb(j),(bb(j)-alpha)/BS(j),re(j)
   !    print*,as(i,j)

      enddo
   
      !print*,ES(2)
    enddo
G2=SUM(GS2)/mc
ahat=SUM(as)/(mc*nproc)
E1=SUM(ESMC(1,:))/mc
E2=SUM(ESMC(2,:))/mc
return
end subroutine


subroutine Sbetas(nproc,ndata,E1,E2,GS,beta)
implicit none
integer nproc,ndata
real*8 E1,E2,GS(2),beta,k,num,dem,u,a,b

k=nproc*(ndata-1)




!IF (GS(1)+GS(2)+E1>0) THEN
!   beta=(4.0*(GS(1)+GS(2)+E1))/(-E2+SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k))
! print*,"1",beta
!ELSE IF (0>GS(1)+GS(2)+E1 .AND. 0>-(E2**2)/(8*k) .AND. E2>0) THEN
!beta=-(4.0*(GS(1)+GS(2)+E1))/(E2+SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k))
!  print*,"2",beta
! ELSE IF (GS(1)+GS(2)+E1==0 .AND. E2<0) THEN
  
!  beta=-E2/k
!   print*,"3",beta
!ELSE 
!  print*,"no hay"

!END IF

num=(8*(k**3)*(E2**2+8*(E1+GS(1)+GS(2)*k-E2*SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k))))
dem=(E2-SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k))**4
!print*,"cond",num/dem,GS(1)+GS(2)+E1,E2,k
beta=-E2/(2.0*k)+SQRT(E2**2+3.0*(GS(1)+GS(2)+E1)*k)/(2.0*k)
!print*,"b1",beta,E2/(2.0*k)+SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k)/(2.0*k),E2/(2.0*k)-SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k)/(2.0*k)
!ELSE 
!call random_number(u)
!beta=1.0+(b-a)*u+a
!END IF
!beta=(4.0*(GS(1)+GS(2)+E1))/(-E2+SQRT(E2**2+8.0*(GS(1)+GS(2)+E1)*k))
!print*,"b2",beta,GS(1)+GS(2)+E1

!print*,
 return 
end subroutine


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
  !print*,alpha,beta

  call normal_01_cdf ( alpha, alpha_cdf )
  call normal_01_cdf ( beta, beta_cdf )
!  print*,alpha_cdf,beta_cdf

call random_number(u)
!print*,"u",u
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
  call normal_01_cdf_inv ( xi_cdf, xi )

  x = mu + s * xi

!print*,mu,s,xi,x

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




