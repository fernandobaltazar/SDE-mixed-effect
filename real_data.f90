
! MAIN PROGRAM
program OU_mean
!Use   mod_tru


 implicit none
      integer ns,ndata,nproc,ndrift,ndiff,num_iter,burn,M_beta,M_all
      real*8 gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2
      real*8 alpha,beta,mu,sigma,Infi,delta
      real*8 u(1000)
  parameter (nproc=311,ns=10,ndata=33937,ndrift=2,ndiff=1,num_iter=1000,M_all=1813910)

  parameter(alpha=0.1,beta=2.0,mu=2.5,sigma=5.0,burn=1,delta=0.1,M_beta=1000,Infi=999999999999999999999999999999.00)
  parameter(alpha_2=-500000,beta_1=1,beta_2=90.0,mu_0=mu,sigma_0=0.001,gamma_1=10.0,gamma_2=10-0)
  integer ndata_all,sd(nproc),i
  real*8 pardrift(ndrift),pardiff(ndiff),redrift(nproc),x,data_all(M_all)
  real*8, allocatable, dimension(:,:) :: datare,bridges_all,time,time1,time2,l
  integer,parameter :: seed = 80533534

allocate(datare(nproc,ndata))
 
allocate(bridges_all(nproc,ndata_all))


OPEN(UNIT=10,FILE='size.txt',STATUS='old') 
OPEN(UNIT=11,FILE='LIF_datavec.txt',STATUS='old') 
  READ(10,*) sd
  READ(11,*) data_all

!    print*,sd
    !print*,data_all(1:3)
   ! CLOSE(1)

call srand(seed)




ndata_all=ns*(ndata-1)+ndata
pardrift(1)=alpha
pardiff(1)=beta







call GenerateData(sd,ndata,M_all,nproc,data_all,datare)

!print*,sd(2),datare(2,1:(sd(2)+1))
!stop



!print*, datare

call MCMC(seed,datare,ndata,sd,ns,nproc,ndata_all,ndrift,pardrift,ndiff,pardiff,delta,num_iter,& 
        burn,M_beta,Infi,alpha,beta,mu,sigma,gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2)
 

print*,"end"
!desallocate(datare) 
end program OU_mean 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
  subroutine GenerateData(sd,ndata,M_all,nproc,data_all,datare)
  implicit none
  integer ndata,nproc,sd(nproc),M_all,ini,end,i
  real*8 datare(nproc,ndata),data_all(M_all)
    
    ini=1
    
    do i=1,(nproc+1)
      

      end=ini+sd(i)-1
    if (i.eq.301) goto 111
      datare(i,1:sd(i))=data_all(ini:end)
111   ini=end+1
      
    enddo
      
  return
  end subroutine 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalvar(seed,x)
implicit none
real*8 x
integer seed

!seed=765432

 call box(seed,x)
!print*,x

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine box(seed,nor)
        integer seed
        real*8 boxmuller,u1,u2,nor
        integer iset
        real*8 fac,gset,rsq,v1,v2
        save iset,gset
        data iset/0/




 1      if (iset.eq.0) then
        u1=rand()
        u2=rand()
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

        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sim_ou_exac(seed,startx,ndata,pardrift,ndrift,pardiff,ndiff,delta,data)
implicit none
integer i,ndata,ndrift,ndiff,seed
real*8 pardrift(ndrift),pardiff(ndiff),delta,alpha,m,sd,data(ndata),var,beta,startx,a
  
  data(1)=startx
  alpha=pardrift(1)
  a=pardrift(2)
  beta=pardiff(1)
  
  do i=2,ndata

    call normalvar(seed,var)  
    sd=beta*sqrt((1-exp(-2*alpha*delta))/(2*alpha))
    m=a/alpha+exp(-alpha*delta)*(data(i-1)-a/alpha)
    data(i)=m+sd*var

  enddo
  
 
  return
end subroutine


       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MCMC(seed,datare,ndata,sd,nsteps,nproc,ndata_all,ndrift,pardrift,ndiff,pardiff,delta,niter,&
                  burn,M_beta,Infi,alpha,beta,mu,sigma,gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2)
implicit none
integer i,ndata,j,nsteps,ndrift,ndiff,numrej,nproc,mc,init,end,ndata_all,niter,burn,k,m,sd(nproc),M_beta,t1_b,t1_s,seed,g1,b1
real*8 datare(nproc,ndata),db,bridges_all(nproc,ndata_all),bridge(nsteps+2),Y(nproc,ndata_all),Infi,data_lam(nproc,ndata)
real*8 alpha,mu,pardrift(ndrift),pardiff(ndiff),beta,sigma,theta_mu,delta,t2_b,t_a,t1_m,t2_m,t2_s,KS(nproc,9)
real*8 alphas(niter),betas(niter),mus(niter),sigmas(niter),redrift2(nproc),KC(nproc,4),KV(nproc,5)
real*8 gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2,mean_re,var_re,ba
character*20 answer,file1,file2,file3,file4



file1="alphasc.txt"

file2="betasc.txt"
file3="musc.txt"
file4="sigmasc.txt"

! delta bridges 
 open(1,file=file1)
open(2,file=file2)

open(3,file=file3)
open(4,file=file4)

    !prior
    
    !call para_prior(alpha,beta,mu,sigma,t_a,t1_b,t2_b,t1_m,t2_m,t1_s,t2_s)
 call prior(seed,gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2,nproc,redrift2,alphas(1),betas(1),mus(1),sigmas(1))
    alphas(1)=18.0
    betas(1)=0.02
!    mus(1)=0.2
!    sigmas(1)=0.5
    print*,"initial parameters"," ","alpha",alphas(1),"beta",betas(1),"mu",mus(1),"sigma",sigmas(1)!,redrift2

 
 
    do i=2,niter
   
  
   ba=MAX(betas(i-1),0.35)
!print*,"iteration",i

    call YSR(seed,datare,ndata,sd,nsteps,nproc,ndata_all,ndrift,alphas(i-1),redrift2,ndiff,ba,delta,Y)
   
! print*,"here" 

    call suf_sta(bridges_all,datare,sd,nproc,ndata,ndata_all,delta,nsteps,Y,KS)

 call draw_alpha_alt(Infi,redrift2,nproc,ndata,betas(i-1),KS,delta,alpha_2,alphas(i))
!     print*,"alpha",alphas(i)

   
call draw_mu(seed,redrift2,nproc,sigmas(i-1),mu_0,sigma_0,mus(i))
!print*,"mu",mus(i)
g1=gamma_1
call draw_sigma(seed,redrift2,nproc,mus(i-1),g1,gamma_2,sigmas(i))
!print*,"sigma",sigmas(i)
b1=beta_2
call draw_beta_normal_new(Infi,redrift2,nproc,ndata,alphas(i),KS,b1,betas(i),K,M_beta)

!betas(i)=1.5

call draw_as(seed,nproc,ndata,alphas(i),betas(i),mus(i),sigmas(i),KS,redrift2)

mean_re=SUM(redrift2)/nproc
!print*,redrift2
!stop
   var_re=sqrt(SUM((redrift2-mean_re)**2)/nproc)!-mean_re**2

print*,"iter",i,"alpha",alphas(i),"beta",betas(i),"mu",mus(i),"sigma",sigmas(i),"mean-RE",mean_re,"var_re",var_re


    enddo


 write(1,*) alphas
 write(2,*) betas

 write(3,*) mus
 write(4,*) sigmas






      endfile(1)
      close(1)
      endfile(2)
      close(2)
      endfile(3)
      close(3)
      endfile(4)
      close(4)


!call system("R CMD BATCH /Users/bladt/Dropbox/SDE_RE/f/Random_Effect/OU/case_2/sta.r sta.out")
!call system("more sta.out")
!call system("open /Applications/Skim.app /Users/bladt/Dropbox/SDE_RE/MCMC/Random_Effect/OU/case_2/esti.pdf")

 return
end subroutine





!prior
subroutine prior(seed,gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2,nproc,redrift2,alpha,beta,mu,sigma)
implicit none
integer i,nproc,theta1_beta,t1_sigma,seed,b1,g1
real*8 gamma_1,gamma_2,mu_0,sigma_0,beta_1,beta_2,alpha_2
real*8 theta_alpha,theta2_beta,theta1_mu,theta2_mu,t2_sigma
real*8 alp(1),bet(1),m(1),sig(1),beta,redrift2(nproc),var,mu,alpha,sigma

!call para_prior(alpha_par,beta_par,mu_par,sigma_par,theta_alpha,theta1_beta,theta2_beta,theta1_mu,theta2_mu,t1_sigma,t2_sigma)

!print*,theta_alpha,theta1_beta,theta2_beta,theta1_mu,theta2_mu,t1_sigma,t2_sigma


!prior alpha
  call rgammas(seed,1,1,alpha_2,alp)
  alpha=alp(1)

!prior beta 
b1=beta_1
!print*,b1,beta_2
call rgammas(seed,1,b1,beta_2,bet)
beta=bet(1)
!print*,bet(1)
!stop
!prior mu
call normalvar(seed,var)
 mu=mu_0+sqrt(sigma_0)*var
 mu=0.02
  
!prior sigma
g1=gamma_1+1
call rgammas(seed,1,g1,gamma_2,sig)
sigma=sqrt(1/sig(1))
sigma=0.01

!print*,sigma,sig(1),t2_sigma,g1

! Generate random effects with initial parameters   
do i=1,nproc
  call normalvar(seed,var)
  redrift2(i)=mu+sigma*var
enddo

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine YSR(seed,data,ndata,sd,nsteps,nproc,ndata_all,ndrift,alpha,redrift2,ndiff,beta,delta,Y)
implicit none
integer i,j,sd(nproc),nproc,ndata_all,nsteps,ndrift,ndiff,ndata,seed
real*8 matbri(nproc,ndata_all),ls(nproc,ndata_all),Y(nproc,ndata_all)
real*8 data(nproc,ndata),alpha,redrift2(nproc),beta,db,delta

db=delta/(nsteps+1)
call mat_bridge(seed,data,ndata,sd,nsteps,nproc,ndata_all,ndrift,alpha,redrift2,ndiff,beta,db,matbri)

call LSR(data,nproc,ndata,nsteps,ndata_all,ndiff,beta,delta,sd,ls)

Y=matbri-ls
   
return
end  

!! complete path with bridges (data,ndata,nsteps,ndrift,pardrift,ndiff,pardiff,del_bri,bridges)

subroutine mat_bridge(seed,datare,ndata,sd,nsteps,nproc,ndata_all,ndrift,alpha,re,ndiff,beta,db,bridges_all)
implicit none
integer i,ndata,sd(nproc),j,nsteps,ndrift,ndiff,numrej,nproc,k,mc,init,end,ndata_all,n,seed
real*8 datare(nproc,ndata),pardrift(ndrift),re(nproc),pardiff(ndiff),db,bridges_all(nproc,ndata_all),bridge(nsteps+2)
real*8 beta,var,delta,alpha,Y(nproc,ndata_all),Z(nproc,ndata_all),ls(nproc,ndata_all),t(nsteps+2),t0,tn,Ypartial(nsteps+2)
real*8 KS(nproc,9)
   
   
    pardiff(1)=beta
    do i=1,nproc
    pardrift(1)=alpha
    pardrift(2)=re(i)
   
   
!print*,"i",i,"drift",redrift(i),"diffusion",rediff(i),"var",var
    do j=1,sd(i)-1
     init=nsteps*(j-1)+j
     end=init+nsteps+1

          call FormerDiffusionBridge(seed,ndrift,pardrift,ndiff,pardiff,db,nsteps,datare(i,j),datare(i,j+1),bridge,numrej)
          bridges_all(i,init:end)=bridge
          


          !Y(i,init:end)=Ypartial
 
    
    enddo
 !   print*,i 
    enddo

!stop


return
end subroutine


 ! Make a bridge from x to y
  subroutine FormerDiffusionBridge(seed,ndrift,pardrift,ndiff,pardiff,delta,numsteps,x,y,bridge,numrej)
  implicit none
  integer i,numrej,j,mp,numsteps,ndrift,ndiff,seed
  real*8 points(numsteps+2),x,y,pardrift(ndrift),pardiff(ndiff),delta,ori_delta,xl,yl
  real*8 points1(numsteps+2),points2(numsteps+2),aux,t0,tn,T(numsteps+2)
  real*8 points3(numsteps+2),bridge(numsteps+2),pardiff_l(ndiff),pardrift_l(ndrift)





ori_delta=delta*(numsteps+1)

 

call sde_lamperti(x,y,ndrift,pardrift,ndiff,pardiff,xl,yl,pardrift_l,pardiff_l)  


t0=0.0
tn=ori_delta
  call sequence(t0,tn,numsteps,T)
 
  numrej=0
!print*,nnd,driftnd,ndiff,pardiff,delta,x,numsteps
1     call diffusion(seed,ndrift,pardrift_l,ndiff,pardiff_l,delta,xl,numsteps,points1)
!print*,"p1",points1


2     call diffusion(seed,ndrift,pardrift_l,ndiff,pardiff_l,delta,yl,numsteps,points2)
!print*,"p2",points2


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
!  print*, numrej,xl,yl
      goto 1

20    do i=1,(mp-1)
       points(i)=points1(i)
      enddo

       do i=1,numsteps+2

          bridge(i)=points(i)
          enddo

 return
  
  end


!! SDE with Lamperti transformation 

subroutine sde_lamperti(x,y,ndrift,pardrift,ndiff,pardiff,xl,yl,pardrift_lam,pardiff_lam)  
implicit none
integer ndrift,ndiff,ndrift_lam,ndiff_lam
real*8 x,y,xl,yl,pardiff(ndiff),pardrift(ndrift),pardrift_lam(ndrift),pardiff_lam(ndiff)

xl=x/pardiff(1)
yl=y/pardiff(1)
pardiff_lam(1)=1.00
pardrift_lam(2)=pardrift(2)/pardiff(1)
pardrift_lam(1)=pardrift(1)

return
end

! sequence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


! Using the Milstein scheme we simulate a diffusion from startx, n steps ahead at
! stepsizes delta. Diffusion depends on external functions
      subroutine diffusion(seed,ndrift,pardrift,ndiff,pardiff,delta,startx,nsteps,points)
  implicit none
  integer nsteps,i,ndrift,ndiff,seed
      real*8 delta,startx,y1,y2,pardrift(ndrift),pardiff(ndiff),y3
     real*8 points(nsteps+2)

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


       call MilsteinStep(seed,delta,points(i-1),y1,y2,y3,points(i))
!print*,'step2',points(i)


    enddo
!print*,points(1)



  return

  end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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







!!!!

      subroutine MilsteinStep(seed,delta,startx,alpha,sigma,sigmax,endx)
      implicit none
      integer seed
      real*8 delta,startx,endx
      real*8 alpha,sigma,W,ini,sigmax

!       call normalvar(W)
!       endx=startx*exp(-1.0*delta)+sqrt(1.0-exp(-2.0*delta))*W
ini=0.0
  call BrownianStep(seed,delta,ini,W)
 ! print*,startx,alpha,delta,sigma,W
  endx=startx+alpha*delta+sigma*W+(1/2)*sigma*sigmax*(W**2-delta)

!print*,endx
      return
  end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine BrownianStep(seed,delta,startx,endx)
  implicit none
  integer seed
  real*8 delta,startx,endx,var

        call normalvar(seed,var)
!print*,'var',var
        endx=startx+sqrt(delta)*var
!print*,'endx',endx
  return
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine LSR(data,nprocess,ndata,nsteps,ndata_all,ndiff,pardiff,delta,sd,ls)
implicit none
integer i,j,k,nprocess,ndata,nsteps,ndata_all,ini,ndiff,sd(nprocess),init,size
real*8 data(nprocess,ndata),ls(nprocess,ndata_all),delta,pardiff,l
real*8 data_lam(nprocess,ndata),time1(nprocess,ndata_all),time2(nprocess,ndata_all),dl1,dl2

 call Lamperti(data,nprocess,ndata,ndiff,pardiff,data_lam)
 call Dif_Times(nsteps,delta,sd,ndata_all,nprocess,time1,time2)
  !print*,time1,time2


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine suf_sta(bridges_all,redata,sd,nproc,ndata,ndata_all,delta,nsteps,Y,KS)
implicit none
integer nproc,ndata_all,i,ndata,j,nsteps,k,sd(nproc),size
real*8 KS(nproc,9),redata(nproc,ndata),bridges_all(nproc,ndata_all),delta,del_bri,Y(nproc,ndata_all)
real*8 l1(nproc,ndata_all)

del_bri=delta/(1.0+nsteps)
 call L1R(redata,nproc,ndata,sd,nsteps,ndata_all,delta,l1)
KS(:,:)=0.00
  do i=1,nproc
size=nsteps*(sd(i)-1)+sd(i)
    KS(i,1)=redata(i,sd(i))-redata(i,1)
    
    
    KS(i,2)=(redata(i,sd(i))**2-redata(i,1)**2)/2
    KS(i,3)=0.0
    KS(i,4)=(sd(i)-1)*delta/2
    KS(i,5)=0.0
    KS(i,6)=0.0
    KS(i,7)=0.0
    KS(i,8)=0.0
    KS(i,9)=0.0

    do j=2,sd(i)
      KS(i,3)=KS(i,3)+(redata(i,j)-redata(i,j-1))**2
    enddo
      KS(i,3)=KS(i,3)/(2*delta)
     
    do k=2,size
        KS(i,5)=KS(i,5)+(del_bri)*(Y(i,k)+Y(i,k-1))/2
        KS(i,6)=KS(i,6)+(del_bri)*(l1(i,k)+l1(i,k-1))/2
        KS(i,7)=KS(i,7)+(del_bri)*(Y(i,k)**2+Y(i,k-1)**2)/2
        KS(i,8)=KS(i,8)+(del_bri)*(Y(i,k)*l1(i,k)+Y(i,k-1)*l1(i,k-1))/2
         KS(i,9)=KS(i,9)+(del_bri)*(l1(i,k)**2+l1(i,k-1)**2)/2
        enddo
      KS(i,7)=KS(i,7)/2 
      KS(i,9)=KS(i,9)/2 
  enddo
  end subroutine





! Lamperti transformation

subroutine Lamperti(data,nprocess,ndata,ndiff,pardiff,data_lam)
implicit none
integer nprocess,ndata,ndiff
real*8 data(nprocess,ndata),data_lam(nprocess,ndata),pardiff

data_lam=data/pardiff
return
end subroutine


 


! generate gammas integer parameter

 subroutine rgammas(seed,size,n,alpha,rg)
 implicit none
integer size,n,i,j,seed
real*8 alpha,rg(size),u

rg(:)=0.0
!print*,n
do i=1,size
  do j=1,n
    u=rand()

    rg(i)=rg(i)-log(u)/alpha

!    print*,i,rg(i)
  enddo
enddo

return

end subroutine


subroutine draw_alpha_alt(Infi,redrift2,nproc,ndata,beta,KS,delta,theta_alpha,alpha)
implicit none

integer nproc,ndata,i,seed
real*8 redrift2(nproc),theta_alpha,alpha,beta,KS(nproc,9),delta,alpha_1,alpha_2,sd,u,al
real*8 pp,xt,cdn,tr,m,var,Infi,alp(1),min
LOGICAL upper
min=0.00

!seed=765432

  alpha_1=SUM(KS(:,7))+(1/beta)*SUM(KS(:,8))+(1/(beta)**2)*SUM(KS(:,9))
  
  
  alpha_2=-theta_alpha-(1/(beta)**2)*SUM(KS(:,2))+SUM(KS(:,4))

  do i=1,nproc
  alpha_2=alpha_2+redrift2(i)*KS(i,5)/beta+redrift2(i)*KS(i,6)/((beta**2))
  enddo
  print*,"alpha1",alpha_1,"alpha2",alpha_2
  m=alpha_2/(2*alpha_1)
  var=1/(2*alpha_1)
  sd=sqrt(var)

  


 call truncated_normal_ab_sample(m,sd,min,Infi,seed,alpha)



    return
    end











subroutine draw_mu(seed,redrift2,nproc,sigma,t1_mu,t2_mu,mu)
implicit none

integer nproc,i,seed
real*8 redrift2(nproc),t1_mu,t2_mu,sigma,mu,A,B,var,m,sd


  A=1/(t2_mu**2)+nproc/(sigma**2)
  B=SUM(redrift2)/(sigma**2)+t1_mu/(t2_mu**2)
!print*,A,B
!stop
 call normalvar(seed,var)
 m=B/A
 sd=sqrt(1/A)
 mu=m+sd*var



    return
    end subroutine




subroutine weight_beta(B,x,y,crit)
implicit none


real*8 B,x,y,crit


crit=B*(sqrt(x)-sqrt(y))

    return
    end subroutine




subroutine draw_beta_normal_new(Infi,redrift2,nproc,ndata,alpha,KS,theta_beta,beta,K,M)
implicit none                
integer nproc,ndata,i,n,M,parar,K,num_acep,seed,theta_beta
real*8 redrift2(nproc),KS(nproc,9),nr,e(1),eta,mt,st,bb,eta_new,eta_old,Infi
real*8 alpha,beta,delta,var(M),A,B,g(M),s,p(M),u,crit,K1,cond,mu,sigma2,sigma,Z,pp,al,acep,x,C(6)
logical upper
real*8 tr,cdn,min
character*20 ans

min=0.00


  A=0.0
  B=0.0

  do i=1,nproc
    A=A-KS(i,1)*redrift2(i)+alpha*KS(i,2)+KS(i,3)+(redrift2(i)**2)*KS(i,4)-alpha*redrift2(i)*KS(i,6)+(alpha**2)*KS(i,9)
    B=B+alpha*redrift2(i)*KS(i,5)-(alpha**2)*KS(i,8)
   
  enddo


  B=B-theta_beta

  print*,"A",A,"B",B
  n=nproc*(ndata-1)
  mu=B/(2*A)
  sigma2=1/(2*A)
  sigma=sqrt(sigma2)
!


num_acep=0

!call cdnormal(tr,cdn)
!print*,"cdn",cdn
!call PHIE(al,pp)
!print*,"pp",pp

st=sigma2
mt=(mu+sqrt(mu**2+4*n*st))/2!+sigma*(pp/(1.0-cdn))
!*(1+al*pp/(1.0-cdn)-(pp/(1.0-pp))**2)
print*,mt,st
!stop

!call normalvar(var(1))
!eta_old=mt+sqrt(st)*var(1) 
!print*,"first",eta_old
 !call trun_normal(1,mt,st,0.0,Infi,eta_tr(1))

 call truncated_normal_ab_sample(mt,sigma,min,Infi,seed,eta_old)

!print*,"second",eta_old
!print*,(mt+sqrt(mt**2+4*n*st))


do i =1,M

!call normalvar(var(1))
!eta_new=(mt+sqrt(mt**2+4*n*st))/2+sqrt(st)*var(1) 

!call trun_normal(1,mt,st,0.0,Infi,eta_tr(1))

call truncated_normal_ab_sample(mt,sigma,min,Infi,seed,eta_new)
 
acep=n*log(eta_new/eta_old)-(1/(2*st))*((eta_old-eta_new)*(mu-sqrt(mu**2+4*n*st)))
call random_number(u) 
!print*,"old",eta_old,"new",eta_new,acep

if(log(u)<acep) then 
  eta_old=eta_new
  num_acep=num_acep+1
 ! print*,"ny",num_acep,eta_old
endif

enddo
print*,"num_acep",num_acep
beta=1/eta_old
!print*,"beta",beta
 

    return
    end subroutine



subroutine draw_sigma(seed,redrift2,nproc,mu,theta1_sigma,theta2_sigma,sigma)
implicit none

integer nproc,i,n,theta1_sigma,seed
real*8 redrift2(nproc),mu,sigma,var(1),C6,theta2_sigma,ms,ms1


  
  C6=0.0
  do i=1,nproc
  C6=C6+(redrift2(i)**2)/2-mu*redrift2(i) 
  enddo
  C6=C6+nproc*(mu**2)/2+theta2_sigma
  n=nproc/2+theta1_sigma
  
  call rgammas(seed,1,n,C6,var(1))

 sigma=sqrt(1/var(1))
ms=sqrt(C6/n)
ms1=sqrt(C6/(n-1))
print*,ms,ms1


    return
    end subroutine





subroutine draw_as(seed,nproc,ndata,alpha,beta,mu,sigma,KS,redrift2)
implicit none
integer nproc,ndata,i,seed
real*8 redrift2(nproc),sigma,beta,KS(nproc,9),var,A(nproc),B(nproc),alpha,mu,m,sd,dd



  
  
  do i=1,nproc
    A(i)=KS(i,4)/(beta**2)+1/(2*sigma**2)
    B(i)=mu/(sigma**2)+KS(i,1)/(beta**2)+(alpha*KS(i,5))/(beta)+alpha*KS(i,6)/(beta**2)
  
  m=B(i)/(2*A(i))
  sd=sqrt(1/2*A(i))

  call normalvar(seed,var)

  redrift2(i)=m+sd*var
 
  enddo
  
  !print*,SUM(A)
  





    return
    end


! Auxiliar Y's 

!!!!!!!!!(delta,datare(i,j),datare(i,j+1),nsteps,pardrift(1),pardrift(2),pardiff(1),bridge,Ypartial)

subroutine bridge_ou_exact(seed,delta,x0,xn,n,alpha,a,sigma,Z,Y)
implicit none
integer i,n,j,seed
real*8 x0,xn,T(n+2),alpha,sigma,X(n+2),Z(n+2),norval,sd,t0,tn,delta,a,Y(n+2)

t0=0.0
tn=delta
!print*,delta

  call sequence(t0,tn,n,T)

  X(1)=x0
  Z(1)=x0

  i=2
  do while(i<=n+2)
    sd=sigma*sqrt((1-exp(-2*alpha*(T(i)-T(i-1))))/(2*alpha))
    call normalvar(seed,norval)
    X(i)=a/alpha+exp(-alpha*(T(i)-T(i-1)))*(X(i-1)-a/alpha)+norval*sd
    i=i+1
  enddo

  j=2
 do while(j<=n+2)
    Z(j)=X(j)+(xn-X(n+2))*(exp(alpha*T(j))-exp(-alpha*T(j)))/(exp(alpha*T(n+2))-exp(-alpha*T(n+2)))
    
    Y(j)=Z(j)/sigma-((1.0-T(j))*(x0/sigma)+T(j)*(xn/sigma))/delta

    j=j+1
  enddo
  Y(1)=Z(1)-((1-T(1))*(x0/sigma))/delta
 return
end subroutine


!!!!!!!!!(delta,datare(i,j),datare(i,j+1),nsteps,pardrift(1),pardrift(2),pardiff(1),bridge,Ypartial)

subroutine bridge_ou_exact_lam(seed,delta,x0,xn,n,alpha,a,sigma,Z,Y)
implicit none
integer i,n,j,seed
real*8 x0,xn,T(n+2),alpha,sigma,X(n+2),Z(n+2),norval,sd,t0,tn,delta,a,Y(n+2),sigmal,al,x0l,xnl

x0l=x0/sigma
xnl=xn/sigma

!print*,"ini",x0,"end",xn

t0=0.0
tn=delta
!print*,delta
al=a/sigma
sigmal=1.0
  call sequence(t0,tn,n,T)

  X(1)=x0l
  Z(1)=x0l

  i=2
  do while(i<=n+2)
    sd=sigmal*sqrt((1-exp(-2*alpha*(T(i)-T(i-1))))/(2*alpha))
    call normalvar(seed,norval)
    X(i)=al/alpha+exp(-alpha*(T(i)-T(i-1)))*(X(i-1)-a/alpha)+norval*sd
    i=i+1
  enddo

  j=2
 do while(j<=n+2)
    Z(j)=X(j)+(xnl-X(n+2))*(exp(alpha*T(j))-exp(-alpha*T(j)))/(exp(alpha*T(n+2))-exp(-alpha*T(n+2)))
    
    Y(j)=Z(j)-((1.0-T(j))*x0l+T(j)*xnl)/delta
!print*,j,Z(j),T(j),x0l,xnl,Y(j)
    j=j+1
  enddo
  Y(1)=Z(1)-((1-T(1))*x0l)/delta
 return
end subroutine





!redata,nproc,ndata,nsteps,ndata_all,delta,l1)

subroutine L1R(data,nprocess,ndata,sd,nsteps,ndata_all,delta,l)
implicit none
integer i,j,k,nprocess,ndata,nsteps,ndata_all,init,sd(nprocess),size
real*8 data(nprocess,ndata),t0,tn,T(nsteps+2),db,l(nprocess,ndata_all),delta
real*8 time1(nprocess,ndata_all),time2(nprocess,ndata_all)

call Dif_Times(nsteps,delta,sd,ndata_all,nprocess,time1,time2)
  !print*,time1,time2

l(:,:)=0.0
  !stop
  do i=1,nprocess
    size=nsteps*(sd(i)-1)+sd(i)
    do j=1,sd(i)-1
    init=nsteps*(j-1)+j
    


    do k=0,nsteps
     !print*,i,k+init
      l(i,init+k)=(time2(i,init+k)*data(i,j)+time1(i,init+k)*data(i,j+1))/delta
     !print*,ls(i,ini+k) 
     
    end do
  end do
  l(i,size)=(time2(i,ndata_all)*data(i,j)+time1(i,ndata_all)*data(i,j+1))/delta
      
  


end do
  
  
  
  return
end subroutine

subroutine MLE(data,n,delta,alpha,beta)
integer n,i
real*8 data(n),alpha,beta,num(n-1),dem(n-1),dif(n-1),s




do i=1,n-1
  num(i)=data(i)*data(i+1)
  dem(i)=data(i)**2

enddo

alpha=-(1/delta)*log(SUM(num)/SUM(dem))
s=0.0
do i=1,n-1
 dif(i)=(data(i+1)-data(i)*exp(-delta*alpha))**2
 s=s+dif(i)
enddo
print*,SUM(dif),s
print*,(1-exp(-2*delta*alpha))*n
beta=sqrt(((2*alpha)*SUM(dif))/(n*(1-exp(-2*delta*alpha))))
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




! end of file normal.f90




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

subroutine trun_normal(seed,n,mean,var,inf,sup,x)
  implicit none
  integer i,n,seed
  real*8 alpha,beta,cdinf,cdsup,xt,pp
  real*8x(n),mean,var,inf,sup,sd,u

!print*,"inputs",n,mean,var,inf,sup
sd=sqrt(var)
  alpha=(inf-mean)/sd
  beta=(sup -mean)/sd
do i=1,n
  
  u=rand()


  call cdnormal(alpha,cdinf)
  call cdnormal(beta,cdsup)
 ! print*,alpha,beta,cdinf,cdsup
  pp=cdinf+u*(cdsup-cdinf)
  if(pp==1.00) then
  pp=0.99999995
  endif  
  call normal_01_cdf_inv(pp,xt)

  !print*,"pp",pp,xt
  x(i)=sd*xt+mean

  !print*,sd*xt
enddo
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

!!! Times

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
 

