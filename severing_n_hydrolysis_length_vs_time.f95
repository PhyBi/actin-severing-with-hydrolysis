!!!!!!!! This Program gives output for length Vs time  for a single filament system    
!!!!!!!! Polymerization at barbed end &  decay & severing from the pointed  end                 
!!!!!!!! Code written by binayak banerjee (binayakbanerjee32@gmail.com)



MODULE numz
    IMPLICIT NONE
    INTEGER,PARAMETER:: DP=KIND(1.0D0)
    REAL(DP),PARAMETER:: pi=3.14159265358979_DP
  END MODULE numz 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parm

    use numz

    implicit none
    integer,parameter:: N=150000,Ttrans=0
    double precision,parameter::h0=0.15d0,h1=1.0d0,GT=0.00d0,GD=0.00d0,s0_boundary=0.00001d0,c=2.50d0,&
            tstop=18000.0d0,tint=1.0d0,r0=11.6d0,ne=0.0d0,s0_bulk=0.0d0,c_actin=1.0d0,nh=3.5d0                                     
    

    ! N= Finite pool size; h0,h1= coarse-grain hydrolysis rate; GT, GD= depolymerization rate base on tip subunit; c= cofilin
    ! concn.; r0= polymerization rate; ne= nucleotide exchange rate in pool(recycling); c_actin= actin concentarion, nh= hill
    ! coefficient; s0_boundary= severing rate when it happen at the junction decorated(NDP+cofilin) and undecrated(NTP) pair;
    ! s0_bulk= severing rate when it happen at the junction of decorated pair; tstop = end time to stop the run

end module parm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program single_filament

    use parm                                                              !Calling the parm module for parameters
    
    implicit none
    integer::l1,NG,F1(N),NT,ND,LT1(N),j,i,code,m1,LD1(N),LD2(N),LBoundary(N),&                                
    p1,idum,n0,n2,m,o,s1,s2,s3,s4,bulk_cnt,boundary_cnt

 
    double precision:: A(6),r1,r2,tau,t,a0,ran2,w,k1,texp,,mean,l2,&                                          
    ti,tf,cms,citr,tprint
    
 
  

    idum=49385



 
    call cpu_time(ti)                                                          !Calculating the CPU time


 
    F1 = 0                                                                    !Initalisation of the filamnet lattice (array)
    l1 = 0                                                                    !Initial length of the filament
    NG = 0                                                                    !Number of ADP subunits in th pool tracker
    t = 0.0d0
    tprint = 8000.0d0                                                         !Time afer which data will be stored  
   

             open( unit = 10, file = "length_vs_time_data.csv")  

    do i =1,l1                                      !Filling the Number of GTP momoners in the filament according to the given intial length
    F1(i)=1
    end do
    citr=0.0d0                                      !Tracker for counting the number of iterations
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

     do                                                                                  !Starting the main loop



        CALL SYSTEM_CLOCK(COUNT=idum)              ! Calling the seed value for random  number generator based on the cpu time
        LT1 = 0                                    ! ATP array: collects the position of ATP(+1) in the filament
        LD1 = 0                                    ! ADP array: collects the position of ADP(-1) in the filament 
        LD2 = 0                                    ! Keep track of consecutive -1
        LBoundary = 0
        NT = 0                                     ! Tracker for counting the number of ATP monomers in the lattice (array) of the filamnet 
        ND = 0                                     ! Tracker for counting the number of ADP monomers in the lattice (array) of the filamnet
        bulk_cnt = 0                               ! bulk_cnt is the no. of possible severing sites i.e -1 -1 boundary   
        boundary_cnt = 0                           ! count boundary nos of +1 and -1 and vice versa
        cap=0

      if( l1 > 0 ) then                                ! Tracking the Number of ATP only when filament length is >0
            do i = 1,l1
                if (F1(i)==1)then
                    NT = NT + 1
                    LT1(NT)=i
                else if (F1(i) == -1) then             !Tracking the Number of ADP only when filament length is >0
                    ND = ND + 1
                    LD1(ND)=i
                end if
            end do
     end if
  
     if( l1 > 1 ) then    
        do i=1,ND-1
         if (LD1(i+1)-LD1(i) .eq. 1) then
                 bulk_cnt = bulk_cnt+1                ! bulk_cnt is the no. of possible severing sites  
                 LD2(bulk_cnt) = LD1(i)               ! Keep track of consecutive -1       
         end if
        end do        
       
                do i=1,l1-1
                  if (F1(i) .ne. F1(i+1)) then
                       boundary_cnt = boundary_cnt + 1
                       LBoundary(boundary_cnt) = i               ! Keep track of  boundary of +1 and -1 and vice versa
                   end if
                 end do
       end if




        A(1) = r0*c_actin                         !r0*(N-l1-NG) (for finite pool)     !Addition
        A(2) = ((h0*c**nh)/(h1 + c**nh))*(NT)                                         !Hydrolysis propensity
        A(4) = s0_bulk*(bulk_cnt)                                                     !Severing propensity  
        A(5) = s0_boundary*boundary_cnt
        A(6) = ne*(NG)                                                                !Nucleotide exchange in Pool


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if (l1 > 0 .and. F1(l1) == 1 ) then          !Conditioning the rates with the length values
            A(3) = GT
            
        else if ( l1 > 0 .and. F1(l1) == -1) then
            A(3) = GD
            
        else if (l1 <= 0) then
            A(2) = 0
            A(3) = 0
            A(4) = 0
            A(5) = 0
            A(1) = r0*c_actin
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        a0 = sum(A)
        cms = 0.0d0
        


        r1 = ran2(idum)
        tau = (-1.0d0)*(1/a0)*DLOG(r1)             !Using Gillsepie to get the case value for selection of the reactions and time
        t = t+tau                                  
 
        if (t .gt. tstop) exit                     ! Data storing at equal interval   
22      if (t .lt. tprint) go to 25
        texp=texp+1
        write(10,*) tprint,l1
        tprint=tprint+tint
        go to 22


25      w = ran2(idum)


        do j=1,6
            code = j
            cms = cms+A(j)
            if (cms >= w*a0) exit
        end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            select case(code)                           !case function
            case(1)
                l1 = l1 + 1                             !Growth 
                if (l1 == 1) then
                F1(l1) = 1                              !Inserting a 1 in the lattice symbolising the presence of  a GTP momner
                else
                F1 = eoshift ( F1, shift = -1)          ! Shifts whole array F1 right by 1   
                F1(1) = 1
                end if
!                print*,'growth'

           case(2)
                k1 = ((NT - 1)*ran2(idum) + 1)           !Calling a random number [1,NT], based upon the total number of GTP monomer in the lattice
                m1 = nint(k1)
                p1 = LT1(m1)                            !Hydrolysis
                F1(p1) = -1                             !Inserting a -1 in the lattice symbolising the presence of  a GDP momner
!                print*,'hydrolysis'

            case(3)                                     !Decay
                if ( F1(l1) == -1) then
                    NG = NG + 1
                end if

                F1(l1) = 0                               !adding a zero in the lattice
                l1 = l1 - 1
                if (A(3) == GD) then
                NG = NG + 1                             !Decay leading to an incrase of GDP monomer in the pool if the decayed subunit was a GDP
                end if
!                print*,'Decay'

            case(4)
                  r2 = ran2(idum)                       ! Severing in bulk
                  m = nint(r2*(bulk_cnt - 1) + 1)       ! Choose random no. betn [1,bulk_cnt]
                  s1 = LD2(m)                           ! s1 is the site where filament breaks  

                  do i=s1+1,l1                            
                  F1(i) = 0
                  end do
        
                  l1 = s1                               ! After severing filament length is s1

                  do i=1,ND-1
                  if  (LD1(i) .eq. s1) then             ! Count no. of ADP right to the severing site     
                      s2 = i
                      exit
                  end if
                  end do

                  NG =  NG + (ND - s2)
                
!                  print*,'severing in bulk'

             case(5)                                                !severing in boundary
                  o = nint(ran2(idum)*(boundary_cnt - 1) + 1)       ! Choose random no. betn [1,boundary_cnt]
                  s3 = LBoundary(o)                                 ! s3 is the site where filament !breaks
                  s4 = 0                                    ! initialisation; counts -1 right to the severing site
                               

                  do i=s3+1,l1    
                  if (F1(i) == -1) then
                         s4 = s4 +1
                         F1(i) = 0
                  else 
                         F1(i) = 0
                  end if 
                  end do
                  
                  l1 = s3
                 
                  NG =  NG +  s4   
!                  print*,'severing in boundary'
  

             case(6)                                                !Nucleotide exchange in pool 
                     NG = NG - 1
!                     print*,'nucleotide exchange in pool'

             case default
                print*,"Error in code !!! ","Code=",code

        end select
        
               
        


      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1      
close(10)

call cpu_time(tf)
end program single_filament

!!!!!!!!!!!!!!!               FUNCTION AND SUBROUTINES USED    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Uniform Random number generators

FUNCTION ran2(idum)
    USE numz
    IMPLICIT NONE
    REAL(DP):: ran2
    !INTEGER,INTENT(inout),OPTIONAL::idum
    INTEGER,INTENT(inout)::idum
    INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
    INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
    INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
    INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
    REAL(DP),PARAMETER::AM=1.0_DP/IM1,EPS=1.2e-7,RNMX=1.0_DP-EPS
    INTEGER::idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    IF (idum<0) THEN
       idum=MAX(-idum,1)
       idum2=idum
        DO j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           IF (idum<0) idum=idum+IM1
           IF (j.LE.NTAB) iv(j)=idum
        ENDDO
        iy=iv(1)
     ENDIF
     k=idum/IQ1
     idum=IA1*(idum-k*IQ1)-k*IR1
     IF (idum<0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
     IF (idum2<0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
     IF(iy.LT.1)iy=iy+IMM1
     ran2=MIN(AM*iy,RNMX)
     RETURN
   END FUNCTION ran2
   
