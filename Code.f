cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program is to calculate the input parameters for all the Tsyganenko   c
c     magnetic field models,including G parameters for T01 and T01s, W           c
c     parameters for TS05.                                                       c
c     By Zhengrui Qin and Richard Denton, 09/01/2007                             c
c     Original code from                                                         c
c     http://virbo.org/svn/virbo/qindenton/                                      c
c           MagParameterProgram-rsw/MagmodelinputONE_corrected.f                 c
c     LANL updated version, 2015 (fix for more recent compilers in 2017)         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc   Choose data format of the OMNIWEB data(http://omniweb.gsfc.nasa.gov/) ccc 
c      write(*,*)'c    IFORMAT = 1 for one minute data                 c'
c      write(*,*)'c    IFORMAT = 5 for five minute data                c'
c      write(*,*)'c    IFORMAT = 60 for hourly data                    c'
c      write(*,*)'cccccccccccccccccccccccccccccccccccccccccccccccccccccc'

      character ychar*4,fchar*24,filein*32,dumchar*19
      real dmy(100)
      common /years/year0,nyrs
      common wst(6)
      integer iyear

      write (*,*) ' please input cadence, starting year, no of years'
      read (*,*) iformat,year0,nyrs
      write (*,111) iformat,year0,year0+nyrs-1
 111  format (' input: cadence',i4,' min  starting year',f6.0,
     & ' ending year',f6.0)

      if(iformat.ne.1.and.iformat.ne.5.and.iformat.ne.60) then
       write(*,*)' cadence must be 1 or 5 or 60; run aborted'
       stop
      endif
      if (iformat.eq.1.or.iformat.eq.5) then
       if (year0.lt.1995) then
        write (*,*)' starting year must be after 1994; run aborted'
        stop
       end if
      end if
      if (iformat.eq.60) then
       if (year0.lt.1963) then
        write (*,*)' starting year must be after 1962; run aborted' 
        stop
       end if
      end if
c     if (year0+nyrs.gt.2016) then
c      write (*,*)' last year cannot be after 2015; run aborted'
c      stop
c     end if

      iyear=year0-1
      write (ychar,'(i4)') iyear
      fchar=ychar//'/QinDenton_'//ychar//'1231_'
      filein=fchar//'1min.txt'
      write (*,*) 'searching ',filein
      open (unit=1,file=filein,err=20,status='old')       ! search for W1-6 pars in previous file
      do j=1,200
       read (1,*,end=20)                                  ! skip header (192 lines)
      end do
 90   read (1,*,err=20,end=30) dumchar,(dmy(i),i=1,31),(wst(j),j=1,6)
       go to 90
 20   filein=fchar//'5min.txt'
      write (*,*) 'searching ',filein
      open (unit=2,file=filein,err=21,status='old')
      do j=1,200
       read (2,*,end=21)                  
      end do
 91   read (2,*,err=21,end=30) dumchar,(dmy(i),i=1,31),(wst(j),j=1,6)
       go to 91
 21   filein=fchar//'hour.txt'
      write (*,*) 'searching ',filein
      open (unit=3,file=filein,err=22,status='old')
      do j=1,200
       read (3,*,end=22)                  
      end do
 92   read (3,*,err=22,end=30) dumchar,(dmy(i),i=1,31),(wst(j),j=1,6)
       go to 92
 22   write (*,*) ' did not find previous WG file, will use averages'
      wst(1) = 0.44    ! at the average values for status 2 in table 3
      wst(2) = 0.42    ! of Qin,Denton,Tsyganenko,Wolf article
      wst(3) = 0.66    ! in Space Weather, vol 5, issue 11
      wst(4) = 0.48    ! http://onlinelibrary.wiley.com/doi/10.1029/2006SW000296/full
      wst(5) = 0.49   
      wst(6) = 0.91

 30   write (*,130) wst(1),wst(2),wst(3),wst(4),wst(5),wst(6)
 130  format ('  initial Ws :',6f7.3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(iformat.eq.1) then
              call kp3_min
              call kpdst_min
              call bsgamma_min
              call av1day_min
              call inter1day_min
              call interpl_min
              call WG_min
      elseif(iformat.eq.5) then
              call kp3_5min
              call kpdst_5min
              call bsgamma_5min
              call av1day_5min
              call inter1day_5min              
              call interpl_5min
              call WG_5min
      else
              call kp3bsgamma
              call av1year
              call av20days
              call inter20days
              call interhour
              call WGhour
      endif
      write(*,*)'     ALL DONE'
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c  *************  FOR 1 MIN DATA ***************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc subroutine kp3  ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c This program calculates the 3 day average of Kp: Kp3.

      subroutine  kp3_min
      integer dst,cap
      character finput*30,fout*30
      real Kp,akp3,Kpa(120),akp,weight,t3day 
      common wst(6)
      common /years/year0,nyrs
           
      finput='1min/kpdst.lst' 
      fout='1min/kpkp3dst.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34

 33   write(*,*) finput//' NOT EXIST'
      goto 1000

 34   open(unit=2,file=fout,status='unknown') 
 
      WRITE(2,*)'YEAR DOY HR   Kp     akp3  DST' 
 
      do I=1,120 
         Kpa(I)=0. 
      enddo 
      J=0      

10    read(1,*,err=20,end=20) iyr4,idoy,ihr,Kp,dst
       IF (iyr4.lt.year0.or.iyr4.ge.year0+nyrs) go to 10  ! Alin
c       IF (iyr4.lt.year0) go to 10  ! Alin
      J=J+1
      if(Kp.eq.99) then
         akp3=99.
      else
         Kp = Kp/10. 
         if (J.lt.120) then 
            cap = J 
         else 
            cap = 120 
         endif 
         do K=1,119 
            Kpa(K)=Kpa(K+1) 
         enddo 
         Kpa(120)=Kp 
C    Calculate akp3 
         
         akp = 0 
         weight = 0 
         t3day = 72         !Time for 3 days in hours 

         do ii=1,cap 
            weight = weight + exp(-((ii-1)/t3day)) 
            akp = akp + exp(-((ii-1)/t3day))*Kpa(121-ii) 
         enddo
 
         akp3= akp/weight 
      endif
   
      WRITE(2,444)iyr4,idoy,ihr,Kp,akp3,dst
 444  FORMAT(I4,1x,I3,1x,I2,1x,2(F6.2,1x),i6) 
      goto 10
 
 20   close(1) 
      CLOSE(2) 
1000  write(*,*)'     Subroutine kp3_min        DONE'
      END  
cccccccccccccc end of subroutine kp3_min  ccccccccccccccccccccccccccccccccccccccccc


ccccccccccccc subroutine kpdstmin   ccccccccccccccccccccccccccccccccccccccccccccccc

c This program interpolate the hourly kp,kp3 and dst to 1 min data.
      subroutine kpdst_min
      common wst(6)
      character finput*30,fout*30
      real kp1,kp31,kp2,kp32,kp,kp3
      integer dst1,dst2,dst 

      finput='1min/kpkp3dst.d' 
      fout='1min/kpkp3dstmin.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34
 33   goto 1000
 34   open(unit=2,file=fout,status='unknown') 
 
      WRITE(2,*)'YEAR DOY HR MIN  Kp     kp3     DST' 
      read(1,*)
      read(1,*,err=20,end=20) iyr1,idoy1,ihr1,kp1,kp31,dst1
      do i=0,29
          write(2,444)iyr1,idoy1,ihr1,i,kp1,kp31,dst1
444   FORMAT(I4,1x,I3,1x,I2,1x,i2,1x,2(F6.2,1x),i6) 
      enddo          

10    read(1,*,err=20,end=20) iyr2,idoy2,ihr2,kp2,kp32,dst2
      do i=30,89
          f2=(i*1.0-29.5)/60.0
          f1=1-f2
          
          if(kp1.eq.99.00 .or. kp2.eq.99.00) then
              kp=99.00
              kp3=99.00
          else
              kp=kp1*f1+kp2*f2
              kp3=kp31*f1+kp32*f2
          endif
 
          if(dst1.eq.99999 .or. dst2.eq.99999) then
              dst=99999
          else
              dst=nint(dst1*f1+dst2*f2)
          endif

          ih=ihr1+i/60
          idoy=idoy1
          if(ih.eq.24) then
             ih=ih-24
             idoy=idoy1+1
          endif
          j=i-i/60*60
          iyr3=iyr1
          if(mod(iyr1,4).eq.0 .and.idoy.gt.366) then
              idoy=1
              iyr3=iyr1+1
          endif
          if(mod(iyr1,4).ne.0 .and.idoy.gt.365) then
              idoy=1
              iyr3=iyr1+1
          endif
          
          write(2,444)iyr3,idoy,ih,j,kp,kp3,dst
      enddo
      
      iyr1=iyr2
      idoy1=idoy2
      ihr1=ihr2
      kp1=kp2
      kp31=kp32
      dst1=dst2
      goto 10
 
 20   do i=30,59
          write(2,444)iyr2,idoy2,ihr2,i,kp2,kp32,dst2
      enddo
    
      close(1,status='delete') 
      CLOSE(2) 
1000  write(*,*)'     Subroutine kpdst_min      DONE'
      END   
cccccccccccc end of subroutine kpdst_min  ccccccccccccccccccccccccccccccccccccccc


cccccccccccc subroutine bsgamma_min  ccccccccccccccccccccccccccccccccccccccccccc

c     this program adds bs**gamma columns    ! runs 10mins for 20years/data, and 5mins for 1year
      subroutine bsgamma_min
      dimension bsg(6),gamma(6),idmy(7),dummy(24)
      character finput*30,fout*30
      real N,V,P
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/
      common wst(6)
      common /years/year0,nyrs
    
      finput='1min/omni_min.asc' 
      fout='1min/model_inputr.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34
 33   write(*,*)finput//' NOT EXIST'
      goto 1000
 34   open(unit=2,file=fout,status='unknown') 
      WRITE(2,400) 'YEAR',' DOY','HR','MIN','ByIMF','BzIMF','V_SW',
     $'Den_P','Pdyn','Bs1g','Bs2g','Bs3g','Bs4g','Bs5g','Bs6g'
400   format(a5,a4,2a3,3a8,a7,a6,6a8)	

10    read(1,40,err=10,end=1000) iyr,idoy,ihr,imin,(idmy(i),i=1,7),d,i1 
     $,(dummy(i),i=1,24) 
       IF (iyr.lt.year0.or.iyr.ge.year0+nyrs) go to 10  ! Alin
         By=dummy(5) 
         Bz=dummy(6)
         N=dummy(13)
         V=dummy(9)
         P=dummy(15)  
         if(Bz.eq.9999.99) then
            do i=1,6
               bsg(i)=9999.99
            enddo
         elseif(Bz.lt.0.0) then
            do i=1,6
               bsg(i)=(-Bz)**gamma(i)
            enddo
         else
            do i=1,6
               bsg(i)=0.0
            enddo
         endif
40    format(2I4,4I3,3I4,2I7,F6.2,I7, 8F8.2,4F8.1,F7.2, 
     $F9.0,F6.2,2F7.2,F6.1,6F8.2) 

	write(2,500)iyr,idoy,ihr,imin,By,Bz,V,N,P,(bsg(j),j=1,6)
	goto 10
500   FORMAT(I5,i4,2i3,2F8.2,F8.1,F7.2,F6.2,6F8.2)
      close(1) 
      close(2) 
	 
1000  write(*,*)'     Subroutine bsgamma_min    DONE'
      end
cccccccccccccccc  end of bsgamma_min  cccccccccccccccccccccccccccccccccccccccc


cccccccccccccccc subroutine av1day  cccccccccccccccccccccccccccccccccccccccccc

c    This program is to take the 1 day average of all parameters.

      subroutine av1day_min

      parameter (nipyr=365,nii=30000,nj=11)
      dimension y(nii),p(nii,nj),w(nii,nj),pv(nj)
      dimension idmy(7),dummy(24),dd(11),igap(3),it1(3),it2(3)
      double precision yrdouble
      character finput*30,fout*30
      common wst(6)
      common /years/year0,nyrs
 
      ni=nipyr*nyrs
      iy1= nint(year0)
      dyr= 1./nipyr
      do  i= 1,ni
         y(i)= (i-0.5)*dyr
      enddo
      do 23002 j= 1,nj
         do 23004 i= 1,ni
            p(i,j)= 0.
            w(i,j)= 0.
23004    continue
23002 continue
     
      finput='1min/model_inputr.d'
      fout='1min/av1day.d'
      open(unit=11,file=finput,status='old')
 34   read(11,20,err=34,end=300)iyr4,idoy,ihr,imin,(pv(j),j=1,nj)
 20   format(I5,i4,2i3,2F8.2,F8.1,F7.2,F6.2,6F8.2)     
        
         icheck=icheck+1
         if( mod(iyr4,4).eq.0 ) then
            doytot= 366.
         else
            doytot= 365.
         endif
         iyr4m= iyr4 - iy1

         yv= iyr4m + ( idoy -1. + (ihr + imin/60.0)/24. )/doytot

         rii= (yv-y(1))/dyr + 1.
         ii= rii
         ip= ii + 1
         fip= rii - ii
         fii= 1. - fip

         do 23008 j= 1,nj

         if( (j.eq.1 .or. j.eq.2).and. pv(j).eq.9999.99 ) go to 23008
            if( j.eq.3 .and. pv(j).eq.99999.9 ) go to 23008
            if( j.eq.4 .and. pv(j).eq.999.99 ) go to 23008
            if( j.eq.5 .and. pv(j).eq.99.99 ) go to 23008
            if( j.ge.6 .and. pv(j).eq.9999.99) go to 23008
            if( ii.ge.1 .and. ii.le.ni ) then
               p(ii,j)= p(ii,j) + fii*pv(j)
               w(ii,j)= w(ii,j) + fii
            endif
            if( ip.ge.1 .and. ip.le.ni ) then
               p(ip,j)= p(ip,j) + fip*pv(j)
               w(ip,j)= w(ip,j) + fip
            endif
23008    continue

      go to 34

 300  close(11,status='delete')

      do 23010 j= 1,nj
         do 23012 i= 1,ni
            if( w(i,j).eq. 0. .or. w(i,j).lt.360 ) then
               p(i,j)= 999.
            else
               p(i,j)= p(i,j)/w(i,j)
            endif
23012    continue
23010 continue
      open(unit=21,file=fout,status='unknown')
      do 23014 i= 1,ni
	   yrdouble=y(i)+iy1
          write(21,2100) yrdouble,(p(i,j),j=1,nj)
 2100       format(f9.4,11(1x,f9.4))
23014 continue
      close(21)
      write(*,*)'     Subroutine av1day_min     DONE'
      end
cccccccccccccccc end of av1day cccccccccccccccccccccccccccccccccccccccc


ccccccccccccccc subroutine inter1day cccccccccccccccccccccccccccccccccc

c     This program is to linearly interpolate across the gaps in the data,calculate bz1-bz6

      subroutine inter1day_min

      parameter (nk=12)     
      character finput*30,fout*30
      dimension av1(nk,36500),gamma(6)
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/
      common wst(6)
      common /years/year0,nyrs

      ni=365*nyrs
      finput='1min/av1day.d'
      fout='1min/interav1day.d'
      
cccc  input the 1 day average
      open(unit=1,file=finput,status='old')
      do i=1,ni
      read(1,*,end=35) (av1(j,i),j=1,nk)
      enddo
 35   close(1,status='delete')

      do  j=1,ni
        do  k=2,nk
        if(av1(k,j).eq.999.) then
            j1=j-1
            j2=j+1
            DO WHILE (av1(k,j1) .eq. 999.)
                j1=j1-1
                if(j1.eq.0) goto  10
            END DO
 10         DO WHILE (av1(k,j2) .eq. 999.)
                j2=j2+1
                if(j2.gt.ni) goto 20
            END DO
 20         if(j1.eq.0) then
                av1(k,j)=av1(k,j2)
            elseif(j2.gt.ni) then
                av1(k,j)=av1(k,j1)
            else
                av1(k,j)=(j-j1)*(av1(k,j2)-av1(k,j1))/(j2-j1)+av1(k,j1)
            endif
        endif
        enddo
      enddo
      
      do j=1,ni
         do k=7,nk
            av1(k,j)=-av1(k,j)**(1.0/gamma(k-6))
         enddo
      enddo

      open(unit=2,file=fout,status='unknown')
      do i=1,ni
         write(2,2100) (av1(j,i),j=1,nk)
 2100    format(f9.4,11(1x,f9.4))
      enddo
      close(2) 
        
      write(*,*)'     Subroutine inter1day_min  DONE'
      end
ccccccccccccccc  end of inter1day cccccccccccccccccccccccccccccccccccccccc


ccccccccccccccc subroutine interpl_min cccccccccccccccccccccccccccccc

c     This program is to interpolate across the gaps in the one-minute data 

      subroutine interpl_min

      parameter(nc=11,nd=386,nmax=12999999)            ! nd not used, nmax < 13,100,000
      dimension av1(nc+1,99999),istatus(nc),i1(nc),i2(nc), 
     $ainter(nc),check1(nc),check2(nc),check(nc)
      double precision tyear,datamin(nc+1,nmax)
      dimension idmy(7),dummy(24),correltime(nc),itime(4,nmax)
      character finput2*30,finput1*30,fout*30
      data correltime/10.59,7.91,64.87,35.95,29.41,
     $7.91,7.91,7.91,7.91,7.91,7.91/
      common wst(6)
      common /years/year0,nyrs
      
      finput1='1min/interav1day.d'
      finput2='1min/omni_min.asc'
      fout='1min/intermin.d'
      n=0

cccc  input the 1 day average ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=1,file=finput1,status='old')
      do i=1,99999
      read(1,*,end=34)(av1(j,i),j=1,12)
      enddo
 34   close(1,status='delete')

      nday=i-1

ccccccc  input the original data   cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=2,file=finput2,status='old')

      do j=1,nmax
444   read(2,40,err=20,end=20) iyr4,idoy,ihr,imin,(idmy(i),i=1,7),d,ii
     $,(dummy(i),i=1,24)
      IF (iyr4.lt.year0.or.iyr4.ge.year0+nyrs) go to 444  ! Alin
         n=n+1
         itime(1,n)=iyr4
         itime(2,n)=idoy
         itime(3,n)=ihr
         itime(4,n)=imin
         call year(iyr4,idoy,ihr,imin,tyear)
         datamin(1,n)=tyear
         datamin(2,n)=dummy(5)  !by
         datamin(3,n)=dummy(6)  !bz
         datamin(4,n)=dummy(9)  !v
         datamin(5,n)=dummy(13) !den
         datamin(6,n)=dummy(15) !pdyn
         datamin(7,n)=dummy(6)
         datamin(8,n)=dummy(6)
         datamin(9,n)=dummy(6)
         datamin(10,n)=dummy(6)
         datamin(11,n)=dummy(6)
         datamin(12,n)=dummy(6)
      enddo
 20   close(2)
40    format(2I4,4I3,3I4,2I7,F6.2,I7, 8F8.2,4F8.1,F7.2,
     $F9.0,F6.2,2F7.2,F6.1,6F8.2)

      ni=n
      
ccccccccccccccccccccc If gap is at the beginning or the end ccccccccccccccc
 43   do j=1,nc
          check1(j)=0
          check2(j)=0
          check(j)=0
      enddo
      do k=2,12
         n=k-1
         if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. datamin(k,1).eq.9999.99)
     &.or. ( k.eq.4 .and. datamin(k,1).eq.99999.9 ) .or. ( k
     &.eq.5 .and. datamin(k,1).eq.999.99 ) .or. ( k.eq.6 .and. 
     $datamin(k,1).eq.99.99) ) then      !if there is a gap at the beginning
                datamin(k,1)=av1(k,1)
                check1(n)=1
         endif
         if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. datamin(k,ni).eq.9999.99)
     &.or. ( k.eq.4 .and. datamin(k,ni).eq.99999.9 ) .or. ( k
     &.eq.5 .and. datamin(k,ni).eq.999.99 ) .or. ( k.eq.6 .and. 
     $datamin(k,ni).eq.99.99) ) then      !if there is a gap at the end
                datamin(k,ni)=av1(k,nday)
                check2(n)=1
         endif
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      open(unit=3,file=fout,status='unknown')
    
      write(3,*) 's=2-->original           s=1-->interpolation 
     $s=0-->1 day average'
      write(3,52) 'Year','Day','Hr','Min','ByIMF','s','BzIMF','s',
     $'V_SW','s','Den_P','s','Pdyn','s','Bz1','s','Bz2','s',
     $'Bz3','s','Bz4','s','Bz5','s','Bz6','s'
52    format(a5,a4,2a3,2(a8,a2),a8,a2,a7,a2,a6,a2,6(a8,a2))
      i00=0

      DO 2002 j=1,ni

cw       if (itime(1,j).gt.i00.and.nyrs.gt.1)
cw     &   write (*,*) ' year',itime(1,j),' started'
cw       i00=itime(1,j)

       do 4000 k=2,12

        n=k-1
        if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. datamin(k,j).ne.9999.99)
     &.or. ( k.eq.4 .and. datamin(k,j).ne.99999.9 ) .or. ( k
     &.eq.5 .and. datamin(k,j).ne.999.99 ) .or. ( k.eq.6 .and. 
     $datamin(k,j).ne.99.99) ) then
            istatus(n)=2
            if(check1(n).eq.1.and.j.eq.1) then
               istatus(n)=0
            endif

            if(check2(n).eq.1.and.j.eq.ni) then
               istatus(n)=0
            endif
            ainter(n)=datamin(k,j)
            check(n)=0
       else

         if(check(n).eq.0)  then !when check(n).eq.0, it means a new gap coming up
            check(n)=1       
            i1(n)=j-1
            i2(n)=j+1
                      
            if( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) then
       	      DO WHILE (datamin(k,i1(n)) .eq. 9999.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datamin(k,i2(n)) .eq. 9999.99)
                 i2(n)=i2(n)+1
              END DO
            endif
         
            if( k.eq.4 ) then
              DO WHILE (datamin(k,i1(n)) .eq. 99999.9)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datamin(k,i2(n)) .eq. 99999.9)
                 i2(n)=i2(n)+1
              END DO
            endif
     
            if( k.eq.5 ) then
              DO WHILE (datamin(k,i1(n)) .eq. 999.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datamin(k,i2(n)) .eq. 999.99)
                 i2(n)=i2(n)+1
              END DO
            endif
    
            if( k.eq.6 ) then
              DO WHILE (datamin(k,i1(n)) .eq. 99.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datamin(k,i2(n)) .eq. 99.99)
                 i2(n)=i2(n)+1
              END DO
            endif
	 endif  
                    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    	   
            timegap=(i2(n)-i1(n)-1)/correltime(n)
           if(timegap .le. 2.0)	then
             ainter(n)=(j-i1(n))*(datamin(k,i2(n))-datamin(k,i1(n)))/
     $(i2(n)-i1(n))+datamin(k,i1(n))
             istatus(n)=1
            endif

           if((timegap .gt. 2.0) .and. (timegap .le. 4.0)) then
              amidtime=(i2(n)+i1(n))/2.0
              if(abs(j-amidtime) .le. correltime(n))  then
                datamin_mid=(datamin(k,i1(n))+datamin(k,i2(n)))/2.0
                ainter(n)=datamin_mid+(j-amidtime)*
     $(datamin(k,i2(n))-datamin(k,i1(n)))/(2*correltime(n))
              else if((j-i1(n)) .lt. correltime(n)) then
                ainter(n)=datamin(k,i1(n))
              else
	        ainter(n)=datamin(k,i2(n))
              endif
              istatus(n)=1
           endif

           if(timegap .gt. 4.0)	then
              ii=int(correltime(n))+1
              if((j-i1(n)) .lt. ii) ainter(n)=datamin(k,i1(n))
       	      if((i2(n)-j) .lt. ii) ainter(n)=datamin(k,i2(n))
              istatus(n)=1
       	      if((j-i1(n)) .ge. ii .and. (i2(n)-j) .ge. ii) then
               w1=max(1-(j-i1(n)-correltime(n))/(2.0*correltime(n)),0.0)
               w2=max(1-(i2(n)-correltime(n)-j)/(2.0*correltime(n)),0.0)
    	         wav=1-w1-w2
                 fav=0.0
	         if(datamin(1,j).lt. av1(1,1)) fav=av1(k,1)
                 if(datamin(1,j).gt. av1(1,nday)) fav=av1(k,nday)
                 do  i3=1,nday-1
                    if( datamin(1,j).ge.av1(1,i3) .and.datamin(1,j).le.
     $av1(1,i3+1))  then
                     fav=(datamin(1,j)-av1(1,i3))*(av1(k,i3+1)
     $-av1(k,i3))/(av1(1,i3+1)-av1(1,i3))+av1(k,i3)
                     goto 35
                    endif
                 enddo	
                                 
 35            ainter(n)=w1*datamin(k,i1(n))+w2*datamin(k,i2(n))+wav*fav
                  if(ainter(n) .eq. fav)  istatus(n)=0
               endif
         endif
        
         if(check1(n).eq.1.and.i1(n).eq.1
     $.and.(i2(n)-j).ge.3*correltime(n)) then
             istatus(n)=0
         endif

         if(check2(n).eq.1.and.i2(n).eq.ni
     $.and.(j-i1(n)).ge.3*correltime(n)) then
             istatus(n)=0
         endif
	          
      endif
                     
4000  continue
      
      write(3,50)itime(1,j),itime(2,j),itime(3,j),itime(4,j),
     $ainter(1),istatus(1),ainter(2),istatus(2),ainter(3),istatus(3),
     $ainter(4),istatus(4),ainter(5),istatus(5),ainter(6),istatus(6),
     $ainter(7),istatus(7),ainter(8),istatus(8),ainter(9),istatus(9),
     $ainter(10),istatus(10),ainter(11),istatus(11)
50    format(i5,i4,2i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,6(f8.2,i2))

2002  Continue

      close(3)
      write(*,*)'     Subroutine interpl_min    DONE'
      end

cccccccccccccccccc end of interpl_min ccccccccccccccccccccccccccccccc      


      subroutine year(iyr4,idoy,ihr,imin,year12)
      double precision year12,dum
      if (mod(iyr4,4).eq.0) then
         dum=(idoy-1+(ihr+imin/60.)/24.)/366.0
         year12=iyr4*1.0+dum
      else
         dum=(idoy-1+(ihr+imin/60.)/24.)/365.0
         year12=iyr4*1.0+dum
      endif
      end


cccccccccccccccccc subroutine WG_min  ccccccccccccccccccccccccccccccccccccc

c    This program is to calculate the W and G parameters
      subroutine WG_min
      common wst(6)
      real kp,kp3
      integer dst
      dimension dat(4,60),idat(4,60),W(6),iW(6),bzn(6),S(6),bsn(6)
      dimension r(6),Wf(6),gamma(6),beta(6),xlumda(6),f(6)
      character str1*8,str2*12,aa      
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/    
      data beta/.8,.18,2.32,1.25,1.6,2.4/
      data xlumda/.39,.46,.39,.42,.41,1.29/
      data r/.39,.7,.031,.58,1.15,.88/
      data f/1.045,1.078,1.017,1.060,1.081,1.280/
      character finput1*30,finput2*30,fout*30,foutput*28
      character datechar*10,timechar*8,ychar*4,mchar*2,dchar*2
            
      do i=1,6
         W(i)=wst(i)/f(i)
         Wf(i)=0.0
      enddo
              
      pid4 = atan2(1.,1.)
      pi2 = 8*pid4
      n=0
      
      finput1='1min/intermin.d'
      finput2='1min/kpkp3dstmin.d'
      fout='1min/WGparametersmin.d'  
          
      open(unit=1,file=finput1,status='old')
      open(unit=2,file=finput2,status='old')
c      open(unit=3,file=fout,status='unknown')
c      write(3,52) 'Year','Day','Hr','Min','ByIMF','BzIMF','V_SW',
c     $'Den_P','Pdyn','G1','G2','G3','8 status','kp','akp3','dst'
c     $,'Bz1','Bz2','Bz3','Bz4','Bz5','Bz6','W1'
c     $,'W2','W3','W4','W5','W6','6 stat'
c52    format(a5,a4,a3,a4,a6,7a7,a10,2a7,a6,6a7,6a9,a8)
     
      read(1,*)
      read(1,*)
      read(2,*)
      idayold=0
       
 10   read(2,*,err=20,end=20)iyr1,idy1,ihr1,imn1,kp,kp3,dst
      read(1,70,err=20,end=20)iyr,idy,ihr,imn,by,iby,bz,ibz,v,iv,den,id,
     $p,ip,bz1,i1,bz2,i2,bz3,i3,bz4,i4,bz5,i5,bz6,i6
       if (iyr1.ne.iyr) write (*,*) 'year mismatch in WGmin input'
 70   format(i5,i4,2i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,6(f8.2,i2))

      call daymonth(iyr,idy,jmo,iday)
      call utc(iyr,jmo,iday,ihr,imn,datechar,timechar,mchar,dchar)
      if (idy.ne.idayold) then
       if (idayold.ne.0) close(3)
       idayold=idy
       write (ychar,'(i4)') iyr
       foutput=ychar//'/QinDenton_'//ychar//mchar//dchar//'_1min'
       open(unit=3,file=foutput,status='unknown')
      end if

      n=n+1

cccccccc  G parameters cccccccccccccccccccccccccccccccccccc
      if(n.le.60) then
           dat(1,n)=by
           dat(2,n)=bz
           dat(3,n)=v
           dat(4,n)=den
           idat(1,n)=iby
           idat(2,n)=ibz
           idat(3,n)=iv
           idat(4,n)=id
      else
           do i=2,60
              dat(1,i-1)=dat(1,i)
              dat(2,i-1)=dat(2,i)
              dat(3,i-1)=dat(3,i)
              dat(4,i-1)=dat(4,i)
              idat(1,i-1)=idat(1,i)
              idat(2,i-1)=idat(2,i)
              idat(3,i-1)=idat(3,i)
              idat(4,i-1)=idat(4,i)
           enddo
           dat(1,60)=by
           dat(2,60)=bz
           dat(3,60)=v
           dat(4,60)=den
           idat(1,60)=iby
           idat(2,60)=ibz
           idat(3,60)=iv
           idat(4,60)=id
      endif
      ncap=min(n,60)
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      iG1=0
      iG2=0
      iG3=0
      
      do i=1,ncap
           bperp = sqrt(dat(1,i)**2+dat(2,i)**2)
           bp40 = bperp/40.
           HH = bp40**2/(1.+bp40)

           if(dat(1,i).eq.0..and.dat(2,i).eq.0.) then
                  theta = 0.
           else
                  theta = atan2(dat(1,i),dat(2,i))
                  if (theta.lt.0) theta = theta + pi2
           endif
           
           if (dat(2,i).lt.0.) then
                  Bs = -dat(2,i)
           else
                  Bs = 0.
           endif
               
           sine = (sin(theta/2.))**3
           sum1 = sum1 + dat(3,i)*HH*sine
           iG1=iG1+idat(1,i)+idat(2,i)+idat(3,i)
           sum2 = sum2 + Bs*dat(3,i)
           iG2=iG2+idat(2,i)+idat(3,i)
           sum3 = sum3 + Bs*dat(3,i)*dat(4,i)
           iG3=iG3+idat(2,i)+idat(3,i)+idat(4,i)
      enddo
      G1 = sum1/(ncap*1.0)
      G2 = sum2/(ncap*1.0)/200.
      G3 = sum3/(ncap*1.0)/2000.
      iG1=nint(iG1/(3.0*ncap))
      iG2=nint(iG2/(2.0*ncap))
      iG3=nint(iG3/(3.0*ncap))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccc W parameters ccccccccccccccccccccccccccccccccccccccccccccc
      bzn(1)=bz1
      bzn(2)=bz2
      bzn(3)=bz3
      bzn(4)=bz4
      bzn(5)=bz5
      bzn(6)=bz6
      do i=1,6
         if(bzn(i).lt.0.0) then
             bsn(i)=-bzn(i)                 
         else
             bsn(i)=0.0
         endif

         S(i)=(den/5)**xlumda(i)*(v/400)**beta(i)*(bsn(i)/5)**gamma(i)
         W(i)=r(i)/60.0*S(i)+W(i)*exp(-r(i)/60.0)
         nf=ibz*iv*id
         
         if(nf.eq.0) then
               Wf(i)=r(i)/60.0*S(i)*0.0+Wf(i)*exp(-r(i)/60.0)	
         elseif(nf.eq.8) then
               Wf(i)=r(i)/60.0*S(i)*2.0+Wf(i)*exp(-r(i)/60.0)
         else
               Wf(i)=r(i)/60.0*S(i)*1.0+Wf(i)*exp(-r(i)/60.0)
         endif
     
         if(W(i).ne.0) then          
               iW(i)=nint(Wf(i)/W(i))
         else
               iW(i)=nint((ibz+iv+id)/3.0)
         endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     

c      if(idy1.eq.idy.and.ihr1.eq.ihr.and.imn1.eq.imn) then
c         write(3,50)iyr,idy,ihr,imn,by,bz,v,den,p,G1,G2,G3,
c     $iby,ibz,iv,id,ip,iG1,iG2,iG3,kp,kp3,dst,(bzn(i),i=1,6),
c     $(W(i)*f(i),i=1,6),(iW(i),i=1,6)
c50    format(I5,I4,2I3,2f7.2,f7.1,2f7.2,3f7.2,i3,7i1,
c     $f7.2,f7.2,i6,6f7.2,6f9.3,i3,5i1)
c      endif
      
      If (idy1.eq.idy.and.ihr1.eq.ihr.and.imn1.eq.imn) then
c       if (dst.ne.99999) 
       write (3,500)
     &  datechar,timechar,iyr,jmo,iday,ihr,imn,
     &  by,bz,v,den,p,G1,G2,G3,
     &  iby,ibz,iv,id,ip,iG1,iG2,iG3,
     &  kp,kp3,dst,(bzn(i),i=1,6),(W(i)*f(i),i=1,6),(iW(i),i=1,6)
      End If
500    format(a10,'T',a8,i5,x,4i3,' 00',x,
     &  2f7.2,f7.1,x,2f7.2,x,3f7.2,x,
     &  8i2,x,
     &  2f6.2,i6,x,6f7.2,x,6f8.3,x,6i2)

      Go To 10

20    close(1,status='delete')
      close(2,status='delete')
      close(3)
      write(*,*)'     Subroutine WG_min         DONE'	     

      end
cccccccccccccccccc end of WG_min cccccccccccccccccccccccccccccccccccccccccccc


c  **********  FOR 5 MIN DATA ********************


cccccc subroutine kp3  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c This program calculates the 3 day average of Kp: Kp3.
     
      subroutine  kp3_5min
      integer dst,cap
      character finput*30,fout*30
      real Kp,akp3,Kpa(120),akp,weight,t3day 
      common wst(6)
      common /years/year0,nyrs
      
      finput='5min/kpdst.lst'              ! Kp, Dst data from 1995 to 2015
      fout='5min/kpkp3dst.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34
 33   write(*,*)finput//' NOT EXIST'
      goto 1000
 34   open(unit=2,file=fout,status='unknown') 
 
      WRITE(2,*)'YEAR DOY HR   Kp     akp3  DST' 
 
      do I=1,120 
         Kpa(I)=0. 
      enddo 
      J=0      

10    read(1,*,err=20,end=20) iyr4,idoy,ihr,Kp,dst
       IF (iyr4.lt.year0) go to 10          ! Alin: load data starting at year0
      J=J+1
      if(Kp.eq.99) then
         akp3=99.
      else
         Kp = Kp/10. 
         if (J.lt.120) then 
            cap = J 
         else 
            cap = 120 
         endif 
         do K=1,119 
            Kpa(K)=Kpa(K+1) 
         enddo 
         Kpa(120)=Kp 
C    Calculate akp3 
         
         akp = 0 
         weight = 0 
         t3day = 72         !Time for 3 days in hours 

         do ii=1,cap 
            weight = weight + exp(-((ii-1)/t3day)) 
            akp = akp + exp(-((ii-1)/t3day))*Kpa(121-ii) 
         enddo
 
         akp3= akp/weight 
      endif
   
      WRITE(2,444)iyr4,idoy,ihr,Kp,akp3,dst
 444  FORMAT(I4,1x,I3,1x,I2,1x,2(F6.2,1x),i6) 
      goto 10
 
 20   close(1) 
      CLOSE(2) 
1000  write(*,*)'     Subroutine kp3_5min       DONE'
      END  
cccccccccccccc end of subroutine kp3_5min  ccccccccccccccccccccccccccccccc


cccccccccccccccc kpdst_5min cccccccccccccccccccccccccccccccccccccccccccccc

c    This program interpolate the houly kp and dst to 5-min data
      subroutine kpdst_5min
      common wst(6)
      character year*4,finput*30,yeard*5,fout*30
      real kp1,kp31,kp2,kp32,kp,kp3
      integer dst1,dst2,dst 
 
      finput='5min/kpkp3dst.d' 
      fout='5min/kpkp3dst5min.d' 
 
      open(unit=1,file=finput,status='old') 
      open(unit=2,file=fout,status='unknown') 
 
      WRITE(2,*)'YEAR DOY HR MIN  Kp     kp3     DST' 
      read(1,*)
      read(1,*,err=20,end=20) iyr1,idoy1,ihr1,kp1,kp31,dst1
      do i=0,25,5
          write(2,444)iyr1,idoy1,ihr1,i,kp1,kp31,dst1
444   FORMAT(I4,1x,I3,1x,I2,1x,i2,1x,2(F6.2,1x),i6) 
      enddo          

10    read(1,*,err=20,end=20) iyr2,idoy2,ihr2,kp2,kp32,dst2
      do i=30,85,5
          f2=(i*1.0-30.0)/60.0
          f1=1-f2

          if(kp1.eq.99.00 .or. kp2.eq.99.00) then
              kp=99.00
              kp3=99.00
          else
              kp=kp1*f1+kp2*f2
              kp3=kp31*f1+kp32*f2
          endif
 
          if(dst1.eq.99999 .or. dst2.eq.99999) then
              dst=99999
          else
              dst=nint(dst1*f1+dst2*f2)
          endif

          ih=ihr1+i/60
          idoy=idoy1
          if(ih.eq.24) then
             ih=ih-24
             idoy=idoy1+1
          endif
          j=i-i/60*60
          iyr3=iyr1
          if(mod(iyr1,4).eq.0 .and.idoy.gt.366) then
              idoy=1
              iyr3=iyr1+1
          endif
          if(mod(iyr1,4).ne.0 .and.idoy.gt.365) then
              idoy=1
              iyr3=iyr1+1
          endif
          write(2,444)iyr3,idoy,ih,j,kp,kp3,dst
      enddo
      
      iyr1=iyr2
      idoy1=idoy2
      ihr1=ihr2
      kp1=kp2
      kp31=kp32
      dst1=dst2
      goto 10
 
 20   do i=30,55,5
          write(2,444)iyr2,idoy2,ihr2,i,kp2,kp32,dst2
      enddo
    
      close(1,status='delete') 
      CLOSE(2) 
1000  write(*,*)'     Subroutine kpdst_5min     DONE'
      END 
cccccccccccccccc end of kpdst_5min ccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccc bsgamma_5min cccccccccccccccccccccccccccccccccccccccccccccccc

c    This program adds the bs**gamma columns
	subroutine bsgamma_5min
        common wst(6)
        dimension bsg(6),gamma(6),idmy(7),dummy(24)
        character finput*30,fout*30,aa
        real N,V,P
	data gamma/0.87,0.67,1.32,1.29,0.69,0.53/

      
      finput='5min/omni_5min.asc'              ! OMNI data from 1981 to 2015
      fout='5min/model_inputr5min.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34
 33   write(*,*)finput//' NOT EXIST'
      goto 1000
 34   open(unit=2,file=fout,status='unknown') 
      WRITE(2,400) 'YEAR',' DOY','HR','MIN','ByIMF','BzIMF','V_SW',
     $'Den_P','Pdyn','Bs1g','Bs2g','Bs3g','Bs4g','Bs5g','Bs6g'
400   format(a5,a4,2a3,3a8,a7,a6,6a8)	
10    read(1,40,err=1000,end=1000) iyr,idoy,ihr,imin,aa,By,Bz,aa,
     $V,aa,N,aa,P,aa
40    format(2I4,2I3,a77,2f8.2,a16,f8.1,a24,f7.2,a9,f6.2,a68) 
         
         if(Bz.eq.9999.99) then
            do i=1,6
               bsg(i)=9999.99
            enddo
         elseif(Bz.lt.0.0) then
            do i=1,6
               bsg(i)=(-Bz)**gamma(i)
            enddo
         else
            do i=1,6
               bsg(i)=0.0
            enddo
         endif

	write(2,500)iyr,idoy,ihr,imin,By,Bz,V,N,P,(bsg(j),j=1,6)
	goto 10
500   FORMAT(I5,i4,2i3,2F8.2,F8.1,F7.2,F6.2,6F8.2)

1000  continue	  
      close(1) 
      close(2) 
      write(*,*)'     Subroutine bsgamma_5min   DONE'
      end
ccccccccccccccc end of bsgamma_5min   ccccccccccccccccccccccccccccc


ccccccccccccccc av1day_5min.f ccccccccccccccccccccccccccccccccccccc

c     This program takes the 1 day average of the data.

      subroutine av1day_5min

      parameter (nipyr=365,nii=30000,nj=11)      
      dimension y(nii),p(nii,nj),w(nii,nj),pv(nj)
      dimension idmy(7),dummy(24),dd(11),igap(3),it1(3),it2(3)
      double precision yrdouble
      character finput*30,fout*30           
      common /years/year0,nyrs
      common wst(6)
 
      ni=nipyr*nyrs
      iy1= nint(year0)              ! was 1995, the 1min data over more than 20 years exceeds RAM 
      dyr= 1./nipyr                   !    also kpdst.lst starts 1995
      do  i= 1,ni
         y(i)= (i-0.5)*dyr
      enddo
      do 23002 j= 1,nj
         do 23004 i= 1,ni
            p(i,j)= 0.
            w(i,j)= 0.
23004    continue
23002 continue
     
      finput='5min/model_inputr5min.d'
      fout='5min/av1day5min.d'

      open(unit=11,file=finput,status='old')

 34   read(11,20,err=34,end=300)iyr4,idoy,ihr,imin,(pv(j),j=1,nj)
 20   format(I5,i4,2i3,2F8.2,F8.1,F7.2,F6.2,6F8.2)     
       
      if( mod(iyr4,4).eq.0 ) then
         doytot= 366.
      else
         doytot= 365.
      endif

      iyr4m= iyr4 - iy1                                        ! can be negative
      yv= iyr4m + ( idoy -1. + (ihr + imin/60.0)/24. )/doytot
      rii= (yv-y(1))/dyr + 1.                                  ! can be negative
      ii= rii                                                  ! integer(rii)
      ip= ii + 1           
      fip= rii - ii
      fii= 1. - fip

      do 23008 j= 1,nj
        if( (j.eq.1 .or. j.eq.2).and. pv(j).eq.9999.99 ) go to 23008
          if( j.eq.3 .and. pv(j).eq.99999.9 ) go to 23008
          if( j.eq.4 .and. pv(j).eq.999.99 ) go to 23008
          if( j.eq.5 .and. pv(j).eq.99.99 ) go to 23008
          if( j.ge.6 .and. pv(j).eq.9999.99) go to 23008
          if( ii.ge.1 .and. ii.le.ni ) then     ! retains data only between year0 and year0+nyrs
             p(ii,j)= p(ii,j) + fii*pv(j)
             w(ii,j)= w(ii,j) + fii
          endif
          if( ip.ge.1 .and. ip.le.ni ) then
             p(ip,j)= p(ip,j) + fip*pv(j)
             w(ip,j)= w(ip,j) + fip
          endif
23008 continue
      go to 34

300   close(11,status='delete')

      do 23010 j= 1,nj
         do 23012 i= 1,ni
            if( w(i,j).eq. 0. .or.w(i,j).lt.72 ) then
               p(i,j)= 999.
            else
               p(i,j)= p(i,j)/w(i,j)
            endif
23012    continue
23010 continue

      open(unit=21,file=fout,status='unknown')
      do 23014 i= 1,ni                                 ! output at year0-(year0+nyrs)
         yrdouble=y(i)+iy1
         write(21,2100) yrdouble,(p(i,j),j=1,nj)
 2100    format(f9.4,11(1x,f9.4))
23014 continue
      close(21)

      write(*,*)'     Subroutine av1day_5min    DONE'

      end
ccccccccccccccccc end of av1day_5min cccccccccccccccccccccccccccccccccccccc


ccccccccccccccccc inter1day_5min cccccccccccccccccccccccccccccccccccccccccc

c    This program linearly interpolates across the gaps in the file,calculate bz1-bz6 

      subroutine inter1day_5min

      parameter (nk=12)                   
      character finput*30,fout*30,finput1*30,finput2*30
      dimension av1(nk,36500),gamma(6)
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/
      common /years/year0,nyrs
      common wst(6)
     
      ni=365*nyrs
      finput='5min/av1day5min.d'
      fout='5min/interav1day5min.d'
    
ccccccc  input the 1 day average
      open(unit=1,file=finput,status='old')
      do i=1,ni
       read(1,*,end=35) (av1(j,i),j=1,nk)        ! loads data at year0-(year0+nyrs)
      enddo
 35   close(1,status='delete')
   
      do  j=1,ni
        do  k=2,nk
        if(av1(k,j).eq.999.) then
            j1=j-1
            j2=j+1
            DO WHILE (av1(k,j1) .eq. 999.)
                j1=j1-1
                if(j1.eq.0) goto  10
            END DO
 10         DO WHILE (av1(k,j2) .eq. 999.)
                j2=j2+1
                if(j2.gt.ni) goto 20
            END DO
 20         if(j1.eq.0) then
                av1(k,j)=av1(k,j2)
            elseif(j2.gt.ni) then
                av1(k,j)=av1(k,j1)
            else
                av1(k,j)=(j-j1)*(av1(k,j2)-av1(k,j1))/(j2-j1)+av1(k,j1)
            endif
        endif
        enddo
      enddo
      
      do j=1,ni
         do k=7,nk
            av1(k,j)=-av1(k,j)**(1.0/gamma(k-6)) 
         enddo
      enddo

      open(unit=2,file=fout,status='unknown')
      do i=1,ni
         write(2,2100) (av1(j,i),j=1,nk)     ! output at year0-(year0+nyrs)
 2100    format(f9.4,11(1x,f9.4))
      enddo
      close(2) 
        
      write(*,*)'     Subroutine inter1day_5min DONE'

      end
ccccccccccccccccccc end of inter1day_5min ccccccccccccccccccccccccccccc


ccccccccccccccccccc interpl_5min       cccccccccccccccccccccccccccccccc

c    This program  interpolates across the gaps in the 5-minute data     

      subroutine interpl_5min

      parameter(nc=11,nmax=5000000)                       ! 105192*nyrs < nmax < 8e6
      dimension av1(nc+1,99999),istatus(nc),i1(nc),i2(nc), 
     $ainter(nc),check1(nc),check2(nc),check(nc)
      double precision tyear,data5min(nc+1,nmax) 
      dimension idmy(7),dummy(24),correltime(nc),itime(4,nmax)
      character finput2*30,finput1*30,fout*30,aa
      data correltime/2.76,2.17,14.55,8.37,6.94,
     $2.17,2.17,2.17,2.17,2.17,2.17/
      common wst(6)
      common /years/year0,nyrs

      finput1='5min/interav1day5min.d'   
      finput2='5min/omni_5min.asc'
      fout='5min/inter5min.d'
      n=0

cccccccc  input the 1 day average ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=1,file=finput1,status='old')
      do i=1,999999
      read(1,*,end=34)(av1(j,i),j=1,12)    ! loads data at year0-(year0+nyrs)
      enddo
 34   close(1,status='delete')
      nday=i-1
      
ccccccc  input the original data   cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      open(unit=2,file=finput2,status='old')

      do i=1,nmax
555   read(2,40,err=20,end=20) iyr4,idoy,ihr,imin,aa,By,Bz,aa,
     $V,aa,den,aa,P,aa
       IF (iyr4.lt.year0.or.iyr4.ge.year0+nyrs) go to 555  ! Alin: load data  at year0-(year0+nyrs)
         n=n+1
         itime(1,n)=iyr4
         itime(2,n)=idoy
         itime(3,n)=ihr
         itime(4,n)=imin
         call year(iyr4,idoy,ihr,imin,tyear)
         data5min(1,n)=tyear
         data5min(2,n)=By  !by
         data5min(3,n)=Bz  !bz
         data5min(4,n)=V   !v 
         data5min(5,n)=den !den
         data5min(6,n)=P   !pdyn
         data5min(7,n)=Bz
         data5min(8,n)=Bz
         data5min(9,n)=Bz
         data5min(10,n)=Bz
         data5min(11,n)=Bz
         data5min(12,n)=Bz
      enddo
 20   close(2)
40    format(2I4,2I3,a77,2f8.2,a16,f8.1,a24,f7.2,a9,f6.2,a68)
      ni=n
     
ccccccccccccccccccccc If gap is at the beginning or the end ccccccccccccccc
 43   do n=1,nc
          check1(n)=0
          check2(n)=0
          check(n)=0
      enddo
      do k=2,12
         n=k-1
         if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. data5min(k,1).eq.9999.99)
     &.or. ( k.eq.4 .and. data5min(k,1).eq.99999.9 ) .or. ( k
     &.eq.5 .and. data5min(k,1).eq.999.99 ) .or. ( k.eq.6 .and. 
     $data5min(k,1).eq.99.99) ) then      !if there is a gap at the beginning
                data5min(k,1)=av1(k,1)
                check1(n)=1
         endif
         if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. data5min(k,ni).eq.9999.99)
     &.or. ( k.eq.4 .and. data5min(k,ni).eq.99999.9 ) .or. ( k
     &.eq.5 .and. data5min(k,ni).eq.999.99 ) .or. ( k.eq.6 .and. 
     $data5min(k,ni).eq.99.99) ) then      !if there is a gap at the end
                data5min(k,ni)=av1(k,nday)
                check2(n)=1
         endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      open(unit=3,file=fout,status='unknown')
      write(3,*) 's=2-->original           s=1-->interpolation 
     $s=0-->1 day average'
      write(3,52) 'Year','Day','Hr','Min','ByIMF','s','BzIMF','s',
     $'V_SW','s','Den_P','s','Pdyn','s','Bz1','s','Bz2','s',
     $'Bz3','s','Bz4','s','Bz5','s','Bz6','s'
52    format(a5,a4,2a3,2(a8,a2),a8,a2,a7,a2,a6,a2,6(a8,a2))
      i00=0

      DO 2002 j=1,ni

cw       if (itime(1,j).gt.i00.and.nyrs.gt.1)
cw     &    write (*,*) ' year',itime(1,j),' started'
       i00=itime(1,j)

       do 4000 k=2,12

        n=k-1

        if( ( ( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) 
     $.and. data5min(k,j).ne.9999.99)
     &.or. ( k.eq.4 .and. data5min(k,j).ne.99999.9 ) .or. ( k
     &.eq.5 .and. data5min(k,j).ne.999.99 ) .or. ( k.eq.6 .and. 
     $data5min(k,j).ne.99.99) ) then
            istatus(n)=2
            if(check1(n).eq.1.and.j.eq.1) then
               istatus(n)=0
            endif

            if(check2(n).eq.1.and.j.eq.ni) then
               istatus(n)=0
            endif
            ainter(n)=data5min(k,j)
            check(n)=0              
       else
         if(check(n).eq.0)  then         !when check(n).eq.0, it means a new gap coming up
            check(n)=1       
            i1(n)=j-1
            i2(n)=j+1
            if( k.eq.2 .or. k.eq.3 .or. k.ge.7 ) then
              DO WHILE (data5min(k,i1(n)) .eq. 9999.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (data5min(k,i2(n)) .eq. 9999.99)
                 i2(n)=i2(n)+1
              END DO
            endif
            if( k.eq.4 ) then
              DO WHILE (data5min(k,i1(n)) .eq. 99999.9)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (data5min(k,i2(n)) .eq. 99999.9)
                 i2(n)=i2(n)+1
              END DO
            endif
            if( k.eq.5 ) then
              DO WHILE (data5min(k,i1(n)) .eq. 999.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (data5min(k,i2(n)) .eq. 999.99)
                 i2(n)=i2(n)+1
              END DO
            endif
            if( k.eq.6 ) then
              DO WHILE (data5min(k,i1(n)) .eq. 99.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (data5min(k,i2(n)) .eq. 99.99)
                 i2(n)=i2(n)+1
              END DO
            endif
        endif  

          timegap=(i2(n)-i1(n)-1)/correltime(n)
           if(timegap .le. 2.0)	then
             ainter(n)=(j-i1(n))*(data5min(k,i2(n))-data5min(k,i1(n)))/
     $(i2(n)-i1(n))+data5min(k,i1(n))
             istatus(n)=1
            endif

           if((timegap .gt. 2.0) .and. (timegap .le. 4.0)) then
              amidtime=(i2(n)+i1(n))/2.0
              if(abs(j-amidtime) .le. correltime(n))  then
                data5min_mid=(data5min(k,i1(n))+data5min(k,i2(n)))/2.0
                ainter(n)=data5min_mid+(j-amidtime)*
     $(data5min(k,i2(n))-data5min(k,i1(n)))/(2*correltime(n))
              else if((j-i1(n)) .lt. correltime(n)) then
                ainter(n)=data5min(k,i1(n))
              else
                ainter(n)=data5min(k,i2(n))
              endif
              istatus(n)=1
           endif

           if(timegap .gt. 4.0)	then
              ii=int(correltime(n))+1
              if((j-i1(n)) .lt. ii) ainter(n)=data5min(k,i1(n))
              if((i2(n)-j) .lt. ii) ainter(n)=data5min(k,i2(n))
              istatus(n)=1
              if((j-i1(n)) .ge. ii .and. (i2(n)-j) .ge. ii) then
               w1=max(1-(j-i1(n)-correltime(n))/(2.0*correltime(n)),0.0)
               w2=max(1-(i2(n)-correltime(n)-j)/(2.0*correltime(n)),0.0)
                 wav=1-w1-w2
                 fav=0.0
                 if(data5min(1,j).lt. av1(1,1)) fav=av1(k,1)
                 if(data5min(1,j).gt. av1(1,nday)) fav=av1(k,nday)
                 do  i3=1,99999
                   if( data5min(1,j).ge.av1(1,i3) .and.data5min(1,j).le.
     $av1(1,i3+1))  then
                     fav=(data5min(1,j)-av1(1,i3))*(av1(k,i3+1)
     $-av1(k,i3))/(av1(1,i3+1)-av1(1,i3))+av1(k,i3)
                     goto 35
                    endif
                 enddo                 
 35          ainter(n)=w1*data5min(k,i1(n))+w2*data5min(k,i2(n))+wav*fav
                  if(ainter(n) .eq. fav)  istatus(n)=0
               endif
         endif
        
         if(check1(n).eq.1.and.i1(n).eq.1
     $.and.(i2(n)-j).ge.3*correltime(n)) then
             istatus(n)=0
         endif

         if(check2(n).eq.1.and.i2(n).eq.nm
     $.and.(j-i1(n)).ge.3*correltime(n)) then
             istatus(n)=0
         endif
          
      endif
                     
4000  continue
      
      write(3,50)itime(1,j),itime(2,j),itime(3,j),itime(4,j),
     $ainter(1),istatus(1),ainter(2),istatus(2),ainter(3),istatus(3),
     $ainter(4),istatus(4),ainter(5),istatus(5),ainter(6),istatus(6),
     $ainter(7),istatus(7),ainter(8),istatus(8),ainter(9),istatus(9),
     $ainter(10),istatus(10),ainter(11),istatus(11)
50    format(i5,i4,2i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,6(f8.2,i2))

2002  Continue

      close(3)       
      write(*,*)'     Subroutine interpl_5min   DONE'	
      end

cccccccccccccccccccccc end of interpl_5min ccccccccccccccccccccccccccccccccccccc 


cccccccccccccccccccccc WG_5min.f ccccccccccccccccccccccccccccccccccccccccccccccc
ccc  Alin: Calculation of W pars is recursive, all but W3 converge quickly

c     This program is to calculate the W and G parameters
      subroutine WG_5min
      common wst(6)
      real kp,kp3
      integer dst
      dimension dat(4,60),idat(4,60),W(6),iW(6),bzn(6),S(6),bsn(6)
      dimension r(6),Wf(6),gamma(6),beta(6),xlumda(6),f(6)
      character str1*8,str2*12,aa      
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/    
      data beta/.8,.18,2.32,1.25,1.6,2.4/
      data xlumda/.39,.46,.39,.42,.41,1.29/
      data r/.39,.7,.031,.58,1.15,.88/
      data f/1.049,1.065,1.051,1.061,1.050,1.251/
      character finput1*30,finput2*30,fout*30,foutput*28
      character datechar*10,timechar*8,ychar*4,mchar*2,dchar*2
     
      do i=1,6
        W(i)=wst(i)/f(i)
        Wf(i)=0.0
      enddo

      pid4 = atan2(1.,1.)
      pi2 = 8*pid4
      n=0
      
      finput1='5min/inter5min.d'
      finput2='5min/kpkp3dst5min.d'
      fout='5min/WGparameters5min.d'      

      open(unit=1,file=finput1,status='old')
      open(unit=2,file=finput2,status='old')
c      open(unit=3,file=fout,status='unknown')
c      write(3,52) 'Year','Day','Hr','Min','ByIMF','BzIMF','V_SW',
c     $'Den_P','Pdyn','G1','G2','G3','8 status','kp','akp3','dst',
c     $'Bz1','Bz2','Bz3','Bz4','Bz5','Bz6',
c     $'W1','W2','W3','W4','W5','W6','6 stat'
c52    format(a5,a4,a3,a4,a6,7a7,a10,2a7,a6,6a7,6a9,a8)

      read(1,*)
      read(1,*)
      read(2,*)
      idayold=0
       
 10   read(2,*,err=20,end=20)iyr1,idy1,ihr1,imn1,kp,kp3,dst
      read(1,70,err=20,end=20)iyr,idy,ihr,imn,by,iby,bz,ibz,v,iv,den,id,
     $p,ip,bz1,i1,bz2,i2,bz3,i3,bz4,i4,bz5,i5,bz6,i6
       if (iyr1.ne.iyr) write (*,*) 'year mismatch in WG_5min input',
     & iyr1, iyr
 70   format(i5,i4,2i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,6(f8.2,i2))
 
      call daymonth(iyr,idy,jmo,iday)
      call utc(iyr,jmo,iday,ihr,imn,datechar,timechar,mchar,dchar)
      if (idy.ne.idayold) then
       if (idayold.ne.0) close(3)
       idayold=idy
       write (ychar,'(i4)') iyr
       foutput=ychar//'/QinDenton_'//ychar//mchar//dchar//'_5min'
       open(unit=3,file=foutput,status='unknown')
      end if

      n=n+1

cccccccc  G parameters cccccccccccccccccccccccccccccccccccc
      if(n.le.12) then
           dat(1,n)=by
           dat(2,n)=bz
           dat(3,n)=v
           dat(4,n)=den
           idat(1,n)=iby
           idat(2,n)=ibz
           idat(3,n)=iv
           idat(4,n)=id
      else
           do i=2,12
              dat(1,i-1)=dat(1,i)
              dat(2,i-1)=dat(2,i)
              dat(3,i-1)=dat(3,i)
              dat(4,i-1)=dat(4,i)
              idat(1,i-1)=idat(1,i)
              idat(2,i-1)=idat(2,i)
              idat(3,i-1)=idat(3,i)
              idat(4,i-1)=idat(4,i)
           enddo
           dat(1,12)=by
           dat(2,12)=bz
           dat(3,12)=v
           dat(4,12)=den
           idat(1,12)=iby
           idat(2,12)=ibz
           idat(3,12)=iv
           idat(4,12)=id
      endif
      ncap=min(n,12)
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      iG1=0
      iG2=0
      iG3=0
      
      do i=1,ncap
           bperp = sqrt(dat(1,i)**2+dat(2,i)**2)
           bp40 = bperp/40.
           HH = bp40**2/(1.+bp40)

           if(dat(1,i).eq.0..and.dat(2,i).eq.0.) then
                  theta = 0.
           else
                  theta = atan2(dat(1,i),dat(2,i))
                  if (theta.lt.0) theta = theta + pi2
           endif
           
           if (dat(2,i).lt.0.) then
                  Bs = -dat(2,i)
           else
                  Bs = 0.
           endif
               
           sine = (sin(theta/2.))**3
           sum1 = sum1 + dat(3,i)*HH*sine
           iG1=iG1+idat(1,i)+idat(2,i)+idat(3,i)
           sum2 = sum2 + Bs*dat(3,i)
           iG2=iG2+idat(2,i)+idat(3,i)
           sum3 = sum3 + Bs*dat(3,i)*dat(4,i)
           iG3=iG3+idat(2,i)+idat(3,i)+idat(4,i)
      enddo
      G1 = sum1/(ncap*1.0)
      G2 = sum2/(ncap*1.0)/200.
      G3 = sum3/(ncap*1.0)/2000.
      iG1=nint(iG1/(3.0*ncap))
      iG2=nint(iG2/(2.0*ncap))
      iG3=nint(iG3/(3.0*ncap))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccc W parameters ccccccccccccccccccccccccccccccccccccccccccccc
      bzn(1)=bz1
      bzn(2)=bz2
      bzn(3)=bz3
      bzn(4)=bz4
      bzn(5)=bz5
      bzn(6)=bz6
      do i=1,6
         if(bzn(i).lt.0.0) then
             bsn(i)=-bzn(i)      !bzn is used for bs temporarily            
         else
             bsn(i)=0.0
         endif

         S(i)=(den/5)**xlumda(i)*(v/400)**beta(i)*(bsn(i)/5)**gamma(i)
         W(i)=r(i)/12.0*S(i)+W(i)*exp(-r(i)/12.0)            ! previous W(i) are used !!!!
         nf=ibz*iv*id
         if(nf.eq.0) then
               Wf(i)=r(i)/12.0*S(i)*0.0+Wf(i)*exp(-r(i)/12.0)	
         elseif(nf.eq.8) then
               Wf(i)=r(i)/12.0*S(i)*2.0+Wf(i)*exp(-r(i)/12.0)
         else
               Wf(i)=r(i)/12.0*S(i)*1.0+Wf(i)*exp(-r(i)/12.0)
         endif
         
         if(W(i).ne.0) then
               iW(i)=nint(Wf(i)/W(i))
         else
               iW(i)=nint((ibz+iv+id)/3.0)
         endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c Alin: matching of years in input files was not inforced (is now), 
c resulting in combining KpDst data starting 1995 with OMNI W/G parameters starting 1981
c day-tag mismatch of leap and non-leap years entries led to some years being skipped 

c      if(idy1.eq.idy.and.ihr1.eq.ihr.and.imn1.eq.imn) then
c         write(3,50)iyr,idy,ihr,imn,by,bz,v,den,p,G1,G2,G3,
c     $iby,ibz,iv,id,ip,iG1,iG2,iG3,kp,kp3,dst,(bzn(i),i=1,6),
c     $(W(i)*f(i),i=1,6),(iW(i),i=1,6)
c50    format(I5,I4,2I3,2f7.2,f7.1,2f7.2,3f7.2,i3,7i1,
c     $f7.2,f7.2,i6,6f7.2,6f9.3,i3,5i1)
c      endif

      If (idy1.eq.idy.and.ihr1.eq.ihr.and.imn1.eq.imn) then
c       if (dst.ne.99999) 
       write (3,500)
     &  datechar,timechar,iyr,jmo,iday,ihr,imn,            
     &  by,bz,v,den,p,G1,G2,G3,
     &  iby,ibz,iv,id,ip,iG1,iG2,iG3,
     &  kp,kp3,dst,(bzn(i),i=1,6),(W(i)*f(i),i=1,6),(iW(i),i=1,6)
      End If
500    format(a10,'T',a8,i5,x,4i3,' 00',x,                
     &  2f7.2,f7.1,x,2f7.2,x,3f7.2,x,
     &  8i2,x,
     &  2f6.2,i6,x,6f7.2,x,6f7.3,x,6i2)
      
      Go To 10

20    close(1,status='delete')
      close(2,status='delete')
      close(3)
      write(*,*)'     Subroutine WG_5min        DONE'
      end
ccccccccccccccccccccccccc end of WG_5min ccccccccccccccccccccccccccccccccccccccc

      Subroutine daymonth(iyear,idy,jmo,iday)

      integer jdm(12)
      data jdm/31,28,31,30,31,30,31,31,30,31,30,31/

      jdm(2)=28
      if (mod(iyear,4).eq.0) jdm(2)=29
      iday=idy
      j=0
 10   j=j+1
      if (iday.le.jdm(j)) then
       jmo=j
      else
       iday=iday-jdm(j)
       go to 10
      end if

      return
      End
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Subroutine utc(iyr,jmo,iday,ihr,imn,datechar,timechar,mo,dy)

      character datechar*10,timechar*8
      character year*4,m*1,mo*2,d*1,dy*2,h*1,hr*2,n*1,nm*2

      write (year,'(i4)') iyr

      if (jmo.lt.10) then
        write (m,'(i1)') jmo
        mo='0'//m
       else
        write (mo,'(i2)') jmo
       end if

       if (iday.lt.10) then
        write (d,'(i1)') iday
        dy='0'//d
       else
        write (dy,'(i2)') iday
       end if

       if (ihr.lt.10) then
        write (h,'(i1)') ihr
        hr='0'//h
       else
        write (hr,'(i2)') ihr
       end if

       if (imn.lt.10) then
        write (n,'(i1)') imn
        nm='0'//n
       else
        write (nm,'(i2)') imn
       end if

       datechar=year//'-'//mo//'-'//dy
       timechar=hr//':'//nm//':00'

      return
      End


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c     ******   FOR 1 HOUR DATA  ********


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c This program calculates the 3 day average of Kp,Kp3,and add bs**gamma columns

      subroutine kp3bsgamma
      integer dst,cap
      character finput*30,fout*30,aa
      real Kp,akp3,Kpa(120),akp,weight,t3day 
      real N,V,P
      dimension bsg(6),gamma(6)
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/
      common wst(6)
      common /years/year0,nyrs
      
      finput='hour/omni2_hour.dat' 
      fout='hour/model_input.d' 
 
      open(unit=1,file=finput,err=33,status='old') 
      goto 34
 33   write(*,*)finput//' NOT EXIST'
      goto 1000
 34   open(unit=2,file=fout,status='unknown') 
 
      WRITE(2,400) 'YEAR',' DOY','HR','ByIMF','BzIMF','V_SW',
     $'Den_P','Pdyn','kp','kp3','dst','Bs1g','Bs2g','Bs3g','Bs4g',
     $'Bs5g','Bs6g'
400   format(a5,a4,a3,3a8,a7,a6,2a7,a6,6a8) 
 
      do I=1,120 
         Kpa(I)=0. 
      enddo 
      J=0      

10    read(1,300,err=20,end=20) iyr4,idoy,ihr,aa,by,bz,aa,N,
     $V,aa,P,aa,Kkp,aa,dst,aa
       IF (iyr4.lt.year0.or.iyr4.ge.year0+nyrs) go to 10  ! Alin
      Kp=Kkp
     
      
 300  format(2i4,i3,a61,2f6.1,a39,f6.1,f6.0,a18,f6.2,a59,i3,a4,i6,a63)
         if(Bz.eq.999.9) then

            do i=1,6
               bsg(i)=9999.99
            enddo
         elseif(Bz.lt.0.0) then
            do i=1,6
               bsg(i)=(-Bz)**gamma(i)
            enddo
         else
            do i=1,6
               bsg(i)=0.0
            enddo
         endif

      J=J+1
      if(Kp.eq.99) then
         akp3=99.
      else
         Kp = Kp/10. 
         if (J.lt.120) then 
            cap = J 
         else 
            cap = 120 
         endif 
         do K=1,119 
            Kpa(K)=Kpa(K+1) 
         enddo 
         Kpa(120)=Kp 
C    Calculate akp3 
         
         akp = 0 
         weight = 0 
         t3day = 72         !Time for 3 days in hours 

         do ii=1,cap 
            weight = weight + exp(-((ii-1)/t3day)) 
            akp = akp + exp(-((ii-1)/t3day))*Kpa(121-ii) 
         enddo
 
         akp3= akp/weight 
      endif
   
      WRITE(2,444)iyr4,idoy,ihr,By,Bz,V,N,p,Kp,akp3,dst,(bsg(j),j=1,6)
 444  FORMAT(I4,1x,I3,1x,I2,1x,2F8.2,F8.1,F7.2,F6.2,2(F6.2,1x),i6,6f8.2) 
      goto 10
 
 20   close(1) 
      CLOSE(2) 
1000  write(*,*)'     Subroutine kp3bsgamma     DONE'
      END 
cccccccccccccccccccc end of kp3bsgamma  cccccccccccccccccccccccccccccccccccccc


cccccccccc  subroutine av1year ccccccccccccccccccccccccccccccccccccccccccccccc

c    This program is to take the 1 year average of the hourly data

      subroutine av1year
      parameter (nipyr=1,nii=1000, nj=11)
      dimension y(nii),p(nii,nj),w(nii,nj),pv(nj)
      dimension idmy(7),dummy(24),dd(11),igap(3),it1(3),it2(3)
      double precision yrdouble
      character finput*30,fout*30
      common wst(6)
      common /years/year0,nyrs
 
      ni=nipyr*nyrs
      iy1=nint(year0)                       ! was 1963
      dyr= 1./nipyr
      do  i= 1,ni
         y(i)= (i-0.5)*dyr
      enddo
      do 23002 j= 1,nj
         do 23004 i= 1,ni
            p(i,j)= 0.
            w(i,j)= 0.
23004    continue
23002 continue
          
      finput='hour/model_input.d' 
      fout='hour/av1year.d'

      open(unit=11,file=finput,status='old')
 34   read(11,20,err=34,end=300)iyr4,idoy,ihr,(pv(j),j=1,5),
     $a,a,idst,(pv(j),j=6,nj)
 20   format(I4,1x,I3,1x,I2,1x,2F8.2,F8.1,F7.2,F6.2,2(F6.2,1x),i6,6f8.2)     
        
         if( mod(iyr4,4).eq.0 ) then
            doytot= 366.
         else
            doytot= 365.
         endif
         iyr4m= iyr4 - iy1

         yv= iyr4m + ( idoy -1. + ihr/24. )/doytot

         rii= (yv-y(1))/dyr + 1.
         ii= rii
         ip= ii + 1
         fip= rii - ii
         fii= 1. - fip

         do 23008 j= 1,nj

         if( (j.eq.1 .or. j.eq.2).and. pv(j).eq.999.9 ) go to 23008
            if( j.eq.3 .and. pv(j).eq.9999.0 ) go to 23008
            if( j.eq.4 .and. pv(j).eq.999.9 ) go to 23008
            if( j.eq.5 .and. pv(j).eq.99.99 ) go to 23008
            if( j.ge.6 .and. pv(j).eq.9999.99) go to 23008
            if( ii.ge.1 .and. ii.le.ni ) then
               p(ii,j)= p(ii,j) + fii*pv(j)
               w(ii,j)= w(ii,j) + fii
            endif
            if( ip.ge.1 .and. ip.le.ni ) then
               p(ip,j)= p(ip,j) + fip*pv(j)
               w(ip,j)= w(ip,j) + fip
            endif
23008    continue
c -- repeat
      go to 34

  300 close(11)
      do 23010 j= 1,nj
         do 23012 i= 1,ni
            if( w(i,j).eq. 0. ) then
               p(i,j)= 0.
            else
               p(i,j)= p(i,j)/w(i,j)
            endif
23012    continue
23010 continue
      
      open(unit=21,file=fout,status='unknown')
      do 23014 i= 1,ni
	   yrdouble=y(i)+iy1
          write(21,2100) yrdouble,(p(i,j),j=1,nj)
 2100       format(f9.4,11(1x,f9.4))
23014 continue
      close(21)
      write(*,*)'     Subroutine av1year        DONE'
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccc subroutine av20days ccccccccccccccccccccccccccccccccccccccccccccccccc

c    this program is to take the 20 days average of the hourly data

      subroutine av20days
      parameter (nipyr=18,nii=10000,nj=11)
      dimension y(nii),p(nii,nj),w(nii,nj),pv(nj)
      dimension idmy(7),dummy(24),dd(11),igap(3),it1(3),it2(3)
      double precision yrdouble
      character finput*30,fout*30
      common wst(6)
      common /years/year0,nyrs
 
      ni=nipyr*nyrs
      iy1=nint(year0)                ! was 1963
      dyr= 1./nipyr
      do  i= 1,ni
         y(i)= (i-0.5)*dyr
      enddo
      do 23002 j= 1,nj
         do 23004 i= 1,ni
            p(i,j)= 0.
            w(i,j)= 0.
23004    continue
23002 continue
          
      finput='hour/model_input.d' 
      fout='hour/av20days.d'
 
      open(unit=11,file=finput,status='old')
 34   read(11,20,err=34,end=300)iyr4,idoy,ihr,(pv(j),j=1,5),
     $a,a,idst,(pv(j),j=6,nj)
20    format(I4,1x,I3,1x,I2,1x,2F8.2,F8.1,F7.2,F6.2,2(F6.2,1x),i6,6f8.2)     
        
         if( mod(iyr4,4).eq.0 ) then
            doytot= 366.
         else
            doytot= 365.
         endif
         iyr4m= iyr4 - iy1

         yv= iyr4m + ( idoy -1. + ihr/24. )/doytot

         rii= (yv-y(1))/dyr + 1.
         ii= rii
         ip= ii + 1
         fip= rii - ii
         fii= 1. - fip

         do 23008 j= 1,nj

            if( (j.eq.1 .or. j.eq.2).and. pv(j).eq.999.9 ) go to 23008
            if( j.eq.3 .and. pv(j).eq.9999.0 ) go to 23008
            if( j.eq.4 .and. pv(j).eq.999.9 ) go to 23008
            if( j.eq.5 .and. pv(j).eq.99.99 ) go to 23008
            if( j.ge.6 .and. pv(j).eq.9999.99) go to 23008
            if( ii.ge.1 .and. ii.le.ni ) then
               p(ii,j)= p(ii,j) + fii*pv(j)
               w(ii,j)= w(ii,j) + fii
            endif
            if( ip.ge.1 .and. ip.le.ni ) then
               p(ip,j)= p(ip,j) + fip*pv(j)
               w(ip,j)= w(ip,j) + fip
            endif
         
23008   continue
      go to 34

 300  close(11)
      do 23010 j= 1,nj
         do 23012 i= 1,ni
            if( w(i,j).eq. 0. .or. w(i,j).lt.120. ) then
               p(i,j)= 999.
            else
               p(i,j)= p(i,j)/w(i,j)
            endif
23012    continue
23010 continue
      open(unit=21,file=fout,status='unknown')
      do 23014 i= 1,ni
	   yrdouble=y(i)+iy1
          write(21,2100) yrdouble,(p(i,j),j=1,nj)
 2100       format(f9.4,11(1x,f9.4))
23014 continue
      close(21)
      write(*,*)'     Subroutine av20days       DONE'
      end
ccccccccccccccc end of av20days   ccccccccccccccccccccccccccccccccccc


cccccccccc subroutine inter20days ccccccccccccccccccccccccccccccccccc

c    This program is to interpolate across the gaps in the 20 day average data and
c    calculate bz1-bz6 

      subroutine inter20days
      common wst(6)
      parameter(nc=11)
      dimension av1y(nc+1,99),istatus(nc),i1(nc),i2(nc), 
     $ainter(nc),check1(nc),check2(nc),check(nc),gamma(6)
      dimension data20d(nc+1,1800)
      double precision tyear
      dimension idmy(7),dummy(24),correltime(nc),itime(4,1800)
      character finput2*30,finput1*30,fout*30
      data correltime/0.48,0.32,0.87,0.72,0.77,
     $0.39,0.39,0.38,0.38,0.39,0.38/
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      finput2='hour/av20days.d'
      finput1='hour/av1year.d'
      fout='hour/inter20days.d'
      
cccccccc  input the 20 day average ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=1,file=finput1,status='old')
      do i=1,99
      read(1,*,end=35)(av1y(j,i),j=1,12)
      enddo
 35   close(1,status='delete')
      nd=i-1
      
ccccccc  input the original data   cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      open(unit=2,file=finput2,status='old')
      
      do n=1,9999
      read(2,40,err=37,end=37) (data20d(i,n),i=1,12)
 40   format(f9.4,11(1x,f9.4))
      enddo
 37   close(2,status='delete')
      nm=n-1
      
ccccccccccccccccccccc If gap is at the beginning or the end ccccccccccccccc
 64   do n=1,nc
          check1(n)=0
          check2(n)=0
          check(n)=0
      enddo
      do k=2,12
         n=k-1
         if(data20d(k,1).eq.999.0) then      !if there is a gap at the beginning
                data20d(k,1)=av1y(k,1)
                check1(n)=1
         endif
         if(data20d(k,nm).eq.999.0) then      !if there is a gap at the end
                data20d(k,nm)=av1y(k,nd)
                check2(n)=1
         endif
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      open(unit=3,file=fout,status='unknown')
      write(3,52) 'Year','ByIMF','BzIMF',
     $'V_SW','Den_P','Pdyn','Bz1','Bz2',
     $'Bz3','Bz4','Bz5','Bz6'
52    format(a9,11(1x,a9))

      do 2002 j=1,nm
        do 4000 k=2,12
        n=k-1
        if(data20d(k,j).ne.999.0) then
            ainter(n)=data20d(k,j)
            check(n)=0
       else

         if(check(n).eq.0)  then
            check(n)=1       
            i1(n)=j-1
            i2(n)=j+1                 
              DO WHILE (data20d(k,i1(n)) .eq. 999.0)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (data20d(k,i2(n)) .eq. 999.0)
                 i2(n)=i2(n)+1
              END DO     
     	 endif  
                	   
            timegap=(i2(n)-i1(n)-1)/correltime(n)
           if(timegap .le. 2.0)	then
             ainter(n)=(j-i1(n))*(data20d(k,i2(n))-data20d(k,i1(n)))/
     $(i2(n)-i1(n))+data20d(k,i1(n))
             istatus(n)=1
            endif

           if((timegap .gt. 2.0) .and. (timegap .le. 4.0)) then
              amidtime=(i2(n)+i1(n))/2.0
              if(abs(j-amidtime) .le. correltime(n))  then
                data20d_mid=(data20d(k,i1(n))+data20d(k,i2(n)))/2.0
                ainter(n)=data20d_mid+(j-amidtime)*
     $(data20d(k,i2(n))-data20d(k,i1(n)))/(2*correltime(n))
              else if((j-i1(n)) .lt. correltime(n)) then
                ainter(n)=data20d(k,i1(n))
              else
	        ainter(n)=data20d(k,i2(n))
              endif
              istatus(n)=1
           endif
           !if(j.le.3)write(*,*)data20d(1,j),timegap
           if(timegap .gt. 4.0)	then
              ii=int(correltime(n))+1
              if((j-i1(n)) .lt. ii) ainter(n)=data20d(k,i1(n))
       	      if((i2(n)-j) .lt. ii) ainter(n)=data20d(k,i2(n))
              istatus(n)=1
       	      if((j-i1(n)) .ge. ii .and. (i2(n)-j) .ge. ii) then
               w1=max(1-(j-i1(n)-correltime(n))/(2.0*correltime(n)),0.0)
               w2=max(1-(i2(n)-correltime(n)-j)/(2.0*correltime(n)),0.0)
    	         wav=1-w1-w2
                 fav=0.0
                 if(data20d(1,j).lt. av1y(1,1)) fav=av1y(k,1)
                 if(data20d(1,j).gt. av1y(1,nd)) fav=av1y(k,nd)
                 do  i3=1,nd-1
                    if( data20d(1,j).ge.av1y(1,i3) .and.data20d(1,j).le.
     $av1y(1,i3+1))  then
                     fav=(data20d(1,j)-av1y(1,i3))*(av1y(k,i3+1)
     $-av1y(k,i3))/(av1y(1,i3+1)-av1y(1,i3))+av1y(k,i3)
                     
                    goto 36
                    endif
                 enddo	 
 36           ainter(n)=w1*data20d(k,i1(n))+w2*data20d(k,i2(n))+wav*fav
                  if(ainter(n) .eq. fav)  istatus(n)=0
               endif
         endif
         if(check1(n).eq.1.and.i1(n).eq.1
     $.and.(i2(n)-j).ge.3*correltime(n)) then
             istatus(n)=0
         endif

         if(check2(n).eq.1.and.i2(n).eq.nm
     $.and.(j-i1(n)).ge.3*correltime(n)) then
             istatus(n)=0
         endif
	          
      endif
                     
4000  continue
         do k=6,11
            ainter(k)=-ainter(k)**(1.0/gamma(k-5))
         enddo
      write(3,50)data20d(1,j),
     $ainter(1),ainter(2),ainter(3),
     $ainter(4),ainter(5),ainter(6),
     $ainter(7),ainter(8),ainter(9),
     $ainter(10),ainter(11)
50    format(f9.4,11(1x,f9.4))
2002	continue					
      close(3)
      write(*,*)'     Subroutine inter20days    DONE'
      end
ccccccccccccccccc  end of inter20days   ccccccccccccccccccccccccccccccccccccccccccc 


ccccccccccc subroutine interhour cccccccccccccccccccccccccccccccccccccccccccccccccc

c     This program is to interpolate across the gaps in the hourly data     

      subroutine interhour
      common wst(6)
      parameter(nc=11)
      dimension av20(nc+1,1800),istatus(nc),i1(nc),i2(nc), 
     $ainter(nc),check1(nc),check2(nc),check(nc)
      double precision tyear,datahr(nc+1,550000)
      dimension dummy(12),idst(550000)
      dimension akp(550000),akp3(550000),correltime(nc),itime(3,550000)
      character finput2*30,finput1*30,yeard*5,fout*30
     $,finput3*30,finput4*30
      data correltime/1.46,0.73,13.50,2.18,1.68,
     $0.73,0.73,0.73,0.73,0.73,0.73/ 
            
      finput2='hour/model_input.d'
      finput1='hour/inter20days.d'
      fout='hour/interhour.d'
      n=0

cccccccc  input the 20 days average ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=1,file=finput1,status='old')
      read(1,*)
      do i=1,9999
      read(1,*,end=34)(av20(j,i),j=1,12)
      enddo
 34   close(1,status='delete')
      nd=i-1  

ccccccc  input the original data   cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      open(unit=2,file=finput2,status='old')
      read(2,*)
      do i=1,999999
          n=n+1
      read(2,40,err=20,end=20)iyr4,idoy,ihr,(dummy(j),j=2,6),
     $akp(n),akp3(n),idst(n),(dummy(j),j=7,12)
40    format(I4,1x,I3,1x,I2,1x,2F8.2,F8.1,F7.2,F6.2,2(F6.2,1x),i6,6f8.2)
        
         itime(1,n)=iyr4
         itime(2,n)=idoy
         itime(3,n)=ihr
         
         call year(iyr4,idoy,ihr,30,tyear)
         datahr(1,n)=tyear
         datahr(2,n)=dummy(2)
         datahr(3,n)=dummy(3)
         datahr(4,n)=dummy(4)
         datahr(5,n)=dummy(5)
         datahr(6,n)=dummy(6)
         datahr(7,n)=dummy(3)
         datahr(8,n)=dummy(3)
         datahr(9,n)=dummy(3)
         datahr(10,n)=dummy(3)
         datahr(11,n)=dummy(3)
         datahr(12,n)=dummy(3)
      enddo
 20   close(2,status='delete')

      nm=n-1
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do n=1,nc
          check1(n)=0
          check2(n)=0
          check(n)=0
      enddo
      do k=2,12
         n=k-1
         if( ( ( k.eq.2 .or. k.eq.3) .and. datahr(k,1).eq.999.9)
     &.or. ( k.eq.4 .and. datahr(k,1).eq.9999.0 ) .or. ( k
     &.eq.5 .and. datahr(k,1).eq.999.9 ) .or. ( k.eq.6 .and. 
     $datahr(k,1).eq.99.99).or.
     $((k.ge.7).and.datahr(k,1).eq.999.9)) then      !if there is a gap at the beginning
                datahr(k,1)=av20(k,1)
                check1(n)=1
         endif
         if( ( ( k.eq.2 .or. k.eq.3).and. datahr(k,nm).eq.999.9)
     &.or. ( k.eq.4 .and. datahr(k,nm).eq.9999. ) .or. ( k
     &.eq.5 .and. datahr(k,nm).eq.999.9 ) .or. ( k.eq.6 .and. 
     $datahr(k,nm).eq.99.99).or.
     $((k.ge.7).and.datahr(k,nm).eq.999.9) ) then      !if there is a gap at the end
                datahr(k,nm)=av20(k,nd)
                check2(n)=1
         endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      open(unit=3,file=fout,status='unknown')
      write(3,*) 's=2-->original           s=1-->interpolation 
     $s=0-->1 day average'
      write(3,52) 'Year','Day','Hr','ByIMF','s','BzIMF','s',
     $'V_SW','s','Den_P','s','Pdyn','s','Kp','Kp3','dst','Bz1','s',
     $'Bz2','s','Bz3','s','Bz4','s','Bz5','s','Bz6','s'
52    format(a5,a4,a3,2(a8,a2),a8,a2,a7,a2,a6,a2,2a6,a6,6(a8,a2))

      do 2002 j=1,nm

        do 4000 k=2,12

        n=k-1
       if( (( k.eq.2 .or. k.eq.3).and. datahr(k,j).ne.999.9).or.
     &(k.eq.4 .and. datahr(k,j).ne.9999. ) .or. 
     $(k.eq.5 .and. datahr(k,j).ne.999.9 ) .or. 
     $(k.eq.6 .and. datahr(k,j).ne.99.99) .or.
     $(k.ge.7 .and. datahr(k,j).ne.999.9)  ) then 
            istatus(n)=2
            
            if(check1(n).eq.1.and.j.eq.1) then
               istatus(n)=0
            endif

            if(check2(n).eq.1.and.j.eq.nm) then
               istatus(n)=0
            endif
            ainter(n)=datahr(k,j)
            check(n)=0
       else
             
         if(check(n).eq.0)  then
            check(n)=1       
            i1(n)=j-1
            i2(n)=j+1
                      
            if( k.eq.2 .or. k.eq.3 ) then
       	      DO WHILE (datahr(k,i1(n)) .eq. 999.9)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datahr(k,i2(n)) .eq. 999.9)
                 i2(n)=i2(n)+1
              END DO
            endif
         
            if( k.eq.4 ) then
              DO WHILE (datahr(k,i1(n)) .eq. 9999.)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datahr(k,i2(n)) .eq. 9999.)
                 i2(n)=i2(n)+1
              END DO
            endif
     
            if( k.eq.5 ) then
              DO WHILE (datahr(k,i1(n)) .eq. 999.9)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datahr(k,i2(n)) .eq. 999.9)
                 i2(n)=i2(n)+1
              END DO
            endif
    
            if( k.eq.6 ) then
              DO WHILE (datahr(k,i1(n)) .eq. 99.99)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datahr(k,i2(n)) .eq. 99.99)
                 i2(n)=i2(n)+1
              END DO
            endif

            if( k.ge.7 ) then
              DO WHILE (datahr(k,i1(n)) .eq. 999.9)
                 i1(n)=i1(n)-1
              END DO
              DO WHILE (datahr(k,i2(n)) .eq. 999.9)
                 i2(n)=i2(n)+1
              END DO
            endif
	 endif  

        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
                	   
            timegap=(i2(n)-i1(n)-1)/correltime(n)
           if(timegap .le. 2.0)	then
             ainter(n)=(j-i1(n))*(datahr(k,i2(n))-datahr(k,i1(n)))/
     $(i2(n)-i1(n))+datahr(k,i1(n))
             istatus(n)=1
            endif

           if((timegap .gt. 2.0) .and. (timegap .le. 4.0)) then
              amidtime=(i2(n)+i1(n))/2.0
              if(abs(j-amidtime) .le. correltime(n))  then
                datahr_mid=(datahr(k,i1(n))+datahr(k,i2(n)))/2.0
                ainter(n)=datahr_mid+(j-amidtime)*
     $(datahr(k,i2(n))-datahr(k,i1(n)))/(2*correltime(n))
              else if((j-i1(n)) .lt. correltime(n)) then
                ainter(n)=datahr(k,i1(n))
              else
	        ainter(n)=datahr(k,i2(n))
              endif
              istatus(n)=1
           endif

           if(timegap .gt. 4.0)	then
              ii=int(correltime(n))+1
              if((j-i1(n)) .lt. ii) ainter(n)=datahr(k,i1(n))
       	      if((i2(n)-j) .lt. ii) ainter(n)=datahr(k,i2(n))
              istatus(n)=1
       	      if((j-i1(n)) .ge. ii .and. (i2(n)-j) .ge. ii) then
               w1=max(1-(j-i1(n)-correltime(n))/(2.0*correltime(n)),0.0)
               w2=max(1-(i2(n)-correltime(n)-j)/(2.0*correltime(n)),0.0)
    	         wav=1-w1-w2
                 fav=0.0
	         if(datahr(1,j).lt. av20(1,1)) fav=av20(k,1)
                 if(datahr(1,j).gt. av20(1,nd)) fav=av20(k,nd)
                 do  i3=1,99999
                    if( datahr(1,j).ge.av20(1,i3) .and.datahr(1,j).le.
     $av20(1,i3+1))  then
                     fav=(datahr(1,j)-av20(1,i3))*(av20(k,i3+1)
     $-av20(k,i3))/(av20(1,i3+1)-av20(1,i3))+av20(k,i3)
                     goto 35
                    endif
                 enddo	                  
 35           ainter(n)=w1*datahr(k,i1(n))+w2*datahr(k,i2(n))+wav*fav
                  if(ainter(n) .eq. fav)  istatus(n)=0
               endif
         endif
        
         if(check1(n).eq.1.and.i1(n).eq.1
     $.and.(i2(n)-j).ge.3*correltime(n)) then
             istatus(n)=0
         endif

         if(check2(n).eq.1.and.i2(n).eq.nm
     $.and.(j-i1(n)).ge.3*correltime(n)) then
             istatus(n)=0
         endif
	          
      endif
                     
4000  continue
      
      write(3,50)itime(1,j),itime(2,j),itime(3,j),
     $ainter(1),istatus(1),ainter(2),istatus(2),ainter(3),istatus(3),
     $ainter(4),istatus(4),ainter(5),istatus(5),
     $akp(j),akp3(j),idst(j),ainter(6),istatus(6),
     $ainter(7),istatus(7),ainter(8),istatus(8),ainter(9),istatus(9),
     $ainter(10),istatus(10),ainter(11),istatus(11)
50    format(i5,i4,i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,
     $2f6.2,i6,6(f8.2,i2))

2002  Continue

      close(3)
      write(*,*)'     Subroutine interhour      DONE'
      end
       
cccccccccccccc end of interhour  ccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccc subroutine WGhour  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     this program is to calculate the W and G parameters

      subroutine WGhour
      common wst(6)
      real kp,kp3,Na,Nb,Ng,Nw
      integer dst
      dimension W(6),iW(6),bzna(6),S(6,4),bznb(6)
      dimension r(6),Wf(6),gamma(6),beta(6),xlumda(6),f(6),bznw(6)
      character str1*8,str2*12  
      data f/1.069748,1.024995,1.247411,1.100432,0.9768723,1.134791/   
      data gamma/0.87,0.67,1.32,1.29,0.69,0.53/    
      data beta/.8,.18,2.32,1.25,1.6,2.4/
      data xlumda/.39,.46,.39,.42,.41,1.29/
      data r/.39,.7,.031,.58,1.15,.88/
      character finput*30,fout*30,aa,foutput*28
      character datechar*10,timechar*8,ychar*4,mchar*2,dchar*2
      
      do i=1,6
        W(i)=wst(i)/f(i)
        Wf(i)=0.0
      enddo
       
      pid4 = atan2(1.,1.)
      pi2 = 8*pid4

      finput='hour/interhour.d'
      fout='hour/WGhour.d'      

      open(unit=1,file=finput,status='old')
c      open(unit=3,file=fout,status='unknown')
c      write(3,52) 'Year','Day','Hr','ByIMF','BzIMF','V_SW',
c     $'Den_P','Pdyn','G1','G2','G3','8 status','kp','akp3','dst'
c     $,'Bz1','Bz2','Bz3','Bz4','Bz5','Bz6','W1'
c     $,'W2','W3','W4','W5','W6','6 stat'
c52    format(a5,a4,a3,a7,7a7,a10,2a7,a6,6a7,6a9,a8)
     
      read(1,*)
      read(1,*)
      Bya = 999.9
      idayold=0
      imn=0
             
 10   read(1,70,err=20,end=20)iyr,idy,ihr,by,iby,bz,ibz,v,iv,den,id,
     $p,ip,kp,kp3,dst,bz1,i1,bz2,i2,bz3,i3,bz4,i4,bz5,i5,bz6,i6
 70   format(i5,i4,i3,2(f8.2,i2),f8.1,i2,f7.2,i2,f6.2,i2,
     $2f6.2,i6,6(f8.2,i2))

      call daymonth(iyr,idy,jmo,iday)
      call utc(iyr,jmo,iday,ihr,imn,datechar,timechar,mchar,dchar)
      if (idy.ne.idayold) then
       if (idayold.ne.0) close(3)
       idayold=idy
       write (ychar,'(i4)') iyr
       foutput=ychar//'/QinDenton_'//ychar//mchar//dchar//'_hour'
       open(unit=3,file=foutput,status='unknown')
      end if

         n=n+1
         Byb = by
         Bzb = bz
         Nb = den
         Vb = v
         bznb(1)=bz1
         bznb(2)=bz2
         bznb(3)=bz3
         bznb(4)=bz4
         bznb(5)=bz5
         bznb(6)=bz6
         if (Bya.eq.999.9) then
               Bya=Byb
               Bza=Bzb
               Na=Nb
               Va=Vb
               bzna(1)=bznb(1)
               bzna(2)=bznb(2)
               bzna(3)=bznb(3)
               bzna(4)=bznb(4)
               bzna(5)=bznb(5)
               bzna(6)=bznb(6) 
         endif

cccccccc  G parameters cccccccccccccccccccccccccccccccccccc
         dBy = (Byb-Bya)/12.
         dBz = (Bzb-Bza)/12.
         dN = (Nb-Na)/12.
         dV = (Vb-Va)/12.
         sum1 = 0.
         sum2 = 0.
         sum3 = 0.
         iG1=0
         iG2=0
         iG3=0
         
            do 45 i=1,12
               Byg = Bya + dBy*i
               Bzg = Bza + dBz*i
               Ng = Na + dN*i
               Vg = Va + dV*i
               bperp = sqrt(Byg**2+Bzg**2)
               bp40 = bperp/40.
               HH = bp40**2/(1.+bp40)
               if(Byg.eq.0..and.Bzg.eq.0.) then
                  theta = 0.
               else
                  theta = atan2(Byg,Bzg)
                  if (theta.lt.0) theta = theta + pi2
               endif
               if (Bzg.lt.0.) then
                  Bs = -Bzg
               else
                  Bs = 0.
               endif
               sine = (sin(theta/2.))**3
               sum1 = sum1 + Vg*HH*sine
               sum2 = sum2 + Bs*Vg
               sum3 = sum3 + Bs*Vg*Ng
 45         continue
            G1 = sum1/12.
            G2 = sum2/12./200.
            G3 = sum3/12./2000.
            iG1=nint((iby+ibz+iv)/3.)
            iG2=nint((ibz+iv)/2.)
            iG3=nint((ibz+iv+id)/3.)
             
cccccccccccccc W parameters ccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,6
      do j=1,4
         Nw=Nb-(Nb-Na)/4.0*(j-1)
         Vw=Vb-(Vb-Va)/4.0*(j-1)
         bznw(i)=bznb(i)-(bznb(i)-bzna(i))/4.0*(j-1) 
         bznw(i)=(abs(bznw(i))-bznw(i))/2.0
        S(i,j)=(Nw/5)**xlumda(i)*(Vw/400)**beta(i)*(bznw(i)/5)**gamma(i)
      enddo
      enddo
     
      do k=1,6  
         W(k)=0.25*r(k)*(S(k,1)+S(k,2)*exp(-0.25*r(k))+S(k,3)
     $*exp(-0.5*r(k))+S(k,4)*exp(-0.75*r(k)))+W(k)*exp(-r(k))
         nf=ibz*iv*id
         if(nf.eq.0) then
           Wf(k)=0.25*r(k)*(S(k,1)+S(k,2)*exp(-0.25*r(k))+S(k,3)
     $*exp(-0.5*r(k))+S(k,4)*exp(-0.75*r(k)))*0.0+Wf(k)*exp(-r(k))
         elseif(nf.eq.8) then
           Wf(k)=0.25*r(k)*(S(k,1)+S(k,2)*exp(-0.25*r(k))+S(k,3)
     $*exp(-0.5*r(k))+S(k,4)*exp(-0.75*r(k)))*2.0+Wf(k)*exp(-r(k))
         else
           Wf(k)=0.25*r(k)*(S(k,1)+S(k,2)*exp(-0.25*r(k))+S(k,3)
     $*exp(-0.5*r(k))+S(k,4)*exp(-0.75*r(k)))*1.0+Wf(k)*exp(-r(k))
         endif
         if(W(k).ne.0) then
               iW(k)=nint(Wf(k)/W(k))
         else
               iW(k)=nint((ibz+iv+id)/3.0)
         endif
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
               Bya=Byb
               Bza=Bzb
               Na=Nb
               Va=Vb
               bzna(1)=bznb(1)
               bzna(2)=bznb(2)
               bzna(3)=bznb(3)
               bzna(4)=bznb(4)
               bzna(5)=bznb(5)
               bzna(6)=bznb(6)
         
c         write(3,50)iyr,idy,ihr,by,bz,v,den,p,G1,G2,G3,
c     $iby,ibz,iv,id,ip,iG1,iG2,iG3,kp,kp3,dst,(bznb(i),i=1,6),
c     $(W(i)*f(i),i=1,6),(iW(i),i=1,6)
c50    format(I5,I4,I3,2f7.2,f7.1,2f7.2,3f7.2,i3,7i1,
c     $f7.2,f7.2,i6,6f7.2,6f9.3,i3,5i1)
 
cx       if (dst.ne.99999) 
        write (3,500)
     &  datechar,timechar,iyr,jmo,iday,ihr,imn,
     &  by,bz,v,den,p,G1,G2,G3,
     &  iby,ibz,iv,id,ip,iG1,iG2,iG3,
     &  kp,kp3,dst,(bznb(i),i=1,6),(W(i)*f(i),i=1,6),(iW(i),i=1,6)
500    format(a10,'T',a8,i5,x,4i3,' 00',x,
     &  2f7.2,f7.1,x,2f7.2,x,3f7.2,x,
     &  8i2,x,
     &  2f6.2,i6,x,6f7.2,x,6f7.3,x,6i2)
      
      Go to 10

20    close(1,status='delete')
      close(3)
     
      write(*,*)'     Subroutine WGhour         DONE'
	     
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
