      PROGRAM getentr

      implicit real*8(a-h,o-z)
      integer L,i,k,R

      real*8 temp(1000),cv(1000),dtemp,en(1000)

      open(unit=20,file='Cv.ed',status='old')
      do i=1,1000
         read(20,*,end=999) temp(i),cv(i)
         L=i
      enddo
 999  continue
      dtemp=temp(2)-temp(1)
      print*,L
      print*,dtemp
      open(34,file='entropy.ed')
      do i=1,L
         do k=1,i
            en(i)=en(i)+cv(k)/temp(k)*dtemp
         enddo
         write(34,*)temp(i),en(i)
      enddo
      end
      
