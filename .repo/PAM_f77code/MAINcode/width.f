      ome1=0.5d0*(ep0+ed0+D+sqrt((delta-D)**2 +4.d0*V*V))
         ome1apx=ed0+V**2/(delta-D)
         
         ome2=0.5d0*(ep0+ed0+D-sqrt((delta-D)**2 +4.d0*V*V))
         ome2apx= ep0+D-V**2/(delta-D)
         
         ome3=0.50d0*(ep0-ed0-D+sqrt((delta+D)**2 +4.d0*V**2))
         ome3apx=ed0+V**2/(delta+D)
         
         ome4=0.50d0*(ep0-ed0-D-sqrt((delta+D)**2 +4.d0*V**2))
         ome4apx=ep0-D-V**2/(delta+D)
         if(1.eq.2)then
            write(56,*) ome1-gzerop,0.d0
            write(56,*) ome1-gzerop,2.d0
            write(57,*) ome1apx-gzerop,0.d0
            write(57,*) ome1apx-gzerop,2.d0
            write(56,*) '            '
            write(57,*) '            '
            write(56,*) ome2-gzerop,0.d0
            write(56,*) ome2-gzerop,2.d0
            write(57,*) ome2apx-gzerop,0.d0
            write(57,*) ome2apx-gzerop,2.d0
            write(56,*) '            '
            write(57,*) '            '
            write(56,*) ome3-gzerop,0.d0
            write(56,*) ome3-gzerop,2.d0
            write(57,*) ome3apx-gzerop,0.d0
            write(57,*) ome3apx-gzerop,2.d0
            write(56,*) '            '
            write(57,*) '            '
            write(56,*) ome4-gzerop,0.d0
            write(56,*) ome4-gzerop,2.d0
            write(57,*) ome4apx-gzerop,0.d0
            write(57,*) ome4apx-gzerop,2.d0
            write(56,*) '            '
            write(57,*) '            '
         endif
         write(58,*) abs(ome2-ome4)
         write(59,*) abs(ome3-ome1)
         write(60,*) ((ome1+ome3)/2.d0 -(ome4+ome2)/2.d0)
