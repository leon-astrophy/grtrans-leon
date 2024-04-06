      SUBROUTINE dhunt(xx,x,jlo)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: jlo
      DOUBLE PRECISION, INTENT(IN) :: x
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xx
      INTEGER :: n,inc,jhi,jm
      LOGICAL :: ascnd
!      write(6,*) 'hunt: ',size(xx)
!      write(6,*) 'hunt: ',x,jlo
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      if (jlo <= 0 .or. jlo > n) then
         jlo=0
         jhi=n+1
         else
            inc=1
            if (x >= xx(jlo) .eqv. ascnd) then
               do
                  jhi=jlo+inc
                  if (jhi > n) then
                     jhi=n+1
                     exit
                     else
                        if (x < xx(jhi) .eqv. ascnd) exit
                        jlo=jhi
                        inc=inc+inc
                        end if
                        end do
                        else
                           jhi=jlo
                           do
                              jlo=jhi-inc
                              if (jlo < 1) then
                                 jlo=0
                                 exit
                                 else
                                    if (x >= xx(jlo) .eqv. ascnd) exit
                                    jhi=jlo
                                    inc=inc+inc
                                 end if
                                 end do
                                 end if
                                 end if
                                 do
                                    if (jhi-jlo <= 1) then
                                       if (x == xx(n)) jlo=n-1
                                       if (x == xx(1)) jlo=1
                                       exit
                                 else
                                       jm=(jhi+jlo)/2
                                    if (x >= xx(jm) .eqv. ascnd) then
                                                jlo=jm
                                    else
                                       jhi=jm
              end if   
            end if                                      
          end do
        END SUBROUTINE dhunt

      SUBROUTINE hunt(xx,x,jlo)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: jlo
      REAL, INTENT(IN) :: x
      REAL, DIMENSION(:), INTENT(IN) :: xx
      INTEGER :: n,inc,jhi,jm
      LOGICAL :: ascnd
!      write(6,*) 'hunt: ',size(xx)
!      write(6,*) 'hunt: ',x,jlo
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      if (jlo <= 0 .or. jlo > n) then
         jlo=0
         jhi=n+1
         else
            inc=1
            if (x >= xx(jlo) .eqv. ascnd) then
               do
                  jhi=jlo+inc
                  if (jhi > n) then
                     jhi=n+1
                     exit
                     else
                        if (x < xx(jhi) .eqv. ascnd) exit
                        jlo=jhi
                        inc=inc+inc
                        end if
                        end do
                        else
                           jhi=jlo
                           do
                              jlo=jhi-inc
                              if (jlo < 1) then
                                 jlo=0
                                 exit
                                 else
                                    if (x >= xx(jlo) .eqv. ascnd) exit
                                    jhi=jlo
                                    inc=inc+inc
                                 end if
                                 end do
                                 end if
                                 end if
                                 do
                                    if (jhi-jlo <= 1) then
                                       if (x == xx(n)) jlo=n-1
                                       if (x == xx(1)) jlo=1
                                       exit
                                 else
                                       jm=(jhi+jlo)/2
                                    if (x >= xx(jm) .eqv. ascnd) then
                                                jlo=jm
                                    else
                                       jhi=jm
              end if   
            end if                                      
          end do
        END SUBROUTINE hunt

