        subroutine printmat(mat,m,n,comment)

c  prints out a m x n matrix

        implicit none
        integer m,n,i,j,indx
        real*8 mat(m,n)
        character*80 comment

c  search for the first blank in the comment - only print up to it
        indx = index(comment,' ')
        write(*,'(a)')comment(1:indx)
        do 10 i=1,m
            write(*,*)(mat(i,j),j=1,n)
10      continue

        print*,' '
c        write(*,'(a)')' press return to continue : '
c        read(*,'(a)')ent

        return
        end

