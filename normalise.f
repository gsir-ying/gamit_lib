       subroutine normalise (vect,norm_vect)

c
c  Purpose: to normalise a 3 x 1 vector (vect) into norm_vect
c
c  IN: vect - vector to be normalised   R*8 (3)
c
c OUT: norm_vect - normalised vector    R*8 (3)
c
       implicit none
c
       real*8 vect(3),norm_vect(3),len,amag3
       integer i

       len = amag3(vect)
       do 10 i=1,3
         norm_vect(i) = vect(i)/len
10     continue

       return
       end
