
C     ********************************************************************
C     *                                                                  *
C     *          GENERATE INITIAL ENSEMBLE FROM A SET OF EOFs            *
C     *                                                                  *
C     ********************************************************************


C -->  N is the number of members to generate;
C -->  dimstate = state dimension;

      INTEGER  N, dimstate
      PARAMETER (N = 10, dimstate = 107376)

      REAL*4   Ens(dimstate,N), stateA(dimstate)

      REAL*4   Minv(N-1,N-1), d(N-1), v(N-1)

      REAL*8   Omega(N-1,N-1)
      
      REAL*8   s, r, ss

      INTEGER nhaz


C -->  Initial Ensemmble is generated with
C
C        Ens = stateA + sqrt(N)*Lb*Omega(t);


C ==>  << Read EOFs outputs: stateA and L >>

      open(1,file='mean_state.bin',
     &     access = 'direct', form = 'unformatted', recl = 4*dimstate)
      read(1,rec=1) stateA

      open(2,file='eof_basis.bin',
     &     access = 'direct', form = 'unformatted', recl = 4*dimstate)
      do i = 1, N-1
         read(2,rec=i) (Ens(j,i),j = 1,dimstate)
      end do



C ==>  << Open output file >>

      open(3,file='ini_ens.bin',
     &     access = 'direct', form = 'unformatted', recl = 4*dimstate) 
      


C ==>  << Generate random orthogonal matrix Omega >>

C -->  Initialize nhaz = 100;

      nhaz = 100


C -->  Compute Omega;

      call randmat(Omega, N-1, nhaz)


C -->  Compute Omega = H(1/sqrt(N))*Omega;

      r = N - sqrt(real(N))

      do j = 1,N-1
         s = 0.0
         do i = 1,N-1
            s = s + Omega(i,j)
         end do
         s = s/r
         do i = 1,N-1
            Omega(i,j) = Omega(i,j) - s
         end do
      end do
      


C ==>  << Compute Omega*M(1/2) = Omega*Minv(1/2)(-1) >>

      do i = 1,N-1
         do j = 1,N-1
            Minv(i,j) = 0.
         end do
         Minv(i,i) = 1.
         d(i) = 1.    
      end do
      
      do i = 1,N-1
         do j = N-1,1,-1
            s = Omega(i,j)
            do k = N-1,j+1,-1
               s = s - Omega(i,k)*Minv(k,j)
            end do
            Omega(i,j) = s/d(j)
         end do
      end do
      


C ==>  << Compute Omega = sqrt(N)*Omega >>

      r = sqrt(real(N))
      
      do i = 1,N-1
         do j = 1,N-1
            Omega(i,j) = r*Omega(i,j)
         end do
      end do



C ==>  << Generate initial Ens >>

      do j = 1,dimstate
         do i = 1,N-1
            v(i) = Ens(j,i)
         end do
         ss = 0.0
         do i = 1,N-1
            s = 0.0
            do k = 1,N-1
               s = s + Omega(i,k)*v(k)
            end do
            Ens(j,i) = stateA(j) + s
            ss = ss + s
         end do
         Ens(j,N) = stateA(j) - ss
      end do
      
      

C ==>  << Generate initial Ens >>

      do i = 1, N
         write(3,rec=i) (Ens(j,i), j = 1, dimstate)
      end do
      
      END




C      **********************************************************************
C      *                                                                    *
C      *                      <<  PROCEDURE RANDMAT  >>                     *
C      *                                                                    *
C      **********************************************************************



C ==>  Cette procedure genere une matrice aleatoire uniforme orthogonale  
C      Womega de dim m*m de la maniere suivante:


C      1- Intialisation:
C      -----------------

C       On part de Womega(1) une matrice 1*1 (un reel) = + ou - 1 avec p = 0.5,
C       il suffit pour cela de tirer un reel x selon la loi uniforme sur [0,1],
C       alors si x < 0.5 on prend Womega(1) = 1 et si x > 0.5, Womega(1) = -1.


C      2- Iteration:  <<  k = 1 --> m >>
C      -------------

C       On genere un vecteur Zk = (z1,...,zk) suivant la loi uniforme sur 
C       la boule unite de Rk. Ceci peut se faire en generant k reels y1,.
C       ..,yk selon la loi normale cetree reduite et puis de prendre 
C                              zi = yi/sqrt(sum(yi*yi))

C       Alors la matrice k*k:
C                         
C                     Womega(k) = [ H(Zk)*Womega(k-1) | Zk ]
C 
C       est une matrice aleatoire uniforme orthogonale de dim k*k. Avec 
C       H(z) la matrice de householder associee au vecteur z.



C ==>  Parametres d entree   =  m     : taille de la matrice Womega a generer;      
C                               ihaz  : entier a choisir:

C ==>  Parametres de sortie  =  Womega : matrice aleatoire uniforme orthogonal de dim m*m;




C ****************<<  Generation d un nombre aleatoire uniforme >>***************    

 
      SUBROUTINE rann(idum,ndum,ran)


      INTEGER   idum, IA, IM, IQ, IR, ndum,k
      PARAMETER (IA=16807,IM=2147483647,IQ=127773,IR=2836,AM=1./IM)
      REAL*8    ran

      k = idum/IQ
      
      ndum = IA*(idum - k*IQ) - IR*k
      if (ndum .lt. 0) ndum = ndum + IM
      ran = AM*ndum

      
      return

      end


C   --------------------------------------------------------------------------------
      

      SUBROUTINE gaussien(idum,ndum,gaussd)


      INTEGER  idum, ndum

      INTEGER  iset
      REAL*8   v1, v2, fac, rsq, gset, x, gaussd

      SAVE     iset, gset
      DATA     iset /0/
      
      
      ndum = idum
      
      if (iset .eq. 0) then

 10      call rann(ndum,ndum,x)
         v1 = 2.0*x - 1.0
         call rann(ndum,ndum,x)
         v2 = 2.0*x - 1.0
         rsq = v1*v1 + v2*v2                    
         if ((rsq .ge. 1.0).or.(rsq .eq. 0)) goto 10     
         fac    = sqrt(- 2.0*log(rsq)/rsq)         
         iset   = 1         
         gset   = v2*fac
         gaussd = v1*fac
      
      else
      
         iset   = 0
         gaussd = gset
      
      end if


      return

      end


C   ---------------------------------------------------------------------------------
     

      SUBROUTINE  randmat(Womega,m,ihaz)



C ==>  << Variables locales >>


      INTEGER   ihaz, m, k, i, np, mp
      REAL*8    Womega(m,m), x 
      REAL*8    r, s



C ==>  << Initialisation: Calcul de Omega(1,1)>>


      call rann(ihaz,np,x)

      if (x .gt. 0.5) then
         Womega(1,1) = 1.0
      else
         Womega(1,1) = -1.0
      end if
      


C ==>  << Construction iteratives de Omega>>


      do k = 2,m


C -->  Generation d un vecteur aleatoire uniforme sur la sphere unite;


         s = 0.0
        
         do i = 1,k
            call gaussien(np,mp,x)
            np = mp
            Womega(k,i) = x
            r           = Womega(k,i)
            s           = s + r**2
         end do     

         s = sqrt(s)
         
         do i = 1,k
            Womega(i,k) = Womega(k,i)/s
         end do
         

C -->  Calcul de la nouvelle matrice Omega;


         r = 1 - Womega(k,k)
         
         do j = 1,k-1
            
            s = 0.0            
            do i = 1,k-1               
               s = s + Womega(i,k)*Womega(i,j) 
            end do
            
            Womega(k,j) = s
            
            s = s/r           
            do i = 1,k-1
               Womega(i,j) = Womega(i,j) - Womega(i,k)*s
            end do
            
         end do

      end do

      ihaz = mp


      return

      END

