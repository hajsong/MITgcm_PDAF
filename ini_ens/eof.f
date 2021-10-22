C     ***********************************************************************
C     *                                                                     *
C     *                        PROGRRAMME PRINCIPAL                         *
C     *                                                                     *
C     ***********************************************************************


C ==>  << Definition des variables >>

C -->   nvars    = nombre des variables dans le vecteur detat;
cbdc               number of variable types in the state vector (T,S,U,V,SSH,.)
C -->   dimstate = dimension des etats;
cbdc               total dimension of state variable (total length)
C -->   nstates  = nombre des etats;
cbdc               number of states to use in making the eofs
C -->   neigen   = nbre des vp a calculer;
cbdc               number of eigenvectors to calculate
C -->   prec     = precision pou le calcul des vp;
cbdc               precision for calculating the eigenvectors
C -->   miter    = nombre maximum d iteration pour le calcul des vp;
cbdc               max number of iterations for calculating the eigenvectors

      INTEGER    nvars, dimstate, nstates, miter, neigen
      INTEGER    dims, dimt, dimu, dimv, dimh

      REAL*4     prec

      PARAMETER (nvars = 5, dimstate = 107376, 
     &           dims = 27405, dimt = 27405, 
     &           dimu = 25459, dimv = 24792, dimh = 2315, 
     &           nstates = 120, neigen = 10, 
     &           prec = 1e-7, miter = 10000)
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C     &           prec = 1e-8, miter = 10000) 
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         
C -->   etat    = vecteur d etat du fichier d entree;
cbdc              state vector as read in for one record
C -->   etat1   = variable de calcul;
cbdc              work array
C -->   etatmoy = etat barycentre des etats;
cbdc              mean state
C -->   etatvar = etat variance des variables;
cbdc              variance of each state variable

      REAL*4    etat(dimstate), etat1(dimstate), 
     &           etatmoy(dimstate), etatvar(dimstate),
     &           etat2(dimstate)

C -->  covmat est un vecteur contenant la matrice de covariance des etats; 

      REAL*4    covmat(nstates*(nstates+1)/2) 


C -->  eigvec est un tableau contenant les vecteurs propres de covmat; 
C -->  eigval est un vecteur contenant les valeurs propres de covmat; 

      REAL*4        eigvec(nstates,neigen), eigval(neigen), 
     &              vec2  (nstates)       , vecr  (nstates) 

      INTEGER       idum, dimvar(5)
cbdc         working variables
      REAL*8        s, tmp, etats(dimstate) 

      
C ==>  << Ouverture des fichiers Input/Ouputs >>  

      open(9,file='eof_vecs.bin',
     &  access = 'direct', form = 'unformatted', recl = 4*dimstate)  
      open(1,file='eof_basis.bin',
     &  access = 'direct', form = 'unformatted', recl = 4*dimstate) 
      open(2,file='mean_state.bin',
     &  access = 'direct', form = 'unformatted', recl = 4*dimstate) 
      open(3,file='eof_eig.bin',  
     &  access = 'direct', form = 'unformatted', recl = 4*1)
      open(4,file='var_state.bin',
     &  access = 'direct', form = 'unformatted', recl = 4*dimstate) 
      open(5,file='var_temp.bin',
     &  access = 'direct', form = 'unformatted', recl = 4*dimstate) 
      open(7,file='out_eof.txt',  
     &  form = 'formatted')


      print*,"" 
      print*,"nbre des etats =", nstates 
      print *, "" 
      

C ==>  << Compute the mean state >>

cbdc initialize mean
      do k = 1, dimstate
         etats(k) = 0.
      end do
   
cbdc read each record and sum
      do nt = 1, nstates
         read(9,rec = nt) etat
         do k = 1, dimstate
            etats(k) = etats(k) + etat(k)
         end do
      end do
      
      do k = 1, dimstate
         etatmoy(k) = etats(k)/nstates
      end do

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cbdc mean state
      write(2,rec=1) etatmoy
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cbdc now etats() is the mean over all nstates realizations

C ==>  << Compute the spatial mean of the variance for each physical variable >>

c dimension of each variable type
      dimvar(1) = dims
      dimvar(2) = dimt
      dimvar(3) = dimu
      dimvar(4) = dimv
      dimvar(5) = dimh

C -->  Compute the variance of each statistical variable;

cbdc initialize variance for each state vector element
      do k = 1, dimstate
         etats(k) = 0.
      end do
      
cbdc read state again, subtract mean, accumulate variance
      do nt = 1, nstates
         read(9,rec = nt) etat
         do k = 1, dimstate
            etats(k) = etats(k) + (etat(k) - etatmoy(k))**2
         end do
      end do

cbdc divide by number of realizations
      do k = 1, dimstate
         etats(k) = etats(k)/nstates
         etat2(k) = etats(k)
      end do

CGG Time mean variance for each point for x,y,z,var)           
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      write(5,rec=1) etat2 
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C -->  Compute the spatial mean of the variance of each physical variable;


         nv = 0

      do iv = 1, nvars
cbdc offset in state vector
C         nv = 0
cbdc total variance over all of this variable type

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         print *,  "dim(", iv, ")", dimvar(iv),nv
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         s = 0.
         do j = 1, dimvar(iv)               
            s = s + etats(nv+j)
         end do
cbdc divide by number of elements of this type to get a mean.
         s = s/dimvar(iv)
cbdc make etatvar(s) have this value for all elements of this type
         do j = 1, dimvar(iv)     
            etatvar(nv+j) = s
         end do
cbdc update offset in state vector
         nv = nv + dimvar(iv)
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         print *,  "dim(", iv, ")", dimvar(iv),nv
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      end do

CGG now etatvar has normalization for each variable

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cbdc variance
      write(4,rec=1) etatvar
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>


C ==>  << Compute the sample covariance matrix of the state vectors >>

C -->  To reduce cost;

CGG reuse etats loops over statevector
      do j = 1,dimstate
         if (etatvar(j).eq.0.) then
cbdc land points are zero
            etats(j) = 0.
         else
cbdc 1/(spatial RMS)
            etats(j) = 1./sqrt(etatvar(j))
         end if
         etat2(j) = etats(j)
      end do

cbdc this will be used to normalize each vector by the spatial (x,y,z)
cbdc average of its variance.
cbdc This will downweight regions with small variance.
cbdc alternatives could be: A) average in x,y, so have different variance
cbdc   for each z, B) average point by point, so each point is equally
cbdc   important, C) read in an externally-chosen normalization for
cbdc   the state vector.

C -->  "etat" contains the last state;
cbdc                            read
cbdc number of unique covariance elements in the realization-realization
cbdc  (usually time-time) covariance
      ii = nstates*(nstates+1)/2
      trace = 0.

cbdc go backwards over states, so that don't have to re-read last one
      do i = nstates, 1, -1
         s = 0.
         do k = 1, dimstate
cbdc subtract mean and divide by spatial RMS
            tmp      = (etat(k) - etatmoy(k))*etats(k)
cbdc save the normalized anomaly
            etat1(k) = tmp
cbdc sum the variance (dot product over all state elements)
            s        = s + tmp**2
         end do

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        tmp         = s/nstates
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cbdc diagonal element
         covmat(ii)  = tmp
cbdc add to trace
         trace       = trace + tmp
         print *,  "Moyenne quadratique etat(", i, ")", tmp
         write(7,*)"Moyenne quadratique etat(", i, ")", tmp
cbdc now do points off the diagonal
cbdc first point in this row
         ii = ii - i
         do j = 1,i-1
cbdc read in the other state
            read(9,rec=j) etat       
            s = 0.
            do k = 1,dimstate
cbdc compute the covariance with the normalized anomaly in etat1()
               s = s + ((etat(k) - etatmoy(k))*etats(k))*etat1(k)
            end do
cbdc off-diagonal covariance

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            covmat(ii+j) = s/nstates
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         end do
      end do

      print *,  ""
      print *,  "Somme des moyennes quadratique ", trace
      print *,  ""

      write(7,*)
      write(7,*)"Somme des moyennes quadratique ", trace
      write(7,*)
       

C ==>  << Computing of the neigen first eigenvectors and eigenvalues of the covariance matrix >>

      idum = 1
      tmp  = 0.

      do ieig = 1, neigen
         do k = 1, nstates
            eigvec(k,ieig) = grand(idum)
         end do
cbdc use ibrahim's routine to calculate just this eigenvector
cbdc (in realization space)
cbdc (this is not a big matrix; this could be done in matlab)
         iter = ieigenpower(covmat, eigvec, nstates, ieig,
     &                      eigvec(1,ieig), vec2, vecr, miter-1,
     &                      ieig*prec, eigval, angle)
cbdc sum of eigenvalues (cumulative)
         tmp = tmp + eigval(ieig)
cbdc put out eigenvalues normalized by the trace
         print *,   "V.p.", ieig, eigval(ieig), " (",
     &        eigval(ieig)/trace, tmp/trace, ")", angle, " iter ", iter
         write(7,*) "V.p.", ieig, eigval(ieig), " (",
     &        eigval(ieig)/trace, tmp/trace, ")", angle, " iter ", iter 
         write(3,rec=ieig)  eigval(ieig)

CGG   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         s = 0.
         do j = 1, nstates
            s = s + eigvec(j,ieig)**2
         end do

         if (abs(s - 1.0) .le. 0.001) then
                print*, "GOOD eigvec", s
                write(7,*) "GOOD eigvec=", s
         else
                print*, "CAUTION: BAD eigvec", s
                write(7,*) "CAUTION: BAD eigvec=", s
         end if

CGG   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      end do


C ==>  << Write results >>

cbdc put out eigenvectors
      do i = 1, neigen
         print*,i
         do j = 1, dimstate
            etats(j) = 0.
         end do
cbdc compute the space eigenvector by dot product
cbdc of the time history of states with the time eigenvector
         do k = 1, nstates   
            read(9,rec = k) etat
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            s = (1./sqrt(real(nstates)))*eigvec(k,i)
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            do j = 1, dimstate
CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C            etats(j) = etats(j) + s*(etat(j) - etatmoy(j))

             etats(j) = etats(j) + s*(etat(j) - etatmoy(j))*etat2(j)   

CGG>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
           end do
         end do

cbdc output normalized space eof
C         write(1,rec=i) (real(etats(j)), j = 1, dimstate)



CGG   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         s = 0.
         do j = 1, dimstate
            s = s + etats(j)**2          
         end do
         if (abs(s - eigval(i)) .le. 0.1) then
                print*, "GOOD eigval", s
                write(7,*) "GOOD eigval=", s
         else
                print*, "CAUTION: BAD eigval", s
                write(7,*) "CAUTION: BAD eigval=", s
         end if
CGG   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


cbdc output NON-NORMALIZED space eof

	do j = 1, dimstate
             etats(j) = etats(j)/etat2(j)   
        end do

        write(1,rec=i) (real(etats(j)), j = 1, dimstate)

      end do

      END




C      **********************************************************************
C      *                                                                    *
C      *                        <<  FONCTION RAN0  >>                       *
C      *                                                                    *
C      **********************************************************************



C ==>  << Cette fonction simule la loi uniforme sur [0,1] >> 

    
      FUNCTION ran0(idum)



C ==>  << Definition des variables >>


      INTEGER idum, IA, IM, IQ, IR
      PARAMETER(IA=16807,IM=2147483647,IQ=127773,IR=2836,AM=1./IM)
      INTEGER k



C ==>  << Programme principale >>


      k=idum/IQ
      idum = IA*(idum - k*IQ) - IR*k
      if (idum .lt. 0) idum = idum + IM
      ran0 = AM*idum

      END



C      **********************************************************************
C      *                                                                    *
C      *                        <<  FONCTION GRAND  >>                      *
C      *                                                                    *
C      **********************************************************************



C ==>  << Cette fonction simule la loi normale sur N(0,1) >> 


      FUNCTION grand(idum)



C ==>  << Definition des variables >>


      INTEGER  idum
      INTEGER  iset
      REAL*4   v1, v2, fac, rsq, gset

      SAVE     iset, gset
      DATA     iset /0/


C ==>  << Programme principale >>


      if (iset .eq. 0) then
 10         v1 = 2.0*ran0(idum) - 1.0
            v2 = 2.0*ran0(idum) - 1.0
            rsq = v1*v1 + v2*v2
         if ((rsq .ge. 1.0).or.(rsq .eq. 0)) goto 10
         fac = sqrt(- 2.0*alog(rsq)/rsq)
         iset = 1
         gset = v2*fac
         grand = v1*fac
      else
         iset = 0
         grand = gset
      endif

      end




C      **********************************************************************
C      *                                                                    *
C      *                      <<  PROCEDURE SMATVEC  >>                     *
C      *                                                                    *
C      **********************************************************************



C ==>  Calcule le produit de mat avec vec et mettre le resultat dans prod;
C      mat est une matrice symetrique d'ordre n, ranger dans un vecteur de
C      taille n*(n+1)/2 commee mat11, mat21, mat22, ..., matn1, ..., matnn,
C      vec et prod sont des vecteurs de taille n (doivent etre distincts)
      


      SUBROUTINE smatvec(mat, vec, prod, n)
    


C ==>  << Definition des variables >>


      REAL*4  mat(*), vec(n), prod(n)

      REAL*8  s



C ==>  << Programme principale >>


      ii = 1
      do 10 i = 1,n
         s = 0
         ij = ii
         do 15 j = 1,i-1
            s = s + mat(ij)*vec(j)
            ij = ij + 1
 15      continue
         do 16 j = i, n
            s = s + mat(ij)*vec(j)
            ij = ij + j
 16      continue
         prod(i) = s
         ii = ii + i
 10   continue

      end




C      **********************************************************************
C      *                                                                    *
C      *                    <<  PROCEDURE IEIGENPOWER  >>                   *
C      *                                                                    *
C      **********************************************************************



C -->  Calcule le i-eme vecteur propre de mat par le methode des puissances;
C      mat est une matrix symetrique d'ordre n, ranger dans un vecteur de
C      taille n*(n+1)/2 commee mat11, mat21, mat22, ..., matn1, ..., matnn;
C      eivecp(i) est un vecteur de taile n: vecteur propre initial, eivecp(0),
C      ..., eivecp(i-1) sont suppose deja calcules et normalises; itermax est
C      le nombre maximum d'iteration, eps est le critere d'arret. La valeur et 
C      le vecteur propre normalise correspondant sont retournes dans eigval(i)
C      et eivec, et angle = 1 - cos(angle entre mat^{-1}eivec et eivec); 
 

  
      FUNCTION ieigenpower(mat, eivecp, n, i, eivec, vec2, vecr,
     &                 itermax, eps, eigval, angle)



C ==>  << Definition des variables >>


      REAL*4  mat(*), eivecp(n,i), eivec(n), eigval(i), vec2(n), vecr(n)

      REAL *8  tmp, a1, a2, b2, b, t, t2



C ==> << Othogonaliser eivec par rapport aux autres vecteurs propres (normalises) >>


      do 100 j = 1,i-1
         tmp = 0
         do 110 k = 1,n
            tmp = tmp + eivec(k)*eivecp(k,j)
 110     continue
         do 120 k = 1,n
            eivec(k) = eivec(k) - tmp*eivecp(k,j)
 120     continue
 100  continue
      b2 = 0
      do 200 k = 1,n
         b2 = b2 + eivec(k)**2
 200  continue
      b = dsqrt(b2)
      do 250 k = 1,n
         eivec(k) = eivec(k)/b
 250  continue

      call smatvec(mat, eivec, vecr, n)

      a1 = 0
      do 300 k = 1,n
         a1 = a1 + vecr(k)*eivec(k)
 300  continue
      do 350 k = 1,n
         vecr(k) = vecr(k) - a1*eivec(k)
 350  continue

      iter = 1



C ==>  << Orthogonalisaton de vect par rapport aux  
C              autres vecteurs propres (normalises) >>


 500     do 510 j = 1,i-1
            tmp = 0
            do 520 k = 1,n
               tmp = tmp + vecr(k)*eivecp(k,j)
 520        continue
            do 530 k = 1,n
               vecr(k) = vecr(k) - tmp*eivecp(k,j)
 530        continue
 510     continue
         b2 = 0
         do 540 k = 1,n
            b2 = b2 + vecr(k)**2
 540     continue
         b = dsqrt(b2)

*         print *, "iter ", iter, " a1 ", a1, " b ", b

c     test d'arret
         angle = b/a1
         if (angle .lt. eps) goto 600
         iter = iter + 1
         if (iter .gt. itermax) goto 600

         do 560 k = 1,n
            vec2(k) = vecr(k)/b
 560     continue
         call smatvec(mat, vec2, vecr, n)
         a2 = 0
         do 570 k = 1,n
            a2 = a2 + vec2(k)*vecr(k)
 570     continue

         t = 2*b/(dabs(a1 - a2) + dsqrt((a1 - a2)**2 + 4*b2))
         if ((a1.le.a2) .and. ((a1.ne. a2).or.(a1+a2.lt.0))) t = -t
         t2 = 1 + t**2
         if (a1**2 .gt. a2**2) then
            a1 = (a1 + t**2*a2 + 2*t*b)/t2 
            t2 = dsqrt(t2)
            tmp = t/t2
            do 580 k = 1,n
               vecr(k) = (vecr(k) - a2*vec2(k) - b*eivec(k))*tmp
               eivec(k) = (eivec(k) + t*vec2(k))/t2
 580        continue
         else
            a1 = (a2 + t**2*a1 - 2*t*b)/t2
            t2 = dsqrt(t2)
            do 590 k = 1,n
               vecr(k) = (vecr(k) - a2*vec2(k) - b*eivec(k))/t2
               eivec(k) = (vec2(k) - t*eivec(k))/t2
 590        continue
         endif
      goto 500

 600  tmp = dsqrt(a1*a1 + b2)
      do 700 k = 1,n
         eivec(k) = (a1*eivec(k) + vecr(k))/tmp
 700  continue
      eigval(i) = a1
      ieigenpower = iter

      end


