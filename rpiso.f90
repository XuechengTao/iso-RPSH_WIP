!##########################################################
!NECESSARY SUBROUTINE FOR ISO-RPSH SIMULATION,
!IF YOU ARE DOING A RPSH-NOKINKS SIMULATION, 
!USE THE SUBROUTINE IN rpshnokinks.f90
!##########################################################
!
!IMPORTANT REMINDER:
![1] flag_nac: CURRENTLY NOT IN USE.
![2] FOR THE SUBROUTINE twoleveliso, MAKE SURE THAT surindex1 IS LOWER THAN surindex2.
![3] DEFAULT MEMORY SPACE FOR THE DIAGONALIZATION SUBROUTINE SUPPORTS nsurfaces LOWER THAN 6, CHECK INTEL-MKL LIBRARY FOR DETAILS REGARDING
!    SUBROUTINE "zheevd" (rpsh.f90), "dsyevd" (rpiso.f90).
!
!NEED TO BE DISCUSSED:
![4] ELECTRONIC WAVEFUNCTIONS / NON-ADIABATIC COUPLING REMIANS UNDETERMINED WITHIN A PHASE FACTOR; 
!CURRENTLY APPLY THE GAUGE: UNIDTA(isur,isur) IS A REAL POSITIVE NUMBER.
!
!##########################################################

subroutine PotentialRead(dimensionality,nsurfaces,natoms,nbeads,flag_nac)
    use potential
    IMPLICIT NONE
    integer,intent(in)::dimensionality,nsurfaces,natoms,nbeads,flag_nac
    real(8)::Viisoa,Vjisoa,Viisod,Vjisod,Kijisod
    real(8)::V(nsurfaces),dV(dimensionality,natoms,nsurfaces)
    real(8),allocatable::dViisoa(:,:),dVjisoa(:,:),dc1isoa(:,:),dViisod(:,:,:),dVjisod(:,:),dKijisod(:,:)
    real(8)::delta,ddelta
    integer::idim,isur,jsur,ksur,iatom,ibead
    real(8),allocatable::eigenVa(:),eigendVa(:,:),constructVd(:,:),diagVd(:,:),unidta(:,:),constructdVd(:,:,:),dVmanybody(:),dmu(:),dmutwobody(:)
    real(8)::work(100),mutwobody,Vmanybody,mutotal,transformationangle,realwork
    integer::iwork(50),info
    
    if (flag_inputrep.eq.0) then
        write(*,*) "THE METHODOLOGY IS UNDER CONSTRUCTION CURRENTLY"
        write(*,*) "FOR AN ab initio CALCULATIONS, "
        write(*,*) "PLEASE FIRST DERIVE THE DIABATIC BASIS AS THE INPUT"
        STOP
        
    else if (flag_inputrep.eq.1) then
        
        if (nsurfaces.eq.1) then

            Viisod=0.0d0
            do ibead=1,nbeads,1
                call PotentialInput(dimensionality,nsurfaces,natoms,ibead,V,dV)
                Viisod=Viisod+V(1)
                if (surface.eq.1) then
                    do idim=1,dimensionality,1
                    do iatom=1,natoms,1
                        derisurface(idim,iatom,ibead)=dV(idim,iatom,1)
                    end do
                    end do
                end if
            end do
        
            potsurface(1)=Viisod/dble(nbeads)

            derivativecoupling=0.0d0
            transformationmatrix(1,1)=1.0d0
            
        else if (nsurfaces.eq.2) then
                
            allocate(dViisoa(dimensionality,nbeads),dVjisoa(dimensionality,nbeads),dc1isoa(dimensionality,nbeads))
            allocate(dViisod(dimensionality,natoms,nbeads),dVjisod(dimensionality,nbeads),dKijisod(dimensionality,nbeads))
            allocate(constructVd(nsurfaces,nsurfaces),unidta(nsurfaces,nsurfaces),constructdVd(nsurfaces,nsurfaces,nbeads))
            allocate(dmu(nbeads),dVmanybody(nbeads),eigenVa(nsurfaces),eigendVa(nsurfaces,nbeads),dmutwobody(nbeads))
   
            call twoleveliso(dimensionality,nsurfaces,nbeads,1,2,Viisoa,Vjisoa,dViisoa,dVjisoa,dc1isoa,Viisod,Vjisod,dViisod,dVjisod,Kijisod,dKijisod)            

            potsurface=0.0d0
            potsurface(1)=Viisoa
            potsurface(2)=Vjisoa
            
            derisurface=0.0d0
            if (surface.eq.1) then
                do ibead=1,nbeads,1
                derisurface(1,ibead)=dViisoa(1,ibead)
                end do
            else if (surface.eq.2) then
                do ibead=1,nbeads,1
                derisurface(1,ibead)=dVjisoa(1,ibead)
                end do
            end if

            derivativecoupling=0.0d0
            do ibead=1,nbeads,1
                derivativecoupling(1,2,1,ibead)=dc1isoa(1,ibead)
                derivativecoupling(2,1,1,ibead)=-derivativecoupling(1,2,1,ibead)
            end do
           
            transformationangle=0.5d0*Datan2(2.0d0*Kijisod,Vjisod-Viisod)
            transformationmatrix(1,1)=Dcos(transformationangle)
            transformationmatrix(1,2)=Dsin(transformationangle)
            transformationmatrix(2,1)=-transformationmatrix(1,2)
            transformationmatrix(2,2)=transformationmatrix(1,1)
                    
            deallocate(dViisoa,dVjisoa,dc1isoa,dViisod,dVjisod,dKijisod)
            deallocate(constructVd,constructdVd,unidta,eigenVa,dVmanybody,dmu,eigendVa,dmutwobody)

        else
            
            potsurface=0.0d0
            derisurface=0.0d0
            derivativecoupling=0.0d0
            
            do isur=1,nsurfaces,1
            do jsur=isur+1,nsurfaces,1
                call twoleveliso(dimensionality,nsurfaces,nbeads,isur,jsur,Viisoa,Vjisoa,dViisoa,dVjisoa,dc1isoa,Viisod,Vjisod,dViisod,dVjisod,Kijisod,dKijisod)            
              
                if (isur.eq.1) then
                    if (jsur.eq.2) then
                        constructVd(isur,isur)=Viisod
                        constructVd(jsur,jsur)=Vjisod
                    else
                        constructVd(jsur,jsur)=Vjisod
                    end if
                end if
                constructVd(isur,jsur)=Kijisod
                constructVd(jsur,isur)=constructVd(isur,jsur)
                
                do ibead=1,nbeads,1
                    constructdVd(isur,isur,ibead)=dViisod(1,ibead)
                    constructdVd(jsur,jsur,ibead)=dVjisod(1,ibead)
                    constructdVd(isur,jsur,ibead)=dKijisod(1,ibead)
                    constructdVd(jsur,isur,ibead)=constructdVd(isur,jsur,ibead)
                end do
            end do
            end do
            
            unidta=constructVd
            call dsyevd('V','U',nsurfaces,unidta,nsurfaces,potsurface,work,100,iwork,50,info)
            do isur=1,nsurfaces,1
                realwork=0.0d0
                do jsur=1,nsurfaces,1
                realwork=realwork+unidta(jsur,isur)*transformationmatrix(jsur,isur)
                end do
                if (realwork.lt.0.0d0) then
                do jsur=1,nsurfaces,1
                unidta(jsur,isur)=-unidta(jsur,isur)
                end do
                end if
            end do
            transformationmatrix=unidta
            
            do ibead=1,nbeads,1
                constructdVd(:,:,ibead)=matmul(constructdVd(:,:,ibead),unidta)
                constructdVd(:,:,ibead)=matmul(transpose(unidta),constructdVd(:,:,ibead))
                do isur=1,nsurfaces,1
                    eigendVa(isur,ibead)=constructdVd(isur,isur,ibead)
                    constructdVd(isur,isur,ibead)=0.0d0
                do jsur=isur+1,nsurfaces,1
                    constructdVd(isur,jsur,ibead)=constructdVd(isur,jsur,ibead)/(potsurface(isur)-potsurface(jsur))
                    constructdVd(jsur,isur,ibead)=-constructdVd(isur,jsur,ibead)
                end do
                end do
                
                derisurface(1,ibead)=eigendVa(surface,ibead)
                
                do isur=1,nsurfaces,1
                do jsur=isur+1,nsurfaces,1
                    derivativecoupling(isur,jsur,1,ibead)=constructdVd(jsur,isur,ibead)
                    derivativecoupling(jsur,isur,1,ibead)=-derivativecoupling(isur,jsur,1,ibead)
                end do 
                end do
            end do
        
            mutwobody=0.0d0
            do isur=1,nsurfaces,1
                mutwobody=mutwobody+Dexp(-inversetemp*potsurface(isur))
            end do
            
            dmutwobody=0.0d0
            do ibead=1,nbeads,1
            do isur=1,nsurfaces,1
                dmutwobody(ibead)=dmutwobody(ibead)-betan*Dexp(-inversetemp*potsurface(isur))*eigendVa(isur,ibead)
            end do
            end do
            
            call mucalc(dimensionality,nsurfaces,nbeads,mutotal,dmu)
            Vmanybody=-Dlog(mutotal/mutwobody)/inversetemp
            potsurface=potsurface+Vmanybody
            
            do ibead=1,nbeads,1
                dVmanybody(ibead)=-dmu(ibead)/inversetemp*mutwobody/mutotal
                dVmanybody(ibead)=dVmanybody(ibead)+dmutwobody(ibead)/inversetemp
                dVmanybody(ibead)=dVmanybody(ibead)/mutwobody*dble(nbeads)
                derisurface(1,ibead)=derisurface(1,ibead)+dVmanybody(ibead)
            end do
            
        end if
        
    else
        STOP "(WRONG flag_inputrep)"
    end if
    
end subroutine

subroutine mucalc(dimensionality,nsurfaces,nbeads,mutotal,dmu)
    use potential
    IMPLICIT NONE
    integer,intent(in)::dimensionality,nsurfaces,nbeads
    real(8),intent(out)::mutotal,dmu(nbeads)
    
    real(8)::V(nsurfaces),dV(nsurfaces),K(nsurfaces,nsurfaces)
    real(8)::constructVd(nsurfaces,nsurfaces),constructdVd(nsurfaces,nsurfaces),unidta(nsurfaces,nsurfaces),dunidta(nsurfaces,nsurfaces)
    real(8)::eigenVa(nsurfaces),eigendVa(nsurfaces),Vaexp(nsurfaces,nsurfaces),dVaexp(nsurfaces,nsurfaces)
    real(8)::Mexp(nsurfaces,nsurfaces,nbeads),dMexp(nsurfaces,nsurfaces,nbeads),holel(nsurfaces,nsurfaces,nbeads),holer(nsurfaces,nsurfaces,nbeads)
    integer::isur,jsur,ibead
    real(8)::work(100),realwork
    integer::iwork(50),info
    
    do ibead=1,nbeads,1
        call PotentialInput(dimensionality,nsurfaces,ibead,V,dV)
        call CouplingInput(dimensionality,nsurfaces,ibead,K)
        
        constructVd=0.0d0
        constructdVd=0.0d0
        do isur=1,nsurfaces,1
            constructVd(isur,isur)=V(isur)
            constructdVd(isur,isur)=dV(isur)
        do jsur=isur+1,nsurfaces,1
            constructVd(isur,jsur)=K(isur,jsur)
            constructVd(jsur,isur)=constructVd(isur,jsur)
            constructdVd(isur,jsur)=K(jsur,isur)
            constructdVd(jsur,isur)=constructdVd(isur,jsur)
        end do
        end do
        unidta=constructVd
        call dsyevd('V','U',nsurfaces,unidta,nsurfaces,eigenVa,work,100,iwork,50,info)
        do isur=1,nsurfaces,1
            realwork=0.0d0
            do jsur=1,nsurfaces,1
            realwork=realwork+unidta(jsur,isur)*beadtransformationmatrix(jsur,isur,ibead)
            end do
            if (realwork.lt.0.0d0) then
                do jsur=1,nsurfaces,1
                unidta(jsur,isur)=-unidta(jsur,isur)
                end do
            end if
        end do
        beadtransformationmatrix(:,:,ibead)=unidta
        
        dunidta=matmul(transpose(unidta),matmul(constructdVd,unidta))
        do isur=1,nsurfaces,1
            eigendVa(isur)=dunidta(isur,isur)
            dunidta(isur,isur)=0.0d0
        do jsur=isur+1,nsurfaces,1
            dunidta(isur,jsur)=dunidta(isur,jsur)/(eigenVa(isur)-eigenVa(jsur))
            dunidta(jsur,isur)=-dunidta(isur,jsur)
        end do
        end do
        dunidta=matmul(dunidta,transpose(unidta))
        dunidta=transpose(dunidta)
        
        Vaexp=0.0d0
        dVaexp=0.0d0
        do isur=1,nsurfaces,1
            Vaexp(isur,isur)=Dexp(-betan*eigenVa(isur))
            dVaexp(isur,isur)=-betan*eigendVa(isur)*Vaexp(isur,isur)
        end do
        
        dMexp(:,:,ibead)=matmul(matmul(unidta,dVaexp),transpose(unidta))
        Mexp(:,:,ibead)=matmul(matmul(dunidta,Vaexp),transpose(unidta))
        dMexp(:,:,ibead)=dMexp(:,:,ibead)+Mexp(:,:,ibead)
        Mexp(:,:,ibead)=matmul(unidta,Vaexp)
        dMexp(:,:,ibead)=dMexp(:,:,ibead)+matmul(Mexp(:,:,ibead),transpose(dunidta))
        Mexp(:,:,ibead)=matmul(Mexp(:,:,ibead),transpose(unidta))
    end do
    
    holel(:,:,1)=Mexp(:,:,1)
    do ibead=2,nbeads-1,1
        holel(:,:,ibead)=matmul(holel(:,:,ibead-1),Mexp(:,:,ibead))
    end do
    holer(:,:,nbeads)=Mexp(:,:,nbeads)
    do ibead=nbeads-1,2,-1
        holer(:,:,ibead)=matmul(Mexp(:,:,ibead),holer(:,:,ibead+1))
    end do
    if (nbeads.gt.1) then
        mutotal=0.0d0
        constructVd=matmul(holel(:,:,nbeads-1),Mexp(:,:,nbeads))
        do isur=1,nsurfaces,1
            mutotal=mutotal+constructVd(isur,isur)
        end do
    else
        mutotal=0.0d0
        do isur=1,nsurfaces,1
            mutotal=mutotal+holel(isur,isur,1)
        end do
    end if
    
    if (nbeads.gt.1) then
        constructVd=matmul(dMexp(:,:,1),holer(:,:,2))
        dmu(1)=0.0d0
        do isur=1,nsurfaces,1
            dmu(1)=dmu(1)+constructVd(isur,isur)
        end do
        do ibead=2,nbeads-1,1
            constructVd=matmul(holel(:,:,ibead-1),dMexp(:,:,ibead))
            constructVd=matmul(constructVd,holer(:,:,ibead+1))
            dmu(ibead)=0.0d0
            do isur=1,nsurfaces,1
                dmu(ibead)=dmu(ibead)+constructVd(isur,isur)
            end do
        end do
        constructVd=matmul(holel(:,:,nbeads-1),dMexp(:,:,nbeads))
        dmu(nbeads)=0.0d0
        do isur=1,nsurfaces,1
            dmu(nbeads)=dmu(nbeads)+constructVd(isur,isur)
        end do
    else
        dmu(1)=0.0d0
        do isur=1,nsurfaces,1
            dmu(1)=dmu(1)+dMexp(isur,isur,1)
        end do
    end if
    
end subroutine

subroutine twoleveliso(dimensionality,nsurfaces,nbeads,surindex1,surindex2,Viisoa,Vjisoa,dViisoa,dVjisoa,dc1isoa,Viisod,Vjisod,dViisod,dVjisod,Kijisod,dKijisod)
    use potential
    IMPLICIT NONE
    integer,intent(in)::dimensionality,nsurfaces,nbeads,surindex1,surindex2
    real(8),intent(out)::Viisoa,Vjisoa,dViisoa(dimensionality,nbeads),dVjisoa(dimensionality,nbeads),dc1isoa(dimensionality,nbeads)
    real(8),intent(out)::Viisod,Vjisod,Kijisod,dViisod(dimensionality,nbeads),dVjisod(dimensionality,nbeads),dKijisod(dimensionality,nbeads)
        
    real(8)::V(nsurfaces),dV(nsurfaces),K(nsurfaces,nsurfaces)
    real(8)::delta,ddelta,Va(2),dVa(2),theta,dtheta,mu,realmatrixtemp(2,2),expsur,realtemp,root
    real(8),allocatable::Mexp(:,:,:),dMexp(:,:,:),holel(:,:,:),holer(:,:,:),dmu(:),droot(:)
    real(8)::unidta(2,2),dunidta(2,2),Vaexp(2,2),dVaexp(2,2)
    integer::idim,isur,ibead

    allocate(Mexp(2,2,nbeads),dMexp(2,2,nbeads),holel(2,2,nbeads),holer(2,2,nbeads),dmu(nbeads),droot(nbeads))
    
    Viisod=0.0d0
    Vjisod=0.0d0
    do ibead=1,nbeads,1
        call PotentialInput(dimensionality,nsurfaces,ibead,V,dV)
        call CouplingInput(dimensionality,nsurfaces,ibead,K)
        
        Viisod=Viisod+V(surindex1)
        Vjisod=Vjisod+V(surindex2)
        dViisod(1,ibead)=dV(surindex1)
        dVjisod(1,ibead)=dV(surindex2)
        
        delta=(V(surindex1)-V(surindex2))**2+4.0d0*K(surindex1,surindex2)**2
        delta=Dsqrt(delta)*0.5d0
        Va(1)=(V(surindex1)+V(surindex2))*0.5d0-delta
        Va(2)=(V(surindex1)+V(surindex2))*0.5d0+delta
        ddelta=(V(surindex1)-V(surindex2))*(dV(surindex1)-dV(surindex2))*0.5d0+2.0d0*K(surindex1,surindex2)*K(surindex2,surindex1)
        ddelta=ddelta*0.5d0/delta
        dVa(1)=(dV(surindex1)+dV(surindex2))*0.5d0-ddelta
        dVa(2)=(dV(surindex1)+dV(surindex2))*0.5d0+ddelta
        theta=0.5d0*Datan2(2.0d0*K(surindex1,surindex2),V(surindex2)-V(surindex1))
        dtheta=(K(surindex2,surindex1)*(V(surindex2)-V(surindex1))-K(surindex1,surindex2)*(dV(surindex2)-dV(surindex1)))
        dtheta=dtheta/(delta*2.0d0)**2
        
        unidta(1,1)=Dcos(theta)
        unidta(1,2)=Dsin(theta)
        unidta(2,1)=-unidta(1,2)
        unidta(2,2)=unidta(1,1)
        dunidta(1,1)=-unidta(1,2)*dtheta
        dunidta(1,2)=unidta(1,1)*dtheta
        dunidta(2,1)=-dunidta(1,2)
        dunidta(2,2)=dunidta(1,1)
        
        Vaexp=0.0d0
        Vaexp(1,1)=Dexp(-betan*Va(1))
        Vaexp(2,2)=Dexp(-betan*Va(2))
        dVaexp=0.0d0
        dVaexp(1,1)=-betan*dVa(1)*Vaexp(1,1)
        dVaexp(2,2)=-betan*dVa(2)*Vaexp(2,2)
        
        dMexp(:,:,ibead)=matmul(matmul(unidta,dVaexp),transpose(unidta))
        Mexp(:,:,ibead)=matmul(matmul(dunidta,Vaexp),transpose(unidta))
        dMexp(:,:,ibead)=dMexp(:,:,ibead)+Mexp(:,:,ibead)
        Mexp(:,:,ibead)=matmul(unidta,Vaexp)
        dMexp(:,:,ibead)=dMexp(:,:,ibead)+matmul(Mexp(:,:,ibead),transpose(dunidta))
        Mexp(:,:,ibead)=matmul(Mexp(:,:,ibead),transpose(unidta))
    end do
    Viisod=Viisod/dble(nbeads)
    Vjisod=Vjisod/dble(nbeads)
    
    holel(:,:,1)=Mexp(:,:,1)
    do ibead=2,nbeads-1,1
        holel(:,:,ibead)=matmul(holel(:,:,ibead-1),Mexp(:,:,ibead))
    end do
    holer(:,:,nbeads)=Mexp(:,:,nbeads)
    do ibead=nbeads-1,2,-1
        holer(:,:,ibead)=matmul(Mexp(:,:,ibead),holer(:,:,ibead+1))
    end do
    if (nbeads.gt.1) then
        realmatrixtemp=matmul(holel(:,:,nbeads-1),Mexp(:,:,nbeads))
        mu=realmatrixtemp(1,1)+realmatrixtemp(2,2)
    else
        mu=holel(1,1,1)+holel(2,2,1)
    end if
    
    if (nbeads.gt.1) then
        realmatrixtemp=matmul(dMexp(:,:,1),holer(:,:,2))
        dmu(1)=realmatrixtemp(1,1)+realmatrixtemp(2,2)
        do ibead=2,nbeads-1,1
            realmatrixtemp=matmul(holel(:,:,ibead-1),dMexp(:,:,ibead))
            realmatrixtemp=matmul(realmatrixtemp,holer(:,:,ibead+1))
			dmu(ibead)=realmatrixtemp(1,1)+realmatrixtemp(2,2)
        end do
        realmatrixtemp=matmul(holel(:,:,nbeads-1),dMexp(:,:,nbeads))
        dmu(nbeads)=realmatrixtemp(1,1)+realmatrixtemp(2,2)
    else
        dmu(1)=dMexp(1,1,1)+dMexp(2,2,1)
    end if
    
    expsur=Dexp(0.5d0*inversetemp*(Viisod+Vjisod))
    root=expsur*mu*0.5d0
    root=Dacosh(root)/inversetemp
    
    Viisoa=0.0d0
    Vjisoa=0.0d0
    if (nsurfaces.eq.2) then
        realtemp=Viisod+Vjisod
        Viisoa=realtemp*0.5d0-root
        Vjisoa=realtemp*0.5d0+root
    end if
    
    dViisoa=0.0d0
    dVjisoa=0.0d0    
    do ibead=1,nbeads,1
        realtemp=dViisod(1,ibead)+dVjisod(1,ibead)
        droot(ibead)=0.25d0*betan*realtemp*mu+0.5d0*dmu(ibead)
        droot(ibead)=droot(ibead)*expsur/betan
        droot(ibead)=droot(ibead)/Dsinh(root*inversetemp)
        if (nsurfaces.eq.2) then
            dViisoa(1,ibead)=(dViisod(1,ibead)+dVjisod(1,ibead))*0.5d0-droot(ibead)
            dVjisoa(1,ibead)=(dViisod(1,ibead)+dVjisod(1,ibead))*0.5d0+droot(ibead)
        end if
    end do
                
    Kijisod=root**2-0.25d0*(Vjisod-Viisod)**2
    if (Kijisod.lt.precision) then
        if (Kijisod.lt.-1.0d-8) then
        write(*,*) "SIGNIFICANT ISOMORPHIC COUPLING NEGLECTED, (Kij^{iso}_d)^2 = ", Kijisod
        end if
        Kijisod=0.0d0
        dKijisod=0.0d0
        dc1isoa=0.0d0
    else
        Kijisod=Dsqrt(Kijisod)
        do ibead=1,nbeads,1
            dKijisod(1,ibead)=root*droot(ibead)-0.25d0*(Vjisod-Viisod)*(dVjisod(1,ibead)-dViisod(1,ibead))
            dKijisod(1,ibead)=dKijisod(1,ibead)/Kijisod
            if (nsurfaces.eq.2) then            
                dc1isoa(1,ibead)=(Vjisod-Viisod)*dKijisod(1,ibead)-Kijisod*(dVjisod(1,ibead)-dViisod(1,ibead))
                dc1isoa(1,ibead)=0.25d0*dc1isoa(1,ibead)/root**2
            end if
        end do  
    end if
    
    deallocate(Mexp,dMexp,holel,holer,dmu,droot)
    
end subroutine
