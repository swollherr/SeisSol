!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stephanie Wollherr
!!
!! @section LICENSE
!! Copyright (c) 2007-2016, SeisSol Group
!! All rights reserved.
!! 
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!! 
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!! 
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!! 
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!! 
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!!
!! @section DESCRIPTION
!! Plasticity module: checks if an element undergoes plastic yielding or not
#include <Initializer/preProcessorMacros.fpp>

MODULE Plasticity_mod
  !---------------------------------------------------------------------------!
  USE TypesDef

#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
#endif
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  !PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Plasticity
     MODULE PROCEDURE Plasticity_3D_high
  END INTERFACE
    INTERFACE Plasticity
     MODULE PROCEDURE Plasticity_3D_DOF
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Plasticity_3D_high
  PUBLIC  :: Plasticity_3D_DOF
  !---------------------------------------------------------------------------!
  CONTAINS


   SUBROUTINE Plasticity_3D_high(DISC, dgvar, DOFStress, nDegFr, nAlignedDegFr, BulkFriction, &
                                 Tv, dt, mu, lambda, parameters , Energy, pstrain, pstrain_gp, &
                                 intGaussP, intGaussW, nVar, nIntGP)


    !-------------------------------------------------------------------------!
  USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization)    :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom
    INTEGER     :: nDegFr
    INTEGER     :: iIntGP, nIntGP
    INTEGER     :: nVar, iPoly
    integer     :: nAlignedDegFr


    REAL        :: stateGP(nVar)                                              ! State in GP
    REAL        :: pstateGP(nVar)                                             ! Primitive state in GP
    REAL        :: Stress_total(1:nIntGP,6)                                   !local stress variable for the yield criterion
    REAL        :: Strain_total(1:nIntGP,6)                                   !local strain variable for elastic strain energy
    REAL        :: Strain_ini(1:nIntGP,6)                                     !local initial strain variable for elastic strain energy
    REAL        :: devStress(1:nIntGP,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nIntGP)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor
    REAL        :: mu, lambda                                                 ! Lame parameters
    REAL        :: tau(1:nIntGP)
    REAL        :: taulim(1:nIntGP)                                           !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv(1:nIntGP)                                           !second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv                                           !plastic parameters
    REAL        :: DOFStress(1:nDegFr,1:6)                                    !initial loading
    REAL        :: dgvar(1:nAlignedDegFr,1:9)                                 !dof's
    REAL        :: dgvar_new(1:nAlignedDegFr,1:6)                             !
    REAL        :: dudt_plastic(1:nDegFr,1:6)                                 !stress change due to plasticity
    REAL        :: dudt_pstrain(1:nDegFr,1:6)                                 !change of plastic strain, dof-wise
    REAL        :: dudt_pstrainGP(1:nIntGP,1:6)                               !change of plastic strain, GP-wise
    REAL        :: pstrain(1:nDegFr,1:7)                                      !plastic strain
    REAL        :: pstrain_gp(1:nIntGP,1:6)                                   !plastic strain at GP
    REAL        :: estrain(1:nDegFr,1:6), estrain_ini(1:nDegFr,1:6)           !total elastic strain
    REAL        :: phi                                                       ! Value of the base function at GP      !
    LOGICAL     :: check
    REAL        :: update(1:nIntGP,6)
    REAL        :: newstateGP(1:6)
    REAL        :: PlasticEnergy_tmp, PValue, EstrainEnergy_tmp   !temp. energies
    REAL        :: Energy(1:2)                                               !plastic and elastic strain energy
    REAL        :: parameters(1:3)                                           !1=volume of the triangle, 2=plastcohesion, 3=density rho
    REAL        :: I1,I1_0,I2,I2_0                                           !first and second invariants of strains
    REAL, POINTER :: IntGaussP(:,:)     =>NULL()
    REAL, POINTER :: IntGaussW(:)       =>NULL()
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, nDegFr, BulkFriction, Tv, dt, mu, lambda, parameters, &
                     intGaussP, intGaussW, nVar, nIntGP

                     !IntGPBaseFunc, MassMatrix
    INTENT(INOUT) :: dgvar, pstrain, pstrain_gp, Energy
    !-------------------------------------------------------------------------!
    dudt_plastic = 0.0D0
    dudt_pstrain = 0.0D0
    dudt_pstrain_gp = 0.0D0
    Stress_total = 0.0D0
    Energy(1:2)  = 0.0D0

    angfric = ATAN(BulkFriction) !angle of friction

    IF (Tv .GT. 0) THEN
       relaxtime = 1.0D0 - EXP(-dt/(Tv)) !Tv: direct input via parameter file; Tv smaller-> stronger plasticity
    ELSE
       relaxtime = 1.0
    ENDIF

    iPoly = DISC%Galerkin%nPoly
    ! Basis func values
    IntGPBaseFunc => DISC%Galerkin%IntGPBaseFunc_Tet(1:nDegFr,1:nIntGP,iPoly)
    ! Mass matrix
    MassMatrix    => DISC%Galerkin%MassMatrix_Tet(1:nDegFr,1:nDegFr,iPoly)


    ! ---[ Calculate trial stress tensor ]---
    ! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz
    ! -> need to be specified throughout the whole medium for every element
    ! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    ! first approach: don't use the mapping of DOFStress and just add it to the get state vector; GP-wise

    DO iIntGP = 1, nIntGP
        !1. get the state at every GP
        stateGP(:) = 0.
        DO iDegFr = 1, nDegFr
           phi = IntGPBaseFunc(iDegFr,iIntGP)
           stateGP(1:nVar) = stateGP(1:nVar) + phi*dgvar(iDegFr,1:nVar)
        ENDDO

        !2. add up all initial loading (constant in an element)
        pstateGP(:) = stateGP(:)
        Stress_total(iIntGP,1:6) = pstateGP(1:6) + DOFStress(1,1:6) !dofstress are in this case just the elementwise initial stresses

        !Difference between high and low order: Stresses are calculated at GP and corresponding strains are also the real values at these points
        !Calculate the total strain from the elastic stress-strain relation
        Strain_total(iIntGP, 1:6) = MATMUL(DISC%Galerkin%Strain_matrix, Stress_total(iIntGP,1:6))
        !Calculate initial strain loading from initial stress loading (elementwise) -> move that outside the routine and calculate beforhand
        Strain_ini(iIntGP,1:6) = MATMUL(DISC%Galerkin%Strain_matrix, DOFStress(1,1:6))
    ENDDO


    !Calculate initial strain loading from initial stress loading (elementwise) -> move that outside the routine and calculate before as it is constant over time
    !Strain_ini(1:6) = MATMUL(DISC%Galerkin%Strain_matrix,DOFStress(1,1:6))

    ! Mean stress, GP-wise
    meanStress(1:nIntGP) = (Stress_total(:,1) + Stress_total(:,2)+ Stress_total(:,3) )/3

    ! Deviatoric stress, GP-wise
    devStress(1:nIntGP,1) = Stress_total(:,1) - meanStress(:)
    devStress(1:nIntGP,2) = Stress_total(:,2) - meanStress(:)
    devStress(1:nIntGP,3) = Stress_total(:,3) - meanStress(:)

    devStress(1:nIntGP,4:6) = Stress_total(:,4:6)


    ! Second invariant of stress deviator
    secInv(1:nIntGP) = 0.5*(devStress(:,1)**2+devStress(:,2)**2+devStress(:,3)**2)+devStress(:,4)**2+devStress(:,5)**2+devStress(:,6)**2

    ! Scalar measure of shear stress
    tau(1:nIntGP)= SQRT(secInv(1:nIntGP))

    ! Yield stress
    taulim(1:nIntGP) = parameters(2)*COS(angfric) - meanStress(1:nIntGP)*SIN(angfric)! minus before sinus is for compressional stress=negative.
    taulim(1:nIntGP) = MAX(0.0, taulim(1:nIntGP))

    check = .false.
    ! Stress deviators are adjusted
    DO iIntGP = 1, nIntGP
       yldfac = 1

       IF (tau(iIntGP) .GT. taulim(iIntGP)) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
           if (check .EQ. .false.) then
           check = .true.
           endif
           yldfac = 1.0D0- (1.0D0 - taulim(iIntGP)/tau(iIntGP))*(relaxtime) !factor by Duan/Day 2008
           ! adjustment of stresses, GP-wise for every variable 1-6

           Stress_total(iIntGP,1) = devStress(iIntGP,1)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,2) = devStress(iIntGP,2)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,3) = devStress(iIntGP,3)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,4:6) = devStress(iIntGP,4:6)*yldfac
       ENDIF

       !update(iIntGP,1:6) = Stress_total(iIntGP, 1:6)-EQN%IniStress(1:6, iElem) !subtract the inital loading; new stress state at every GP
       update(iIntGP,1:6) = Stress_total(iIntGP, 1:6)-DOFStress(1,1:6) !subtract the inital loading; new stress state at every GP

       !update dudt_pstrain GP-wise
       dudt_pstrain_gp(iIntGP,1:6) = ((1-yldfac)/(2.0*mu))*devStress(iIntGP,1:6)

    ENDDO !adjustment over all GP-points


   IF (check .EQ. .false.) THEN
      dudt_plastic(1:nDegFr,1:6) = 0.0
      !
   ELSE !back projection is necessary because at least one GP needs to be adjusted
    !muss vorher einmal allokiert werden in dg_Setup
   dgvar_new = 0.0

    ! back projection to the DOFs
    DO iIntGP = 1, nIntGP

         DO iDegFr = 1, nDegFr
           phi = IntGPBaseFunc(iDegFr,iIntGP) !basis function number idegfr at point iIntGp
           newstateGP = update(iIntGP,1:6) !new value at that GP (=old value if no adjustment)
           !dgvar kann hier nicht benutzt werden, da ungleich Null, geht nur am Anfang, wenn alles Null ist
          dgvar_new(iDegFr,1:6) =  dgvar_new(iDegFr,1:6) + IntGaussW(iIntGP)*newstateGP(1:6)*phi
        ENDDO
     ENDDO !nIntGP
     !divide by diagonal mass matrix entries
     DO iDegFr = 1, nDegFr
            dgvar_new(iDegFr,:) = dgvar_new(iDegFr,:) / MassMatrix(iDegFr,iDegFr)
     ENDDO

     !Update dof-wise
     dudt_plastic(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6)- dgvar_new(1:nDegFr,1:6)
     dudt_pstrain(1:nDegFr,1:6) = (1/(2.0*mu))*dudt_plastic(1:nDegFr,1:6) !for the plastic strain just take the first dof


   ENDIF !check = .true.

    dgvar(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6) - dudt_plastic(1:nDegFr,1:6)

!============ LOW-ORDER
    !update plastic strain
    !pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6) !plastic strain tensor
    !accumulated plastic strain
    !pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                 !+ dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

    !calculate energies
    !PlasticEnergy_tmp = Stress_total(1,1)*dudt_pstrain(1) + Stress_total(1,2)*dudt_pstrain(2) + Stress_total(1,3)*dudt_pstrain(3) + 2.0*Stress_total(1,4)*dudt_pstrain(4) &
                      !+ 2.0*Stress_total(1,5)*dudt_pstrain(5) + 2.0*Stress_total(1,6)*dudt_pstrain(6)

    !estrain(1:6) = Strain_total(1,1:6) - dudt_pstrain(1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    !estrain_ini(1:6) = Strain_ini(1:6)

    !first and second invariants of the total elastic strain
    !I1 = estrain(1) + estrain(2) + estrain(3)
    !I2 = estrain(1)**2 + estrain(2)**2 + estrain(3)**2 + 2.0*estrain(4)**2 + 2.0*estrain(5)**2 + 2.0*estrain(6)**2
    !first and second invariants of the initial strain loading
    !I1_0 = estrain_ini(1) + estrain_ini(2) + estrain_ini(3)
    !I2_0 = estrain_ini(1)**2 + estrain_ini(2)**2 + estrain_ini(3)**2 + 2.0*estrain_ini(4)**2 + 2.0*estrain_ini(5)**2 + 2.0*estrain_ini(6)**2
    !Elastic strain energy
    !subtracted the initial elastic strain
    !EstrainEnergy_tmp = 0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0)

    !Energy(1) = PlasticEnergy_tmp*parameters(1) !multiplied by volume to get integral over element
    !Energy(2) = EstrainEnergy_tmp*parameters(1)
!============== HIGH-ORDER

    !update plastic strain, degFr-wise but only first degfr is outputted
    pstrain(1:nDegFr,1:6) = pstrain(1:nDegFr,1:6) + dudt_pstrain(1:nDegFr,1:6) !plastic strain tensor
    !accumulated plastic strain
    pstrain(1:nDegFr,7) = pstrain(1:nDegFr,7)+ dt*sqrt(0.5*(dudt_pstrain(1:nDegFr,1)**2 + dudt_pstrain(1:nDegFr,2)**2 &
                                                   + dudt_pstrain(1:nDegFr,3)**2)+ dudt_pstrain(1:nDegFr,4)**2 + dudt_pstrain(1:nDegFr,5)**2 + dudt_pstrain(1:nDegFr,6)**2)

    !update plastic strain but GP-wise
    pstrain_gp(1:nIntGP,1:6) = pstrain_gp(1:nIntGP,1:6) + dudt_pstrain_gp(1:nIntGP,1:6)

    !calculate energies here the strains are GP wise
    !estrain(1:nDegFr,1:6) = Strain_total(1:nDegFr,1:6) - dudt_pstrain(1:nDegFr,1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    estrain(1:nIntGP,1:6) = Strain_total(1:nIntGP,1:6) -pstrain_gp(1:nIntGP,1:6) ! !total elastic strain, if no plastic yielding -> elastic strain = total strain

    estrain_ini(1:nIntGP,1:6) = Strain_ini(1:nIntGP,1:6)

    EstrainEnergy_tmp = 0.0
    PlasticEnergy_tmp = 0.0

    !Caluclated GP-wise and integrated with gaussian quadrature
    DO iIntGP=1,nIntGP
        I1 = 0.0
        I2 = 0.0
        I1_0 = 0.0
        I2_0 = 0.0
        PValue = 0.0

        !first and second invariants of the total elastic strain
        I1 = estrain(iIntGP,1) + estrain(iIntGP,2) + estrain(iIntGP,3)
        I2 = estrain(iIntGP,1)**2 + estrain(iIntGP,2)**2 + estrain(iIntGP,3)**2 + 2.0D0*estrain(iIntGP,4)**2 &
             + 2.0D0*estrain(iIntGP,5)**2 + 2.0D0*estrain(iIntGP,6)**2
        !first and second invariants of the initial strain loading
        I1_0 = estrain_ini(iIntGP,1) + estrain_ini(iIntGP,2) + estrain_ini(iIntGP,3)
        I2_0 = estrain_ini(iIntGP,1)**2 + estrain_ini(iIntGP,2)**2 + estrain_ini(iIntGP,3)**2 + 2.0D0*estrain_ini(iIntGP,4)**2 &
               + 2.0D0*estrain_ini(iIntGP,5)**2 + 2.0D0*estrain_ini(iIntGP,6)**2

        PValue = Stress_total(iIntGP,1)*dudt_pstrain_gp(iIntGP,1) + Stress_total(iIntGP,2)*dudt_pstrain_gp(iIntGP,2) &
                 + Stress_total(iIntGP,3)*dudt_pstrain_gp(iIntGP,3) + 2.0D0*Stress_total(iIntGP,4)*dudt_pstrain_gp(iIntGP,4) &
                 + 2.0D0*Stress_total(iIntGP,5)*dudt_pstrain_gp(iIntGP,5) + 2.0D0*Stress_total(iIntGP,6)*dudt_pstrain_gp(iIntGP,6)

        !Elastic strain energy
        !subtracted the initial elastic strain
        EstrainEnergy_tmp = EstrainEnergy_tmp + IntGaussW(iIntGP)*(0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0))
        !Plastic energy
        PlasticEnergy_tmp = PlasticEnergy_tmp + IntGaussW(iIntGP)*(PValue)


    ENDDO

    Energy(1) = PlasticEnergy_tmp*6.0D0*parameters(1) !multiplied by |J| for the reference element
    Energy(2) = EstrainEnergy_tmp*6.0D0*parameters(1)



 END SUBROUTINE Plasticity_3D_high


!yldfac is only caluclated from the first DOF, and all DOF's are adjusted by the same coefficient

  SUBROUTINE Plasticity_3D_DOF(DISC,dgvar, DOFStress, nDegFr, nAlignedDegFr, BulkFriction, Tv, dt, mu,lambda, parameters , Energy, pstrain)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    !TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    !TYPE(tUnstructMesh)      :: MESH
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: nDegFr
    integer     :: nAlignedDegFr

    REAL        :: Stress(1:nDegFr,6)                                         !local stress variable for the yield criterion
    REAL        :: Strain_ini(1:nDegFr,6)                                            !local initial strain variable for elastic strain energy
    REAL        :: Strain_total(1:nDegFr,6)                                          !local strain variable for elastic strain energy
    REAL        :: devStress(1:nDegFr,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nDegFr)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor  
    REAL        :: mu, lambda                                                 !Lame parameters mu and lambda
    REAL        :: tau,taulim                                                 !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv                                                    !secInv=second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv
    REAL        :: DOFStress(1:nDegFr,1:6)
    REAL        :: dgvar(1:nAlignedDegFr,1:6)
    REAL        :: dudt_pstrain(1:nDegFr,1:6)
    REAL        :: pstrain(1:nDegFr,1:7)
    REAL        :: estrain(1:nDegFr,1:6), estrain_ini(1:nDegFr,1:6)                             !total elastic strain
    REAL        :: PlasticEnergy_tmp, EstrainEnergy_tmp
    REAL        :: Energy(1:2), massMatrix
    REAL        :: parameters(1:3)
    REAL        :: I1,I1_0,I2,I2_0
    INTEGER     :: LocnPoly
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, DOFStress, nDegFr, BulkFriction, Tv, dt, mu, lambda, parameters
    INTENT(INOUT) :: dgvar, pstrain, Energy
    !-------------------------------------------------------------------------!

    dudt_pstrain = 0.0
    LocnPoly     = DISC%Galerkin%nPoly
    Energy(1:2)  = 0.0


    angfric = ATAN(BulkFriction) !angle of friction

    IF (Tv .GT. 0) THEN
       relaxtime = 1.0D0 - EXP(-dt/(Tv)) !Tv: direct input via parameter file; Tv smaller-> stronger plasticity
    ELSE
       relaxtime = 1.0
    ENDIF


! ---[ Calculate trial stress tensor ]---
! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz -> need to be specified throughout the whole medium for every element 
! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    Stress(:,1:6)= dgvar(1:nDegFr,1:6)  + DOFStress(:,1:6)   !act.Stress + initial stress_xx

    !Calculate the total strain from the elastic stress-strain relation
    DO iDegFr = 1,nDegFR
       Strain_total(iDegFr,1:6) = MATMUL(DISC%Galerkin%Strain_matrix, Stress(iDegFr,1:6))

       !Calculate initial strain loading from initial stress loading (elementwise) -> move that outside the routine and calculate beforhand
        Strain_ini(iDegFr,1:6) = MATMUL(DISC%Galerkin%Strain_matrix, DOFStress(iDegFr,1:6))
    ENDDO

! ---[ Calculate trial yield stress ]---

    ! Mean stress
    meanStress(1:nDegFr) = (Stress(1:nDegFr,1) + Stress(1:nDegFr,2)+ Stress(1:nDegFr,3) )*(1.0D0/3.0D0)

    ! Deviatoric stress
    devStress(1:nDegFr,1) = Stress(1:nDegFr,1) - meanStress(1:nDegFr)
    devStress(1:nDegFr,2) = Stress(1:nDegFr,2) - meanStress(1:nDegFr)
    devStress(1:nDegFr,3) = Stress(1:nDegFr,3) - meanStress(1:nDegFr)

    devStress(1:nDegFr,4:6) = Stress(1:nDegFr,4:6)

    ! Second invariant of stress deviator
    secInv = 0.5*(devStress(1,1)**2+devStress(1,2)**2+devStress(1,3)**2)+devStress(1,4)**2+devStress(1,5)**2+devStress(1,6)**2 

    ! Scalar measure of shear stress
    tau= SQRT(secInv)

    ! Yield stress   
    taulim = parameters(2)*COS(angfric) - meanStress(1)*SIN(angfric)! minus before sinus is for compressional stress=negative.
    taulim = MAX(0.0, taulim)


    ! Stress deviators are adjusted

    IF (tau .GT. taulim) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
       yldfac = 1.0D0- (1.0D0 - taulim/tau)*(relaxtime) !factor by Duan/Day


       Stress(1:nDegFr,1) = devStress(1:nDegFr,1)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,2) = devStress(1:nDegFr,2)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,3) = devStress(1:nDegFr,3)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,4:6) = devStress(1:nDegFr,4:6)*yldfac


       !----Change the dofs-----
       dgvar(1:nDegFr,1:6) = Stress(1:nDegFr,1:6) - DOFStress(1:nDegFr,1:6)
       !dividing by 2*mu due to tensor convention
       dudt_pstrain(1:nDegFr,1:6) = ((1-yldfac)/(2*mu))*devStress(1:nDegFr, 1:6) !only the first dof is considered for plastic strain tensor


        
    ENDIF !yield criterion check

!================== LOW-ORDER
    !update plastic strain
    !pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6) !plastic strain tensor

    !accumulated plastic strain
    !pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                                                   !+ dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

    !calculate energies
    !PlasticEnergy_tmp = Stress(1,1)*dudt_pstrain(1) + Stress(1,2)*dudt_pstrain(2) + Stress(1,3)*dudt_pstrain(3) + 2.0*Stress(1,4)*dudt_pstrain(4) &
                      !+ 2.0*Stress(1,5)*dudt_pstrain(5) + 2.0*Stress(1,6)*dudt_pstrain(6)

    !estrain(1:6) = Strain_total(1:6) - dudt_pstrain(1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    !estrain_ini(1:6) = Strain_ini(1:6)

    !first and second invariants of the total elastic strain
    !I1 = estrain(1) + estrain(2) + estrain(3)
    !I2 = estrain(1)**2 + estrain(2)**2 + estrain(3)**2 + 2.0*estrain(4)**2 + 2.0*estrain(5)**2 + 2.0*estrain(6)**2
    !first and second invariants of the initial strain loading
    !I1_0 = estrain_ini(1) + estrain_ini(2) + estrain_ini(3)
    !I2_0 = estrain_ini(1)**2 + estrain_ini(2)**2 + estrain_ini(3)**2 + 2.0*estrain_ini(4)**2 + 2.0*estrain_ini(5)**2 + 2.0*estrain_ini(6)**2
    !Elastic strain energy
    !subtracted the initial elastic strain
    !EstrainEnergy_tmp = 0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0)

    !Energy(1) = PlasticEnergy_tmp*parameters(1) !multiplied by volume to get integral over element
    !Energy(2) = EstrainEnergy_tmp*parameters(1)

!================= HIGH_ORDER
    !update plastic strain, only first degFr is outputted
    pstrain(1:nDegFr,1:6) = pstrain(1:nDegFr,1:6) + dudt_pstrain(1:nDegFr,1:6) !plastic strain tensor
    !accumulated plastic strain
    pstrain(1:nDegFr,7) = pstrain(1:nDegFr,7)+ dt*sqrt(0.5*(dudt_pstrain(1:nDegFr,1)**2 + dudt_pstrain(1:nDegFr,2)**2 &
                                                   + dudt_pstrain(1:nDegFr,3)**2)+ dudt_pstrain(1:nDegFr,4)**2 + dudt_pstrain(1:nDegFr,5)**2 + dudt_pstrain(1:nDegFr,6)**2)

    !calculate energies
    !estrain(1:nDegFr,1:6) = Strain_total(1:nDegFr,1:6) - dudt_pstrain(1:nDegFr,1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    estrain(1:nDegFr,1:6) = Strain_total(1:nDegFr,1:6) -pstrain(1:nDegFr,1:6) ! !total elastic strain, if no plastic yielding -> elastic strain = total strain

    estrain_ini(1:nDegFr,1:6) = Strain_ini(1:nDegFr,1:6)


    EstrainEnergy_tmp = 0.0
    PlasticEnergy_tmp = 0.0

    DO iDegFr=1,nDegFr
        I1=0.0
        I2=0.0
        I1_0=0.0
        I2_0=0.0
        massMatrix = DISC%Galerkin%MassMatrix_Tet(iDegFr,iDegFr,LocnPoly)
        !first and second invariants of the total elastic strain
        I1 = estrain(iDegFr,1) + estrain(iDegFr,2) + estrain(iDegFr,3)
        I2 = estrain(iDegFr,1)**2 + estrain(iDegFr,2)**2 + estrain(iDegFr,3)**2 + 2.0D0*estrain(iDegFr,4)**2 &
             + 2.0D0*estrain(iDegFr,5)**2 + 2.0D0*estrain(iDegFr,6)**2
        !first and second invariants of the initial strain loading
        I1_0 = estrain_ini(iDegFr,1) + estrain_ini(iDegFr,2) + estrain_ini(iDegFr,3)
        I2_0 = estrain_ini(iDegFr,1)**2 + estrain_ini(iDegFr,2)**2 + estrain_ini(iDegFr,3)**2 + 2.0D0*estrain_ini(iDegFr,4)**2 &
               + 2.0D0*estrain_ini(iDegFr,5)**2 + 2.0D0*estrain_ini(iDegFr,6)**2
        !Elastic strain energy
        !subtracted the initial elastic strain
        EstrainEnergy_tmp = EstrainEnergy_tmp + massMatrix*(0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0))

        PlasticEnergy_tmp = PlasticEnergy_tmp + massMatrix*(Stress(iDegFr,1)*dudt_pstrain(iDegFr,1) + Stress(iDegFr,2)*dudt_pstrain(iDegFr,2) + Stress(iDegFr,3)*dudt_pstrain(iDegFr,3) &
                            + 2.0D0*Stress(iDegFr,4)*dudt_pstrain(iDegFr,4) + 2.0D0*Stress(iDegFr,5)*dudt_pstrain(iDegFr,5) + 2.0D0*Stress(iDegFr,6)*dudt_pstrain(iDegFr,6))
    ENDDO

    Energy(1) = PlasticEnergy_tmp*6.0D0*parameters(1) !multiplied by volume to get integral over element
    Energy(2) = EstrainEnergy_tmp*6.0D0*parameters(1)


 END SUBROUTINE Plasticity_3D_DOF


END MODULE Plasticity_mod
