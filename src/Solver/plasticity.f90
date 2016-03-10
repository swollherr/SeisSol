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


   SUBROUTINE Plasticity_3D_high(dgvar, DOFStress, nDegFr, nAlignedDegFr, BulkFriction, Tv, PlastCo, dt, mu, pstrain, intGaussP, intGaussW, &
                                    !IntGPBaseFunc, MassMatrix, &
                                    DISC, nVar, nIntGP)


    !-------------------------------------------------------------------------!
  USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization)    :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: nDegFr
    INTEGER     :: iIntGP, nIntGP
    INTEGER     :: nVar, iPoly
    integer     :: nAlignedDegFr


    REAL        :: stateGP(nVar)     ! State in GP            !
    REAL        :: pstateGP(nVar)    ! Primitive state in GP  !
    REAL        :: Stress_total(nIntGP,6)                                         !local stress variable for the yield criterion
    REAL        :: devStress(1:nIntGP,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nIntGP)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor
    REAL        :: mu                                                         ! Lame parameter mu
    REAL        :: LocnVar                                          !
    REAL        :: tau(1:nIntGP)
    REAL        :: taulim(1:nIntGP)                                                !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv(1:nIntGP)                                               !secInv=second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv, PlastCo
    REAL        :: DOFStress(1:nDegFr,1:6)
    REAL        :: dgvar(1:nAlignedDegFr,1:6)
    REAL        :: dgvar_new(1:nAlignedDegFr,1:6)
    REAL        :: dudt_plastic(1:nDegFr,1:6)
    REAL        :: dudt_pstrain(1:6)
    REAL        :: pstrain(1:7)
    REAL, POINTER :: IntGaussP(:,:)     =>NULL()
    REAL, POINTER :: IntGaussW(:)       =>NULL()
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()
    REAL         :: phi                                                            ! Value of the base function at GP      !
    LOGICAL      :: check
    REAL         :: update(1:nIntGP,6)
    REAL         :: newstateGP(1:6)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: nDegFr, BulkFriction, Tv, PlastCo, dt, mu, intGaussP, intGaussW, DISC, nVar, nIntGP
                     !IntGPBaseFunc, MassMatrix
    INTENT(INOUT) :: dgvar, pstrain
    !-------------------------------------------------------------------------!
    dudt_plastic = 0.0
    dudt_pstrain = 0.0
    Stress_total = 0.0

    angfric = ATAN(BulkFriction) !angle of friction
    relaxtime = Tv !direct input via parameter file
    !dt/(Tv) !Tv=dx/V_s with dx=min(dx);  Tv smaller-> stronger plasticity
    iPoly = DISC%Galerkin%nPoly
    ! Basis func values
    IntGPBaseFunc => DISC%Galerkin%IntGPBaseFunc_Tet(1:nDegFr,1:nIntGP,iPoly)
    ! Mass matrix
    MassMatrix    => DISC%Galerkin%MassMatrix_Tet(1:nDegFr,1:nDegFr,iPoly)


! ---[ Calculate trial stress tensor ]---
! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz -> need to be specified throughout the whole medium for every element
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
    ENDDO

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
    taulim(1:nIntGP) = PlastCo*COS(angfric) - meanStress(1:nIntGP)*SIN(angfric)! minus before sinus is for compressional stress=negative.
    taulim(1:nIntGP) = MAX(0.0, taulim(1:nIntGP))

    check = .false.
    ! Stress deviators are adjusted
    DO iIntGP = 1, nIntGP

       IF (tau(iIntGP) .GT. taulim(iIntGP)) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
           if (check .EQ. .false.) then
           check = .true.
           endif
           yldfac = 1.0D0- (1.0D0 - taulim(iIntGP)/tau(iIntGP))*(1.0D0 - EXP(-relaxtime)) !factor by Duan/Day 2008
           ! adjustment of stresses, GP-wise for every variable 1-6

           Stress_total(iIntGP,1) = devStress(iIntGP,1)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,2) = devStress(iIntGP,2)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,3) = devStress(iIntGP,3)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,4:6) = devStress(iIntGP,4:6)*yldfac
       ENDIF

       !update(iIntGP,1:6) = Stress_total(iIntGP, 1:6)-EQN%IniStress(1:6, iElem) !subtract the inital loading; new stress state at every GP
       update(iIntGP,1:6) = Stress_total(iIntGP, 1:6)-DOFStress(1,1:6) !subtract the inital loading; new stress state at every GP
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

     !Update
      dudt_plastic(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6)- dgvar_new(1:nDegFr,1:6)
      dudt_pstrain(1:6) = (1/mu)*dudt_plastic(1,1:6) !for the plastic strain just take the first dof
     !mu or 2*mu?
   ENDIF !check = .true.

   dgvar(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6) - dudt_plastic(1:nDegFr,1:6)

    !update plastic strain
    pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6) !plastic strain tensor
    !accumulated plastic strain
    pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                                                   + dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)


 END SUBROUTINE Plasticity_3D_high


!yldfac is only caluclated from the first DOF, and all DOF's are adjusted by the same coefficient
  SUBROUTINE Plasticity_3D_DOF(dgvar, DOFStress, nDegFr, nAlignedDegFr, BulkFriction, Tv, dt, mu, parameters , Energy, pstrain)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    !TYPE(tEquations)         :: EQN
    !TYPE(tDiscretization)    :: DISC
    !TYPE(tUnstructMesh)      :: MESH
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: nDegFr
    integer     :: nAlignedDegFr

    REAL        :: Stress(1:nDegFr,6)                                         !local stress variable for the yield criterion
    REAL        :: devStress(1:nDegFr,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nDegFr)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor  
    REAL        :: mu                                                         ! Lame parameter mu
    REAL        :: LocnVar                                          !
    REAL        :: tau,taulim                                                 !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv                                                     !secInv=second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv
    REAL        :: DOFStress(1:nDegFr,1:6)
    REAL        :: dgvar(1:nAlignedDegFr,1:6)
    REAL        :: dudt_pstrain(1:6)
    REAL        :: pstrain(1:7)
    REAL        :: PlasticEnergy_tmp, Energy(1:2)
    REAL        :: parameters(1:2)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DOFStress, nDegFr, BulkFriction, Tv, dt, mu, parameters
    INTENT(INOUT) :: dgvar, pstrain, Energy
    !-------------------------------------------------------------------------!

    dudt_pstrain = 0.0



    angfric = ATAN(BulkFriction) !angle of friction
    relaxtime = dt/(Tv) !Tv=dx/V_s with dx=min(dx);  Tv smaller-> stronger plasticity


! ---[ Calculate trial stress tensor ]---
! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz -> need to be specified throughout the whole medium for every element 
! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    Stress(:,1:6)= dgvar(1:nDegFr,1:6)  + DOFStress(:,1:6)   !act.Stress + initial stress_xx

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
       yldfac = 1.0D0- (1.0D0 - taulim/tau)*(1.0D0 - EXP(-relaxtime)) !factor by Duan/Day


       Stress(1:nDegFr,1) = devStress(1:nDegFr,1)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,2) = devStress(1:nDegFr,2)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,3) = devStress(1:nDegFr,3)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,4:6) = devStress(1:nDegFr,4:6)*yldfac


       !----Change the dofs-----
       dgvar(1:nDegFr,1:6) = Stress(1:nDegFr,1:6) - DOFStress(1:nDegFr,1:6)


       dudt_pstrain(1:6) = ((1-yldfac)/mu)*devStress(1, 1:6) !only the first dof is considered for plastic strain tensor

        
    ENDIF !yield criterion check

    !update plastic strain
    pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6) !plastic strain tensor
    !calculate increment of dissipated plastic energy for this element
    !take stress or dgvar?
    PlasticEnergy_tmp = dgvar(1,1)*dudt_pstrain(1) + dgvar(1,2)*dudt_pstrain(2) + dgvar(1,3)*dudt_pstrain(3) + 2*dgvar(1,4)*dudt_pstrain(4) &
                      + 2*dgvar(1,5)*dudt_pstrain(5) + 2*dgvar(1,6)*dudt_pstrain(6)
    Energy(1) = PlasticEnergy_tmp*parameters(1) !volume

    !accumulated plastic strain
    pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                                                   + dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

 END SUBROUTINE Plasticity_3D_DOF

END MODULE Plasticity_mod
