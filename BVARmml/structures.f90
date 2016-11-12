module structures

!type Matrix
!    double precision, allocatable :: M(:,:)
!end type

integer, parameter :: varnamelength = 8, tokenlength = 100

type BreakInfo
    integer nr, rsize   ! number of regimes, min regime size
    integer, allocatable :: bdates(:)
    character(len=3), allocatable :: btypes(:)
    logical fixed
end type

type ModelInfo
    character*4 modtype
    character*varnamelength, allocatable :: dependentvar(:), deterministicvar(:), exogvar(:)
    integer lags
    type(BreakInfo) breaks
end type

type DataArray
    double precision, allocatable :: datamat(:,:)
    character*varnamelength, allocatable :: varnames(:)
end type

type MCMCInfo
    integer reps, burn, thin
    character(len=tokenlength) dumpfile
end type

type LittermanSpec
    double precision, allocatable :: tightness(:), firstlagmean(:), &
                                     deterministicmean(:,:), deterministicvar(:,:)
end type

type PriorInfo
    character*4 priortype
    type(LittermanSpec) Litterman
    integer priordf
end type

type VARData
	! Data matrices and OLS estimates
	! T = number of obs
	! m = number of variables dimension of model
	! p = lag length
	! q = number of deterministic vars
	! k = number parameters in one equation
	double precision, allocatable :: Y(:,:), Z(:,:), D(:,:), W(:,:), B(:,:), &
	                                 ZPZ(:,:), ZPY(:,:), &
									 Gamma_OLS(:,:), vecGamma_OLS(:), ESS(:,:)
	character(len=varnamelength), allocatable :: Ylabels(:)
	character(len=varnamelength+3), allocatable :: Zlabels(:)
	integer Nobs, firstobs, lastobs, m, p, q, k
end type

type SampleInfo
    integer first, last
end type

type priorND
	! Normal-Diffuse prior
	! prior mean and variance of regression parameters
	! prior variance as full matrix or vector if diagonal
	double precision, allocatable :: Gamma(:,:), & ! prior mean
	                                 varGamma(:,:), & ! prior variance
	                                 invvarGamma(:,:), & ! inverse of prior variance
	                                 diagvarGamma(:), & ! prior var if diagonal
	                                 invdiagvarGamma(:), &  ! invers if variance is diagonal
									 invvarGammaGamma(:), &  ! Sigma^-1 * vec(Gamma)
									 priorPsi(:,:), priorVecPsi(:)
    integer priordf	
    logical diagPsi								 
end type

type priorNW
	! Normal-Wishart prior
	! Prior mean and variance matrices, for conditional normal and inverse Wishart
	! variance matrices may be full or diagonal
	double precision, allocatable :: priorGamma(:,:), priorOmega(:,:), priorVecOmega(:), &
									 priorPsi(:,:), priorVecPsi(:)
	integer priordf
	logical diagOmega, diagPsi
end type



type priorSS
    ! prior for Steady State model
    ! prior mean and variance of parameters for steady state (Lambda)
    ! and dynamics (Gamma), Normal
    ! error variance (Psi) is inverse Wishart
    double precision, allocatable :: Gamma(:,:), diagvarGamma(:), invdiagvarGamma(:), &
                                     invvarGammaGamma(:), &
                                     Lambda(:,:), diagvarLambda(:), &
                                     invdiagvarLambda(:), invvarLambdaLambda(:), Psi(:,:)
     integer priordf
end type

type Breakprior
    integer minregime
    double precision breakprob
end type

type VARprior
	character*4 priortype
	type(priorND) :: ND
	type(priorNW) :: NW
	type(priorSS) :: SS
	type(Breakprior) :: breaks
end type

type NDsampler
    ! *_0 is prior quantities
    ! *_1 is posterior quantities
    
	! state of sampler
	double precision, allocatable :: invPsi_Sampled(:,:,:) ! inverse of current Psi
	! temporary arrays
	double precision, allocatable :: varGamma_1(:,:), SGSG(:), Gamma_1(:), &
									 Gamma_gen(:), err(:,:), ESS(:,:), &
									 Psi_root(:,:), S_chol(:,:)
end type

type NIWposterior

     double precision, allocatable :: varGamma(:,:), vecGamma(:), &
									 S_chol(:,:)
								
END TYPE


type NWposterior

     double precision, allocatable :: Omega_invchol(:,:),Omega(:,:), vecGamma(:), &
									 S_chol(:,:)
	integer postdf							
END TYPE


type NWsampler
	double precision, allocatable :: Gamma(:,:,:),Psi(:,:,:),Gammapost(:,:),Omegapost(:,:),Psipost(:,:),psipostchol(:),omegapostchol(:)
end type

type SSsampler
    type(BreakInfo) :: breaks
    double precision, allocatable :: Lambda_Sampled(:,:,:), invPsi_Sampled(:,:,:), &
                                     YGamma(:,:), ZGamma(:,:), YLambda(:,:), &
                                     XTX(:,:), SGSG(:), BB(:,:), BBPsi(:,:), &
                                     F(:,:), SLSL(:), &
                                     Gamma_1(:), varGamma_1(:,:), &
                                     Gamma_gen(:), &
                                     Lambda_1(:), varLambda_1(:,:), Lambda_gen(:), &
                                     err(:,:), ESS(:,:), &
									 Psi_root(:,:), S_chol(:,:) 
end type

type SamplerState
	! keep state and precalculated quantities for samplers
	type(NDsampler) :: ND
	type(NWsampler) :: NW
	type(SSsampler) :: SS
	integer firstobs, lastobs
end type

! Start values for samples

!type SS_start
!    double precision, allocatable :: Psi(:,:), Lambda(:,:)
!end type

type StartVal
!    type(SS_start) :: SS
    double precision, allocatable :: Psi(:,:), Lambda(:,:)
    type(BreakInfo) :: breaks
end type

type RMSEs
    real*8, allocatable :: bma(:,:,:,:), top(:,:,:,:), med(:,:,:,:), ar(:,:,:,:), rw(:,:,:,:), rmean(:,:,:,:)
end type

type RMSEss
    real*8, allocatable :: bma(:,:), top(:,:), med(:,:), ar(:,:), rw(:,:), rmean(:,:)
end type


type mats
    real*8, allocatable :: bma(:,:,:,:), top(:,:,:,:), med(:,:,:,:), ar(:,:,:,:), rw(:,:,:,:), rmean(:,:,:,:),main(:,:,:,:)
end type

type weights
    real*8, allocatable  :: priorval(:,:),dims(:,:),var(:,:),sumposts(:,:),marg(:,:,:),pred(:,:,:)
    integer, allocatable :: topmod(:,:,:) ,ranksum(:,:)
end type 

type mods
    integer n,p,q,num,no_mods_del,true
    character*12 :: posteval
    real*8, allocatable  :: lnprior(:,:),priorsize(:)
    integer, allocatable :: models(:,:),evalvec(:),bma_comb(:),rmse_obs(:)
end type 

type dgp
    integer N,H,p,m_syst,m_spur,m,q,mcburn,mcreps,simreps,present_top,dims(2)
    real*8, allocatable :: param_syst(:,:,:), param_spur(:,:,:),var_err(:)
    character*4  :: dgps 
    CHARACTER(len=12),allocatable :: Ser_name(:)
end type 



end module