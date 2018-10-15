c
c	Array sizes for neighbourhood algorithm routines
c
c	PARAMETER		MEANING
c                                       The NA routines use the
c                                       following include file to define
c                                       all parameters controlling
c                                       memory required by the arrays.
c
c                                       The following is a description
c                                       of each parameter in na_param.in
c                                       that you are allowed to change to suit
c                                       your application. They are
c                                       checked against input parameters
c                                       read in from file na.in and an 
c                                       error is reported if any values
c                                       are too small.
c                                       
c       PARAMETER               MEANING
c
c       nd_max                  Maximum number of variables, i.e. the 
c                               dimension of the parameter space.
c
c       nit_max                 Maximum number of iterations
c
c       nsample_max             Maximum number of models generated
c                               in each sample
c
c       nh_max                  Maximum size of header of NAD output file
c                               (leave unchanged if you choose not to add
c                                any material to the NAD header. See manual)
c
c	maxseq			Maximum number of random sequences generated
c				by quasi random number generator. This
c				value can be set to 1 if the quasi
c				random number generator is not used.
c
c-----------------------------------------------------------------------
c
c				The following parameters are fixed
c				and should not be changed.
c
c	nmod_max		Maximum number of models to be generated
c			        (determined by nsample_max and nit_max)
c
c	nsleep_max		Maximum number of samples skipped over 
c				(Currently fixed at 1. Do not change)
c
c-----------------------------------------------------------------------
c
        parameter       (nsample_max=1000,
     &                   nit_max=5000,
     &                   nd_max=25,
     &                   nh_max=1000,
     &                   nsleep_max=1,
c    &                   maxseq=nd_max*nsample_max,
     &                   maxseq=50,
     &                   nmod_max = nsample_max*(nit_max+1))
