/* model.h
	 Definitions and prototypes for model.c
 */

#ifndef MODEL_H
#define MODEL_H

#define STRICT

/* Useful constants */
#define PI 3.1415926535
#define	SECS_PER_MA		3.1556736e13
#define	KELVINS_AT_0C	273.15
#define	U238YR	1.55125e-10
#define	U238MA	1.55125e-4
#define U238SEC 4.91575e-18
#define SQRT2PI	2.50662827463

/* Annealing models */
#define NUM_ANNEALING_MODELS 5
#define KETCHAM_ET_AL	0
#define	DAI       	  1
#define LASLETT_DUR   2
#define CROWLEY_DUR   3
#define CROWLEY_F_AP  4

/* Kinetic parameters supported */
#define NUM_KINETIC_TYPES 4
#define	ETCH_PIT_LENGTH		0
#define CL_PFU						1
#define OH_PFU						2
#define	CL_WT_PCT					3

/* Ways to treat initial track length */
#define NUM_L0_MODELS     2
#define L0_FROM_KPAR      0
#define L0_FROM_USER      1
#define L0_USER_MIN       14.0
#define L0_USER_MAX       20.0
#define L0_DEF_SIZE       NUM_KINETIC_TYPES * 2

/* Options for age to use for fitting */
#define NUM_AGE_TYPES		2
#define POOLED					0
#define WEIGHTED_MEAN		1
#define	CENTRAL					2

/* Options for length statistic */
#define NUM_LEN_STAT_TYPES	2
#define KS_TEST	0
#define KUIPER	1

/* Bounds for standard length reduction */
#define STD_LEN_RED_MIN		0.85
#define STD_LEN_RED_MAX		1.0

/* Ways to calculated kinetic value for modeling data */
#define DP_MIDPT  0
#define DP_WTMEAN 1

/* Array sizes, etc. */
#define	PDF_NUMSD					4
#define	PDF_ELTSPERSTDEV 	200
#define	PDF_ARRAYSIZE 		PDF_NUMSD*PDF_ELTSPERSTDEV+1

#define NUM_PDF_PTS    20		/* 81 */
#define PDF_BIN_SIZE	1.0  /*0.25 */

/*/#define	MAX_NUM_TIME_STEPS 170 */
#define	MAX_NUM_TIME_STEPS 5000
#define	MAX_NUM_DATA			 250
#define	MAX_NUM_AGES			 100

typedef struct  ageType {
	int     numSpon;  /* Number of spontaneous tracks */
	int     numInd;   /* Number of induced tracks */
	double  age;      /* Single-grain age (Ma) */
	double  ageErr;   /* Single-grain age error (1 standard deviation) */
  double  kPar;     /* Kinetic parameter value for grain */
} ageRec, *agePtr;

typedef struct	tlType {
	double	length;   /* Fission track length (microns) */
  double  pLength;  /* Track length projected to c-axis (microns) */
	double	kPar;     /* Kinetic parameter of grain in which track measured */
} tlRec, *tlPtr;

typedef struct  ttPathType {
	double  time;     /* Time of nodal point (in Ma for spec, sec for model) */
	double  temp;     /* Temperature at point (in deg C for spec, K for model) */
} ttPathRec, *ttPathPtr;
typedef ttPathRec *PttPath;

/* StatsType -- For track length distributions */
typedef struct statsType {
	double	mean,stDev,skewness,kurtosis,stdError;
} statsRec, *statsPtr;

typedef struct annealModelType {
	double c0,c1,c2,c3,a,b,lMin;
} annealModelRec, *annealModelPtr;

typedef struct annealConvType {
	double kPar,a0,a1;
} annealConvRec, *annealConvPtr;

typedef struct l0DefType {
	double m,b;
} l0Def, *l0DefPtr;

typedef struct tableType {
	double par,value;
} table, *tablePtr;

/* Prototypes */
void	InitPDFAxisSimple(double pdfAxis[],int *numPDFPts);

int  ForwardModel(
                    int        numTTDefs,         /* Number of nodes in t-T path definition */
                    double    stdLengthReduction, /* Length reduction in age standard */
                    double    kinPar,             /* Kinetic parameter value */
                    double    pctPerTimeStep,     /* Maximum % of model duration per interpolated time step */
                    int       annealModel,        /* Annealing model (code) */
                    int       doProject,          /* Do c-axis projection (1=yes, 0=no) */
                    int       usedCf,             /* Was Cf-irradiation used for analysis? (1=yes, 0=no) */
                    int       kinParType,         /* Kinetic parameter type (code) */
                    int       l0model,            /* Initial length model used (code) */
                    double    l0user,             /* User-defined initial length */
                    double    cdf[],              /* Calculated track length cdf */
                    int       numPDFPts,          /* Number of points in pdf and cdf arrays */
                    double    pdfAxis[],          /* PDF and CDF x-axis values */
                    double    pdf[],              /* Calculated track length pdf */
                    double    *oldestModelAge,    /* Calculated age of oldest fission track observed */
                    double    *ftModelAge,        /* Calculated fission-track age */
                    int       *numPopulations,   /* Number of interpolated time steps used in which f-t lengths were non-zero */
                    ttPathPtr    tTDef);             /* Time-temperature path definition, youngest age first */

#endif
