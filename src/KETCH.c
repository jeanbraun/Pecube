/* model.c
   Implementation of fission track model in C.
 */

#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "ketch.h"
int kinpar;
double kinetic_par;
int l0model;


float main_he(int ntime,float He_time[],float He_temp[]);
/* Variables declared globally */

/* Annealing model parameters */
  annealModelRec modKetchamEtAl = {-19.844,0.38951,-51.253,-7.6423,-0.12327,-11.988,0.0};

  annealModelRec modTILm = {-1.66965,0.0000241755,-12.4864,0.000843004,0.675508,4.16615,0.0};
  annealModelRec modTILc = {-2.36910,0.0000603834,-8.65794,0.000972676,0.404700,1.65355,9.0};
/* This is Crowley Durango model from Willett, 1992, in turn from where???
  annealModelRec modCrowDur = {-2.959,0.00008336,-17.0585,0.0005382,0.51,2.97,0.0}; 
*/
  annealModelRec modCrowDur = {-3.202,0.00009367,-19.6328,0.0004200,0.49,3.00,0.0};
  annealModelRec modCrowFAp = {-1.508,0.00002076,-10.3227,0.0009967,0.76,4.30,0.0};
  annealModelRec modCrowSrAp = {-1.123,0.00001055,-5.0085,0.001195,0.97,4.16,0.0};
  annealModelRec modLasDur  = {-4.87,0.000168,-28.12,0.0,0.35,2.7,0.0};
  int numConvDAILmDpar = 3;
  annealConvRec convDAILmDpar[3] = {{4.58,0.00000,1.00000},  /* B2 */
                               {2.43,0.5750,0.38374},  /* Fish Canyon */
										 {1.65,0.78477,0.21224}}; /* Renfrew */
  int numConvDAILmCl = 3;
  annealConvRec convDAILmCl[3] = {{2.95,0.00000,1.00000},  /* B2 */
                               {0.81,0.5750,0.38374},  /* Fish Canyon */
                             {0.03,0.78477,0.21224}}; /* Renfrew */

/* Initial track length (l0) function definition array */
  l0Def l0vsKparDefault[L0_DEF_SIZE] = {{0.283,15.63},  /* dpar mean */
                                 {0.205,16.10},  /* dpar projected */
                                 {0.544,16.18},  /* Cl(pfu) mean */
                                 {0.407,16.49},  /* Cl(pfu) projected */
                                 {0.000,16.18},  /* OH(pfu) mean */
                                 {0.000,16.57},  /* OH(pfu) projected */
                                 {0.13824,16.288},  /* Cl(wt%) mean used only for DAI legacy model */
                                 {0.17317,16.495}}; /* Cl(wt%) proj used only for DAI legacy model */

  l0Def l0vsKpar[L0_DEF_SIZE];

 
/* Tables */
#define NUM_RMEAN_TO_RCPAR 18
  table   meanToCpar[NUM_RMEAN_TO_RCPAR] = {{1.0,1.0},{0.95881,0.970588},{0.917508,0.941176},
                                            {0.876082,0.911765},{0.834515,0.882353},
                                            {0.792789,0.852941},{0.750883,0.823529},
                                            {0.708769,0.794118},{0.666416,0.764706},
                                            {0.644024,0.735294},{0.592261,0.705882},
                                            {0.543816,0.676471},{0.499281,0.647059},
                                            {0.460204,0.617647},{0.428183,0.588235},
                                            {0.403557,0.558824},{0.385398,0.529412},
                                            {0.0,0.0}};
#define NUM_RMEAN_TO_RDEN 18
  table   lengthToDensity[NUM_RMEAN_TO_RDEN] = {{1.0,1.0},{0.95881,0.956529},{0.917508,0.912921},
                                            {0.876082,0.869161},{0.834515,0.82523},
                                            {0.792789,0.781104},{0.750883,0.736758},
                                            {0.708769,0.692159},{0.666416,0.647269},
                                            {0.644024,0.541049},{0.592261,0.442989},
                                            {0.543816,0.358767},{0.499281,0.286844},
                                            {0.460204,0.225362},{0.428183,0.172461},
                                            {0.403557,0.126281},{0.385398,0.085095},
                                            {0.0,0.0}};

/* ----------------------- Utility routines ------------------------ */

/*  InitialTrackLength
    Returns the initial track length for the population based on the apatite
    kinetics, using data from experiment H0 by W.D.Carlson and R.A.Donelick
    (UT Austin), or a constant value specified by the user.
 */
double  InitialTrackLength(double kineticPar,
                           int kineticParType,
                           int doProject,
                           int l0model,
                           double l0user)
{ int index ;
  if (l0model == L0_FROM_USER) return(l0user);
/* else... */   
    index = 2*kineticParType + doProject;
    /*  return l0vsKpar[index].m*kineticPar + l0vsKpar[index].b; KG */
  return l0vsKparDefault[index].m*kineticPar + l0vsKparDefault[index].b;
}


/* FindValue -- A routine that looks up a value in a table; the table
   is stored in a structure containing the lookup parameters and
   corresponding values.  The lookup parameters must be sorted in
   descending order.  If the parameter is not matched exactly, but is
   bounded by lookup values, a linear interpolation is used.  If the
   parameter is beyond either endpoint of the lookup values, the
   value returned is equal to the value for the appropriate endpoint.
 */
double FindValue(double par,tablePtr table,int numTable)
{
  int i;
  double frac;

  if (par >= table[0].par) return(table[0].value);
  if (par <= table[numTable-1].par) return(table[numTable-1].value);
  for (i=1;par < table[i].par;i++) ;
  frac = (par - table[i].par)/(table[i-1].par - table[i].par);
  return(table[i-1].value*frac + table[i].value*(1.0-frac));
}

/*  ObservationalBias
    Basically, the probability of a population being observed relative to
    the probability of the longest population of tracks.  For mean length
    models, this the probability is simply the reduced length (Line segment
    theory; Parker and Cowan, 1976; Laslett et al., 1982).

    Updated 9/13/99 -- Actually, this should take into account loss of some
    tracks to total annealing.  The closest thing we have so far is the
    relationship between observed length and observed density.
 */
double  ObservationalBias(double redLength)
{
  return(FindValue(redLength,lengthToDensity,NUM_RMEAN_TO_RDEN));
}


/* ---------------------- T-t Path Definition --------------------- */

/*  InterpolateTTPathKet
    Takes the time-temperature path specification and subdivides it for
    calculation in isothermal intervals.
    Does it based on model of Ketcham et al., in review.
    It is calibrated to facilitate 0.5% accuracy for end-member F-apatite by
    having a maximum temperature step of 3.5 degrees C when the model
    temperature is within 10 C of the total annealing temperature.  Before this
    cutoff the maximum temperature step required is 8 C.  If the overall model
    time steps are too large, these more distant requirements may not be met.
    
 */
#define NEAR_ANNEAL_CUTOFF_KET     10.0
#define MAX_TEMP_STEP_NEAR_TA_KET  3.5
#define MAX_TEMP_STEP_KET           8
int  InterpolateTTPathKet(int       numTTDefs,
							 ttPathPtr  tTDef,
							 int       *numTTNodes,
							 ttPathPtr tTPath,
							 double    pctPerTimeStep)
{
  int      dN,n,i;
  double  rate,absRate;   /* Rate of temperature change (K/m.y.) */
  double  maxTMult;       /* Used to help find max temp in a time step */
  double  maxTemp;        /* Maximum temperature for tiem step (K) */
  double  nearAnnealTemp; /* Temperature to limit temp steps at (K) */
  double  timeStep;       /* Size of individual time step (m.y.) */
  double  defTimeStep;    /* Overall default time step (m.y.) */
  double  tempPerTimeStep;/* Temperature change per default time step (K) */
  double  currDefTimeStep;/* Default time step for the current path segment (m.y.) */
  double  altTimeStep;    /* Alternative time step for high temps (m.y.) */
  double  endTemp;        /* Temperature at end of current t-T segment */



/* Initialize path, in case it starts out below maximum annealing temp */
  *numTTNodes = 1;
  altTimeStep=0.0;
  tTPath[0].temp = tTDef[numTTDefs-1].temp + KELVINS_AT_0C;
  tTPath[0].time = tTDef[numTTDefs-1].time;
  defTimeStep = tTDef[numTTDefs-1].time*pctPerTimeStep/100.0;
  for (dN=numTTDefs-1;dN>0;dN--) {
/* Calculate rate and total annealing temperature for this t-T segment */
	 rate = (tTDef[dN].temp-tTDef[dN-1].temp)/(tTDef[dN].time-tTDef[dN-1].time+0.0001);
	 absRate = fabs(rate);
	 tempPerTimeStep = absRate*defTimeStep;
	 currDefTimeStep = (tempPerTimeStep <= MAX_TEMP_STEP_KET) ? defTimeStep : MAX_TEMP_STEP_KET/absRate;
	 maxTMult = rate > 0 ? 0 : -1;
	 endTemp = tTDef[dN-1].temp + KELVINS_AT_0C;
/* Calculate Ta; if rate is low, set nearAnnealTemp arbitrarily high */
	 if (absRate < 0.1)
		nearAnnealTemp = 1000.0;
	 else {
		nearAnnealTemp = 3.7767*pow(absRate,0.019837) - NEAR_ANNEAL_CUTOFF_KET;
		altTimeStep = MAX_TEMP_STEP_NEAR_TA_KET/absRate;
	 }
	 while (tTPath[*numTTNodes-1].time > tTDef[dN-1].time) {
/* Make sure we haven't specified too many nodes */
		if (*numTTNodes + 1 > MAX_NUM_TIME_STEPS) return(0);
		maxTemp = tTPath[*numTTNodes-1].temp + defTimeStep*rate*maxTMult;
/* If heating, make sure maxTemp not higher than end of segment */
		if ((rate < 0) && (maxTemp > endTemp)) maxTemp = endTemp;
		timeStep = currDefTimeStep;
		if (maxTemp > nearAnnealTemp)
		  if (altTimeStep < defTimeStep) timeStep = altTimeStep;
/* Check to see if this is final step for this segment.  A small factor
	is added to account for the possibility of roundoff. NOTE: This factor must
	be significantly shorter than any time step. */
		if (timeStep+0.001 > tTPath[*numTTNodes-1].time - tTDef[dN-1].time) {
		  tTPath[*numTTNodes].time = tTDef[dN-1].time;
		  tTPath[*numTTNodes].temp = endTemp;
		  /* printf("TT %d %f %f %f %f \n",*numTTNodes,tTPath[*numTTNodes].temp,rate,timeStep,tTPath[*numTTNodes].time);  */
		}
		else {
		  tTPath[*numTTNodes].time = tTPath[*numTTNodes-1].time - timeStep;
		  tTPath[*numTTNodes].temp = tTPath[*numTTNodes-1].temp - rate*timeStep;
/*        printf("TR %d %f %f %f %f \n",*numTTNodes,tTPath[*numTTNodes].temp,rate,timeStep,tTPath[*numTTNodes].time);  */
		}
		(*numTTNodes)++;
	 }
  }
/* Convert Ma to seconds */
  for (n=0; n < *numTTNodes; n++)
	 tTPath[n].time *= SECS_PER_MA;




  return(1);
}


/*  InterpolateTTPath
	 Takes the time-temperature path specification and subdivides it for
	 calculation in isothermal intervals.
 */
#define B_1 435.34
#define B_2 0.013793

#define NEAR_ANNEAL_CUTOFF     6
#define MAX_TEMP_STEP_NEAR_TA  1
#define MAX_TEMP_STEP  5
int  InterpolateTTPath(int       numTTDefs,
							 ttPathPtr  tTDef,
							 int       *numTTNodes,
							 ttPathPtr tTPath,
							 double    pctPerTimeStep)
{
  int      dN,n;
  double  rate,absRate;   /* Rate of temperature change (K/m.y.) */
  double  maxTMult;       /* Used to help find max temp in a time step */
  double  maxTemp;        /* Maximum temperature for tiem step (K) */
  double  nearAnnealTemp; /* Temperature to limit temp steps at (K) */
  double  timeStep;       /* Size of individual time step (m.y.) */
  double  defTimeStep;    /* Default time step (m.y.) */
  double  altTimeStep=0.0;    /* Alternative time step for high temps (m.y.) */
  double  endTemp;        /* Temperature at end of current t-T segment */

/* Initialize path, in case it starts out below maximum annealing temp */
  *numTTNodes = 1;
  tTPath[0].temp = tTDef[numTTDefs-1].temp + KELVINS_AT_0C;
  tTPath[0].time = tTDef[numTTDefs-1].time;
  defTimeStep = tTDef[numTTDefs-1].time*pctPerTimeStep/100.0;
  for (dN=numTTDefs-1;dN>0;dN--) {
/* Calculate rate and total annealing temperature for this t-T segment */
    rate = (tTDef[dN].temp-tTDef[dN-1].temp)/(tTDef[dN].time-tTDef[dN-1].time);
    absRate = fabs(rate);
    maxTMult = rate > 0 ? 0 : -1;
    endTemp = tTDef[dN-1].temp + KELVINS_AT_0C;
/* Calculate Ta; if rate is low, set nearAnnealTemp arbitrarily high */
    if (absRate < 0.1)
		nearAnnealTemp = 1000.0;
	 else {
      nearAnnealTemp = (B_1*pow(absRate,B_2)) - NEAR_ANNEAL_CUTOFF;
      altTimeStep = MAX_TEMP_STEP_NEAR_TA/absRate;
/* The next line would be used if we weren't using the near-Ta cutoff */
/*     altTimeStep = (absRate > 0.1) ? MAX_TEMP_STEP/absRate : defTimeStep;      */
    }
    while (tTPath[*numTTNodes-1].time > tTDef[dN-1].time) {
/* Make sure we haven't specified too many nodes */
      if (*numTTNodes + 1 > MAX_NUM_TIME_STEPS) return(0);
      maxTemp = tTPath[*numTTNodes-1].temp + defTimeStep*rate*maxTMult;
/* If heating, make sure maxTemp not higher than end of segment */
      if ((rate < 0) && (maxTemp > endTemp)) maxTemp = endTemp;
		timeStep = defTimeStep;
      if (maxTemp > nearAnnealTemp)
        if (altTimeStep < defTimeStep) timeStep = altTimeStep;
/* Check to see if this is final step for this segment.  A small factor
   is added to account for the possibility of roundoff. NOTE: This factor must
   be significantly shorter than any time step. */
      if (timeStep+0.001 > tTPath[*numTTNodes-1].time - tTDef[dN-1].time) {
        tTPath[*numTTNodes].time = tTDef[dN-1].time;
        tTPath[*numTTNodes].temp = endTemp;
      }
      else {
        tTPath[*numTTNodes].time = tTPath[*numTTNodes-1].time - timeStep;
        tTPath[*numTTNodes].temp = tTPath[*numTTNodes-1].temp - rate*timeStep;
      }
      (*numTTNodes)++;
    }
  }
/* Convert Ma to seconds */
  for (n=0; n < *numTTNodes; n++)
    tTPath[n].time *= SECS_PER_MA;
  return(1);
}



/* --------------------------- Modeling --------------------------- */

/*  InitPDFAxisSimple
    Compiles the x-axis values for which the pdf and cdf functions will be
    calculated.  Currently assumes a uniform, PDF_BIN_SIZE micron spacing. 
 */
void  InitPDFAxisSimple(double pdfAxis[],int *numPDFPts)
{
  int  i;

  *numPDFPts = NUM_PDF_PTS;
  for (i=0; i<NUM_PDF_PTS; i++) pdfAxis[i] = i*PDF_BIN_SIZE;
}


/*  ReducedStdev
	 Calculates the reduced standard deviation of a track population length
	 from the reduced mean length.  Based on Carlson and Donelick (unpub.
    data).
 */
double  ReducedStdev(double redLength,int doProject)
{
  if (doProject) return(0.1081-0.1642*redLength+0.1052*redLength*redLength);
  else return(0.4572-0.8815*redLength+0.4947*redLength*redLength);
}


/*  ChooseAnnealingModel
    Selects which annealng model the user has specified, and arranges the
	 appropriate constants.
 */
void  ChooseAnnealingModel(int    annealingModel,
                           int    kinParType,
                           int    *numConv,
                           annealModelRec *annMod,
                           annealConvPtr *annConv)
{
  switch (annealingModel) {
    case DAI:
      *annMod = modTILm;
      if (kinParType == ETCH_PIT_LENGTH) {
        *numConv = numConvDAILmDpar;
        *annConv = convDAILmDpar;
      }
      else {  /* kinParType = CL_WT_PCT */
        *numConv = numConvDAILmCl;
        *annConv = convDAILmCl;
      }
      break;
    case LASLETT_DUR:
      *annMod = modLasDur;
      *numConv = 0;
      break;
    case CROWLEY_DUR:
      *annMod = modCrowDur;
      *numConv = 0;
      break;
    case CROWLEY_F_AP:
      *annMod = modCrowFAp;
      *numConv = 0;
      break;
  }
}

/* CalcModelLengthsKet
   Calculates the model track length distribution for a given time=
   temperature history based on the calibration of Ketcham et al. (in
   review).  The length calculated is the reduced modeled mean
   c-axis-parallel length (Rcmod).
 */                                               
#define MIN_OBS_RCMOD  0.55
void CalcModelLengthsKet(ttPathPtr    tTPath,
                         int          numTTNodes,
                         double        redLength[],
                         double        kinPar,
								 int          kinParType,
                         int          *firstTTNode)
{
  int      node, nodeB;
  double  equivTime;
  double  timeInt,x1,x2,x3;
  double  totAnnealLen;
  double  equivTotAnnLen;
  double  rmr0,k;  /* Apatite-apatite conversion factorsq */
  double  calc;
  double  tempCalc;
  rmr0 = 0;
/* Calculate the rmr0-k values for the kinetic parameter given */
  switch (kinParType) {
    case ETCH_PIT_LENGTH:
      if (kinPar <= 1.75) rmr0 = 0.84;
      else if (kinPar >= 4.58) rmr0 = 0.0;
      else rmr0 = 1.0-exp(0.647*(kinPar-1.75)-1.834);
      break;
    case CL_WT_PCT:
      /* Just convert the kinetic parameter to Cl apfu
         Note that this invalidates kinPar for the rest of the routine */
      kinPar = kinPar * 0.2978;
    case CL_PFU:
      calc = fabs(kinPar-1.0);
      if (calc <= 0.130) rmr0 = 0.0;
      else rmr0 = 1.0-exp(2.107*(1.0-calc)-1.834);
      break;
    case OH_PFU:
      calc = fabs(kinPar-1.0);
      rmr0 = 0.84*(1.0-pow(1.0-calc,4.5));
      break;
  }
  k = 1-rmr0;

  totAnnealLen = MIN_OBS_RCMOD;
/*/ equivTotAnnLen is the length of the more resistant apatite at the length of
//  total annealing for the less resistant apatite we're modeling.
//  In the future, if this routine is adapted to solve for many different apatite
//  kinetic populations at once, we would use the rmr0 and k values for the most
//  resistant apatite being modeled.
*/
  equivTotAnnLen = pow(totAnnealLen,1.0/k)*(1.0-rmr0)+rmr0;

  equivTime = 0.0;
  tempCalc = log(1.0/((tTPath[numTTNodes-2].temp + tTPath[numTTNodes-1].temp)/2.0));
  for (node = numTTNodes-2; node >= 0; node--) {
    timeInt = tTPath[node].time - tTPath[node+1].time + equivTime;

	 x1 = (log(timeInt) - modKetchamEtAl.c2)/(tempCalc - modKetchamEtAl.c3);
    x2 = 1.0 + modKetchamEtAl.a * (modKetchamEtAl.c0 + modKetchamEtAl.c1 * x1);
    redLength[node] = pow(x2,1.0/modKetchamEtAl.a);
    x3 = 1.0 - modKetchamEtAl.b * redLength[node];
    redLength[node] = (x3 < 0) ? 0.0 : pow(x3, 1.0/modKetchamEtAl.b);

    if (redLength[node] < equivTotAnnLen)
      redLength[node] = 0.0;

/* Check to see if we've reached the end of the length distribution
   If so, we then do the kinetic conversion. */
    if ((redLength[node] == 0.0) || (node == 0)) {
      *firstTTNode = (node ? node+1 : node);
      for (nodeB = *firstTTNode; nodeB < numTTNodes-1; nodeB++) {
          if (redLength[nodeB] <= rmr0) {
            redLength[nodeB] = 0.0;
            *firstTTNode = nodeB;
          }
          else {
            redLength[nodeB] = pow((redLength[nodeB] - rmr0)/(1.0 - rmr0),k);
            if (redLength[nodeB] < totAnnealLen) {
              redLength[nodeB] = 0.0;
              *firstTTNode = nodeB;
            }
          }
      }
      return;
    }

/* Update tiq for this time step */
    if (redLength[node] < 0.999) {
      tempCalc = log(1.0/((tTPath[node-1].temp + tTPath[node].temp)/2.0));
      equivTime = pow((1.0-pow(redLength[node],modKetchamEtAl.b))/modKetchamEtAl.b,modKetchamEtAl.a);
      equivTime = ((equivTime - 1.0)/modKetchamEtAl.a - modKetchamEtAl.c0)/modKetchamEtAl.c1;
      equivTime = exp(equivTime*(tempCalc-modKetchamEtAl.c3)+modKetchamEtAl.c2);
    }
  }  
}

/*  CalcModelLengths
    Calculates the model track length distribution for the given time-
    temperature history.  For each T-t segment, it finds the reduced
    mean and standard deviation for the population of track lengths based
    on the model of Laslett et al. (1987).  The tiq calculation uses
    Goswami et al. (1984) and Duddy et al. (1988).
 */
#define MIN_OBS_RM  0.4095
void  CalcModelLengths(ttPathPtr      tTPath,
                       int            numTTNodes,
                       double          redLength[],
                       double         kinPar,
                       annealModelRec annMod,
                       annealConvPtr  annConv,
                       int            numConv,
							  int            *firstTTNode)
{
  int      node, nodeB;
  double  deltaTimes[MAX_NUM_TIME_STEPS];    
  double  meanTemps[MAX_NUM_TIME_STEPS];     
  double  equivTime;
  double  timeInt,x1,x2,x3;
  double  totAnnealLen;
  double  frac;
  double  a0=0.0;
  double  a1 = 0.0;
  int     c0;

/* Find position along kinetic line, calculate relative zero point */
  if (numConv) {
    if (kinPar >= annConv[0].kPar) c0 = -1;
    else for (c0=numConv-1;c0 && (kinPar > annConv[c0].kPar);c0--) ;
    if (c0 == -1) {
      a0 = annConv[0].a0;
      a1 = annConv[0].a1;
    }
    else if (c0 == numConv-1) {
      a0 = annConv[numConv-1].a0;
      a1 = annConv[numConv-1].a1;
    }
    else {
      frac = (kinPar-annConv[c0+1].kPar)/(annConv[c0].kPar-annConv[c0+1].kPar);
      a0 = annConv[c0].a0*frac + annConv[c0+1].a0*(1-frac);
      a1 = annConv[c0].a1*frac + annConv[c0+1].a1*(1-frac);
    }
  }
  totAnnealLen = MIN_OBS_RM;
  for (node = 0; node < numTTNodes-1; node++) {
    deltaTimes[node] = tTPath[node].time - tTPath[node+1].time;  
    meanTemps[node] = (tTPath[node].temp + tTPath[node+1].temp)/2.0;
  }

  equivTime = 0.0;
  for (node = numTTNodes-2; node >= 0; node--) {
    timeInt = deltaTimes[node] + equivTime;
    x1 = (log(timeInt) - annMod.c2)/(1.0 / meanTemps[node] - annMod.c3);
    x2 = 1.0 + annMod.a * (annMod.c0 + annMod.c1 * x1);
    redLength[node] = pow(x2,1.0/annMod.a);
    x3 = 1.0 - annMod.b * redLength[node];
    redLength[node] = (x3 < 0) ? 0.0 : pow(x3, 1.0/annMod.b);
    if (redLength[node] < totAnnealLen)
      redLength[node] = 0.0;

/* Check to see if we've reached the end of the length distribution
   If so, we then do the kinetic conversion. */
    if ((redLength[node] == 0.0) || (node == 0)) {
      *firstTTNode = (node ? node+1 : node);
//		printf(" FIRSAT 2 = %d\n",*firstTTNode);
		if (numConv)
		  for (nodeB = *firstTTNode; nodeB < numTTNodes-1; nodeB++) {
			 if (redLength[nodeB] <= a0) {
				redLength[nodeB] = 0.0;
            *firstTTNode = nodeB;
//				printf("FIRAZT 3 \n",*firstTTNode );
			 }
          else {
            redLength[nodeB] = pow((redLength[nodeB] - a0)/(1 - a0),a1);
            if (redLength[nodeB] < totAnnealLen) redLength[nodeB] = 0.0;
          }
        }

      return;
    }

/* Update tiq for this time step */
    if (redLength[node] < 0.999) {
      equivTime = pow((1.0-pow(redLength[node],annMod.b))/annMod.b,annMod.a);
      equivTime = ((equivTime - 1.0)/annMod.a - annMod.c0)/annMod.c1;
      equivTime = exp(equivTime*(1.0/meanTemps[node-1]-annMod.c3)+annMod.c2);
    }
  }  
}


/*  AgeCorrectionKet
    Does the conversion from length to density for the Ketcham et al., 1999 model.
    The routine is placed "way up here" because it will also be used to estimate
    bias for population summing.

    Assumes we're passing in a c-axis-projected length
 */
double AgeCorrectionKet(double cparlen)
{
  if (cparlen > 0.757) return(1.600*cparlen-0.599);
  if (cparlen >= MIN_OBS_RCMOD) return(9.205*cparlen*cparlen-9.157*cparlen+2.269);
  return(0.0);
}



/*  SumPopulationsKet
    Sums the individual model track length populations into an overall
    population, and normalizes.  Takes care of conversion from projected to mean
    lengths and finding the standard deviation of the population distribution.
 */
/* MIN_LENGTH -- the minimum observable length */
#define MIN_LENGTH 2.15

void  SumPopulationsKet(int    numPDFPts,
                      int      numTTNodes,
                      int      firstTTNode,
                      int     doProject,
                      int      usedCf,
							 ttPathPtr tTPath,
                      double  pdfAxis[],
                      double  pdf[],
                      double  cdf[],
                      double  initLength,
							 double  redLength[])
{
  int      i,j;
  double  weight,rLen,rStDev,obsBias,rmLen,calc,z;
  double   wt1,wt2;

/* Sum curves for pdf */
  for (i=0; i < numPDFPts; i++) pdf[i] = 0.0;

  wt1 = exp(U238SEC*tTPath[firstTTNode].time)/U238SEC;
  for (j=firstTTNode; j < numTTNodes-1; j++) {
    wt2 = exp(U238SEC*tTPath[j+1].time)/U238SEC;
    weight = wt1-wt2;
    wt1 = wt2;
    rmLen = usedCf ? 1.396*redLength[j]-0.4017 : -1.499*redLength[j]*redLength[j]+4.150*redLength[j]-1.656;
    rLen = doProject ? redLength[j] : rmLen;
    rStDev = ReducedStdev(rLen,doProject);
    obsBias = AgeCorrectionKet(redLength[j]);
    calc = weight*obsBias/(rStDev*SQRT2PI);
    if (rLen > 0) {
      for (i=0; i < numPDFPts; i++) {
        if (pdfAxis[i] >= MIN_LENGTH) {
          z = (rLen-pdfAxis[i]/initLength)/rStDev;
          if (z <= PDF_NUMSD) pdf[i] += calc*exp(-(z*z)/2.0);
        }
      }
    }
  }
/* Calculate cdfs. */
  cdf[0] = pdf[0];
  for (i=1; i < numPDFPts; i++)
    cdf[i] = cdf[i-1]+((pdf[i]+pdf[i-1])/2.0)*(pdfAxis[i]-pdfAxis[i-1]);
/* Normalize */
  if (cdf[numPDFPts-1] > 0.0)  /* Some non-zero lengths */
    for (i=0; i < numPDFPts; i++)  {
      pdf[i] = pdf[i]/cdf[numPDFPts-1];
      cdf[i] = cdf[i]/cdf[numPDFPts-1];
    }
}


/*  SumPopulations
    Sums the individual model track length populations into an overall
    population, and normalizes.  Takes care of conversion from mean to projected
    lengths and finding the standard deviation of the population distribution.
 */
void  SumPopulations(  int      numPDFPts,
                      int      numTTNodes,
                      int      firstTTNode,
							 int     doProject,
                      ttPathPtr tTPath,
                      double  pdfAxis[],
                      double  pdf[],
                      double  cdf[],
                      double  initLength,
                      double  redLength[])
{
  int      i,j;
  double  weight,rLen,rStDev,obsBias,calc,z;
  double  wt1, wt2;

  for (i=0; i < numPDFPts; i++) pdf[i] = 0.0;
/* Sum curves for pdf */
  wt1 = exp(U238SEC*tTPath[firstTTNode].time)/U238SEC;
  for (j=firstTTNode; j < numTTNodes-1; j++) {
    wt2 = exp(U238SEC*tTPath[j+1].time)/U238SEC;
    weight = wt1-wt2;
    wt1 = wt2;
    rLen = doProject ? FindValue(redLength[j],meanToCpar,NUM_RMEAN_TO_RCPAR) : redLength[j];
    rStDev = ReducedStdev(redLength[j],doProject);
    obsBias = ObservationalBias(redLength[j]);
    calc = weight*obsBias/(rStDev*SQRT2PI);
    if (rLen > 0) {
      for (i=0; i < numPDFPts; i++) {
        if (pdfAxis[i] >= MIN_LENGTH) {
          z = (rLen-pdfAxis[i]/initLength)/rStDev;
            if (z <= PDF_NUMSD) pdf[i] += calc*exp(-(z*z)/2.0);
        }
      }
    }
  }
/* Sum the cdf. */
  cdf[0] = pdf[0];
  for (i=1; i < numPDFPts; i++)
    cdf[i] = cdf[i-1]+((pdf[i]+pdf[i-1])/2.0)*(pdfAxis[i]-pdfAxis[i-1]);
/* Normalize */
  if (cdf[numPDFPts-1] > 0.0)  /* Some non-zero lengths */
    for (i=0; i < numPDFPts; i++)  {
      pdf[i] = pdf[i]/cdf[numPDFPts-1];
      cdf[i] = cdf[i]/cdf[numPDFPts-1];
    }
}


/*  CalcModelStats
    Calculates the descriptive statistics of the model track length
    distribution.
 */
void  CalcModelStats( int       numPopulations,
                      int       numPDFPts,
                      double    pdfAxis[],
                      double    pdf[],
                      statsPtr  stats)
{
  double  areas[MAX_NUM_DATA];
  double  sumAreas;
  int      i;
  double  weight;

  stats->stDev = 0.0;
  stats->skewness = 0.0;
  stats->kurtosis = 0.0;

  if (numPopulations == 0) {
    stats->mean = 0.0;
    stats->stdError = 0.0;
    return;
  }
  areas[0] = 0.0;
  sumAreas = areas[0];
  stats->mean = 0.0;

  for (i=1; i<numPDFPts; i++) {
    areas[i] = ((pdf[i]+pdf[i-1])/2.0)*(pdfAxis[i]-pdfAxis[i-1]);
    sumAreas += areas[i];
    stats->mean += areas[i]*(pdfAxis[i]+pdfAxis[i-1])/2.0;
  }
  if (sumAreas == 0.0) {
    stats->stdError = 0.0;
    return;
  }
  stats->mean /= sumAreas;

  for (i=0; i<numPDFPts; i++)
    stats->stDev += (areas[i]/sumAreas)*pow(pdfAxis[i] - stats->mean,2.0);

  stats->stDev = sqrt(stats->stDev);

  for (i=0; i<numPDFPts; i++) {
    weight = areas[i]/sumAreas;
    stats->skewness += weight*pow((pdfAxis[i] - stats->mean)/stats->stDev,3.0);
    stats->kurtosis += weight*pow((pdfAxis[i] - stats->mean)/stats->stDev,4.0);
  }

  stats->kurtosis -= 3.0;

  stats->stdError = stats->stDev / sqrt(numPopulations);
}


/*  CalcModelAgesKet
    Calculates the estimated age which would be measured from the model
    fission track population.  Each time interval is added in, and
    corrected by the amount of track length reduction (causing the age
    to appear smaller).
    This version adapted to the model of Ketcham et al., 1999,
    assuming we're using c-axis-parallel lengths.

    9/17/99 Updated to make it use the midpoint length during a time step,
    rather than the endpoint one, for calculating the age correction.
 */
void  CalcModelAgesKet(ttPathPtr  tTPath,
                    double  redLength[],
                    int      numTTNodes,
                    int      firstNode,
                    double  *oldestModelAge,
                    double  *ftModelAge,
                    double  stdLengthReduction)
{
  int     node;
  double  midLength;

  *oldestModelAge = tTPath[firstNode].time/SECS_PER_MA;
  for (*ftModelAge=0.0, node=firstNode; node < numTTNodes-2; node++) {
/* Correct each time interval for length reduction */
    midLength = (redLength[node]+redLength[node+1])/2.0;
    *ftModelAge += AgeCorrectionKet(midLength)*(tTPath[node].time-tTPath[node+1].time);
  }
  *ftModelAge += AgeCorrectionKet(redLength[numTTNodes-2])*(tTPath[node].time-tTPath[node+1].time);

/* Account for length reduction in length standard, convert to Ma */
  *ftModelAge /= (stdLengthReduction*SECS_PER_MA);
}


/*  AgeCorrection
    Estimates the correction in the fission track age caused by length
    reduction over a time interval.  If orientation is ignored, the
    appropriate answer is to follow the relationship of track length
    reduction to track density reduction (Green, 1988; Willett, 1992).
    If we project track lengths, we can simply use the reduced track
    length.
 */
double AgeCorrection(double redLength)
{
  return(FindValue(redLength,lengthToDensity,NUM_RMEAN_TO_RDEN));
}


/*  CalcModelAges
    Calculates the estimated age which would be measured from the model
    fission track population.  Each time interval added in, and
    corrected by the amount of density reduction expected given the
    length reduction.

    9/17/99 Updated to use mid-time-step length, rather than end-time-step,
    for calculating expected density reduction.
 */
void  CalcModelAges(ttPathPtr  tTPath,
                    double  redLength[],
                    int      numTTNodes,
						  int      firstNode,
                    double  *oldestModelAge,
                    double  *ftModelAge,
                    double  stdLengthReduction)
{
  int     node;
  double  midLength;

  *oldestModelAge = tTPath[firstNode].time/SECS_PER_MA;
  for (*ftModelAge=0.0, node=firstNode; node < numTTNodes-2; node++) {
/* Correct each time interval for length reduction */
    midLength = (redLength[node]+redLength[node+1])/2.0;
    *ftModelAge += AgeCorrection(midLength)*(tTPath[node].time-tTPath[node+1].time);
  }
  *ftModelAge += AgeCorrection(redLength[numTTNodes-2])*(tTPath[node].time-tTPath[node+1].time);
/* Account for length reduction in length standard, convert to Ma */
  *ftModelAge /= (stdLengthReduction*SECS_PER_MA);
}


/*  ForwardModel
    This is the basic forward model routine.  Input consists of a time-
    temperature path, some kinetically influenced parameters, and the
    lengths at which pdf's should be calculated.  Output consists of
    the resulting pdf and cdf, track length distribution statistics, and
    fission-track age.
 */
 int  ForwardModel( int       numTTDefs,
                    double    stdLengthReduction,
                    double    kinPar,
                    double    pctPerTimeStep,
                    int       annealModel,
                    int       doProject,
                    int       usedCf,
                    int       kinParType,
                    int       l0model,
                    double    l0user,
                    double    cdf[],
                    int        numPDFPts,
                    double    pdfAxis[],
                    double    pdf[],
                    double    *oldestModelAge,
                    double    *ftModelAge,
                    int       *numPopulations,
                    ttPathPtr    tTDef
                    )
{
  double  redLength[MAX_NUM_TIME_STEPS]; 
  double  initLength;
  int     firstTTNode,numTTNodes;
  annealModelRec annMod;
  annealConvPtr  annConv;
  int     numConv;
  ttPathRec  tTPath[MAX_NUM_TIME_STEPS];  /* Interpolated time-temperature path */

  initLength = InitialTrackLength(kinPar,kinParType,doProject,l0model,l0user);
  /* KG */ if(initLength < 1) initLength = l0user;
//printf("initLenght= %f",initLength);
  if (annealModel == KETCHAM_ET_AL) {
    if (InterpolateTTPathKet(numTTDefs,tTDef,&numTTNodes,tTPath,pctPerTimeStep)) {
      CalcModelLengthsKet(tTPath,numTTNodes,redLength,kinPar,kinParType,&firstTTNode);
      SumPopulationsKet(numPDFPts,numTTNodes,firstTTNode,doProject,usedCf,tTPath,pdfAxis,
                      pdf,cdf,initLength,redLength);  
      CalcModelAgesKet(&tTPath[0],redLength,numTTNodes,firstTTNode,oldestModelAge,
							 ftModelAge,stdLengthReduction);
		*numPopulations = numTTNodes - firstTTNode;
	 }
	 else
		*numPopulations = 0;
  }
  else {
	 if (InterpolateTTPath(numTTDefs,tTDef,&numTTNodes,tTPath,pctPerTimeStep)) {
		ChooseAnnealingModel(annealModel,kinParType,&numConv,&annMod,&annConv);

		CalcModelLengths(tTPath,numTTNodes,redLength,kinPar,annMod,annConv,numConv,&firstTTNode);
		SumPopulations(numPDFPts,numTTNodes,firstTTNode,doProject,tTPath,pdfAxis,
							 pdf,cdf,initLength,redLength);
		CalcModelAges(tTPath,redLength,numTTNodes,firstTTNode,oldestModelAge,
							 ftModelAge,stdLengthReduction);
		*numPopulations = numTTNodes - firstTTNode;
	 }
	 else
		*numPopulations = 0;
  }
  return (*numPopulations);
}




void ketch_main_(int *ntime,float ketchtime[],float ketchtemp[],double *alo,double *final_age,double *oldest_age,double *fmean,double fdist[])
{
 // ttPathRec kerryTt[22];
 ttPathRec kerryTt[*ntime];
#define nbins  200
  int numPDFPts = nbins;
  double pdfAxis[nbins];
  double cdf[nbins];
  double pdf[nbins];
  double oldestModelAge;
  double ftModelAge;
  double meanlength;
  double lengthred;
  int numPopulations;
  int i,ktime;
  int annealmodel,doProject,usedCf;
//  int l0model;
  double lgeol = 14.56;
//  double lgeol;
//  lgeol=14.56;
  FILE *fg;
  

  for(i=0;i<numPDFPts;i++) pdfAxis[i] = (double)(i*1.0+0.5)*20.0/numPDFPts;

//for (i=0;i<*ntime;i++)printf("time( %i )=%f and temp( %i )=%f \n",i,ketchtime[i],i,ketchtemp[i]);
fg = fopen("Tt_test.dat","wt");
for (i=0;i<*ntime;i++) fprintf(fg,"%f %f \n",ketchtime[i],ketchtemp[i]);
fclose(fg);



  ktime = *ntime;
  for(i=0;i<*ntime;i++)  {
	  kerryTt[i].time = ketchtime[i];//+0.001*i; /*ad hoc modification to ensure times are different....) */
	  kerryTt[i].temp = ketchtemp[i];
//	  printf(" KT %3d %6.2f %6.2f \n",i,ketchtime[i],ketchtemp[i]);
  }

  lengthred = lgeol/(*alo);
  annealmodel = KETCHAM_ET_AL; /* LASLETT_DUR */
//  kinpar = ETCH_PIT_LENGTH;
//  kinetic_par=1.5;
  //over ride
kinpar = CL_WT_PCT;
kinetic_par = 0.0; //value of cl. weight pc
//  l0model = L0_FROM_KPAR;
l0model = L0_FROM_USER; /*mod VKP - made this active again*/

  doProject = 0;
  usedCf = 0;

//  lengthred = lgeol/alo;

//  printf("before forward...\n");

//  printf("ktime = %3d\n",ktime);
//  printf("lgeol = %8.3f\n",lgeol);
//  printf("alo = %8.3f\n",*alo);
//  printf("lengthred = %8.3f\n",lengthred);
//  printf("kinpar = %3d\n",kinpar);
//  printf("kinetic_par = %f\n",kinetic_par);
//  printf("annealmodel = %3d\n",annealmodel);
//  printf("doProject = %3d\n",doProject);
//  printf("usedCf = %3d\n",usedCf);

  ForwardModel(ktime,lengthred,kinetic_par,(double)1.0,annealmodel,doProject,usedCf,kinpar,l0model,*alo,cdf,numPDFPts,pdfAxis,pdf,&oldestModelAge,&ftModelAge,&numPopulations,kerryTt);

  meanlength = 0.0;
  for(i=0;i<numPDFPts;i++) {
		meanlength += pdfAxis[i]*pdf[i];
		fdist[i+1] = pdf[i];
  } 
  meanlength *=20.0/(double)numPDFPts;
  *fmean = meanlength ;
  *final_age = ftModelAge;
  *oldest_age = oldestModelAge;
  /* printf(" OUT 1.5  %f %f\n",ftModelAge,meanlength);  */

  return ;


}
