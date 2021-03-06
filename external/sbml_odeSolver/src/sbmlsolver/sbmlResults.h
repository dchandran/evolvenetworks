/*
  Last changed Time-stamp: <2005-12-15 20:34:37 raim>
  $Id: sbmlResults.h,v 1.10 2005/12/15 19:54:06 raimc Exp $
*/
/* 
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. The software and
 * documentation provided hereunder is on an "as is" basis, and the
 * authors have no obligations to provide maintenance, support,
 * updates, enhancements or modifications.  In no event shall the
 * authors be liable to any party for direct, indirect, special,
 * incidental or consequential damages, including lost profits, arising
 * out of the use of this software and its documentation, even if the
 * authors have been advised of the possibility of such damage.  See
 * the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * The original code contained here was initially developed by:
 *
 *     Rainer Machne
 *
 * Contributor(s):
 *
 */

#ifndef _SBMLRESULTS_H_
#define _SBMLRESULTS_H_

#include "sbmlsolver/exportdefs.h"


#ifdef __cplusplus
extern "C" {
#endif


  typedef struct timeCourse timeCourse_t ;
  /** A simple structure containing a variable name,
     and its time courses generated by simulation routines.
  */
  struct timeCourse {
    int timepoints;       /**< number of timepoints, including initial
			       conditions */
    char *name;           /**< variable name */
    double *values;       /**< variable time course */
    double **sensitivity; /**< sensitivity time courses */
  } ;

  typedef struct timeCourseArray timeCourseArray_t ;
  /** A simple structure containing num_val time courses */
  struct timeCourseArray {
    int num_val;        /**< number of time courses  */
    timeCourse_t **tc;  /**< time courses */
  } ;


  typedef struct _SBMLResults SBMLResults_t;
  /** A simple structure that contains time courses - represented
      by the timeCourseArray structure - for
      SBML structures, such as species, non-constant compartments
      and parameters and reaction fluxes.
  */
  struct _SBMLResults {

    timeCourse_t *time;                /**< the time points */

    /* concentration and variable parameter and compartment time series */
    timeCourseArray_t *species;       /**< time courses for all species */  
    timeCourseArray_t *compartments;  /**< time courses for all non-constant
					 compartments */
    timeCourseArray_t *parameters;    /**< time courses for all non-constant
					 global parameters */    
    timeCourseArray_t *fluxes;        /**< time courses of reaction fluxes */

    int nsens;
    /** parameters IDs for which sensitivities have been calculated */
    char **param;
    
  } ;



  typedef struct _SBMLResultsMatrix SBMLResultsMatrix_t;
  /** A matrix of _SBMLResults used for batch integration with
      parameter variation via varySettings */
  struct _SBMLResultsMatrix {
    SBMLResults_t ***results;
    int i;
    int j;
  } ;

  SBML_ODESOLVER_API timeCourse_t *SBMLResults_getTime(SBMLResults_t *);  
  SBML_ODESOLVER_API timeCourse_t *SBMLResults_getTimeCourse(SBMLResults_t *, const char *);
  SBML_ODESOLVER_API int SBMLResults_getNout(SBMLResults_t *);
  SBML_ODESOLVER_API int SBMLResults_getNumSens(SBMLResults_t *);
  SBML_ODESOLVER_API const char *SBMLResults_getSensParam(SBMLResults_t *, int);
  SBML_ODESOLVER_API timeCourse_t *Compartment_getTimeCourse(Compartment_t *, SBMLResults_t *);
  SBML_ODESOLVER_API timeCourse_t *Species_getTimeCourse(Species_t *, SBMLResults_t *);
  SBML_ODESOLVER_API timeCourse_t *Parameter_getTimeCourse(Parameter_t *, SBMLResults_t *);
  SBML_ODESOLVER_API const char*TimeCourse_getName(timeCourse_t *);
  SBML_ODESOLVER_API int TimeCourse_getNumValues(timeCourse_t *);
  SBML_ODESOLVER_API double TimeCourse_getValue(timeCourse_t *, int);
  SBML_ODESOLVER_API double TimeCourse_getSensitivity(timeCourse_t *, int, int);
  SBML_ODESOLVER_API void SBMLResults_dump(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResults_dumpSpecies(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResults_dumpCompartments(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResults_dumpParameters(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResults_dumpFluxes(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResults_free(SBMLResults_t *);
  SBML_ODESOLVER_API void SBMLResultsMatrix_free(SBMLResultsMatrix_t *);
  SBML_ODESOLVER_API SBMLResults_t *SBMLResultsMatrix_getResults(SBMLResultsMatrix_t *, int i, int j);

 
#ifdef __cplusplus
}
#endif

/* not part of the API */
SBMLResults_t *SBMLResults_create(Model_t *, int timepoints);
SBMLResultsMatrix_t *SBMLResultsMatrix_allocate(int values, int timepoints);

#endif

/* End of file */
