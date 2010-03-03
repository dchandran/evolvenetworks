#ifndef DESIGNFEEFFORWARDGRN_H
#define DESIGNFEEFFORWARDGRN_H

#include <stdlib.h>
#include <stdio.h>

/*! 
* \brief set the number of inputs in the system
* \param int number of inputs
* \ingroup ffndesign
*/
void setNumInputs(int);

/*! 
* \brief set the number of margins in each input, i.e. number of times the output value changes
* \param int * number of margins for each input
* \ingroup ffndesign
*/
void setGridSizes(int *);

/*! 
* \brief the input->output mapping function
* \param double * array of inputs
* \param int output = 1 or 0 (high or low)
* \ingroup ffndesign
*/
typedef int (*FFN_IO_FUNC)(double*);

/*! 
* \brief set the intput->output mapping function
* \param FFNIOFUNC function
* \ingroup ffndesign
*/
void setTargetFunction(FFN_IO_FUNC);

/*! 
* \brief generates the graphviz + parameters output file
* \param FILE* file for graphviz output
* \param FILE* file for list of parameter values
* \ingroup ffndesign
*/
void generateGraphFile(FILE*, FILE*);

#endif DESIGNFEEFFORWARDGRN_H
