/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef MAIN_H_
#define MAIN_H_

#include <iostream>
#include <stdlib.h>
#include <getopt.h>
#include <omp.h>
#include <locale>

#include "inputReader/readLoader.h"
#include "economyGraph/hashTable.h"
#include "economyGraph/economyGraph.h"
#include "overlapGraph/overlapGraph.h"
#include "overlapGraph/simplification.h"
#include "overlapGraph/copycountEstimator.h"
#include "matePair/matePair.h"
#include "matePair/mergeContigs.h"
#include "matePair/scaffolding.h"
using namespace std;

//command-line option variables
uint8_t minStep, maxStep;
uint16_t minOverlap;
string fileInput, listInput, outputDir, prefixName, inputPrefix;
bool saveAll, debugging;
int lowThreshold, highThreshold;
//global variables
time_t time_start;
ofstream logStream;
uint64_t genomeSize=0, averageReadLength=0;

void printUsage();
void printListOfArgs();
void defaultVals();
void initializeEverything();
void closeEverything();
void parseArgs(int argc, char *argv[]);
bool checkRequired();

#endif /* MAIN_H_ */
