/*****************************************************************************
 * Author		:
 * Date			:
 * Description	:
 ****************************************************************************/

#ifndef UTILS_H_
#define UTILS_H_

#include <stdint.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

using namespace std;
#define OPEN_FILE "Could not open file"
#define MEM_ALLOC "Malloc failed"
#define OUT_INDEX "Out of Index"
#define OTHERS "ERROR"

string trim(string str, string pattern=" \n\r\t");
string trimBack(string str, string pattern=" \n\r\t");
void printError(string error, string mssg);
int makeDirectory(string path);
int deleteFile(string path);
string reverseComplement(string &str);
uint8_t *charsToBytes(string &read);
string bytesToChars(uint8_t *readInt, int readLength);
int isGoodRead(string &read, int minOvlp);
uint64_t* get64Bit2Int(uint8_t *read, uint16_t start, uint16_t length);
uint64_t get64BitInt(uint8_t *read, uint16_t start, uint16_t length);
uint8_t  reverseEdgeType(uint8_t type);
int stringCompareInBytes(uint8_t *read1, uint16_t read1Length, uint8_t *read2, uint16_t read2Length);

#endif /* UTILS_H_ */
