/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef READLOADER_H_
#define READLOADER_H_

#include <string>
#include <algorithm>
#include <parallel/algorithm>
#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <zlib.h>
#include "../utils.h"
#include "inputReader.h"
using namespace std;

class Read
{
public:
	uint8_t *readInt;
	uint8_t *readReverseInt;
	uint16_t frequency;
	uint16_t length;
	bool operator < (const Read &r2) const;
	bool operator == (const Read &r2) const;
};

class ReadLoader
{
public:
	uint64_t numberOfUniqueReads;
	uint64_t numberOfReads;
	uint16_t minOverlap;
	uint64_t totalBP;
private:
	uint64_t totalReads;
	uint64_t sizeOfList;
	Read *readsList;
public:
	ReadLoader(uint16_t minOvlp);
	~ReadLoader();
	void loadFromList(string listPath);
	void readDatasetInBytes(string mateFile1, string mateFile2="");
	void insertReadIntoList(string &read1);
	void organizeReads();
	void swapReads(Read &read1, Read &read2);
	void saveReadsInFile(string path);
	void loadReadsFromFile(string path);
	Read* getRead(uint64_t readId);
	int64_t getIdOfRead(string read);
	void shrinkReadsList();
};

#endif /* READLOADER_H_ */
