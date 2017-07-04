/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "../inputReader/readLoader.h"
#include "../utils.h"

class HashElement
{
public:
	uint64_t readId:40;
	uint64_t type:2;
};

class HashTable
{
	friend class EconomyGraph;
	uint16_t hashStringLength;
	uint64_t hashThreshold;
	uint64_t longHash;
	uint64_t sizeOfHashTable;
	uint64_t precomputeHash;
	uint64_t hashMiss;
	uint64_t sizeOfHashTablePrevious;
	HashElement **hashTableList;
	uint16_t minOverlap;
	ReadLoader *loaderObj;
public:
	HashTable(uint16_t minOvlp, ReadLoader *loader1);
	~HashTable();
	void hashPrefixesAndSuffix();
	uint64_t hashTableInsert(uint64_t *value, uint64_t readNumber, uint8_t type);
	int64_t hashTableSearch(uint64_t *value, uint64_t &numberOfHashMiss);
	uint64_t getHashValue(uint64_t *value);
	uint64_t findNextPrime(uint64_t value);
	void saveHashTableInFile(string path);
	void loadHashTableFromFile(string path);
	uint64_t findPreviousPrime(uint64_t value);
};

#endif /* HASHTABLE_H_ */
