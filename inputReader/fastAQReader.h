/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#ifndef FASTAQREADER_H_
#define FASTAQREADER_H_

#include <iostream>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "../utils.h"

using namespace std;

typedef struct
{
	string name;
	string comment;
	string sequence;
	string quality;
} read_t;

class FastAQReader
{
	KSEQ_INIT(gzFile, gzread)
	gzFile fp;
	kseq_t *seq;
public:
	FastAQReader(string fileName);
	~FastAQReader();
	bool getNextRead(read_t &read);
};

#endif /* FASTAQREADER_H_ */
