/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	: This class facilitates reading from FASTA/Q files, both interleaved mate-pairs and separated files mate-pairs
 *
 ****************************************************************************/

#ifndef INPUTREADER_H_
#define INPUTREADER_H_

#include <string>
#include "fastAQReader.h"
using namespace std;

class InputReader
{
	string mateFile1;
	string mateFile2;
	FastAQReader *mateReader1;
	FastAQReader *mateReader2;
public:
	read_t read;
	InputReader(string file1, string file2);
	~InputReader();
	bool getNextRead(uint64_t readId);
};

#endif /* INPUTREADER_H_ */
