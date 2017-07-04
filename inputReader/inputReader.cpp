/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	: This class facilitates reading from FASTA/Q files, both interleaved mate-pairs and separated files mate-pairs
 *
 ****************************************************************************/

#include "inputReader.h"

InputReader::InputReader(string file1, string file2)
{
	mateFile1 = file1;
	mateFile2 = file2;
	mateReader1 = new FastAQReader(mateFile1);
	if(mateFile2!="")
		mateReader2 = new FastAQReader(mateFile2);
}

InputReader::~InputReader()
{
	delete mateReader1;
	if(mateFile2!="")
		delete mateReader2;
}

bool InputReader::getNextRead(uint64_t readId)
{
	if(readId%2==0)
	{
		if(mateReader1->getNextRead(read))
			return true;
		return false;
	}
	else
	{
		if(mateFile2!="")
		{
			if(mateReader2->getNextRead(read))
				return true;
			return false;
		}
		else
		{
			if(mateReader1->getNextRead(read))
				return true;
			return false;
		}
	}
}
