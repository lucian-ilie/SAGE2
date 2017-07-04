/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "fastAQReader.h"

FastAQReader::FastAQReader(string fileName)
{
	if((fp=gzopen(fileName.c_str(), "r"))==NULL)
		printError(OPEN_FILE, fileName);
	seq = kseq_init(fp);
}

bool FastAQReader::getNextRead(read_t &read)
{
	int resVal;
	resVal = kseq_read(seq);
	if(resVal<0)
		return false;
	read.name = seq->name.s;
	seq->name.s[0]='\0';
	read.sequence = seq->seq.s;
	seq->seq.s[0]='\0';
	if(seq->comment.s!=NULL)
	{
		read.comment = seq->comment.s;
		seq->comment.s[0]='\0';
	}
	else
	{
		read.comment = "";
	}
	if(seq->qual.s!=NULL)
	{
		read.quality = seq->qual.s;
		seq->qual.s[0]='\0';
	}
	else
	{
		read.quality = "";
	}
	return true;
}

FastAQReader::~FastAQReader()
{
	gzclose(fp);
	kseq_destroy(seq);
}
