/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 ****************************************************************************/

#include "readLoader.h"

extern ofstream logStream;
extern uint64_t averageReadLength;
bool Read::operator < (const Read &r2) const
{
	int compRes = stringCompareInBytes(this->readInt, this->length, r2.readInt, r2.length);
	if(compRes < 0)
		return true;
	else
		return false;
}

bool Read::operator == (const Read &r2) const
{
	int compRes = stringCompareInBytes(this->readInt, this->length, r2.readInt, r2.length);
	if(compRes == 0)
		return true;
	else
		return false;
}

ostream& operator<<(ostream &os, const Read &read1)
{
	string read, reverseRead;
	read = bytesToChars(read1.readInt, read1.length);
	reverseRead = bytesToChars(read1.readReverseInt, read1.length);
	os << read1.frequency << "\t" << read1.length << "\t" << read << "\t" << reverseRead;
	return os;
}

istream& operator>>(istream &is, Read &read1)
{
	int freq, len;
	string read, reverseRead;
	is >> freq >> len >> read >> reverseRead;
	read1.frequency = freq;
	read1.length = len;
	read1.readInt = charsToBytes(read);
	read1.readReverseInt = charsToBytes(reverseRead);
	return is;
}

ReadLoader::ReadLoader(uint16_t minOvlp)
{
	numberOfUniqueReads=0;
	totalReads = 0;
	numberOfReads = 0;
	sizeOfList = 0;
	totalBP = 0;
	readsList = NULL;
	minOverlap = minOvlp;
}

ReadLoader::~ReadLoader()
{
	uint64_t i;
	#pragma omp parallel for
	for(i=1; i<=numberOfUniqueReads; i++)
	{
		free(readsList[i].readInt);
		free(readsList[i].readReverseInt);
	}
	free(readsList);
}

void ReadLoader::loadFromList(string listPath)
{
	logStream<< "In function loadFromList().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);

	string line;
	ifstream fin(listPath.c_str());
	if(fin.is_open()==false)
		printError(OPEN_FILE, listPath);
	string var, val, var1, val1;
	uint32_t mateFile=0;
	size_t pos;
	while(getline(fin, line)>0)
	{
		if(line.empty() || line[0]=='#')
			continue;
		pos = line.find('=');
		if(pos==line.npos)
		{
			cout<< "[ERROR] List of input files in a wrong format!\n\n";
			exit(0);
		}
		var = trim(line.substr(0, pos));
		val = trim(line.substr(pos+1));
		if(mateFile%2==0 && var=="f1")
		{
			var1 = var;
			val1 = val;
		}
		else if(mateFile%2==0 && var=="f")
		{
			readDatasetInBytes(val);
			mateFile++;
		}
		else if(mateFile%2==1 && var=="f2")
		{
			readDatasetInBytes(val1, val);
		}
		else
		{
			cout<< "[ERROR] List of input files in a wrong format!\n\n";
			exit(0);
		}
		mateFile++;
	}

	averageReadLength = totalBP/numberOfReads;
	
	logStream<< setw(20) << "Number of datasets: " << mateFile/2 << "\n";
	logStream<< setw(20) << "Total reads: " << totalReads << "\n";
	logStream<< setw(20) << "Good reads: " << numberOfReads << "\n";
	logStream<< setw(20) << "Bad reads: " << totalReads-numberOfReads << "\n";
	logStream<< "\t" << setw(21) << "Average read length: " << averageReadLength << "\n";
	logStream<< setw(20) << "Size of array: " << sizeOfList << "\n";

	logStream<< "Function loadFromList() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

void ReadLoader::readDatasetInBytes(string mateFile1, string mateFile2)
{
	logStream<< "In function readDatasetInBytes().\n";
	time_t seconds_s=time(NULL);
	string read1;
	uint64_t readInFile = 0;
	uint64_t goodInFile = 0;
	uint64_t smallRead = 0;
	
	InputReader *readerObj = new InputReader(mateFile1, mateFile2);
	logStream << "Reading from file: " << mateFile1 << " " << (mateFile2!=""? "& "+mateFile2 : "" ) << "\n";
	logStream.flush();
	while(readerObj->getNextRead(readInFile))
	{
		read1 = readerObj->read.sequence;
		if(read1.size() <= minOverlap)
			smallRead++;
		else
		{
			if(isGoodRead(read1, minOverlap)) 
			{
				insertReadIntoList(read1);
				goodInFile++;
			}
		}
		totalReads++;
		readInFile++; 
	}
	averageReadLength = totalBP/numberOfReads;
	delete readerObj;

	logStream<< "\t" << setw(21) << "Total reads in file: " << readInFile << "\n";
	logStream<< "\t" << setw(21) << "Good reads in file: " << goodInFile << "\n";
	logStream<< "\t" << setw(21) << "Bad reads in file: " << readInFile-goodInFile << "\n";
	logStream<< "\t" << setw(21) << "Reads shorter than minOverlap: " << smallRead << "\n";
	logStream<< "\t" << setw(21) << "Average read length: " << averageReadLength << "\n";
	logStream<< "\t" << setw(21) << "Size of array: " << sizeOfList << "\n";
	

	logStream<< "Function readDatasetInBytes() in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

/* ============================================================================================
   Insert a good read from the file to the list of reads.
   ============================================================================================ */
void ReadLoader::insertReadIntoList(string &read1)
{
	string read1Reverse;
	if(readsList==NULL)
	{
		sizeOfList=1000000;
		if((readsList=(Read *)malloc(sizeOfList*sizeof(Read)))==NULL)
			printError(MEM_ALLOC, "readsList");
	}
	if(numberOfReads>=sizeOfList-10) 
	{
		sizeOfList=sizeOfList+1000000;
		if((readsList=(Read *)realloc(readsList,sizeOfList*sizeof(Read)))==NULL)
			printError(MEM_ALLOC,"reallocating readsList failed");
	}
	read1Reverse=reverseComplement(read1);
	if(read1.compare(read1Reverse)<0) 
	{
		numberOfReads++;
		readsList[numberOfReads].frequency = 0;
		readsList[numberOfReads].length = read1.size();
		totalBP += readsList[numberOfReads].length;
		readsList[numberOfReads].readInt = charsToBytes(read1);
		readsList[numberOfReads].readReverseInt = NULL;
	}
	else
	{
		numberOfReads++;
		readsList[numberOfReads].frequency = 0;
		readsList[numberOfReads].length = read1Reverse.size();
		totalBP += readsList[numberOfReads].length;
		readsList[numberOfReads].readInt = charsToBytes(read1Reverse);
		readsList[numberOfReads].readReverseInt = NULL;
	}
}

void ReadLoader::organizeReads()
{
	logStream<< "\nIn function organizeReads().\n";
	logStream.flush();
	time_t seconds_s=time(NULL);
	uint64_t i, j;
	__gnu_parallel::sort(readsList+1, readsList+numberOfReads+1);
	logStream<< "\tQuicksort reads finished in " << time(NULL)-seconds_s << " sec.\n";
	
	seconds_s=time(NULL);
	for(i=1, j=1; i<=numberOfReads; i++) 
	{
		if(!(readsList[j]==readsList[i]))
		{
			j++;
			if(i!=j)
				swapReads(readsList[j], readsList[i]);
		}
		readsList[j].frequency++;
	}
	numberOfUniqueReads=j;
	
	for(i=numberOfUniqueReads+1; i<=numberOfReads; i++)
		free(readsList[i].readInt);
	
	if(numberOfUniqueReads+1<sizeOfList)
	{
		if((readsList=(Read *)realloc(readsList, (numberOfUniqueReads+1)*sizeof(Read)))==NULL)
			printError(MEM_ALLOC,"reallocating readsList failed");
	}
	logStream<< "\tRemoving duplicate reads finished in " << time(NULL)-seconds_s << " sec.\n";
	
	seconds_s=time(NULL);
	string read1Reverse, read1;
	#pragma omp parallel for default(none) private(i, read1, read1Reverse)
	for(i=1; i<=numberOfUniqueReads; i++)
	{
		read1 = bytesToChars(readsList[i].readInt, readsList[i].length);
		read1Reverse = reverseComplement(read1);
		readsList[i].readReverseInt = charsToBytes(read1Reverse);
	}

	logStream<< "\tNumber of unique reads: " << numberOfUniqueReads << "\n";
	logStream<< "\tComputing reverse complements finished in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

void ReadLoader::swapReads(Read &read1, Read &read2)
{
	Read tmp;
	tmp = read1;
	read1 = read2;
	read2 = tmp;
}

void ReadLoader::saveReadsInFile(string path)
{
	logStream<< "\nIn function saveReadsInFile().\n";
	time_t seconds_s=time(NULL);
	logStream<< "\tSaving in the file : " << path << "\n";
	logStream.flush();
	ofstream fout(path.c_str());
	if(fout.is_open()==false)
		printError(OPEN_FILE, path);
	fout<< numberOfUniqueReads << "\n";
	for(uint64_t i=1; i<=numberOfUniqueReads; i++)
	{
		fout<< readsList[i] << "\n";
	}
	fout.close();
	logStream<< "Function saveReadsInFile in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

void ReadLoader::loadReadsFromFile(string path)
{
	logStream<< "\nIn function loadReadsFromFile().\n";
	time_t seconds_s=time(NULL);
	uint64_t i;
	logStream<< "\tLoading from the file : " << path << endl;
	ifstream fin(path.c_str());
	if(fin.is_open()==false)
		printError(OPEN_FILE, path);
	fin>> numberOfUniqueReads;
	if((readsList=(Read *)malloc((numberOfUniqueReads+1)*sizeof(Read)))==NULL)
		printError(MEM_ALLOC, "readsList");
	for(i=1; i<=numberOfUniqueReads; i++)
		fin >> readsList[i];
	fin.close();
	logStream<< "\tNumber of unique reads: " << numberOfUniqueReads << "\n";
	logStream<< "Function loadReadsFromFile in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}

Read* ReadLoader::getRead(uint64_t readId)
{
	if(readId<1 || readId>numberOfUniqueReads)
		printError(OUT_INDEX, "readsList");
	return (readsList+readId);
}

/* ============================================================================================
   This function finds ID of a read from the sorted list by binary search.
   ============================================================================================ */
int64_t ReadLoader::getIdOfRead(string read)
{
	int64_t lb=1, ub=numberOfUniqueReads, mid;
	int comp_result, flag;
	uint8_t *readInteger;
	string read_r=reverseComplement(read);
	if(read.compare(read_r)<0)
	{
		readInteger=charsToBytes(read);
		flag=1;
	}
	else
	{
		readInteger=charsToBytes(read_r);
		flag=-1;
	}
	while(lb<=ub)
	{
		mid=(ub+lb)>>1;
		comp_result=stringCompareInBytes(readInteger, read.size(), readsList[mid].readInt, readsList[mid].length);
		if(comp_result==0)
		{
			free(readInteger);
			return mid*flag;
		}
		else if(comp_result>0)
			lb=mid+1;
		else
			ub=mid-1;
	}
	free(readInteger);
	logStream<< "\tReturning 0 in getIdOfRead()\n";
	logStream.flush();
	return 0;
}

void ReadLoader::shrinkReadsList()
{
	logStream<< "\nIn function shrinkReadsList().\n";
	time_t seconds_s=time(NULL);
	uint64_t i;
	Read *readsList_tmp;
	if((readsList_tmp=(Read *)malloc((numberOfUniqueReads+1)*sizeof(Read)))==NULL)
		printError(MEM_ALLOC, "readsList_tmp");
	for(i=1; i<=numberOfUniqueReads; i++)
	{
		readsList_tmp[i].frequency = readsList[i].frequency;
		readsList_tmp[i].length = readsList[i].length;
		string read = bytesToChars(readsList[i].readInt, readsList[i].length);
		readsList_tmp[i].readInt = charsToBytes(read);
		free(readsList[i].readInt);
		string reverseRead = bytesToChars(readsList[i].readReverseInt, readsList[i].length);
		readsList_tmp[i].readReverseInt = charsToBytes(reverseRead);
		free(readsList[i].readReverseInt);
	}
	free(readsList);
	readsList = readsList_tmp;
	logStream<< "Function shrinkReadsList in " << time(NULL)-seconds_s << " sec.\n";
	logStream.flush();
}
