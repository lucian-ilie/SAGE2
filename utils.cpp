/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 * Changes		: 
 ****************************************************************************/

#include "utils.h"

extern ofstream logStream;

string trim(string str, string pattern)
{
    size_t s = str.find_first_not_of(pattern);
    size_t e = str.find_last_not_of(pattern);

    if((string::npos == s) || (string::npos == e))
        return "";
    else
        return str.substr(s, e-s+1);
}

string trimBack(string str, string pattern)
{
    size_t e = str.find_last_not_of(pattern);

    if(string::npos == e)
        return "";
    else
        return str.substr(0, e+1);
}

/* ==================================================================================
   This function prints an error message
   ================================================================================== */
void printError(string error, string mssg)
{
	logStream<< error << " : " << mssg << "!\n";
	exit(EXIT_FAILURE);
}

/* ============================================================================================
   Make a directory called directoryName.
   ============================================================================================ */
int makeDirectory(string path)
{
	if(path != "")
	{
		int systemRet;
		string command = "";
		command = "mkdir -p " + path;
		systemRet = system(command.c_str());
		return systemRet;
	}
	else
		return 0;
}

int deleteFile(string path)
{
	int systemRet;
	string command = "";
	command = "rm -rf " + path;
	systemRet = system(command.c_str());
	return systemRet;
}

/* ============================================================================================
   This function returns the reverse complement of a given string.
   If the string is GAAAAAACCCCCCG
   It will return   CGGGGGGTTTTTTC
   ============================================================================================ */
string reverseComplement(string &str)
{
	int stringLength1=str.size();
	string retStr(stringLength1, 'N');
	int i, j=0;
	for(i=stringLength1-1; i>=0; i--)
	{
		if(str[i]=='A')
			retStr[j]='T';
		else if(str[i]=='C')
			retStr[j]='G';
		else if(str[i]=='G')
			retStr[j]='C';
		else if(str[i]=='T')
			retStr[j]='A';
		j++;
	}
	return retStr;
}

/* ============================================================================================
   This function computes the bit representation of the reads.
   ============================================================================================ */
uint8_t *charsToBytes(string &read)
{
	int readLength = read.size();
	int readArrayLength=(int)ceil((double)readLength/4.0);
	int i,j=0,shift=2*(4*readArrayLength-readLength);
	uint8_t x, *returnArray;
	if((returnArray=(uint8_t *)malloc((readArrayLength)*sizeof(uint8_t)))==NULL)
		printError(MEM_ALLOC, "returnArray");
	for(i=0;i<(int)readArrayLength;i++)
		returnArray[i]=0;
	for(i=0;i<(int)(readLength);i++)
	{
		j=i/4;
		x=0;
		if(read[i]=='A' || read[i]=='a')	x=0;
		else if(read[i]=='C' || read[i]=='c')	x=1;
		else if(read[i]=='G' || read[i]=='g')	x=2;
		else if(read[i]=='T' || read[i]=='t')	x=3;

		returnArray[j]=(returnArray[j]<<2)|x;
	}
	returnArray[j]=returnArray[j]<<shift;
	return returnArray;
}

/* ============================================================================================
   This function returns the read from the bit representation.
   ============================================================================================ */
string bytesToChars(uint8_t *readInt, int readLength)
{
	string array(readLength, 'N');
	int i,x;
	for(i=0; i<readLength; i++)
	{
		x=(readInt[i>>2]>>(8-2*(i%4+1)))&0X03;
		if(x==0) array[i]='A';
		else if(x==1) array[i]='C';
		else if(x==2) array[i]='G';
		else if(x==3) array[i]='T';
	}
	return array;
}

/* ============================================================================================
   Checks if a read is good or not.
   Good read = No base pairs other than A,C,G and T. Returns 0 if it is a bad read, returns 1 
   if it is a good read.
   ============================================================================================ */
int isGoodRead(string &read, int minOvlp)
{
	int readLength = read.size();
	if(readLength <= minOvlp)// very short read
	{
		return 0;
	}
	int i=0;
	for(i=0;i<readLength;i++)
	{
		if(read[i]=='A' || read[i]=='C' || read[i]=='G' || read[i]=='T')
			continue;
		if(read[i]=='a') read[i]='A';
		else if(read[i]=='c') read[i]='C';
		else if(read[i]=='g') read[i]='G';
		else if(read[i]=='t') read[i]='T';
		else break;
	}
	if(i==readLength)
		return 1;
	else
		return 0;
}

/* ============================================================================================
   Returns 64bit integers that represent a substring of minOverlap.
   ============================================================================================ */
uint64_t* get64Bit2Int(uint8_t *read, uint16_t start, uint16_t length)
{
	uint64_t *number;
	if((number=(uint64_t *)malloc(2*sizeof(uint64_t)))==NULL)
		printError(MEM_ALLOC, "number");
	number[0]=0;number[1]=0;
	if(length<=32)
	{
		number[1]=get64BitInt(read, start, length);
	}
	else
	{
		number[0]=get64BitInt(read, start, length-32);
		number[1]=get64BitInt(read, start+length-32, 32);
	}
	return number;
}

uint64_t get64BitInt(uint8_t *read, uint16_t start, uint16_t length)
{
	uint64_t number=0;
	int byte, fraction1=(start&0X3)<<1, fraction2=((start+length)&0X3)<<1;
	if(start>>2==(start+length)>>2) // Start and End in the same byte.
	{
		number=(read[start>>2]&(0XFF>>fraction1))>>(8-fraction2);
		return number;
	}
	for(byte=start>>2;byte<(start+length)>>2;byte++)
	{
		if(byte==start>>2)
			number=read[byte]&(0XFF>>fraction1);
		else
			number=(number<<8)|read[byte];
	}
	number=(number<<fraction2)|(read[byte]>>(8-fraction2));
	return number;
}

/* ============================================================================================
   This function returns the reverse type of an edge.
   ============================================================================================ */
uint8_t  reverseEdgeType(uint8_t type)
{
	if(type==0)
		return 3;
	else if(type==3)
		return 0;
	return type;
}

/* ============================================================================================
   Compare two strings in bytes.
   ============================================================================================ */
int stringCompareInBytes(uint8_t *read1, uint16_t read1Length, uint8_t *read2, uint16_t read2Length)
{
	int i;
	int read1ArrayLength, read2ArrayLength;
	read1ArrayLength = (int)ceil((double)read1Length/4.0);
	read2ArrayLength = (int)ceil((double)read2Length/4.0);
	for(i=0; i<read1ArrayLength && i<read2ArrayLength; i++)
	{
		if(read1[i]<read2[i])
			return -1;
		if(read1[i]>read2[i])
			return 1;
	}
	if(read1Length<read2Length)
		return -1;
	else if(read1Length>read2Length)
		return 1;
	return 0;
}
