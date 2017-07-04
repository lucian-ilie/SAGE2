/*****************************************************************************
 * Author		: Mike Molnar and Ehsan Haghshenas
 * Date			: August 27, 2015
 * Description	:
 * Changes		: 
 *
 ****************************************************************************/

#include "main.h"

int main(int argc, char *argv[])
{
	//set default values and initializations
	defaultVals();
	//parse the command line arguments
	parseArgs(argc, argv);
	//check if all the required options are set
	if(!checkRequired())
	{
		cout<< "\n";
		printUsage();
		cout<< "(For more information run ./SAGE2 -h)\n\n";
		exit(0);
	}
	initializeEverything();

	ReadLoader *loaderObj=NULL;
	HashTable *hashObj=NULL;
	EconomyGraph *economyObj=NULL;
	OverlapGraph *graphObj=NULL;
	CopyCountEstimator *estimateObj=NULL;
	MatePair *mateObj=NULL;
	ContigExtender *mergeObj=NULL;
	ScaffoldMaker *scaffObj=NULL;

	//step 1
	if(minStep<=1 && maxStep>=1)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 1\n";
		logStream<< "                                          organizing reads\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();
		loaderObj = new ReadLoader(minOverlap);
		if(listInput!="")
			loaderObj->loadFromList(listInput);
		else
			loaderObj->readDatasetInBytes(fileInput);
		loaderObj->organizeReads();

		if(maxStep==1)
		{
			loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			delete loaderObj;
		}
		else if(saveAll==true)
		{
			loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
		}

	}
	
	if(minStep<=2 && maxStep>=2)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 2\n";
		logStream<< "                                         building hash table\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();
		if(minStep==2)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");
		}

		hashObj = new HashTable(minOverlap, loaderObj);
		hashObj->hashPrefixesAndSuffix(); 
		if(maxStep==2)
		{
			if(minStep==1 && saveAll==false)
				loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			hashObj->saveHashTableInFile(outputDir+prefixName+".hashTable");
			delete hashObj;
			delete loaderObj;
		}
		else if(saveAll==true)
		{
			hashObj->saveHashTableInFile(outputDir+prefixName+".hashTable");
		}
	}

	if(minStep<=3 && maxStep>=3)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 3\n";
		logStream<< "                                      building overlap graph\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();
		if(minStep==3)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");

			hashObj = new HashTable(minOverlap, loaderObj);
			hashObj->loadHashTableFromFile(outputDir+inputPrefix+".hashTable");
		}

		economyObj = new EconomyGraph(minOverlap, hashObj);
		economyObj->buildInitialOverlapGraph();

		economyObj->buildOverlapGraphEconomy();
		delete hashObj; 

		economyObj->sortEconomyGraph(); 

		graphObj = new OverlapGraph(economyObj, loaderObj);
		graphObj->convertGraph();
		delete economyObj;

		if(maxStep==3)
		{
			if(minStep==1 && saveAll==false)
				loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph3");
			delete graphObj;
			delete loaderObj;
		}
		else if(saveAll==true)
		{
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph3");
		}
	}

	if(minStep<=4 && maxStep>=4)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 4\n";
		logStream<< "                                       simplify overlap graph\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();
		if(minStep==4)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");

			graphObj = new OverlapGraph(loaderObj);
			graphObj->loadOverlapGraphFromFile(outputDir+inputPrefix+".graph3");
		}

		int threshold=0, closeValue=10;
		contractCompositePaths(graphObj, loaderObj);
		removeDeadEnds(graphObj, loaderObj, threshold);
		removeBubbles(graphObj, loaderObj, closeValue);
		contractCompositePaths(graphObj, loaderObj);
		
		uint64_t removedEdges = 0;
		
		while(1)
		{
			removedEdges = 0;
			removedEdges += removeDeadEnds(graphObj, loaderObj, threshold);
			removedEdges += removeBubbles(graphObj, loaderObj, closeValue);
			removedEdges += contractCompositePaths(graphObj, loaderObj);

			if(removedEdges==0)	
				break;
			if(closeValue<50)
				closeValue += 10;			
			if(threshold<3)
				threshold++;
		}

		if(maxStep==4)
		{
			if(minStep==1 && saveAll==false)
				loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			delete graphObj;
			delete loaderObj;
		}
	}

	if(minStep<=5 && maxStep>=5)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 5\n";
		logStream<< "                    copy count estimation (min cost flow) & single end assembly\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();

		if(minStep==5)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");

			graphObj = new OverlapGraph(loaderObj);
			graphObj->loadOverlapGraphFromFile(outputDir+inputPrefix+".graph4");
		}

		graphObj->checkIfAllReadsPresent(); 
		estimateObj = new CopyCountEstimator(graphObj, loaderObj);
		genomeSize = estimateObj->genomeSizeEstimation();
		estimateObj->computeMinCostFlow(outputDir);
		delete estimateObj;

		if(maxStep==5)
		{
			if(minStep==1 && saveAll==false)
				loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph5");
			delete graphObj;
			delete loaderObj;
		}
		else if(saveAll==true)
		{
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph5");
		}
	}

	if(minStep<=6 && maxStep>=6)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 6\n";
		logStream<< "                     Map mate pairs & merge single end assemblies to get contigs\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();

		if(minStep==6)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");

			graphObj = new OverlapGraph(loaderObj);
			graphObj->loadOverlapGraphFromFile(outputDir+inputPrefix+".graph5");
		}

		mateObj = new MatePair(graphObj, loaderObj);
		if(listInput!="")
			mateObj->mapMatePairsFromList(listInput);
		else
			mateObj->mapMatePairs(fileInput, "", 1);
		mateObj->meanSdEstimation();

		mergeObj = new ContigExtender(graphObj, mateObj, loaderObj);
		mergeObj->extendContigs(minOverlap);
		delete mergeObj;
		graphObj->printGraph(outputDir+prefixName+"_contig");

		if(maxStep==6)
		{
			if(minStep==1 && saveAll==false)
				loaderObj->saveReadsInFile(outputDir+prefixName+".reads");
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph6");
			delete mateObj;
			delete graphObj;
			delete loaderObj;
		}
		else if(saveAll==true)
		{
			graphObj->saveOverlapGraphInFile(outputDir+prefixName+".graph6");
		}
	}

	if(minStep<=7 && maxStep>=7)
	{
		logStream<< "***********************************************************************************************************\n";
		logStream<< "                                               STEP 7\n";
		logStream<< "                                  merge contigs to get scaffolds\n";
		logStream<< "***********************************************************************************************************\n";
		logStream.flush();
		if(minStep==7)
		{
			loaderObj = new ReadLoader(minOverlap);
			loaderObj->loadReadsFromFile(outputDir+inputPrefix+".reads");

			graphObj = new OverlapGraph(loaderObj);
			graphObj->loadOverlapGraphFromFile(outputDir+inputPrefix+".graph6");
			
			mateObj = new MatePair(graphObj, loaderObj);
			if(listInput!="")
				mateObj->mapMatePairsFromList(listInput);
			else
				mateObj->mapMatePairs(fileInput, "", 1);
			mateObj->meanSdEstimation();
		}

		scaffObj = new ScaffoldMaker(graphObj, mateObj, loaderObj);
		scaffObj->makeScaffolds(minOverlap);

		delete scaffObj;
		graphObj->printGraph(outputDir+prefixName+"_scaffold", true);

		delete mateObj;
		delete graphObj;
		delete loaderObj;
	}
	logStream<< "\n*************************************************** END ***************************************************\n";
	logStream.close();
	return 0;
}

void printUsage()
{
	cout<< "USAGE:\n";
	cout<< "\t./SAGE2 [options] -f <inputFile> -k <minOverlap>\n\n";
}

void printListOfArgs()
{
	cout<< "Full list of options:\n";

	cout<< "\t-f|--fileInput <string>\n";
	cout<< "\t\tInput file in fasta/fastq format (interleaved). Mutually\n\t\texclusive with option -l|--listInput. [Required]\n\n";

	cout<< "\t-l|--listInput <string>\n";
	cout<< "\t\tA list of input files in fasta/fastq format. Mutually exclusive\n\t\twith option -f|--fileInput. (for more information consult with\n\t\tthe README file). [Required]\n\n";

	cout<< "\t-k|--minOverlap <int>\n";
	cout<< "\t\tThe minimum length of an overlap to be considered for graph\n\t\tbuilding. [Required]\n\n";

	cout<< "\t-o|--outputDir <string>\n";
	cout<< "\t\tPath to the output directory. [Default: null]\n\n";

	cout<< "\t-p|--prefix <string>\n";
	cout<< "\t\tPrefix of all output files in the output directory.\n\t\t[Default: \"untitled\"]\n\n";

	cout<< "\t-i|--inputPrefix <string>\n";
	cout<< "\t\tPrefix of all input files in the output directory.\n\t\t[Default: value of prefix]\n\n";

	cout<< "\t-m|--minStep <int>\n";
	cout<< "\t\tThe first step which will be run by the program (for more\n\t\tinformation consult with the README file). [Default: 1]\n\n";

	cout<< "\t-M|--maxStep <int>\n";
	cout<< "\t\tThe last step which will be run by the program (for more\n\t\tinformation consult with the README file). [Default: 7]\n\n";

	cout<< "\t-s|--saveAll\n";
	cout<< "\t\tSaves the result of all intermediate steps.\n\n";

	cout<< "\t-d|--debug\n";
	cout<< "\t\tOutputs and saves more results for debugging.\n\n";

	cout<< "\t-h|--help\n";
	cout<< "\t\tPrints the usage and full list of options.\n\n";
}

void defaultVals()
{
	minStep = 1;
	maxStep = 7;
	saveAll = false;
	debugging = false;
	minOverlap = 0;
	fileInput = "";
	listInput = "";
	outputDir = "";
	prefixName = "untitled";
	inputPrefix = "";
	lowThreshold = 1;
	highThreshold = 5;
}

void initializeEverything()
{
	makeDirectory(outputDir);
	logStream.open((outputDir+prefixName+".log").c_str());
	logStream.imbue(std::locale(""));
	logStream<< std::fixed;
	logStream<< setprecision(2);
	logStream<< "***********************************************************************************************************\n";
	logStream<< "\tEXECUTING PROGRAM: SAGE2 version 1.0\n";
	if(listInput!="")
		logStream<< "\t  INPUT LIST PATH: " << listInput << "\n";
	else
		logStream<< "\t  INPUT FILE PATH: " << fileInput << "\n";
	logStream<< "\t OUTPUT DIRECTORY: " << outputDir << "\n";
	logStream<< "\t    OUTPUT PREFIX: " << prefixName << "\n";
	logStream<< "\t  MINIMUM OVERLAP: " << minOverlap << "\n";
	logStream<< "\t       START STEP: " << (int)minStep << "\n";
	logStream<< "\t         END STEP: " << (int)maxStep << "\n";
	logStream<< "\t   SAVE ALL FILES: " << (saveAll? "TRUE" : "FALSE") << "\n";
	logStream<< "***********************************************************************************************************\n\n";
	logStream.flush();
}

void parseArgs(int argc, char *argv[])
{
	int c;
	int option_index = 0;
	bool fFlag=false, lFlag=false;

	static struct option long_options[] =
	{
		{"help",			no_argument,		0,	'h'},
		{"fileInput",		required_argument,	0,	'f'},
		{"minOverlap",		required_argument,	0,	'k'},
		{"listInput",		required_argument,	0,	'l'},
		{"outputDir",		required_argument,	0,	'o'},
		{"prefix",			required_argument,	0,	'p'},
		{"inputPrefix",		required_argument,	0,	'i'},
		{"minStep",			required_argument,	0,	'm'},
		{"maxStep",			required_argument,	0,	'M'},
		{"saveAll",			no_argument,		0,	's'},
		{"debug",			no_argument,		0,	'd'},
		{0, 0, 0, 0}
	};

	while((c = getopt_long(argc, argv, "hf:l:k:o:p:i:m:M:sd", long_options, &option_index)) != -1)
	{
		switch(c)
		{
		case 0:
			cout<< "option --" << long_options[option_index].name;
			if(optarg)
				cout<< " with argument " << optarg;
			cout<< "\n";
			break;

		case 'h':
			cout<< "\n";
			printUsage();
			printListOfArgs();
			exit(0);
			break;

		case 'f':
			fileInput = optarg;
			fFlag = true;
			break;

		case 'l':
			listInput = optarg;
			lFlag = true;
			break;

		case 'k':
			minOverlap = atoi(optarg);
			break;

		case 'o':
			outputDir = optarg;
			if(outputDir!="")
				outputDir = trimBack(outputDir, "/") + "/";
			break;

		case 'p':
			prefixName = optarg;
			break;

		case 'i':
			inputPrefix = optarg;
			break;

		case 'm':
			minStep = atoi(optarg);
			if(minStep<1)
				minStep = 1;
			break;

		case 'M':
			maxStep = atoi(optarg);
			if(maxStep>7)
				maxStep = 7;
			break;

		case 's':
			saveAll = true;
			break;

		case 'd':
			debugging = true;
			break;

		case '?':
			cout<< "\n";
			exit(0);
			break;

		default:
			cout<< "[ERROR] Wrong command line arguments!\n\n";
			exit(0);
		}
	}

	if (optind < argc)
	{
		cout<< "[WARNING] There are some non-option arguments: ";
		while (optind < argc)
			cout<< argv[optind++] << " ";
		cout<< "\n";
	}

	if(fFlag && lFlag)
	{
		cout<< "[ERROR] Options -f|--fileInput and -l|--listInput are mutually exclusive!\n\n";
		exit(0);
	}
}

bool checkRequired()
{
	bool allSet = true;
	if(fileInput=="" && listInput=="")
	{
		cout<< "[ERROR] One of the options -f|--fileInput or -l|--listInput is required.\n";
		allSet = false;
	}
	if(minOverlap==0)
	{
		cout<< "[ERROR] Option -k|--minOverlap is required.\n";
		allSet = false;
	}
	if(maxStep<minStep)
	{
		cout<< "[ERROR] maxStep should not be smaller than minStep!\n\n";
		allSet = false;
	}
	if(inputPrefix=="")
	{
		inputPrefix = prefixName;
	}
	return allSet;
}
