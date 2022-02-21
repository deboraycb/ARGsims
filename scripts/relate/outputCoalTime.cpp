// Programmer: April Wei
// Use c++11
// Start: Jun 17, 2020
// Last edit: Jun 17, 2020
/*
Usage: This program aims at outputing the quantile of the msPrime simulated coalescence time compare to the Relate MCMC samples

Compile with: 
	g++ -Wall -std=c++11 -g outputCoalTime.cpp -o outputCoalTime

Run with (example): ./outputCoalTime testCoalTimeQuantile.txt coalTimeMSPrime/coalTimeInd0Ind1Sample4Replicate1.txt \
relateFiles/relateTreeGenomeSample4Replicate1/genome1000MCMCSample4Replicate1.anc \
relateFiles/relateTreeGenomeSample4Replicate1/genome1000MCMCSample4Replicate1.mut 1000 0 1

Check out the output:
less -N testMCMCCounts.txt |awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' | less
*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <algorithm>

using namespace std;


class NodeClass
{
    private:
        int parentNode;
        double *branchLength;
        double *coalTimeAtThisNode;
        //int* listOfPopAtThisNode;
        int* listOfChildren;
        int numOfChildenFound;
        int childIndex1;
        int childIndex2;
    public:
        void readNode(string line, int numMCMC);
        int getParentNode();
        int* getListOfChildren();
        int getNumOfChildenFound();
        void setCoalTime(double* inBranchLength, double* inTime, int numMCMC);
        double* getCoalTime();
        double* getBranchLength();
        void assignSelfPop(int inPop);
        void addChildren(int index, int inChild);
        void addnumOfChildenFound();
        ~NodeClass();
};



int calcCoalescnceTime(NodeClass* nodeList, int numHaplotypes, int numMCMC, int node1Index, int node2Index)
{
    int numNodes = 2*numHaplotypes - 1;
    int indexParent;
    int iterNodes = 0;
    while(iterNodes < numNodes -1)// the last node is the root and it does not have a parent node;
    {
        indexParent = nodeList[iterNodes].getParentNode();
        if (nodeList[indexParent].getNumOfChildenFound() == 0)
        {
            nodeList[indexParent].addnumOfChildenFound();
            nodeList[indexParent].addChildren(0,iterNodes);
            nodeList[indexParent].setCoalTime(nodeList[iterNodes].getBranchLength(), nodeList[iterNodes].getCoalTime(), numMCMC);
        }
        else
        {
            nodeList[indexParent].addChildren(1,iterNodes);
            nodeList[indexParent].addnumOfChildenFound();
        }
        iterNodes++;
    }
    int ancestorOfNodes[numNodes] = {0};
    while (node1Index < numNodes - 1)
    {
        node1Index = nodeList[node1Index].getParentNode();
        ancestorOfNodes[node1Index] = 1;
    }
    while (true)
    {
        indexParent = nodeList[node2Index].getParentNode();
        if (ancestorOfNodes[indexParent] == 1)
        {
            return indexParent;
        }
        else
        {
            node2Index = indexParent;
        }
    }
}

void outputCoalTime(NodeClass* nodeList, int indexCommonAncestor, int numMCMC, ofstream& outFile, int twoNe)
{
    double* coalTime = nodeList[indexCommonAncestor].getCoalTime();
    for (int iter = 0; iter < numMCMC; iter++)
    {
        outFile << '\t' << coalTime[iter]/twoNe;
    }
    outFile << '\n';
}

bool readMutFile(ifstream& mutFile, int* coordinateInMut)
{
    string line;
    int iter = 0;
    getline(mutFile,line);
    string tempStr;
    int position;
    while (getline(mutFile, line))
    {
        stringstream ssin(line);
        getline(ssin, tempStr, ';');
        iter = stoi(tempStr);
        getline(ssin, tempStr, ';');
        position = stoi(tempStr);
        coordinateInMut[iter] = position;
    }
    return true;
}

bool readMSPrimeFile(ifstream& msPrimeFile, double* msPrimeCoorStart, double* msPrimeCoorEnd, double* msPrimeCoalTime)
{
    string line;
    string tempStr;
    int iter = 0;
    while (getline(msPrimeFile, line))
    {
        stringstream ssin(line);
        getline(ssin, tempStr, '\t');
        msPrimeCoorStart[iter] = stod(tempStr);
        getline(ssin, tempStr, '\t');
        msPrimeCoorEnd[iter] = stod(tempStr);
        getline(ssin, tempStr);
        msPrimeCoalTime[iter] = stod(tempStr);
        iter++;
    }
    return true;
}

bool  readTreeFileHeader(ifstream& inFile, int& numHaplotypes, int& numTrees, int& numMCMC)
{
    string line;
    string tempStr;
    for (int iter = 0; iter < 3; iter++)
    {
        getline(inFile, line);
        stringstream ssin(line);
        ssin >> tempStr;
        ssin >> tempStr;
        if (iter == 0)
        {
            numHaplotypes = stoi(tempStr);
        }
        else if (iter == 1)
        {
            numTrees = stoi(tempStr);
        }
        else
        {
            numMCMC = stoi(tempStr);
        }
    }
    cout << "numHaplotypes = " << numHaplotypes << ", numTrees = " << numTrees 
         << ", numMCMC = " << numMCMC << endl;
    return true;
}


void printMyInfo()
{
    cout << "\n\n";
    cout << "----------------------------------------------------------------\n";
    cout << "---------You are using the program written by April Wei---------\n";
    cout << "For questions and comments, please email: aprilwei@berkeley.edu.\n";
    cout << "----------------------------------------------------------------\n\n";
}

void printExampleInput()
{
    cout << "The input arguments are wrong!" << endl;
    cout << "Here is the example input: ./outputCoalTime <OutputFileName>  <InputMSPrimeFileName> <InputRelateMCMCFileName> <InputMutFileName> EffPopSize(int) node1Index node2Index" << endl;
}


int main(int argc, char *argv [])
{
    ifstream inFile;
    string inFileName;
    ofstream outFile;
    string outFileName;
    string msPrimeFileName;
    ifstream msPrimeFile;
    string mutFileName;
    ifstream mutFile;
    if (argc != 8)
    {
        printExampleInput();
    }
    else
    {
        printMyInfo();
        //----------collect commandline input---------
        outFileName = argv[1];
        outFile.open(outFileName.c_str());
        msPrimeFileName = argv[2];//the msPrimeFile that contains the true coalescence time
        msPrimeFile.open(msPrimeFileName.c_str());
        inFileName = argv[3];
        inFile.open(inFileName.c_str());
        mutFileName = argv[4];
        mutFile.open(mutFileName.c_str());//should replace this by .mut file
        int Ne =  stoi(argv[5]);
        int node1Index =  stoi(argv[6]);
        int node2Index =  stoi(argv[7]);
        cout << "\nparameters read: \n"
             << "InputRelateMCMCFileName: " << inFileName << '\n'
             << "InputMSPrimeFileName: " << msPrimeFileName << '\n'
             << "outFileName: " << outFileName << '\n'
             << "InputMutFileName: " << mutFileName << '\n'
             << "Ne (diploid)= " << Ne << '\n'
             << "node1Index = " << node1Index << '\n' 
             << "node2Index = " << node2Index << '\n';
 //------------------------------------------
        int numHaplotypes;
        int numTrees;
        int numMCMC;
        if (readTreeFileHeader(inFile, numHaplotypes, numTrees, numMCMC))
        {
            cout << "Tree file header read\n";
            int* coordinateInMut;
            int numLines = count(istreambuf_iterator<char>(mutFile),istreambuf_iterator<char>(), '\n');
            coordinateInMut = new int [numLines - 1];
            mutFile.close();
            mutFile.open(mutFileName.c_str());
            double * msPrimeCoalTimeArray;
            double * msPrimeCoorEndArray;
            double * msPrimeCoorStartArray;
            int numMSPrimeLines = count(istreambuf_iterator<char>(msPrimeFile),istreambuf_iterator<char>(), '\n');
            msPrimeCoorEndArray = new double [numMSPrimeLines];
            msPrimeCoorStartArray = new double [numMSPrimeLines];
            msPrimeCoalTimeArray = new double [numMSPrimeLines];
            msPrimeFile.close();
            msPrimeFile.open(msPrimeFileName.c_str());
            int numNodes = 2*numHaplotypes - 1;
            if (readMutFile(mutFile, coordinateInMut) && readMSPrimeFile(msPrimeFile, msPrimeCoorStartArray, msPrimeCoorEndArray, msPrimeCoalTimeArray))
            {
                cout << "mut file is read\n msPrimeFile is read \n";
                string line;
                NodeClass* nodeList;
        //read in the first tree
                getline(inFile, line);
                stringstream inTreeStream(line);
                string tempStr;
                inTreeStream >> tempStr;
                nodeList = new NodeClass [numNodes];
                for (int iter = 0; iter < numNodes; iter++)//create nodeList;
                {
                    getline(inTreeStream, tempStr, ')');
                    nodeList[iter].readNode(tempStr, numMCMC);
                }
        //create the parameters needed.
                double coorRelateTreeEnd;
                int iterMSPrimeTree = 0;
                double outputRegionStart = 0;
                int indexOfCommonAncestor;
                while (getline(inFile, line))
                {
            //create the stringstream of the next line
                    stringstream inTreeStream(line);
                    string tempStr;
                    inTreeStream >> tempStr;
                    coorRelateTreeEnd = (double) coordinateInMut[stoi(tempStr)];
            //process the coalescence time and output this tree;
                    indexOfCommonAncestor = calcCoalescnceTime(nodeList, numHaplotypes, numMCMC, node1Index, node2Index);
                    if (msPrimeCoorEndArray[iterMSPrimeTree] > coorRelateTreeEnd) //move to the next relate tree in the outer while loop
                    {
                        outFile << coorRelateTreeEnd - outputRegionStart << '\t' << msPrimeCoalTimeArray[iterMSPrimeTree];
                        outputCoalTime(nodeList, indexOfCommonAncestor, numMCMC, outFile, 2 * Ne);
                        outputRegionStart = coorRelateTreeEnd;
                    }
                    else
                    {
                        while (msPrimeCoorEndArray[iterMSPrimeTree] < coorRelateTreeEnd)
                        {
                            outFile << msPrimeCoorEndArray[iterMSPrimeTree] - outputRegionStart  << '\t' << msPrimeCoalTimeArray[iterMSPrimeTree];
                            outputCoalTime(nodeList, indexOfCommonAncestor, numMCMC, outFile, 2 * Ne);
                            outputRegionStart = msPrimeCoorEndArray[iterMSPrimeTree];
                            iterMSPrimeTree++;
                        }
                    //finish outputing the remaining region of the msPrimeTree overlapping with thisRelateTree 
                        outFile << coorRelateTreeEnd - outputRegionStart  << '\t' << msPrimeCoalTimeArray[iterMSPrimeTree];
                        outputCoalTime(nodeList, indexOfCommonAncestor, numMCMC, outFile, 2 * Ne);
                        outputRegionStart = coorRelateTreeEnd;
                    }
            //delete this nodeList in the end to free up the memory
                    delete [] nodeList;
            //read in the nodeList of this already read line.
                    nodeList = new NodeClass [numNodes];
                    for (int iter = 0; iter < numNodes; iter++)//create nodeList;
                    {
                        getline(inTreeStream, tempStr, ')');
                        nodeList[iter].readNode(tempStr, numMCMC);
                    }
                }
                indexOfCommonAncestor = calcCoalescnceTime(nodeList, numHaplotypes, numMCMC, node1Index, node2Index);
                coorRelateTreeEnd = msPrimeCoorEndArray[numMSPrimeLines - 1];
                while (msPrimeCoorEndArray[iterMSPrimeTree] < coorRelateTreeEnd)
                {
                    outFile << msPrimeCoorEndArray[iterMSPrimeTree] - outputRegionStart  << '\t' << msPrimeCoalTimeArray[iterMSPrimeTree];
                    outputCoalTime(nodeList, indexOfCommonAncestor, numMCMC, outFile, 2 * Ne);
                    outputRegionStart = msPrimeCoorEndArray[iterMSPrimeTree];
                    iterMSPrimeTree++;
                }
                outFile << coorRelateTreeEnd - outputRegionStart  << '\t' << msPrimeCoalTimeArray[iterMSPrimeTree];
                outputCoalTime(nodeList, indexOfCommonAncestor, numMCMC, outFile, 2 * Ne);
            }
            else
            {
                cout << "cannot read mutFile\n";
                return (1);
            }
        }
        else
        {
            cout << "cannot read tree file heade\n";
            return (1);
        }   
    }
    return (0);
}


void NodeClass::readNode(string line, int numMCMC)
{
    stringstream inNodeStream(line);
    string tempStr;
    //read parentNode and branchLength
    getline(inNodeStream, tempStr, ':');
    parentNode = stoi(tempStr);
    getline(inNodeStream, tempStr, '(');
    branchLength = new double [numMCMC];
    coalTimeAtThisNode = new double [numMCMC];
    for (int iter = 0; iter < numMCMC; iter++)
    {
        inNodeStream >> tempStr;
        branchLength[iter] = stod(tempStr);
        coalTimeAtThisNode[iter] = 0;
    }
    //initiate the parameter yet to be determined
    numOfChildenFound = 0;
    listOfChildren = new int [2];
    listOfChildren[0] = 0;
    listOfChildren[1] = 0;
}


NodeClass::~NodeClass() 
{ 
    delete [] branchLength;
    delete [] coalTimeAtThisNode;
    delete [] listOfChildren;
} 
void NodeClass::addnumOfChildenFound()
{
    numOfChildenFound++;
}

int NodeClass::getNumOfChildenFound()
{
    return numOfChildenFound;
}

int* NodeClass::getListOfChildren()
{
    return listOfChildren;
}

int NodeClass::getParentNode()
{
    return parentNode;
}

void NodeClass::setCoalTime(double* inBranchLength, double* inTime, int numMCMC)
{
    for (int iter = 0; iter < numMCMC; iter++)
    {
        coalTimeAtThisNode[iter] = inBranchLength[iter] + inTime[iter];
    }
}

double* NodeClass::getCoalTime()
{
    return coalTimeAtThisNode;
}


double* NodeClass::getBranchLength()
{
    return branchLength;
}

void NodeClass::addChildren(int index, int inChild)
{
    listOfChildren[index] = inChild;
}


