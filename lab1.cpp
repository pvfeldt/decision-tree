#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>  
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/imgproc.hpp>
using namespace cv;
using namespace std;

//read dataset;
//in:file name, the number of data lines;
//out:data read and stored in 2D char array DNAReading(column number=361);
char** readFile(string fileName,int lines)
{
    ifstream file(fileName, ios::binary);
    char** DNAReading = new char*[lines];
    for (int i = 0; i < lines; i++)
    {
        DNAReading[i] = new char[361];
    }
    for (int i = 0; i < lines; i++)
    {
        file.getline(DNAReading[i], 361 * sizeof(int));
    }
    return DNAReading;
}

//translate the binary lines to ATCG lines;
//in:2D char array DNAReading,the line number;
//out:2D char array DNATranslated(column number=60);
char** translate(char** dna, int lines)
{
    char** DNATranslated = new char* [lines];
    for (int i = 0; i < lines; i++)
    {
        DNATranslated[i] = new char[61];
    }
    for (int i = 0; i < lines; i++)
    {
        for (int j = 0; j < 60; j++)
        {
            if (dna[i][j * 6] == '0' && dna[i][j * 6 + 2] == '0' && dna[i][j * 6 + 4] == '0')
                DNATranslated[i][j] = 'T';
            else if (dna[i][j * 6] == '0' && dna[i][j * 6 + 2] == '0' && dna[i][j * 6 + 4] == '1')
                DNATranslated[i][j] = 'G';
            else if (dna[i][j * 6] == '0' && dna[i][j * 6 + 2] == '1' && dna[i][j * 6 + 4] == '0')
                DNATranslated[i][j] = 'C';
            else if (dna[i][j * 6] == '1' && dna[i][j * 6 + 2] == '0' && dna[i][j * 6 + 4] == '0')
                DNATranslated[i][j] = 'A';
        }
        DNATranslated[i][60] = dna[i][360];
    }
    return DNATranslated;
}

//optimized log2() function;avoid log2() function outputing NAN;
//in:a number
//out:if number=0, return 0;if number>0,return log2(number);
double log(double num)
{
    if (num > 0)
        return log2(num);
    else return 0;
}

//count the entropy of type ei,ie and n(total);
//in:2D char array DNAReading,the line number;
//out:entropy(total);
double entropyTotal(char** dnat,int lines)
{
    int ei = 0;
    int ie = 0;
    int n = 0;
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][60] == '1')
            ei++;
        else if (dnat[i][60] == '2')
            ie++;
        else if (dnat[i][60] == '3')
            n++;
    }
    int sum = ei + ie + n;
    double eiRate = (double)ei / sum;
    double ieRate = (double)ie / sum;
    double nRate = (double)n / sum;
    double entropy = -eiRate * log(eiRate) - ieRate * log(ieRate) - nRate * log(nRate);
    return entropy;
}

//count the information gain of each column;
//in:2D char array DNAReading,2D char array DNATranslated,line number,column index;
//out:information gain of each column;
double informationGainPerColumn(char** dnat, int lines, int columnindex)
{
    int eiA = 0;
    int ieA = 0;
    int nA = 0;
    int eiT = 0;
    int ieT = 0;
    int nT = 0;
    int eiC = 0;
    int ieC = 0;
    int nC = 0;
    int eiG = 0;
    int ieG = 0;
    int nG = 0;
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][60] == '1')
        {
            if (dnat[i][columnindex] == 'A')
                eiA++;
            else if (dnat[i][columnindex] == 'G')
                eiG++;
            else if (dnat[i][columnindex] == 'T')
                eiT++;
            else if (dnat[i][columnindex] == 'C')
                eiC++;
        }
        else if (dnat[i][60] == '2')
        {
            if (dnat[i][columnindex] == 'A')
                ieA++;
            else if (dnat[i][columnindex] == 'G')
                ieG++;
            else if (dnat[i][columnindex] == 'T')
                ieT++;
            else if (dnat[i][columnindex] == 'C')
                ieC++;
        }
        else if (dnat[i][60] == '3')
        {
            if (dnat[i][columnindex] == 'A')
                nA++;
            else if (dnat[i][columnindex] == 'G')
                nG++;
            else if (dnat[i][columnindex] == 'T')
                nT++;
            else if (dnat[i][columnindex] == 'C')
                nC++;
        }
    }
    int ASum = eiA + ieA + nA;
    int TSum = eiT + ieT + nT;
    int CSum = eiC + ieC + nC;
    int GSum = eiG + ieG + nG;
    int sum = ASum + TSum + CSum + GSum;
    double eiARate = 0;
    double ieARate = 0;
    double nARate = 0; 
    double eiGRate = 0;
    double ieGRate = 0;
    double nGRate = 0;
    double eiTRate = 0;
    double ieTRate = 0;
    double nTRate = 0;
    double eiCRate = 0;
    double ieCRate = 0;
    double nCRate = 0;
    double ARate = 0;
    double GRate = 0;
    double TRate = 0;
    double CRate = 0;
    eiARate = (double)eiA / ASum;
    ieARate = (double)ieA / ASum;
    nARate = (double)nA / ASum;
    if (ASum == 0)
    {
        eiARate = 0;
        ieARate = 0;
        nARate = 0;
    }
    eiGRate = (double)eiG / GSum;
    ieGRate = (double)ieG / GSum;
    nGRate = (double)nG / GSum;
    if (GSum == 0)
    {
        eiGRate = 0;
        ieGRate = 0;
        nGRate = 0;
    }
    eiTRate = (double)eiT / TSum;
    ieTRate = (double)ieT / TSum;
    nTRate = (double)nT / TSum;
    if (TSum == 0)
    {
        eiTRate = 0;
        ieTRate = 0;
        nTRate = 0;
    }
    eiCRate = (double)eiC / CSum;
    ieCRate = (double)ieC / CSum;
    nCRate = (double)nC / CSum;
    if (CSum == 0)
    {
        eiCRate = 0;
        ieCRate = 0;
        nCRate = 0;
    }
    
    ARate = (double)ASum / sum;
    GRate = (double)GSum / sum;
    TRate = (double)TSum / sum;
    CRate = (double)CSum / sum;
    double AEntropy = -eiARate * log(eiARate) - ieARate * log(ieARate) - nARate * log(nARate);
    double GEntropy = -eiGRate * log(eiGRate) - ieGRate * log(ieGRate) - nGRate * log(nGRate);
    double TEntropy = -eiTRate * log(eiTRate) - ieTRate * log(ieTRate) - nTRate * log(nTRate);
    double CEntropy = -eiCRate * log(eiCRate) - ieCRate * log(ieCRate) - nCRate * log(nCRate);
    double informationGain = entropyTotal(dnat, lines) - ARate * AEntropy - GRate * GEntropy - TRate * TEntropy - CRate * CEntropy;
    return informationGain;
}

//find the index of the maximum of information gain of all the columns;
//in:2D char array DNATranslated,line number;
//out:index of the maximum information gain;
int findMaxInformationGain(char** dnat, int lines)
{
    double max = 0;
    int maxColumnIndex = 0;
    for (int i = 0; i < 60; i++)
    {
        double IG = informationGainPerColumn(dnat, lines, i);
        if (informationGainPerColumn(dnat, lines, i) > max)
        {
            max = informationGainPerColumn(dnat, lines, i);
            maxColumnIndex = i;
        }
          
    }
    return maxColumnIndex;
}

//get the updated 2D char array with certain feature from upper level node('A''T''G''C');
//in:2D char array DNATranslated,line number,the upper level node(the max information gain index),certain feature('A''T''G''C');
//out:2D char array DNAUpdated;
char** DNAProcess(char** dnat,int lines,int maxindex,char feature)
{
    vector<int>index;
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][maxindex] == feature)
        {
            index.push_back(i);
        }
    }
    int updatedLines = index.size();
    char** DNAUpdated = new char* [updatedLines];
    for (int i = 0; i < updatedLines; i++)
    {
        DNAUpdated[i] = new char[61];
    }
    for (int i = 0; i < updatedLines; i++)
    {
        for (int j = 0; j < 61; j++)
        {
            DNAUpdated[i][j] = dnat[index[i]][j];
        }
    }
    return DNAUpdated;
}

//the other version to get the updated 2D char array with certain feature from upper level node('A''T''G''C'),for those lines>1,features different and can't separate;
//in:2D array DNATranslated,line number,the upper level node(the max information gain index),certain feature('A''T''G''C');
//out:2D single-line array with certain feature(that appears most);
char** DNAProcess_(char** dnat, int lines, int maxindex, char feature)
{
    int num[3] = { 0,0,0 };
    char numIndex = ' ';
    int maxNum = 0;
    int newIndex = 0;
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][60] == '1')
            num[0]++;
        else if (dnat[i][60] == '2')
            num[1]++;
        else if (dnat[i][60] == '3')
            num[2]++;
    }
    maxNum = num[0];
    numIndex = '1';
    for (int i = 0; i < 3; i++)
    {
        if (num[1] > maxNum)
        {
            maxNum = num[1];
            numIndex = '2';
        }

        else if (num[2] > maxNum)
        {
            maxNum = num[2];
            numIndex = '3';
        }
    }
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][60] == numIndex)
        {
            newIndex = i;
            break;
        }
    }
    char** DNAUpdated_ = new char* [1];
    DNAUpdated_[0] = new char[61];
    for (int i = 0; i < 61; i++)
    {
        DNAUpdated_[0][i] = dnat[newIndex][i];
    }
    return DNAUpdated_;
}

//get the line number of the updated 2D char array DNAUpdated;
//in:2D char array DNATranslated,line number,the upper level node(the max information gain index),certain feature('A''T''G''C');
//out:the updated line number;
int DNAProcessLine(char** dnat, int lines, int maxindex, char feature)
{
    vector<int>index;
    for (int i = 0; i < lines; i++)
    {
        if (dnat[i][maxindex] == feature)
        {
            index.push_back(i);
        }
    }
    int updatedLineNumber = index.size();
    if (updatedLineNumber == lines)
        return 1;
    return updatedLineNumber;
}

//build structure tree;
struct treeNode
{
    int index;
    int layer;
    char tag;
    treeNode* root;
    treeNode* ALeaf;
    treeNode* GLeaf;
    treeNode* TLeaf;
    treeNode* CLeaf;
    treeNode* classification;
    treeNode() :index(0), layer(0), tag(' '), root(NULL), ALeaf(NULL), GLeaf(NULL), TLeaf(NULL), CLeaf(NULL), classification() {}
};



//a global 1D int array for further drawing of the tree;
vector<int> countTrainTreeLayer;

//a int to char function;
//in:number n;
//out:the char form of n;
string intToChar(int n)
{
    char* num = new char[10];
    _itoa_s(n, num, sizeof(num), 10);
    string s = num;
    delete[] num;
    return s;
}

//build tree structure;
//in:2D char array DNATranslated,line number,tree node,layer number;
void tree(char** dnat, int lines, treeNode* t, int layer)
{
    if (lines == 0 || t->tag=='r') return;
    else if (lines == 1)
    {
        cout << "classification:" << dnat[0][60] << " " << "layer:" << t->layer << " " << "tag:" << t->tag << endl;
        countTrainTreeLayer.push_back(layer);
        t->classification = new treeNode();
        t->classification->index = atoi(&dnat[0][60]);
        t->classification->layer = t->layer;
        t->classification->tag = '0';
        return;
    }
    else if (lines > 1)
    {
        int count = 0;
        for (int i = 0; i < lines; i++)
        {
            if (dnat[i][60] == dnat[0][60])
                count++;
        }
        if (count == lines)
        {
            cout << "classification:" << dnat[0][60] << " " << "layer:" << t->layer << " " << "tag:" << t->tag << endl;
            countTrainTreeLayer.push_back(layer);
            t->classification = new treeNode();
            t->classification->index = atoi(&dnat[0][60]);
            t->classification->layer = t->layer;
            t->classification->tag = '0';
            return;
        }
        
      
    }
    int maxindex = findMaxInformationGain(dnat, lines);
    t->index = maxindex;
    t->layer = layer;
    cout << "node index:"<<t->index << " "<<"node layer:" << t->layer << endl;

    char** AUpdated = DNAProcess(dnat, lines, maxindex, 'A');
    int AUpdatedLines = DNAProcessLine(dnat, lines, maxindex, 'A');
    int ALayer = layer + 1;
    t->ALeaf = new treeNode();
    t->ALeaf->index = findMaxInformationGain(AUpdated, AUpdatedLines);
    t->ALeaf->layer = ALayer;
    t->ALeaf->tag = 'a';
    cout << "A->" << endl;
    int AUpdatedLines_ = DNAProcessLine(AUpdated, AUpdatedLines, maxindex, 'A');
    if (AUpdatedLines==AUpdatedLines_)
    {
        int countClass = 0;
        for (int i = 0; i < AUpdatedLines; i++)
        {
            if (AUpdated[i][60] == AUpdated[0][60])
                countClass++;
        }
        if (countClass != AUpdatedLines)
        {
            AUpdated = DNAProcess_(AUpdated, AUpdatedLines, maxindex, 'A');
            AUpdatedLines = 1;
        }
    }
    if (AUpdatedLines == 0)
    {
        t->ALeaf->tag = 'r';
    }
    tree(AUpdated, AUpdatedLines, t->ALeaf, ALayer);

    char** GUpdated = DNAProcess(dnat, lines, maxindex, 'G');
    int GUpdatedLines = DNAProcessLine(dnat, lines, maxindex, 'G');
    int GLayer = layer + 1;
    t->GLeaf = new treeNode();
    t->GLeaf->index = findMaxInformationGain(GUpdated, GUpdatedLines);
    t->GLeaf->layer = GLayer;
    t->GLeaf->tag = 'g';
    cout << "G->" << endl; 
    int GUpdatedLines_ = DNAProcessLine(GUpdated, GUpdatedLines, maxindex, 'G');
    if (GUpdatedLines == GUpdatedLines_)
    {
        int countClass = 0;
        for (int i = 0; i < GUpdatedLines; i++)
        {
            if (GUpdated[i][60] == GUpdated[0][60])
                countClass++;
        }
        if (countClass != GUpdatedLines)
        {
            GUpdated = DNAProcess_(GUpdated, GUpdatedLines, maxindex, 'G');
            GUpdatedLines = 1;
        }
    }
    if (GUpdatedLines == 0)
    {
        t->GLeaf->tag='r';
    }
    tree(GUpdated, GUpdatedLines, t->GLeaf, GLayer);

    char** TUpdated = DNAProcess(dnat, lines, maxindex, 'T');
    int TUpdatedLines = DNAProcessLine(dnat, lines, maxindex, 'T');
    int TLayer = layer + 1;
    t->TLeaf = new treeNode();
    t->TLeaf->index = findMaxInformationGain(TUpdated, TUpdatedLines);
    t->TLeaf->layer = TLayer;
    t->TLeaf->tag = 't';
    cout << "T->" << endl;
    int TUpdatedLines_ = DNAProcessLine(TUpdated, TUpdatedLines, maxindex, 'T');
    if (TUpdatedLines == TUpdatedLines_)
    {
        int countClass = 0;
        for (int i = 0; i < TUpdatedLines; i++)
        {
            if (TUpdated[i][60] == TUpdated[0][60])
                countClass++;
        }
        if (countClass != TUpdatedLines)
        {
            TUpdated = DNAProcess_(TUpdated, TUpdatedLines, maxindex, 'T');
            TUpdatedLines = 1;
        }
    }
    if (TUpdatedLines == 0)
    {
        t->TLeaf->tag = 'r';
    }
    tree(TUpdated, TUpdatedLines, t->TLeaf, TLayer);

    char** CUpdated = DNAProcess(dnat, lines, maxindex, 'C');
    int CUpdatedLines = DNAProcessLine(dnat, lines, maxindex, 'C');
    int CLayer = layer + 1;
    t->CLeaf = new treeNode();
    t->CLeaf->index = findMaxInformationGain(CUpdated, CUpdatedLines);
    t->CLeaf->layer = CLayer;
    t->CLeaf->tag = 'c';
    cout << "C->" << endl;
    int CUpdatedLines_ = DNAProcessLine(CUpdated, CUpdatedLines, maxindex, 'C');
    if (CUpdatedLines == CUpdatedLines_)
    {
        int countClass = 0;
        for (int i = 0; i < CUpdatedLines; i++)
        {
            if (CUpdated[i][60] == CUpdated[0][60])
                countClass++;
        }
        if (countClass != CUpdatedLines)
        {
            CUpdated = DNAProcess_(CUpdated, CUpdatedLines, maxindex, 'C');
            CUpdatedLines = 1;
        }
    }
    if (CUpdatedLines == 0)
    {
        t->CLeaf->tag = 'r';
    }
    tree(CUpdated, CUpdatedLines, t->CLeaf, CLayer);
}

//whole process of creating a decision tree;
//in:file name,lines;
//out:tree node;
treeNode* createDecisionTree(string filename,int lines)
{
    char** DNAReading = readFile(filename,lines);
    char** DNATranslated = translate(DNAReading, lines);
    int maxIndex = findMaxInformationGain(DNATranslated, lines);
    int layerDefault = 0;
    treeNode* t = new treeNode();
    t->index = maxIndex;
    t->layer = layerDefault;
    tree(DNATranslated, lines, t, layerDefault);
    return t;
}

//count the maximum layer number of the tree;
//in:1D array countTrainTreeLayer or countTestTreeLayer;
//out:the maximum layer number;
int maxLayer(vector<int> cl)
{
    int max = 0;
    int clsize = cl.size();
    for (int i = 0; i < clsize; i++)
    {
        if (cl[i] > max)
            max = cl[i];
    }
    return max;
}

//load in image "untitled.jpg";
Mat image = imread("untitled.jpg");

//draw a tree;
//in:tree node,starting point x,starting point y,string output(initial:null),1D array countTrainTreeLayer;
void drawTree(treeNode* t, double pointX, double pointY, string output, vector<int> ctl)
{
    double radius =0.15;//r;
    int maxTrainLayerCount = maxLayer(ctl);
    double previousX = pointX;
    double previousY = pointY;
    if (t == NULL) return;
    else
    {
        if (t->tag == 'a')
        {
            double L = pow(4, maxTrainLayerCount) * 2 * radius - pow(4, maxTrainLayerCount - t->layer) * 2 * radius;
            double L0 = L * 3 / (pow(4, t->layer) - 1);
            if (t->ALeaf == NULL && t->GLeaf == NULL && t->TLeaf == NULL && t->CLeaf == NULL)
            {
                L = pow(4, maxTrainLayerCount) * 2 * radius;
                L0 = L / pow(4, t->layer - 1) - 2 * radius;
            }
            pointX += -L0 / 2;
            pointY += 1000;
        }
        else if (t->tag == 'g')
        {
            double L = pow(4, maxTrainLayerCount) * 2 * radius - pow(4, maxTrainLayerCount - t->layer) * 2 * radius;
            double L0 = L * 3 / (pow(4, t->layer) - 1);
            if (t->ALeaf == NULL && t->GLeaf == NULL && t->TLeaf == NULL && t->CLeaf == NULL)
            {
                L = pow(4, maxTrainLayerCount) * 2 * radius;
                L0 = L / pow(4, t->layer - 1) - 2 * radius;
            }
            pointX += -L0 / 6;
            pointY += 1000;
        }
        else if (t->tag == 't')
        {
            double L = pow(4, maxTrainLayerCount) * 2 * radius - pow(4, maxTrainLayerCount - t->layer) * 2 * radius;
            double L0 = L * 3 / (pow(4, t->layer) - 1);
            if (t->ALeaf == NULL && t->GLeaf == NULL && t->TLeaf == NULL && t->CLeaf == NULL)
            {
                L = pow(4, maxTrainLayerCount) * 2 * radius;
                L0 = L / pow(4, t->layer - 1) - 2 * radius;
            }
            pointX += L0 / 6;
            pointY += 1000;
        }
        else if (t->tag == 'c')
        {
            double L = pow(4, maxTrainLayerCount) * 2 * radius - pow(4, maxTrainLayerCount - t->layer) * 2 * radius;
            double L0 = L * 3 / (pow(4, t->layer) - 1);
            if (t->ALeaf == NULL && t->GLeaf == NULL && t->TLeaf == NULL && t->CLeaf == NULL)
            {
                L = pow(4, maxTrainLayerCount) * 2 * radius;
                L0 = L / pow(4, t->layer - 1) - 2 * radius;
            }
            pointX += L0 / 2;
            pointY += 1000;
        }
        else if (t->tag == ' ')
        {
            pointX = 10000;
            pointY = 1000;
        }
        output = intToChar(t->index);
        if (t->classification != NULL)
            output = intToChar(t->classification->index);
        line(image, Point(previousX, previousY), Point(pointX, pointY), Scalar(0, 0, 0), 20);
        circle(image, Point(previousX, previousY), radius, Scalar(0, 0, 0), 20);
        circle(image, Point(previousX, previousY), radius - 0.1, Scalar(255, 255, 255), -1);
        circle(image, Point(pointX, pointY), radius, Scalar(0, 0, 0), 20);
        circle(image, Point(pointX, pointY), radius - 0.1, Scalar(255, 255, 255), -1);
        drawTree(t->ALeaf, pointX, pointY, output, ctl);
        drawTree(t->GLeaf, pointX, pointY, output, ctl);
        drawTree(t->TLeaf, pointX, pointY, output, ctl);
        drawTree(t->CLeaf, pointX, pointY, output, ctl);
        putText(image, &output[0], Point(pointX - 100, pointY + 100), 0, 10, Scalar(0, 0, 0), 20);
    }
}

//draw a tree totally;
//in:tree node,1D array countTrainTreeLayer;
void drawTreeTotal(treeNode* t, vector<int> ctl)
{
    namedWindow("decision tree",CV_WINDOW_NORMAL);
    double circlePointX = 10000;
    double circlePointY = 1000;
    drawTree(t,circlePointX,circlePointY,"",ctl);
    imshow("decision tree", image);
}

//traverse the trained tree;
//in:a single line from 2D array DNATranslated(test set),tree node;
//out:if correct,return 1;if incorrect,return 0;
int traverseTree(char* dnatsingleline, treeNode* t)
{
    if (t != NULL)
    {
      if (t->classification != NULL)
        {
            if (t->classification->index == atoi(&dnatsingleline[60]))
                return 1;
            else
                return 0;
        }
        if (dnatsingleline[t->index] == 'A')
        {
            traverseTree(dnatsingleline, t->ALeaf);
        }
        else if (dnatsingleline[t->index] == 'G')
        {
            traverseTree(dnatsingleline, t->GLeaf);
        }
        else if (dnatsingleline[t->index] == 'T')
        {
            traverseTree(dnatsingleline, t->TLeaf);
        }
        else if (dnatsingleline[t->index] == 'C')
        {
            traverseTree(dnatsingleline, t->CLeaf);
        }
    }
}

//count accuracy of the tree model;
//in:file name,the line number,tree number;
//out:the result of accuracy;
double accuracy(string filename, int lines,treeNode* t)
{
    char** DNAReadingTest = readFile(filename, lines);
    char** DNATranslatedTest = translate(DNAReadingTest, lines);
    vector<int> traverseResult;
    int traverse = 0;
    for (int i = 0; i < lines; i++)
    {
        traverse = traverseTree(DNATranslatedTest[i], t);
        traverseResult.push_back(traverse);
    }
    int count = 0;
    for (int i = 0; i < lines; i++)
    {
        if (traverseResult[i] == 1)
            count++;
    }
    double accuracyResult = (double)count / lines;
    return accuracyResult;
}
int main()
{
    cout << "training set:" << endl;
    treeNode* treeTrain=createDecisionTree("dna(1).data", 2000);
    drawTreeTotal(treeTrain, countTrainTreeLayer);
    cout << "test set accuracy:" << endl;
    cout << accuracy("dna(1).test", 1186, treeTrain) * 100 << "%" << endl;
    delete[] treeTrain;
    waitKey(0);
    system("pause");
    return 0;
}

