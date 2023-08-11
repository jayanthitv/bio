
/*
 *  This code evolves error correcting codes over the Hamming metric
 *  using the Conway crossover operator.  In selects those data from
 *  and input file of example DNA sequences.  This version has been
 *  modified to permit command line arguments for the input data file,
 *  the minimum distance, the random material rate, and the population
 *  size.
 *
 *  By collaborator request, this copy is also hacked to create codes
 *  that share no elements at all.
 *
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

using namespace std;

#include "stat.h"

//random number seed
#define RNS 73214899

//DNA sequence length
#define dimMAX 1000
int dim=1;

//Fitness size - number of sequences used to evaluate fitness
#define fez 100
char trd[fez][dimMAX];  //trial sequences used to evaluate fitness

static char datafyle[1000];  //the data file name

//minimum distance and its square
double mind;

//run controls - number of runs, mating events per run and rerporting interval
#define runs 1
#define mevs 100000
#define RI 1000
//running report
#define verbose 1

//population specifiers, size, random material rate, upper bound on code size
#define popsizeMAX 1000
int name;
int popsize;              //this holds the actual popsize
int RMR;
#define MAX 1000
//conway crossover buffer size
#define CCB 2500

//Conway variables
char pot[CCB][dimMAX];  //buffer for crossover
char res[CCB][dimMAX];  //result for crossover
int Cz;  //resulting crossover size

//population variables
char pop[popsizeMAX][MAX][dimMAX];    //holds population members
int csize[popsizeMAX];             //current code sizes
double fit[popsizeMAX];            //holds fitness value
int dx[popsizeMAX];                //sorting index

//Data variables - maximum number of data points
#define maxp 100000
static char *D[maxp];              //data store
int used[maxp];                    //is this data point used?
int ndp;                           //number of actual data points

void readdata();  //read in the data used as the "random points"
void initalg();   //initialize the algorithm
int okay(char *pt,double code[MAX][dimMAX],int z);  //is a new point okay? 
void rpoint(char *pt);                    //generate a random point
void cpoint(char *a,char *b);             //copy a point
void spoint(char *a,char *b);             //copy a point
void create(char gene[MAX][dimMAX],int &z);  //make an initial population member
double HD(char *a,char *b);               //Hamming distance between points
void initpop();                           //intialize a population
void ccross(int a,int b);                 //Conway crossover
void matingevent();                       //perform a mating event
void report(ostream &aus);                //make a statistical report
void reportbest(ostream &aus);            //report the best structure

//datafyle

int main(int argc,char **argv){//main program

int run,mev;         //run and mating event look counters
fstream stat,crit;   //statistics and structure reporting
char fn[60];         //file name construction buffer

  if(argc!=7){//check for command line arguments
    cout << "EvolverCL  <datafile>  <string length>  <population size>  <minimum distance>  <random material rate>" << endl;
    return(0);  //Remind the user if they are missing
  }

         /*******Pull in the command line arguments************/
  strcpy(datafyle,argv[1]);
  dim=atoi(argv[2]);       //put the data file name where it is expected
  name=atoi(argv[3]);  //name of subcluster
  popsize=atoi(argv[4]);     //get the population size
  if(popsize>popsizeMAX){//safety check
    cout << "The maximum population size is " << popsizeMAX << " please adjust this to use a population size of " << popsize << endl;
  }
  //Get the minimum distance
  mind=atof(argv[5]);
  //Get the random material rate
  RMR=atoi(argv[6]);
                   /*****Command line read******/

  initalg();      //initialize the algorithm

  char str1[60] = "best.pcld_";
  strcat(str1,argv[2]); 
  char separator[] = "_";
  strcat(str1, separator);
  strcat(str1,argv[3]);
  crit.open(str1,ios::out);  //open the good structure channel
  for(run=0;run<runs;run++){//loop over runs
    snprintf(fn, sizeof(fn), "run%03d.dat", run);   //create statistics filename
    stat.open(fn,ios::out);  //open the statistics channel
    initpop();      //initialize a population
    report(stat);  //make an intial report
    for(mev=0;mev<mevs;mev++){//loop over mating events
      matingevent();  //do a mating event
      if((mev+1)%RI==0){//make a report if its time
        if(verbose==1)cout << run << " " << (mev+1)/RI << " ";
        report(stat);  //to the report
      }
    }
    stat.close();  //close the statistics channel
    reportbest(crit);  //report the best structure
  }

  crit.close(); //close the reporting channel

  return(0);  //keep the system happy

}

void readdata(){//read in the data used as the "random points"

char buf[2000];    //character input buffer
fstream inp;       //input stream
int k;             //character index
int i;             //loop index

  for(i=0;i<maxp;i++)used[i]=0;      //nothing is used initially
  inp.open(datafyle,ios::in);        //open the data file
  ndp=0;                             //zero the number of data points
  inp.getline(buf,1999);             //get the labels
  while((strlen(buf)>3)&&(ndp<maxp)){//loop over data records
    D[ndp]=new char[dimMAX];            //get storage for the point
    inp.getline(D[ndp],1999);        //get the next data record
    //for(i=0;i<dim;i++)cout << D[ndp][i];cout << endl; //DEBUG
    ndp++;                           //record that the point exists
    inp.getline(buf,1999);           //get the next label
  }
  inp.close();
  cout << "Data file contains " << ndp << " items." << endl;
}

void initalg(){//initialize the algorithm

  srand48(RNS); //seed the random number generator
  readdata();   //read in the data

}

double HD(char *a,char *b){//find the hamming distance between points

double delta,accu;    //distance scratch variables
int i;                //loop index

  accu=0.0;  //zero the accumulator
  for(i=0;i<dim;i++){//loop over the coordinates
    if(a[i]!=b[i])accu+=1.0;  //accumulate the Hamming distance
  }
  return(accu);  //report the Hamming distance
}

int okay(char *pt,char code[MAX][dimMAX],int z){//is a new point okay? 

int i,j;              //loop index variables

  for(i=0;i<z;i++){//loop over poitns currently in code
    if(HD(code[i],pt)<mind)return(0); //too close
  }

  return(1); //survived all checks

}

void rpoint(char *pt){//generate a random point

int rp;   //random position

  do {//random point not marked as used already
    rp=lrand48()%ndp;  //select a random data point
  }while(used[rp]==1); //check for used...
  for(int i=0;i<dim;i++)pt[i]=D[rp][i];  //and copy it

}

void cpoint(char *a,char *b){//copy a point

  for(int i=0;i<dim;i++)a[i]=b[i];  //copy the values

}

void spoint(char *a,char *b){//swap points

char sw; //swap variables

 for(int i=0;i<dim;i++){sw=a[i];a[i]=b[i];b[i]=sw;}  //swap the values

}

//This routine implements the Conway algorithm to make new population members
void create(char gene[MAX][dimMAX],int &z){//make an initial population member

char pt[dimMAX];  //tentative new point
int i,j;         //loop index variables

  z=0;
  for(i=0;i<MAX;i++){//loop over samples
    rpoint(pt); //make a random point
    if(okay(pt,gene,z)){//found an acceptable point
      for(j=0;j<dim;j++)gene[z][j]=pt[j];  //transfer the point
      z++;  //register its existence
    }
  }
}

void initpop(){//intialize a population

int i;  //loop index variable

  for(i=0;i<popsize;i++){//loop over population
    create(pop[i],csize[i]);    //go make the population member
    fit[i]=((double)csize[i]);  //convert the fitness to a real number
    dx[i]=i;                    //refresh the sorting index
    cout << fit[i] << " ";
  }
  cout << endl;
}

void ccross(int a,int b){//perform Conway crossover of popmembers a and b

int i,j,k;       //index variables
char pt[dimMAX];  //new point

  //Assemble the buffer
  k=0;  //number of things in the pot buffer
  for(i=0;i<csize[a];i++){//copy the first guy
    cpoint(pot[k++],pop[a][i]);
  }
  for(i=0;i<csize[b];i++){//copy the first guy
    cpoint(pot[k++],pop[b][i]);
  }
  for(i=0;i<RMR;i++){//add the random material
    rpoint(pt);
    cpoint(pot[k++],pt);
  }

  //Shuffle the buffer
  for(i=0;i<k;i++){//loop over points in buffer
    j=lrand48()%k; //select random coordinate
    if(i!=j)spoint(pot[j],pot[i]);  //swap them if they are different
  }
  
  //Ready for lexicode algorithm
  Cz=0;  //zero the size
  for(i=0;i<k;i++)if(okay(pot[i],res,Cz)){//if the point fits
    cpoint(res[Cz],pot[i]);  //copy it in
    Cz++;  //and record its existence
  }

}

void matingevent(){//perform a mating event

int i;   //index variable

  tselect(fit,dx,3,popsize);  //select the three guys
  ccross(dx[1],dx[2]);        //perform conway crossover
  if(Cz>=csize[dx[0]]){//result at least as good?
    csize[dx[0]]=Cz;  //update size
    for(i=0;i<Cz;i++)cpoint(pop[dx[0]][i],res[i]); //copy the code
    fit[dx[0]]=((double)Cz); //update fitness
  }
 
}

void report(ostream &aus){//make a statistical report

dset D;

  D.add(fit,popsize);  //build data set
  if(verbose)cout << D.Rmu() << " " << D.RCI95() << " " 
                  << D.Rsg() << " " << D.Rmax() << endl;
  aus << D.Rmu() << " " << D.RCI95() << " "                //print out
      << D.Rsg() << " " << D.Rmax() << endl;

}

int equal(char *a,char *b){//are these two words equal?

int i;   //loop index

 for(i=0;i<dim;i++)if(a[i]!=b[i])return(0);  //Mismatch? Not equal!
 return(1);  //No mismatch, equal!
}

void reportbest(ostream &aus){//report the best structure

int i,j,k,b;

  b=0;  //initialize best pointer
  for(i=1;i<popsize;i++)if(fit[i]>fit[b])b=i;
  aus << fit[b] << " -fitness" << endl;
  k=0;
  for(i=0;i<csize[b];i++){//loop over vectors
    aus << pop[b][i][0];
    for(j=1;j<dim;j++)aus << pop[b][i][j];
    aus << endl;
    //scan and mark as used
    for(k=0;k<ndp;k++)if(equal(pop[b][i],D[k])){//used in code?
      cout << pop[b][i][0];
      for(j=1;j<dim;j++)cout << pop[b][i][j];
      cout << " -used!" << endl;
      used[k]=1;
    }
  }
}
