#include <iostream>
#include <stdio.h>
#include<stdlib.h>
#include<fstream>
#include<math.h>
using namespace std;



double absolute(double a)  //  We create a function to take absolute value of the expressions
{
	if (a<0)
		return -1*a;
	else
		return a;
}



double findmax(double *xvector,int size){   // findmax function will be used in power iteration process, it finds the max value of the vector given

	double temp=-1;
		for(int i=0;i<size;i++)
		{
	 	
			if(absolute(xvector[i])>temp)
				temp=absolute(xvector[i]);
		
		}

return temp;
}



double *dividebymax(double *vector,int size,double max)  // this function normalize the vector given, by dividing it the max of last iteration
{	
	for(int i=0;i<size;i++)
		{
	 	
			vector[i]/=max;
		
		}

return vector;
}



double dotproduct (double *vect,double *vect2,int size) // It calculates the dotproduct we will use it in householder transformation
{	
	
	double sum=0;
	for (int i=0;i<size;i++)
	{
		 
	 sum+=vect[i]*vect2[i];
	}	 
	
	
return sum;
} 


	 
double** vectormultiply (double *vect,int size,double coeff)  // this function computes nxn matrix from nx1 vector namely;  V * transpose(V)
{	
	
	double sum=0;
	double **newvect=new double*[size];
	for (int i=0;i<size;i++)
		newvect[i]=new double[size];
	
	for (int i=0;i<size;i++)
	{
		 
		for (int j=0;j<size;j++)
			{
				sum=0;
				sum+= vect[i] * vect[j];
				newvect[i][j]=sum*coeff;
			}
		
	}	 
	
	
return newvect;
} 



double **MatrixSubtraction (double **Matrix1,double **Matrix2,int size )   // we subtract first matrix given from the second
{	
	
	double sum=0;
	double **newmatrix=new double*[size];
	for (int i=0;i<size;i++)
		newmatrix[i]=new double[size];
	
	for (int i=0;i<size;i++)
	{
		 
		for (int j=0;j<size;j++)
			{
				 
				newmatrix[i][j]=Matrix1[i][j]-Matrix2[i][j];
			}
		
	}	 
	
	
return newmatrix;	
} 



double** transpose(double** matrix,int size  )  // we take the transpose of given matrix
{
	 
	double** tempmatrix ;
	for(int i=0;i<size;i++)
	{
		tempmatrix[i]=new double[size];
	}
	
	
	for(int i=0;i<size;i++)
		{	  	
			for(int j=0;j<size;j++)
				{
				tempmatrix[i][j]=1;
				}
		}
	
	for(int i=0;i<size;i++)
		{	  	
			for(int j=0;j<size;j++)
				{
				tempmatrix[j][i]=matrix[i][j];
				}
		}
	

return tempmatrix;
}



double** matrixmultiplication(double** matrix1,double** matrix2,int size)   //we calculate multiplication of given 2 matrices
{
	double** tempmatrix=new double*[size],sum;
	for(int i=0;i<size;i++)
	{
		tempmatrix[i]=new double[size];
	}
	
	for(int i=0;i<size;i++)
	{
		
		for(int j=0;j<size;j++)
			{
		     sum=0;
			for(int k=0;k<size;k++)
				{
					sum+=*&matrix1[i][k]*(*&matrix2[k][j]);
	 
					  
				}
				tempmatrix[i][j]=sum;
	
			}
   }
   
return tempmatrix;	
}





class Matrix{      // now we create class named matrix to calculate power iteration
	
	
	public:
	
	double **matrix,*xvector,*xnew,eigenvalue,*eigenvector,**H,**HT;
	int size;
	
	
	Matrix(int n,double **m)  // constructor: we take matrix and size of it
	{
	 size=n;
	 
	 
	 matrix=new double*[n];
	 for(int i=0;i<n;i++){
	  
	 	matrix[i]=new double[n];
	 }
	
	 matrix=m;
	}
	
 
	
	double* initializeXvector()   // Initiliazer: it creates a vector with the largest value is 1 others is 0
	{
		double *xvector=new double[size];
		xvector[0]=1;
		return xvector;
	}
	

	
	double* multiply(double *newxvector) // Multiply and get Xnew, it computes Xnew=A*X
	{	
		
		double sum=0;
		xnew=new double[size];
		for (int i=0;i<size;i++)
		{
			sum=0;
			for (int j=0;j<size;j++)
				{
					
					sum+=matrix[i][j]*newxvector[j];
					
				}
			xnew[i]=sum;
		}	 
		
		
	return xnew;
		
	}
	int PowerMethod(double tolerans)   // this is the main method, we initialize and continue to multiply the vector and matrix until convergence occurs
	{
		int negative=0,times=0;
		double difference=0,max=0,tempmax=0;
		double *xvector=new double[size];
		double *X0=new double[size];
		double *tempx=new double[size],*tempx2=new double[size];
		
		X0=initializeXvector(); // initialize some vector
		tempx=X0;
	 
		do
		{
		tempmax=max;
		xvector=multiply(tempx); // A*X
		max=findmax(xvector,size); // find max, candidate of eigenvector
		tempx=dividebymax(xvector,size,max); //normalize it
		 
		if(tempx[0]*tempx2[0]<0)   // this where we check is eigenvalue is negative or not, we track signs every iteration and check
			negative++;
			times++;
		tempx2=*&tempx;   
		 
		difference=absolute(tempmax-max); // we check convergence here if it is less then tolerans value iteration stops
		} while(difference>tolerans);
		
		eigenvalue=max;
		eigenvector=tempx;
		if(negative>times/2) // if in half of the iteration signs changes that means we have negative eigenvalue
			negative=1;
	return negative;
    }
	
	 
	double **HouseholdTransform(double *eigenvector)  // we need H, househodler matrix in order to use deflation process
	{
	
		double total=0,*vect,**matr,**identity, coeff ;
		
		identity=new double*[size];
	 
		for(int i=0;i<size;i++)  // we create identity matrix
		{
		  
		 	identity[i]=new double[size];
	 
		}
		 
		for(int i=0;i<size;i++) // assigns 1's to identity
		{
		  for(int j=0;j<size;j++)
			{
				if(i==j)
				 	identity[i][j]=1;
				else 
					identity[i][j]=0;	 
			}
		}
		 

		for (int i=0;i<size;i++)
			{
			total+=eigenvector[i]*eigenvector[i];
			}
		total=sqrt(total);
		vect=eigenvector;
		vect[0]=vect[0]-total; // this is the Hx1=ce1  process 
		
		coeff=2/dotproduct(vect,vect,size);
		
		H=MatrixSubtraction(identity,vectormultiply(vect,size,coeff),size);// we computed the H and it is stored inside class so no return
	 
	 
	
	
	
		
	}

double** deflation()  // in deflation part we get rid of from the largest eigen value and compute second largest one
	{
		
		double **result=new double*[size-1],**temp,**AandH,**HT;
		 
	 
		for(int i=0;i<size-1;i++)
		{
		 	result[i]=new double[size-1];
		} 
		
		HT=transpose(H,size);                        //  transpose(H)
		AandH=matrixmultiplication(matrix,HT,size);  //  A*H'
	 	temp=matrixmultiplication(H,AandH,size);     //  H*A*H'
	 
	    for(int i=0;i<size;i++)
		{
		 	for(int j=0;j<size;j++)
				{
				 	 if(j>0 and i>0)
				 	 	result[i-1][j-1]=temp[i][j];   // We take the B matrix here, it is other than first column and row
				} 
			 
		} 		
		return result;

	}
};





int main( int argc, char *argv[] )  {
 
   
   if( argc >4 )                                // we check input arguments 
   {
      cout<<"Too many arguments supplied";
   }
   else if (argc<4)
   {
      cout<<"more argument expected";
   }

 
 
 	double tolerans;
 	tolerans= atof(argv[2]);

	
	ifstream myfile;  /////// we open file to read it with ifstream
   	myfile.open(argv[1]);  //// we open

 
	int n; /// n will be dimension of n by n matrix
	double a; // we create a variable to store elements of matrix
///////////////////////////// we should compute the size of matrix in order to see dimensions of matrix we take the position data 
///////////////////////////// we count the elements in matrix
	while(myfile>>a)	
	{
		n++;
		
	}
	myfile.close();  //////////////// we close file
	n=sqrt(n);/// square root of number of elements in square matrix gives us a dimension of matrix
	//now n is size of array
	///////////////////////////////////////we now create 2D array with dynamic memory allocation
    int i,j;
	double **arr=new double*[n];  //// we use pointer of pointer to create 2d array
	for(i=0;i<n;i++)
    {
    	arr[i]=new double[n];
	}

	int size=n;
	
	myfile.open(argv[1]);//////////////////// we open file again to take the elements in it
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)		// we use 2 for loop for 2d array and temporarily store elements in a first then assign it to specific matrix position
			{	
			myfile>>a;
			arr[i][j]=a;
			}	
	}
	
	myfile.close();////////////// we took the data from file and closed it
		
    //Power Iteration For 1st Eigenvalue
	int negative=0,negativ2=0;
 	
	Matrix MatrixGiven (size,arr); // create class object of given matrix
    negative=MatrixGiven.PowerMethod(tolerans);  // call power method find eigenvalue
  	
	if(negative==1)  // check negativity
  		MatrixGiven.eigenvalue=MatrixGiven.eigenvalue*-1;	
	
	ofstream my_X_file(argv[3]);	  // write the results into output file
	my_X_file<<"Eigenvalue#1: ";
	my_X_file<<round(MatrixGiven.eigenvalue*100)/100<<endl;
	for(i=0;i<n;i++)
	{
		my_X_file<<round(MatrixGiven.eigenvector[i]*100)/100<<endl;
	}
 	
	 
	//Power Iteration For 2nd Eigenvalue and Deflation
	MatrixGiven.HouseholdTransform(MatrixGiven.eigenvector); // call householder and compute H
	double **B; 
    B=MatrixGiven.deflation();  	 // call deflation and compute matrix B
	
	Matrix MatrixB(size-1,B); // create class object of matrix B
	negativ2=MatrixB.PowerMethod(tolerans); // call power method second time and find 2nd dominant eigenvalue
	if(negativ2==1) // check negativity
  		MatrixB.eigenvalue=MatrixB.eigenvalue*-1;
 

	// write results into file
	my_X_file<<"Eigenvalue#2: ";
	my_X_file<<round(MatrixB.eigenvalue*100)/100<<endl;


}
