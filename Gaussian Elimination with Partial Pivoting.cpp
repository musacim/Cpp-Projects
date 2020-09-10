#include <iostream>
#include<fstream>
#include <new>
#include <cmath>
#include <stdlib.h>
using namespace std;
double AbsVal(double a)	////////////////////we have absolute value finder function here
{
	if (a>0)
	{return a;}
	else
	return a*-1;
}

int main(int argc, char *argv[])
{
	
////////////////////////////////////////we will start by getting the name of file from the user for matrix A
	
	
	string filename;    //we create string variable to store name of file 
	filename=argv[1];   
		
	
	ifstream myfile;  /////// we open file to read it with ifstream
    
	myfile.open(filename.c_str());  //// we open

 
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


	
	myfile.open(filename.c_str());//////////////////// we open file again to take the elements in it
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)		// we use 2 for loop for 2d array and temporarily store elements in a first then assign it to specific matrix position
			{	
			myfile>>a;
			arr[i][j]=a;
			}	
	}
	
	myfile.close();////////////// we took the data from file and closed it
		
	




			//----------------------------------------   condition number part     --------------------------
		if(n==2) //// we check the matrix is it 2 by 2?, ýf yes we will compute condition numbers
		{
			double inv_arr[2][2],det; /// we will find the invers of matrix first
			
			inv_arr[0][0]=arr[1][1];	/// we completed step one for inverse calculation manually because just n=2 case is wanted;
			inv_arr[1][1]=arr[0][0];
			inv_arr[0][1]=-arr[0][1];
			inv_arr[1][0]=-arr[1][0];
			
			det=(arr[0][0]*arr[1][1])-(arr[0][1]*arr[1][0]);//// we need to calculate determinant, we will divide all elements by determinant in next step
			
			for(i=0;i<2;i++)			//// we completed inverse of 2 by 2 matrix
			{
				for(j=0;j<2;j++)
				{
				inv_arr[i][j]/=det;
				}
			}	
			
			
			
			
			
				/// now we have inverse matrix, so we will compute norms of them at 1 and infinity;
					
				// we first find norm at 1 for matrix and then for the inverse of it,next repeat process at infinity	
					
				double row1,row2,col1,col2,norm_col,norm_row,norm_col_inv,norm_row_inv;
		
				col1=AbsVal(arr[0][0])+AbsVal(arr[1][0]);/// we compute norm of matrix at 1
				col2=AbsVal(arr[0][1])+AbsVal(arr[1][1]);
				if(col1>=col2)
				{
					norm_col=col1;
				}
				else
				{ norm_col=col2;}
				
				
				col1=AbsVal(inv_arr[0][0])+AbsVal(inv_arr[1][0]);/// we compute norm of inverse matrix at 1
				col2=AbsVal(inv_arr[0][1])+AbsVal(inv_arr[1][1]);
				if(col1>=col2)
				{
					norm_col_inv=col1;
				}
				else
				{ norm_col_inv=col2;}
				
				cout<<"\nCondition number of Matrix at 1 is: "<<norm_col*norm_col_inv<<endl;
				
				
				///// now we will compute norm at infinity
				
				
				row1=AbsVal(arr[0][0])+AbsVal(arr[0][1]);/// we compute norm of matrix at infinity
				row2=AbsVal(arr[1][0])+AbsVal(arr[1][1]);
				if(row1>=row2)
				{
					norm_row=row1;
				}
				else
				{ norm_row=row2;}
				
				
				row1=AbsVal(inv_arr[0][0])+AbsVal(inv_arr[0][1]);/// we compute norm of inverse matrix at infinity
				row2=AbsVal(inv_arr[1][0])+AbsVal(inv_arr[1][1]);
				if(row1>=row2)
				{
					norm_row_inv=row1;
				}
				else
				{ norm_row_inv=row2;}
				
				cout<<"Condition number of Matrix at infinity is: "<<norm_row*norm_row_inv<<endl; //result is multiple of this two according to formula
				cout<<"\n";
				
				
				
				
				
				
		}
			
			


////////////////////////////////////////////////////////////we now take matrix b
	string filename2;
	filename2=argv[2];
		
	int m =0;// we use m variable when counting elements in vector
	double store_element_b;

	ifstream myfile2;
    
	myfile2.open(filename2.c_str());//////////////////////open it
	
	double *array_b = new (nothrow)double[n];////////////////we create an array with dynamic memory allocation we used n because ýt should be n, otherwise it is not
															// possible to multiply this vector with n by n matrix

	while(myfile2>>store_element_b)/////////////////////// we take the data from file
	{
		
			array_b[m]=store_element_b;
		
		m++;
		
		
	}
	myfile2.close();//////////we close it


////////////////////////////////////////////Forward Substitution Part//////////////////////////////////////


		int pivot_num=0;		///we will sort the matrix n-1 times so I have to use pivot_num to follow number of calculations
		//we will make n-1 step calculation so we should create for loop 
	for(pivot_num =0;pivot_num<n-1;pivot_num ++)
    {
	      
				
					
			////////////////////////////////////////////////////  first I will find the max of the firts element among rows
					
					double max=-1;
					int index_row_of_max;
						
					for (i=pivot_num;i<n;i++)
					{
						if (AbsVal(arr[i][pivot_num]) >max)         ////    AbsVal()  is a function that I defined, it makes the element non-negative
						{
						
						max=AbsVal(arr[i][pivot_num]);
						index_row_of_max=i;
						}
					}
					
					
			     
			//////////////////////////////////////////////////////// Then change the order of the rows that max element will be at the top
																	
								
					
														
				 			
							 
							 
						
					double temp,temp2;	 				
					
					for(i=0;i<n;i++) //////////////////// we changed the order of rows of Matrix A
					{
						temp=arr[pivot_num][i];
						arr[pivot_num][i]=arr[index_row_of_max][i];
						arr[index_row_of_max][i]=temp;
					}
									
																				
					
			////////////////////////////////////////////////////////we changed the order of rows of Vector b
					
						temp2=array_b [pivot_num];
						array_b[pivot_num]=array_b[index_row_of_max];
						array_b[index_row_of_max]=temp2;
				
				
									////////////////////////////////	we divide elements of first column other than current pivot element by max of first column then	
									///////////////////////////////   multiply with first row then substract it from the other and repeat process for the rest of rows
									//////////////////////////////    this is forward substitiuon part of partial pivoting
				
			
			
					double pivot_element,coefficient,temp3,temp4;
					pivot_element=arr[pivot_num ][pivot_num];
						
									///// we initiliaze the i with piv_num+1 because elements below the pivot will be processed
					for(i=pivot_num+1;i<n;i++)	
					{
							
						coefficient=arr[i][pivot_num]/pivot_element;    ///// we multiply the first row with coefficients and then substract from the other rows
						for(j=0;j<n;j++)
						{
						
						temp3=arr[pivot_num][j] *coefficient;
						arr[i][j]-=temp3;
								
								
							
						}
						temp4=array_b[pivot_num]*coefficient;  //now we find the b matrix according to calculations we made previously, it is out of j loop since it has one column				
						array_b[i]-=temp4;
						
			 	   }
							
							
															
																
					

							
	}
		double control_num;			
		int count=0;///////////////////////////we are checking singularity here, if square matrix is full rank then it is nonsingular
		for(i=0;i<n;i++)
		{
		
			
			control_num=round(arr[i][i]*10000)/10000;///// here we set the precision 4 digits by using round func.
		
			
			if (control_num!=0 && isnan(arr[i][i])==0 && isinf(arr[i][i])==0) /// we also check if nan or inf occurs during calculation
			{
			count++;}   ////// we count nonzero diagonals here then if it is not equal to dimension of matrix n,it is singular
		}
		if (count!=n)
		{
			cout<<"Error:\nThis matrix is singular, system is stopped. "<<endl;//////error message is printed here
			exit(0);   ///////we quit
		}
							
						

	
			
	
			
/////////////////////////////////////////////////Back Substitution /////////////////////////////
		
			
			double X[n],temp5=0;		///we create empty X vector
			
			X[n-1]=array_b[n-1]/arr[n-1][n-1]; ////we can calculate easily the last element of X 
			
			for(i=n-2;i>=0;i--)//////////// we have back substitiuon algorithm here
			{
				
				
				temp5=array_b[i];
				for(j=i+1;j<n;j++)
				
				{
					
					
					temp5-=arr[i][j]*X[j]; //// we are solving the equation here for each varible
				}
				
					X[i]=temp5/arr[i][i];
				
			
			}
			
		
	
			cout<<"\nX vector is: "<<endl;/// print out the x vector
			for(i=0;i<n;i++)
			{
				cout<<X[i]<<endl;
			}
		
		
			
			
			
			
			
			
			ofstream my_X_file("X.txt");		//////we create a file named "X" then write X vector in it, it is created into working directiory folder
			
			for(i=0;i<n;i++)
			{
				my_X_file<<X[i]<<"\n";
			}
		
			
			
			
			
			
			
			
		
		  delete []arr;			////we finished, so clear the pointers
		  delete []array_b;
			
		



	return 0;
}
