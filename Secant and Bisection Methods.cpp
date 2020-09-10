#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

double abs(double a)  // we create absolute value function
{
	if(a<0)
		return a*-1;
	else
		return a;
}

double power(double number,int n) // we create function for taking power of numbers
{
	double sum;
	if (n==0)
		sum=1;
	else
	{
	sum=number;
	
	for (int i=1;i<n;i++)
	{
		sum*=number;
	}
	
	}
	return sum;
	
}
// we create a class which contains the methods
class Algorithms{
	
	public:
	double x1,x0,tol,*arg,sendx0,sendx1;
	int order,sendcount;
	Algorithms(int n,double *arguments)   // constructor takes inputs and the number of given parameters
	{
		
	 	arg=arguments;
	    tol=arg[n-1];
	    x1 =arg[n-2];
	    x0 =arg[n-3];
		order=n-5;	
 		 
	}
	
    
 	
	double f (double x)   // this is a function, we calculate the f(x) here with given coefficients
	{
			double result=0;
			
			for (int i=0;i<=order;i++)
			{
			
			result +=arg[i+1] * power(x,order-i) ;
		    
			}
			
			return result;
	}
	
	double Secant(bool enable_hybrid=false)  // secant algorithm, takes enable and if hybrid enable is "on" 
											  //we take x values from bisection method after 2 iteration is calculated in bisection
	{
		double Xnew,difference,X1,X0;
		if(enable_hybrid==false)  // hybrid enable is checked then, use initial x values accordingly
		{
			X0=x0;
			X1=x1;
		}
		else
		{
			X0=sendx0;  // we take inclass variable here 
			X1=sendx1;
		}
		int count_h;
		int count_s=0;  // we will count iterations, if hybrid is enable we will add extra iterations from bisection method
		
		do{
		
		count_s++;
	
		Xnew= X1 -(f(X1)* ((X1-X0)/(f(X1)-f(X0)) )); // this is the main algorithm we calculate xnew at every iteration until converges
		difference=Xnew-X1;
	 
		X0=X1;
		X1=Xnew;
		
		}while(abs(difference)>tol);  // check difference betwween last two xnew and stop iteration if needed
		count_h=count_s+sendcount;
	 	if(enable_hybrid==false) // we check enable again and print out the iteration number accordingly
			cout<<"Root is "<<Xnew<<endl<<"# of iterations is "<<count_s<<" for SECANT"<<endl;
		else
			cout<<"Root is "<<Xnew<<endl<<"# of iterations is "<<count_h<<" for HYBRID"<<endl;
	}
	
	double Bisection(bool enable_hybrid=false )  // bisection algorithm we send x values to secant if enable hybrid is "on" , if not calculate them and print root
	{
		double Xnew,temp,difference,X1,X0;
	
		int count_b=0;
		X0=x0;
		X1=x1;
		do
		{
		temp=Xnew;
		count_b++;
		
		if (count_b==3 and enable_hybrid==true) // we check iteration number here, if it reached three that means we have two,so break the loop and send information to secant
		{   
		    
			sendx0=X0;
			sendx1=X1;
			sendcount=2;
			return 0;

		}
	    
	    	
		Xnew=(X0+X1)/2;  // bisection algorithm divides into 2 at every iteration, and check the sign of multiplication of two options
		if(f(Xnew)*f(X0)<0)
		{
			X1=Xnew;
			X0=X0;
			 
		}
		else
		{
			X0=Xnew;
			X1=X1;
		 
		}
		difference=abs(temp-Xnew);
		
		if (difference<tol and count_b<3 and enable_hybrid==true)//in hybrid,we make first two iteration in bisection, we should check whether it is finished befor 2 iteration
			sendx0=X0;
			sendx1=X1;
			sendcount=count_b;
		}while(difference>tol); // we stop iterating if tolerance provided then print out the root and # of iteration
	 	if(enable_hybrid==false)
	 	
			cout<<"Root is "<<Xnew<<endl<<"# of iterations is "<<count_b<<" for BISECTION"<<endl;
		 
		
	}

	double Hybrid() // hybrid algorithm, send enable as "true" first bisection calculates fist 2 iteration then secant take new x values and finish the process
	{
		
		
		Bisection(true);
		Secant(true);
		
		
	}
};
int main(int argc, char** argv) {
	
	int n,i;
    n=argc;  // number of arguments we have
	double *arguments=new double[n];   // we create dynamically allocated array for inputs
 	
 	for(i=0;i<n;i++)
    {
    	arguments[i]=atof(argv[i]);  // we store them as a number
    	
	 
	}
	
	Algorithms Algos(n,arguments); // create a class which includes methods in it, then call these methods one by one
	
	Algos.Bisection();
	cout<<endl;
	Algos.Secant();
	cout<<endl;
	Algos.Hybrid();
	
	return 0;
}
