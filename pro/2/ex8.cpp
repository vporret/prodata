
#include <iostream> 

using namespace std; 

// Recursive function

int g_c_d(int a, int b) 
    
{ 
        
    if (a <= 0 || b <= 0) // a parameter is not positive
        return 0; // return 0
   
    // if equal
    if (a == b) 
        return a; 
   
    // if a is bigger
    
    if (a > b) 
        return g_c_d(a-b, b); 
    return g_c_d(a, b-a); 
} 

// Non-recursive function

int g_c_d2(int c, int d)     	
{
   int  number;                
   while (d != 0) 
   {      
      number = c % d;          
      c = d;            
      d = number;
   }                      
   return c;              
} 

// Test function : 

int main() 
{ 
    int a = 155, b = 5; 
    cout<<"GCD '(recursive)' of "<<a<<" and "<<b<<" is "<< g_c_d(a, b)<<"\n"; 
    int c = 155, d = 5; 
    cout<<"GCD2 '(non-recursive) of "<<a<<" and "<<b<<" is "<< g_c_d2(a, b); 
    return 0; 
    
} 

