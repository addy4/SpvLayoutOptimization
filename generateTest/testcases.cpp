#include</Users/PremBhatia1/stdc.h>
#include<stdlib.h> 
#include<stdio.h> 
#include<string.h> 
#include<iostream> 
#include<time.h> 
#include<math.h> 
#include<cmath> 

using namespace std ; 

void generate(int X, int Y, int P, int Avg,int **arr, int total)
{
    int counts = 0 ; 
    srand(time(0)) ; 
    int block = (int)sqrt((double)P) ; 
    block = 5 ;  
    int Boxes[block][block] ; 
    // divided the area into Boxes -> #block by #block . As of now, each element Boxes has O panels 
    for(int i = 0 ; i < block ; i++)
    {
        for(int j = 0 ; j < block ; j++)
        {
            Boxes[i][j] = 0 ; 
        }
    }
    //Limit = No. of Clusters * Average 
    int Limit = total ; 
    int panels = 0 ; 
    // For mapping a panel to Boxes entry 
    double H = (double)block/X ; /* For determining the X index of Boxes */ 
    double V = (double)block/Y ; /* For determining the Y index of Boxes */
    //cout << "H = " << H << " " << "V = " << V << endl ; 
    while(panels < Limit)
    {
        int row = rand()%Y ; 
        int col = rand()%X ; 
        //if a row,col is such that arr[row][col] != 1 means that we have not made it 1, it is being made 1 when we are printing it as a panel
        //once it is made 1 , it will not be made 1 in any of the next iterations 
        //whenever a row/col is printed , panels++ and arr[row][col] = 1 and thus same row col 
        //cannot be made 1 in the next iteration due to the check arr[row][col] != 1  
        if(arr[row][col] == 0  /* && Boxes[(int)(row*H)][(int)(col*V)] < Avg */) 
        {
            //cout << col << " " << row << " Belongs to " << (int)(col*H) << " " << (int)(row*V) << endl ;
            //cout << col << " " << row << endl ; 
            //Boxes[(int)(row*H)][(int)(col*V)]++ ;  
            arr[row][col] = 1 ;
            cout << col << " " << row << endl ;  
            panels++ ; 
        } 
       //panels++ ; 
    }

    /*
    for(int i = 0 ; i < block ; i++)
    {
        for(int j = 0 ; j < block ; j++)
        {
            cout << Boxes[i][j] << " " << endl ; 
        }
    }
    */ 

    /* 
    for(int ii = 0 ; ii < Y ; ii++)
    {
        for(int jj = 0 ; jj < X ; jj++)
        {
            cout << arr[ii][jj] << " " ;
            if(arr[ii][jj] == 1)
            {
                counts++ ; 
            }
        }
        cout << endl ; 
    }
    */  
}

int main(int argc, char** argv)
{
    int X,Y,P,Avg,Total ; 
    X = 75 ; //length of field
    Y = 91 ; //width of field
    P = 12 ; //#medians
    Total = 315 ; //#spv arrays = N
    Avg = Total/P + 10000 ;    
    int** arr = (int**)malloc(Y*sizeof(int*)) ; 
    for(int l = 0 ; l < Y ; l++)
    {
        arr[l] = (int*)calloc(X,sizeof(int)) ; 
    }
    cout << X << " " << Y << endl ; 
    cout << Total << endl ; 
    cout << P << endl ; 
    generate(X,Y,P,Avg,arr,Total) ; 
}
