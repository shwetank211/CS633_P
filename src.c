// Including necessary libraries
#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include<math.h>

int main(int argc, char* argv[])
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    int myrank,size_total;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&size_total);

    // Extracting parameters from command line arguments
    int N=sqrt((double)atoi(argv[2]));
    int t=atoi(argv[3]);
    int seed=atoi(argv[4]);
    int stencil_type=atoi(argv[5]);

    // Initializing the array 'old' with random values
    double* old;
    old=(double*)malloc(N*N*sizeof(double));
    srand(seed*(myrank+10));
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            old[i*N+j]=abs(rand()+(i*rand()+j*myrank))/100;
        }
    }
    // Determining the process position in the grid
    int Pi,Pj;
    int nr=atoi(argv[1]),nc=12/nr;
    Pi = myrank/nc;
    Pj = myrank%nc;
    int rsize=nc,csize=nr;
    int lflag=0,rflag=0,uflag=0,dflag = 0;
    if(Pi>0)
        uflag=1;
    if(Pi<csize-1)
        dflag=1;
    if(Pj>0)
        lflag=1;
    if(Pj<rsize-1)
        rflag=1;

    double sTime,eTime,maxTime;// Declaring variables for timing
    if(stencil_type==9)//9 point stencil
    {
        double left[N],right[N],up[N],down[N];//arrays for storing boundary values in each direction
        double send_left[2*N],send_right[2*N],send_up[2*N],send_down[2*N];//arrays for sending boundary values
        double recv_left[2*N],recv_right[2*N],recv_up[2*N],recv_down[2*N];//arrays for receiving boundary values
        double left1[N],right1[N],up1[N],down1[N];//arrays for storing penultimate boundary values in each direction
        for(int i=0;i<N;i++)//initializing boundary value arrays to zero
        {
            left[i]=0;
            right[i]=0;
            up[i]=0;
            down[i]=0;
            left1[i]=0;
            right1[i]=0;
            up1[i]=0;
            down1[i]=0;
        }
        double* new;//allocating memory for the new array and initializing it with values from old
        new=(double*)malloc(N*N*sizeof(double));
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                new[i*N+j]=old[i*N+j];
            }
        }
        int* den;//allocating memory for the denominator matrix den and initializing it with 9s
        den=(int*)malloc(N*N*sizeof(int));//implicit type conversion will be used in division and space is optmised
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                den[i*N+j]=9;
            }
        }
        //adjusting denominator matrix for boundary conditions
        if(!uflag)
        {
            for(int i=0;i<N;i++)
            {
                den[0*N+i]-=2;
                den[1*N+i]--;
            }
        }
        if(!dflag)
        {
            for(int i=0;i<N;i++)
            {
                den[(N-1)*N+i]-=2;
                den[(N-2)*N+i]--;
            }
        }
        if(!lflag)
        {
            for(int i=0;i<N;i++)
            {
                den[i*N+0]-=2;
                den[i*N+1]--;
            }
        }
        if(!rflag)
        {
            for(int i=0;i<N;i++)
            {
                den[i*N+N-1]-=2;
                den[i*N+N-2]--;
            }
        }
        //timing the code
        sTime = MPI_Wtime();
        //loop for iterations
        while(t--)
        {
            int position=0;//initializing position for packing boundary values
            if(uflag) //if there is process above
            {
                for(int i=0;i<N;i++) //sending row 0
                {
                    MPI_Pack(old+i,1,MPI_DOUBLE,send_up,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                    MPI_Pack(old+i+N,1,MPI_DOUBLE,send_up,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_up,position,MPI_PACKED,myrank-nc,myrank,MPI_COMM_WORLD,&request);
            }
            if(dflag) //if there is a process below  
            {
                position=0;
                for(int i=0;i<N;i++) //sending row (N-1)
                {
                    MPI_Pack(old+(N-1)*N+i,1,MPI_DOUBLE,send_down,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                    MPI_Pack(old+(N-2)*N+i,1,MPI_DOUBLE,send_down,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_down,position,MPI_PACKED,myrank+nc,myrank,MPI_COMM_WORLD,&request);  
            }
            if(lflag) //if there is a process on left
            {
                position=0;
                for(int i=0;i<N;i++) //sending column 0
                {
                    MPI_Pack(old+i*N,1,MPI_DOUBLE,send_left,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                    MPI_Pack(old+i*N+1,1,MPI_DOUBLE,send_left,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_left,position,MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&request);
            }
            if(rflag) //if there is a process on right
            {
                position=0;
                for(int i=0;i<N;i++) //sending column (N-1)
                {
                    MPI_Pack(old+(i*N+N-1),1,MPI_DOUBLE,send_right,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                    MPI_Pack(old+(i*N+N-2),1,MPI_DOUBLE,send_right,2*N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_right,position,MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&request);          
            }
            //receiving data from neighbouring processes
            if(uflag) //receiving from above process
            {
                position=0;
                MPI_Recv(recv_up,2*N*sizeof(double),MPI_PACKED,myrank-nc,myrank-nc,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_up,2*N*sizeof(double),&position,up+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    MPI_Unpack(recv_up,2*N*sizeof(double),&position,up1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(dflag) //receiving from below process
            {
                position=0;
                MPI_Recv(recv_down,2*N*sizeof(double),MPI_PACKED,myrank+nc,myrank+nc,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_down,2*N*sizeof(double),&position,down+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    MPI_Unpack(recv_down,2*N*sizeof(double),&position,down1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(lflag) //receiving from left process
            {
                position=0;
                MPI_Recv(recv_left,2*N*sizeof(double),MPI_PACKED,myrank-1,myrank-1,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_left,2*N*sizeof(double),&position,left+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    MPI_Unpack(recv_left,2*N*sizeof(double),&position,left1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(rflag) //receiving from right process
            {
                position=0;
                MPI_Recv(recv_right,2*N*sizeof(double),MPI_PACKED,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_right,2*N*sizeof(double),&position,right+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                    MPI_Unpack(recv_right,2*N*sizeof(double),&position,right1+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            //updating interior points using the 9-point stencil computation
            for(int i=2;i<N-2;i++)
            {
                for(int j=2;j<N-2;j++)
                {
                    new[i*N+j]=(old[(i-1)*N+j]+old[(i+1)*N+j]+old[i*N+j-1]+old[i*N+j+1]+old[(i-2)*N+j]+old[(i+2)*N+j]+old[i*N+j-2]+old[i*N+j+2]+old[i*N+j])/9;
                }
            }
            //updating corner points using the boundary conditions
            new[0*N+0]=(up[0]+left[0]+up1[0]+left1[0]+old[0*N+1]+old[0*N+2]+old[1*N+0]+old[2*N+0]+old[0*N+0])/den[0*N+0];
            new[0*N+N-1]=(up[N-1]+up1[N-1]+right[0]+right1[0]+old[0*N+N-2]+old[1*N+N-1]+old[0*N+N-3]+old[2*N+N-1]+old[0*N+N-1])/den[0*N+N-1];
            new[(N-1)*N+0]=(down[0]+left[N-1]+old[(N-1)*N+1]+old[(N-2)*N+0]+down1[0]+left1[N-1]+old[(N-1)*N+2]+old[(N-3)*N+0]+old[(N-1)*N+0])/den[(N-1)*N+0];
            new[(N-1)*N+N-1]=(down[N-1]+right[N-1]+old[(N-2)*N+N-1]+old[(N-1)*N+N-2]+down1[N-1]+right1[N-1]+old[(N-3)*N+N-1]+old[(N-1)*N+N-3]+old[(N-1)*N+N-1])/den[(N-1)*N+N-1];
            new[0*N+1]=(up[1]+left[0]+up1[1]+old[0*N+0]+old[0*N+2]+old[0*N+3]+old[1*N+1]+old[2*N+1]+old[0*N+1])/den[0*N+1];
            new[0*N+N-2]=(up[N-2]+up1[N-2]+right[0]+old[0*N+N-3]+old[0*N+N-4]+old[1*N+N-2]+old[2*N+N-2]+old[0*N+N-2]+old[0*N+N-1])/den[0*N+N-2];
            new[(N-1)*N+1]=(down[1]+left[N-1]+old[(N-1)*N+2]+old[(N-2)*N+1]+down1[1]+old[(N-1)*N+0]+old[(N-1)*N+3]+old[(N-3)*N+1]+old[(N-1)*N+1])/den[(N-1)*N+1];
            new[(N-1)*N+N-2]=(down[N-2]+right[N-1]+old[(N-2)*N+N-2]+old[(N-1)*N+N-3]+down1[N-2]+old[(N-1)*N+N-1]+old[(N-3)*N+N-2]+old[(N-1)*N+N-4]+old[(N-1)*N+N-2])/den[(N-1)*N+N-2];                
            new[1*N+0]=(up[0]+left[1]+left1[1]+old[1*N+1]+old[0*N+0]+old[2*N+0]+old[3*N+0]+old[1*N+2]+old[1*N+0])/den[1*N+0];
            new[1*N+N-1]=(up[N-1]+right[1]+right1[1]+old[1*N+N-2]+old[0*N+N-1]+old[2*N+N-1]+old[3*N+N-1]+old[1*N+N-3]+old[1*N+N-1])/den[1*N+N-1];
            new[(N-2)*N+0]=(down[0]+left[N-2]+left1[N-2]+old[(N-2)*N+1]+old[(N-1)*N+0]+old[(N-3)*N+0]+old[(N-4)*N+0]+old[(N-2)*N+2]+old[(N-2)*N+0])/den[(N-2)*N+0];
            new[(N-2)*N+N-1]=(down[N-1]+right[N-2]+right1[N-2]+old[(N-2)*N+N-2]+old[(N-1)*N+N-1]+old[(N-3)*N+N-1]+old[(N-4)*N+N-1]+old[(N-2)*N+N-3]+old[(N-2)*N+N-1])/den[(N-2)*N+N-1];
            new[1*N+1]=(up[1]+left[1]+old[1]+old[1*N+0]+old[1*N+2]+old[1*N+3]+old[2*N+1]+old[3*N+1]+old[1*N+1])/den[1*N+1];
            new[1*N+N-2]=(up[N-2]+old[N-2]+right[1]+old[1*N+N-3]+old[1*N+N-4]+old[2*N+N-2]+old[3*N+N-2]+old[1*N+N-1]+old[1*N+N-2])/den[1*N+N-2];
            new[(N-2)*N+1]=(down[1]+old[(N-1)*N+1]+left[N-2]+old[(N-2)*N+0]+old[(N-2)*N+2]+old[(N-2)*N+3]+old[(N-3)*N+1]+old[(N-4)*N+1]+old[(N-2)*N+1])/den[(N-2)*N+1];
            new[(N-2)*N+N-2]=(down[N-2]+old[(N-1)*N+N-2]+right[N-2]+old[(N-2)*N+N-3]+old[(N-2)*N+N-4]+old[(N-3)*N+N-2]+old[(N-4)*N+N-2]+old[(N-2)*N+N-2]+old[(N-2)*N+N-1])/den[(N-2)*N+N-2];
            //updating boundary points (excluding corners)
            for(int i=2;i<N-2;i++)
            {
                new[0*N+i]= (old[0*N+i+1]+old[0*N+i-1]+old[1*N+i]+up[i]+old[0*N+i+2]+old[0*N+i-2]+old[2*N+i]+up1[i]+old[0*N+i])/den[0*N+i];
                new[(N-1)*N+i]= (old[(N-1)*N+i+1]+old[(N-1)*N+i-1]+old[(N-2)*N+i]+down[i]+old[(N-1)*N+i+2]+old[(N-1)*N+i-2]+old[(N-3)*N+i]+down1[i]+old[(N-1)*N+i])/den[(N-1)*N+i];
                new[i*N+0]= (old[(i+1)*N+0]+old[(i-1)*N+0]+old[i*N+1]+left[i]+old[(i+2)*N+0]+old[(i-2)*N+0]+old[i*N+2]+left1[i]+old[i*N+0])/den[i*N+0];
                new[i*N+N-1]=(old[(i+1)*N+N-1]+old[(i-1)*N+N-1]+old[i*N+N-2]+right[i]+old[(i+2)*N+N-1]+old[(i-2)*N+N-1]+old[i*N+N-3]+right1[i]+old[i*N+N-1])/den[i*N+N-1];
            }
            //updating penultimate boundary points (excluding corners)
            for(int i=2;i<N-2;i++)
            {
                new[1*N+i]= (old[1*N+i+1]+old[1*N+i-1]+old[2*N+i]+up[i]+old[1*N+i+2]+old[1*N+i-2]+old[3*N+i]+old[0*N+i]+old[1*N+i])/den[1*N+i];
                new[(N-2)*N+i]= (old[(N-2)*N+i+1]+old[(N-2)*N+i-1]+old[(N-1)*N+i]+down[i]+old[(N-2)*N+i+2]+old[(N-2)*N+i-2]+old[(N-3)*N+i]+old[(N-4)*N+i]+old[(N-2)*N+i])/den[(N-2)*N+i];
                new[i*N+1]= (old[(i+1)*N+1]+old[(i-1)*N+1]+old[i*N+2]+left[i]+old[(i+2)*N+1]+old[(i-2)*N+1]+old[i*N+3]+old[i*N+0]+old[i*N+1])/den[i*N+1];
                new[i*N+N-2]=(old[(i+1)*N+N-2]+old[(i-1)*N+N-2]+old[i*N+N-3]+right[i]+old[(i+2)*N+N-2]+old[(i-2)*N+N-2]+old[i*N+N-4]+old[i*N+N-1]+old[i*N+N-2])/den[i*N+N-2];
            }
            //updating old array with new array    
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    old[i*N+j]=new[i*N+j];
                }
            }
        }
        eTime=MPI_Wtime();
    }
    else //5 point stencil
    {
        double left[N],right[N],up[N],down[N];//arrays for storing boundary values in each direction
        double send_left[N],send_right[N],send_up[N],send_down[N];//arrays for sending boundary values
        double recv_left[N],recv_right[N],recv_up[N],recv_down[N];//arrays for receiving boundary values
        for(int i=0;i<N;i++)//initializing boundary value arrays to zero
        {
            left[i]=0;
            right[i]=0;
            up[i]=0;
            down[i]=0;
        }
        //allocating memory for the new array and initializing it with values from old
        double* new;
        new=(double*)malloc(N*N*sizeof(double));
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                new[i*N+j]=old[i*N+j];
            }
        }
        //allocating memory for the denominator matrix den and initializing it with 5s
        int* den;
        den=(int*)malloc(N*N*sizeof(int));
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                den[i*N+j]=5;
            }
        }
        //adjusting denominator matrix for boundary conditions
        if(!uflag)
        {
            for(int i=0;i<N;i++)
            {
                den[0*N+i]--;
            }
        }
        if(!dflag)
        {
            for(int i=0;i<N;i++)
            {
                den[(N-1)*N+i]--;
            }
        }
        if(!lflag)
        {
            for(int i=0;i<N;i++)
            {
                den[i*N+0]--;
            }
        }
        if(!rflag)
        {
            for(int i=0;i<N;i++)
            {
                den[i*N+N-1]--;
            }
        } 
        //timing the code
        sTime=MPI_Wtime();   
        while(t--)//loop for iterations
        {
            //initializing position for packing boundary values
            int position=0;
            if(uflag)//if there is process above
            {
                for(int i=0;i<N;i++) //sending row 0
                {
                    MPI_Pack(old+i,1,MPI_DOUBLE,send_up,N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_up,position,MPI_PACKED,myrank-nc,myrank,MPI_COMM_WORLD,&request);
            }
            if(dflag) //if there is a process below  
            {
                position=0;
                for(int i=0;i<N;i++) //sending row (N-1)
                {
                    MPI_Pack(old+(N-1)*N+i,1,MPI_DOUBLE,send_down,N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_down,position,MPI_PACKED,myrank+nc,myrank,MPI_COMM_WORLD,&request);
                
            }
            if(lflag) //if there is a process on left
            {
                position=0;
                for(int i=0;i<N;i++) //sending column 0
                {
                    MPI_Pack(old+i*N,1,MPI_DOUBLE,send_left,N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_left,position,MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&request);
            }
            if(rflag) //if there is a process on right
            {
                position=0;
                for(int i=0;i<N;i++) //sending column (n-1)
                {
                    MPI_Pack(old+(i*N+N-1),1,MPI_DOUBLE,send_right,N*sizeof(double),&position,MPI_COMM_WORLD);
                }
                MPI_Isend(send_right,position,MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&request);
            }
            //receiving data from neighbouring processes
            if(uflag)
            {
                position=0;
                MPI_Recv(recv_up,N*sizeof(double),MPI_PACKED,myrank-nc,myrank-nc,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_up,N*sizeof(double),&position,up+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(dflag)
            {
                position=0;
                MPI_Recv(recv_down,N*sizeof(double),MPI_PACKED,myrank+nc,myrank+nc,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_down,N*sizeof(double),&position,down+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(lflag)
            {
                position=0;
                MPI_Recv(recv_left,N*sizeof(double),MPI_PACKED,myrank-1,myrank-1,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_left,N*sizeof(double),&position,left+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(rflag)
            {
                position=0;
                MPI_Recv(recv_right,N*sizeof(double),MPI_PACKED,myrank+1,myrank+1,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_right,N*sizeof(double),&position,right+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            // MPI_Barrier(MPI_COMM_WORLD);
            //updating interior points using the 5-point stencil computation
            for(int i=1;i<N-1;i++)
            {
                for(int j=1;j<N-1;j++)
                {
                    new[i*N+j]=(old[(i-1)*N+j]+old[(i+1)*N+j]+old[i*N+j-1]+old[i*N+j+1]+old[i*N+j])/5;
                }
            }    
            //updating corner points using the boundary conditions
            new[0*N+0]= (up[0]+left[0]+old[0*N+1]+old[1*N+0]+old[0*N+0])/den[0*N+0];
            new[0*N+N-1]= (up[N-1]+right[0]+old[0*N+N-2]+old[1*N+N-1]+old[0*N+N-1])/den[0*N+N-1];
            new[(N-1)*N+0]= (down[0]+left[N-1]+old[(N-1)*N+1]+old[(N-2)*N+0]+old[(N-1)*N+0])/den[(N-1)*N+0];
            new[(N-1)*N+N-1]= (down[N-1]+right[N-1]+old[(N-2)*N+N-1]+old[(N-1)*N+N-2]+old[(N-1)*N+N-1])/den[(N-1)*N+N-1];        
            //updating boundary points (excluding corners)
            for(int i=1;i<N-1;i++)
            {
                new[0*N+i]= (old[0*N+i+1]+old[0*N+i-1]+old[1*N+i]+up[i]+old[0*N+i])/den[0*N+i];
                new[(N-1)*N+i]= (old[(N-1)*N+i+1]+old[(N-1)*N+i-1]+old[(N-2)*N+i]+down[i]+old[(N-1)*N+i])/den[(N-1)*N+i];
                new[i*N+0]= (old[(i+1)*N+0]+old[(i-1)*N+0]+old[i*N+1]+left[i]+old[i*N+0])/den[i*N+0];
                new[i*N+N-1]=(old[(i+1)*N+N-1]+old[(i-1)*N+N-1]+old[i*N+N-2]+right[i]+old[i*N+N-1])/den[i*N+N-1];
            }
            //updating old array with new array
            for(int i=0;i<N;i++)
            {
                for(int j=0;j<N;j++)
                {
                    old[i*N+j]=new[i*N+j];
                }
            }
        }
        eTime=MPI_Wtime();
    }
    //calculating the time taken by the code
    double time = eTime-sTime;
    //reducing the time taken by all processes to get the maximum time taken
    MPI_Reduce(&time,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(myrank==0)
    {
		printf("%lf\n",maxTime);    	
    }
    MPI_Finalize();
    return 0;
}