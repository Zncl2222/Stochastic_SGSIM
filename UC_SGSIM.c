#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <io.h>

//global parameter

int start=0;
int end=150;
int mlen=150;
int len=10;
int C0=1;
double a=17.32;
int nR=5;
int seed=1516;
int neigh=0;
int loglevel=0;


int* model_x;
double u;
double*u_array,*dist_temp,*data_temp,*weights,*dist_temp2,*flatten_temp,*distcov_temp,*sampled,*z;
double**deter_temp,**fac_temp,** b_temp,**inverse,**b_temp2,**b,**datacov_temp,**z_array,**inv_temp,**pdist_temp;


FILE *output;
char fhead[]="R";
char ftail[]=".txt";
char path[]="./Realizations/";
char number1[10];


double variance(double* array)
{
    double mean=0;
    double var=0;
    for (int i=0;i<mlen;i++)
    {
        mean=mean+array[i];
    }

    mean=mean/mlen;

    for (int i=0;i<mlen;i++)
    {
        var=var+(pow(array[i]-mean,2));
    }
    var=var/mlen;

    return var;
}


int* randompath(int* rpath,int length)
{ 
  int rtemp;
  for(int i=length-1;i;i--)
  {
    int rindex= rand()%i;
    
    do
    {
      rtemp=rpath[rindex];
      
      rpath[rindex]=rpath[i];
      
      rpath[i]=rtemp;

    }while(0);

  }

  return rpath;
}


void LUdecomposition(double** a ,double* b, double* x , int n) 
{
   int i = 0, j = 0, k = 0;
   double l[10][10],u[10][10];
   double c[10];
 
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (j < i)
         l[j][i] = 0;
         else {
            l[j][i] = a[j][i];
            for (k = 0; k < i; k++) {
               l[j][i] = l[j][i] - l[j][k] * u[k][i];
            }
         }
      }
      for (j = 0; j < n; j++) {
         if (j < i)
         u[i][j] = 0;
         else if (j == i)
         u[i][j] = 1;
         else {
            u[i][j] = a[i][j] / l[i][i];
            for (k = 0; k < i; k++) {
               u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
            }
         }
      }
   }

    // Solve L(Ux)=b, assume Ux=c 
    c[0]=b[0]/l[0][0];
    
    for (i=1;i<n;i++)
    {
        c[i]=b[i];
        
        for (j=0;j<i;j++)
        {
            c[i]=c[i]-l[i][j]*c[j];
            
        }
        c[i]=c[i]/l[i][i];
    }

    x[n]=c[n];

    for (i=n-1;i>=0;i--)
    {
        x[i]=c[i];
        for (j=i+1;j<n;j++)
        {
            x[i]-=u[i][j]*x[j];
        }
    }


}

int printarray_pdist(double**x,int n_dim,char* c )
{   
    printf(c);
    for (int i=0;i<n_dim;i++)
    {
        for (int j=0;j<n_dim;j++)
        {
            printf("%lf  ",x[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    return 0;
}

void matrixform(double* x,double**matrix,int n_dim)
{
 
    for (int i=0;i<n_dim;i++)
    {
        for (int j=0;j<n_dim;j++)
        {
            matrix[i][j]=x[n_dim*i+j];
        }

    }

}
double** sort2d(double** x, int n_dim)
{   
    double tempp;
    double tempp2;
    double tempp3;

      for(int i = 0; i < n_dim; i++) {

       for(int j = i; j < n_dim; j++) {

           if( x[j][2] < x[i][2] ) {

                tempp = x[j][0];
                tempp2= x[j][1];
                tempp3= x[j][2];

                x[j][0] = x[i][0];
                x[j][1] = x[i][1];
                x[j][2] = x[i][2];

                x[i][0] = tempp;
                x[i][1] = tempp2;
                x[i][2] = tempp3;
           }
       }
   }
   return x;
}


int* arange(int x)
{
    int* space;

    space=(int*)malloc(x*sizeof(int));

    for (int i=0;i<x;i++)
    {
        space[i]=i;
    }
    return space;
}


float random_normal()
{   
    
    float x= rand()/(float)RAND_MAX;
    float y= rand()/(float)RAND_MAX;
    float z=sqrt(-2*log(x))*cos(2*M_PI*y);
    return z;
}

void pdist(double* x, double** c,int n_dim)
{


    for (int i=0;i<n_dim;i++)
    {
        for (int j=0;j<n_dim;j++)
        {
            c[i][j]=fabs(x[j]-x[i]);
        }
    }

}

void Cov_model(double *x, double* cov,int n_dim)
{

    for(int i=0;i<n_dim*n_dim;i++)
    {   
            cov[i]=C0-(C0*(1-exp(-3*pow(x[i],2)/pow(a,2))));  
    }


}

void Cov_model2(double **x,double* cov,int n_dim)
{

    for(int i=0;i<n_dim;i++)
    {
        for (int j=0;j<n_dim;j++)
        {
            cov[n_dim*i+j]=C0-(C0*(1-exp(-3*pow(x[i][j],2)/pow(a,2))));
        }
        
    }


}



void flatten(double** x, int n_dim)
{
    double *flat;

    flat=(double*)malloc(n_dim*n_dim*sizeof(double));


    for (int i=0;i<n_dim;i++)
    {
        for (int j=0;j<n_dim;j++)
        {
            flat[n_dim*i+j]=x[i][j];
        }
    }

}


int Print_Log1(double* a,double* b,double* c, int curr, int n_dim)
{
   
    printf("\nU=%d, current=%d\n\n",(int)u,curr);
    printf("Location = ");
    for (int j=0;j<n_dim;j++)
    {
        printf("%lf   ",a[j]);
    }
    printf("\nData =    ");
    for (int j=0;j<n_dim;j++)
    {
        printf("%lf   ",b[j]);
    }
    printf("\nDistdiff = ");
    for (int j=0;j<n_dim;j++)
    {
        printf("%lf   ",c[j]);
    }
    printf("\n");

    return 0;
}

int Print_Log2(double** a,double** b,double* c,double* d,double z_temp, double fix_temp, int n_dim)
{
    printf("\nPdist = \n");
    for (int j=0;j<n_dim;j++)
    {   
        for (int k=0;k<n_dim;k++)
        {
            printf("%lf ",a[j][k]);
        }
        printf("\n");
    }
    printf("\n");

    printf("DataCov = \n");
    for (int j=0;j<n_dim;j++)
    {   
        for (int k=0;k<n_dim;k++)
        {
            printf("%lf ",b[j][k]);
        }
        printf("\n");
        
    }
    printf("\n");

    printf("Distcov = \n");
    for (int j=0;j<n_dim;j++)
    {
        printf("%lf   ",c[j]);
    }
    printf("\n");

    printf("Weights = \n");
    for (int j=0;j<n_dim;j++)
    {
        printf("%f ",d[j]);
    }

    printf("\n\nFix= %f ",fix_temp);
    printf("\nEstimate=%lf\n",z_temp);
    
    
    return 0;
}

int main(void)
{   
    printf("One dimensional UnConditional Sequential Gaussian SIMulation (UCSGSIM)\n");

    printf("Please input the parameters for simulation\n");

    printf("Length of model\n");

    scanf(" %d",&end);

    printf("Number of realizations\n");

    scanf(" %d",&nR);

    printf("Correlation range for kriging\n");

    scanf(" %f",&a);

    printf("Initial random seed number\n");

    scanf(" %d",&seed);

    printf("Log level (0,1,2,3)?  (defaule is 0)\n");

    scanf(" %d", &loglevel);

    mlen=end+1;

    const char *Folder="./Realizations/";


    _mkdir(Folder);



    char fullname[strlen(path)+strlen(ftail)+strlen(fhead)+strlen(number1)];
    
    srand(seed);

    double z0=0;
    double svar=0;
    float fix=0;
    
    clock_t start, end;
    start=clock();

    /* Memory allocation*/

    pdist_temp=(double**)malloc(6*sizeof(double));
    dist_temp=(double*)malloc(6*sizeof(double));
    dist_temp2=(double*)malloc(6*sizeof(double));
    data_temp=(double*)malloc(6*sizeof(double));
    z_array=(double**)malloc(mlen*sizeof(double*));
    inv_temp=(double**)malloc(6*sizeof(double*));
    weights=(double*)malloc(6*sizeof(double));
    sampled=(double*)malloc(mlen*sizeof(double));
    flatten_temp=(double*)malloc(6*6*sizeof(double));
    datacov_temp=(double**)malloc(6*sizeof(double));
    distcov_temp=(double*)malloc(6*sizeof(double));

    deter_temp=(double**)malloc(6*sizeof(double*));
    fac_temp=(double**)malloc(6*sizeof(double*));
    b_temp=(double**)malloc(6*sizeof(double*));
    inverse=(double**)malloc(6*sizeof(double*));
    b_temp2=(double**)malloc(6*sizeof(double*));
    b=(double**)malloc(6*sizeof(double*));

    for (int i=0;i<mlen;i++)
    {
        z_array[i]=(double*)malloc(3*sizeof(double));
    }

    for (int i=0;i<6;i++)
    {
        pdist_temp[i]=(double*)malloc(6*sizeof(double));
        inv_temp[i]=(double*)malloc(6*sizeof(double)); 
        datacov_temp[i]=(double*)malloc(6*sizeof(double));
        deter_temp[i]=(double*)malloc(6*sizeof(double));
        fac_temp[i]=(double*)malloc(6*sizeof(double));
        inverse[i]=(double*)malloc(6*sizeof(double));
        b_temp2[i]=(double*)malloc(6*sizeof(double));
        b_temp[i]=(double*)malloc(6*sizeof(double));
        b[i]=(double*)malloc(6*sizeof(double));
    }


    z=(double*)malloc(mlen*sizeof(double));
        
    u_array=(double*)malloc(mlen*sizeof(double));

    model_x=arange(mlen);

        
    for (int n=0;n<nR;n++)
    {   
         

    printf("\n\nRealizations: %d\n\n",n);
    for (int i=0;i<mlen;i++)
    {
        z[i]=0;
        u_array[i]=-1;
        sampled[i]=0;
    }

    

    model_x=randompath(model_x,mlen);


    int neigh=0;
    int currentlen=0;
    int judge=0;

    // Simple Kriging
    
    for (int i=0;i<mlen;i++)
    {   
        u=model_x[i];
        
        if (neigh==0)
        {
            z[(int)u]=random_normal();
            sampled[i]=u;
        }

        else
        { 
          
            for (int j=0; j<currentlen;j++)
            {   
                u_array[j]=fabs(sampled[j]-u);
                
            }

            
            int close=0;

            for (int j=0;j<currentlen;j++)
            {
                if(u_array[j]<a*1.732){close++;}
            }
            
            if (close==0)
            {
                z[(int)u]=random_normal();
                sampled[i]=u;
                
            }

            else
            {

                for (int j=0;j<currentlen;j++)
                {
                    z_array[j][0]=sampled[j];
                    z_array[j][1]=z[(int)(sampled[j])]; 
                    z_array[j][2]=u_array[j];
                }

                if (neigh>=2)
                {
                    z_array=sort2d(z_array,currentlen);
                }

                
                for (int j=0;j<neigh;j++)
                {
                    dist_temp[j]=z_array[j][0];
                    data_temp[j]=z_array[j][1];
                    dist_temp2[j]=z_array[j][2];
                }
   
                // data cov matrix
                pdist(dist_temp,pdist_temp,neigh);
                
                Cov_model2(pdist_temp,flatten_temp,neigh);

                matrixform(flatten_temp,datacov_temp,neigh);
                
                // dist cov matrix
                Cov_model(dist_temp2,distcov_temp,neigh);
  
                if (neigh>=1)
                {   
                    LUdecomposition(datacov_temp,distcov_temp,weights,neigh);
                }

                z0=0;
                svar=0;
                fix=0;

                for (int j=0;j<neigh;j++)
                {   
                    z0=z0+data_temp[j]*weights[j];
                    svar=svar+distcov_temp[j]*weights[j];
                }

                svar=1-svar;
                if(svar<0){svar=0;}
                fix=random_normal()*pow(svar,0.5);

                z[(int)u]=z0+fix;
                seed++;

                if (loglevel==1)
                {
                    Print_Log1(dist_temp,data_temp,dist_temp2,currentlen,neigh);
                }
                else if (loglevel==2)
                {
                    Print_Log1(dist_temp,data_temp,dist_temp2,currentlen,neigh);
                    Print_Log2(pdist_temp,datacov_temp,distcov_temp,weights,z0,fix,neigh);
                }
                
                if (isfinite(z[(int)u])==0)
                {
                    judge++;
                }

            }
            

        }
    
        if (neigh<5){neigh+=1;}
        sampled[i]=u;
        currentlen++;
        seed++;

    }   
        
        sprintf(number1,"%d", n);
        memset(fullname,'\0',strlen(path)+strlen(ftail)+strlen(fhead)+strlen(number1));
        strcat(fullname,path);
        strcat(fullname,fhead);
        strcat(fullname,number1);
        strcat(fullname,ftail);   
        
        if(0.35<variance(z)<10 && judge==0)
        {
        
            output= fopen (fullname, "w");
            
            for(int g=0;g<mlen;g++)
            {
                fprintf(output,"%d\t%f\n",g,z[g]);
            }

            fclose (output); 
        }

        else{n--;}
        seed++;
    }
       

    end=clock();
    double totaltime=(double)(end-start)*1e-3;
    printf("\n1D case finised\n");
    printf("Time:%.4fs\n",totaltime);
    

    system("PAUSE");
    return 0;
}