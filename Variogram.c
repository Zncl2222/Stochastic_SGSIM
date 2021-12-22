#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <io.h>

//global parameter

double mlen,nR,seed;
double a;
double gmean,gstd;
int hs=35;
int steps=1;
char vario_fhead[]="Variogram";
char fhead[]="Realizations";
char ftail[]=".txt";
char number1[15];
char n1[]="000";
char n2[]="00";
char n3[]="0";

double**pdist_temp;
double* variogram_array;
FILE* vario_output;

int parameter_read(char* Path)
{   
    double parameter[20];
    int count=0;

    FILE* para=fopen(Path,"r");
    
    if(para==NULL)
    {
        printf("ERROR,no data\n");
        return 0;
    }
    
   const unsigned MAX_LENGTH = 1000;
   char buffer[MAX_LENGTH];

   while (fgets(buffer, MAX_LENGTH, para))
   {
       parameter[count]=atof(buffer);
       count++;
   }

    mlen=parameter[0];
    a=parameter[1];
    nR=parameter[3];
    seed=parameter[4];
    gmean=parameter[6];
    gstd=parameter[7];

    return 0;
}

double* data_read(char* Path)
{   
    int count=0;
    double** data_temp;
    double* data;

    data=(double*)malloc(mlen*sizeof(double));
    data_temp=(double**)malloc(mlen*sizeof(double));

    for (int i=0;i<mlen;i++)
    {
        data_temp[i]=(double*)malloc(2*sizeof(double));
    }

    FILE* para=fopen(Path,"r");
    
    if(para==NULL)
    {
        printf("ERROR,no data\n");
        return 0;
    }
    
    for (int i=0;i<mlen;i++)
    {   
        for (int j=0;j<2;j++)
        {
            fscanf(para,"%lf",&data_temp[i][j]);
        }
            
    }

    fclose(para);

    for (int i=0;i<mlen;i++)
    {
        data[i]=data_temp[i][1];
    }

    return data;
}

double* d_arange(int x)
{
    double* space;

    space=(double*)malloc(x*sizeof(double));

    for (int i=0;i<x;i++)
    {
        space[i]=i;
    }
    return space;
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

void variogram(double* array, double* v ,int hs, int steps)
{ 
    double Z_temp;
    double* temp;
    double count;

    pdist_temp=(double**)malloc(mlen*sizeof(double));
    

    for (int i=0;i<mlen;i++)
    {
        pdist_temp[i]=(double*)malloc(mlen*sizeof(double));
    }

    temp=d_arange(mlen);

    pdist(temp,pdist_temp,mlen);

    for (int i=0;i<hs;i+=steps)
    {   
        Z_temp=0;
        count=0;

        for (int j=0;j<mlen;j++)
        {
            for (int k=j+1;k<mlen;k++)
            {   
                
                if (pdist_temp[j][k]>=i-steps && pdist_temp[j][k]<=i+steps)
                {
                    Z_temp=Z_temp+pow((array[j]-array[k]),2);
                    count+=1;
                    //printf("Test %d\n",i);
                }
            }
        }

        if (Z_temp>=1e-6)
        {   
            v[i]=Z_temp/(2*count);
        }
        
    }

};

int main()
{
    char Data_path[100];
    char Parameter_path[strlen(Data_path)+50];
    char Folder_v[100];
    char Variogram[]="\\Variogram\\";
    char fullname[strlen(Data_path)+strlen(ftail)+strlen(fhead)+strlen(number1)];

    double* temp_data;
    
    printf("Input the Realizations path\n");
    scanf(" %s",&Data_path);

    memset(Parameter_path,'\0',strlen(Data_path)+50);
    strcat(Parameter_path,Data_path);
    strcat(Parameter_path,"\\ParameterSettings.txt");

    printf("Data = %s\n",Data_path);
    printf("Parameter = %s\n",Parameter_path);

    parameter_read(Parameter_path);

    printf("Model lengtg = %f\n",mlen);
    printf("Effective range = %f\n",a);
    printf("Number of realizations = %f\n",nR);
    printf("Random seed = %f\n",seed);
    printf("Mean vaule = %f\n",gmean);
    printf("Standard deviation = %f\n",gstd);

    printf("Input the len of lag distance\n");
    scanf(" %d",&hs);
    printf("Input the steps of lag distance\n");
    scanf(" %d",&steps);

    strcat(Folder_v,Data_path);
    strcat(Folder_v,Variogram);
    strcat(Data_path,"\\");
    _mkdir(Folder_v);

    variogram_array=(double*)malloc(hs*sizeof(double));

    for (int i=0;i<nR;i++)
    {   
        sprintf(number1,"%d", i);
        memset(fullname,'\0',strlen(Data_path)+strlen(ftail)+strlen(fhead)+5+strlen(number1));
        strcat(fullname,Data_path);
        strcat(fullname,fhead);
        if (i<10)
        {   
            strcat(fullname,n1);
        }
        else if (i<100 && i>=10)
        {   
            strcat(fullname,n2);
        }
        else if (i<1000 && i>=100)
        {   
            strcat(fullname,n3);
        }

        strcat(fullname,number1);
        strcat(fullname,ftail);   
        
        temp_data=data_read(fullname);

        for (int j=0;j<mlen;j++)
        {
            temp_data[j]=(temp_data[j]-gmean)/gstd;
        }

        variogram(temp_data,variogram_array,hs,steps);

        memset(fullname,'\0',strlen(Data_path)+strlen(ftail)+strlen(fhead)+12+strlen(number1));
        strcat(fullname,Folder_v);
        strcat(fullname,vario_fhead);

        if (i<10)
        {   
            strcat(fullname,n1);
        }
        else if (i<100 && i>=10)
        {   
            strcat(fullname,n2);
        }
        else if (i<1000 && i>=100)
        {   
            strcat(fullname,n3);
        }
                    
        strcat(fullname,number1);
        strcat(fullname,ftail);

        vario_output=fopen(fullname,"w");

        for (int j=0;j<hs;j+=steps)
        {
            fprintf(vario_output,"%d\t%.10lf\n",j,variogram_array[j]);
        }

        fclose(vario_output);

        printf("Variogram %d \n ",i);

    }

    
    system("PAUSE");

    return 0;
}