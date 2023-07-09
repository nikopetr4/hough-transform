#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>



void read_rawimage(fname, length, width, image)
    char *fname; unsigned long length; unsigned long width; unsigned char image[length][width];
{
        short i;
        FILE *file;

        file=fopen(fname,"r");
        for (i=0; i<length; i++)
                        fread(image[i], 1, width, file);
        fclose(file);
}

void write_rawimage(fname, length, width, image)
    char *fname; unsigned long length; unsigned long width; unsigned char image[length][width];
{
        short i;
        FILE *file;

        file=fopen(fname,"w");
        for (i=0; i<length; i++)
                       fwrite(image[i], 1, width, file);


        fclose (file);
}
int find_edge(int x, int y, int rho, int theta)
{
	float PI=3.14159265359;
    float rad = PI / 180;
    int a;
    a = ((x*sin(rad*(theta))) - (y*cos(rad*(theta)))+rho);
    
    if (a < 1 && a > -1)
        {return 1;}
    else
        {return 0;}
}
//φιλτρο
void copy_in_2_out_img (length, width, inimg, outimg)
    unsigned long length, width;
    unsigned char inimg[length][width], outimg[length][width];

	{ 	

	
	int i, j, x, y, sobelx, sobely, sgmval;
    int sgm_max = 0;

    double start,end,Ts;
    double start2,end2,Tp;

	start =omp_get_wtime();

    int max_rho ;
	max_rho = sqrt(length*length +width*width);
    int bins[2*max_rho][180];
    int rho,theta;
	float PI=3.14159265359;
    float rad = PI / 180;

	int a;
	
	for(i=0;i<length;i++)
		for(j=0;j<width;j++){
			outimg[i][j]=inimg[i][j];
			inimg[i][j]=0;
		}


 

		for(i=0;i<2*max_rho;i++)
			for(j=0;j<180;j++){
				bins[i][j]=0;
			}


    int K,L;    
    for(y = 1; y < length-1; y++)
    {
        for(x = 1; x < width-1; x++)
        {
            /* dE/dx */
            sobelx = (-(outimg[y+1][x-1] + 2*outimg[y][x-1] + outimg[y-1][x-1]) + outimg[y+1][x+1] + 2*outimg[y][x+1] + outimg[y-1][x+1]);
                
            /* dE/dy */
           sobely = (-(outimg[y+1][x-1] + 2*outimg[y+1][x] + outimg[y+1][x+1]) + outimg[y-1][x-1] + 2*outimg[y-1][x] + outimg[y-1][x+1]);           
            /* SGM */
    
            K=sobelx*sobelx;
            L=sobely*sobely;
        
            sgmval = sqrt(K+L);

            inimg[y][x] = sgmval;   

            if(sgmval > sgm_max)
                { sgm_max = sgmval;}
       }
    }
	//printf("max is %d \n",sgm_max);
       
        for (i=0; i<length; i++)
                for (j=0; j<width; j++){
                   		
                        if (outimg[i][j]>0) {
                                for (theta = 0; theta < 180; theta++) {
                                		rho = j * cos(theta*rad) + i * sin(theta*rad); 
                						if(inimg[i][j] > 0.1*sgm_max && inimg[i][j] <0.2*sgm_max){
											//printf("SGM is %d \n",inimg[i][j]);
                							bins[max_rho+rho][theta]++;
               								}
                                   		}
								}																			                                                                    
                }	
        

	//normalize bins
	int maxbins=0;
    for(i = 0; i < 2*max_rho; i++)
    {
        for(j = 0; j < 180; j++)
        {
            if(bins[i][j] > maxbins)
            {
            	maxbins=bins[i][j];
            }
        }
    } 
	int l;
    /* Threshold bins */
    for(i = 0; i < 2*max_rho; i++)
    {
        for(j = 0; j < 180; j++)
        
        {	
            if(bins[i][j]> 0.78*maxbins)
            {
            	l=i-max_rho;
            	//printf("[%d,%d]=%d\n",l,j,bins[i][j]);

            for(y=0;y<length;y++)
            	for(x=0;x<width;x++){
					inimg[y][x]=outimg[y][x];
					   		
             		if(find_edge(x, y, l-(length/2), j))
					 {
            			outimg[y][x]=255;
					 }
				}
               
            }
        }
    }
	end=omp_get_wtime();
	Ts=end-start;
	printf("Time for serial took %f sec\n",Ts);
	





    
	int P;
	for(P=2;P<65;P*=2){
	start2=omp_get_wtime();
	
	
	#pragma omp parallel for private(i,j) shared(outimg) num_threads(P) collapse(2)
	for(i=0;i<length;i++)
	for(j=0;j<width;j++){
		outimg[i][j]=inimg[i][j];
		inimg[i][j]=0;
	}


 
	#pragma omp parallel for private(i,j) shared(bins) num_threads(P) collapse(2)
	for(i=0;i<2*max_rho;i++)
		for(j=0;j<180;j++){
			bins[i][j]=0;
		}


    int K,L;    
    #pragma omp parallel for private(y,x) num_threads(P) collapse(2)
    for(y = 1; y < length-1; y++)
    {
        for(x = 1; x < width-1; x++)
        {
            /* dE/dx */
            sobelx = (-(outimg[y+1][x-1] + 2*outimg[y][x-1] + outimg[y-1][x-1]) + outimg[y+1][x+1] + 2*outimg[y][x+1] + outimg[y-1][x+1]);
                
            /* dE/dy */
           sobely = (-(outimg[y+1][x-1] + 2*outimg[y+1][x] + outimg[y+1][x+1]) + outimg[y-1][x-1] + 2*outimg[y-1][x] + outimg[y-1][x+1]);           
            /* SGM */
    
            K=sobelx*sobelx;
            L=sobely*sobely;
        
            sgmval = sqrt(K+L);

            inimg[y][x] = sgmval;   

            if(sgmval > sgm_max)
                { sgm_max = sgmval;}
       }
    }
//printf("max is %d \n",sgm_max);
       #pragma omp parallel for private(i,j,theta) num_threads(P) collapse(2)
        for (i=0; i<length; i++)
                for (j=0; j<width; j++){
                   		
                        if (outimg[i][j]>0) {
                                for (theta = 0; theta < 180; theta++) {
                                		rho = j * cos(theta*rad) + i * sin(theta*rad); 
                						if(inimg[i][j] > 0.1*sgm_max && inimg[i][j] < 0.2*sgm_max ){
											//printf("SGM is %d \n",inimg[i][j]);
                							bins[max_rho+rho][theta]++;
               								}
                                   		}
								}																			                                                                    
                }	
        

//normalize bins
	int maxbins=0;
	#pragma omp parallel for private(i,j) num_threads(P) collapse(2)
    for(i = 0; i < 2*max_rho; i++)
    {
        for(j = 0; j < 180; j++)
        {
            if(bins[i][j] > maxbins)
            {
            	maxbins=bins[i][j];
            }
        }
    } 
	int l;
    /* Threshold bins */
    #pragma omp parallel for private(i,j,x,y) shared(outimg) num_threads(P) collapse(2)
    for(i = 0; i < 2*max_rho; i++)
    {
        for(j = 0; j < 180; j++)
        
        {	
            if(bins[i][j]> 0.78*maxbins)
            {
            	l=i-max_rho;
            	//printf("[%d,%d]=%d\n",l,j,bins[i][j]);

            for(y=0;y<length;y++)
            	for(x=0;x<width;x++){
					inimg[y][x]=outimg[y][x];            		
             		if(find_edge(x, y, l-(length/2), j))
					 {
            			outimg[y][x]=255;
					 }
				}
               
            }
        }
    }
	

	end2=omp_get_wtime();
	Tp=end2-start2;
	printf("Time for %d threads took %f sec\n",P,Tp);
}

	return;	      
}

	





// τελος φιλτρου 
void main(int argc, char *argv[]) {
    char infname[50], outfname[50];
    unsigned char **inimg, **outimg;
    unsigned long height, width, i, j;

    if (argc<4) {printf("usage is: %s inimg height width [outimg]\n", argv[0]);exit(-1);}
    strcpy(infname,argv[1]);
    height=(unsigned long)atoi(argv[2]);
    width=(unsigned long)atoi(argv[3]);
    strcpy(outfname,argv[4]);
    
    inimg = (unsigned char **) malloc (height*width);
    read_rawimage(infname, height, width, inimg);



    outimg = (unsigned char **) malloc (height*width);
      
  

  

    copy_in_2_out_img (height, width, inimg, outimg);
   




  
    write_rawimage(outfname, height, width, outimg);
    free(outimg);
		
    free(inimg);
}

	
