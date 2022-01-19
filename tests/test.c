#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <dirent.h>
#include <sys/types.h>
#include <errno.h>
#include "lina.h"

#define check assert

//Print the matrix A with size m by n
static void pmatrix(FILE *fp, double *A, int m, int n);


typedef struct dot_test{

	double *A;
	double *B;
	double *C;
	int m;
	int n;
	int l;

}dot_test;

typedef struct add_test{

	double *A;
	double *B;
	double *C;
	int m;
	int n;

}add_test;

typedef struct scale_test{

	double *A;
	double *B;
	double s;
	int m;
	int n;

}scale_test;

typedef struct transpose_test{

	double *A;
	double *B;
	int m;
	int n;

}transpose_test;

#define PATH "./tests/"

int main()
{	
	//Defining pointers to the test structures
	add_test *add_tests;
	dot_test *dot_tests;
	scale_test *scale_tests;
	transpose_test *transpose_tests;

	//Number of tests for each lina functions
	int n_dot_tests, n_add_tests, n_scale_tests, n_transpose_tests;
	
	//Opening dir stream
	DIR *dir = opendir(PATH);
	struct dirent *ep;
	check(dir != NULL);

	//Loading all the tests from files
	while (ep = readdir(dir))
	{
		if(ep->d_type != DT_DIR || !strcmp(ep->d_name, ".") || !strcmp(ep->d_name, ".."))
			continue;
		
		if(!strcmp(ep->d_name,"add"))
		{

			//Implement loading lina_add test matrices
			char sub_path[256] = PATH;
			strcat(sub_path,ep->d_name);
			strcat(sub_path,"/");
			
			DIR *sub_dir = opendir(sub_path);
			check(sub_dir != NULL);

			struct dirent *sub_ep;

			unsigned int count = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;
				count += 1;
			}

			closedir(sub_dir);

			sub_dir = opendir(sub_path);
			check(sub_dir != NULL);
			
			add_tests = malloc(sizeof(add_test)*count);
			n_add_tests = count;
			int i = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;

				char file_pos[256];
				strcat(file_pos,sub_path);
				strcat(file_pos,sub_ep->d_name);

				FILE *fp;
				fp = fopen(file_pos,"r");
				check(fp != NULL);

				int m,n;
				char *error;

				add_tests[i].A = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(add_tests[i].A  != NULL);

				add_tests[i].B = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(add_tests[i].B  != NULL);

				add_tests[i].C = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(add_tests[i].C  != NULL);

				add_tests[i].m = m;
				add_tests[i].n = n;

				i += 1;
				fclose(fp);

			}

			closedir(sub_dir);

		}
		else if (!strcmp(ep->d_name,"dot"))
		{
			//Implement loading lina_dot test matrices
			char sub_path[256] = PATH;
			strcat(sub_path,ep->d_name);
			strcat(sub_path,"/");
			
			DIR *sub_dir = opendir(sub_path);
			check(sub_dir != NULL);

			struct dirent *sub_ep;

			unsigned int count = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;
				count += 1;
			}

			closedir(sub_dir);

			sub_dir = opendir(sub_path);
			check(sub_dir != NULL);
			dot_tests = malloc(sizeof(dot_test)*count);
			n_dot_tests = count;

			int i = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;

				char file_pos[256];
				strcpy(file_pos,sub_path);
				strcat(file_pos,sub_ep->d_name);

				FILE *fp;
				fp = fopen(file_pos,"r");
				check(fp != NULL);

				int m,n,l;
				char *error;

				dot_tests[i].A = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(dot_tests[i].A  != NULL);

				dot_tests[i].B = lina_loadMatrixFromStream(fp,&l,&n,&error);
				check(dot_tests[i].B  != NULL);

				dot_tests[i].C = lina_loadMatrixFromStream(fp,&l,&m,&error);
				check(dot_tests[i].C  != NULL);

				dot_tests[i].m = m;
				dot_tests[i].n = n;
				dot_tests[i].l = l;

				i += 1;
				fclose(fp);

			}

			closedir(sub_dir);


		}
		else if (!strcmp(ep->d_name,"scale"))
		{
			//Implement loading lina_scale test matrices
			char sub_path[256] = PATH;
			strcat(sub_path,ep->d_name);
			strcat(sub_path,"/");
			
			DIR *sub_dir = opendir(sub_path);
			check(sub_dir != NULL);

			struct dirent *sub_ep;

			unsigned int count = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;
				count += 1;
			}

			closedir(sub_dir);

			sub_dir = opendir(sub_path);
			check(sub_dir != NULL);

			scale_tests = malloc(sizeof(scale_test)*count);
			
			n_scale_tests = count;

			int i = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;

				char file_pos[256];
				strcpy(file_pos,sub_path);

				strcat(file_pos,sub_ep->d_name);

				FILE *fp;
				fp = fopen(file_pos,"r");
				check(fp != NULL);

				int m,n;
				int useless1,useless2;
				double *scale;
				char *error;

				scale_tests[i].A = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(scale_tests[i].A  != NULL);

				scale_tests[i].B = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(scale_tests[i].B  != NULL);

				scale = lina_loadMatrixFromStream(fp,&useless1,&useless2,&error);
				check(scale != NULL);

				scale_tests[i].m = m;
				scale_tests[i].n = n;
				scale_tests[i].s = scale[0];
				free(scale);


				i += 1;
				fclose(fp);

			}

			closedir(sub_dir);


		}
		else if (!strcmp(ep->d_name,"transpose"))
		{
			//Implement loading lina_transpose test matrices
			char sub_path[256] = PATH;
			strcat(sub_path,ep->d_name);
			strcat(sub_path,"/");
			
			DIR *sub_dir = opendir(sub_path);
			check(sub_dir != NULL);

			struct dirent *sub_ep;

			unsigned int count = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;
				count += 1;
			}

			closedir(sub_dir);

			sub_dir = opendir(sub_path);
			check(sub_dir != NULL);
			
			transpose_tests = malloc(sizeof(transpose_test)*count);
			n_transpose_tests = count;

			int i = 0;

			while (sub_ep = readdir(sub_dir))
			{
				if(sub_ep->d_type == DT_DIR)
					continue;

				char file_pos[256];
				strcpy(file_pos,sub_path);
				strcat(file_pos,sub_ep->d_name);

				FILE *fp;
				fp = fopen(file_pos,"r");
				check(fp != NULL);

				int m,n;
				char *error;
				
				transpose_tests[i].A = lina_loadMatrixFromStream(fp,&n,&m,&error);
				check(transpose_tests[i].A  != NULL);

				transpose_tests[i].B = lina_loadMatrixFromStream(fp,&m,&n,&error);
				check(transpose_tests[i].B  != NULL);

				transpose_tests[i].m = m;
				transpose_tests[i].n = n;

				i += 1;
				fclose(fp);

			}

			closedir(sub_dir);


		}

	}
	
	closedir(dir);

	//Starting the lina_add tests
	{
		int passed_tests = 0;
		fprintf(stdout,"Starting tests on lina_add():\n");
		for(int i=0;i<n_add_tests;i++){
			
			double *C = (double*) malloc(sizeof(*C)*add_tests[i].m * add_tests[i].n);
			check(C != NULL);

			lina_add(add_tests[i].A, add_tests[i].B, C,add_tests[i].m,add_tests[i].n);
			
			if( !memcmp(add_tests[i].C, C, sizeof(*C)*add_tests[i].m * add_tests[i].n) )
				passed_tests += 1;
			else{
				
				fprintf(stderr,"----------------------------------------------------\n");
				fprintf(stderr,"Test on lina_add() failed on the following matrices:\n");
				pmatrix(stderr,add_tests[i].A, add_tests[i].m, add_tests[i].n);
				fprintf(stderr,"+\n");
				pmatrix(stderr,add_tests[i].B, add_tests[i].m, add_tests[i].n);
				fprintf(stderr,"lina_add() gives following output:\n");
				pmatrix(stderr,C, add_tests[i].m, add_tests[i].n);
				fprintf(stderr,"instead of:\n");
				pmatrix(stderr,add_tests[i].C, add_tests[i].m, add_tests[i].n);
				fprintf(stderr,"----------------------------------------------------\n");


			}
			free(C);

		}
		if(n_add_tests != 0)
			fprintf(stdout, "Test on lina_add() finished: %d out of %d tests were succesfull\n",passed_tests,n_add_tests);
		else
			fprintf(stdout, "There are no tests for lina_add() function.\n");
	}

	//Starting the lina_dot tests
	{
		int passed_tests = 0;
		fprintf(stdout,"\nStarting tests on lina_dot():\n");
		for(int i=0;i<n_dot_tests;i++){
			
			double *C = (double*) malloc(sizeof(*C)*dot_tests[i].m * dot_tests[i].l);
			check(C != NULL);

			lina_dot(dot_tests[i].A, dot_tests[i].B, C, dot_tests[i].m, dot_tests[i].n, dot_tests[i].l);
			
			if( !memcmp(dot_tests[i].C, C, sizeof(*C)*dot_tests[i].m * dot_tests[i].l) )
				passed_tests += 1;
			else{
				
				fprintf(stderr,"----------------------------------------------------\n");
				fprintf(stderr,"Test on lina_dot() failed on the following matrices:\n");
				pmatrix(stderr,dot_tests[i].A, dot_tests[i].m, dot_tests[i].n);
				fprintf(stderr,"*\n");
				pmatrix(stderr,dot_tests[i].B, dot_tests[i].n, dot_tests[i].l);
				fprintf(stderr,"lina_dot() gives following output:\n");
				pmatrix(stderr,C, dot_tests[i].m, dot_tests[i].l);
				fprintf(stderr,"instead of:\n");
				pmatrix(stderr,dot_tests[i].C, dot_tests[i].m, dot_tests[i].l);
				fprintf(stderr,"----------------------------------------------------\n");


			}
			free(C);

		}
		if(n_dot_tests != 0)
			fprintf(stdout, "Test on lina_dot() finished: %d out of %d tests were succesfull\n",passed_tests,n_dot_tests);
		else
			fprintf(stdout, "There are no tests for lina_dot() function.\n");
	}

	//Starting the lina_transpose tests
	{
		int passed_tests = 0;
		fprintf(stdout,"\nStarting tests on lina_transpose():\n");
		for(int i=0;i<n_transpose_tests;i++){
			
			double *C = (double*) malloc(sizeof(*C)*transpose_tests[i].m * transpose_tests[i].n);
			check(C != NULL);

			lina_transpose(transpose_tests[i].A, C, transpose_tests[i].m, transpose_tests[i].n);
			
			if( !memcmp(transpose_tests[i].B, C, sizeof(*C)*transpose_tests[i].m * transpose_tests[i].n) )
				passed_tests += 1;
			else{
				
				fprintf(stderr,"----------------------------------------------------\n");
				fprintf(stderr,"Test on lina_transpose() failed on the following matrices:\n");
				pmatrix(stderr,transpose_tests[i].A, transpose_tests[i].m, transpose_tests[i].n);
				fprintf(stderr,"lina_transpose() gives following output:\n");
				pmatrix(stderr,C, transpose_tests[i].n, transpose_tests[i].m);
				fprintf(stderr,"instead of:\n");
				pmatrix(stderr,transpose_tests[i].B, transpose_tests[i].n, transpose_tests[i].m);
				fprintf(stderr,"----------------------------------------------------\n");


			}
			free(C);

		}
		if(n_transpose_tests != 0)
			fprintf(stdout, "Test on lina_transpose() finished: %d out of %d tests were succesfull\n",passed_tests,n_transpose_tests);
		else
			fprintf(stdout, "There are no tests for lina_transpose() function.\n");
	}

	//Starting the lina_scale tests
	{
		int passed_tests = 0;
		fprintf(stdout,"\nStarting tests on lina_scale():\n");
		for(int i=0;i<n_scale_tests;i++){
			
			double *C = (double*) malloc(sizeof(*C)*scale_tests[i].m * scale_tests[i].n);
			check(C != NULL);

			lina_scale(scale_tests[i].A, C, scale_tests[i].s, scale_tests[i].m, scale_tests[i].n);
			
			if( !memcmp(scale_tests[i].B, C, sizeof(*C)*scale_tests[i].m * scale_tests[i].n) )
				passed_tests += 1;
			else{
				
				fprintf(stderr,"----------------------------------------------------\n");
				fprintf(stderr,"Test on lina_scale() failed on the following matrices:\n");
				pmatrix(stderr,scale_tests[i].A, scale_tests[i].m, scale_tests[i].n);
				fprintf(stderr,"lina_scale() gives following output:\n");
				pmatrix(stderr, C, scale_tests[i].m, scale_tests[i].n);
				fprintf(stderr,"instead of:\n");
				pmatrix(stderr,scale_tests[i].B, scale_tests[i].m, scale_tests[i].n);
				fprintf(stderr,"----------------------------------------------------\n");


			}
			free(C);

		}
		if(n_scale_tests != 0)
			fprintf(stdout, "Test on lina_scale() finished: %d out of %d tests were succesfull\n",passed_tests,n_scale_tests);
		else
			fprintf(stdout, "There are no tests for lina_scale() function.\n");
	}

	//Freeing the memory of all the heap variables
	{
		//Freeing add_tests memory
		for(int i=0;i< n_add_tests;i++)
		{
			free(add_tests[i].A);
			free(add_tests[i].B);
			free(add_tests[i].C);
		}

		//Freeing dot_tests memory
		for(int i=0;i< n_dot_tests;i++)
		{
			free(dot_tests[i].A);
			free(dot_tests[i].B);
			free(dot_tests[i].C);
		}

		//Freeing scale_tests memory
		for(int i=0;i< n_scale_tests;i++)
		{
			free(scale_tests[i].A);
			free(scale_tests[i].B);
		}

		//Freeing transpose_tests memory
		for(int i=0;i< n_transpose_tests;i++)
		{
			free(transpose_tests[i].A);
			free(transpose_tests[i].B);
		}
	}

	return 0;
}

static void pmatrix(FILE *fp, double *A,int m,int n){

	for(int i = 0; i<m; i++){
		fprintf(fp, "    | ");
		for(int j = 0; j< n; j++)
			fprintf(fp, "%g ",A[i*n + j]);
		
		fprintf(fp, "|\n");

	}

	fprintf(fp, "\n");

}