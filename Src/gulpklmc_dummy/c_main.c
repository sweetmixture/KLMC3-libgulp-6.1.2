#include<mpi.h>
#include<stdio.h>
#include<string.h>

extern void gulpklmc( MPI_Comm*, char*, int*, int* );
//subroutine gulpklmc( MPI_comm_klmc, klmc_task_iopath, klmc_task_id, klmc_worker_id ) bind (C,name="gulpklmc")

int main()
{
	int rank,size;
	char iopath[512];

	memset(iopath,' ',sizeof(iopath));
        //sprintf(iopath,"/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/gulpklmc/samplerun");
	//printf("C iopath> %s\n",iopath);

	int task_id = 1;
	int worker_id = 1;

	MPI_Comm base_comm;

	MPI_Init(NULL,NULL);
	base_comm = MPI_COMM_WORLD;

	MPI_Comm_size(base_comm,&size);
	MPI_Comm_rank(base_comm,&rank);

	// 1. CALL GULP
	memset(iopath,' ',sizeof(iopath));
        sprintf(iopath,"/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/gulpklmc/sample1");
	printf("C iopath> %s\n",iopath);
        if( rank == 0 ) printf("C main> before calling gulpklmc\n");
	printf("===========================================================================\n");
	gulpklmc( &base_comm, iopath, &task_id, &worker_id );
        if( rank == 0 ) printf("C main> after calling gulpklmc\n");

	// 2. CALL GULP
	memset(iopath,' ',sizeof(iopath));
        sprintf(iopath,"/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/gulpklmc/sample2");
	printf("C iopath> %s\n",iopath);
        if( rank == 0 ) printf("C main> before calling gulpklmc\n");
	printf("===========================================================================\n");
	gulpklmc( &base_comm, iopath, &task_id, &worker_id );
        if( rank == 0 ) printf("C main> after calling gulpklmc\n");

	MPI_Finalize();
	return 0;
}
