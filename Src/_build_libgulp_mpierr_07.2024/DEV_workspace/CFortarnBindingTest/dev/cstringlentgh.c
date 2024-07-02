#include<stdio.h>
#include<string.h>

int main()
{
        int rank,size;
        char iopath[512];

        memset(iopath,' ',sizeof(iopath));
        //sprintf(iopath,"/work/e05/e05/wkjee/Software/gulp-6.1.2/Src/Custom/path_test");
        //sprintf(iopath,"/work/e05/e05/wkjee/Software/klmc_demo/test1");
        sprintf(iopath,"/mnt/lustre/a2fs-work2/work/e05/e05/wkjee/Software/master-worker-C-Fortran-Python/mpi/C/master_worker/archer2_klmcgulp/build/dev_multi_2/A8");
	printf("length> %d\n",strlen(iopath));
        printf("iopath> %s\n",iopath);

        return 0;
}

