//-*-C++-*-

#include "../include_linalg.h"

void cuda_organizer(int status);

void cuda_thread_synchronize();




//#include "test_reading.h"

#include "test_GEMM.h"
#include "test_GEMD.h"
#include "test_BENNET.h"
#include "test_TRSV.h"

#include "test_copy_from.h"

//#include "test_GEINV.h"
//#include "test_operator.h"



int main(int argc,char *argv[])
{
    //testing_cuda_init();
    cuda_organizer(0);

    std::cout << "\n\n\t hello world !\n\n";


    test_GEMM();

    test_copy_from();

    //test_GEMD();
    
    //test_BENNET();

    test_TRSV();

/*
    test_reading();

    test_GEINV();

    test_operator();
*/
    //void testing_cuda_finalize();

    cuda_organizer(1);

    return 0;
}
