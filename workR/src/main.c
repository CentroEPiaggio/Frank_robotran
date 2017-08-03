   /**
    *
    *   Universite catholique de Louvain
    *   CEREM : Centre for research in mechatronics
    *   http://www.robotran.be  
    *   Contact : info@robotran.be
    *
    *
    *   MBsysC main script template for simple model:
    *   -----------------------------------------------
    *    This template loads the data file *.mbs and execute:
    *      - the coordinate partitioning module
    *      - the direct dynamic module (time integration of
    *        equations of motion).
    *    It may be adapted and completed by the user.
    * 
    *    (c) Universite catholique de Louvain
    *
    * To turn this file as a C++ file, just change its extension to .cc (or .cpp).
    * If you plan to use some C++ files, it is usually better that the main is compiled as a C++ function.
    * Currently, most compilers do not require this, but it is a safer approach to port your code to other computers.
    */

#include <stdio.h>
#include "mbs_data.h"
#include "mbs_dirdyn.h"
#include "mbs_part.h"
#include "realtime.h"
#include "mbs_set.h"
#include "mbs_load_xml.h"
#include "cmake_config.h"

int main(int argc, char const *argv[])
{
    
    MbsData *mbs_data;
    
    MbsPart *mbs_part;
    MbsDirdyn *mbs_dirdyn;

    printf("Starting Frank_segway MBS project!\n");


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    /*                     LOADING                               *
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    
    printf("Loading the Frank_segway data file !\n");
    mbs_data = mbs_load(PROJECT_SOURCE_DIR"/../dataR/Frank_segway.mbs", BUILD_PATH);
    printf("*.mbs file loaded!\n");
    
    // Set the number of user constraints
    int N_usr_c = 5;
    mbs_set_nb_userc(mbs_data, N_usr_c);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    /*              COORDINATE PARTITIONING                      *
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    mbs_part = mbs_new_part(mbs_data);
    mbs_part->options->rowperm = 1;
    mbs_part->options->verbose = 1;
    mbs_run_part(mbs_part, mbs_data);
    
    mbs_delete_part(mbs_part);


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    /*                   DIRECT DYNANMICS                        *
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    mbs_dirdyn = mbs_new_dirdyn(mbs_data);

    // dirdyn options (see documentations for additional options)
    mbs_dirdyn->options->dt0 = 1e-3;
    mbs_dirdyn->options->tf  = 10.0;
    mbs_dirdyn->options->save2file = 1;
    mbs_dirdyn->options->respath  = PROJECT_SOURCE_DIR"/../resultsR";
    mbs_dirdyn->options->animpath = PROJECT_SOURCE_DIR"/../animationR";
    //mbs_dirdyn->options->realtime = 1;
    
    mbs_run_dirdyn(mbs_dirdyn, mbs_data);
    
    mbs_delete_dirdyn(mbs_dirdyn, mbs_data);
    
    
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    /*                   CLOSING OPERATIONS                      *
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    mbs_delete_data(mbs_data);
    
    return 0;
}

