/** ---------------------------
  * Robotran - MBsysC
  * 
  * Template file for equilibrium module
  * 
  * This files enable the user to call custom code just 
  * after loading the project. 
  * The call is done inside the mbs_load function.
  * 
  * (c) Universite catholique de Louvain
  *     
  */

#include "math.h"

#include "mbs_data.h"

#include "mbs_set.h"

/*! \brief user own initialization functions
 *
 * \param[in,out] mbs_data data structure of the model
 *
 */
void user_load_post(MbsData *mbs_data)
{
	
    // Set the number of user constraints
    int N_usr_c = 5;
    mbs_set_nb_userc(mbs_data, N_usr_c);	
}
