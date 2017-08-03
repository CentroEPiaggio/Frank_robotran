/* ---------------------------
 * Robotran - MBsysC
 * 
 * Template file for direct dynamics module
 * 
 * This files enable the user to call custom at
 * specific places in the time simulation. It is a template
 * file that can be edited by the user.
 * 
 * (c) Universite catholique de Louvain
 *     
 */

#define _USE_MATH_DEFINES  // to use M_PI in Windows
#include "math.h"

#include "mbs_data.h"
#include "mbs_dirdyn_struct.h"
//
#include "set_output.h"
#include "user_all_id.h"
#include "mbs_project_interface.h"

#include "user_IO.h"
//
/*! \brief user own initialization functions
 *
 * \param[in,out] mbs_data data structure of the model
 * \param[in,out] mbs_dd general structure of the direct dynamic module (for advance users)
 *
 * For beginners, it is advised to only use the MbsData structure.
 * The field MbsDirdyn is provided for more advance users.
 */
void user_dirdyn_init(MbsData *mbs_data, MbsDirdyn *mbs_dd)
{
   /* MbsSensor psens[1];                        // Creation of a pointer to a sensor struct.
    allocate_sensor(psens, mbs_data->njoint);  // Allocate the Jacobian at the correct dimension
    init_sensor(psens, mbs_data->njoint);      // Initialize all value to zero
    mbs_sensor(psens, mbs_data, Sensor_right_wheel_id); // Compute the sensor
    //printf("%f", psens->P[1]);                  // Print the X position (for example)
    //free_sensor(psens);                        // Free the memory (always better) */
}

/*! \brief user own loop functions
 *
 * \param[in,out] mbs_data data structure of the model
 * \param[in,out] mbs_dd general structure of the direct dynamic module (for advance users)
 *
 * For beginners, it is advised to only use the MbsData structure.
 * The field MbsDirdyn is provided for more advance users.
 */
void user_dirdyn_loop(MbsData *mbs_data, MbsDirdyn *mbs_dd)
{
    int id = Sensor_right_wheel_id;
    
    // retrieve the pointer to the sensor structure defined in mbs_aux
    MbsSensor *PtrSensor = mbs_dd->mbs_aux->psens;

    // compute the sensor (position, velocity...)
    mbs_sensor(PtrSensor, mbs_data, id);

    // save the x, y, z positions
    set_output(PtrSensor->P[1], "X_Position");
    set_output(PtrSensor->P[2], "Y_Position");
    set_output(PtrSensor->P[3], "Z_Position");
    
    //*mbs_data->user_IO->Sens_R_wheel = *PtrSensor->P;
    
//     MbsSensor psens[1];                        // Creation of a pointer to a sensor struct.
//     allocate_sensor(psens, mbs_data->njoint);  // Allocate the Jacobian at the correct dimension
//     init_sensor(psens, mbs_data->njoint);      // Initialize all value to zero
//     mbs_sensor(psens, mbs_data, Sensor_right_wheel_id); // Compute the sensor
//     //printf("%f", psens->P[1]);                  // Print the X position (for example)
//     
//     mbs_data->user_IO->Sens_R_wheel = psens->P;
//     free_sensor(psens);
//     //Sens = mbs_data->user_IO->Sensor_right_wheel;
}

/*! \brief user own finishing functions
 *
 * \param[in,out] mbs_data data structure of the model
 * \param[in,out] mbs_dd general structure of the direct dynamic module (for advance users)
 *
 * For beginners, it is advised to only use the MbsData structure.
 * The field MbsDirdyn is provided for more advance users.
 */
void user_dirdyn_finish(MbsData *mbs_data, MbsDirdyn *mbs_dd)
{

}
