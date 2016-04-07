/*!
 *  \author Peter Staar
 *  \author Andrei Plamada
 */
#ifndef PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_INCLUDE_MATERIAL_HAMILTONIANS_H
#define PHYS_LIBRARY_PARAMETERS_MODELS_MATERIAL_HAMILTONIANS_INCLUDE_MATERIAL_HAMILTONIANS_H
enum    material_name {NiO_symmetric,NiO_unsymmetric, CuO2, CuO2_1band, SrVO3};
using material_name_type=material_name ;

#include "material_lattice_template.h"
#include "material_interaction_template.h"

#include "template_specialization_NiO/material_lattice_NiO.h"
#include "template_specialization_CuO2/material_lattice_CuO2.h"
#include "template_specialization_CuO2_1band/material_lattice_CuO2_1band.h"
#include "template_specialization_SrVO3/material_lattice_SrVO3.h"

#endif