/*
 * Copyright 2009, 2010 Jan Gasthaus (j.gasthaus@gatsby.ucl.ac.uk)
 * 
 * This file is part of libPLUMP.
 * 
 * libPLUMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * libPLUMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with libPLUMP.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "libplump/random.h"

#include <gsl/gsl_rng.h>

namespace gatsby { namespace libplump {

gsl_rng* global_rng = 0;

void init_rng() {
       const gsl_rng_type * T;
       gsl_rng_env_setup();
     
       T = gsl_rng_default;
       global_rng = gsl_rng_alloc (T);
}

/**
 * Free the global RNG.
 */
void free_rng() {
       gsl_rng_free (global_rng);
}

}}
