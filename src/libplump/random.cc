/*
 * Copyright 2008-2016 Jan Gasthaus
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
