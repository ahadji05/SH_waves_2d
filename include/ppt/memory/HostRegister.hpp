/**
 * @file
 *
 * @author  Andreas Hadjigeorgiou, The Cyprus Institute,
 *          Personal-site: https://ahadji05.github.io,
 *          E-mail: a.hadjigeorgiou@cyi.ac.cy
 *
 * @copyright 2022 CaSToRC (The Cyprus Institute), Delphi Consortium (TU Delft)
 *
 * @version 1.0
 *
 * @section LICENCE
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PPT_HOST_REGISTER_HPP
#define PPT_HOST_REGISTER_HPP

#ifdef PPT_ENABLE_CUDA_BACKEND
#include "cuda_runtime.h"
#else
#include "sys/mman.h"
#endif

namespace ppt {

class HostRegister {

    public:

        static int mem_lock( void *m, size_t nbytes ){
            #ifdef PPT_ENABLE_CUDA_BACKEND
                cudaHostRegister(m, nbytes, cudaHostRegisterPortable);
            #else
                mlock(m, nbytes);
            #endif
            return 0;
        }

        static int mem_unlock( void *m, [[maybe_unused]] size_t nbytes ){
            #ifdef PPT_ENABLE_CUDA_BACKEND
                cudaHostUnregister(m);
            #else
                munlock(m, nbytes);
            #endif
            return 0;
        }
};

} // namespace ppt

#endif