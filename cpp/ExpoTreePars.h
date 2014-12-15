/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the ETH Zurich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __EXPOTREEPARS_H__
#define __EXPOTREEPARS_H__

#include <cstring>
#include <iostream>
#include <iomanip>

class ExpoTreePars {
  public:
    ExpoTreePars(int nchains, int npar, int max_shifts);
    ExpoTreePars(const ExpoTreePars&);
    virtual ~ExpoTreePars();

    inline double* t(int par, int chain = 0) { return _theta + chain*dim + par*(max_shifts+1); }
    inline const double* t(int par, int chain = 0) const { return _theta + chain*dim + par*(max_shifts+1); }

    inline double* s(int ival, int chain = 0) { return _shifts + chain*max_shifts + ival; }
    inline const double* s(int ival, int chain = 0) const { return _shifts + chain*max_shifts + ival; }

    inline int* st(int ival, int chain = 0) { return _shift_types + chain*max_shifts + ival; }
    inline const int* st(int ival, int chain = 0) const { return _shift_types + chain*max_shifts + ival; }

    inline int* ns(int chain = 0) { return _nshifts + chain; }
    inline const int* ns(int chain = 0) const { return _nshifts + chain; }

    inline double& theta(int par, int ival, int chain = 0) { return *(t(par,chain)+ival); }
    inline double theta(int par, int ival, int chain = 0) const { return *(t(par,chain)+ival); }

    inline double& shift(int ival, int chain = 0) { return *(s(ival,chain)); }
    inline double shift(int ival, int chain = 0) const { return *(s(ival,chain)); }

    inline int& shift_type(int ival, int chain = 0) { return *(st(ival,chain)); }
    inline int shift_type(int ival, int chain = 0) const { return *(st(ival,chain)); }

    inline int& nshifts(int chain = 0) { return *(ns(chain)); }
    inline int nshifts(int chain = 0) const { return *(ns(chain)); }

    inline int num_chains() const { return nchains; }

    friend std::ostream& operator<<(std::ostream& out, const ExpoTreePars& etp);

  protected:
    int nchains;
    int npar;
    int max_shifts;
    double* _theta;
    double* _shifts;
    int* _shift_types;
    int* _nshifts;

  private:
    int dim;
};

#endif
