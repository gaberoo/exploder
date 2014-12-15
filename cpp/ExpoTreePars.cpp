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

#include "ExpoTreePars.h"

ExpoTreePars::ExpoTreePars(int nc, int np, int ms) 
  : nchains(nc), npar(np), max_shifts(ms)
{
  dim = np*(ms+1);
  _theta = new double[nchains*npar*(max_shifts+1)];
  _shifts = new double[nchains*max_shifts];
  _shift_types = new int[nchains*max_shifts];
  _nshifts = new int[nchains];
}

ExpoTreePars::ExpoTreePars(const ExpoTreePars& etp)
  : nchains(etp.nchains), npar(etp.npar), max_shifts(etp.max_shifts),
    dim(etp.dim)
{
  _theta = new double[nchains*npar*(max_shifts+1)];
  _shifts = new double[nchains*max_shifts];
  _shift_types = new int[nchains*max_shifts];
  _nshifts = new int[nchains];

  memcpy(_theta,etp._theta,nchains*npar*(max_shifts+1)*sizeof(double));
  memcpy(_shifts,etp._shifts,nchains*max_shifts*sizeof(double));
  memcpy(_shift_types,etp._shift_types,nchains*max_shifts*sizeof(double));
  memcpy(_nshifts,etp._nshifts,nchains*sizeof(int));
}

ExpoTreePars::~ExpoTreePars() 
{
  delete[] _theta;
  delete[] _shifts;
  delete[] _shift_types;
  delete[] _nshifts;
}

std::ostream& operator<<(std::ostream& out, const ExpoTreePars& etp) {
  out << "Number of chains:    " << etp.num_chains() << std::endl;
  out << "Maximum shifts:      " << etp.max_shifts << std::endl;
  out << "Number of paramters: " << etp.npar << std::endl;

  for (int chain(0); chain < etp.num_chains(); ++chain) {
    for (int ival(0); ival <= etp.nshifts(chain) && ival <= etp.max_shifts; ++ival) {
      for (int par(0); par < etp.npar; ++par) {
        out << std::setw(12) << etp.theta(par,ival,chain) << " ";
      }
      if (ival > 0) out << " | " << etp.shift(ival-1,chain);
      out << std::endl;
    }
    out << std::endl;
  }

  return out;
}


