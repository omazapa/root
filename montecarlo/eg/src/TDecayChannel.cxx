// @(#)root/eg:$Id$
// Author: P.Murat   15/02/2001

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//-----------------------------------------------------------------------------
//  Feb 16 2001 P.Murat: description of the decay channel
////////////////////////////////////////////////////////////////////////////////

#include "TDecayChannel.h"

ClassImp(TDecayChannel)

////////////////////////////////////////////////////////////////////////////////
///default constructor

TDecayChannel::TDecayChannel()
{
   fNumber            = 0;
   fMatrixElementCode = 0;
   fBranchingRatio    = 0;
}

////////////////////////////////////////////////////////////////////////////////
///constructor

TDecayChannel::TDecayChannel(Int_t    NumberD,
                             Int_t    MatrixElementType,
                             Double_t BRatio,
                             Int_t    NumberDaughters,
                             Int_t*   DaughterCode)
{
   fNumber            = NumberD;
   fMatrixElementCode = MatrixElementType;
   fBranchingRatio    = BRatio;
   fDaughters.Set(NumberDaughters,DaughterCode);
}

////////////////////////////////////////////////////////////////////////////////

TDecayChannel::~TDecayChannel() {
}

