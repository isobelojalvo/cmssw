

#include "DataFormats/L1Trigger/interface/L1PFTau.h"

using std::ostream;
using std::endl;
using std::hex;
using std::dec;

typedef std::vector<L1PFTau> L1PFTauCollection;

// default constructor
L1PFTau::L1PFTau() : m_data(0) { }


// destructor
L1PFTau::~L1PFTau() { }


// print to stream
ostream& operator << (ostream& os, const L1PFTau& clus) {
  os << "L1PFTau:";
  os << " Et=" << clus.et();
  return os;
}

