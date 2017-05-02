

#include "DataFormats/L1CaloTrigger/interface/L1CaloCluster.h"

using std::ostream;
using std::endl;
using std::hex;
using std::dec;

// default constructor
L1CaloCluster::L1CaloCluster() : m_data(0) { }



// destructor
L1CaloCluster::~L1CaloCluster() { }


// print to stream
ostream& operator << (ostream& os, const L1CaloCluster& clus) {
  os << "L1CaloCluster:";
  os << " Et=" << clus.et();
  return os;
}

