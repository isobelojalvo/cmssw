

#include "DataFormats/L1CaloTrigger/interface/L1PFObject.h"

using std::ostream;
using std::endl;
using std::hex;
using std::dec;

// default constructor
L1PFObject::L1PFObject() : m_data(0) { }



// destructor
L1PFObject::~L1PFObject() { }


// print to stream
ostream& operator << (ostream& os, const L1PFObject& clus) {
  os << "L1PFObject:";
  os << " Et=" << clus.et();
  return os;
}

