
#include <vector>
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace DataFormats_L1CaloTrigger {
  struct dictionary {
    L1CaloEmCollection em;
    L1CaloRegionCollection rgn;
    L1CaloClusterCollection cluster;
    L1PFObjectCollection l1pf;
    
    edm::Wrapper<L1CaloEmCollection> w_em;
    edm::Wrapper<L1CaloRegionCollection> w_rgn;
    edm::Wrapper<L1CaloClusterCollection> w_cluster;
    edm::Wrapper<L1PFObjectCollection> w_l1pf;
  };
}
