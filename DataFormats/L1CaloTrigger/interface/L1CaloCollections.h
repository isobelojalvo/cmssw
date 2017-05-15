
#ifndef L1CALOCOLLECTIONS_H
#define L1CALOCOLLECTIONS_H

#include <vector>

#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCluster.h"
#include "DataFormats/L1CaloTrigger/interface/L1PFObject.h"

typedef std::vector<L1CaloEmCand> L1CaloEmCollection;
typedef std::vector<L1CaloRegion> L1CaloRegionCollection;
typedef std::vector<L1CaloCluster> L1CaloClusterCollection;
typedef std::vector<L1PFObject> L1PFObjectCollection;


#endif
