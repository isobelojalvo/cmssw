/// ////////////////////////////////////////
/// Stacked Tracker Simulations          ///
/// ////////////////////////////////////////

#include "DataFormats/Common/interface/Wrapper.h"


/*********************/
/** L1 CALO TRIGGER **/
/*********************/

#include "DataFormats/CaloClusters/interface/L1CaloCluster.h"

namespace DataFormats_CaloClustersger {
  L1CaloCluster cluster;
  edm::Wrapper<L1CaloCluster> w_cluster;

  std::vector<L1CaloCluster> caloClusterCollection;
  edm::Wrapper<std::vector<L1CaloCluster> > w_caloClusterCollection;

};

namespace {
  namespace {

    L1CaloCluster cluster;
    edm::Wrapper<L1CaloCluster> w_cluster;
    
    std::vector<L1CaloCluster> caloClusterCollection;
    edm::Wrapper<std::vector<L1CaloCluster> > w_caloClusterCollection;
    
  }
}

