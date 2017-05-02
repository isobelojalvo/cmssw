#ifndef L1CALOCLUSTER_H
#define L1CALOCLUSTER_H

#include <ostream>

class L1CaloCluster
{
public:

  /// default constructor
  L1CaloCluster();


  
  /// destructor
  ~L1CaloCluster();



  // get/set methods for the data

  /// reset the data content (not position id!)
  void reset() { m_data = 0; }

  /// get raw data
  uint16_t raw() const { return m_data; }

  /// get Et
  unsigned et() const { return (0); }

  /// set data
  void setRawData(uint32_t data) { m_data = data; }


  /// equality operator, including rank, feature bits, and position
  //int operator==(const L1CaloCluster& c) const { return ((m_data==c.raw() && m_id==c.id()) || (this->empty() && c.empty())); }

  /// inequality operator
  //int operator!=(const L1CaloCluster& c) const { return !(*this == c); }

  /// is there any information in the candidate
  bool empty() const { return (m_data == 0); }

  /// print to stream
  friend std::ostream& operator << (std::ostream& os, const L1CaloCluster& reg);

 private:

  uint16_t m_data;

};


#endif /*L1CALOCLUSTER_H*/
