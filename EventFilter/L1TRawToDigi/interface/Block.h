#ifndef Block_h
#define Block_h

namespace l1t {
   typedef uint32_t BlockId;

   class BlockHeader {
      public:
         BlockHeader(unsigned int id, unsigned int size) : data_(((id & ID_mask) << ID_shift) | ((size & size_mask) << size_shift)) {};
         BlockHeader(const uint32_t *data) : data_(data[0]) {};

         inline unsigned int getID() const { return (data_ >> ID_shift) & ID_mask; };
         inline unsigned int getSize() const { return (data_ >> size_shift) & size_mask; };

         inline uint32_t raw() const { return data_; };

      private:
         static const unsigned int ID_shift = 24;
         static const unsigned int ID_mask = 0xff;
         static const unsigned int size_shift = 16;
         static const unsigned int size_mask = 0xff;

         uint32_t data_;
   };

   class Block {
      public:
         Block(const uint32_t *data) : header_(0) {};
         Block(unsigned int id, const std::vector<uint32_t>& payload) : header_(id, payload.size()), payload_(payload) {};

         inline unsigned int getSize() const { return payload_.size() + 1; };

         BlockHeader header() const { return header_; };
         std::vector<uint32_t> payload() const { return payload_; };

      private:
         BlockHeader header_;
         std::vector<uint32_t> payload_;
   };

   typedef std::vector<Block> Blocks;
}

#endif
