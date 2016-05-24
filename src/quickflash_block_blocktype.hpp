// C++ header file quickflash_block_blocktype.hpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Block type description.
*/


#ifndef QUICKFLASH_BLOCK_BLOCKTYPE_HPP
#define QUICKFLASH_BLOCK_BLOCKTYPE_HPP


#include <vector>
#include <string>


namespace QuickFlash
{
namespace Block
{
// Node ID processing

enum BlockType { Leaf, Branch };

extern const char Leaf_BlockType_Label[];
extern const char Branch_BlockType_Label[];

extern const unsigned int Leaf_ID;

BlockType convert_block_type(const int block_type_id);

void convert_block_type(const std::vector<int> & block_type_id,
			std::vector<BlockType> & block_type);

void convert_block_type(const std::vector<unsigned int> & block_type_id,
			std::vector<BlockType> & block_type);

void get_block_type_label(const BlockType block_type, 
			  std::string & block_type_label);


}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_BLOCKTYPE_HPP
