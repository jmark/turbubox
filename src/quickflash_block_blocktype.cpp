// C++ program file quickflash_block_blocktype.cpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Block type description.
*/


#include "quickflash_block_blocktype.hpp"
#include <vector>
#include <string>


namespace QuickFlash
{
namespace Block
{

// BlockType Labels

const char Leaf_BlockType_Label[] = "Leaf";
const char Branch_BlockType_Label[] = "Branch";

// Node ID processing

const unsigned int Leaf_ID = 1;


BlockType convert_block_type(const int block_type_id)
  {
  BlockType block_type = Leaf;

  if (block_type_id != static_cast<int>(Leaf_ID))
    block_type = Branch;

  return block_type;
  }


void convert_block_type(const std::vector<int> & block_type_id,
			std::vector<BlockType> & block_type)
  {
  const unsigned int num_blocks = block_type_id.size();

  block_type.resize(num_blocks);

  for (unsigned int i = 0; i < num_blocks; i++)
    block_type[i] = convert_block_type(block_type_id[i]);
  }


void convert_block_type(const std::vector<unsigned int> & block_type_id,
			std::vector<BlockType> & block_type)
  {
  const unsigned int num_blocks = block_type_id.size();

  block_type.resize(num_blocks);

  for (unsigned int i = 0; i < num_blocks; i++)
    block_type[i] = convert_block_type(static_cast<int>(block_type_id[i]));
  }


void get_block_type_label(const BlockType block_type, 
			  std::string & block_type_label)
  {
  switch (block_type)
    {
    case Leaf :
      block_type_label = Leaf_BlockType_Label;
      break;

    case Branch :
      block_type_label = Branch_BlockType_Label;
      break;
    }
  }


}  // End namespace Block
}  // End namespace QuickFlash
