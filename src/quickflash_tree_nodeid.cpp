// C++ program file quickflash_tree_nodeid.cpp

/*
  By Nathan C. Hearn
     October 12, 2006

  Node ID vector for use with identifying tree nodes.

  NOTE: Since level zero (the root node) has no entries, levels and indexes
  are off by one, i.e., the information for level i is stored at index (i-1).
  Think of level zero as the basement, level one is at index zero,
  level two is at index 1, level three is at index 2, and so on.
*/


#include "quickflash_tree_nodeid.hpp"


namespace QuickFlash
{
namespace Tree
{

bool NodeID::branch_member(const NodeID & child_id) const
  {
  bool is_member = true;

  const unsigned int child_level = child_id.get_level();
  const unsigned int my_level = get_level();


  if (child_level < my_level)
    is_member = false;  // child_id further up tree trunk
  else
    {
    // Verify id's up to my level

    unsigned int i = 0;

    while ((i < my_level) && is_member)
      {
      if (child_id[i] != id[i])
	is_member = false;

      i++;
      }
    }

  return is_member;
  }


void NodeID::invert_level_bit(const unsigned int level,
			      const unsigned int bit_index)
  {
  // ERROR CHECKING???

  const unsigned int level_id = id[level];

  const unsigned int bit_rep = 1 << bit_index;

  if ((level_id & bit_rep) > 0)
    id[level] = level_id - bit_rep;
  else
    id[level] = level_id + bit_rep;
  }


void NodeID::invert_bits(const unsigned int bit_index)
  {
  const unsigned int num_levels = id.size();

  const unsigned int bit_rep = 1 << bit_index;

  for (unsigned int level = 0; level < num_levels; level++)
    {
    const unsigned int level_id = id[level];

    if ((level_id & bit_rep) > 0)
      id[level] = level_id - bit_rep;
    else
      id[level] = level_id + bit_rep;
    }
  }


}  // End namespace Tree
}  // End namespace QuickFlash
