// C++ program file quickflash_tree_block.cpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Root structure for set of QuickFlash tree objects.
*/


#include "quickflash_tree_block.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Tree
{

// Class TreeBlockNodeID

int TreeBlockNodeID::get_branch(ConstTreeIter & tree_iter) const
  {
  int ret_val = -1;  // Error state

  if (block_ptr != 0)
    {
    tree_iter.reset(*block_ptr);

    ret_val = tree_iter.branch(tree_index, nodeid);
    }

  return ret_val;
  }


int TreeBlockNodeID::get_node(ConstTreeIter & tree_iter) const
  {
  int ret_val = -1;  // Error state

  if (block_ptr != 0)
    {
    tree_iter.reset(*block_ptr);

    ret_val = tree_iter.node(tree_index, nodeid);
    }

  return ret_val;
  }


// Class TreeIter

int TreeIter::next() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.next() >= 0)
      ret_val = 0;  // OK
    else
      {
      treeindex++;

      if (treeindex < blockptr->get_num_cells())
	{
	(*blockptr)[treeindex].get_node_iter(nodeiter);

	if (nodeiter.valid_node())
	  ret_val = 0;
	}
      else
	nodeiter.set_invalid();
      }
    }

  return ret_val;
  }


int TreeIter::prev() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.prev() >= 0)
      ret_val = 0;  // OK
    else if (treeindex > 0)
      {
      treeindex--;

      (*blockptr)[treeindex].get_node_iter(nodeiter);

      if (nodeiter.goto_last() >= 0)
	ret_val = 0;
      }
    else
      nodeiter.set_invalid();
    }

  return ret_val;
  }


int TreeIter::next_branch() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.next_skipBranch() >= 0)
      ret_val = 0;  // OK
    else
      {
      treeindex++;

      if (treeindex < blockptr->get_num_cells())
	{
	(*blockptr)[treeindex].get_node_iter(nodeiter);

	if (nodeiter.valid_node())
	  ret_val = 0;
	}
      else
	nodeiter.set_invalid();
      }
    }

  return ret_val;
  }


int TreeIter::first() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    if (blockptr->get_num_cells() > 0)
      {
      treeindex = 0;

      (*blockptr)[0].get_node_iter(nodeiter);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

  return ret_val;
  }


int TreeIter::last() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    const unsigned int num_cells = blockptr->get_num_cells();

    if (num_cells > 0)
      {
      (*blockptr)[num_cells - 1].get_node_iter(nodeiter);

      if (nodeiter.goto_last() >= 0)
	ret_val = 0;
      }
    }

  return ret_val;
  }


int TreeIter::local_branch(const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    ret_val = nodeiter.goto_branch(node_id);

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int TreeIter::local_node(const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    ret_val = nodeiter.goto_node(node_id);

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int TreeIter::branch(const unsigned int tree_block_index, 
		     const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    if (tree_block_index < blockptr->get_num_cells())
      {
      treeindex = tree_block_index;

      (*blockptr)[tree_block_index].get_node_iter(nodeiter);

      nodeiter.goto_branch(node_id);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int TreeIter::node(const unsigned int tree_block_index, 
		   const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    if (tree_block_index < blockptr->get_num_cells())
      {
      treeindex = tree_block_index;

      (*blockptr)[tree_block_index].get_node_iter(nodeiter);

      nodeiter.goto_node(node_id);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


// Class ConstTreeIter

int ConstTreeIter::next() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.next() >= 0)
      ret_val = 0;  // OK
    else
      {
      treeindex++;

      if (treeindex < blockptr->get_num_cells())
	{
	(*blockptr)[treeindex].get_node_iter(nodeiter);

	if (nodeiter.valid_node())
	  ret_val = 0;
	}
      else
	nodeiter.set_invalid();
      }
    }

  return ret_val;
  }


int ConstTreeIter::prev() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.prev() >= 0)
      ret_val = 0;  // OK
    else if (treeindex > 0)
      {
      treeindex--;

      (*blockptr)[treeindex].get_node_iter(nodeiter);

      if (nodeiter.goto_last() >= 0)
	ret_val = 0;
      }
    else
      nodeiter.set_invalid();
    }

  return ret_val;
  }


int ConstTreeIter::next_branch() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    if (nodeiter.next_skipBranch() >= 0)
      ret_val = 0;  // OK
    else
      {
      treeindex++;

      if (treeindex < blockptr->get_num_cells())
	{
	(*blockptr)[treeindex].get_node_iter(nodeiter);

	if (nodeiter.valid_node())
	  ret_val = 0;
	}
      else
	nodeiter.set_invalid();
      }
    }

  return ret_val;
  }


int ConstTreeIter::first() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    if (blockptr->get_num_cells() > 0)
      {
      treeindex = 0;

      (*blockptr)[0].get_node_iter(nodeiter);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

  return ret_val;
  }


int ConstTreeIter::last() const
  {
  int ret_val = -1;

  if (blockptr != 0)
    {
    const unsigned int num_cells = blockptr->get_num_cells();

    if (num_cells > 0)
      {
      (*blockptr)[num_cells - 1].get_node_iter(nodeiter);

      if (nodeiter.goto_last() >= 0)
	ret_val = 0;
      }
    }

  return ret_val;
  }


int ConstTreeIter::local_branch(const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    ret_val = nodeiter.goto_branch(node_id);

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int ConstTreeIter::local_node(const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    ret_val = nodeiter.goto_node(node_id);

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int ConstTreeIter::branch(const unsigned int tree_block_index, 
			  const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    if (tree_block_index < blockptr->get_num_cells())
      {
      treeindex = tree_block_index;

      (*blockptr)[tree_block_index].get_node_iter(nodeiter);

      nodeiter.goto_branch(node_id);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


int ConstTreeIter::node(const unsigned int tree_block_index, 
			const NodeID & node_id) const
  {
  int ret_val = -1;  // Error state

  if (blockptr != 0)
    {
    if (tree_block_index < blockptr->get_num_cells())
      {
      treeindex = tree_block_index;

      (*blockptr)[tree_block_index].get_node_iter(nodeiter);

      nodeiter.goto_node(node_id);

      if (nodeiter.valid_node())
	ret_val = 0;
      }

    if (ret_val < 0)
      set_invalid();
    }

  return ret_val;
  }


}  // End namespace Tree
}  // End namespace QuickFlash
