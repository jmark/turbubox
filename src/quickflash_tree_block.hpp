// C++ header file quickflash_tree_block.hpp

/*
  By Nathan C. Hearn
     October 25, 2006

  Root structure and iterators for set of QuickFlash tree objects.
*/


#ifndef QUICKFLASH_TREE_BLOCK_HPP
#define QUICKFLASH_TREE_BLOCK_HPP


#include "quickflash_tree_node.hpp"
#include "quickflash_tree_nodedata.hpp"
#include "quickflash_block_domain.hpp"
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_tree_nodeid.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Tree
{

typedef Node<NodeData> TreeNode;
typedef NodeIter<NodeData> TreeNodeIter;
typedef ConstNodeIter<NodeData> ConstTreeNodeIter;

typedef Block::Domain<TreeNode> TreeBlock;


// Forward declarations

class TreeIter;
class ConstTreeIter;


// Class TreeBlockNodeID

class TreeBlockNodeID
  {
  public :

    TreeBlockNodeID() :
      block_ptr(0), tree_index(0), nodeid()
      { }

    TreeBlockNodeID(const TreeBlock & tree_block,
		    const unsigned int tree_block_index,
		    const NodeID & node_id) :
      block_ptr(0), tree_index(0), nodeid()
      { reset(tree_block, tree_block_index, node_id); }

    TreeBlockNodeID(const TreeBlock * tree_block,
		    const unsigned int tree_block_index,
		    const NodeID & node_id) :
      block_ptr(0), tree_index(0), nodeid()
      { reset(tree_block, tree_block_index, node_id); }

    TreeBlockNodeID(const TreeBlockNodeID & source) :
      block_ptr(0), tree_index(0), nodeid()
      { reset(source); }

    ~TreeBlockNodeID() { }

    TreeBlockNodeID & operator=(const TreeBlockNodeID & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      block_ptr = 0;
      tree_index = 0;

      nodeid.reset();
      }

    void reset(const TreeBlock & tree_block,
	       const unsigned int tree_block_index, const NodeID & node_id)
      {
      block_ptr = &tree_block;
      tree_index = tree_block_index;

      nodeid = node_id;
      }
      
    void reset(const TreeBlock * tree_block,
	       const unsigned int tree_block_index, const NodeID & node_id)
      {
      block_ptr = tree_block;
      tree_index = tree_block_index;

      nodeid = node_id;
      }

    void reset(const TreeBlockNodeID & source)
      {
      if (&source != this)
	{
	block_ptr = source.block_ptr;
	tree_index = source.tree_index;

	nodeid = source.nodeid;
	}
      }

    void invert_id(const unsigned int axis)
      {
      if (block_ptr == 0)
	throw Except("No tree block specified", __FILE__, __LINE__);

      tree_index = block_ptr->invert_cell_index(axis, tree_index);

      nodeid.invert_bits(axis);
      }

    unsigned int get_tree_block_index() const { return tree_index; }

    const NodeID & get_node_id() const { return nodeid; }

    unsigned int get_level() const { return nodeid.get_level(); }

    unsigned int get_refinement_level() const
      {
      // Flash refinement levels are one more than the node level

      const unsigned int refine_level = nodeid.get_level() + 1;

      return refine_level;
      }

    int get_branch(ConstTreeIter & tree_iter) const;

    int get_node(ConstTreeIter & tree_iter) const;

  private :

    const TreeBlock * block_ptr;

    unsigned int tree_index;

    NodeID nodeid;
  };


// Class TreeIter

class TreeIter
  {
  friend class ConstTreeIter;

  public:
    TreeIter() : blockptr(0), treeindex(0), nodeiter() { reset(); }

    TreeIter(TreeBlock & block) :
      blockptr(0), treeindex(0), nodeiter()
      { reset(block); }

    TreeIter(const TreeIter & source) :
      blockptr(0), treeindex(0), nodeiter()
      { reset(source); }

    ~TreeIter() { }

    TreeIter & operator=(const TreeIter & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      blockptr = 0;
      treeindex = 0;
      nodeiter.reset();
      }

    void reset(TreeBlock & block)
      {
      blockptr = &block;
      treeindex = 0;

      if (block.get_num_cells() > 0)
	block[0].get_node_iter(nodeiter);
      else
	nodeiter.reset();
      }

    void reset(const TreeIter & source)
      {
      if (&source != this)
	{
	blockptr = source.blockptr;
	treeindex = source.treeindex;
	nodeiter = source.nodeiter;
	}
      }

    NodeData & get_node_data() { return nodeiter.get_data(); }
    const NodeData & get_node_data() const { return nodeiter.get_data(); }

    const Block::BlockInfo & get_block_info() const 
      { return nodeiter.get_data().get_block_info(); }

    unsigned int get_block_index() const 
      { return nodeiter.get_data().get_block_index(); }

    bool valid_node() const
      {
      bool valid = false;

      if (blockptr != 0)
	valid = nodeiter.valid_node();

      return valid;
      }

    void set_invalid() const { nodeiter.set_invalid(); }

    bool is_leaf() const
      {
      bool leaf_state = false;

      if (blockptr != 0)
	if (!(nodeiter.children_present()))
	  leaf_state = true;

      return leaf_state;
      }

    int next() const;
    int prev() const;

    int next_branch() const;

    int first() const;
    int last() const;

    int local_branch(const NodeID & node_id) const;
    int local_node(const NodeID & node_id) const;

    int branch(const unsigned int tree_block_index, const NodeID & node_id)
      const;
    int node(const unsigned int tree_block_index, const NodeID & node_id)
      const;

    int branch(const TreeBlockNodeID & id) const
      { return branch(id.get_tree_block_index(), id.get_node_id()); }

    int node(const TreeBlockNodeID & id) const
      { return node(id.get_tree_block_index(), id.get_node_id()); }

    const TreeBlock & get_tree_block() const { return *blockptr; }

    unsigned int get_tree_block_index() const { return treeindex; }

    const NodeID & get_node_id() const { return nodeiter.get_nodeID(); }

    unsigned int get_level() const { return get_node_id().get_level(); }

    unsigned int get_refinement_level() const
      { 
      // Flash refinement levels are one more than the node level

      const unsigned int refine_level = get_node_id().get_level() + 1;

      return refine_level;
      }

    void get_tree_block_nodeid(TreeBlockNodeID & id) const
      { id.reset(blockptr, treeindex, nodeiter.get_nodeID()); }

    TreeBlockNodeID get_tree_block_nodeid() const
      { return TreeBlockNodeID(blockptr, treeindex, nodeiter.get_nodeID()); }

    const TreeNodeIter & get_node_iter() const { return nodeiter; }

  protected:
    mutable TreeBlock * blockptr;

    mutable unsigned int treeindex;

    mutable TreeNodeIter nodeiter;
  };


// Class ConstTreeIter

class ConstTreeIter
  {
  public:
    ConstTreeIter() : blockptr(0), treeindex(0), nodeiter() { reset(); }

    ConstTreeIter(const TreeBlock & block) :
      blockptr(0), treeindex(0), nodeiter()
      { reset(block); }

    ConstTreeIter(const TreeIter & source) :
      blockptr(0), treeindex(0), nodeiter()
      { reset(source); }

    ConstTreeIter(const ConstTreeIter & source) :
      blockptr(0), treeindex(0), nodeiter()
      { reset(source); }

    ~ConstTreeIter() { }

    ConstTreeIter & operator=(const TreeIter & source)
      {
      reset(source);
      return *this;
      }

    ConstTreeIter & operator=(const ConstTreeIter & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      blockptr = 0;
      treeindex = 0;
      nodeiter.reset();
      }

    void reset(const TreeBlock & block)
      {
      blockptr = &block;
      treeindex = 0;

      if (block.get_num_cells() > 0)
	block[0].get_node_iter(nodeiter);
      else
	nodeiter.reset();
      }

    void reset(const TreeIter & source)
      {
      blockptr = source.blockptr;
      treeindex = source.treeindex;
      nodeiter = source.nodeiter;
      }

    void reset(const ConstTreeIter & source)
      {
      if (&source != this)
	{
	blockptr = source.blockptr;
	treeindex = source.treeindex;
	nodeiter = source.nodeiter;
	}
      }

    const NodeData & get_node_data() const { return nodeiter.get_data(); }

    const Block::BlockInfo & get_block_info() const 
      { return nodeiter.get_data().get_block_info(); }

    unsigned int get_block_index() const 
      { return nodeiter.get_data().get_block_index(); }

    bool valid_node() const
      {
      bool valid = false;

      if (blockptr != 0)
	valid = nodeiter.valid_node();

      return valid;
      }

    void set_invalid() const { nodeiter.set_invalid(); }

    bool is_leaf() const
      {
      bool leaf_state = false;

      if (blockptr != 0)
	if (!(nodeiter.children_present()))
	  leaf_state = true;

      return leaf_state;
      }

    int next() const;
    int prev() const;

    int next_branch() const;

    int first() const;
    int last() const;

    int local_branch(const NodeID & node_id) const;
    int local_node(const NodeID & node_id) const;

    int branch(const unsigned int tree_block_index, const NodeID & node_id)
      const;
    int node(const unsigned int tree_block_index, const NodeID & node_id)
      const;

    int branch(const TreeBlockNodeID & id) const
      { return branch(id.get_tree_block_index(), id.get_node_id()); }

    int node(const TreeBlockNodeID & id) const
      { return node(id.get_tree_block_index(), id.get_node_id()); }

    const TreeBlock & get_tree_block() const { return *blockptr; }

    unsigned int get_tree_block_index() const { return treeindex; }

    const NodeID & get_node_id() const { return nodeiter.get_nodeID(); }

    unsigned int get_level() const { return get_node_id().get_level(); }

    unsigned int get_refinement_level() const
      { 
      // Flash refinement levels are one more than the node level

      const unsigned int refine_level = get_node_id().get_level() + 1;

      return refine_level;
      }

    void get_tree_block_nodeid(TreeBlockNodeID & id) const
      { id.reset(blockptr, treeindex, nodeiter.get_nodeID()); }

    TreeBlockNodeID get_tree_block_nodeid() const
      { return TreeBlockNodeID(blockptr, treeindex, nodeiter.get_nodeID()); }

    const ConstTreeNodeIter & get_node_iter() const { return nodeiter; }

  protected:
    const TreeBlock * blockptr;

    mutable unsigned int treeindex;

    mutable ConstTreeNodeIter nodeiter;
  };


}  // End namespace Tree
}  // End namespace QuickFlash


#endif  // QUICKFLASH_TREE_BLOCK_HPP
