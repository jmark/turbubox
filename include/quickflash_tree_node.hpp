// C++ header file quickflash_tree_node.hpp

/*
  By Nathan C. Hearn
     October 25, 2006

  Template class for tree nodes with a data container.

  Tree data D must have both default and copy-style reset functions.

  NOTE: In NodeIter do NOT let the branch pointer be null if root is not null.
*/


#ifndef QUICKFLASH_TREE_NODE_HPP
#define QUICKFLASH_TREE_NODE_HPP

#include <vector>
#include <list>
#include "quickflash_tree_nodeid.hpp"


namespace QuickFlash
{
namespace Tree
{

// Forward declarations

template<class D>
class Node;

template<class D>
class NodeIter;

template<class D>
class ConstNodeIter;


// Node Class Definition

template<class D>
class Node
  {
  public:
    typedef NodeIter<D> Iter;
    typedef ConstNodeIter<D> ConstIter;

  private:
    typedef typename std::list<Node<D> *> NodePtrList;

  public:
    Node() :
      nodeID(), rootNode(true), parent(0), prevSibling(0), nextSibling(0), 
      children(0), data()
      { reset(); }
    Node(const Node<D> & source) :
      nodeID(), rootNode(true), parent(0), prevSibling(0), nextSibling(0), 
      children(0), data(source.data)
      { reset(source); }
    ~Node() { delete_children(); }

    Node<D> & operator=(const Node<D> & source)
      {
      reset(source);
      return *this;
      }

    void reset();
    void reset(const Node<D> & source);

    const NodeID & get_nodeID() const { return nodeID; }
    
    unsigned int get_level() const { return nodeID.get_level(); }

    LevelID_Type get_levelID() const { return nodeID.get_level_id(); }

    LevelID_Type get_levelID(const unsigned int level) const
      { return nodeID.get_level_id(level); }

    void get_node_iter(Iter & node_iter)
      { node_iter.reset(this); }

    void get_node_iter(ConstIter & node_iter) const
      { node_iter.reset(this); }

    Iter get_node_iter() { return Iter(this); }
    ConstIter get_node_iter() const { return ConstIter(this); }

    bool is_root() const { return rootNode; }

    bool children_present() const { return !(children.empty()); }

    bool child_present(const LevelID_Type child_id) const
      { return (get_child(child_id) != 0); }

    bool branch_member(const NodeID & child_id) const
      { return nodeID.branch_member(child_id); }

    bool branch_member(const Node<D> * child_ptr) const
      { return branch_member(child_ptr->get_nodeID()); }

    bool branch_member(const Node<D> & child_node) const
      { return branch_member(child_node.get_nodeID()); }

    Node<D> * get_nodePtr() { return this; }
    const Node<D> * get_nodePtr() const { return this; }
    
    Node<D> * get_parent() { return parent; }
    const Node<D> * get_parent() const { return parent; }

    Node<D> * get_prevSibling() { return prevSibling; }
    const Node<D> * get_prevSibling() const { return prevSibling; }

    Node<D> * get_nextSibling() { return nextSibling; }
    const Node<D> * get_nextSibling() const { return nextSibling; }

    Node<D> * get_nextNode();
    const Node<D> * get_nextNode() const;

    Node<D> * get_prevNode();
    const Node<D> * get_prevNode() const;

    Node<D> * get_nextNode_skipChildren();
    const Node<D> * get_nextNode_skipChildren() const;

    Node<D> * get_lastNode_branch();
    const Node<D> * get_lastNode_branch() const;

    Node<D> * get_firstNode_branchLevel(const unsigned int level);
    const Node<D> * get_firstNode_branchLevel(const unsigned int level) const;

    Node<D> * get_rootNode();
    const Node<D> * get_rootNode() const;

    Node<D> * get_lastNode() { return get_rootNode()->get_lastNode_branch(); }
    const Node<D> * get_lastNode() const
      { return get_rootNode()->get_lastNode_branch(); }

    Node<D> * get_nextNode_level() { return get_nextNode_level(get_level()); }
    const Node<D> * get_nextNode_level() const
      { return get_nextNode_level(get_level()); }

    Node<D> * get_nextNode_level(const unsigned int level);
    const Node<D> * get_nextNode_level(const unsigned int level) const;

    Node<D> * get_prevNode_level() { return get_prevNode_level(get_level()); }
    const Node<D> * get_prevNode_level() const
      { return get_prevNode_level(get_level()); }

    Node<D> * get_prevNode_level(const unsigned int level);
    const Node<D> * get_prevNode_level(const unsigned int level) const;

    Node<D> * get_node(const NodeID & id_vect);
    const Node<D> * get_node(const NodeID & id_vect) const;

    Node<D> * get_branch(const NodeID & id_vect);
    const Node<D> * get_branch(const NodeID & id_vect) const;

    Node<D> * get_firstChild();
    const Node<D> * get_firstChild() const;

    Node<D> * get_lastChild();
    const Node<D> * get_lastChild() const;

    Node<D> * get_child(const LevelID_Type child_id);
    const Node<D> * get_child(const LevelID_Type child_id) const;

    Node<D> * create_child(const LevelID_Type child_id);

    void delete_child(const LevelID_Type child_id);

    void delete_children();

    void set_data(const D & dataSource) { data = dataSource; }

    D & get_data() { return data; }
    const D & get_data() const { return data; }

  private:
    NodeID nodeID;

    bool rootNode;

    Node<D> * parent;

    Node<D> * prevSibling;
    Node<D> * nextSibling;
    
    NodePtrList children;

    // std::map<LevelID_Type, Node<D> *> child_map;

    D data;
  };


// NodeIter Class Definition

template<class D>
class NodeIter
  {
  friend class ConstNodeIter<D>;

  public:
    typedef Node<D> TreeNode;

  public:
    NodeIter() : root(0), branch(0), current(0) { reset(); }
    NodeIter(TreeNode & branch_node) :
      root(0), branch(0), current(0)
      { reset(&branch_node); }
    NodeIter(TreeNode * branch_node) :
      root(0), branch(0), current(0)
      { reset(branch_node); }
    NodeIter(const NodeIter<D> & source) :
      root(0), branch(0), current(0)
      { reset(source); }
    ~NodeIter() { }

    NodeIter<D> & operator=(const NodeIter<D> & source)
      {
      reset(source);
      return *this;
      }

    NodeIter<D> & operator=(TreeNode & branch_node)
      {
      reset(&branch_node);
      return *this;
      }

    NodeIter<D> & operator=(TreeNode * branch_node)
      {
      reset(branch_node);
      return *this;
      }

    void reset()
      {
      root = 0;
      branch = 0;
      current = 0;
      }

    void reset(const NodeIter<D> & source);
    void reset(TreeNode * branch_node);

    bool valid_node() const { return (current != 0); }

    bool valid_tree() const 
      { return ((root != 0) && (branch != 0) && (current != 0)); }

    void set_invalid() const { current = 0; }

    const NodeID & get_nodeID() const { return current->get_nodeID(); }

    int next() const;
    int prev() const;

    int next_skipBranch() const;

    int goto_node(const NodeID & nodeID) const;

    int goto_branch(const NodeID & nodeID) const;

    int level_next() const;
    int level_prev() const;

    int level_next(const unsigned int level) const;
    int level_prev(const unsigned int level) const;

    int set_branch() const;
    int set_branch(const NodeID & nodeID) const;

    void goto_root() const { current = root; }
    int goto_last() const;

    int branch_next() const;
    int branch_prev() const;

    void branch_root() const { current = branch; }
    int branch_last() const;

    int set_data(const D & source_data);

    D & get_data() { return current->get_data(); }
    const D & get_data() const { return current->get_data(); }

    TreeNode & get_current() { return *current; }
    const TreeNode & get_current() const { return *current; }

    bool children_present() const
      {
      bool present = false;

      if (current != 0)
	if (current->children_present())
	  present = true;

      return present;
      }

  private:
    mutable TreeNode * root;
    mutable TreeNode * branch;
    mutable TreeNode * current;
  }; 


// ConstNodeIter Class Definition

template<class D>
class ConstNodeIter
  {
  public:
    typedef Node<D> TreeNode;

  public:
    ConstNodeIter() : root(0), branch(0), current(0) { reset(); }
    ConstNodeIter(const TreeNode & branch_node) :
      root(0), branch(0), current(0)
      { reset(&branch_node); }
    ConstNodeIter(const TreeNode * branch_node) :
      root(0), branch(0), current(0)
      { reset(branch_node); }
    ConstNodeIter(const NodeIter<D> & source) :
      root(0), branch(0), current(0)
      { reset(source); }
    ConstNodeIter(const ConstNodeIter<D> & source) :
      root(0), branch(0), current(0)
      { reset(source); }
    ~ConstNodeIter() { }

    ConstNodeIter<D> & operator=(const ConstNodeIter<D> & source)
      {
      reset(source);
      return *this;
      }

    ConstNodeIter<D> & operator=(const NodeIter<D> & source)
      {
      reset(source);
      return *this;
      }

    ConstNodeIter<D> & operator=(const TreeNode & branch_node)
      {
      reset(&branch_node);
      return *this;
      }

    ConstNodeIter<D> & operator=(const TreeNode * branch_node)
      {
      reset(branch_node);
      return *this;
      }

    void reset()
      {
      root = 0;
      branch = 0;
      current = 0;
      }

    void reset(const ConstNodeIter<D> & source);
    void reset(const NodeIter<D> & source);
    void reset(const TreeNode * branch_node);

    bool valid_node() const { return (current != 0); }

    bool valid_tree() const 
      { return ((root != 0) && (branch != 0) && (current != 0)); }

    void set_invalid() const { current = 0; }

    const NodeID & get_nodeID() const { return current->get_nodeID(); }

    int next() const;
    int prev() const;

    int next_skipBranch() const;

    int goto_node(const NodeID & nodeID) const;

    int goto_branch(const NodeID & nodeID) const;

    int level_next() const;
    int level_prev() const;

    int level_next(const unsigned int level) const;
    int level_prev(const unsigned int level) const;

    int set_branch() const;
    int set_branch(const NodeID & nodeID) const;

    void goto_root() const { current = root; }
    int goto_last() const;

    int branch_next() const;
    int branch_prev() const;

    void branch_root() const { current = branch; }
    int branch_last() const;

    const D & get_data() const { return current->get_data(); }

    const TreeNode & get_current() const { return *current; }

    bool children_present() const
      {
      bool present = false;

      if (current != 0)
	if (current->children_present())
	  present = true;

      return present;
      }

  private:
    mutable const TreeNode * root;
    mutable const TreeNode * branch;
    mutable const TreeNode * current;
  }; 


}  // End namespace Tree
}  // End namespace QuickFlash


// Include the function definitions

#include "quickflash_tree_node.tcc"


#endif  // QUICKFLASH_TREE_NODE_HPP
