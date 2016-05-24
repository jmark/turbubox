// C++ auxiliary template file quickflash_tree_node.tcc

/*
  By Nathan C. Hearn
     March 6, 2007

  For inclusion only by quickflash_tree_node.hpp.
*/


#ifndef QUICKFLASH_TREE_NODE_TCC
#define QUICKFLASH_TREE_NODE_TCC


#include "quickflash_tree_nodeid.hpp"



namespace QuickFlash
{
namespace Tree
{

// Node Member Definitions

template<class D>
void Node<D>::reset()
  {
  nodeID.reset();

  rootNode = true;

  parent = 0;

  prevSibling = 0;
  nextSibling = 0;

  delete_children();

  data.reset();
  }


template<class D>
void Node<D>::reset(const Node<D> & source)
  {
  if(&source != this)
    {
    // Copy only the data, not any tree structure, treating this node as root

    reset();

    data.reset(source.data);
    }
  }


template<class D>
Node<D> * Node<D>::get_nextNode()
  {
  /*
    Finds the next node in a particular tree-traversal walk.  The algorithm
    for finding the next node is as follows:

    1. Return a pointer to the first child in children.
    2. If no children exist, return a pointer to the next sibling of this
       node.
    3. If the call to getNextSibling() is null, call getNextSibling on my
       parent node.  Return the node that is found.
    4. If no node is found, repeat step 3 for the parent of my parent (and
       so on) until either a node is found or the root node is reached.
       (The root node has no siblings.)  Return the first node found, or else 
       return a null pointer.

    If one starts with the root node, this function can be used to run through
    all of the other nodes in the tree, one by one.
  */

  Node<D> * nextNode = 0;

  // Pointer to first child (if one exists)

  if (!(children.empty()))
    nextNode = children.front()->get_nodePtr();
  else
    {
    // Find sibling (or parent with sibling)

    Node<D> * nodePtr = this;
    bool foundNode = false;

    while ((!foundNode) && (!(nodePtr->is_root())))
      {
      Node<D> * testNodePtr = nodePtr->get_nextSibling();

      // If I have a sibling, use it; else, ask my parent if it has a sibling

      if (testNodePtr != 0)
	{
	nodePtr = testNodePtr;
	foundNode = true;
	}
      else
	nodePtr = nodePtr->parent;
      }

    if (foundNode)
      nextNode = nodePtr;
    }

  return nextNode;
  }


template<class D>
  const Node<D> * Node<D>::get_nextNode() const
  {
  /*
    Finds the next node in a particular tree-traversal walk.  The algorithm
    for finding the next node is as follows:

    1. Return a pointer to the first child in children.
    2. If no children exist, return a pointer to the next sibling of this
       node.
    3. If the call to getNextSibling() is null, call getNextSibling on my
       parent node.  Return the node that is found.
    4. If no node is found, repeat step 3 for the parent of my parent (and
       so on) until either a node is found or the root node is reached.
       (The root node has no siblings.)  Return the first node found, or else 
       return a null pointer.

    If one starts with the root node, this function can be used to run through
    all of the other nodes in the tree, one by one.
  */

  const Node<D> * nextNode = 0;

  // Pointer to first child (if one exists)

  if (!(children.empty()))
    nextNode = children.front()->get_nodePtr();
  else
    {
    // Find sibling (or parent with sibling)

    const Node<D> * nodePtr = this;
    bool foundNode = false;

    while((!foundNode) && (!(nodePtr->is_root())))
      {
      const Node<D> * testNodePtr = nodePtr->get_nextSibling();

      // If I have a sibling, use it; else, ask my parent if it has a sibling

      if(testNodePtr != 0)
	{
	nodePtr = testNodePtr;
	foundNode = true;
	}
      else
	nodePtr = nodePtr->parent;
      }

    if (foundNode)
      nextNode = nodePtr;
    }

  return nextNode;
  }


template<class D>
Node<D> * Node<D>::get_prevNode()
  {
  /*
    Finds the previous node in a particular tree-traversal walk.  The algorithm
    for finding the next node is as follows:

    1. If I am root, return null.
    2. Return a pointer to the last node in the branch of the previous sibling.
    3. If no previous siblings exist, return a pointer to the parent of this 
       node.

    If one starts with the last node in the tree, this function can be used to
    run through all of the other nodes in the tree, one by one, ending with
    the root node.
  */

  Node<D> * prevNode = 0;

  if (!rootNode)
    {
    // Pointer to previous sibling, if one exists

    Node<D> * nodePtr = get_prevSibling();

    if (nodePtr != 0)
      prevNode = nodePtr->get_lastNode_branch();
    else
      prevNode = parent;
    }

  return prevNode;
  }


template<class D>
const Node<D> * Node<D>::get_prevNode() const
  {
  /*
    Finds the previous node in a particular tree-traversal walk.  The algorithm
    for finding the next node is as follows:

    1. If I am root, return null.
    2. Return a pointer to the last node in the branch of the previous sibling.
    3. If no previous siblings exist, return a pointer to the parent of this 
       node.

    If one starts with the last node in the tree, this function can be used to
    run through all of the other nodes in the tree, one by one, ending with
    the root node.
  */

  const Node<D> * prevNode = 0;

  if (!rootNode)
    {
    // Pointer to previous sibling, if one exists

    const Node<D> * nodePtr = get_prevSibling();

    if (nodePtr != 0)
      prevNode = nodePtr->get_lastNode_branch();
    else
      prevNode = parent;
    }

  return prevNode;
  }


template<class D>
Node<D> * Node<D>::get_nextNode_skipChildren()
  {
  /*
    Returns a pointer to the next node, skipping the children of the current
    node.
  */

  Node<D> * nodePtr = this;
  Node<D> * nextSib = get_nextSibling();

  while (!(nodePtr->is_root()) && (nextSib == 0))
    {
    nodePtr = nodePtr->get_parent();
    nextSib = nodePtr->get_nextSibling();
    }

  return nextSib;
  }


template<class D>
const Node<D> * Node<D>::get_nextNode_skipChildren() const
  {
  /*
    Returns a pointer to the next node, skipping the children of the current
    node.
  */

  const Node<D> * nodePtr = this;
  const Node<D> * nextSib = get_nextSibling();

  while (!(nodePtr->is_root()) && (nextSib == 0))
    {
    nodePtr = nodePtr->get_parent();
    nextSib = nodePtr->get_nextSibling();
    }

  return nextSib;
  }


template<class D>
Node<D> * Node<D>::get_lastNode_branch()
  {
  /*
    Returns a pointer to the "last" node in this node's branch.  Starting
    with this node, the function follows down the path of last children,
    returning a pointer to the last node reached that has no children.
  */

  Node<D> * nodePtr = this;

  while (nodePtr->children_present())
    nodePtr = nodePtr->get_lastChild();

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_lastNode_branch() const
  {
  /*
    Returns a pointer to the "last" node in this node's branch.  Starting
    with this node, the function follows down the path of last children,
    returning a pointer to the last node reached that has no children.
  */

  const Node<D> * nodePtr = this;

  while (nodePtr->children_present())
    nodePtr = nodePtr->get_lastChild();

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_firstNode_branchLevel(const unsigned int level)
  {
  /*
    Returns a pointer to the first node in this node's branch that has the 
    specified level.  If none is found, a null pointer is returned.
  */

  Node<D> * nodePtr = this;

  const unsigned int current_level = get_level();

  if (level < current_level)
    nodePtr = 0;
  else if (level > current_level)
    {
    bool found_it = false;

    while ((nodePtr != 0) && (!found_it))
      {
      const unsigned int node_level = nodePtr->get_level();

      if (node_level == level)
	found_it = true;
      else if (node_level <= current_level)
	nodePtr = 0;
      else
	nodePtr = nodePtr->get_nextNode();
      }
    }

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_firstNode_branchLevel(const unsigned int level)
  const
  {
  /*
    Returns a pointer to the first node in this node's branch that has the 
    specified level.  If none is found, a null pointer is returned.
  */

  const Node<D> * nodePtr = this;

  const unsigned int current_level = get_level();

  if (level < current_level)
    nodePtr = 0;
  else if (level > current_level)
    {
    bool found_it = false;

    while ((nodePtr != 0) && (!found_it))
      {
      const unsigned int node_level = nodePtr->get_level();

      if (node_level == level)
	found_it = true;
      else if (node_level <= current_level)
	nodePtr = 0;
      else
	nodePtr = nodePtr->get_nextNode();
      }
    }

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_rootNode()
  {
  Node<D> * nodePtr = this;

  while (!(nodePtr->is_root()))
    nodePtr = nodePtr->get_parent();

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_rootNode() const
  {
  const Node<D> * nodePtr = this;

  while (!(nodePtr->is_root()))
    nodePtr = nodePtr->get_parent();

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_nextNode_level(const unsigned int level)
  {
  /*
    Returns the next node in the tree (after the current one) that is at
    the level denoted by the argument.
  */

  Node<D> * nodePtr = 0;

  if (level > 0)
    {
    bool found_it = false;

    nodePtr = get_nextNode();

    while ((nodePtr != 0) && (!found_it))
      {
      const unsigned int node_level = nodePtr->get_level();

      if (node_level == level)
	found_it = true;
      else if (node_level < level)
	{
	Node<D> * testPtr = nodePtr->get_firstNode_branchLevel(level);

	if (testPtr != 0)
	  {
	  nodePtr = testPtr;
	  found_it = true;
	  }
	}

      if (!found_it)
	nodePtr = nodePtr->get_nextNode_skipChildren();
      }
    }

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_nextNode_level(const unsigned int level) const
  {
  /*
    Returns the next node in the tree (after the current one) that is at
    the level denoted by the argument.
  */

  const Node<D> * nodePtr = 0;

  if (level > 0)
    {
    bool found_it = false;

    nodePtr = get_nextNode();

    while ((nodePtr != 0) && (!found_it))
      {
      const unsigned int node_level = nodePtr->get_level();

      if (node_level == level)
	found_it = true;
      else if (node_level < level)
	{
	const Node<D> * testPtr = nodePtr->get_firstNode_branchLevel(level);

	if (testPtr != 0)
	  {
	  nodePtr = testPtr;
	  found_it = true;
	  }
	}

      if (!found_it)
	nodePtr = nodePtr->get_nextNode_skipChildren();
      }
    }

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_prevNode_level(const unsigned int level)
  {
  /*
    Returns the previous node in the tree (before the current one) that is at
    the level denoted by the argument.
  */

  Node<D> * nodePtr = get_prevNode();

  bool found_it = false;

  while ((nodePtr != 0) && (!found_it))
    {
    const unsigned int current_level = nodePtr->get_level();

    if (current_level < level)
      nodePtr = nodePtr->get_prevNode();
    else
      {
      const unsigned int delta = current_level - level;

      for (int i = 0; i < delta; i++)
	nodePtr = nodePtr->get_parent();

      found_it = true;
      }
    }

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_prevNode_level(const unsigned int level) const
  {
  /*
    Returns the previous node in the tree (before the current one) that is at
    the level denoted by the argument.
  */

  const Node<D> * nodePtr = get_prevNode();

  bool found_it = false;

  while ((nodePtr != 0) && (!found_it))
    {
    const unsigned int current_level = nodePtr->get_level();

    if (current_level < level)
      nodePtr = nodePtr->get_prevNode();
    else
      {
      const unsigned int delta = current_level - level;

      for (int i = 0; i < delta; i++)
	nodePtr = nodePtr->get_parent();

      found_it = true;
      }
    }

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_node(const NodeID & id_vect)
  {
  Node<D> * nodePtr = this;

  if (!(branch_member(id_vect)))
    nodePtr = get_rootNode();

  const unsigned int vect_level = id_vect.get_level();

  unsigned int current_level = nodePtr->get_level();

  while ((nodePtr != 0) && (current_level < vect_level))
    nodePtr = nodePtr->get_child(id_vect[current_level++]);

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_node(const NodeID & id_vect) const
  {
  const Node<D> * nodePtr = this;

  if (!(branch_member(id_vect)))
    nodePtr = get_rootNode();

  const unsigned int vect_level = id_vect.get_level();

  unsigned int current_level = nodePtr->get_level();

  while ((nodePtr != 0) && (current_level < vect_level))
    nodePtr = nodePtr->get_child(id_vect[current_level++]);

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_branch(const NodeID & id_vect)
  {
  Node<D> * nodePtr = this;

  if (!(branch_member(id_vect)))
    nodePtr = get_rootNode();

  const unsigned int vect_level = id_vect.get_level();

  bool found_it = false;

  unsigned int current_level = nodePtr->get_level();

  while (!found_it)
    {
    Node<D> * testPtr = nodePtr->get_child(id_vect[current_level++]);

    if (testPtr != 0)
      {
      nodePtr = testPtr;

      if (current_level >= vect_level)
	found_it = true;
      }
    else
      found_it = true;
    }

  return nodePtr;
  }


template<class D>
const Node<D> * Node<D>::get_branch(const NodeID & id_vect) const
  {
  const Node<D> * nodePtr = this;

  if (!(branch_member(id_vect)))
    nodePtr = get_rootNode();

  const unsigned int vect_level = id_vect.get_level();

  bool found_it = false;

  unsigned int current_level = nodePtr->get_level();

  while (!found_it)
    {
    const Node<D> * testPtr = nodePtr->get_child(id_vect[current_level++]);

    if (testPtr != 0)
      {
      nodePtr = testPtr;

      if (current_level >= vect_level)
	found_it = true;
      }
    else
      found_it = true;
    }

  return nodePtr;
  }


template<class D>
Node<D> * Node<D>::get_firstChild()
  {
  Node<D> * childPtr = 0;

  if(children_present())
    childPtr = children.front();

  return childPtr;
  }


template<class D>
const Node<D> * Node<D>::get_firstChild() const
  {
  const Node<D> * childPtr = 0;

  if(children_present())
    childPtr = children.front();

  return childPtr;
  }


template<class D>
Node<D> * Node<D>::get_lastChild()
  {
  Node<D> * childPtr = 0;

  if(children_present())
    childPtr = children.back();

  return childPtr;
  }


template<class D>
const Node<D> * Node<D>::get_lastChild() const
  {
  const Node<D> * childPtr = 0;

  if(children_present())
    childPtr = children.back();

  return childPtr;
  }


template<class D>
Node<D> * Node<D>::get_child(const LevelID_Type child_id)
  {
  // Return null if child is not present

  Node<D> * childPtr = 0;

  if(children_present())
    {
    childPtr = get_firstChild();

    bool foundChild = false;

    const unsigned int childID_index = get_level();

    while((childPtr != 0) && (!foundChild))
      {
      if(childPtr->get_nodeID()[childID_index] == child_id)
	foundChild = true;
      else
	childPtr = childPtr->get_nextSibling();
      }
    }

  return childPtr;
  }


template<class D>
const Node<D> * Node<D>::get_child(const LevelID_Type child_id) const
  {
  // Return null if child is not present

  const Node<D> * childPtr = 0;

  if(children_present())
    {
    childPtr = get_firstChild();

    bool foundChild = false;

    const unsigned int childID_index = get_level();

    while((childPtr != 0) && (!foundChild))
      {
      if(childPtr->get_nodeID()[childID_index] == child_id)
	foundChild = true;
      else
	childPtr = childPtr->get_nextSibling();
      }
    }

  return childPtr;
  }


template<class D>
Node<D> * Node<D>::create_child(const LevelID_Type child_id)
  {
  // Return null if child already exists

  Node<D> * childPtr = 0;

  const unsigned int childID_index = get_level();

  if(children_present())
    {
    // Find the sibling that comes after the one to be created

    typename NodePtrList::iterator sibNode = children.begin();
    const typename NodePtrList::iterator endSibNode = children.end();

    bool foundSib = false;

    LevelID_Type sibID = 0;

    while((sibNode != endSibNode) && (!foundSib))
      {
      Node<D> * sibNodePtr = *sibNode;

      sibID = sibNodePtr->get_nodeID()[childID_index];

      if(sibID < child_id)
	++sibNode;
      else
	foundSib = true;
      }

    // Only create child if one with childID does not exist

    if(sibID != child_id)
      {
      Node<D> * prevSibPtr = 0;
      Node<D> * nextSibPtr = 0;

      if(foundSib)
	{
	nextSibPtr = *sibNode;
	prevSibPtr = nextSibPtr->get_prevSibling();
	}
      else
	prevSibPtr = children.back();

      // Insert the new child before sibNode

      childPtr = new Node<D>;

      children.insert(sibNode, childPtr);

      // Set the child's siblings

      childPtr->prevSibling = prevSibPtr;
      childPtr->nextSibling = nextSibPtr;

      // Reset the siblings pointers

      if(prevSibPtr != 0)
	prevSibPtr->nextSibling = childPtr;

      if(nextSibPtr != 0)
	nextSibPtr->prevSibling = childPtr;
      }
    }
  else
    {
    childPtr = new Node<D>;

    children.insert(children.end(), childPtr);
    }

  // If child created, set up parent, nodeID, and non-root status

  if(childPtr != 0)
    {
    childPtr->rootNode = false;

    childPtr->parent = this;

    // Build the node ID vector

    const unsigned int childLevel = childID_index + 1;

    NodeID & child_nodeID = childPtr->nodeID;

    child_nodeID.resize(childLevel);

    // Node ID is my ID with childID appended to end

    for(unsigned int i = 0; i < childID_index; i++)
      child_nodeID[i] = nodeID[i];

    child_nodeID[childID_index] = child_id;
    }

  return childPtr;
  }


template<class D>
void Node<D>::delete_child(const LevelID_Type child_id)
  {
  // Find the child

  typename NodePtrList::iterator child_iter = children.begin();
  const typename NodePtrList::iterator end_child_iter = children.end();

  bool found_child = false;

  while ((child_iter != end_child_iter) && (!found_child))
    {
    if ((*child_iter)->get_levelID() == child_id)
      found_child = true;
    else
      ++child_iter;
    }

  if (found_child)
    {
    Node<D> * childNode = *child_iter;

    // Reset the sibling pointers of childNode's siblings

    Node<D> * prev_sib = childNode->get_prevSibling();
    Node<D> * next_sib = childNode->get_nextSibling();

    if (prev_sib != 0)
      prev_sib->nextSibling = next_sib;

    if (next_sib != 0)
      next_sib->prevSibling = prev_sib;

    // Delete childNode

    children.erase(child_iter);
    }
  }


template<class D>
void Node<D>::delete_children()
  {
  while(!(children.empty()))
    {
    delete children.back();
    children.pop_back();
    }
  }


// NodeIter Member Definitions

template<class D>
void NodeIter<D>::reset(const NodeIter<D> & source)
  {
  if (&source != this)
    {
    root = source.root;
    branch = source.branch;
    current = source.current;
    }
  }


template<class D>
void NodeIter<D>::reset(TreeNode * branch_node)
  {
  if (branch_node == 0)
    reset();
  else
    {
    root = branch_node->get_rootNode();
    branch = branch_node;
    current = branch_node;
    }
  }


template<class D>
int NodeIter<D>::next() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::prev() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::next_skipBranch() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_skipChildren();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::goto_node(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    current = root->get_node(nodeID);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::goto_branch(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    current = root->get_branch(nodeID);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::level_next() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_level();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::level_prev() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode_level();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::level_next(const unsigned int level) const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_level(level);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::level_prev(const unsigned int level) const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode_level(level);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int NodeIter<D>::set_branch() const
  {
  int ret_val = -1;  // Invalid branch

  branch = current;

  if (branch != 0)
    ret_val = 0;  // Valid branch
  else
    branch = root;  // Don't let branch be invalid if root is valid

  return ret_val;
  }


template<class D>
int NodeIter<D>::set_branch(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid branch

  // Set the branch pointer to a SPECIFIC node ONLY if the node exists

  if (root != 0)
    {
    // Don't replace branch with an invalid node

    TreeNode test_branch = root->get_node(nodeID);

    if (test_branch != 0)
      {
      branch = test_branch;
      ret_val = 0;  // Valid branch
      }
    }

  return ret_val;
  }


template<class D>
int NodeIter<D>::goto_last() const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    {
    current = root->get_lastNode();
    ret_val = 0;  // Valid node
    }

  return ret_val;
  }


template<class D>
int NodeIter<D>::branch_next() const
  {
  int ret_val = -1;  // Valid node

  if (branch != 0)
    {
    next();

    if (current != 0)
      {
      if (branch->branch_member(current))
	ret_val = 0;  // Valid node
      else
	current = 0;  // Node not in branch
      }
    }

  return ret_val;
  }


template<class D>
int NodeIter<D>::branch_prev() const
  {
  int ret_val = -1;  // Valid node

  if (branch != 0)
    {
    prev();

    if (current != 0)
      {
      if (branch->branch_member(current))
	ret_val = 0;  // Valid node
      else
	current = 0;  // Node not in branch
      }
    }

  return ret_val;
  }


template<class D>
int NodeIter<D>::branch_last() const
  {
  int ret_val = -1;  // Invalid node

  if (branch != 0)
    {
    current = branch->get_lastNode_branch();
    ret_val = 0;  // Valid node
    }

  return ret_val;
  }


template<class D>
int NodeIter<D>::set_data(const D & source_data)
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    {
    current->set_data(source_data);
    ret_val = 0;  // Stored data in valid node
    }

  return ret_val;
  }


// ConstNodeIter Member Definitions

template<class D>
void ConstNodeIter<D>::reset(const ConstNodeIter<D> & source)
  {
  if (&source != this)
    {
    root = source.root;
    branch = source.branch;
    current = source.current;
    }
  }


template<class D>
void ConstNodeIter<D>::reset(const NodeIter<D> & source)
  {
  root = source.root;
  branch = source.branch;
  current = source.current;
  }


template<class D>
void ConstNodeIter<D>::reset(const TreeNode * branch_node)
  {
  if (branch_node == 0)
    reset();
  else
    {
    root = branch_node->get_rootNode();
    branch = branch_node;
    current = branch_node;
    }
  }


template<class D>
int ConstNodeIter<D>::next() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::prev() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::next_skipBranch() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_skipChildren();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::goto_node(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    current = root->get_node(nodeID);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::goto_branch(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    current = root->get_branch(nodeID);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::level_next() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_level();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::level_prev() const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode_level();

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::level_next(const unsigned int level) const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_nextNode_level(level);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::level_prev(const unsigned int level) const
  {
  int ret_val = -1;  // Invalid node

  if (current != 0)
    current = current->get_prevNode_level(level);

  if (current != 0)
    ret_val = 0;  // Valid node

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::set_branch() const
  {
  int ret_val = -1;  // Invalid branch

  branch = current;

  if (branch != 0)
    ret_val = 0;  // Valid branch
  else
    branch = root;  // Don't let branch be invalid if root is valid

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::set_branch(const NodeID & nodeID) const
  {
  int ret_val = -1;  // Invalid branch

  // Set the branch pointer to a SPECIFIC node ONLY if the node exists

  if (root != 0)
    {
    // Don't replace branch with an invalid node

    TreeNode test_branch = root->get_node(nodeID);

    if (test_branch != 0)
      {
      branch = test_branch;
      ret_val = 0;  // Valid branch
      }
    }

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::goto_last() const
  {
  int ret_val = -1;  // Invalid node

  if (root != 0)
    {
    current = root->get_lastNode();
    ret_val = 0;  // Valid node
    }

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::branch_next() const
  {
  int ret_val = -1;  // Valid node

  if (branch != 0)
    {
    next();

    if (current != 0)
      {
      if (branch->branch_member(current))
	ret_val = 0;  // Valid node
      else
	current = 0;  // Node not in branch
      }
    }

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::branch_prev() const
  {
  int ret_val = -1;  // Valid node

  if (branch != 0)
    {
    prev();

    if (current != 0)
      {
      if (branch->branch_member(current))
	ret_val = 0;  // Valid node
      else
	current = 0;  // Node not in branch
      }
    }

  return ret_val;
  }


template<class D>
int ConstNodeIter<D>::branch_last() const
  {
  int ret_val = -1;  // Invalid node

  if (branch != 0)
    {
    current = branch->get_lastNode_branch();
    ret_val = 0;  // Valid node
    }

  return ret_val;
  }


}  // End namespace Tree
}  // End namespace QuickFlash


#endif  // QUICKFLASH_TREE_NODE_TCC
