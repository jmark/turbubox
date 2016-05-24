// C++ header file quickflash_tree_nodeid.hpp

/*
  By Nathan C. Hearn
     October 12, 2006

  Node ID vector for use with identifying tree nodes.

  NOTE: Since level zero (the root node) has no entries, levels and indexes
  are off by one, i.e., the information for level i is stored at index (i-1).
  Think of level zero as the basement, level one is at index zero,
  level two is at index 1, level three is at index 2, and so on.
*/


#ifndef QUICKFLASH_TREE_NODEID_HPP
#define QUICKFLASH_TREE_NODEID_HPP


#include <vector>

namespace QuickFlash
{
namespace Tree
{

typedef unsigned int LevelID_Type;

class NodeID
  {
  private:
    typedef std::vector<LevelID_Type> IDVect;

  public:
    NodeID() : id() { }
    NodeID(const IDVect & id_vector) : id(id_vector) { }
    NodeID(const NodeID & source) : id(source.id) { }
    ~NodeID() { }

    NodeID & operator=(const NodeID & source)
      {
      reset(source);
      return *this;
      }

    LevelID_Type & operator[](const unsigned int index) 
      { return id[index]; }
    const LevelID_Type & operator[](const unsigned int index) const
      { return id[index]; }

    unsigned int get_level() const { return id.size(); }
    unsigned int get_numLevels() const { return (id.size() + 1); }

    LevelID_Type get_level_id() const
      { return get_level_id(get_level()); }

    LevelID_Type get_level_id(const unsigned int level) const
      {
      LevelID_Type levelID = 0;

      if (level > 0)
	{
	const unsigned int index = level - 1;
	levelID = id[index];
	}

      return levelID;
      }

    bool branch_member(const NodeID & child_id) const;

    bool in_branch(const NodeID & parent_id) const
      { return parent_id.branch_member(*this); }

    void reset() { id.clear(); }

    void reset(const IDVect & id_vector) { id = id_vector; }

    void reset(const NodeID & source)
      {
      if (&source != this)
	id = source.id;
      }

    void add_level(const LevelID_Type level_id_elem) 
      { id.push_back(level_id_elem); }

    void resize(const unsigned int new_level, 
		const LevelID_Type default_value=0)
      { id.resize(new_level, default_value); }

    void invert_level_bit(const unsigned int level,
			  const unsigned int bit_index);

    void invert_bits(const unsigned int bit_index);

  private:
    std::vector<LevelID_Type> id;
  };


}  // End namespace Tree
}  // End namespace QuickFlash


#endif  // QUICKFLASH_TREE_NODEID_HPP
