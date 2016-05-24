// C++ header file quickflash_registry.hpp

/*
  By Nathan C. Hearn
     September 19, 2006

  Template for building registries.
*/


#ifndef QUICKFLASH_REGISTRY_HPP
#define QUICKFLASH_REGISTRY_HPP

#include <vector>
#include <map>


namespace QuickFlash
{
namespace Registry
{

template<typename KeyType, typename DataType>
class Registry
  {
  public :
    typedef typename std::map<KeyType, unsigned int> KeyIndexMap;

  public :

    Registry() : reg(), reg_index_map() { }

    ~Registry() { reset(); }

    void reset();

    DataType & operator[](const unsigned int reg_index)
      { return *(reg[reg_index]); }
    const DataType & operator[](const unsigned int reg_index) const
      { return *(reg[reg_index]); }

    DataType & operator[](const KeyType & key)
      { return *(reg[get_reg_index(key)]); }
    const DataType & operator[](const KeyType & key) const
      { return *(reg[get_reg_index(key)]); }

    unsigned int get_num_entries() const { return reg.size(); }

    bool key_present(const KeyType & key) const;

    void get_key_list(std::vector<KeyType> & key_list) const;

    unsigned int get_reg_index(const KeyType & key) const;

    unsigned int add_new(const KeyType & key);
    unsigned int add_new(const KeyType & key, const DataType & original);

  private :

    std::vector<DataType *> reg;

    KeyIndexMap reg_index_map;
  };


template<typename KeyType, typename DataType>
void Registry<KeyType, DataType>::reset()
  {
  // Clean out the map

  reg_index_map.clear();  // I hope someone is not pointing to these pointers

  // Delete the objects in the registry

  while (!(reg.empty()))
    {
    delete reg.back();
    reg.back() = 0;

    reg.pop_back();
    }
  }


template<typename KeyType, typename DataType>
bool Registry<KeyType, DataType>::key_present(const KeyType & key) const
  {
  bool present = false;

  const typename KeyIndexMap::const_iterator map_iter 
    = reg_index_map.find(key);

  if (map_iter != reg_index_map.end())
    present = true;

  return present;
  }


template<typename KeyType, typename DataType>
void Registry<KeyType, DataType>::get_key_list(std::vector<KeyType> & key_list)
  const
  {
  const unsigned int list_size = reg_index_map.size();

  key_list.clear();
  key_list.reserve(list_size);

  const typename KeyIndexMap::const_iterator map_iter = reg_index_map.begin();

  while (map_iter != reg_index_map.end())
    {
    key_list.push_back(map_iter->first);
    ++map_iter;
    }
  }


template<typename KeyType, typename DataType>
unsigned int Registry<KeyType, DataType>::get_reg_index(const KeyType & key) 
  const
  {
  const typename KeyIndexMap::const_iterator map_iter 
    = reg_index_map.find(key);

  if (map_iter == reg_index_map.end())
    throw Except("Key not found", __FILE__, __LINE__);

  return map_iter->second;
  }


template<typename KeyType, typename DataType>
unsigned int Registry<KeyType, DataType>::add_new(const KeyType & key)
  {
  // See if the key already exists

  if (key_present(key))
    throw Except("Key already present", __FILE__, __LINE__);

  // Create the new object

  const unsigned int new_index = reg.size();

  reg.resize((new_index + 1), 0);  // Creates a null vector at reg[new_index]

  reg[new_index] = new DataType;

  // Add to the registry

  reg_index_map[key] = new_index;

  return new_index;
  }


template<typename KeyType, typename DataType>
unsigned int Registry<KeyType, DataType>::add_new(const KeyType & key, 
						  const DataType & original)
  {
  // See if the key already exists

  if (key_present(key))
    throw Except("Key already present", __FILE__, __LINE__);

  // Create the new object

  const unsigned int new_index = reg.size();

  reg.resize((new_index + 1), 0);  // Creates a null vector at reg[new_index]

  reg[new_index] = new DataType(original);

  // Add to the registry

  reg_index_map[key] = new_index;

  return new_index;
  }


}  // End namespace Registry
}  // End namespace QuickFlash


#endif  // QUICKFLASH_REGISTRY_HPP
