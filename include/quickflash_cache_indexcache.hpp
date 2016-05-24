// C++ header file quickflash_cache_indexcache.hpp

/*
  By Nathan C. Hearn
     October 26, 2006

  IndexCache class template.

  To do: Add exception handling for non-initialized access???
*/


#ifndef QUICKFLASH_CACHE_INDEXCACHE_HPP
#define QUICKFLASH_CACHE_INDEXCACHE_HPP


#include <vector>
#include <list>
#include <memory>  // std::pair
#include <algorithm>  // std::sort
#include "quickflash_cache_item.hpp"


namespace QuickFlash
{
namespace Cache
{

// Class IndexCache

template<typename Item_Type>
class IndexCache
  {
  private:
    typedef Item<Item_Type> CacheItem;
    typedef std::list<CacheItem> CacheItemList;

  public:
    IndexCache() : 
      initialized(false), max_size(0), clean_size(0), current_size(0),
      access_counters(), create_counters(), item_iter_list(), item_list(), 
      clean_counter(0)
      { reset(); }

    IndexCache(const unsigned int num_items, 
	       const unsigned int max_cache_size=0,
	       const unsigned int clean_cache_size=0) :
      initialized(false), max_size(0), clean_size(0), current_size(0),
      access_counters(), create_counters(), item_iter_list(), item_list(), 
      clean_counter(0)
      { reset(num_items, max_cache_size, clean_cache_size); }

    IndexCache(const IndexCache<Item_Type> & source) :
      initialized(false), max_size(0), clean_size(0), current_size(0),
      access_counters(), create_counters(), item_iter_list(), item_list(), 
      clean_counter(0)
      { reset(source); }

    ~IndexCache() { reset(); }

    IndexCache<Item_Type> & operator=(const IndexCache<Item_Type> & source)
      {
      reset(source);
      return *this;
      }

    void reset() 
      {
      flush();

      max_size = 0;
      clean_size = 0;

      current_size = 0;

      access_counters.clear();
      create_counters.clear();

      item_iter_list.clear();

      clean_counter = 0;

      initialized = false;
      }

    void reset(const unsigned int num_items, 
	       const unsigned int max_cache_size=0,
	       const unsigned int clean_cache_size=0)
      {
      reset();

      set_max_cache_size(max_cache_size, clean_cache_size);

      access_counters.resize(num_items, 0);
      create_counters.resize(num_items, 0);

      item_iter_list.resize(num_items, item_list.end());

      clean_counter = 0;

      initialized = true;
      }

    void reset(const IndexCache<Item_Type> & source);

    void reset_access_counters()
      {
      const unsigned int num_items = get_num_items();

      for (unsigned int index = 0; index < num_items; index++)
	access_counters[index] = 0;
      }

    void reset_create_counters()
      {
      const unsigned int num_items = get_num_items();

      for (unsigned int index = 0; index < num_items; index++)
	create_counters[index] = 0;
      }

    bool is_initialized() const { return initialized; }

    void set_max_cache_size(const unsigned int max_cache_size,
			    const unsigned int clean_cache_size=0)
      { 
      max_size = max_cache_size;
      clean_size = clean_cache_size;

      if ((clean_size > max_size) || (clean_cache_size == 0))
	clean_size = max_size;

      if (max_size > 0)
	clean_cache(clean_size, max_size);
      }

    unsigned int get_clean_counter() const { return clean_counter; }

    unsigned int get_max_cache_size() const { return max_size; }

    unsigned int get_num_items() const { return item_iter_list.size(); }

    unsigned int get_num_cache_items() const { return current_size; }

    void add_data(const unsigned int index, const Item_Type & source);

    void remove_data(const unsigned int index);

    void flush(const unsigned int keep_item_count=0)
      { clean_cache(keep_item_count); }

    bool data_present(const unsigned int index) const
      { return (item_iter_list[index] != item_list.end()); }

    Item_Type & get_item(const unsigned int index) 
      {
      access_counters[index]++;
      promote_item(index);
      return item_iter_list[index]->get_data();
      }

    const Item_Type & get_item(const unsigned int index) const
      {
      access_counters[index]++;
      promote_item(index);
      return item_iter_list[index]->get_data();
      }

    Item_Type & operator[](const unsigned int index) 
      { return get_item(index); }

    const Item_Type & operator[](const unsigned int index) const
      { return get_item(index); }

    unsigned int get_access_counter(const unsigned int index) const
      { return access_counters[index]; }

    unsigned int get_create_counter(const unsigned int index) const
      { return create_counters[index]; }

    unsigned int get_num_created_items() const;
    unsigned int get_num_multi_created_items() const;      

  private:
    void clean_cache(const unsigned int keep_item_count,
		     const unsigned int max_item_count=0);

    void promote_item(const unsigned int index) const;

  private:
    bool initialized;

    unsigned int max_size;
    unsigned int clean_size;

    unsigned int current_size;

    mutable std::vector<unsigned int> access_counters;
    std::vector<unsigned int> create_counters;

    std::vector<typename CacheItemList::iterator> item_iter_list;

    mutable CacheItemList item_list;

    unsigned int clean_counter;
  };
    


}  // End namespace Cache
}  // End namespace QuickFlash


// Include the function definitions

#include "quickflash_cache_indexcache.tcc"


#endif  // QUICKFLASH_CACHE_INDEXCACHE_HPP
