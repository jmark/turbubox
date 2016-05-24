// C++ auxiliary template file quickflash_cache_indexcache.tcc

/*
  By Nathan C. Hearn
     March 6, 2007

  For inclusion only by quickflash_cache_indexcache.hpp.
*/


#ifndef QUICKFLASH_CACHE_INDEXCACHE_TCC
#define QUICKFLASH_CACHE_INDEXCACHE_TCC

#include "quickflash_cache_indexcache.hpp"

namespace QuickFlash
{
namespace Cache
{

template<typename Item_Type>
void IndexCache<Item_Type>::reset(const IndexCache<Item_Type> & source)
  {
  // NOTE: Does not copy access or creation counters

  if (&source != this)
    {
    if (source.is_initialized())
      {
      const unsigned int num_items = source.get_num_items();

      reset(num_items);

      max_size = 0;
      clean_size = 0;

      for (unsigned int index = 0; index < num_items; index++)
	if (source.data_present(index))
	  add_data(index, source[index]);

      max_size = source.max_size;
      clean_size = source.clean_size;

      if (max_size > 0)
	clean_cache(clean_size, max_size);

      initialized = true;
      }
    else
      reset();
    }
  }


template<typename Item_Type>
void IndexCache<Item_Type>::add_data(const unsigned int index, 
				     const Item_Type & source)
  {
  // NOTE: Creation counters only updated if data is not currently present

  if (!data_present(index))
    {
    if (max_size > 0)
      clean_cache((clean_size - 1), max_size); // Will make space for this item

    // Place the item at the front of the list

    CacheItem empty_item;

    item_list.push_front(empty_item);
    item_list.front().reset(index, source);  // The index is the item ID

    // Make the iterator at index reference the new item

    item_iter_list[index] = item_list.begin();

    // Add to the creation counter

    create_counters[index]++;

    current_size++;
    }
  else
    item_iter_list[index]->set_data(source);
  }


template<typename Item_Type>
void IndexCache<Item_Type>::remove_data(const unsigned int index)
  {
  // Deletes the item with ID index from item_list, and resets its iterator

  if (data_present(index))
    {
    item_list.erase(item_iter_list[index]);
    item_iter_list[index] = item_list.end();

    current_size--;
    }
  }


template<typename Item_Type>
void IndexCache<Item_Type>::clean_cache(const unsigned int keep_item_count,
					const unsigned int max_item_count)
  {
  // Only operate if the cache is too full

  bool cache_needs_cleaning = false;

  if (current_size > keep_item_count)
    {
    if (max_item_count > 0)
      {
      if (current_size > max_item_count)
	cache_needs_cleaning = true;
      }
    else
      cache_needs_cleaning = true;
    }

  if (cache_needs_cleaning)
    {
    clean_counter++;

    // Reduce the size of item_list to clean_size by deleting items at the 
    // end of item_list (and setting their iterators to the list end)

    const unsigned int delete_count = current_size - keep_item_count;

    for (unsigned int counter = 0; counter < delete_count; counter++)
      {
      // Get the index of the last item in item_list

      const unsigned int item_index = item_list.back().get_id();

      // Delete the item

      remove_data(item_index);
      }
    }
  }


template<typename Item_Type>
void IndexCache<Item_Type>::promote_item(const unsigned int index) const
  {
  // Move the iterator for the item at index to the front of the cache list

  item_list.splice(item_list.begin(), item_list, item_iter_list[index]);
  }


template<typename Item_Type>
unsigned int IndexCache<Item_Type>::get_num_created_items() const
  {
  // Determines the number of items that were created at least once

  unsigned int ret_val = 0;

  const unsigned int num_items = get_num_items();

  for (unsigned int index = 0; index < num_items; index++)
    if (create_counters[index] > 0)
      ret_val++;

  return ret_val;
  }


template<typename Item_Type>
unsigned int IndexCache<Item_Type>::get_num_multi_created_items() const
  {
  // Determines the number of items that were created more than once

  unsigned int ret_val = 0;

  const unsigned int num_items = get_num_items();

  for (unsigned int index = 0; index < num_items; index++)
    if (create_counters[index] > 1)
      ret_val++;

  return ret_val;
  }


}  // End namespace Cache
}  // End namespace QuickFlash


#endif  // QUICKFLASH_CACHE_INDEXCACHE_TCC
