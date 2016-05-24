// C++ header file quickflash_cache_item.hpp

/*
  By Nathan C. Hearn
     October 26, 2006

  Template for cache storage items.
*/


#ifndef QUICKFLASH_CACHE_ITEM_HPP
#define QUICKFLASH_CACHE_ITEM_HPP


namespace QuickFlash
{
namespace Cache
{

// Class Item

template<typename Item_Type>
class Item
  {
  public:
    Item() : initialized(false), id(0), data() { }

    Item(const unsigned int item_id, const Item_Type & item_data_source) :
      initialized(false), id(0), data()
      { reset(item_id, item_data_source); }

    Item(const Item<Item_Type> & source) :
      initialized(false), id(0), data()
      { reset(source); }

    ~Item() { }

    Item<Item_Type> & operator=(const Item<Item_Type> & source)
      {
      reset(source);
      return *this;
      }

    bool set_up() const { return initialized; }

    unsigned int get_id() const { return id; }

    void reset()
      {
      initialized = false;

      id = 0;
      data.reset();
      }
      
    void reset(const unsigned int item_id, const Item_Type & item_data_source)
      {
      initialized = true;

      id = item_id;
      data = item_data_source;
      }

    void reset(const Item<Item_Type> & source)
      {
      if (&source != this)
	{
	initialized = source.initialized;

	id = source.id;
	data = source.data;
	}
      }

    void set_data(const Item_Type & item_data_source)
      {
      // NOTE: Does not change initialized status!!  (No ID if not initialized)

      data = item_data_source;
      }

    Item_Type & get_data() { return data; }
    const Item_Type & get_data() const { return data; }

    Item_Type * get_data_ptr() { return (&data); }
    const Item_Type * get_data_ptr() const { return (&data); }

  private:
    bool initialized;

    unsigned int id;
    Item_Type data;
  };
	    

}  // End namespace Cache
}  // End namespace QuickFlash


#endif  // QUICKFLASH_CACHE_ITEM_HPP
