// C++ program file quickflash_counters.cpp

/*
  By Nathan C. Hearn
     October 30, 2006

  Definitions for counter classes.
*/


#include "quickflash_counters.hpp"
#include <vector>


namespace QuickFlash
{

// Class VectorCounter

VectorCounter::VectorCounter() :
  counter_dims(), current_state(), num_states(0), counter_in_bounds(false)
  { reset(); }


VectorCounter::VectorCounter(const std::vector<unsigned int> & limits) :
  counter_dims(), current_state(), num_states(0), counter_in_bounds(false)
  { reset(limits); }


VectorCounter::VectorCounter(const VectorCounter & source) :
  counter_dims(), current_state(), num_states(0), counter_in_bounds(false)
  { reset(source); }


void VectorCounter::reset()
  {
  counter_dims.resize(0);
  current_state.resize(0);

  num_states = 0;

  counter_in_bounds = false;
  }


void VectorCounter::reset(const std::vector<unsigned int> & limits)
  {
  const unsigned int dims = limits.size();

  counter_dims = limits;

  current_state.resize(dims);

  for (unsigned int i = 0; i < dims; i++)
    current_state[i] = 0;

  unsigned int num_counter_states = 1;

  for (unsigned int i = 0; i < dims; i++)
    num_counter_states *= limits[i];

  if ((num_counter_states < 1) || (dims < 1))
    reset();
  else
    {
    counter_in_bounds = true;
    num_states = num_counter_states;
    }
  }


void VectorCounter::reset(const VectorCounter & source)
  {
  if (&source != this)
    {
    counter_dims = source.counter_dims;
    current_state = source.current_state;

    num_states = source.num_states;

    counter_in_bounds = source.counter_in_bounds;
    }
  }


int VectorCounter::increment()
  {
  int ret_val = -1;  // Error state

  // Increment only if counter in bounds

  if (counter_in_bounds)
    {
    /*
      Start with least-sigificant element (relativeCounterState[dims - 1]), 
      and increment it.  If it goes out of bounds, set it to zero; move to
      next more-significant element and repeat.  If most significant element
      goes out of bounds, return an empty vector and set counterInBounds to
      false.
    */

    bool increment_completed = false;

    const unsigned int dims = get_dims();

    unsigned int elem = dims - 1;

    while((!increment_completed) && (counter_in_bounds))
      {
      current_state[elem]++;

      if(current_state[elem] < counter_dims[elem])
	{
	increment_completed = true;
	ret_val = 0;  // OK
	}
      else
	{
	current_state[elem] = 0;  // Roll the odometer over...

	if (elem > 0)
	  elem--;  // Next more-sigificant element
	else
	  counter_in_bounds = false;  // Counter is out of bounds!
	}
      }
    }

  return ret_val;
  }


void VectorCounter::reset_counter()
  {
  // Set all elements to zero

  const unsigned int dims = get_dims();

  if (dims > 0)
    {
    for(unsigned int i = 0; i < dims; i++)
      current_state[i] = 0;

    counter_in_bounds = true;
    }
  }


unsigned int VectorCounter::get_scalar_index() const
  {
  unsigned int scalar_index = 0;

  if (counter_in_bounds)
    {
    /*
      Like a C array, let the lowest order elements have the highest 
      precedence.  This stepping is compatible with the increment operation.
    */

    const unsigned int dims = get_dims();

    scalar_index = current_state[0];

    for(unsigned int i = 1; i < dims; i++)
      {
      scalar_index *= counter_dims[i];
      scalar_index += current_state[i];
      }
    }
  else
    scalar_index = num_states;  // Error condition!!!

  return scalar_index;
  }


unsigned int VectorCounter::get_scalar_index_reverse() const
  {
  unsigned int scalar_index = 0;

  if (counter_in_bounds)
    {
    /*
      UNLIKE a C array, let the HIGHEST order elements have the highest 
      precedence.  This stepping is NOT compatible with the increment 
      operation.

      However, this is the way Fortran arrays are traditionally stored.
      Using this function along with the increment function will perform
      a Fortran to C conversion.
    */

    const unsigned int dims = get_dims();

    scalar_index = current_state[dims - 1];

    for(unsigned int i = (dims - 1); i > 0; i--)
      {
      scalar_index *= counter_dims[i - 1];
      scalar_index += current_state[i - 1];
      }
    }
  else
    scalar_index = num_states;  // Error condition!!!

  return scalar_index;
  }


}  // End namespace QuickFlash
